/*
 * tbborrt_server.cpp
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2022 Matthew (matthewoots at gmail.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 * 
 * 
 */

#define KNRM  "\033[0m"
#define KRED  "\033[31m"
#define KGRN  "\033[32m"
#define KYEL  "\033[33m"
#define KBLU  "\033[34m"
#define KMAG  "\033[35m"
#define KCYN  "\033[36m"
#define KWHT  "\033[37m"

#include <tbborrt_server.h>

namespace tbborrt_server
{
    vector<Eigen::Vector3d> tbborrt_server_node::find_path(
        vector<Eigen::Vector3d> previous_input, std::pair<Eigen::Vector3d, Eigen::Vector3d> start_end)
    {
        std::lock_guard<std::mutex> octree_lock(octree_mutex);

        Eigen::Vector3d translational_difference = start_end.second - start_end.first;
        double translational_difference_distance = translational_difference.norm();
        Eigen::Vector3d translational_difference_vec = 
            translational_difference / translational_difference.norm();
        Eigen::Quaterniond q = quaternion_from_pitch_yaw(
            Vector3d(1,0,0), translational_difference_vec);

        // Setup global_to_vector_transform transform
        global_to_vector_transform = Affine3d::Identity(); 
        
        global_to_vector_transform.rotate(q);
        global_to_vector_transform.translate(-start_end.first);

        // Initialize variables for RRT (reset)
        nodes.clear();

        start_node.position = start_end.first;
        start_node.parent = NULL;
        start_node.children.clear();

        nodes.push_back(&start_node);
        end_node.position = start_end.second;
        end_node.parent = NULL;
        end_node.children.clear();

        // The buffer for the xyz search area determined by _sensor_range
        double buffer_factor = 2.0;
        _sensor_buffer = buffer_factor * _sensor_range;

        Eigen::Vector3d transformed_start = transform_vector_with_affine(
            start_end.first, global_to_vector_transform);
        Eigen::Vector3d transformed_end = transform_vector_with_affine(
            start_end.second, global_to_vector_transform);

        std::cout << "transformed_start = " << KBLU << transformed_start.transpose() << KNRM << 
            " transformed_end = " << KBLU << transformed_end.transpose() << KNRM << std::endl;

        _goal_within_sensory_bounds = translational_difference_distance - _sensor_buffer < 0;
        
        time_point<std::chrono::system_clock> fail_timer = system_clock::now();
        
        reached = false;

        while(duration<double>(system_clock::now() - fail_timer).count() < _runtime_error.second)
        {
            time_point<std::chrono::system_clock> fail_sub_timer = system_clock::now();
            iteration = 0;

            while(!reached)
            {
                time_point<std::chrono::system_clock> search_timer = system_clock::now();
                search_single_node();

                /** @brief Debug message **/
                // std::cout << "    search time = " << duration<double>(system_clock::now() - search_timer).count()*1000 << "ms" << std::endl;
                
                if (duration<double>(system_clock::now() - fail_sub_timer).count() > _runtime_error.first)
                {
                    break;
                }
                iteration++;
            }

            /** @brief Debug message **/
            std::cout << "Sub-iterations taken = " << iteration << std::endl;
            
            if (reached)
                break;
        }
        std::cout << (reached ? "Successful" : "Unsuccessful") << " search complete after " << 
            duration<double>(system_clock::now() - fail_timer).count()*1000 << "ms" << std::endl;
        
        if (!reached)
            return std::vector<Vector3d>();

        std::cout << "Search complete with " << nodes.size() << " nodes" << std::endl;

        return path_extraction();
    }

    void tbborrt_server_node::search_single_node()
    {
        std::mt19937 generator(dev());
        
        // Setup bounds
        std::uniform_real_distribution<double> dis_middle(-_sensor_buffer, _sensor_buffer);

        Node* step_node = new Node;

        
        Eigen::Vector3d transformed_random_vector, random_vector;

        if (_occupied_points > 0) // When there are points in the cloud
        {
            while (1)
            {
                pcl::PointXYZ point;

                transformed_random_vector = 
                    Eigen::Vector3d(dis_middle(generator), dis_middle(generator), dis_middle(generator));
                random_vector = transform_vector_with_affine(
                    transformed_random_vector, global_to_vector_transform.inverse());             

                // Check octree boundary
                if (!point_within_octree(random_vector))
                    break;
                
                point.x = random_vector.x();
                point.y = random_vector.y();
                point.z = random_vector.z();
                    
                if (!_octree.isVoxelOccupiedAtPoint(point))
                    break;
            }
        }
        else // When there is no points in the cloud
        {
            transformed_random_vector = 
                Eigen::Vector3d(dis_middle(generator), dis_middle(generator), dis_middle(generator));
            random_vector = transform_vector_with_affine(
                transformed_random_vector, global_to_vector_transform.inverse());
        }

        double transformed_random_vector_distance = transformed_random_vector.norm();

        step_node->position = random_vector;

        int index = get_nearest_node(*step_node, start_node);

        // Clamp z axis
        random_vector.z() = max(min(random_vector.z(), _height_constrain.second), _height_constrain.first);

        for (int i = 0; i < (int)_no_fly_zone.size(); i++)
        {
            // x_min, x_max, y_min, y_max in original frame
            double x_min = _no_fly_zone[i][0], x_max = _no_fly_zone[i][1];
            double y_min = _no_fly_zone[i][2], y_max = _no_fly_zone[i][3];
            
            // Reject point if it is in no fly zone
            if (random_vector.x() <= x_max && random_vector.x() >= x_min &&
                random_vector.y() <= y_max && random_vector.y() >= y_min)
                return;
                                
        }

        /** @brief Debug message **/
        // std::cout << "Random_node = " << random_vector.transpose() << std::endl;

        bool flag = check_line_validity(
            step_node->position, nodes[index]->position);
        
        if(!flag)
            return;
        
        // Add the new node into the list and push_back data on children and parent
        step_node->parent = nodes[index];
        nodes.push_back(step_node);
        nodes[index]->children.push_back(step_node);

        if ((transformed_random_vector_distance > _sensor_buffer && random_vector.x() > 0) ||
            check_line_validity(end_node.position, step_node->position))
        {
            reached = true;
            end_node.parent = step_node;
            nodes.push_back(&end_node);
            (nodes[nodes.size()-1]->children).push_back(&end_node);
            return;
        }
        else
            return;

        iteration++;
    }

    bool tbborrt_server_node::check_line_validity(
        Eigen::Vector3d p, Eigen::Vector3d q)
    {
        if (_occupied_points == 0)
            return true;

        /** @brief Method 1 **/
        // float precision = (float)_resolution/2;
        // pcl::PointCloud<pcl::PointXYZ>::VectorType voxels_in_line_search;
        // Eigen::Vector3f p_f = Eigen::Vector3f((float)p.x(), (float)p.y(), (float)p.z());
        // Eigen::Vector3f q_f = Eigen::Vector3f((float)q.x(), (float)q.y(), (float)q.z());
        
        // Eigen::Vector3d direction = (q-p) / (q-p).norm();

        // // Make sure that doing getApproxIntersectedVoxelCentersBySegment is within the octree boundary
        // if (!point_within_octree(p))
        // {

        // }

        // if (!point_within_octree(q))
        // {
            
        // }

        // int voxels = (int)_octree.getApproxIntersectedVoxelCentersBySegment(p_f, q_f, voxels_in_line_search, precision);
        // int intersects = 0;

        // /** @brief Debug message **/
        // // std::cout << "    voxel_size = " << voxels_in_line_search.size() << std::endl;
        
        // for (int i = 0; i < voxels_in_line_search.size(); i++)
        // {
        //     // Check octree boundary
        //     if (!point_within_octree(Eigen::Vector3d(
        //         voxels_in_line_search[i].x, 
        //         voxels_in_line_search[i].y, 
        //         voxels_in_line_search[i].z)))
        //         continue;

        //     if (_octree.isVoxelOccupiedAtPoint(voxels_in_line_search[i]))
        //         return false;
        // }
        /** @brief End of Method 1 **/
        

        /** @brief Method 2 **/
        pcl::PointXYZ search_point, end_point, obstacle_point;
        double accumulated_distance = 0;
        search_point.x = p.x();
        search_point.y = p.y();
        search_point.z = p.z();

        end_point.x = p.x();
        end_point.y = q.y();
        end_point.z = q.z();

        double distance_target = (q-p).norm();
        Eigen::Vector3d direction_vector = (q-p) / (q-p).norm();

        while (accumulated_distance - distance_target < 0)
        {
            // K nearest neighbor search
            int K = 1;

            std::vector<int> pointIdxNKNSearch;
            std::vector<float> pointNKNSquaredDistance;

            time_point<std::chrono::system_clock> timer = system_clock::now();
            if (_octree.nearestKSearch (search_point, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
            {
                obstacle_point = (*_store_cloud)[pointIdxNKNSearch[0]];
                float distance = sqrtf(pointNKNSquaredDistance[0]);
                Eigen::Vector3d distance_travelled = direction_vector * distance_target;
                
                search_point.x = distance_travelled.x() + search_point.x;
                search_point.y = distance_travelled.y() + search_point.y;
                search_point.z = distance_travelled.z() + search_point.z;
                 
                accumulated_distance += (double)distance; 

                Eigen::Vector3d d = 
                    Eigen::Vector3d(obstacle_point.x, obstacle_point.y, obstacle_point.z) -
                    Eigen::Vector3d(search_point.x, search_point.y, search_point.z);
                std::cout << "    d.norm() = " << d.norm() << " " 
                    << " accumulated_distance/distance_target = " << accumulated_distance << " / "
                    << distance_target <<  " time-taken = " 
                    << duration<double>(system_clock::now() - timer).count()*1000 << "ms" << std::endl;
                if (d.norm() < _resolution * 1.5)
                    return false;
            }
        }
        /** @brief End of Method 2 **/

        return true;
    }

    // [get_nearest_node] is responsible for finding the nearest node in the tree 
    // for a particular random node. 
    inline int tbborrt_server_node::get_nearest_node(Node random, Node base_node)
    {
        // We give dist a default value if total node is 1 it will fall back on this
        double dist = separation(base_node.position, random.position);
        double min_dist = dist;
        int linking_node = 0;

        for(int i = 0; i < (int)nodes.size(); i++)
        {
            // Other nodes than start node
            dist = separation(nodes[i]->position, random.position);
            // Evaluate distance
            if(dist < min_dist)
            {
                min_dist = dist;
                linking_node = i;
            }
        }
        return linking_node;
    }
}