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
    vector<Eigen::Vector3d> tbborrt_server_node::find_path(vector<Eigen::Vector3d> previous_input)
    {
        std::lock_guard<std::mutex> octree_lock(octree_mutex);

        // Check if we have valid start and end positions given
        if (param.s_e.first.norm() <= 1E-6 && param.s_e.second.norm() <= 1E-6)
        {
            std::cout << KRED << "Failure finding path return vector.size() = 0" << KNRM << std::endl;
            return std::vector<Vector3d>();
        }

        Eigen::Vector3d translational_difference = param.s_e.second - param.s_e.first;
        double translational_difference_distance = translational_difference.norm();
        Eigen::Vector3d translational_difference_vec = translational_difference.normalized();
        
        double bearing = tan(translational_difference_vec.y() / translational_difference_vec.x());

        // Initialize variables for RRT (reset)
        nodes.clear();

        start_node.position = param.s_e.first;
        start_node.parent = NULL;
        start_node.children.clear();

        nodes.push_back(&start_node);
        end_node.position = param.s_e.second;
        end_node.parent = NULL;
        end_node.children.clear();
        
        time_point<std::chrono::system_clock> fail_timer = system_clock::now();
        
        reached = false;

        while(duration<double>(system_clock::now() - fail_timer).count() < param.r_e.second)
        {
            time_point<std::chrono::system_clock> fail_sub_timer = system_clock::now();
            iteration = 0;

            while(!reached)
            {
                // time_point<std::chrono::system_clock> search_timer = system_clock::now();
                search_single_node();

                /** @brief Debug message **/
                // std::cout << "    search time = " << duration<double>(system_clock::now() - search_timer).count()*1000 << "ms" << std::endl;
                
                if (duration<double>(system_clock::now() - fail_sub_timer).count() > param.r_e.first)
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
        {
            std::cout << KRED << "Failure finding path return vector.size() = 0" << KNRM << std::endl;
            return std::vector<Vector3d>();
        }

        std::cout << "Search complete with " << nodes.size() << " nodes" << std::endl;

        return path_extraction();
    }

    void tbborrt_server_node::search_single_node()
    {
        std::mt19937 generator(dev());
        
        // Setup bounds
        std::uniform_real_distribution<double> dis(-1.0, 1.0);
        std::uniform_real_distribution<double> dis_off(0.2, 1.0);
        std::uniform_real_distribution<double> dis_height(param.h_c.first, param.h_c.second);

        Node* step_node = new Node;

        
        Eigen::Vector3d transformed_random_vector, random_vector;

        if (_occupied_points > 0) // When there are points in the cloud
        {
            while (1)
            {
                // https://mathworld.wolfram.com/SphericalCoordinates.html
                // https://karthikkaranth.me/blog/generating-random-points-in-a-sphere/
                pcl::PointXYZ point;
                double theta = dis(generator) * M_PI;
                double phi = abs(dis(generator)) * M_PI;
                double r = dis_off(generator) * param.s_b;
                double sin_theta = sin(theta); 
                double cos_theta = cos(theta);
                double sin_phi = sin(phi); 
                double cos_phi = cos(phi);
                // Make sure to clamp the height
                random_vector = 
                    Eigen::Vector3d(start_node.position.x() + (r * sin_phi * cos_theta), 
                    start_node.position.y() + (r * sin_phi * sin_theta), 
                    min(max(start_node.position.z() + (r * cos_phi), param.h_c.first), param.h_c.second));

                // Check octree boundary, if it is exit from this loop
                if (!point_within_octree(random_vector))
                    break;
                
                point.x = random_vector.x();
                point.y = random_vector.y();
                point.z = random_vector.z();

                // Check whether voxel is occupied
                if (!_octree.isVoxelOccupiedAtPoint(point))
                    break;
            }
        }
        else // When there is no points in the cloud
        {
            // https://mathworld.wolfram.com/SphericalCoordinates.html
            // https://karthikkaranth.me/blog/generating-random-points-in-a-sphere/
            double theta = dis(generator) * M_PI;
            double phi = abs(dis(generator)) * M_PI;
            double r = dis_off(generator) * param.s_b;
            double sin_theta = sin(theta); 
            double cos_theta = cos(theta);
            double sin_phi = sin(phi); 
            double cos_phi = cos(phi);
            // Make sure to clamp the height
            random_vector = 
                Eigen::Vector3d(start_node.position.x() + (r * sin_phi * cos_theta), 
                start_node.position.y() + (r * sin_phi * sin_theta), 
                min(max(start_node.position.z() + (r * cos_phi), param.h_c.first), param.h_c.second));
        }

        step_node->position = random_vector;

        int index = get_nearest_node(*step_node, start_node);

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
            nodes[index]->position, step_node->position);
        
        if(!flag)
            return;
        
        // Add the new node into the list and push_back data on children and parent
        step_node->parent = nodes[index];
        nodes.push_back(step_node);
        nodes[index]->children.push_back(step_node);

        if (check_line_validity(step_node->position, end_node.position))
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
        /** @brief Method 1 **/
        // Get the translational difference p to q
        Eigen::Vector3d t_d = q - p;
        // Get the translational vector p to q
        Eigen::Vector3d t_d_pq = t_d.normalized();
        // Get the translational vector q to p
        Eigen::Vector3d t_d_qp = -t_d_pq;
        // Get the translational norm
        double t_n = t_d.norm();

        Eigen::Vector3d p_fd = p;
        Eigen::Vector3d q_fd = q;
        double dist_counter = 0.0;
        double step = param.r * 2;
        while (!point_within_octree(p_fd))
        {
            if (dist_counter > t_n)
                return true;
            Eigen::Vector3d vector_step = t_d_pq * step;
            // Find new q_f
            p_fd += vector_step;
            dist_counter += step;
        }

        while (!point_within_octree(q_fd))
        {
            if (dist_counter > t_n)
                return true;
            Eigen::Vector3d vector_step = t_d_qp * step;
            // Find new q_f
            q_fd += vector_step;
            dist_counter += step;
        }

        if ((q_fd - p_fd).norm() < step)
            return true;

        // pcl::PointCloud<pcl::PointXYZ>::VectorType voxels_in_line_search;
        // int voxels = (int)_octree.getApproxIntersectedVoxelCentersBySegment(
        //         p_fd, q_fd, voxels_in_line_search, (float)step);
        
        return check_approx_intersection_by_segment(
            p_fd, q_fd, (float)step);

        /** @brief End of Method 1 **/
        

        /** @brief Method 2 **/
        // pcl::PointXYZ search_point, end_point, obstacle_point;
        // double accumulated_distance = 0;
        // search_point.x = p.x();
        // search_point.y = p.y();
        // search_point.z = p.z();

        // end_point.x = p.x();
        // end_point.y = q.y();
        // end_point.z = q.z();

        // double distance_target = (q-p).norm();
        // Eigen::Vector3d direction_vector = (q-p) / (q-p).norm();
        // int tries = 0, max_tries = 30;

        // while (accumulated_distance - distance_target < 0 && tries < max_tries)
        // {
        //     tries++;
            
        //     Eigen::Vector3d distance_travelled;

        //     time_point<std::chrono::system_clock> timer = system_clock::now();
            
        //     int index;
        //     float sq_dist;
            
        //     _octree.approxNearestSearch(search_point, index, sq_dist);

        //     // Get the point in the pointcloud
        //     obstacle_point = (*_store_cloud)[index];
        //     float distance = sqrtf(sq_dist);

        //     // Add the distance travelled to the accumulated distance
        //     distance_travelled = direction_vector * distance;

        //     accumulated_distance += (double)distance; 
        //     std::cout << "    tries = " << tries
        //         << " distance = " << distance << " " 
        //         << " accumulated_distance/distance_target = " << accumulated_distance << " / "
        //         << distance_target <<  " time-taken = " 
        //         << duration<double>(system_clock::now() - timer).count()*1000 << "ms" << std::endl;
            
        //     if (distance < (float)param.r * 1.5f)
        //     {
        //         accumulated_distance = pow(10,6);
        //         std::cout << KRED << "Return on violation" << KNRM << std::endl;
        //         return false;
        //     }

        //     search_point.x = distance_travelled.x() + search_point.x;
        //     search_point.y = distance_travelled.y() + search_point.y;
        //     search_point.z = distance_travelled.z() + search_point.z;
        // }

        // return true;

        /** @brief End of Method 2 **/

        
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

    // Edit function from:
    // https://pointclouds.org/documentation/octree__pointcloud_8hpp_source.html#l00269
    // Definition at line 269 of file octree_pointcloud.hpp
    bool tbborrt_server_node::check_approx_intersection_by_segment(
        const Eigen::Vector3d origin, const Eigen::Vector3d end, float precision)
    {
        Eigen::Vector3d direction = end - origin;
        double norm = direction.norm();
        direction.normalize();

        const double step_size = param.r * precision;
        // Ensure we get at least one step for the first voxel.
        const auto nsteps = std::max<std::size_t>(1, norm / step_size);
        
        pcl::octree::OctreeKey prev_key;
        
        bool bkeyDefined = false;
        
        // Walk along the line segment with small steps.
        for (std::size_t i = 0; i < nsteps; ++i) {
            Eigen::Vector3d p = origin + (direction * step_size * static_cast<float>(i));
        
            pcl::PointXYZ octree_p;
            octree_p.x = p.x();
            octree_p.y = p.y();
            octree_p.z = p.z();
        
            pcl::octree::OctreeKey key;
            gen_octree_key_for_point(octree_p, key);
        
            // Not a new key, still the same voxel.
            if ((key == prev_key) && (bkeyDefined))
            continue;
        
            prev_key = key;
            bkeyDefined = true;

            pcl::PointXYZ point;
            gen_leaf_node_center_from_octree_key(key, point);
            if (_octree.isVoxelOccupiedAtPoint(point))
                return false;
        }
        
        pcl::octree::OctreeKey end_key;
        pcl::PointXYZ end_p;
        end_p.x = end.x();
        end_p.y = end.y();
        end_p.z = end.z();

        gen_octree_key_for_point(end_p, end_key);
        if (!(end_key == prev_key)) {
            pcl::PointXYZ point;
            gen_leaf_node_center_from_octree_key(end_key, point);
            if (_octree.isVoxelOccupiedAtPoint(point))
                return false;
        }
        
        return true;
    }
}