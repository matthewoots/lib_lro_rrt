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
        
        global_to_vector_transform.rotate(q.inverse());
        global_to_vector_transform.translate(-start_end.first);

        // Initialize variables for RRT
        nodes.clear();
        Node start_node;
        Node end_node;

        start_node.position = start_end.first;
        start_node.parent = NULL;

        nodes.push_back(&start_node);
        end_node.position = start_end.second;

        // The buffer for the xyz search area determined by _sensor_range
        double buffer_factor = 2.0;
        _sensor_buffer = buffer_factor * _sensor_range;

        Eigen::Vector3d transformed_start = transform_vector_with_affine(
            start_end.first, global_to_vector_transform);
        Eigen::Vector3d transformed_end = transform_vector_with_affine(
            start_end.second, global_to_vector_transform);

        std::cout << "transformed_start = " << transformed_start.transpose() << 
            " transformed_end = " << transformed_end.transpose() << std::endl;

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
                search_single_node(start_node, end_node);
                // std::cout << "    search time = " << duration<double>(system_clock::now() - search_timer).count()*1000 << "ms" << std::endl;
                if (duration<double>(system_clock::now() - fail_sub_timer).count() > _runtime_error.first)
                {
                    break;
                }
            }
            
            if (reached)
                break;
        }
        std::cout << (reached ? "Successful" : "Unsuccessful") << " search complete after " << 
            duration<double>(system_clock::now() - fail_timer).count() << "s" << std::endl;
        
        if (!reached)
            return std::vector<Vector3d>();

        std::cout << "Search complete with " << nodes.size() << " nodes" << std::endl;

        return path_extraction(start_node, end_node);
    }

    void tbborrt_server_node::search_single_node(Node start, Node end)
    {
        std::mt19937 generator(dev());
        
        // Setup bounds
        std::uniform_real_distribution<double> dis_middle(-_sensor_buffer, _sensor_buffer);

        Node* step_node = new Node;
        Eigen::Vector3d transformed_random_vector = 
            Eigen::Vector3d(dis_middle(generator), dis_middle(generator), dis_middle(generator));

        double transformed_random_vector_distance = transformed_random_vector.norm();
        
        Eigen::Vector3d random_vector = transform_vector_with_affine(
            transformed_random_vector, global_to_vector_transform.inverse());

        step_node->position = random_vector;

        int index = get_nearest_node(*step_node, start);

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
            check_line_validity(end.position, step_node->position))
        {
            reached = true;
            end.parent = step_node;
            nodes.push_back(&end);
            (nodes[nodes.size()-1]->children).push_back(&end);
            return;
        }
        else
            return;

        iteration++;
    }

    bool tbborrt_server_node::check_line_validity(
        Eigen::Vector3d p, Eigen::Vector3d q)
    {
        float precision = (float)_resolution;
        pcl::PointCloud<pcl::PointXYZ>::VectorType voxels_in_line_search;
        Eigen::Vector3f p_f = Eigen::Vector3f((float)p.x(), (float)p.y(), (float)p.z());
        Eigen::Vector3f q_f = Eigen::Vector3f((float)q.x(), (float)q.y(), (float)q.z());
        int voxels = (int)_octree.getApproxIntersectedVoxelCentersBySegment(p_f, q_f, voxels_in_line_search, precision);
        int intersects = 0;
        for (int i = 0; i < voxels_in_line_search.size(); i++)
        {
            if (_octree.isVoxelOccupiedAtPoint(voxels_in_line_search[i]))
                return false;
        }

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