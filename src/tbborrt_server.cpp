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
    inline vector<Eigen::Vector3d> tbborrt_server_node::find_path(
        vector<Eigen::Vector3d> previous_input, std::pair<Eigen::Vector3d, Eigen::Vector3d> start_end)
    {
        std::lock_guard<std::mutex> octree_lock(octree_mutex);

        Eigen::Vector3d translational_difference = start_end.second - start_end.first;
        Eigen::Vector3d translational_difference_vec = 
            translational_difference / translational_difference.norm();
        Eigen::Quaterniond q = quaternion_from_pitch_yaw(
            Vector3d(1,0,0), translational_difference_vec);

        // Setup global_to_vector_transform transform
        global_to_vector_transform = Affine3d::Identity(); 
        global_to_vector_transform.translate(translational_difference);
        global_to_vector_transform.rotate(q.inverse());

        // Initialize variables for RRT
        nodes.clear();
        Node start_node, end_node;

        start_node.position = start_end.first;
        start_node.parent = NULL;

        nodes.push_back(&start_node);
        end_node.position = start_end.second;

        // The buffer for the xyz search area determined by _sensor_range
        double buffer_factor = 2.0;
        double sensor_buffer = buffer_factor * _sensor_range;
        
        time_point<std::chrono::system_clock> fail_timer_start = system_clock::now();
        
        bool error = false;
        while(!reached)
        {
            search_single_node();
            if (duration<double>(system_clock::now() - fail_timer_start).count() > _runtime_error.second)
            {
                error = true;
                break;
            }
        }

        if (error)
            return std::vector<Vector3d>();

        return path_extraction(start_node, end_node);
    }

    inline void tbborrt_server_node::search_single_node()
    {
        
    }

    inline bool tbborrt_server_node::check_line_validity(
        Eigen::Vector3d p, Eigen::Vector3d q)
    {
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