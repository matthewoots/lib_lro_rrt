/*
 * lro_rrt_server.cpp
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

#include <lro_rrt_server.h>

namespace lro_rrt_server
{
    /** *********************************
     * **********************************
     * **********************************
     * public functions 
     * **********************************
     * **********************************
     * **********************************
     * **/

    /** @brief Main run module for the RRT server **/ 
    bool lro_rrt_server_node::get_path(
        vector<Eigen::Vector3d> previous_input,
        vector<Eigen::Vector3d> &output)
    {
        std::lock_guard<std::mutex> octree_lock(octree_mutex);

        output.clear();

        // Check if we have valid start and end positions given
        if (search_param.s_e.first.norm() <= 1E-6 && search_param.s_e.second.norm() <= 1E-6)
        {
            std::cout << KRED << "fail to find path, return false" << KNRM << std::endl;
            return false;
        }

        // Initialize variables for RRT (reset)
        nodes.clear();

        start_node = new Node;
        start_node->position = search_param.s_e.first;
        start_node->parent = NULL;
        start_node->children.clear();
        start_node->cost_from_start = 0.0;

        nodes.push_back(start_node);

        end_node = new Node;
        end_node->position = search_param.s_e.second;
        end_node->parent = NULL;
        end_node->cost_from_start = FLT_MAX;
        end_node->children.clear();
        
        sampler = {};
        sampler.set_constrained_circle_parameters(
            search_param.s_e, param.s_d_n, param.s_b, param.s_l_h, param.s_l_v);

        time_point<std::chrono::system_clock> fail_timer = system_clock::now();
        
        reached = false;

        int iteration = 0;
        // When there are points in the cloud
        if (search_param.o_p != 0) 
        {
            while (duration<double>(system_clock::now() - fail_timer).count() < param.r_e.second)
            {
                time_point<std::chrono::system_clock> fail_sub_timer = system_clock::now();
                iteration = 0;

                kd_free(kd_tree);
                kd_tree = kd_create(3);
                int v = kd_insert3(
                    kd_tree, start_node->position.x(), 
                    start_node->position.y(), start_node->position.z(),
                    start_node);
                
                // while(!reached)
                while(1)
                {
                    // time_point<std::chrono::system_clock> search_timer = system_clock::now();
                    query_single_node();

                    /** @brief Debug message **/
                    // std::cout << "sub-search time(" << KBLU <<
                    //    duration<double>(system_clock::now() - search_timer).count()*1000 << 
                    //    KNRM << "ms" << std::endl;
                    
                    if (duration<double>(system_clock::now() - fail_sub_timer).count() > param.r_e.first)
                        break;
                    iteration++;
                }

                /** @brief Debug message **/
                std::cout << "sub_iterations(" << KBLU <<
                    iteration << KNRM << ")" << std::endl;
                
                if (reached)
                    break;
            }
        }
        // When there are no points in the cloud
        else
        {
            reached = true;
            end_node->parent = start_node;
            nodes.push_back(end_node);
            // (end_node->children).push_back(end_node);
        }
        std::cout << (reached ? "Successful" : "Unsuccessful") << " search complete after " << 
            duration<double>(system_clock::now() - fail_timer).count()*1000 << "ms" << std::endl;
        
        if (!reached)
        {
            std::cout << KRED << "fail to find path, return false" << KNRM << std::endl;
            return false;
        }
        
        nodes.push_back(end_node);

        output = path_extraction();

        std::cout << "intermediate_nodes(" << KBLU << nodes.size() - 2 << KNRM
            ") iterations(" << KBLU << iteration << KNRM << ") path_size(" << 
            KBLU << output.size()-2 << KNRM << ")" << std::endl;

        return true;
    }

    /** @brief Check whether the line between the pair of points is obstacle free **/
    bool lro_rrt_server_node::check_line_validity(Eigen::Vector3d p, Eigen::Vector3d q)
    {
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
        double step = param.r * 0.9;

        // time_point<std::chrono::system_clock> t_b_t = system_clock::now();

        while (!point_within_aabb(
            p_fd, search_param.mn_b, search_param.mx_b))
        {
            if (dist_counter > t_n)
                return true;
            Eigen::Vector3d vector_step = t_d_pq * step;
            // Find new p_f
            p_fd += vector_step;
            dist_counter += step;
        }

        while (!point_within_aabb(
            q_fd, search_param.mn_b, search_param.mx_b))
        {
            if (dist_counter > t_n)
                return true;
            Eigen::Vector3d vector_step = t_d_qp * step;
            // Find new q_f
            q_fd += vector_step;
            dist_counter += step;
        }

        // bool l1 = false;
        // if (check_line_box(search_param.mn_b, search_param.mx_b, 
        //     p, q, p_fd))
        //     p_fd += t_d_pq * param.r;
        // else
        //     l1 = true;
        
        // if (check_line_box(search_param.mn_b, search_param.mx_b, 
        //     q, p, q_fd))
        //     q_fd += t_d_qp * param.r;
        // else
        //     if (l1) return true;

        // std::cout << "touch_boundary_time = " << KGRN <<
        //     duration<double>(system_clock::now() - 
        //     t_b_t).count()*1000 << "ms" << KNRM << std::endl;

        if ((q_fd - p_fd).norm() < step)
            return true;

        // pcl::PointCloud<pcl::PointXYZ>::VectorType voxels_in_line_search;
        // int voxels = (int)_octree.getApproxIntersectedVoxelCentersBySegment(
        //         p_fd, q_fd, voxels_in_line_search, (float)step);
        
        Eigen::Vector3d intersect;
        return check_approx_intersection_by_segment(
            p_fd, q_fd, (float)step, intersect);
        
    }

    /** @brief Setup the parameters for the search **/ 
    void lro_rrt_server_node::set_parameters(parameters parameter)
    {
        // Clean up the previous data
        param = parameters();
        // Write to the parameter data
        param = parameter;
        // Set the resolution for the octree
        _octree.setResolution(param.r);
        // The buffer for the xyz search area determined by the sensor range
        param.s_b = param.s_bf * param.s_r;
    }

    /** @brief Setup no fly zone limits for the search **/ 
    void lro_rrt_server_node::set_no_fly_zone(vector<Eigen::Vector4d> no_fly_zone)
    {
        // Clean up the previous data
        search_param.n_f_z.clear();
        // Write to the no fly zone data
        search_param.n_f_z = no_fly_zone;
    }

    /** @brief Update the pose and the octree, since the octree is centered around the pose due to perception range **/ 
    void lro_rrt_server_node::update_pose_and_octree(
        pcl::PointCloud<pcl::PointXYZ>::Ptr obs_pcl, Eigen::Vector3d p, Eigen::Vector3d q)
    {
        std::lock_guard<std::mutex> octree_lock(octree_mutex);

        search_param.s_e.first = p;
        search_param.s_e.second = q;

        _octree.deleteTree();
        
        _octree.setInputCloud(obs_pcl);
        _octree.addPointsFromInputCloud();

        p_c = obs_pcl;
        
        /** @brief Debug message **/
        // std::cout << "pointcloud size (" << obs_pcl->points.size() << ")" << std::endl;
        
        search_param.o_p = _octree.getOccupiedVoxelCenters(occupied_voxels);
        
        _octree.getBoundingBox(
            search_param.mn_b.x(), search_param.mn_b.y(), search_param.mn_b.z(),
            search_param.mx_b.x(), search_param.mx_b.y(), search_param.mx_b.z());

        // Altered from:
        // pcl::octree::OctreePointCloud<PointT, LeafContainerT, BranchContainerT, OctreeT>::getKeyBitSize()
        const float min_value = std::numeric_limits<float>::epsilon();
        
        search_param.m_k[0] =
            static_cast<uint8_t>(std::ceil(
            (search_param.mx_b.x() - search_param.mn_b.x() - min_value) / param.r));
        search_param.m_k[1] =
            static_cast<uint8_t>(std::ceil(
            (search_param.mx_b.y() - search_param.mn_b.y() - min_value) / param.r));
        search_param.m_k[2] =
            static_cast<uint8_t>(std::ceil(
            (search_param.mx_b.z() - search_param.mn_b.z() - min_value) / param.r));
        
        /** @brief Debug message **/
        // std::cout << "Minimum Boundary = " << KBLU << search_param.mn_b.transpose() << KNRM << " " << 
        //     "Maximum Boundary = " << KBLU << search_param.mx_b.transpose() << KNRM << std::endl;
    }

    // Edit function from:
    // https://pointclouds.org/documentation/octree__pointcloud_8hpp_source.html#l00269
    // Definition at line 269 of file octree_pointcloud.hpp
    bool lro_rrt_server_node::check_approx_intersection_by_segment(
        const Eigen::Vector3d origin, const Eigen::Vector3d end, 
        float precision, Eigen::Vector3d& intersect)
    {
        Eigen::Vector3d direction = end - origin;
        double norm = direction.norm();
        direction.normalize();

        const double step_size = precision;
        // Ensure we get at least one step for the first voxel.
        const auto nsteps = std::max<std::size_t>(1, norm / step_size);
        
        pcl::octree::OctreeKey prev_key;

        vector<Eigen::Vector3i> index;
        // Setup the neighbouring boxes
        for (int i = -1; i <= 1; i++)
            for (int j = -1; j <= 1; j++)
                for (int k = -1; k <= 1; k++)
                    index.push_back(Eigen::Vector3i(i, j, k));
        
        // Walk along the line segment with small steps, To include the last point + 1
        for (std::size_t i = 0; i < nsteps + 1; i++) 
        {
            Eigen::Vector3d p = origin + (direction * step_size * static_cast<float>(i));
        
            pcl::PointXYZ octree_p(p.x(), p.y(), p.z());
        
            pcl::octree::OctreeKey main_key;
            gen_octree_key_for_point(octree_p, main_key);
            // Not a new key, still the same voxel.
            if ((main_key == prev_key))
                continue;
        
            prev_key = main_key;

            // Check for surrounding voxels if the radius touches them
            for (Eigen::Vector3i &offset : index)
            {
                pcl::PointXYZ point;
                pcl::octree::OctreeKey query_key(
                    main_key.x + offset.x(),
                    main_key.y + offset.y(),
                    main_key.z + offset.z());
                gen_leaf_node_center_from_octree_key(query_key, point);
                
                Eigen::Vector3d point_eigen(point.x, point.y, point.z);
                if ((point_eigen - p).norm() > param.r * 1.75)
                    continue;

                if (_octree.isVoxelOccupiedAtPoint(point))
                {
                    intersect = Eigen::Vector3d(point.x, point.y, point.z);
                    return false;
                }
            }
        }
        
        return true;
    }

    /** *********************************
     * **********************************
     * **********************************
     * private functions 
     * **********************************
     * **********************************
     * **********************************
     * **/

    /** @brief get_nearest_node is responsible for finding the nearest node in the tree
     * for a particular random node **/ 
    inline int lro_rrt_server_node::get_nearest_node(Node random, Node base_node)
    {
        // We give dist a default value if total node is 1 it will fall back on this
        double dist = (base_node.position - random.position).norm();
        double min_dist = dist;
        int linking_node = 0;

        for(int i = 0; i < (int)nodes.size(); i++)
        {
            // Other nodes than start node
            dist = (nodes[i]->position - random.position).norm();
            // Evaluate distance
            if(dist < min_dist)
            {
                min_dist = dist;
                linking_node = i;
            }
        }
        return linking_node;
    }

    /** @brief Adds in a random node and check for the validity **/
    void lro_rrt_server_node::query_single_node()
    {
        double eps = 0.00001;

        Node* step_node = new Node;
        Eigen::Vector3d random_vector;
        while (1)
        {
            random_vector = sampler.get_rand_point_in_circle();
            // Make sure to clamp the height
            random_vector.z() = 
                min(max(random_vector.z(), param.h_c.first), param.h_c.second);

            // Check octree boundary, if it is exit from this loop
            if (!point_within_aabb(
                random_vector, search_param.mn_b, search_param.mx_b))
                break;
            
            pcl::PointXYZ point(random_vector.x(), random_vector.y(), random_vector.z());

            // Check whether voxel is occupied
            if (!_octree.isVoxelOccupiedAtPoint(point))
                break;
        }

        for (int i = 0; i < (int)search_param.n_f_z.size(); i++)
        {
            // x_min, x_max, y_min, y_max in original frame
            double x_min = search_param.n_f_z[i][0], x_max = search_param.n_f_z[i][1];
            double y_min = search_param.n_f_z[i][2], y_max = search_param.n_f_z[i][3];
            
            // Reject point if it is in no fly zone
            if (random_vector.x() <= x_max && random_vector.x() >= x_min &&
                random_vector.y() <= y_max && random_vector.y() >= y_min)
                return;
        }

        step_node->position = random_vector;
        
        struct kdres *nn = kd_nearest3(
            kd_tree, random_vector.x(), random_vector.y(), random_vector.z());
        if (nn == nullptr)
        {
          std::cout << KRED << "nearest query error" << KNRM << std::endl;
          return;
        }

        Node* nearest_node = (Node*)kd_res_item_data(nn);
        // int index = nearest_node->index;
        kd_res_free(nn);
        // int index = get_nearest_node(*step_node, *start_node);
        double min_distance_to_node = 
            (nearest_node->position - random_vector).norm();
        double min_distance_from_start = 
            nearest_node->cost_from_start + min_distance_to_node + 0.0001;

        struct kdres *neighbours;
        neighbours = kd_nearest_range3(
            kd_tree, step_node->position.x(), 
            step_node->position.y(), step_node->position.z(),
            min_distance_to_node * 2);
        
        if (neighbours == nullptr)
            return;
        
        std::vector<node_status_check> neighbour_nodes;
        double accepted_neighbour_count = 20, count = 0, rejected_count = 0;
        neighbour_nodes.reserve(accepted_neighbour_count);

        while (!kd_res_end(neighbours) && count < accepted_neighbour_count)
        {
            double pos[3];
            node_status_check status_node;
            Node *c_n = (Node*)kd_res_item(neighbours, pos);
            status_node.node = c_n;
            status_node.is_checked = status_node.is_valid = false;
            neighbour_nodes.push_back(status_node);
            // store range query result so that we dont need to query again for rewire;
            kd_res_next(neighbours); // go to next in kd tree range query result
            count++;
        }

        Node *optimal_node = new Node;
        for (auto &current_node : neighbour_nodes)
        {
            double distance_to_node = 
                (current_node.node->position - random_vector).norm();

            double total_distance_from_start = 
                current_node.node->cost_from_start + distance_to_node;

            current_node.is_checked = true;
            if (total_distance_from_start >= min_distance_from_start ||
                !check_line_validity(current_node.node->position, step_node->position))
            {
                rejected_count++;
                continue;
            }
            current_node.is_valid = true;
            min_distance_from_start = total_distance_from_start;
            optimal_node = current_node.node;
        }

        if (rejected_count == (int)neighbour_nodes.size())
            return;

        /** @brief Debug message **/
        // std::cout << "Random_node = " << random_vector.transpose() << std::endl;

        // if(!check_line_validity(
        //     nodes[index]->position, step_node->position))
        //     return;
        
        // Add the new node into the list and push_back data on children and parent
        // step_node->parent = nodes[index];
        step_node->parent = optimal_node;
        step_node->cost_from_start = min_distance_from_start;
        step_node->cost_from_parent = (optimal_node->position - random_vector).norm();
        optimal_node->children.push_back(step_node);
        nodes.push_back(step_node);
        // nodes[optimal_node->index]->children.push_back(step_node);
        
        kd_insert3(
            kd_tree, step_node->position.x(), 
            step_node->position.y(), step_node->position.z(), 
            step_node);

        if (check_line_validity(step_node->position, end_node->position))
        {
            double dist_to_goal = 
                (step_node->position - end_node->position).norm();
            bool is_better_path = end_node->cost_from_start > 
                dist_to_goal + step_node->cost_from_start;
            
            if (is_better_path)
            {
                change_node_parent(end_node, step_node, dist_to_goal);
                std::vector<Eigen::Vector3d> path = path_extraction();
                reached = true;
            }
            return;
        }
        else
            return;

        // Conduct rewiring
        for (auto &current_node : neighbour_nodes)
        {
            double dist_to_potential_child = 
                (step_node->position, current_node.node->position).norm();
            bool not_consistent = 
                step_node->cost_from_start + dist_to_potential_child < 
                current_node.node->cost_from_start ? true : false;
            bool promising = 
                step_node->cost_from_start + dist_to_potential_child + 
                (current_node.node->position + step_node->position).norm() < 
                end_node->cost_from_start ? true : false;

            if (not_consistent && promising)
            {
                bool connected(false);
                if (current_node.is_checked)
                    connected = current_node.is_valid;
                else
                    connected = check_line_validity(
                        step_node->position, current_node.node->position);
            
            if (connected)
            {
                double best_cost_before_rewire = end_node->cost_from_start;
                change_node_parent(
                    current_node.node, step_node, 
                    dist_to_potential_child);
            }
          }
        }
    }

    /** @brief Reorder then shorten the path by finding shortest obstacle free path through the nodes **/
    std::vector<Eigen::Vector3d> 
        lro_rrt_server_node::path_extraction()
    {
        Node up, down;
        down = *end_node;
        up = *end_node->parent;
        std::vector<Eigen::Vector3d> path;

        while(1)
        {
            path.push_back(down.position);
            if(up.parent == NULL)
                break;
            up = *(up.parent);
            down = *(down.parent);
        }
        path.push_back(search_param.s_e.first);

        std::vector<Eigen::Vector3d> reordered_path = 
            get_reorder_path(path);
            
        return reordered_path;

        // for (int i = 0; i < (int)reordered_path.size(); i++)
        //     std::cout << KCYN << reordered_path[i].transpose() << KNRM << std::endl;
    }

    /** @brief Reorder the path since the output of RRT is inverted **/
    std::vector<Eigen::Vector3d> 
    lro_rrt_server_node::get_reorder_path(std::vector<Eigen::Vector3d> path)
    {
        std::vector<Eigen::Vector3d> reordered_path;
        for (int i = (int)path.size()-1; i >= 0; i--)
            reordered_path.push_back(path[i]);            

        return reordered_path;
    }

    /** @brief Shorten the RRT path by trimming the nodes **/
    // std::vector<Eigen::Vector3d> 
    // lro_rrt_server_node::get_shorten_path(std::vector<Eigen::Vector3d> path)
    // {
    //     std::vector<Eigen::Vector3d> shortened_path;
    //     shortened_path.push_back(path[0]);

    //     for (int i = 1; i < (int)path.size(); i++)
    //     {
    //         if (!check_line_validity(
    //             shortened_path[(int)shortened_path.size()-1], path[i]))
    //         {        
    //             shortened_path.push_back(path[i-1]);
    //         }
    //     }   
    //     shortened_path.push_back(path[path.size()-1]);      

    //     return shortened_path;
    // }

    /** @brief Check if the point is within the AABB **/
    bool lro_rrt_server_node::point_within_aabb(
        Eigen::Vector3d point, Eigen::Vector3d min,
        Eigen::Vector3d max)
    {
        // Check octree boundary
        if (point.x() < max.x() - param.r/2 && 
            point.x() > min.x() + param.r/2 &&
            point.y() < max.y() - param.r/2 && 
            point.y() > min.y() + param.r/2 &&
            point.z() < max.z() - param.r/2 && 
            point.z() > min.z() + param.r/2)
            return true;
        else
            return false;
    }

    /** @brief Constrain angle to between -pi to pi **/
    double lro_rrt_server_node::constrain_between_180(double x)
    {
        x = fmod(x + M_PI,2*M_PI);
        if (x < 0)
            x += 2*M_PI;
        return x - M_PI;
    }

}