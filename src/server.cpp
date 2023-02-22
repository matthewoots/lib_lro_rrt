/*
 * server.cpp
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

#include "server.h"

namespace lro
{
    /** 
     * **********************************
     * **********************************
     * **********************************
     * public functions 
     * **********************************
     * **********************************
     * **********************************
    **/

    /** 
     * @brief get_search_path
     * Main run function of the lro module 
    **/ 
    bool server::get_search_path(
        Eigen::Vector3d &start, Eigen::Vector3d &end,
        vector<Eigen::Vector3d> &output, bool sample_tree)
    {
        std::lock_guard<std::mutex> octree_lock(octree_mutex);

        output.clear(); 

        // Check if we have valid start and end positions given
        if ((start - end).norm() <= 1E-6)
        {
            std::cout << KRED << "[too close] fail to find path, return false" << KNRM << std::endl;
            return false;
        }

        bool valid_p1 = true, valid_p2 = true;
        if (is_point_within_octree(start))
        {
            pcl::PointXYZ point1(
                start.x(), start.y(), start.z());
            valid_p1 = !_octree.isVoxelOccupiedAtPoint(point1);
        }
        
        if (is_point_within_octree(end))
        {
            pcl::PointXYZ point2(
                end.x(), end.y(), end.z());
            valid_p2 = !_octree.isVoxelOccupiedAtPoint(point2);
        }

        // Check whether voxel is occupied
        if (!valid_p1 || !valid_p2)
        {
            std::cout << KRED << "[occupied] fail to find path, return false" << KNRM << std::endl;
            return false;
        }

        start_node = new Node;
        start_node->position = start;
        start_node->parent = NULL;
        start_node->children.clear();
        start_node->cost_from_start = 0.0;

        nodes.push_back(start_node);

        end_node = new Node;
        end_node->position = end;
        end_node->parent = NULL;
        end_node->cost_from_start = FLT_MAX;
        end_node->children.clear();

        // std::cout << "start: " << start.transpose() << 
        //     " end: " << end.transpose() << std::endl;
        
        sampler = {};
        sampler.set_constrained_circle_parameters(
            std::make_pair(start, end), _scaled_min_dist, 
            _sensor_range, _search_limit_hfov_list, 
            _search_limit_vfov_list);

        time_point<std::chrono::system_clock> fail_timer = system_clock::now();
        
        reached = false;

        int iteration = 0;
        // When there are points in the cloud
        if (search_param.o_p != 0) 
        {
            time_point<std::chrono::system_clock> fail_sub_timer = system_clock::now();
            iteration = 0;

            kd_tree = kd_create(3);
            int v = kd_insert3(
                kd_tree, start_node->position.x(), 
                start_node->position.y(), start_node->position.z(),
                start_node);
            
            double time_offset = 0.0;
            bool set_offset = false;
            
            while(1)
            {
                // time_point<std::chrono::system_clock> search_timer = system_clock::now();
                query_single_node();

                /** @brief Debug message **/
                // std::cout << "sub-search time(" << KBLU <<
                //    duration<double>(system_clock::now() - search_timer).count()*1000 << 
                //    KNRM << "ms" << std::endl;
                
                if (reached && !set_offset)
                {
                    time_offset = 
                        _runtime_error - _refinement_time -
                        duration<double>(system_clock::now() - fail_sub_timer).count();
                    set_offset = true;
                }
                
                if (duration<double>(system_clock::now() - 
                    fail_sub_timer).count() + time_offset > _runtime_error)
                    break;
                
                iteration++;
            }

            /** @brief Debug message **/
            std::cout << "[lro] iterations(" << KBLU <<
                iteration << KNRM << ")" << std::endl;
        }
        // When there are no points in the cloud
        else
        {
            reached = true;
            end_node->parent = start_node;
            nodes.push_back(end_node);

            std::cout << "[lro] " << (reached ? "successful" : "unsuccessful") << " search complete after " << 
                duration<double>(system_clock::now() - fail_timer).count()*1000 << "ms" << std::endl;

            output = extract_final_path(end_node, start);

            std::cout << "[lro] " << "intermediate_nodes(" << KBLU << nodes.size() - 2 << KNRM
                ") iterations(" << KBLU << iteration << KNRM << ") path_size(" << 
                KBLU << output.size()-2 << KNRM << ")" << std::endl;

            return true;
        }

        std::cout << "[lro] " << (reached ? "Successful" : "Unsuccessful") << " search complete after " << 
            duration<double>(system_clock::now() - fail_timer).count()*1000 << "ms" << std::endl;

        if (!reached)
        {
            Node *safe_node = new Node;
            std::cout << "[lro] " << KRED << "[use safe path] fail to find path, return false" << KNRM << std::endl;
            // safe_node = get_safe_point_in_tree(start_node, 4.0);
            safe_node = nullptr;

            if (safe_node != nullptr)
                output = extract_final_path(safe_node, start);
            else
                std::cout << "[lro] " << KRED << 
                    "start_node does not have children" << KNRM << std::endl;
                
            std::cout << "[lro] " << KRED << "extract last safe path" << KNRM << std::endl;

            return false;
        }

        nodes.push_back(end_node);

        output = extract_final_path(end_node, start);

        if (sample_tree)
            sample_whole_tree(start_node);

        std::cout << "[lro] " << "intermediate_nodes(" << KBLU << nodes.size() - 2 << KNRM
            ") iterations(" << KBLU << iteration << KNRM << ") path_size(" << 
            KBLU << output.size()-2 << KNRM << ")" << std::endl;

        kd_free(kd_tree);

        return true;
    }

    /** 
     * @brief get_line_validity
     * Check whether the line between the pair of points is 
     * obstacle free 
    **/
    bool server::get_line_validity(Eigen::Vector3d p, Eigen::Vector3d q)
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
        double step = _resolution * 0.9;

        // time_point<std::chrono::system_clock> t_b_t = system_clock::now();

        while (!is_point_within_octree(p_fd))
        {
            if (dist_counter > t_n)
                return true;
            Eigen::Vector3d vector_step = t_d_pq * step;
            // Find new p_f
            p_fd += vector_step;
            dist_counter += step;
        }

        while (!is_point_within_octree(q_fd))
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
        
        Eigen::Vector3d intersect;
        return check_approx_intersection_by_segment(
            p_fd, q_fd, intersect);
        
    }

    /** 
     * @brief get_tree
     * Get the full spanning tree found by the search
    **/
    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> server::get_tree()
    {
        return edges;
    }

    /** 
     * @brief set_no_fly_zone
     * Setup no fly zone limits for the search 
    **/
    void server::set_no_fly_zone(vector<Eigen::Vector4d> no_fly_zone)
    {
        // Clean up the previous data
        search_param.n_f_z.clear();
        // Write to the no fly zone data
        search_param.n_f_z = no_fly_zone;
    }

    /** 
     * @brief update_octree
     * Update the octree, the octree is centered 
     * around the pose due to perception range 
    **/ 
    void server::update_octree(
        pcl::PointCloud<pcl::PointXYZ>::Ptr obs_pcl)
    {
        std::lock_guard<std::mutex> octree_lock(octree_mutex);

        _octree.deleteTree();
        
        _octree.setInputCloud(obs_pcl);
        _octree.addPointsFromInputCloud();
        
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
            (search_param.mx_b.x() - search_param.mn_b.x() - min_value) / _resolution));
        search_param.m_k[1] =
            static_cast<uint8_t>(std::ceil(
            (search_param.mx_b.y() - search_param.mn_b.y() - min_value) / _resolution));
        search_param.m_k[2] =
            static_cast<uint8_t>(std::ceil(
            (search_param.mx_b.z() - search_param.mn_b.z() - min_value) / _resolution));
        
        /** @brief Debug message **/
        // std::cout << "Minimum Boundary = " << KBLU << search_param.mn_b.transpose() << KNRM << " " << 
        //     "Maximum Boundary = " << KBLU << search_param.mx_b.transpose() << KNRM << std::endl;
    }

    /** 
     * @brief check_approx_intersection_by_segment
     * Checks for intersected voxels that contain pointclouds 
     * Edit function from:
     * https://pointclouds.org/documentation/octree__pointcloud_8hpp_source.html#l00269
     * Definition at line 269 of file octree_pointcloud.hpp
    **/ 
    bool server::check_approx_intersection_by_segment(
        const Eigen::Vector3d origin, const Eigen::Vector3d end, 
        Eigen::Vector3d& intersect)
    {
        pcl::octree::OctreeKey origin_key, end_key;
        pcl::PointXYZ octree_origin(
            origin.x(), origin.y(), origin.z());
        pcl::PointXYZ octree_end(
            end.x(), end.y(), end.z());
        gen_octree_key_for_point(octree_origin, origin_key);
        gen_octree_key_for_point(octree_end, end_key);

        Eigen::Vector3i d(
            end_key.x - origin_key.x,
            end_key.y - origin_key.y,
            end_key.z - origin_key.z
        );

        std::vector<Eigen::Vector3i> idx_list;
        get_bresenham_3d_from_origin(
            d.x(), d.y(), d.z(), idx_list);
        
        for (int i = 0; i < (int)idx_list.size(); i++) 
        {            
            pcl::PointXYZ point;
            // Fall back to find the root of the expanded volume
            pcl::octree::OctreeKey query_key(
                origin_key.x + idx_list[i].x(),
                origin_key.y + idx_list[i].y(),
                origin_key.z + idx_list[i].z());
            gen_leaf_node_center_from_octree_key(query_key, point);

            for (Eigen::Vector3i &offset : offset_list)
            {
                pcl::octree::OctreeKey query_point(
                    query_key.x + offset.x(),
                    query_key.y + offset.y(),
                    query_key.z + offset.z()
                );

                if (!is_index_within_octree(query_point))
                    continue;
                
                if (i > 0)
                {
                    if (abs(idx_list[i].x() + offset.x() - idx_list[i-1].x()) <= EXPANSION &&
                    abs(idx_list[i].y() + offset.y() - idx_list[i-1].y()) <= EXPANSION &&
                    abs(idx_list[i].z() + offset.z() - idx_list[i-1].z()) <= EXPANSION)
                        continue;
                }

                pcl::PointXYZ neighbour_point;
                gen_leaf_node_center_from_octree_key(
                    query_point, neighbour_point);
                if (_octree.isVoxelOccupiedAtPoint(neighbour_point))
                {
                    intersect = Eigen::Vector3d(
                        point.x, point.y, point.z);
                    return false;
                }
            }
        }
        
        return true;
    }

    /** 
     * @brief get_path_validity
     * Check the pairs of points in the path and return whether 
     * the path is still valid
    **/ 
    bool server::get_path_validity(
        std::vector<Eigen::Vector3d> path)
    {
        // If it contains just 1 node which is its current point
        if (path.size() < 2)
            return false;

        // Check to see whether the new control point and the previous inputs
        // have any pointclouds lying inside
        for (size_t i = 0; i < path.size()-1; i++)
        {
            if (!get_line_validity(path[i], path[i+1]))
                return false;
        }
        
        return true;
    }

    /** 
     * @brief create_expanded_list
    **/ 
    std::vector<Eigen::Vector3i> 
        server::create_expanded_list(int extra)
    {
        std::vector<Eigen::Vector3i> out;
        // Setup the neighbouring boxes
        for (int i = -extra; i <= extra; i++)
            for (int j = -extra; j <= extra; j++)
                for (int k = -extra; k <= extra; k++)
                    out.push_back(Eigen::Vector3i(i, j, k));

        return out;
    }
    

    /** 
     * @brief check_trajectory_collision
    **/ 
    bool server::check_trajectory_collision(
        std::vector<Eigen::Vector3d> traj, int &index)
    {
        if (traj.empty())
            return true;
        
        pcl::octree::OctreeKey previous_key;
        for (size_t i = 0; i < traj.size(); i++)
        {
            pcl::octree::OctreeKey octree_key;
            pcl::PointXYZ octree_point(
                traj[i].x(), traj[i].y(), traj[i].z());
                
            gen_octree_key_for_point(octree_point, octree_key);

            if (previous_key == octree_key)
                continue;
            else
                previous_key = octree_key;

            for (Eigen::Vector3i &offset : offset_list)
            {
                pcl::octree::OctreeKey query_point(
                    octree_key.x + offset.x(),
                    octree_key.y + offset.y(),
                    octree_key.z + offset.z()
                );

                if (!is_index_within_octree(query_point))
                    continue;
                
                if (i > 0)
                {
                    if (abs(traj[i].x() + offset.x() - traj[i-1].x()) <= EXPANSION &&
                    abs(traj[i].y() + offset.y() - traj[i-1].y()) <= EXPANSION &&
                    abs(traj[i].z() + offset.z() - traj[i-1].z()) <= EXPANSION)
                        continue;
                }

                pcl::PointXYZ neighbour_point;
                gen_leaf_node_center_from_octree_key(
                    query_point, neighbour_point);
                if (_octree.isVoxelOccupiedAtPoint(neighbour_point))
                {
                    index = i;
                    return true;
                }

                // std::cout << "point (" << neighbour_point.x << 
                //     " " << neighbour_point.y << " " << neighbour_point.z << ")" << std::endl;
            }
        }
        return false;
    }

    /** 
     * @brief gen_leaf_node_center_from_octree_key
     * Edited from the protected function for octree
     * void pcl::octree::OctreePointCloud<PointT, LeafContainerT, BranchContainerT, OctreeT>::
     * genLeafNodeCenterFromOctreeKey(const OctreeKey& key, PointT& point) const
    **/
    void server::gen_leaf_node_center_from_octree_key(
        const pcl::octree::OctreeKey key, pcl::PointXYZ& point)
    {
        // define point to leaf node voxel center
        point.x = static_cast<float>(
            (static_cast<double>(key.x) + 0.5f) * 
            _resolution + search_param.mn_b.x());
        point.y = static_cast<float>(
            (static_cast<double>(key.y) + 0.5f) * 
            _resolution + search_param.mn_b.y());
        point.z =static_cast<float>(
            (static_cast<double>(key.z) + 0.5f) * 
            _resolution + search_param.mn_b.z());
    }

    /** 
     * @brief change_node_parent
     * Swap and update the parent node
    **/
    void server::change_node_parent(
        Node* &node, Node* &parent, 
        const double &cost_from_parent)
    {
        // Remove node from its parent's children list
        if (node->parent)
        {
            for (auto i = node->parent->children.begin(); i != node->parent->children.end(); ++i) 
            {
                if (*i == node)
                {
                    node->parent->children.erase(i);
                    break;
                }
            }
        }
        node->parent = parent;
        node->cost_from_parent = cost_from_parent;
        node->cost_from_start = parent->cost_from_start + cost_from_parent;
        parent->children.push_back(node);

        // for all its descedants, change the cost_from_start and tau_from_start;
        Node *descendant(node);
        std::queue<Node*> Q;
        Q.push(descendant);
        while (!Q.empty())
        {
            descendant = Q.front();
            Q.pop();
            for (const auto &leafptr : descendant->children)
            {
                leafptr->cost_from_start = leafptr->cost_from_parent + descendant->cost_from_start;
                Q.push(leafptr);
            }
        }
    }

    /** 
     * **********************************
     * **********************************
     * **********************************
     * private functions 
     * **********************************
     * **********************************
     * **********************************
    **/

    /** 
     * @brief query_single_node
     * Query a random node found by the sampler and 
     * check for its validity 
    **/
    void server::query_single_node()
    {
        Node* step_node = new Node;
        Eigen::Vector3d random_vector;
        while (1)
        {
            if (!reached)
                random_vector = sampler.get_rand_point_in_circle();
            else
                random_vector = sampler.get_rand_point_in_circle();
            
            Node* nearest_node = 
                get_nearest_neighbour(random_vector);

            Eigen::Vector3d v = 
                (random_vector - nearest_node->position).normalized(); 
            
            random_vector = 
                nearest_node->position + v * 1.5;

            // Make sure to clamp the height
            random_vector.z() = 
                min(max(random_vector.z(), _height_range.first), 
                _height_range.second);

            // Check octree boundary, if it is exit from this loop
            if (!is_point_within_octree(random_vector))
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
        
        Node* nearest_node = 
            get_nearest_neighbour(random_vector);
        
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
        double accepted_neighbour_count = 10, count = 0, rejected_count = 0;
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
                !get_line_validity(current_node.node->position, step_node->position))
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

        // if(!get_line_validity(
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

        if (get_line_validity(step_node->position, end_node->position))
        {
            double dist_to_goal = 
                (step_node->position - end_node->position).norm();
            bool is_better_path = end_node->cost_from_start > 
                dist_to_goal + step_node->cost_from_start;
            
            if (is_better_path)
            {
                change_node_parent(end_node, step_node, dist_to_goal);
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
                (step_node->position - current_node.node->position).norm();
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
                    connected = get_line_validity(
                        step_node->position, current_node.node->position);
            
            if (connected)
            {
                // double best_cost_before_rewire = end_node->cost_from_start;
                change_node_parent(
                    current_node.node, step_node, 
                    dist_to_potential_child);
            }
          }
        }
        kd_res_free(neighbours);

    }

    /** 
     * @brief extract_final_path
     * Find the shortest obstacle free path through the nodes 
     * when the goal node is reached
    **/
    std::vector<Eigen::Vector3d> 
        server::extract_final_path(Node *n, Eigen::Vector3d &start)
    {
        Node up, down;
        down = *n;
        up = *n->parent;
        std::vector<Eigen::Vector3d> path;

        while(1)
        {
            path.push_back(down.position);
            if(up.parent == NULL)
                break;
            up = *(up.parent);
            down = *(down.parent);
        }
        path.push_back(start);

        std::vector<Eigen::Vector3d> reordered_path;
        for (int i = (int)path.size()-1; i >= 0; i--)
            reordered_path.push_back(path[i]); 
            
        return reordered_path;
    }

    void server::sample_whole_tree(Node* &root)
    {
        vertices.clear();
        edges.clear();

        // whatever dfs or bfs
        Node *node(root);
        std::queue<Node*> Q;
        Q.push(node);
        while (!Q.empty())
        {
            node = Q.front();
            Q.pop();
            for (const auto &leafptr : node->children)
            {
                vertices.push_back(leafptr->position);
                edges.emplace_back(
                    std::make_pair(node->position, leafptr->position));
                Q.push(leafptr);
            }
        }
    }

    Node* server::get_safe_point_in_tree(
        Node* root, double distance)
    {
        if (root == nullptr)
            return nullptr;
        
        if ((root->children).empty())
        {
            std::cout << KGRN << "no children" << KNRM << std::endl;
            return nullptr;
        }
        
        // whatever dfs or bfs
        Node *node(root);
        std::queue<Node*> Q;
        Q.push(node);

        int count = 0;
        while (!Q.empty())
        {
            count++;
            std::cout << KGRN << count << KNRM << std::endl;
            if (node == nullptr)
            {
                Q.pop();
                continue;
            }
            else
            {
                node = Q.front();
                Q.pop();
            }
            
            if(!node->children.empty())
            {
                printf("children %d\n", (int)node->children.size());
                for (const auto &leafptr : node->children)
                {
                    std::cout << KGRN << "leaf of " << count << KNRM << std::endl;            
                    if (leafptr->cost_from_start > distance)
                    {
                        std::cout << KGRN << "found safe leaf" << KNRM << std::endl;
                        return leafptr;
                    }
                    Q.push(leafptr);
                }
            }
        }
    }

    /** 
     * @brief is_point_within_octree
     * Check if the point is within the octree 
    **/
    bool server::is_point_within_octree(
        Eigen::Vector3d point)
    {
        // Check octree boundary
        if (point.x() < search_param.mx_b.x() - _resolution/2 && 
            point.x() > search_param.mn_b.x() + _resolution/2 &&
            point.y() < search_param.mx_b.y() - _resolution/2 && 
            point.y() > search_param.mn_b.y() + _resolution/2 &&
            point.z() < search_param.mx_b.z() - _resolution/2 && 
            point.z() > search_param.mn_b.z() + _resolution/2)
            return true;
        else
            return false;
    }

    /** 
     * @brief is_index_within_octree
     * Check if the index is within the octree 
    **/
    bool server::is_index_within_octree(
        pcl::octree::OctreeKey idx)
    {
        // Check octree boundary
        if (idx.x < 0 && 
            idx.x > search_param.m_k[0] &&
            idx.y < 0 && 
            idx.y > search_param.m_k[1] &&
            idx.z < 0 &&
            idx.z > search_param.m_k[2])
            return false;
        else
            return true;
    }

    /** 
     * @brief gen_octree_key_for_point
     * Edited from the protected function for octree
     * void pcl::octree::OctreePointCloud<PointT, LeafContainerT, BranchContainerT, OctreeT>::
     * genOctreeKeyforPoint(const PointT& point_arg, OctreeKey& key_arg) const
    **/
    void server::gen_octree_key_for_point(
        const pcl::PointXYZ point_arg, pcl::octree::OctreeKey& key_arg)
    {
        // calculate integer key for point coordinates
        key_arg.x = static_cast<uint8_t>((
            point_arg.x - search_param.mn_b.x()) / _resolution);
        key_arg.y = static_cast<uint8_t>((
            point_arg.y - search_param.mn_b.y()) / _resolution);
        key_arg.z = static_cast<uint8_t>((
            point_arg.z - search_param.mn_b.z()) / _resolution);
        
        // assert(key_arg.x <= search_param.m_k[0]);
        // assert(key_arg.y <= search_param.m_k[1]);
        // assert(key_arg.z <= search_param.m_k[2]);
    }

}