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
    vector<Eigen::Vector3d> lro_rrt_server_node::find_path(vector<Eigen::Vector3d> previous_input)
    {
        std::lock_guard<std::mutex> octree_lock(octree_mutex);

        // Check if we have valid start and end positions given
        if (search_param.s_e.first.norm() <= 1E-6 && search_param.s_e.second.norm() <= 1E-6)
        {
            std::cout << KRED << "Failure finding path return vector.size() = 0" << KNRM << std::endl;
            return std::vector<Vector3d>();
        }

        Eigen::Vector3d translational_difference = search_param.s_e.second - search_param.s_e.first;
        double translational_difference_distance = translational_difference.norm();
        Eigen::Vector3d translational_difference_vec = translational_difference.normalized();
        
        double bearing = tan(translational_difference_vec.y() / translational_difference_vec.x());

        // Initialize variables for RRT (reset)
        nodes.clear();

        start_node.position = search_param.s_e.first;
        start_node.parent = NULL;
        start_node.children.clear();

        nodes.push_back(&start_node);
        end_node.position = search_param.s_e.second;
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
        
        /** @brief Debug message **/
        // std::cout << "pointcloud size (" << obs_pcl->points.size() << ")" << std::endl;
        
        pcl::PointCloud<pcl::PointXYZ>::VectorType _occupied_voxels;
        search_param.o_p = _octree.getOccupiedVoxelCenters(_occupied_voxels);
        
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
        const Eigen::Vector3d origin, const Eigen::Vector3d end, float precision, Eigen::Vector3d& intersect)
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
            {
                intersect.x() = point.x;
                intersect.y() = point.y;
                intersect.z() = point.z;
                return false;
            }
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
            {
                intersect.x() = point.x;
                intersect.y() = point.y;
                intersect.z() = point.z;
                return false;
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
    void lro_rrt_server_node::search_single_node()
    {
        std::mt19937 generator(dev());
        
        // Setup bounds
        std::uniform_real_distribution<double> dis(-1.0, 1.0);
        std::uniform_real_distribution<double> dis_off(0.2, 1.0);
        std::uniform_real_distribution<double> dis_height(param.h_c.first, param.h_c.second);

        Node* step_node = new Node;

        
        Eigen::Vector3d transformed_random_vector, random_vector;

        if (search_param.o_p > 0) // When there are points in the cloud
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

    /** @brief Reorder then shorten the path by finding shortest obstacle free path through the nodes **/
    std::vector<Eigen::Vector3d> 
    lro_rrt_server_node::path_extraction()
    {
        Node up, down;
        down = end_node;
        up = *(end_node.parent);
        std::vector<Eigen::Vector3d> path;

        while(1)
        {
            path.push_back(down.position);
            if(up.parent == NULL)
                break;
            up = *(up.parent);
            down = *(down.parent);
        }

        std::vector<Eigen::Vector3d> reordered_path = get_reorder_path(path);
        std::vector<Eigen::Vector3d> shortened_path = get_shorten_path(reordered_path);

        for (int i = 0; i < (int)reordered_path.size(); i++)
            std::cout << KCYN << reordered_path[i].transpose() << KNRM << std::endl;

        return shortened_path;
    }

    /** @brief Reorder the path since the output of RRT is inverted **/
    std::vector<Eigen::Vector3d> 
    lro_rrt_server_node::get_reorder_path(std::vector<Eigen::Vector3d> path)
    {
        std::vector<Eigen::Vector3d> reordered_path;
        reordered_path.push_back(start_node.position);
        for (int i = (int)path.size()-1; i >= 0; i--)
            reordered_path.push_back(path[i]);            

        return reordered_path;
    }

    /** @brief Shorten the RRT path by trimming the nodes **/
    std::vector<Eigen::Vector3d> 
    lro_rrt_server_node::get_shorten_path(std::vector<Eigen::Vector3d> path)
    {
        std::vector<Eigen::Vector3d> shortened_path;
        shortened_path.push_back(path[0]);

        for (int i = 1; i < (int)path.size(); i++)
        {
            if (!check_line_validity(
                shortened_path[(int)shortened_path.size()-1], path[i]))
            {        
                shortened_path.push_back(path[i-1]);
            }
        }   
        shortened_path.push_back(path[path.size()-1]);      

        return shortened_path;
    }

    /** @brief Check if the point is within the octree **/
    inline bool lro_rrt_server_node::point_within_octree(Eigen::Vector3d point)
    {
        // Check octree boundary
        if (point.x() < search_param.mx_b.x() - param.r/2 && 
            point.x() > search_param.mn_b.x() + param.r/2 &&
            point.y() < search_param.mx_b.y() - param.r/2 && 
            point.y() > search_param.mn_b.y() + param.r/2 &&
            point.z() < search_param.mx_b.z() - param.r/2 && 
            point.z() > search_param.mn_b.z() + param.r/2)
            return true;
        else
            return false;
    }

}