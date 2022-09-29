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
        Eigen::Quaterniond q = quaternion_from_pitch_yaw(
            Vector3d(1,0,0), translational_difference_vec);

        // Setup global_to_vector_transform transform
        global_to_vector_transform = Eigen::Affine3d::Identity(); 
        // global_to_vector_transform.linear() = q.inverse().toRotationMatrix();
        global_to_vector_transform.rotate(q.inverse());
        global_to_vector_transform.translate(-param.s_e.first);

        // Eigen::Affine3d test_transform = Eigen::Affine3d::Identity(); 
        // test_transform.linear() = q.toRotationMatrix();
        // test_transform.translate(param.s_e.first);
        // test_transform.rotate(q);

        // std::cout << "affineMatrixTest1 = " << std::endl << 
        //     KBLU << global_to_vector_transform.inverse().rotation() << std::endl <<
        //     KBLU << global_to_vector_transform.inverse().translation() << std::endl <<
        //     KNRM << "affineMatrixTest2 = "  << std::endl << 
        //     KBLU << test_transform.rotation() << KNRM << std::endl <<
        //     KBLU << test_transform.translation() << KNRM << std::endl;

        // Initialize variables for RRT (reset)
        nodes.clear();

        start_node.position = param.s_e.first;
        start_node.parent = NULL;
        start_node.children.clear();

        nodes.push_back(&start_node);
        end_node.position = param.s_e.second;
        end_node.parent = NULL;
        end_node.children.clear();

        // Eigen::Vector3d transformed_start = transform_vector_with_affine(
        //     param.s_e.first, global_to_vector_transform);
        // Eigen::Vector3d transformed_end = transform_vector_with_affine(
        //     param.s_e.second, global_to_vector_transform);

        /** @brief Debug message **/
        // std::cout << "transformed_start = " << KBLU << transformed_start.transpose() << KNRM << 
        //     " transformed_end = " << KBLU << transformed_end.transpose() << KNRM << std::endl;

        // Eigen::Vector3d original_start = transform_vector_with_affine(
        //     transformed_start, global_to_vector_transform.inverse());
        // Eigen::Vector3d original_end = transform_vector_with_affine(
        //     transformed_end, global_to_vector_transform.inverse());

        /** @brief Debug message **/
        // std::cout << "original_start = " << KBLU << original_start.transpose() << KNRM << 
        //     " original_end = " << KBLU << original_end.transpose() << KNRM << std::endl;

        _goal_within_sensory_bounds = 
            translational_difference_distance - param.s_b < 0;
        
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
        std::uniform_real_distribution<double> dis_middle(-param.s_b, param.s_b);
        std::uniform_real_distribution<double> dis_height(param.h_c.first, param.h_c.second);

        Node* step_node = new Node;

        
        Eigen::Vector3d transformed_random_vector, random_vector;

        if (_occupied_points > 0) // When there are points in the cloud
        {
            while (1)
            {
                pcl::PointXYZ point;

                // param.s_b/2 is needed to add an offset so that the random point doesnt lie so far behind the agent
                transformed_random_vector = Eigen::Vector3d(
                    min(dis_middle(generator) + param.s_b/2, param.s_b), 
                    dis_middle(generator), dis_middle(generator));
                random_vector = transform_vector_with_affine(
                    transformed_random_vector, global_to_vector_transform.inverse());             

                // Clamp z axis
                random_vector.z() = dis_height(generator);

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
                Eigen::Vector3d(dis_middle(generator) + param.s_b/2, dis_middle(generator), dis_middle(generator));
            random_vector = transform_vector_with_affine(
                transformed_random_vector, global_to_vector_transform.inverse());
            // Clamp z axis
            random_vector.z() = dis_height(generator);
        }

        double transformed_random_vector_distance = transformed_random_vector.norm();

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
        std::cout << "Random_node = " << random_vector.transpose() << std::endl;

        bool flag = check_line_validity(
            nodes[index]->position, step_node->position);
        
        if(!flag)
            return;
        
        // Add the new node into the list and push_back data on children and parent
        step_node->parent = nodes[index];
        nodes.push_back(step_node);
        nodes[index]->children.push_back(step_node);

        if ((transformed_random_vector_distance > param.s_b && transformed_random_vector.x() > 0) ||
            check_line_validity(step_node->position, end_node.position))
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

        /** @brief Debug message **/
        // std::cout << "p = " << KGRN << p.transpose() << KNRM << 
        //     " q = " << KGRN << q.transpose() << KNRM << std::endl;
            
        /** @brief Method 1 **/
        Eigen::Vector3d translational_difference = q - p;

        Eigen::Vector3d translational_difference_vec = 
            translational_difference / translational_difference.norm();
        Eigen::Quaterniond quat = quaternion_from_pitch_yaw(
            Vector3d(1,0,0), translational_difference_vec);

        // Setup vector transform
        Eigen::Affine3d vector_transform = Eigen::Affine3d::Identity(); 
        vector_transform.rotate(quat.inverse());
        vector_transform.translate(-p);
        
        Eigen::Vector3d t_p = transform_vector_with_affine(p, vector_transform);        
        Eigen::Vector3d t_q = transform_vector_with_affine(q, vector_transform);

        /** @brief Debug message **/
        // std::cout << "t_p = " << KGRN << t_p.transpose() << KNRM << 
        //     " t_q = " << KGRN << t_q.transpose() << KNRM << std::endl;

        vector<std::pair<Eigen::Vector3d,Eigen::Vector3d>> query_pairs;

        int expansion = 0;

        // Create a transformed discretized plane
        for (int i = 0; i < 1 + expansion * 2; i++)
        {
            for (int j = 0; j < 1 + expansion * 2; j++)
            {
                for (int k = 0; k < 1 + expansion * 2; k++)
                {
                    std::pair<Eigen::Vector3d,Eigen::Vector3d> point_pair;
                    point_pair.first = Eigen::Vector3d(
                        t_p.x() - param.r * expansion + i*param.r,
                        t_p.y() - param.r * expansion + j*param.r,
                        t_p.z() - param.r * expansion + k*param.r
                    );
                    point_pair.second = Eigen::Vector3d(
                        t_q.x() - param.r * expansion + i*param.r,
                        t_q.y() - param.r * expansion + j*param.r,
                        t_q.z() - param.r * expansion + k*param.r
                    );

                    /** @brief Debug message **/
                    // std::cout << "param.r = " << param.r <<
                    //     " point_pair.first = " << KBLU << point_pair.first.transpose() << KNRM << 
                    //     " point_pair.second = " << KBLU << point_pair.second.transpose() << KNRM << std::endl;

                    query_pairs.push_back(point_pair);
                }
            }
        }


        // Expand the voxel search to encompass a larger box-like area
        for (int i = 0; i < (int)query_pairs.size(); i++)
        {
            float precision = (float)param.r;
            
            /** @brief Debug message **/
            // std::cout << "t_p1 = " << KCYN << query_pairs[i].first.transpose() << KNRM << 
            //     " t_q1 = " << KCYN << query_pairs[i].second.transpose() << KNRM << std::endl;

            Eigen::Vector3d o_p = transform_vector_with_affine(
                query_pairs[i].first, vector_transform.inverse());
            Eigen::Vector3d o_q = transform_vector_with_affine(
                query_pairs[i].second, vector_transform.inverse());
            if (!point_within_octree(o_p) && !point_within_octree(o_q))
                continue;

            /** @brief Debug message **/
            // std::cout << "o_p = " << KRED << o_p.transpose() << KNRM << 
            //     " o_q = " << KRED << o_q.transpose() << KNRM << std::endl;

            pcl::PointCloud<pcl::PointXYZ>::VectorType voxels_in_line_search;
            Eigen::Vector3f p_f = Eigen::Vector3f(
                (float)o_p.x(), (float)o_p.y(), (float)o_p.z());
            // Eigen::Vector3f q_f = Eigen::Vector3f(
            //     (float)o_q.x(), (float)o_q.y(), (float)o_q.z());

            double distance_target = (o_q-o_p).norm();
            Eigen::Vector3d direction = (o_q-o_p).normalized();
            double distance = precision * 2;
            Eigen::Vector3d extension = distance * direction;
            double accumulated_distance = 0;

            Eigen::Vector3f q_f;
            
            // if o_p is inside the octree
            if (point_within_octree(o_p))
            {
                std::cout << KGRN << "o_p inside octree" << KNRM << std::endl;
                Eigen::Vector3d query = o_p;
                while (accumulated_distance - distance_target < 0)
                {
                    query = query + extension;
                    /** @brief Debug message **/
                    // std::cout << "query = " << query.transpose() << std::endl;

                    if (!point_within_octree(query))
                    {
                        query = query - extension;
                        q_f = Eigen::Vector3f(
                            (float)query.x(), (float)query.y(), 
                            (float)query.z());
                        break;
                    }

                    accumulated_distance += distance;
                }

                if (accumulated_distance > distance_target)
                    q_f = Eigen::Vector3f(
                        (float)o_q.x(), (float)o_q.y(), (float)o_q.z());
            }
            // else we have to find an suitable o_p inside the octree to use 
            // octree.getApproxIntersectedVoxelCentersBySegment
            else
            {
                std::cout << KRED << "o_p outside octree" << KNRM << std::endl;
                Eigen::Vector3d query = o_p;
                bool p_f_found, q_f_found;
                // move forward into the octree
                while (accumulated_distance - distance_target < 0)
                {
                    query = query + extension;
                    if (point_within_octree(query) && !p_f_found)
                    {
                        p_f = Eigen::Vector3f(
                            (float)query.x(), (float)query.y(), 
                            (float)query.z());
                        p_f_found = true;
                    }

                    if (!point_within_octree(query) && p_f_found)
                    {
                        query = query - extension;
                        q_f = Eigen::Vector3f(
                            (float)query.x(), (float)query.y(), 
                            (float)query.z());
                        q_f_found = true;
                        break;
                    }

                    accumulated_distance += distance;
                }

                // means q_f and p_f are both outside of the octree
                if (!q_f_found)
                {
                    std::cout << KGRN << "q_f and p_f are both outside" << KNRM << std::endl;
                    continue;
                }
            }
            
            std::cout << "p_f " << p_f.transpose() << " q_f " << q_f.transpose() << std::endl;
            std::cout << "getApproxIntersectedVoxelCentersBySegment" << std::endl;
            int voxels = (int)_octree.getApproxIntersectedVoxelCentersBySegment(
                p_f, q_f, voxels_in_line_search, precision);

            /** @brief Debug message **/
            std::cout << "    voxel_size = " << voxels_in_line_search.size() << std::endl;
            
            for (int j = 0; j < voxels; j++)
            {
                // Check octree boundary
                if (!point_within_octree(Eigen::Vector3d(
                    voxels_in_line_search[j].x, 
                    voxels_in_line_search[j].y, 
                    voxels_in_line_search[j].z)))
                    continue;

                if (_octree.isVoxelOccupiedAtPoint(voxels_in_line_search[j]))
                    return false;
            }
        }
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
        /** @brief End of Method 2 **/

        // std::cout << KGRN << "Out of loop" << KNRM << std::endl;
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