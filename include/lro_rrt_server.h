/*
 * lro_rrt_server.h
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
/* 
* Some Documentations
* https://pointclouds.org/documentation/classpcl_1_1octree_1_1_octree_point_cloud.html
* https://pcl.readthedocs.io/projects/tutorials/en/latest/octree.html
* https://github.com/otherlab/pcl/blob/master/test/octree/test_octree.cpp 
*/

#ifndef LRO_RRT_SERVER_H
#define LRO_RRT_SERVER_H

#include <iostream>
#include <cmath>
#include <random>
#include <math.h>
#include <float.h>
#include <queue>
#include <vector>
#include <chrono>
#include <mutex>
#include <Eigen/Core>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

// #include <pcl/octree/octree.h>
#include <pcl/octree/octree_search.h>
#include "kdtree.h"
#include "lro_rrt_sampler.h"
#include "lro_rrt_helper.h"
#include "lro_rrt_struct.h"

using namespace Eigen;
using namespace std;
using namespace std::chrono;

#define KNRM  "\033[0m"
#define KRED  "\033[31m"
#define KGRN  "\033[32m"
#define KYEL  "\033[33m"
#define KBLU  "\033[34m"
#define KMAG  "\033[35m"
#define KCYN  "\033[36m"
#define KWHT  "\033[37m"

#define EXPANSION 1

typedef time_point<std::chrono::system_clock> t_p_sc; // giving a typename

namespace lro_rrt_server
{

    class lro_rrt_server_node
    {
        public:

            Node *start_node, *end_node;

            std::vector<Node*> nodes;
            kdtree *kd_tree;

            std::vector<Eigen::Vector3d> global_path;
            std::vector<Eigen::Vector3d> vertices;
            std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> edges;

            pcl::PointCloud<pcl::PointXYZ>::VectorType occupied_voxels;

            /** 
             * @brief lro_rrt_server_node
             * Constructor of the RRT node
            **/ 
            lro_rrt_server_node()
            {
                offset_list = create_expanded_list(EXPANSION);
            }

            /** 
             * @brief ~lro_rrt_server_node
             * Destructor of the RRT node
            **/ 
            ~lro_rrt_server_node()
            {
                _octree.deleteTree();
            }

            /** 
             * @brief get_path
             * Main run function of the rrt module 
            **/ 
            bool get_path(
                vector<Eigen::Vector3d> &output, bool sample_tree);

            /** 
             * @brief get_line_validity
             * Check whether the line between the pair of points is 
             * obstacle free 
            **/
            bool get_line_validity(
                Eigen::Vector3d p, Eigen::Vector3d q);

            /** 
             * @brief initialized and deinitialized
             * To check the initialization of the rrt parameters 
            **/ 
            bool initialized() {return init;}
            void deinitialized() {init = false;}

            /** 
             * @brief set_parameters
             * Setup the parameters for the rrt 
            **/ 
            void set_parameters(parameters parameter);

            /** 
             * @brief get_receding_path
            **/ 
            void get_receding_path(
                Eigen::Vector3d point, double distance,
                std::vector<Eigen::Vector3d> &output)
            {
                output.clear();

                double min_dist = FLT_MAX;
                int index = 0;

                if (global_path.empty())
                    return;
                
                for (size_t i = 1; i < global_path.size(); i++)
                {
                    double line_dist;
                    Eigen::Vector3d v;
                    get_nearest_distance_to_line(
                        point, global_path[i-1], global_path[i], 
                        line_dist, v);

                    if (line_dist < min_dist)
                    {
                        index = i;
                        min_dist = line_dist;
                    }
                }

                printf("index %d\n", index);

                std::vector<Eigen::Vector3d> tmp;

                double travel_dist = 0.0, previous;
                tmp.push_back(point);
                for (int i = index; i < (int)global_path.size(); i++)
                {
                    previous = travel_dist;
                    travel_dist += 
                        (tmp.back() - global_path[i]).norm();

                    printf("%d travel_dist %lf/%lf\n", i, 
                        travel_dist, distance);
                    if (travel_dist <= distance)
                        tmp.push_back(global_path[i]);
                    else
                    {
                        Eigen::Vector3d dir_vector = 
                            (global_path[i] - tmp.back()).normalized();
                        double leftover = 
                            abs(previous - distance);
                        
                        tmp.push_back(
                            tmp.back() + dir_vector * leftover);
                        
                        break;
                    }
                }

                get_discretized_path(tmp, output);
            }

            /** 
             * @brief set_no_fly_zone
             * Setup no fly zone limits for the search 
            **/ 
            void set_no_fly_zone(vector<Eigen::Vector4d> no_fly_zone);

            /** 
             * @brief update_pose_goal
             * Update the pose and goal
            **/
            void update_pose_goal(
                Eigen::Vector3d p, Eigen::Vector3d q);

            /** 
             * @brief update_octree
             * Update the octree, the octree is centered 
             * around the pose due to perception range 
            **/ 
            void update_octree(
                pcl::PointCloud<pcl::PointXYZ>::Ptr obs_pcl);
            
            /** 
             * @brief check_approx_intersection_by_segment
             * Checks for intersected voxels that contain pointclouds 
             * Edit function from:
             * https://pointclouds.org/documentation/octree__pointcloud_8hpp_source.html#l00269
             * Definition at line 269 of file octree_pointcloud.hpp
            **/ 
            bool check_approx_intersection_by_segment(
                const Eigen::Vector3d origin, const Eigen::Vector3d end, 
                Eigen::Vector3d& intersect);

            /** 
             * @brief get_path_validity
             * Check the pairs of points in the path and return whether 
             * the path is still valid
            **/ 
            bool get_path_validity(
                std::vector<Eigen::Vector3d> path);

            /** 
             * @brief gen_leaf_node_center_from_octree_key
             * Edited from the protected function for octree
             * void pcl::octree::OctreePointCloud<PointT, LeafContainerT, BranchContainerT, OctreeT>::
             * genLeafNodeCenterFromOctreeKey(const OctreeKey& key, PointT& point) const
            **/
            void gen_leaf_node_center_from_octree_key(
                const pcl::octree::OctreeKey key, pcl::PointXYZ& point);

            /** 
             * @brief change_node_parent
             * Swap and update the parent node
            **/
            void change_node_parent(
                Node* &node, Node* &parent, 
                const double &cost_from_parent);

            void sample_whole_tree(Node* &root);

            Node* get_safe_point_in_tree(
                Node* root, double distance);

            /** 
             * @brief check_trajectory_collision
            **/ 
            bool check_trajectory_collision(
                std::vector<Eigen::Vector3d> traj, int &index);

        private:

            std::mutex octree_mutex; 

            /** @param param fixed/static parameters for the search **/
            lro_rrt_server::parameters param;

            /** @param search_param dynamic parameters that vary by instances **/
            lro_rrt_server::search_parameters search_param;

            bool reached = false;
            bool init = false;

            vector<Eigen::Vector3i> offset_list;
            
            /** @param sampler initialize the sampling node **/
            lro_rrt_server::sampler_node sampler;
            
            /** @param _octree pcl converted octree class **/
            pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> _octree = 
                decltype(_octree)(0.1);
            
            /** 
             * @brief create_expanded_list
            **/ 
            std::vector<Eigen::Vector3i> 
                create_expanded_list(int extra);

            /** 
             * @brief query_single_node
             * Query a random node found by the sampler and 
             * check for its validity 
            **/
            void query_single_node();

            /** 
             * @brief extract_final_path
             * Find the shortest obstacle free path through the nodes 
             * when the goal node is reached
            **/
            std::vector<Eigen::Vector3d> 
                extract_final_path(Node *n);

            /** 
             * @brief is_point_within_octree
             * Check if the point is within the octree 
            **/
            bool is_point_within_octree(
                Eigen::Vector3d point);
            
            /** 
             * @brief is_index_within_octree
             * Check if the index is within the octree 
            **/
            bool is_index_within_octree(
                pcl::octree::OctreeKey idx);

            /** 
             * @brief gen_octree_key_for_point
             * Edited from the protected function for octree
             * void pcl::octree::OctreePointCloud<PointT, LeafContainerT, BranchContainerT, OctreeT>::
             * genOctreeKeyforPoint(const PointT& point_arg, OctreeKey& key_arg) const
            **/
            void gen_octree_key_for_point(
                const pcl::PointXYZ point_arg, pcl::octree::OctreeKey& key_arg);
    
            Node* get_nearest_neighbour(Eigen::Vector3d p)
            {
                struct kdres *nn = kd_nearest3(
                kd_tree, p.x(), p.y(), p.z());
                if (nn == nullptr)
                {
                    std::cout << KRED << "nearest query error" << KNRM << std::endl;
                    return nullptr;
                }

                Node* nearest_node = (Node*)kd_res_item_data(nn);
                kd_res_free(nn);

                return nearest_node;
            }
    };
}

#endif