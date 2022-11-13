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

namespace lro_rrt_server
{

    class lro_rrt_server_node
    {
        public:

            Node *start_node, *end_node;

            std::vector<Node*> nodes;
            kdtree *kd_tree;

            pcl::PointCloud<pcl::PointXYZ>::VectorType occupied_voxels;

            /** 
             * @brief lro_rrt_server_node
             * Constructor of the RRT node
            **/ 
            lro_rrt_server_node(){}

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
            bool get_path(vector<Eigen::Vector3d> &output);

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
             * @brief set_no_fly_zone
             * Setup no fly zone limits for the search 
            **/ 
            void set_no_fly_zone(vector<Eigen::Vector4d> no_fly_zone);

            /** 
             * @brief update_pose_and_octree
             * Update the pose and the octree, since the octree is centered 
             * around the pose due to perception range 
            **/ 
            void update_pose_and_octree(
                pcl::PointCloud<pcl::PointXYZ>::Ptr obs_pcl, Eigen::Vector3d p, Eigen::Vector3d q);
            
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
                std::vector<Eigen::Vector3d> path)
            {
                // If it contains just 1 node which is its current point
                if (path.size() == 1)
                    return false;

                // Check to see whether the new control point and the previous inputs
                // have any pointclouds lying inside
                int last_safe_idx = -1;
                for (int i = 0; i < path.size()-1; i++)
                {
                    if (!get_line_validity(
                        path[i], path[i+1]))
                    {
                        last_safe_idx = i;
                        break;
                    }
                }

                if (last_safe_idx >= 0)
                    return false;
                else
                    return true;
            }

            /** 
             * @brief gen_leaf_node_center_from_octree_key
             * Edited from the protected function for octree
             * void pcl::octree::OctreePointCloud<PointT, LeafContainerT, BranchContainerT, OctreeT>::
             * genLeafNodeCenterFromOctreeKey(const OctreeKey& key, PointT& point) const
            **/
            void gen_leaf_node_center_from_octree_key(
                const pcl::octree::OctreeKey key, pcl::PointXYZ& point)
            {
                // define point to leaf node voxel center
                point.x = static_cast<float>(
                    (static_cast<double>(key.x) + 0.5f) * 
                    param.r + search_param.mn_b.x());
                point.y = static_cast<float>(
                    (static_cast<double>(key.y) + 0.5f) * 
                    param.r + search_param.mn_b.y());
                point.z =static_cast<float>(
                    (static_cast<double>(key.z) + 0.5f) * 
                    param.r + search_param.mn_b.z());
            }

            

            void extract_point_cloud_within_boundary(
                Eigen::Vector3d c, double r, 
                pcl::PointCloud<pcl::PointXYZ>::Ptr &pc)
            {
                pc->points.clear();
                if (search_param.o_p == 0)
                    return;

                double rr = pow(r, 2);
                for (auto &v : occupied_voxels)
                {
                    double xx = pow(v.x - c.x(), 2);
                    double yy = pow(v.y - c.y(), 2);
                    double zz = pow(v.z - c.z(), 2);
                    if (xx + yy + zz < rr)
                        pc->points.push_back(v);
                }
            }

            void get_estimated_center_of_point(
                Eigen::Vector3d p, Eigen::Vector3d &est)
            {
                Eigen::Vector3d dist = (search_param.mn_b - p);
                Eigen::Vector3d v = (search_param.mn_b - p).normalized();
                int nx = (int)round(abs(dist.x())/param.r);
                int ny = (int)round(abs(dist.y())/param.r);
                int nz = (int)round(abs(dist.z())/param.r);
                est = Eigen::Vector3d(
                    (nx + 0.5f)*param.r + search_param.mn_b.x(),
                    (ny + 0.5f)*param.r + search_param.mn_b.y(),
                    (nz + 0.5f)*param.r + search_param.mn_b.z()
                );
            }

            /** 
             * @brief change_node_parent
             * Swap and update the parent node
            **/
            void change_node_parent(
                Node* &node, Node* &parent, 
                const double &cost_from_parent);

        private:

            std::mutex octree_mutex; 

            /** @param param fixed/static parameters for the search **/
            lro_rrt_server::parameters param;

            /** @param search_param dynamic parameters that vary by instances **/
            lro_rrt_server::search_parameters search_param;

            bool reached = false;
            bool init = false;
            
            /** @param sampler initialize the sampling node **/
            lro_rrt_server::sampler_node sampler;
            
            /** @param _octree pcl converted octree class **/
            pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> _octree = 
                decltype(_octree)(0.1);
            
            /** @param p_c store pointcloud **/
            pcl::PointCloud<pcl::PointXYZ>::Ptr p_c;

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
             * @brief gen_octree_key_for_point
             * Edited from the protected function for octree
             * void pcl::octree::OctreePointCloud<PointT, LeafContainerT, BranchContainerT, OctreeT>::
             * genOctreeKeyforPoint(const PointT& point_arg, OctreeKey& key_arg) const
            **/
            void gen_octree_key_for_point(
                const pcl::PointXYZ point_arg, pcl::octree::OctreeKey& key_arg);
    };
}

#endif