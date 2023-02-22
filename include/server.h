/*
 * server.h
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

#ifndef LRO_SERVER_H
#define LRO_SERVER_H

#include <cmath>
#include <math.h>
#include <float.h>
#include <queue>
#include <vector>
#include <chrono>
#include <mutex>
#include <Eigen/Core>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

#include "kdtree.h"
#include "sampler.h"
#include "helper.h"
#include "struct.h"

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

typedef std::chrono::time_point<std::chrono::system_clock> t_p_sc; // giving a typename

namespace lro
{
    class server
    {
        public:

            Node *start_node, *end_node;

            std::vector<Node*> nodes;
            kdtree *kd_tree;

            std::vector<Eigen::Vector3d> vertices;

            pcl::PointCloud<pcl::PointXYZ>::VectorType occupied_voxels;

            /** 
             * @brief server
             * Constructor of the lro server node
            **/ 
            server(double sensor_range, 
                double resolution, 
                double scaled_min_dist,
                double runtime_error,
                double refinement_time,
                std::pair<double, double> slh, 
                std::pair<double, double> slv, 
                std::pair<double, double> height_range) : 
                _sensor_range(sensor_range), 
                _resolution(resolution), 
                _scaled_min_dist(scaled_min_dist),
                _runtime_error(runtime_error),
                _refinement_time(refinement_time),
                _search_limit_hfov_list(slh), 
                _search_limit_vfov_list(slv),
                _height_range(height_range)
            {
                _octree.setResolution(_resolution);
                offset_list = create_expanded_list(EXPANSION);
            }

            /** 
             * @brief ~server
             * Destructor of the lro server node
            **/ 
            ~server()
            {
                _octree.deleteTree();
            }

            /** 
             * @brief get_search_path
             * Main run function of the lro module 
            **/ 
            bool get_search_path(
                Eigen::Vector3d &start, Eigen::Vector3d &end,
                vector<Eigen::Vector3d> &output, bool sample_tree);

            /** 
             * @brief get_tree
             * Get the full spanning tree found by the search
            **/
            std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> get_tree();

            /** 
             * @brief set_no_fly_zone
             * Setup no fly zone limits for the search 
            **/ 
            void set_no_fly_zone(
                std::vector<Eigen::Vector4d> no_fly_zone);

            /** 
             * @brief update_octree
             * Update the octree, the octree is centered 
             * around the pose due to perception range 
            **/ 
            void update_octree(
                pcl::PointCloud<pcl::PointXYZ>::Ptr obs_pcl);

            /** 
             * @brief check_trajectory_collision
            **/ 
            bool check_trajectory_collision(
                std::vector<Eigen::Vector3d> traj, int &index);

        private:

            std::mutex octree_mutex; 

            double _sensor_range;
            double _resolution;
            double _scaled_min_dist;
            double _runtime_error;
            double _refinement_time;
            std::pair<double, double> _search_limit_hfov_list;
            std::pair<double, double> _search_limit_vfov_list;
            std::pair<double, double> _height_range;

            /** @param search_param dynamic parameters that vary by instances **/
            lro::search_parameters search_param;

            std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> edges;

            bool reached = false;
            bool init = false;

            vector<Eigen::Vector3i> offset_list;
            
            /** @param sampler initialize the sampling node **/
            lro::sampler sampler;
            
            /** @param _octree pcl converted octree class **/
            pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> _octree = 
                decltype(_octree)(0.1);
            
            /** 
             * @brief create_expanded_list
            **/ 
            std::vector<Eigen::Vector3i> 
                create_expanded_list(int extra);

            /** 
             * @brief get_line_validity
             * Check whether the line between the pair of points is 
             * obstacle free 
            **/
            bool get_line_validity(
                Eigen::Vector3d p, Eigen::Vector3d q);

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
                extract_final_path(Node *n, Eigen::Vector3d &start);

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