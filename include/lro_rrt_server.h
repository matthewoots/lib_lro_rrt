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
#include <vector>
#include <chrono>
#include <mutex>
#include <Eigen/Core>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

// #include <pcl/octree/octree.h>
#include <pcl/octree/octree_search.h>

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

            /** 
            * @brief Parameters for the RRT module performance
            * @param height_constrain (h_c) The minimum and maximum height of the search
            * @param protected_zone = Any protected zone to avoid considering waypoints inside
            * @param runtime_error,first() = _sub_runtime_error = The timeout for each search in the module
            * @param runtime_error,second() = _runtime_error = The overall timeout before we close the program 
            * @param octree = pcl converted octree class
            **/
            struct parameters 
            {
                // @runtime_error : consist of _sub_runtime_error and _runtime_error
                std::pair<double,double> r_e; 
                // @height_constrain : consist of _min_height and _max_height
                std::pair<double,double> h_c; 
                double s_i; // @search_interval
                double s_r; // @sensor_range
                double p_z; // @protected_zone
                std::pair<double,double> s_l_h; // @search_limit_hfov
                std::pair<double,double> s_l_v; // @search_limit_vfov
                double s_d_n; // @scaled_min_dist_from_node
                double r; // @resolution
                double m_s; // @map_size
                double s_bf; // @sensor_buffer_factor
                double s_b; // @sensor_buffer
            };

            struct search_parameters 
            {
                std::pair<Eigen::Vector3d, Eigen::Vector3d> s_e; // start and end pair
                vector<Eigen::Vector4d> n_f_z; // no fly zones (minx maxx miny maxy)
                int o_p; // occupied points in the octree
                Eigen::Vector3d mn_b; // minimum boundary for the octree
                Eigen::Vector3d mx_b; // maximum boundary for the octree
                std::uint8_t m_k[3]; // max key for the octree
            };

            /** @brief Constructor of the RRT node**/ 
            lro_rrt_server_node(){}

            /** @brief Destructor of the RRT node**/ 
            ~lro_rrt_server_node(){ _octree.deleteTree();}

            /** @brief Main run module for the RRT server **/ 
            vector<Eigen::Vector3d> find_path(
                vector<Eigen::Vector3d> previous_input);

            /** @brief Check whether the line between the pair of points is obstacle free **/
            bool check_line_validity(
                Eigen::Vector3d p, Eigen::Vector3d q);

            /** @brief To check the initialization of the RRT parameters **/ 
            bool initialized() {return init;}
            void deinitialized() {init = false;}

            /** @brief Setup the parameters for the search **/ 
            void set_parameters(parameters parameter);

            /** @brief Setup no fly zone limits for the search **/ 
            void set_no_fly_zone(vector<Eigen::Vector4d> no_fly_zone);

            /** @brief Update the pose and the octree, since the octree is centered around the pose due to perception range **/ 
            void update_pose_and_octree(
                pcl::PointCloud<pcl::PointXYZ>::Ptr obs_pcl, Eigen::Vector3d p, Eigen::Vector3d q);
            
            /** @brief Used inside check_line_validity and it checks for intersected voxels that contain pointclouds **/ 
            bool check_approx_intersection_by_segment(
                const Eigen::Vector3d origin, const Eigen::Vector3d end, float precision, 
                Eigen::Vector3d& intersect, string mode);

            /** @brief Edited from the protected function for octree
             * void pcl::octree::OctreePointCloud<PointT, LeafContainerT, BranchContainerT, OctreeT>::
             * genOctreeKeyforPoint(const PointT& point_arg, OctreeKey& key_arg) const
            **/
            void gen_octree_key_for_point(
                const pcl::PointXYZ point_arg, pcl::octree::OctreeKey& key_arg)
            {
                // calculate integer key for point coordinates
                key_arg.x = static_cast<uint8_t>((point_arg.x - search_param.mn_b.x()) / param.r);
                key_arg.y = static_cast<uint8_t>((point_arg.y - search_param.mn_b.y()) / param.r);
                key_arg.z = static_cast<uint8_t>((point_arg.z - search_param.mn_b.z()) / param.r);
                
                assert(key_arg.x <= search_param.m_k[0]);
                assert(key_arg.y <= search_param.m_k[1]);
                assert(key_arg.z <= search_param.m_k[2]);
            }

            /** @brief Edited from the protected function for octree
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

            void nearest_distance_to_line(
                Eigen::Vector3d p, Eigen::Vector3d s, Eigen::Vector3d e, double &d, Eigen::Vector3d &q)
            {
                Eigen::Vector3d s_e = e - s;
                Eigen::Vector3d s_p = p - s;
                double s_e_d = s_e.norm();
                double t = (s_e.x() * s_p.x() + s_e.y() * 
                    s_p.y() + s_e.z() * s_p.z()) / s_e_d;
                if(t < 0.0)
                    t = 0;
                if(t > 1.0)
                    t = 1;

                d = (t * s_e).norm();
                q = s + t * s_e;
            }

        private:

            struct Node 
            {
                vector<Node *> children;
                Node *parent;
                Eigen::Vector3d position;
            };

            std::mutex octree_mutex; 

            lro_rrt_server::lro_rrt_server_node::parameters param;

            lro_rrt_server::lro_rrt_server_node::search_parameters search_param;

            Node start_node, end_node;

            vector<Node*> nodes;

            bool reached = false;
            bool set_resolution = false;
            bool init = false;

            int iteration;

            double bearing;

            std::random_device dev;
            
            /** @param _octree pcl converted octree class **/
            pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> _octree  = decltype(_octree)(0.1);
            
            pcl::PointCloud<pcl::PointXYZ>::Ptr p_c;

            /** @brief get_nearest_node is responsible for finding the nearest node in the tree
            * for a particular random node **/ 
            inline int get_nearest_node(Node random, Node base_node);

            /** @brief Adds in a random node and check for the validity **/
            void search_single_node();

            // bool check_valid_hfov(double current_bearing, double transformed_bearing, double fov);

            /** @brief Reorder then shorten the path by finding shortest obstacle free path through the nodes **/
            std::vector<Eigen::Vector3d> path_extraction();

            /** @brief Reorder the path since the output of RRT is inverted **/
            std::vector<Eigen::Vector3d> get_reorder_path(std::vector<Eigen::Vector3d> path);

            /** @brief Shorten the RRT path by trimming the nodes **/
            std::vector<Eigen::Vector3d> get_shorten_path(std::vector<Eigen::Vector3d> path);

            /** @brief Check if the point is within the octree **/
            inline bool point_within_octree(Eigen::Vector3d point);

            /** @brief Constrain angle to between -pi to pi **/
            double constrain_between_180(double x);

            inline Eigen::Quaterniond quaternion_from_pitch_yaw(
                Eigen::Vector3d v1, Eigen::Vector3d v2)
            {
                Eigen::Quaterniond q;

                // https://github.com/toji/gl-matrix/blob/f0583ef53e94bc7e78b78c8a24f09ed5e2f7a20c/src/gl-matrix/quat.js#L54
                // Eigen::Vector3d xUnitVec3 = Eigen::Vector3d(1,0,0);
                // Eigen::Vector3d yUnitVec3 = Eigen::Vector3d(0,1,0);
                // double dot = v1.x()*v2.x() + 
                //     v1.y()*v2.y() + v1.z()*v2.z();

                // Eigen::Quaterniond q;

                // if (dot < -0.999999)
                // {
                //     Eigen::Vector3d tmpvec3 = xUnitVec3.cross(v1);
                //     if (tmpvec3.norm() < 0.000001)
                //         tmpvec3 = yUnitVec3.cross(v1);
                //     Eigen::Vector3d axis = tmpvec3.normalized();
                //     Eigen::Matrix3d m;
                //     m = AngleAxisd(axis.x()*M_PI, Vector3d::UnitX())
                //         * AngleAxisd(axis.y()*M_PI, Vector3d::UnitY())
                //         * AngleAxisd(axis.z()*M_PI, Vector3d::UnitZ());

                //     Eigen::Quaterniond n_q(m);
                //     q = n_q;
                //     q.normalized();
                // }
                // else if (dot > 0.999999) {
                //     q.x() = 0;
                //     q.y() = 0;
                //     q.z() = 0;
                //     q.w() = 1;
                // }
                // else 
                // {
                //     Eigen::Vector3d tmpvec3 = v1.cross(v2);
                //     q.x() = tmpvec3.x();
                //     q.y() = tmpvec3.y();
                //     q.z() = tmpvec3.z();
                //     q.w() = 1 + dot;
                //     q.normalized();
                // }

                // https://stackoverflow.com/questions/1171849/finding-quaternion-representing-the-rotation-from-one-vector-to-another
                // dot_product check
                // Eigen::Quaterniond q = Eigen::Quaterniond::Identity();
                // if (v1.x()*v2.x() + 
                //     v1.y()*v2.y() + v1.z()*v2.z() > 0.999999)
                //     return q;

                // if (v1.x()*v2.x() + 
                //     v1.y()*v2.y() + v1.z()*v2.z() < -0.999999)
                //     return q;

                // Eigen::Vector3d a = v1.cross(v2);
                // q.vec() = Vector3d(a.x(), a.y(), a.z());
                // q.w() = sqrt(pow(v1.norm(),2) * pow(v2.norm(),2)) + 
                //     v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
                // q.normalize();

                // q = Eigen::Quaterniond::FromTwoVectors(v1.normalized(), v2.normalized());

                return q;
            }

    };
}

#endif