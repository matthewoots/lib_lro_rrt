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
    class sampler_node
    {
        public:

            /** 
            * @brief Parameters for the RRT module performance
            * @param r Radius of the ellipsoid
            * @param c Center of the ellipsoid 
            * @param t Transform and oriented to start-end vector
            * @param s Scale on the xyz wrt to body axis
            **/
            struct ellipsoid
            {
                double r;
                Eigen::Affine3d t;
                Eigen::Vector3d c;
                Eigen::Vector3d s;
            };

            struct circle
            {
                std::uniform_real_distribution<double> hfov;
                std::uniform_real_distribution<double> vfov;
                std::uniform_real_distribution<double> range;
                Eigen::Vector3d c;
                double b;
            };

            sampler_node()
            {
                std::random_device dev;
                generator = std::mt19937_64(dev());
                uniform_random = std::uniform_real_distribution<double>(-1.0, 1.0);
                skewed_uniform_random = std::uniform_real_distribution<double>(-0.0, 1.0);
            }
            
            ellipsoid get_ellipsoid_parameters(
                std::pair<Eigen::Vector3d, Eigen::Vector3d> s_e,
                double p_z)
            {
                ellipsoid e;
                Eigen::Vector3d direction = s_e.second - s_e.first;
                e.c = s_e.second + direction/2.0;
                e.r = direction.norm() + 2 * p_z;

                direction.normalized();
                double yaw = atan2(direction.y(), direction.x());
                Eigen::Vector2d h_xy = Eigen::Vector2d(direction.x(), direction.y());
                double length_h_xy = h_xy.norm();
                double pitch = atan2(direction.z(), length_h_xy);
                e.t = get_affine_transform(
                    e.c, Eigen::Vector3d(0.0, pitch, yaw));

                return e;
            }

            void set_constrained_circle_parameters(
                std::pair<Eigen::Vector3d, Eigen::Vector3d> s_e,
                double minimum, double range,
                std::pair<double,double> hfov, 
                std::pair<double,double> vfov)
            {
                cir = {};
                cir.hfov = std::uniform_real_distribution<double>(hfov.first, hfov.second);
                cir.vfov = std::uniform_real_distribution<double>(vfov.first, vfov.second);
                cir.range = std::uniform_real_distribution<double>(minimum * range, 1.0 * range);
                
                cir.c = s_e.first;
                Eigen::Vector3d direction = s_e.second - s_e.first;
                direction.normalized();
                cir.b = atan2(direction.y(), direction.x());
            }

            Eigen::Vector3d get_rand_point_in_circle()
            {
                // https://mathworld.wolfram.com/SphericalCoordinates.html
                // https://karthikkaranth.me/blog/generating-random-points-in-a-sphere/
                double a = -M_PI + cir.b;
                Eigen::Matrix3d yaw;
                yaw << cos(a), -sin(a), 0,
                    sin(a), cos(a), 0,
                    0, 0, 1;
                double theta = cir.hfov(generator) * 2.0*M_PI;
                double phi = cir.vfov(generator) * M_PI;
                double r = cir.range(generator);
                double sin_theta = sin(theta); 
                double cos_theta = cos(theta);
                double sin_phi = sin(phi); 
                double cos_phi = cos(phi);
                
                Eigen::Vector3d point = Eigen::Vector3d(
                    (r * sin_phi * cos_theta),
                    (r * sin_phi * sin_theta),
                    (r * cos_phi));
                
                Eigen::Vector3d rot_point = yaw * point;
                
                return cir.c + rot_point;
            }
        
        private:
            
            std::mt19937_64 generator;
            std::uniform_real_distribution<double> uniform_random;
            std::uniform_real_distribution<double> skewed_uniform_random;
            std::vector<ellipsoid> elps;
            circle cir;

            /**
             * @brief get_affine_transform
             * @param pos Translational position
             * @param rpy Euler angles
             * @param (Return) Affine3d matrix
            **/
            Eigen::Affine3d get_affine_transform(
                Eigen::Vector3d pos, Eigen::Vector3d rpy)
            {
                Eigen::Affine3d affine;

                Eigen::Vector3d orientated_rpy = 
                    Eigen::Vector3d(0.0, rpy.y(), -rpy.z());
                // Get rotation matrix from RPY
                // https://stackoverflow.com/a/21414609
                Eigen::AngleAxisd rollAngle(orientated_rpy.x(), Eigen::Vector3d::UnitX());
                Eigen::AngleAxisd pitchAngle(orientated_rpy.y(), Eigen::Vector3d::UnitY());
                Eigen::AngleAxisd yawAngle(orientated_rpy.z(), Eigen::Vector3d::UnitZ());
                Eigen::Quaternion<double> q = rollAngle * pitchAngle * yawAngle;
                Eigen::Matrix3d rot = q.matrix();

                Eigen::Vector3d rot_pos = rot * -pos;

                affine.translation() = rot_pos;
                affine.linear() = rot;

                return affine;
            }
    };

    class lro_rrt_server_node
    {
        public:

            struct Node 
            {
                EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	            Node() : parent(NULL), cost_from_start(FLT_MAX), cost_from_parent(0.0){};
                std::vector<Node*> children;
                Node *parent;
                // int index;
                Eigen::Vector3d position;
                double cost_from_start;
	            double cost_from_parent; 
            };

            struct node_status_check 
            {
                Node *node;
                bool is_checked;
	            bool is_valid;
            };

            /** 
            * @brief Global parameters for the RRT setup
            * @param r_e
            * @param h_c
            * @param s_i
            * @param s_r
            * @param p_z
            * @param s_l_h
            * @param s_l_v
            * @param s_d_n
            * @param r
            * @param m_s
            * @param s_bf
            * @param s_b
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

            /** 
            * @brief Local parameters for an RRT instance
            * @param s_e start and end pair
            * @param n_f_z no fly zones (minx maxx miny maxy)
            * @param o_p occupied points in the octree
            * @param mn_b minimum boundary for the octree
            * @param mx_b maximum boundary for the octree
            * @param m_k max key for the octree
            **/
            struct search_parameters 
            {
                std::pair<Eigen::Vector3d, Eigen::Vector3d> s_e; // start and end pair
                vector<Eigen::Vector4d> n_f_z; // no fly zones (minx maxx miny maxy)
                int o_p; // occupied points in the octree
                Eigen::Vector3d mn_b; // minimum boundary for the octree
                Eigen::Vector3d mx_b; // maximum boundary for the octree
                std::uint8_t m_k[3]; // max key for the octree
            };

            Node *start_node, *end_node;

            std::vector<Node*> nodes;
            kdtree *kd_tree;

            pcl::PointCloud<pcl::PointXYZ>::VectorType occupied_voxels;

            /** @brief Constructor of the RRT node**/ 
            lro_rrt_server_node(){}

            /** @brief Destructor of the RRT node**/ 
            ~lro_rrt_server_node(){ _octree.deleteTree();}

            /** @brief Main run module for the RRT server **/ 
            bool get_path(
                vector<Eigen::Vector3d> previous_input,
                vector<Eigen::Vector3d> &output);

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
                Eigen::Vector3d& intersect);

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

            void get_nearest_distance_to_line(
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

            // http://www.3dkingdoms.com/weekly/weekly.php?a=3
            bool inline get_intersection(
                double fDst1, double fDst2, Eigen::Vector3d P1, 
                Eigen::Vector3d P2, Eigen::Vector3d &Hit) 
            {
                if ((fDst1 * fDst2) >= 0.0) return false;
                if (abs(fDst1 - fDst2) < 0.0001) return false; 
                Hit = P1 + (P2 - P1) * (-fDst1 / (fDst2 - fDst1));
                return true;
            }

            bool inline in_box(
                Eigen::Vector3d Hit, Eigen::Vector3d B1, 
                Eigen::Vector3d B2, const int Axis) 
            {
                if ( Axis==1 && Hit.z() > B1.z() && Hit.z() < B2.z() && Hit.y() > B1.y() && Hit.y() < B2.y()) return true;
                if ( Axis==2 && Hit.z() > B1.z() && Hit.z() < B2.z() && Hit.x() > B1.x() && Hit.x() < B2.x()) return true;
                if ( Axis==3 && Hit.x() > B1.x() && Hit.x() < B2.x() && Hit.y() > B1.y() && Hit.y() < B2.y()) return true;
                return false;
            }

            // The box in this article is Axis-Aligned and so can be defined by only two 3D points:
            // B1 - the smallest values of X, Y, Z
            // B2 - the largest values of X, Y, Z 
            // returns true if line (L1, L2) intersects with the box (B1, B2)
            // returns intersection point in Hit
            bool check_line_box(
                Eigen::Vector3d B1, Eigen::Vector3d B2, 
                Eigen::Vector3d L1, Eigen::Vector3d L2, 
                Eigen::Vector3d &Hit)
            {
                if (L2.x() < B1.x() && L1.x() < B1.x()) return false;
                if (L2.x() > B2.x() && L1.x() > B2.x()) return false;
                if (L2.y() < B1.y() && L1.y() < B1.y()) return false;
                if (L2.y() > B2.y() && L1.y() > B2.y()) return false;
                if (L2.z() < B1.z() && L1.z() < B1.z()) return false;
                if (L2.z() > B2.z() && L1.z() > B2.z()) return false;
                if (L1.x() > B1.x() && L1.x() < B2.x() &&
                    L1.y() > B1.y() && L1.y() < B2.y() &&
                    L1.z() > B1.z() && L1.z() < B2.z()) 
                    {Hit = L1; 
                    return true;}
                if ( (get_intersection( L1.x()-B1.x(), L2.x()-B1.x(), L1, L2, Hit) && in_box( Hit, B1, B2, 1 ))
                || (get_intersection( L1.y()-B1.y(), L2.y()-B1.y(), L1, L2, Hit) && in_box( Hit, B1, B2, 2 )) 
                || (get_intersection( L1.z()-B1.z(), L2.z()-B1.z(), L1, L2, Hit) && in_box( Hit, B1, B2, 3 )) 
                || (get_intersection( L1.x()-B2.x(), L2.x()-B2.x(), L1, L2, Hit) && in_box( Hit, B1, B2, 1 )) 
                || (get_intersection( L1.y()-B2.y(), L2.y()-B2.y(), L1, L2, Hit) && in_box( Hit, B1, B2, 2 )) 
                || (get_intersection( L1.z()-B2.z(), L2.z()-B2.z(), L1, L2, Hit) && in_box( Hit, B1, B2, 3 )))
                    return true;

                return false;
            }

            void change_node_parent(
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

            void extract_point_cloud_within_boundary(
                Eigen::Vector3d min, Eigen::Vector3d max,
                pcl::PointCloud<pcl::PointXYZ>::Ptr &pc)
            {
                pc->points.clear();
                if (search_param.o_p == 0)
                    return;

                for (auto &v : occupied_voxels)
                {
                    if (point_within_aabb(
                        Eigen::Vector3d(v.x, v.y, v.z), min, max))
                        pc->points.push_back(v);
                }
            }

        private:

            std::mutex octree_mutex; 

            /** @param param fixed/static parameters for the search **/
            lro_rrt_server::lro_rrt_server_node::parameters param;
            /** @param search_param dynamic parameters that vary by instances **/
            lro_rrt_server::lro_rrt_server_node::search_parameters search_param;

            bool reached = false;
            bool set_resolution = false;
            bool init = false;

            lro_rrt_server::sampler_node sampler;
            
            /** @param _octree pcl converted octree class **/
            pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> _octree  = decltype(_octree)(0.1);
            /** @param p_c store pointcloud **/
            pcl::PointCloud<pcl::PointXYZ>::Ptr p_c;

            /** @brief get_nearest_node 
             * responsible for finding the nearest node in the tree
             * for a particular random node 
            **/ 
            inline int get_nearest_node(
                Node random, Node base_node);

            /** @brief query_single_node
             * Adds in a random node and check for the validity 
            **/
            void query_single_node();

            /** @brief path_extraction
             * Reorder then shorten the path by finding shortest obstacle free path through the nodes 
            **/
            std::vector<Eigen::Vector3d> path_extraction();

            /** @brief Reorder the path since the output of RRT is inverted **/
            std::vector<Eigen::Vector3d> get_reorder_path(std::vector<Eigen::Vector3d> path);

            /** @brief Shorten the RRT path by trimming the nodes **/
            // std::vector<Eigen::Vector3d> get_shorten_path(std::vector<Eigen::Vector3d> path);

            /** @brief Check if the point is within the octree **/
            bool point_within_aabb(
                Eigen::Vector3d point, Eigen::Vector3d min,
                Eigen::Vector3d max);

            /** @brief Constrain angle to between -pi to pi **/
            double constrain_between_180(double x);

    };
}

#endif