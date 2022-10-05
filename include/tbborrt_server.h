/*
 * tbborrt_server.h
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

#ifndef TBBORRT_SERVER_H
#define TBBORRT_SERVER_H

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

namespace tbborrt_server
{

    class tbborrt_server_node
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
                std::pair<Eigen::Vector3d, Eigen::Vector3d> s_e; // start and end pair
                // @runtime_error : consist of _sub_runtime_error and _runtime_error
                std::pair<double,double> r_e; 
                // @height_constrain : consist of _min_height and _max_height
                std::pair<double,double> h_c; 
                double s_i; // @search_interval
                double s_r; // @sensor_range
                double p_z; // @protected_zone
                double r; // @resolution
                double m_s; // @map_size
                double s_bf; // @sensor_buffer_factor
                double s_b; // @sensor_buffer
            };

            /** @brief Constructor of the rrt_server node**/ 
            tbborrt_server_node(){}

            /** @brief Destructor of the rrt_server node**/ 
            ~tbborrt_server_node()
            {
                _octree.deleteTree();
            }

            bool check_line_validity(Eigen::Vector3d p, Eigen::Vector3d q);

            /** @brief Main run module for the rrt_server node
            * @param previous_input = The previous input found that is reusable in the search
            **/ 
            vector<Eigen::Vector3d> find_path(vector<Eigen::Vector3d> previous_input);

            /** @brief Setup the parameters for the search **/ 
            void set_parameters(
                parameters parameter, vector<Eigen::Vector4d> no_fly_zone)
            {
                param = parameters();
                // Write over the previous data
                param = parameter;

                _no_fly_zone.clear();
                _no_fly_zone = no_fly_zone;
                if (!set_resolution)
                {
                    set_resolution = true;
                    _octree.setResolution(param.r);
                }

                // The buffer for the xyz search area determined by the sensor range
                param.s_b = param.s_bf * param.s_r;
            }

            void update_pose_and_octree(
                pcl::PointCloud<pcl::PointXYZ>::Ptr obs_pcl, 
                Eigen::Vector3d p, Eigen::Vector3d q)
            {
                std::lock_guard<std::mutex> octree_lock(octree_mutex);

                param.s_e.first = p;

                // Get boundary defined by p and q that is including buffer
                double max_x_arg = p.x() + param.s_b;
                double min_x_arg = p.x() - param.s_b;
                
                double max_y_arg = p.y() + param.s_b;
                double min_y_arg = p.y() - param.s_b;

                double max_z_arg = param.h_c.second;
                double min_z_arg = param.h_c.first;

                _octree.deleteTree();
                _octree.defineBoundingBox(
                    min_x_arg, min_y_arg, min_z_arg,
                    max_x_arg, max_y_arg, max_z_arg);
                _octree.setInputCloud(obs_pcl);
                _octree.addPointsFromInputCloud();
                /** @brief Debug message **/
                // std::cout << "Pointcloud size = " << KBLU << 
                //     obs_pcl->points.size() << KNRM << std::endl;
                
                _store_cloud = obs_pcl;
                pcl::PointCloud<pcl::PointXYZ>::VectorType _occupied_voxels;
                _occupied_points = _octree.getOccupiedVoxelCenters(_occupied_voxels);
                
                _octree.getBoundingBox(
                    min_bnd.x(), min_bnd.y(), min_bnd.z(),
                    max_bnd.x(), max_bnd.y(), max_bnd.z());
                /** @brief Debug message **/
                // std::cout << "Minimum Boundary = " << KBLU << min_bnd.transpose() << KNRM << " " << 
                //     "Maximum Boundary = " << KBLU << max_bnd.transpose() << KNRM << std::endl;
            }

        private:

            std::mutex octree_mutex, search_mutex; 

            struct Node 
            {
                vector<Node *> children;
                Node *parent;
                Eigen::Vector3d position;
            };

            tbborrt_server::tbborrt_server_node::parameters param;

            Node start_node;
            Node end_node;

            vector<Node*> nodes;
            bool reached = false;
            int iteration;
            std::random_device dev;
            bool set_resolution = false;
            
            /** @param _octree = pcl converted octree class **/
            pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> _octree  = decltype(_octree)(0.1);
            int _occupied_points;
            pcl::PointCloud<pcl::PointXYZ>::Ptr _store_cloud;
            Eigen::Vector3d min_bnd, max_bnd;

            
            /** @brief Parameters for the local map expansion size in the RRT module **/
            bool _goal_within_sensory_bounds;
            vector<Eigen::Vector4d> _no_fly_zone;

            double bearing;
            
            inline int get_nearest_node(Node random, Node base_node);

            inline double separation(Eigen::Vector3d p, Eigen::Vector3d q) {return (p - q).norm();}

            void search_single_node();

            bool check_valid_hfov(double current_bearing, double transformed_bearing, double fov);

            inline Eigen::Quaterniond quaternion_from_pitch_yaw(
                Eigen::Vector3d v1, Eigen::Vector3d v2)
            {

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

                Eigen::Quaterniond q = 
                    Eigen::Quaterniond::FromTwoVectors(v1.normalized(), v2.normalized());

                return q;
            }

            /** @brief Transform pose according to the translation and q given in affine form*/
            inline Eigen::Vector3d transform_vector_with_affine(
                Eigen::Vector3d p, Eigen::Affine3d t)
            {
                Eigen::Affine3d pose(Affine3d::Identity());
                pose.translation() = p;

                Eigen::Affine3d transformed_pose = t * pose;   

                return transformed_pose.translation();
            }

            std::vector<Eigen::Vector3d> path_extraction()
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

            inline std::vector<Eigen::Vector3d> get_reorder_path(
                std::vector<Eigen::Vector3d> path)
            {
                std::vector<Eigen::Vector3d> reordered_path;
                reordered_path.push_back(start_node.position);
                for (int i = (int)path.size()-1; i >= 0; i--)
                    reordered_path.push_back(path[i]);            

                return reordered_path;
            }

            inline std::vector<Eigen::Vector3d> get_shorten_path(
                std::vector<Eigen::Vector3d> path)
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

            inline bool point_within_octree(Eigen::Vector3d point)
            {
                // Check octree boundary
                if (point.x() < max_bnd.x() - param.r/2 && 
                    point.x() > min_bnd.x() + param.r/2 &&
                    point.y() < max_bnd.y() - param.r/2 && 
                    point.y() > min_bnd.y() + param.r/2 &&
                    point.z() < max_bnd.z() - param.r/2 && 
                    point.z() > min_bnd.z() + param.r/2)
                    return true;
                else
                    return false;
            } 

    };
}

#endif