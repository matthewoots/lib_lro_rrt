/*
 * lro_rrt_sampler.h
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

#ifndef LRO_RRT_SAMPLER_H
#define LRO_RRT_SAMPLER_H

#include <cmath>
#include <random>
#include <math.h>
#include <queue>
#include <vector>
#include <Eigen/Core>

using namespace Eigen;
using namespace std;

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
}

#endif