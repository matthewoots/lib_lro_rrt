/*
 * lro_rrt_helper.cpp
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
#include "lro_rrt_helper.h"

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
    /**
     * @brief get_intersection
     * Get the intersection point between the line and the plane
     * http://www.3dkingdoms.com/weekly/weekly.php?a=3
    **/
    bool inline get_intersection(
        double fDst1, double fDst2, Eigen::Vector3d P1, 
        Eigen::Vector3d P2, Eigen::Vector3d &Hit) 
    {
        if ((fDst1 * fDst2) >= 0.0) return false;
        if (abs(fDst1 - fDst2) < 0.0001) return false; 
        Hit = P1 + (P2 - P1) * (-fDst1 / (fDst2 - fDst1));
        return true;
    }

    /**
     * @brief in_box
     * Whether hit point is in bounding box
     * http://www.3dkingdoms.com/weekly/weekly.php?a=3
    **/
    bool inline in_box(
        Eigen::Vector3d Hit, Eigen::Vector3d B1, 
        Eigen::Vector3d B2, const int Axis) 
    {
        if ( Axis==1 && Hit.z() > B1.z() && Hit.z() < B2.z() && Hit.y() > B1.y() && Hit.y() < B2.y()) return true;
        if ( Axis==2 && Hit.z() > B1.z() && Hit.z() < B2.z() && Hit.x() > B1.x() && Hit.x() < B2.x()) return true;
        if ( Axis==3 && Hit.x() > B1.x() && Hit.x() < B2.x() && Hit.y() > B1.y() && Hit.y() < B2.y()) return true;
        return false;
    }

    /**
     * @brief check_line_box
     * The box in this article is Axis-Aligned and so can be defined by only two 3D points:
     * B1 - the smallest values of X, Y, Z
     * B2 - the largest values of X, Y, Z 
     * returns true if line (L1, L2) intersects with the box (B1, B2)
     * returns intersection point in Hit
     * http://www.3dkingdoms.com/weekly/weekly.php?a=3
    **/
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

    /** 
     * @brief get_discretized_path
     * Discretize the path according their individual legs
    **/
    void get_discretized_path(
        std::vector<Eigen::Vector3d> input, 
        std::vector<Eigen::Vector3d> &output)
    {
        output.clear();
        std::vector<Eigen::Vector3d> dir_vector;
        std::vector<double> norm_vector;
        std::vector<int> seg_vector;
        
        for (int i = 1; i < (int)input.size(); i++)
        {
            Eigen::Vector3d difference = input[i] - input[i-1];
            
            norm_vector.push_back(difference.norm());
            seg_vector.push_back((int)round(difference.norm() / (1.25)));
            dir_vector.push_back(difference.normalized());
        }

        for (int i = 0; i < (int)seg_vector.size(); i++)
        {
            double n = norm_vector[i] / (double)seg_vector[i];
            for (int j = 0; j < seg_vector[i]; j++)
            {
                output.push_back(input[i] + dir_vector[i] * n * j);
            }
        }
        
        output.push_back(input.back());
    }

    /** 
     * @brief get_nearest_distance_to_line
     * Nearest distance from a point to a line, defined by the start
     * and end point
    **/
    void get_nearest_distance_to_line(
        Eigen::Vector3d p, Eigen::Vector3d s, 
        Eigen::Vector3d e, double &d, Eigen::Vector3d &q)
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

    /** 
     * @brief get_nearest_node 
     * Responsible for finding the nearest node in the tree
     * for a particular random node 
    **/ 
    inline int get_nearest_node(
        Node random, Node base_node, std::vector<Node*> nodes)
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

    /** 
     * @brief constrain_between_180
     * Constrain angle to between -pi to pi 
    **/
    double constrain_between_180(double x)
    {
        x = fmod(x + M_PI,2*M_PI);
        if (x < 0)
            x += 2*M_PI;
        return x - M_PI;
    }

    /**
     * @brief get_bresenham_3d_from_origin
     * Modified from:
     * https://gist.github.com/yamamushi/5823518
     */
    void get_bresenham_3d_from_origin( 
        const int x2, const int y2, const int z2, 
        std::vector<Eigen::Vector3i> &idx)
    {
        idx.clear();

        int i, dx, dy, dz, l, m, n, x_inc, y_inc, z_inc, err_1, err_2, dx2, dy2, dz2;
        Eigen::Vector3i point(0, 0, 0);

        dx = x2 - 0;
        dy = y2 - 0;
        dz = z2 - 0;
        x_inc = (dx < 0) ? -1 : 1;
        l = abs(dx);
        y_inc = (dy < 0) ? -1 : 1;
        m = abs(dy);
        z_inc = (dz < 0) ? -1 : 1;
        n = abs(dz);
        dx2 = l << 1;
        dy2 = m << 1;
        dz2 = n << 1;
        
        if ((l >= m) && (l >= n)) {
            err_1 = dy2 - l;
            err_2 = dz2 - l;
            for (i = 0; i < l; i++) 
            {
                idx.push_back(point);
                if (err_1 > 0) {
                    point[1] += y_inc;
                    err_1 -= dx2;
                }
                if (err_2 > 0) {
                    point[2] += z_inc;
                    err_2 -= dx2;
                }
                err_1 += dy2;
                err_2 += dz2;
                point[0] += x_inc;
            }
        } else if ((m >= l) && (m >= n)) {
            err_1 = dx2 - m;
            err_2 = dz2 - m;
            for (i = 0; i < m; i++) {
                idx.push_back(point);
                if (err_1 > 0) {
                    point[0] += x_inc;
                    err_1 -= dy2;
                }
                if (err_2 > 0) {
                    point[2] += z_inc;
                    err_2 -= dy2;
                }
                err_1 += dx2;
                err_2 += dz2;
                point[1] += y_inc;
            }
        } else {
            err_1 = dy2 - n;
            err_2 = dx2 - n;
            for (i = 0; i < n; i++) {
                idx.push_back(point);
                if (err_1 > 0) {
                    point[1] += y_inc;
                    err_1 -= dz2;
                }
                if (err_2 > 0) {
                    point[0] += x_inc;
                    err_2 -= dz2;
                }
                err_1 += dy2;
                err_2 += dx2;
                point[2] += z_inc;
            }
        }
        idx.push_back(point);
    }

}