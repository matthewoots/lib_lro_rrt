/*
 * lro_rrt_helper.h
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

#ifndef LRO_RRT_HELPER_H
#define LRO_RRT_HELPER_H

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
#include "lro_rrt_struct.h"

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
        Eigen::Vector3d P2, Eigen::Vector3d &Hit);

    /**
     * @brief in_box
     * Whether hit point is in bounding box
     * http://www.3dkingdoms.com/weekly/weekly.php?a=3
    **/
    bool inline in_box(
        Eigen::Vector3d Hit, Eigen::Vector3d B1, 
        Eigen::Vector3d B2, const int Axis);
    
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
        Eigen::Vector3d &Hit);

    /** 
     * @brief get_discretized_path
     * Discretize the path according their individual legs
    **/
    void get_discretized_path(
        std::vector<Eigen::Vector3d> input, 
        std::vector<Eigen::Vector3d> &output);

    /** 
     * @brief get_nearest_distance_to_line
     * Nearest distance from a point to a line, defined by the start
     * and end point
    **/
    void get_nearest_distance_to_line(
        Eigen::Vector3d p, Eigen::Vector3d s, 
        Eigen::Vector3d e, double &d, Eigen::Vector3d &q);
    
    /** 
     * @brief get_nearest_node 
     * Responsible for finding the nearest node in the tree
     * for a particular random node 
    **/ 
    inline int get_nearest_node(
        Node random, Node base_node, std::vector<Node*> nodes);

    /** 
     * @brief constrain_between_180
     * Constrain angle to between -pi to pi 
    **/
    double constrain_between_180(double x);

    /**
     * @brief get_bresenham_3d_from_origin
     * Modified from:
     * https://gist.github.com/yamamushi/5823518
     */
    void get_bresenham_3d_from_origin( 
        const int x2, const int y2, const int z2, 
        std::vector<Eigen::Vector3i> &idx);
}

#endif