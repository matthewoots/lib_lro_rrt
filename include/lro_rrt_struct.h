/*
 * lro_rrt_struct.h
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

#ifndef LRO_RRT_STRUCT_H
#define LRO_RRT_STRUCT_H

#include <float.h>
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
        double r_t; // @refinement_time
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
}

#endif

