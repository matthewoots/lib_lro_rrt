/*
 * test.cpp
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

#include "tbborrt_server.h"
#include <cmath>
#include <random>
#include <ctime>

#define KNRM  "\033[0m"
#define KRED  "\033[31m"
#define KGRN  "\033[32m"
#define KYEL  "\033[33m"
#define KBLU  "\033[34m"
#define KMAG  "\033[35m"
#define KCYN  "\033[36m"
#define KWHT  "\033[37m"

using namespace tbborrt_server;

int main()
{
    srand ((unsigned int) time (NULL));

    time_point<std::chrono::system_clock> time = system_clock::now();

    double time_duration = duration<double>(system_clock::now() - time).count();

    std::random_device dev;
    std::mt19937 generator(dev());
    
    double map_size = 6.0;
    std::pair<double,double> height_constrain{0.0, 5.0};
    std::pair<double,double> runtime_error{0.025, 0.1};
    double sensor_range = 3.0;
    double protected_zone = 0.5;
    vector<Eigen::Vector4d> no_fly_zone;
    vector<Eigen::Vector3d> previous_input;
    pcl::PointCloud<pcl::PointXYZ>::Ptr obs_pcl (new pcl::PointCloud<pcl::PointXYZ>());

    double resolution = protected_zone;
    tbborrt_server::tbborrt_server_node rrt(resolution);

    bool fail = false;

    while (!fail)
    {
        // Generate pointcloud data and produce lines as obstacles
        int obstacle_count = 10;
        double line_size = 1.0;
        obs_pcl->width = obstacle_count;
        obs_pcl->height = 1;
        obs_pcl->points.resize (obs_pcl->width * obs_pcl->height);

        for (std::size_t i = 0; i < obs_pcl->size (); ++i)
        {
            (*obs_pcl)[i].x = 2 * map_size * rand () / (RAND_MAX + 1.0f) - map_size;
            (*obs_pcl)[i].y = 2 * map_size * rand () / (RAND_MAX + 1.0f) - map_size;
            (*obs_pcl)[i].z = height_constrain.second * rand () / (RAND_MAX + 1.0f);
        }

        int divisions = (int)(line_size / (protected_zone/2.0));
        double divisions_interval = line_size / divisions;
        for (int i = 0; i < obstacle_count; ++i)
        {
            for (int j = 0; j < divisions; j++)
            {
                pcl::PointXYZ point;
                point.x = (*obs_pcl)[i].x - line_size/2 + j*divisions_interval;
                point.y = (*obs_pcl)[i].y - line_size/2 + j*divisions_interval;
                point.z = (*obs_pcl)[i].z - line_size/2 + j*divisions_interval;
                obs_pcl->points.push_back(point);
            }
        }

        std::cout << "obstacle size = " << obs_pcl->points.size() << std::endl;

        std::uniform_real_distribution<double> dis_middle(-map_size, map_size);
        std::uniform_real_distribution<double> dis_height(height_constrain.first, height_constrain.second);

        rrt.set_parameters(protected_zone, no_fly_zone, 
            runtime_error, height_constrain, sensor_range);

        rrt.update_octree(obs_pcl);

        Eigen::Vector3d start = 
            Eigen::Vector3d(dis_middle(generator), dis_middle(generator), dis_height(generator));
        Eigen::Vector3d end = 
            Eigen::Vector3d(dis_middle(generator), dis_middle(generator), dis_height(generator));
        
        std::cout << "start_position = " << KGRN << start.transpose() << KNRM << " " <<
                    "end_position = " << KGRN << end.transpose() << KNRM << " " <<
                    "distance = " << KGRN << (start-end).norm() << KNRM << std::endl;
        
        std::pair<Eigen::Vector3d, Eigen::Vector3d> start_end;
        start_end.first = start;
        start_end.second = end;
        vector<Eigen::Vector3d> search = rrt.find_path(previous_input, start_end);
        if ((int)search.size() > 0)
        {
            fail = false;
            std::cout << KGRN << "Success! Size of " << search.size() << " found" << KNRM << std::endl;
        }
        else
        {
            fail = true;
            std::cout << KRED << "Failure! Size of 0 found" << KNRM << std::endl;
        }

        std::cout << "Duration = " << 
            KGRN << time_duration << KNRM << "s" << std::endl;
    }

    return 0;
}
