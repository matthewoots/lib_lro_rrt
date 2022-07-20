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
#include <pcl/filters/crop_box.h>
#include <cmath>
#include <random>
#include <thread>

#define KNRM  "\033[0m"
#define KRED  "\033[31m"
#define KGRN  "\033[32m"
#define KYEL  "\033[33m"
#define KBLU  "\033[34m"
#define KMAG  "\033[35m"
#define KCYN  "\033[36m"
#define KWHT  "\033[37m"

using namespace tbborrt_server;
using namespace std::this_thread; // sleep_for, sleep_until
using namespace std::chrono; // nanoseconds, system_clock, seconds

int main()
{
    std::random_device dev;
    std::mt19937 generator(dev());
    
    double map_size = 20.0;
    std::pair<double,double> height_constrain{0.0, 5.0};
    std::pair<double,double> runtime_error{0.020, 0.10};
    double sensor_range = 4.0;
    double protected_zone = 0.3;
    double pointcloud_resolution = 0.2;
    vector<Eigen::Vector4d> no_fly_zone;
    vector<Eigen::Vector3d> previous_input;
    pcl::PointCloud<pcl::PointXYZ>::Ptr obs_pcl (new pcl::PointCloud<pcl::PointXYZ>());

    std::uniform_real_distribution<double> dis_middle(-map_size, map_size);
    std::uniform_real_distribution<double> dis_height(height_constrain.first, height_constrain.second);

    double resolution = protected_zone;
    tbborrt_server::tbborrt_server_node rrt(resolution);

    // Generate pointcloud data and produce AA cubes as obstacles
    int obstacle_count = 10;
    double cube_size = 5.0;
    int divisions = (int)(cube_size / pointcloud_resolution);
    int point_counts = pow(divisions, 3);
    double divisions_interval = cube_size / (double)divisions;
    std::cout << "divisions = " << divisions << " " << 
        "point_counts = " << point_counts << " " << 
        "divisions_interval = " << divisions_interval << std::endl;

    obs_pcl->width = obstacle_count * point_counts;
    obs_pcl->height = 1;
    obs_pcl->points.resize (obs_pcl->width * obs_pcl->height);

    // obstacle vector
    vector<Eigen::Vector3d> cube_pos;

    for (int i = 0; i < obstacle_count; i++)
    {
        cube_pos.push_back(
            Eigen::Vector3d(dis_middle(generator), 
            dis_middle(generator), 
            dis_height(generator))
        );
    }
    
    for (int i = 0; i < obstacle_count; i++)
    {
        for (int j = 0; j < divisions; j++)
        {
            for (int k = 0; k < divisions; k++)
            {
                for (int l = 0; l < divisions; l++)
                {
                    pcl::PointXYZ point;
                    point.x = cube_pos[i].x() - cube_size/2 + j*divisions_interval;
                    point.y = cube_pos[i].y() - cube_size/2 + k*divisions_interval;
                    point.z = cube_pos[i].z() - cube_size/2 + l*divisions_interval;
                    obs_pcl->points.push_back(point);
                }
            }
        }
    }

    bool fail = false;
    int iteration = 0;

    while (!fail)
    {
        std::cout << "Iteration " << KBLU << iteration << KNRM << std::endl;

        rrt.set_parameters(protected_zone, no_fly_zone, 
            runtime_error, height_constrain, sensor_range);

        bool start_end_not_valid = true;
        Eigen::Vector3d start, end;

        // To prevent random points being close to any obstacles / inside any obstacles
        while (start_end_not_valid)
        {
            start = 
                Eigen::Vector3d(dis_middle(generator), dis_middle(generator), dis_height(generator));
            end = 
                Eigen::Vector3d(dis_middle(generator), dis_middle(generator), dis_height(generator));
            
            int count = 0;
            for (std::size_t i = 0; i < obs_pcl->points.size(); ++i)
            {
                count = (int)i;
                Eigen::Vector3d start_diff = Eigen::Vector3d(
                    obs_pcl->points[(int)i].x - start.x(), 
                    obs_pcl->points[(int)i].y - start.y(), 
                    obs_pcl->points[(int)i].z - start.z());
                Eigen::Vector3d end_diff = Eigen::Vector3d(
                    obs_pcl->points[(int)i].x - end.x(), 
                    obs_pcl->points[(int)i].y - end.y(), 
                    obs_pcl->points[(int)i].z - end.z());
                if (start_diff.norm() < protected_zone * 1.5 || end_diff.norm() < protected_zone * 1.5 )
                    break;
            }

            if (count == (int)obs_pcl->points.size() - 1)
                break;
        }
        
        std::cout << "start_position = " << KBLU << start.transpose() << KNRM << " " <<
                    "end_position = " << KBLU << end.transpose() << KNRM << " " <<
                    "distance = " << KBLU << (start-end).norm() << KNRM << std::endl;
        
        // Crop to simulate local sensing
        pcl::PointCloud<pcl::PointXYZ>::Ptr output(
            new pcl::PointCloud<pcl::PointXYZ>);
        
        Eigen::Vector3d dimension = Eigen::Vector3d(sensor_range, sensor_range, sensor_range);

        Eigen::Vector3d min = start - dimension;
        Eigen::Vector3d max = start + dimension;

        pcl::CropBox<pcl::PointXYZ> box_filter;
        box_filter.setMin(Eigen::Vector4f(min.x(), min.y(), min.z(), 1.0));
        box_filter.setMax(Eigen::Vector4f(max.x(), max.y(), max.z(), 1.0));

        box_filter.setInputCloud(obs_pcl);
        box_filter.filter(*output);

        std::cout << "local/original size = " << KBLU << 
            output->points.size() << KNRM << " / " << 
            KBLU << obs_pcl->points.size() << KNRM << std::endl;

        rrt.update_octree(output);


        time_point<std::chrono::system_clock> time = system_clock::now();
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

        double time_duration = duration<double>(system_clock::now() - time).count();

        std::cout << "Duration = " << 
            KGRN << time_duration*1000 << KNRM << "ms" << std::endl;

        /** @brief Debug message **/
        for (int i = 0; i < search.size(); i++)
            std::cout << "    node " << i << " = " << search[i].transpose() << std::endl;
        
        iteration++;
        std::cout << std::endl;

        sleep_for(seconds(1));
    }

    return 0;
}
