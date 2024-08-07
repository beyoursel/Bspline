#include <iostream>
#include <vector>
#include <algorithm>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/common/common.h>
#include <chrono>
#include <boost/filesystem.hpp>
#include <omp.h>
#include <bspline.h>
#include <queue>
#include <limits>
#include <cmath>
#include <iomanip>
#include "bspline_api.h"


#define M_PI 3.14159265358979323846
#define HOME_LAT_ 31.8351078
#define HOME_LON_ 118.7806222

float home_x_;
float home_y_;

using namespace std;


void SetHome(double latitude, double longitude) {

  home_x_ = longitude / (360.0 / (40030173.0 * cos(latitude * M_PI / 180.0)));
  home_y_ = latitude / 0.00000899;
}


double M2Lat(double d) {
  // 1米对应的纬度变化约为1/111111.0度
  return HOME_LAT_ + d * 0.00000899;    
}


double M2Lon(double d)
{
  // 根据当前纬度计算1米对应的经度变化
  double latitudeInRadians = HOME_LAT_ * M_PI / 180.0;
  double metersPerDegreeLongitude = 111320.0 * cos(latitudeInRadians);
  return HOME_LON_ + d / metersPerDegreeLongitude;
}


double Lat2M(double latitude)
{
  return latitude / 0.00000899 - home_y_;
}

double Lon2M(double longitude)
{
  return longitude / (360.0 / (40030173.0 * cos(HOME_LAT_ * M_PI / 180.0))) - home_x_;
}



void VisualizePointCloudV2(pcl::PointCloud<pcl::PointXYZ>::Ptr downsample_cloud,
    const std::vector<Point>&  control_points)
{
    pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
    viewer->setBackgroundColor(0, 0, 0);

    pcl::PointCloud<pcl::PointXYZ>::Ptr controlCloud(new pcl::PointCloud<pcl::PointXYZ>);

    for (const auto& point : control_points)
    {
        controlCloud->points.push_back(pcl::PointXYZ(point.x, point.y, point.z));
    }

    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud_color_handler(downsample_cloud, 255, 255, 255);
    viewer->addPointCloud(downsample_cloud, cloud_color_handler, "downsample_points");

    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> control_color_handler(controlCloud, 255, 0, 0);
    viewer->addPointCloud(controlCloud, control_color_handler, "control_points");

    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "downsample_points");
    // viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4, "fitted_points");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "control_points");

    viewer->addCoordinateSystem(1.0);
    viewer->initCameraParameters();

    while (!viewer->wasStopped())
    {
        viewer->spinOnce(100);
    }
}


void VisualizePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr downsample_cloud,
    const std::vector<Point>& vertices,
    pcl::PointCloud<pcl::PointXYZ>::Ptr  control_points)
{
    pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
    viewer->setBackgroundColor(0, 0, 0);

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    // pcl::PointCloud<pcl::PointXYZ>::Ptr controlCloud(new pcl::PointCloud<pcl::PointXYZ>);

    for (const auto& vertex : vertices)
    {
        cloud->points.push_back(pcl::PointXYZ(vertex.x, vertex.y, vertex.z));
    }

    // for (const auto& point : control_points)
    // {
    //     controlCloud->points.push_back(pcl::PointXYZ(point[0], point[1], point[2]));
    // }


    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud_color_handler(cloud, 255, 255, 255);
    viewer->addPointCloud(downsample_cloud, cloud_color_handler, "downsample_points");

    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> downsample_color_handler(cloud, 0, 255, 0);
    viewer->addPointCloud(cloud, downsample_color_handler, "fitted_points");

    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> control_color_handler(control_points, 255, 0, 0);
    viewer->addPointCloud(control_points, control_color_handler, "control_points");

    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "control_points");

    viewer->addCoordinateSystem(1.0);
    viewer->initCameraParameters();

    while (!viewer->wasStopped())
    {
        viewer->spinOnce(100);
    }
}

float GetAverageHeight(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud) {

    float temp_height = 0.0;
    for (auto& point: cloud->points) {
        temp_height += point.z;
    }
    return temp_height / cloud->points.size();
}


void VoxelDownSample(pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud , pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered, float voxel_size) {
    // Create the filtering object
    pcl::VoxelGrid<pcl::PointXYZ> sor;

    sor.setInputCloud(input_cloud);
    sor.setLeafSize(voxel_size, voxel_size, voxel_size);  // voxel_size     
    // sor.setLeafSize(0.5f, 0.5f, 0.5f);  // voxel_size 
    sor.filter(*cloud_filtered);
}

void StatisticalRemoveOutlier(pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud , pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered, int mean_k, float std_thresh) {
    // Create the filtering object
    pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor;
    sor.setInputCloud (input_cloud);
    sor.setMeanK(mean_k); // default: 10 越大越严格
    sor.setStddevMulThresh (std_thresh); // default: 2.0 越小越严格
    sor.filter(*cloud_filtered);
}

double time_inc(std::chrono::high_resolution_clock::time_point &t_end,
                std::chrono::high_resolution_clock::time_point &t_begin) {
  
  return std::chrono::duration_cast<std::chrono::duration<double>>(t_end -t_begin).count() * 1000;
}


int main(int argc, char** argv) {


    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_pcd_file>" << std::endl;
        return -1;
    }

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    if (pcl::io::loadPCDFile<pcl::PointXYZ>(argv[1], *cloud) == -1) {
        PCL_ERROR("Couldn't read file %s \n", argv[1]);
        return -1;
    }


    // auto start = std::chrono::high_resolution_clock::now();
    pcl::PointCloud<pcl::PointXYZ>::Ptr filtered_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr filtered_ground(new pcl::PointCloud<pcl::PointXYZ>);

    // PassthroghFilter(cloud, cloud);

    VoxelDownSample(cloud, filtered_cloud, 1.0);

    StatisticalRemoveOutlier(filtered_cloud, filtered_ground, 50, 1.0);

    std::cout << "cloud size is: " << filtered_ground->points.size() << std::endl;


    auto start = std::chrono::high_resolution_clock::now();
    BspFitting bsp_fit(filtered_ground, 5.0, 3);



    Point ptc = bsp_fit.GetBsplinePoint(0, 0); // test

    auto end = std::chrono::high_resolution_clock::now();
    double duration = time_inc(end, start);
    std::cout << "Bspline-fitting time: " << duration << " milliseconds" << std::endl;  

    std::cout << "x: " << ptc.x << std::endl;
    std::cout << "y: " << ptc.y << std::endl;

    // -----------------test--------------------
    // Point test_point_;

    // SetHome(HOME_LAT_, HOME_LON_);


    // double x_intr;
    // double y_intr;


    // std::vector<std::vector<double>> LonandLat = {{118.7806222, 31.8351078},
    //                                                 {118.7810634, 31.8347478},
    //                                                 {118.7814473, 31.8344676},
    //                                                 {118.7814108, 31.8348904},
    //                                                 {118.7806701, 31.8351548},
    //                                                 {118.7815429, 31.8352704},
    //                                                 {118.7811188, 31.8351964},
    //                                                 {118.7820834, 31.8346503}
    //                                             };

    // for (int i = 0; i < LonandLat.size(); i++) {
        
    //     double longti = LonandLat[i][0];
    //     double lati = LonandLat[i][1];

    //     x_intr = Lon2M(longti);
    //     y_intr = Lat2M(lati);

    //     test_point_ = bsp_fit.GetBsplinePoint(x_intr,  y_intr);

    //     std::cout << std::fixed << std::setprecision(8) << "x: " << x_intr << " " << "y: " << y_intr << std::endl;
    //     std::cout << std::fixed << std::setprecision(8) << "latitude: " << lati << " " << "longitude: " << longti << " " << "height: " << test_point_.z + 21.48 <<  std::endl;
    // }


    // std::vector<std::vector<Point>> ctr_points_set = bsp_fit.ctr_points_; // 需要修改ctr_points_的权限,目前为private

    // // 提取控制点
    // std::vector<Point> control_points;
    // for (const auto& row : ctr_points_set)
    // {
    //     for (const auto& point : row)
    //     {
    //         control_points.push_back(point);
    //     }
    // }


    // // 可视化拟合点和控制点
    // VisualizePointCloudV2(filtered_ground, control_points);

    return 0;
}
