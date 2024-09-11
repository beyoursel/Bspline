#include "bspline_surface.h"
#include <pcl/io/pcd_io.h>
#include <pcl/visualization/pcl_visualizer.h>

void VisualizePointCloudV2(pcl::PointCloud<pcl::PointXYZ>::Ptr downsample_cloud,
                           pcl::PointCloud<pcl::PointXYZ>::Ptr ctp_cloud,
                           pcl::PointCloud<pcl::PointXYZ>::Ptr sur_cloud) {
  pcl::visualization::PCLVisualizer::Ptr viewer(
      new pcl::visualization::PCLVisualizer("3D Viewer"));
  viewer->setBackgroundColor(0, 0, 0);

  pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud_color_handler(
      downsample_cloud, 255, 255, 255);
  viewer->addPointCloud(downsample_cloud, cloud_color_handler,
                        "downsample_points");

  pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ>
      control_color_handler(ctp_cloud, 255, 0, 0);
  viewer->addPointCloud(ctp_cloud, control_color_handler, "ctp_cloud");

  // pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> surf_color_handler(
  //     sur_cloud, 0, 255, 0);
  // viewer->addPointCloud(sur_cloud, surf_color_handler, "sur_cloud");

  viewer->setPointCloudRenderingProperties(
      pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "downsample_points");
  // viewer->setPointCloudRenderingProperties(
  //     pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "sur_cloud");
  viewer->setPointCloudRenderingProperties(
      pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "ctp_cloud");

  viewer->addCoordinateSystem(1.0);
  viewer->initCameraParameters();

  while (!viewer->wasStopped()) {
    viewer->spinOnce(100);
  }
}

int main(int argc, char **argv) {

  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <input_pcd_file>" << std::endl;
    return -1;
  }

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
  if (pcl::io::loadPCDFile<pcl::PointXYZ>(argv[1], *cloud) == -1) {
    PCL_ERROR("Couldn't read file %s \n", argv[1]);
    return -1;
  }

  pcl::PointCloud<pcl::PointXYZ>::Ptr filtered_cloud(new pcl::PointCloud<pcl::PointXYZ>);
  pcl::PointCloud<pcl::PointXYZ>::Ptr filtered_ground(new pcl::PointCloud<pcl::PointXYZ>);

  VoxelDownSample(cloud, filtered_cloud, 1.0);

  StatisticalRemoveOutlier(filtered_cloud, filtered_ground, 50, 1.0);

  std::cout << "cloud size is: " << filtered_ground->points.size() << std::endl;


  PointCloud bsp_cloud;

  for (auto& ptc: filtered_ground->points) {
    PointT bsp_ptc;
    bsp_ptc.position.x() = ptc.x;
    bsp_ptc.position.y() = ptc.y;
    bsp_ptc.position.z() = ptc.z;
    bsp_cloud.push_back(bsp_ptc);
  }


  auto start = std::chrono::high_resolution_clock::now();

  BspSurface bsp_fit(3, 10.0);
  bsp_fit.SetData(bsp_cloud, 0.1);

  auto end = std::chrono::high_resolution_clock::now();

  double ptc = bsp_fit.GetHeight(137.99, -50.81); // test

  auto end1 = std::chrono::high_resolution_clock::now();

  double duration = time_inc(end, start);
  std::cout << "Bspline-fitting Initialization: " << duration << " milliseconds"
            << std::endl;

  double duration1 = time_inc(end1, end);
  std::cout << "Query Bspline-fitting Point: " << duration1 << " milliseconds"
            << std::endl;

  // std::cout << "x: " << ptc.x << std::endl;
  // std::cout << "y: " << ptc.y << std::endl;
  std::cout << "z: " << ptc << std::endl;

  // -----------------visualization--------------------
  pcl::PointCloud<pcl::PointXYZ>::Ptr surf_cloud(new pcl::PointCloud<pcl::PointXYZ>);
  PointCloud bsp_surf_cloud;
  bsp_fit.GetSurface(bsp_surf_cloud);

  for (auto& ptc: bsp_surf_cloud) {
    pcl::PointXYZ temp;
    temp.x = ptc.position.x();
    temp.y = ptc.position.y();
    temp.z = ptc.position.z();    
    surf_cloud->points.push_back(temp);
  }
  pcl::PointCloud<pcl::PointXYZ>::Ptr ctr_cloud(new pcl::PointCloud<pcl::PointXYZ>);

  PointCloud bsp_ctr_cloud = bsp_fit.GetCtrlPts();
  for (auto& ptc: bsp_ctr_cloud) {
    pcl::PointXYZ temp;
    temp.x = ptc.position.x();
    temp.y = ptc.position.y();
    temp.z = ptc.position.z();    
    ctr_cloud->points.push_back(temp);
  }
  
  VisualizePointCloudV2(filtered_ground, ctr_cloud, surf_cloud);

  return 0;
}