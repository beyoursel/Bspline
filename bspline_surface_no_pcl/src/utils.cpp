#include "utils.h"


void getMinMax3D(const PointCloud& input, PointT& min_pt, PointT& max_pt) {
    // initial max&min
    min_pt.position = Eigen::Vector3f(std::numeric_limits<float>::max(),
                                      std::numeric_limits<float>::max(),
                                      std::numeric_limits<float>::max());
    max_pt.position = Eigen::Vector3f(std::numeric_limits<float>::lowest(),
                                      std::numeric_limits<float>::lowest(),
                                      std::numeric_limits<float>::lowest());
    
    
    for (const auto& point : input) {
        min_pt.position = min_pt.position.cwiseMin(point.position);  
        max_pt.position = max_pt.position.cwiseMax(point.position);
    }
}


double SquareError(double sx, double sy, double fx, double fy) {
  return std::sqrt((fx - sx) * (fx - sx) + (fy - sy) * (fy - sy));
}

double time_inc(std::chrono::high_resolution_clock::time_point &t_end,
                std::chrono::high_resolution_clock::time_point &t_begin) {

  return std::chrono::duration_cast<std::chrono::duration<double>>(t_end -
                                                                   t_begin)
             .count() *
         1000;
}

void VoxelDownSample(const pcl::PointCloud<pcl::PointXYZ>::Ptr &input_cloud,
                     pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud_filtered, float voxel_size) {
  // Create the filtering object
  pcl::VoxelGrid<pcl::PointXYZ> sor;

  sor.setInputCloud(input_cloud);
  sor.setLeafSize(voxel_size, voxel_size, voxel_size); // voxel_size
  // sor.setLeafSize(0.5f, 0.5f, 0.5f);  // voxel_size
  sor.filter(*cloud_filtered);
}

void StatisticalRemoveOutlier(const pcl::PointCloud<pcl::PointXYZ>::Ptr &input_cloud,
                              pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud_filtered, int mean_k,
                              float std_thresh) {
  // Create the filtering object
  pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor;
  sor.setInputCloud(input_cloud);
  sor.setMeanK(mean_k);               // default: 10 越大越严格
  sor.setStddevMulThresh(std_thresh); // default: 2.0 越小越严格
  sor.filter(*cloud_filtered);
}