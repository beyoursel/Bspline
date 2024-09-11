/**
 * \file utils.h
 * \author LeTao (1747524097@qq.com)
 * \brief 工具函数
 * \version 0.1
 * \date 2024-09-11
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef UTILS_H
#define UTILS_H

#include <chrono>
#include <cmath>
#include <iostream>
#include <pcl/common/common.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#include <queue>
#include "nanoflann.hpp"

struct PointT {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Eigen::Vector3f position;
  float vel;
  float rcs;

  bool operator<(const PointT &other) const { return position.x() < other.position.x(); }
};


typedef std::vector<PointT, Eigen::aligned_allocator<PointT>> PointCloud;

// kd_tree adapter
struct PointCloudAdapter {
    const PointCloud& pts;

    PointCloudAdapter(const PointCloud& pointCloud) : pts(pointCloud) {}

    inline size_t kdtree_get_point_count() const { return pts.size(); }

    inline float kdtree_get_pt(const size_t idx, const size_t dim) const {
        if (dim == 0) return pts[idx].position.x();
        else if (dim == 1) return pts[idx].position.y();
        else return pts[idx].position.z();
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};


/**
 * \brief 计算点云数据集的最远最小点
 */

void getMinMax3D(const PointCloud &input, PointT &min_pt, PointT &max_pt);

/**
 * \brief 计算平方差
 *
 * \param sx 查询点x
 * \param sy 查询点y
 * \param fx 拟合点x
 * \param fy 拟合点y
 * \return 返回误差平方
 */
double SquareError(double sx, double sy, double fx, double fy);

/**
 * \brief 计时
 *
 * \param t_end
 * \param t_begin
 * \return double
 */
double time_inc(std::chrono::high_resolution_clock::time_point &t_end,
                std::chrono::high_resolution_clock::time_point &t_begin);

/**
 * \brief 体素滤波
 *
 * \param input_cloud
 * \param cloud_filtered
 * \param voxel_size 体素大小
 */
void VoxelDownSample(const pcl::PointCloud<pcl::PointXYZ>::Ptr &input_cloud,
                     pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud_filtered,
                     float voxel_size);

/**
 * \brief 去离群点
 *
 * \param input_cloud
 * \param cloud_filtered
 * \param mean_k
 * \param std_thresh
 */
void StatisticalRemoveOutlier(
    const pcl::PointCloud<pcl::PointXYZ>::Ptr &input_cloud,
    pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud_filtered, int mean_k,
    float std_thresh);

#endif // UTILS_H