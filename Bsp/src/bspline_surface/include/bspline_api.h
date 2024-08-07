#ifndef BSPLINE_API_H
#define BSPLINE_API_H

#include "bspline.h"
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/common.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <vector>
#include <queue>


constexpr double MAX_DOUBLE = std::numeric_limits<double>::max();
typedef pcl::PointCloud<pcl::PointXYZ>::Ptr PointCloud;


struct Grid_T {
    std::vector<Point> points;
};


class BspFitting
{
public:
    BspFitting(PointCloud ptc, double grid_size, int k): k_(k), grid_size_(grid_size) {
        
        pcl::getMinMax3D(*ptc, min_pt_, max_pt_);
        x_range_ = max_pt_.x - min_pt_.x;
        y_range_ = max_pt_.y - min_pt_.y;

        k_ex_ = 2; // 在control_points外围增加一圈


        GetControlPoint(ptc); // 根据输入点云获得control_points
        bsp_ = BspSurface(ctr_points_, k_); // 初始化BspSurface实例
    }

    Point GetBsplinePoint(double x, double y);
    double InterpolatePoint(const std::vector<Point>& points, double x, double y, int k = 3);
    void InterpolatePointKd(std::vector<std::vector<Point>>& cn_points, std::vector<Point>& ptc_lists, int k_kd = 3);



private:
    BspSurface bsp_;
    pcl::PointXYZ min_pt_, max_pt_;
    double grid_size_;
    double x_range_, y_range_;
    int k_, k_ex_;
    std::vector<std::vector<Point>> ctr_points_;
    void GetControlPoint(PointCloud cloud);
    
};


#endif
