#ifndef BSPLINE_API_H
#define BSPLINE_API_H

#include "bspline.h"
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/common.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <vector>
#include <queue>
#include <chrono>
#include <algorithm>
#include <cmath>


#define M_PI 3.14159265358979323846
#define HOME_LAT_ 31.8351078
#define HOME_LON_ 118.7806222
#define HOME_LATITUDE 21.48


constexpr double MAX_DOUBLE = std::numeric_limits<double>::infinity();
typedef pcl::PointCloud<pcl::PointXYZ>::Ptr PointCloud;



struct Grid_T {
    std::vector<Point> points;
};


class BspFitting
{
public:
    std::vector<std::vector<Point>> cn_point_pub_; 

/**
 *  \brief Construct a new BspFitting object
 *  \param PointCloud 输入点云
 *  \param grid_size 栅格尺寸
 *  \param k 阶数
 */
    BspFitting(PointCloud ptc, double grid_size, int k);

/**
 *  \brief 查询(x, y)位置的z
 *  \param x
 *  \param y
 *  \return Point
 */
    Point GetBsplinePoint(double x, double y);

/**
 *  \brief 输入经纬度，查询对应位置的z
 *  \param lon
 *  \param lat
 *  \return z
 */
    double GetBsplinePointLonLat(double lon, double lat);


private:
    BspSurface bsp_; // BspSurface实例
    pcl::PointXYZ min_pt_, max_pt_; // 输入点云最近、最远点
    double grid_size_; // 栅格尺寸
    double x_range_, y_range_; // 输入点云范围
    int k_, k_ex_; // 阶数，控制点扩展维度
    double home_x_;
    double home_y_;

/**
 *  \brief 计算XY平面上两点距离
 */
    double DistanceXY(const Point& p1, const Point& p2);
    std::vector<std::vector<Point>> ctr_points_;

/**
 *  \brief 根据输入点云计算控制点
 *  \param cloud
 */
    void GetControlPoint(PointCloud cloud);

/**
 *  \brief Kd-tree版本控制点插值
 *  \param ctr_ptc 控制点
 *  \param ptc_lists 由点云得到的控制点
 *  \param k_kd 最近邻点个数
 */
    void InterpolatePointKd(std::vector<std::vector<Point>>& ctr_ptc, std::vector<Point>& ptc_lists, int k_kd = 3);

/**
 *  \brief 将纬度转化为ENU下局部坐标y
 *  \param latitude 纬度
 *  \return 返回y
 */
    double Lat2M(double latitude);

/**
 *  \brief 将经度转化为ENU下局部坐标x
 *  \param longtitude 经度
 *  \return 返回x
 */  
    double Lon2M(double longitude);

/**
 *  \brief 设置home点
 *  \param lati_home home点的纬度
 *  \param lon_home home点的经度
 */
    void SetHome(double lati_home, double lon_home);

};

#endif
