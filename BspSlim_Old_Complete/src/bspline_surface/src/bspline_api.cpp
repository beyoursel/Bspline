#include "bspline_api.h"


double BspFitting::Lat2M(double latitude)
{
  return latitude / 0.00000899 - home_y_;
}


double BspFitting::Lon2M(double longitude)
{
  return longitude / (360.0 / (40030173.0 * cos(HOME_LAT_ * M_PI / 180.0))) - home_x_;
}


void BspFitting::SetHome(double latitude, double longitude) {

  home_x_ = longitude / (360.0 / (40030173.0 * cos(latitude * M_PI / 180.0)));
  home_y_ = latitude / 0.00000899;
}


double BspFitting::DistanceXY(const Point& p1, const Point& p2) {
    return std::sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}


BspFitting::BspFitting(PointCloud ptc, double grid_size, int k)
{
    k_ = k;
    grid_size_ = grid_size;
    pcl::getMinMax3D(*ptc, min_pt_, max_pt_);
    x_range_ = max_pt_.x - min_pt_.x;
    y_range_ = max_pt_.y - min_pt_.y;

    // std::cout << "x_range is: " << min_pt_.x << " " << max_pt_.x << std::endl;
    // std::cout << "y_range is: " << min_pt_.y << " " << max_pt_.y << std::endl;    
    // k_ex_ = 1;

    GetControlPoint(ptc);
    bsp_ = BspSurface(ctr_points_, k_);
    cn_point_pub_ = bsp_.m_cn_point_; // for vis
    k_point_pub_ = bsp_.knot_cn_point_; // for_vis
}


Point BspFitting::GetBsplinePoint(double x, double y) {

    if ((x < min_pt_.x || x > max_pt_.x || y < min_pt_.y || y > max_pt_.y)) { // 修改成knot_range
        std::cerr << "the query point is out of the range" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return bsp_.GetFittingPoint(x, y);
}


double BspFitting::GetBsplinePointLonLat(double lon, double lat) {

    SetHome(HOME_LAT_, HOME_LON_); 
    double x_intr = Lon2M(lon);
    double y_intr = Lat2M(lat);

    Point query_p = GetBsplinePoint(x_intr, y_intr);

    return query_p.z + HOME_LATITUDE;
}

void BspFitting::GetControlPoint(PointCloud cloud) {

    auto start = std::chrono::high_resolution_clock::now();
    int M = static_cast<int>(x_range_ / grid_size_) + 1;
    int N = static_cast<int>(y_range_ / grid_size_) + 1;   

    std::vector<std::vector<Grid_T>> grids(M, std::vector<Grid_T>(N));

    for (size_t i = 0; i < cloud->points.size(); ++i) {
        const auto& point = cloud->points[i];
        int grid_x = static_cast<int>((point.x - min_pt_.x) / grid_size_);
        int grid_y = static_cast<int>((point.y - min_pt_.y) / grid_size_);

        if (grid_x >= 0 && grid_x < M && grid_y >= 0 && grid_y < N) {

            grids[grid_x][grid_y].points.push_back({point.x, point.y, point.z});
        }
    }

    std::vector<std::vector<Point>> ctr_ptc(M, std::vector<Point>(N));

    std::vector<Point> ptc_list;

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            auto& grid_points = grids[i][j].points;

            if (!grid_points.empty()) {


                std::sort(grid_points.begin(), grid_points.end(), [](Point a, Point b) { return a.z < b.z; });

                size_t retain_count = static_cast<size_t>(grid_points.size() * 0.1);

                if (retain_count < 1) {
                    retain_count = 1;
                }

                // one point at the bottom of 10%
                ctr_ptc[i][j] =  Point((i + 0.5) * grid_size_ + min_pt_.x, (j + 0.5) * grid_size_ + min_pt_.y, grid_points[retain_count-1].z);
                ptc_list.push_back(ctr_ptc[i][j]);

            } else {
                ctr_ptc[i][j] = Point((i + 0.5) * grid_size_ + min_pt_.x, (j + 0.5) * grid_size_ + min_pt_.y, MAX_DOUBLE);
            }

        }
    }


    InterpolatePointKd(ctr_ptc, ptc_list);

    ctr_points_ = ctr_ptc;
}


// void BspFitting::GetControlPoint(PointCloud cloud) {

//     auto start = std::chrono::high_resolution_clock::now();
//     int M = static_cast<int>(x_range_ / grid_size_) + 2 * k_ex_;
//     int N = static_cast<int>(y_range_ / grid_size_) + 2 * k_ex_;   

//     double min_limit_x = min_pt_.x - k_ex_ * grid_size_;
//     double min_limit_y = min_pt_.y - k_ex_ * grid_size_;

//     std::vector<std::vector<Grid_T>> grids(M, std::vector<Grid_T>(N));

//     for (size_t i = 0; i < cloud->points.size(); ++i) {
//         const auto& point = cloud->points[i];
//         int grid_x = static_cast<int>((point.x - min_pt_.x) / grid_size_) + k_ex_;
//         int grid_y = static_cast<int>((point.y - min_pt_.y) / grid_size_) + k_ex_;

//         if (grid_x >= 0 && grid_x < M && grid_y >= 0 && grid_y < N) {

//             grids[grid_x][grid_y].points.push_back({point.x, point.y, point.z});
//         }
//     }

//     std::vector<std::vector<Point>> ctr_ptc(M, std::vector<Point>(N));

//     std::vector<Point> ptc_list;

//     for (int i = 0; i < M; ++i) {
//         for (int j = 0; j < N; ++j) {
//             auto& grid_points = grids[i][j].points;

//             if (!grid_points.empty()) {


//                 std::sort(grid_points.begin(), grid_points.end(), [](Point a, Point b) { return a.z < b.z; });

//                 size_t retain_count = static_cast<size_t>(grid_points.size() * 0.1);

//                 if (retain_count < 1) {
//                     retain_count = 1;
//                 }

//                 // one point at the bottom of 10%
//                 ctr_ptc[i][j] =  Point((i + 0.5) * grid_size_ + min_limit_x, (j + 0.5) * grid_size_ + min_limit_y, grid_points[retain_count-1].z);
//                 ptc_list.push_back(ctr_ptc[i][j]);

//             } else {
//                 ctr_ptc[i][j] = Point((i + 0.5) * grid_size_ + min_limit_x, (j + 0.5) * grid_size_ + min_limit_y, MAX_DOUBLE);
//             }

//         }
//     }


//     InterpolatePointKd(ctr_ptc, ptc_list);

//     ctr_points_ = ctr_ptc;
// }


void BspFitting::InterpolatePointKd(std::vector<std::vector<Point>>& cn_points, std::vector<Point>& ptc_lists, int k_kd) {

    PointCloud cloud_2d_presudo(new pcl::PointCloud<pcl::PointXYZ>);
    for (auto point: ptc_lists) {
        cloud_2d_presudo->points.push_back(pcl::PointXYZ(point.x, point.y, 0.0f));
    }

    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(cloud_2d_presudo);

    std::vector<int> point_idx_knsearch(k_kd);
    std::vector<float> point_kdn_sq_dis(k_kd);
    for (int i = 0; i < cn_points.size(); ++i) {
        for (int j = 0; j < cn_points[0].size(); ++j) {
            pcl::PointXYZ search_point;
            if (std::isinf(cn_points[i][j].z)) {
                search_point.x = cn_points[i][j].x;
                search_point.y = cn_points[i][j].y;
                search_point.z = 0.0;

                if (kdtree.nearestKSearch(search_point, k_kd, point_idx_knsearch, point_kdn_sq_dis) > 0) {
                    double height_sum = 0;
                    for (size_t idx = 0; idx < point_idx_knsearch.size(); ++idx) {
                      height_sum += ptc_lists[point_idx_knsearch[idx]].z;
                    }
                    cn_points[i][j].z = height_sum / k_kd;
                }
            }
        }
    }
}