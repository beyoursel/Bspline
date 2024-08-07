#include "bspline_api.h"



Point BspFitting::GetBsplinePoint(double x, double y) {

        if ((x < min_pt_.x || x > max_pt_.x || y < min_pt_.y || y > max_pt_.y)) {
            std::cerr << "the query point is out of the range" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        return bsp_.GetFittingPoint(x, y);
    }


void BspFitting::GetControlPoint(PointCloud cloud) {


    int M = static_cast<int>(x_range_ / grid_size_) + 2 + 2 * k_ex_;
    int N = static_cast<int>(y_range_ / grid_size_) + 2 + 2 * k_ex_;   

    double min_limit_x = min_pt_.x - k_ex_ * grid_size_;
    double min_limit_y = min_pt_.y - k_ex_ * grid_size_;

    std::vector<std::vector<Grid_T>> grids(M, std::vector<Grid_T>(N));

    for (size_t i = 0; i < cloud->points.size(); ++i) {
        const auto& point = cloud->points[i];
        int grid_x = static_cast<int>((point.x - min_pt_.x) / grid_size_) + k_ex_;
        int grid_y = static_cast<int>((point.y - min_pt_.y) / grid_size_) + k_ex_;

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
                ctr_ptc[i][j] =  Point((i + 0.5) * grid_size_ + min_limit_x, (j + 0.5) * grid_size_ + min_limit_y, grid_points[retain_count-1].z);
                ptc_list.push_back(ctr_ptc[i][j]);
                // mean height of 10% pointcloud at bottom
                // double height_sum = 0.0;
                // for (auto it = grid_points.begin(); it != grid_points.begin() + retain_count; ++it) {
                //     height_sum += it->z;
                // }
                // ctr_ptc[i][j] = Point((i + 0.5) * grid_size_ + min_pt_.x, (j + 0.5) * grid_size_ + min_pt_.y, height_sum/retain_count));
            } else {
                ctr_ptc[i][j] = Point((i + 0.5) * grid_size_ + min_limit_x, (j + 0.5) * grid_size_ + min_limit_y, MAX_DOUBLE);
            }

        }
    }

    // interpolate the empty grid(z-value is maxDoudble) with neareast ground points
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            if (ctr_ptc[i][j].z == MAX_DOUBLE) {
                double x_em = ctr_ptc[i][j].x; // convert to real coordiantes
                double y_em = ctr_ptc[i][j].y;
                ctr_ptc[i][j].z = InterpolatePoint(ptc_list, x_em, y_em);
            }
        }
    }

    ctr_points_ = ctr_ptc;
}


double BspFitting::InterpolatePoint(const std::vector<Point>& points, double x, double y, int k) {
    if (points.empty()) {
        throw std::invalid_argument("Point list is empty.");
    }

    std::priority_queue<std::pair<double, Point>> pq;

    // search the nearest k points
    for (const auto& pt : points) {
        double dist = std::sqrt((pt.x - x) * (pt.x - x) + (pt.y - y) * (pt.y - y));
        pq.push(std::make_pair(dist, pt));
        if (pq.size() > k) { 
            pq.pop();
        } 
    }

    double sum_z = 0.0;
    int count = 0;
    while (!pq.empty()) {
        auto top = pq.top();
        pq.pop();
        sum_z += top.second.z;
        count++;
    }
    return sum_z / count;
}