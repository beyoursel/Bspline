#include "bspline_api.h"


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

    k_ex_ = 1;

    GetControlPoint(ptc);
    bsp_ = BspSurface(ctr_points_, k_);
}


Point BspFitting::GetBsplinePoint(double x, double y) {

        if ((x < min_pt_.x || x > max_pt_.x || y < min_pt_.y || y > max_pt_.y)) {
            std::cerr << "the query point is out of the range" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        return bsp_.GetFittingPoint(x, y);
    }


void BspFitting::GetControlPoint(PointCloud cloud) {

    auto start = std::chrono::high_resolution_clock::now();
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

    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() * 1000;
    std::cout << "Assign PointCloud to Grid Time: " << duration << " milliseconds" << std::endl;  

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

            } else {
                ctr_ptc[i][j] = Point((i + 0.5) * grid_size_ + min_limit_x, (j + 0.5) * grid_size_ + min_limit_y, MAX_DOUBLE);
            }

        }
    }

    std::cout << "valid ctr-points from ground pointcloud: " << ptc_list.size() << std::endl;

    auto end1 = std::chrono::high_resolution_clock::now();
    double duration1 = std::chrono::duration_cast<std::chrono::duration<double>>(end1 - end).count() * 1000;
    std::cout << "Generate Ground Points Time: " << duration1 << " milliseconds" << std::endl;  

    // version 1 kd-tree

    // InterpolatePointKd(ctr_ptc, ptc_list);

    // version 2 adaptive radius search

    InterpolatePointAdaptRadius(ctr_ptc);

    // version 3
    // InterpolatePointAdaptRadiusQueue(ctr_ptc);

    auto end2 = std::chrono::high_resolution_clock::now();
    double duration2 = std::chrono::duration_cast<std::chrono::duration<double>>(end2 - end1).count() * 1000;
    std::cout << "Interpolate Empty Position Time: " << duration2 << " milliseconds" << std::endl;  

    ctr_points_ = ctr_ptc;
}


void BspFitting::InterpolatePointAdaptRadiusQueue(std::vector<std::vector<Point>>& temp, int k_nei) {

    std::vector<std::vector<Point>> ctr_ptc = temp;
    int M = ctr_ptc.size();
    int N = ctr_ptc[0].size();

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            if (std::isinf(ctr_ptc[i][j].z)) {
                std::priority_queue<std::pair<double, double>, std::vector<std::pair<double, double>>, Compare> pq;
                int searchRadius = 3;

                // 搜索整个半径区域内的所有点
                while (pq.size() < 3 && searchRadius < std::max(M, N)) {
                    for (int di = -searchRadius; di <= searchRadius; ++di) {
                        for (int dj = -searchRadius; dj <= searchRadius; ++dj) {
                            int ni = i + di;
                            int nj = j + dj;

                            if (ni >= 0 && ni < M && nj >= 0 && nj < N && !std::isinf(ctr_ptc[ni][nj].z)) {
                                double distSq = DistanceXY(ctr_ptc[i][j], ctr_ptc[ni][nj]);
                                pq.emplace(distSq, ctr_ptc[ni][nj].z);
                            }
                        }
                    }
                    searchRadius++;
                }

                // 检查找到的所有点，保留最小的k_nei个点
                std::vector<double> topKNearestZ;
                while (!pq.empty() && topKNearestZ.size() < k_nei) {
                    topKNearestZ.push_back(pq.top().second);
                    pq.pop();
                }

                // 计算k_nei个最近点的z值的平均值
                double sumz = 0.0;
                for (const double z : topKNearestZ) {
                    sumz += z;
                }
                temp[i][j].z = sumz / k_nei;
            }
        }
    }
}



void BspFitting::InterpolatePointAdaptRadius(std::vector<std::vector<Point>>& temp, int k_nei) {
   
    std::vector<std::vector<Point>> ctr_ptc = temp;
    int M = ctr_ptc.size();
    int N = ctr_ptc[0].size();

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            if (std::isinf(ctr_ptc[i][j].z)) {
                std::vector<Point> nearpoint;
                std::vector<std::pair<double, double>> nearby_h; 

                int searchRadius = 3; 

                // adaptive expand search radius
                while (nearpoint.size() < 3 && searchRadius < std::max(M, N)) {
                    for (int di = -searchRadius; di <= searchRadius; ++di) {
                        for (int dj = -searchRadius; dj <= searchRadius; ++dj) {
                            int ni = i + di;
                            int nj = j + dj;

                            if (ni >= 0 && ni < M && nj >= 0 && nj < N && !std::isinf(ctr_ptc[ni][nj].z)) {
                                nearpoint.push_back(ctr_ptc[ni][nj]);
                            }
                        }
                    }
                    searchRadius = searchRadius + 1; 
                }

                for (auto ptc: nearpoint) {
                    double dis = DistanceXY(ctr_ptc[i][j], ptc);
                    nearby_h.emplace_back(dis, ptc.z);
                }
                std::sort(nearby_h.begin(), nearby_h.end());

                double sumz = 0.0;

                for (int i = 0; i < k_nei; i++) {
                    sumz += nearby_h[i].second;
                }
                temp[i][j].z = sumz / k_nei;
            }
        }
    }
}


void BspFitting::InterpolatePointKd(std::vector<std::vector<Point>>& cn_points, std::vector<Point>& ptc_lists, int k_kd) {

    PointCloud cloud_2d_presudo(new pcl::PointCloud<pcl::PointXYZ>);
    for (auto point: ptc_lists) {
        cloud_2d_presudo->points.push_back(pcl::PointXYZ(point.x, point.y, 0.0f));
    }

    auto start_kd = std::chrono::high_resolution_clock::now();
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(cloud_2d_presudo);

    auto end_kd = std::chrono::high_resolution_clock::now();
    double duration_kd = std::chrono::duration_cast<std::chrono::duration<double>>(end_kd - start_kd).count() * 1000;
    std::cout << "Build Kd_tree Time: " << duration_kd << " milliseconds" << std::endl;  

    std::vector<int> point_idx_knsearch(k_kd);
    std::vector<float> point_kdn_sq_dis(k_kd);
    for (int i = 0; i < cn_points.size(); ++i) {
        for (int j = 0; j < cn_points[0].size(); ++j) {
            pcl::PointXYZ search_point;
            if (std::isinf(cn_points[i][j].z)) {
                search_point.x = cn_points[i][j].x; // convert to real coordiantes
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