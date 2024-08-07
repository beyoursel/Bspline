#include <iostream>
#include <vector>
#include <algorithm>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/common/common.h>
#include <chrono>
#include <boost/filesystem.hpp>
#include <omp.h>
#include <bspline.h>
#include <queue>
#include <limits>
#include <cmath>
#include <iomanip>


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


// 将pcl::PointCloud<pcl::PointXYZ>::Ptr转换为std::vector<Point>
std::vector<Point> PclPointCloudToVector(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud) {
    std::vector<Point> points;

    for (const auto& pt : cloud->points) {
        Point p;
        p.x = pt.x;
        p.y = pt.y;
        p.z = pt.z;
        points.push_back(p);
    }

    return points;
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


double InterpolatePointKdTree(const pcl::KdTreeFLANN<Point>& kdtree, const pcl::PointCloud<Point>::Ptr& cloud, double x, double y, int k = 3) {
    std::vector<int> pointIdxKNNSearch(k);
    std::vector<float> pointKNNSquaredDistance(k);

    Point searchPoint;
    searchPoint.x = x;
    searchPoint.y = y;

    if (kdtree.nearestKSearch(searchPoint, k, pointIdxKNNSearch, pointKNNSquaredDistance) > 0) {
        double sum_z = 0.0;
        for (int idx : pointIdxKNNSearch) {
            sum_z += cloud->points[idx].z;
        }
        return sum_z / k;
    } else {
        throw std::runtime_error("KdTree search failed.");
    }
}


// 插值函数：使用KNN（k近邻）算法来填充空栅格
double InterpolatePoint(const std::vector<Point>& points, double x, double y, int k = 3) {
    if (points.empty()) {
        throw std::invalid_argument("Point list is empty.");
    }
    // std::priority_queue是一个最大堆，队列顶部(pq.top())存储的是队列中最大的元素
    std::priority_queue<std::pair<double, Point>> pq; // std::pair的规则是先比较第一个元素，若第一个元素相同则比较第二个元素

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


struct Grid {
    std::vector<Point> points;
};


bool compareHeight(const Point& a, const Point& b) {
    return a.z < b.z;
}


void StatisticalRemoveOutlier(pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud , pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered, int mean_k, float std_thresh) {
      // Create the filtering object
    pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor;
    sor.setInputCloud (input_cloud);
    sor.setMeanK(mean_k); // default: 10 越大越严格
    sor.setStddevMulThresh (std_thresh); // default: 2.0 越小越严格
    sor.filter(*cloud_filtered);
}


void PassthroghFilter(pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud , pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered) {
    // Create the filtering object
    pcl::PassThrough<pcl::PointXYZ> pass;
    pass.setInputCloud (input_cloud);
    pass.setFilterFieldName ("z");
    pass.setFilterLimits (-5, 20);
    //pass.setNegative (true);
    pass.filter (*cloud_filtered);
    // ROS_INFO("the number of passthrough ptc is %ld", cloud_filtered->size()); 
    std::cout << "the number of passthrough ptc is "  << cloud_filtered->size() << std::endl;
}


void VoxelDownSample(pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud , pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered, float voxel_size) {
    // Create the filtering object
    pcl::VoxelGrid<pcl::PointXYZ> sor;

    sor.setInputCloud(input_cloud);
    sor.setLeafSize(voxel_size, voxel_size, voxel_size);  // voxel_size     
    // sor.setLeafSize(0.5f, 0.5f, 0.5f);  // voxel_size 
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

    // PassthroghFilter(cloud, cloud);

    VoxelDownSample(cloud, filtered_cloud, 1.0);

    pcl::PointXYZ min_pt, max_pt;
    pcl::getMinMax3D(*filtered_cloud, min_pt, max_pt); //


    double x_range = max_pt.x - min_pt.x;
    double y_range = max_pt.y - min_pt.y;
    // grid_size
    double grid_width = 5.0;
    double grid_height = 5.0;

    int M = static_cast<int>(x_range / grid_width) + 2;
    int N = static_cast<int>(y_range / grid_height) + 2;
     

    std::vector<std::vector<Grid>> grids(M, std::vector<Grid>(N));

    for (size_t i = 0; i < filtered_cloud->points.size(); ++i) {
        const auto& point = filtered_cloud->points[i];
        int grid_x = static_cast<int>((point.x - min_pt.x) / grid_width);
        int grid_y = static_cast<int>((point.y - min_pt.y) / grid_height);

        if (grid_x >= 0 && grid_x < M && grid_y >= 0 && grid_y < N) {

            grids[grid_x][grid_y].points.push_back({point.x, point.y, point.z});
        }
    }

    pcl::PointCloud<pcl::PointXYZ>::Ptr ground_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr ground_filtered(new pcl::PointCloud<pcl::PointXYZ>);    
    pcl::PointCloud<pcl::PointXYZ>::Ptr tree_top(new pcl::PointCloud<pcl::PointXYZ>);

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            auto& grid_points = grids[i][j].points;
            pcl::PointCloud<pcl::PointXYZ>::Ptr grid_cloud_outliers(new pcl::PointCloud<pcl::PointXYZ>);
            if (!grid_points.empty()) {


                std::sort(grid_points.begin(), grid_points.end(), compareHeight);

                size_t retain_count = static_cast<size_t>(grid_points.size() * 0.1);
                size_t range_count = static_cast<size_t>(grid_points.size() * 0.05);

                if (retain_count < 1) {
                    retain_count = 1;
                }

                size_t retain_top = static_cast<size_t>(grid_points.size() * 0.07);
                if (retain_top < 1) {
                    retain_top = 1;
                }


                // one point at the bottom of 10%
                ground_cloud->points.push_back(pcl::PointXYZ((i + 0.5) * grid_width + min_pt.x, (j + 0.5) * grid_height + min_pt.y, grid_points[retain_count-1].z));

                for (auto it = grid_points.begin()+range_count+range_count; it != grid_points.end(); ++it) {
                    grid_cloud_outliers->points.push_back(pcl::PointXYZ(it->x, it->y, it->z));
                }

                for (auto it = grid_points.begin(); it != grid_points.begin() + retain_count - range_count; ++it) {
                    grid_cloud_outliers->points.push_back(pcl::PointXYZ(it->x, it->y, it->z));
                }

                if (grid_cloud_outliers->points.size() > 50) {
                    *tree_top += *grid_cloud_outliers;
                }

                // mean height of 10% pointcloud at bottom
                // float height_sum = 0.0;
                // for (auto it = grid_points.begin(); it != grid_points.begin() + retain_count; ++it) {
                //     height_sum += it->z;
                // }
                // ground_cloud->points.push_back(pcl::PointXYZ((i + 0.5) * grid_width + min_pt.x, (j + 0.5) * grid_height + min_pt.y, height_sum/retain_count));
            }
        }
    }

    std::cout << "top number: " << tree_top->points.size() << std::endl;

    std::cout << "before statistical: " << ground_cloud->points.size() << std::endl;
    StatisticalRemoveOutlier(ground_cloud, ground_filtered, 50, 0.35);
    std::cout << "after statistical: " << ground_filtered->points.size() << std::endl;


    // obtain the height difference of bspline-fitting surface
    pcl::PointXYZ min_pt_g, max_pt_g;
    pcl::getMinMax3D(*ground_filtered, min_pt_g, max_pt_g);


    // print the range of ground point
    // std::cout << "min_pt_g.x is: " << min_pt_g.x << std::endl;
    // std::cout << "max_pt_g.x is: " << max_pt_g.x << std::endl;
    // std::cout << "min_pt_g.y is: " << min_pt_g.y << std::endl;
    // std::cout << "max_pt_g.y is: " << max_pt_g.y << std::endl;


    // print the range of raw pointcloud
    std::cout << "min_pt.x is: " << min_pt.x << std::endl;
    std::cout << "max_pt.x is: " << max_pt.x << std::endl;
    std::cout << "min_pt.y is: " << min_pt.y << std::endl;
    std::cout << "max_pt.y is: " << max_pt.y << std::endl;

    // ***********************************version 1*************************************
    // cover all raw pointcloud

    // std::vector<std::vector<Eigen::Vector3d>> cn_point(M+2, std::vector<Eigen::Vector3d>(N+2));
    // double maxDouble = std::numeric_limits<double>::max();


    // // create control points grid
    // for (size_t i = 0; i < M+2; i++) {
    //     for (size_t j = 0; j < N+2; j++) {
    //         // cn_point[i][j] = Eigen::Vector3d(min_pt.x + i * grid_width, min_pt.y + j * grid_height, low_peak);
    //         cn_point[i][j](0) = min_pt.x + i * grid_width - grid_width;
    //         cn_point[i][j](1) = min_pt.y + j * grid_height - grid_height;
    //         cn_point[i][j](2) = maxDouble;
    //     }
    // }

    // for (const auto& point : ground_filtered->points) {
    //     float x_ptc = static_cast<float>(point.x);
    //     float y_ptc = static_cast<float>(point.y);
    //     cn_point[(x_ptc - min_pt.x) / grid_width + 1][(y_ptc - min_pt.y) / grid_height + 1](2) = point.z;
    // }

    // std::vector<Point> ptc = PclPointCloudToVector(ground_filtered);
    // // // 插值填充空栅格
    // for (int i = 0; i < M+2; ++i) {
    //     for (int j = 0; j < N+2; ++j) {
    //         if (cn_point[i][j](2) == maxDouble) {
    //             double x_em = min_pt.x + i * grid_width - grid_width;
    //             double y_em = min_pt.y + j * grid_height - grid_height;
    //             cn_point[i][j](2) = interpolatePoint(ptc, x_em, y_em);
    //         }
    //     }
    // }
     
    // ***********************************version 2********************************************** 
    // ***cover all the extracted ground point
    // ***create control points grid

    double x_range_g = max_pt_g.x - min_pt_g.x;
    double y_range_g = max_pt_g.y - min_pt_g.y;
    int M_g = static_cast<int>(x_range_g / grid_width) + 1 + 2; // x_range_g is int-type
    int N_g = static_cast<int>(y_range_g / grid_height) + 1 + 2;

    std::vector<std::vector<Point>> cn_point(M_g, std::vector<Point>(N_g));
    double maxDouble = std::numeric_limits<double>::max();

    std::cout << x_range_g << " " << y_range_g << std::endl;

    for (size_t i = 0; i < M_g; i++) {
        for (size_t j = 0; j < N_g; j++) {
            // cn_point[i][j] = Eigen::Vector3d(min_pt.x + i * grid_width, min_pt.y + j * grid_height, low_peak);
            cn_point[i][j].x = min_pt_g.x + i * grid_width - grid_width * 1;
            cn_point[i][j].y = min_pt_g.y + j * grid_height - grid_height * 1;
            cn_point[i][j].z = maxDouble;
        }
    }

    for (const auto& point : ground_filtered->points) {
        float x_ptc = static_cast<float>(point.x);
        float y_ptc = static_cast<float>(point.y);
        cn_point[(x_ptc - min_pt_g.x) / grid_width + 1][(y_ptc - min_pt_g.y) / grid_height + 1].z = point.z;
    }


    // version 1
    auto start_1 = std::chrono::high_resolution_clock::now();

    // std::vector<Point> ptc = PclPointCloudToVector(ground_filtered); // convert ground_cloud to vector<point>


    // interpolate the empty grid(z-value is maxDoudble) with neareast ground points
    // for (int i = 0; i < M_g; ++i) {
    //     for (int j = 0; j < N_g; ++j) {
    //         if (cn_point[i][j].z == maxDouble) {
    //             double x_em = min_pt_g.x + i * grid_width - grid_width * 1; // convert to real coordiantes
    //             double y_em = min_pt_g.y + j * grid_height - grid_height * 1;
    //             cn_point[i][j].z = InterpolatePoint(ptc, x_em, y_em);
    //         }
    //     }
    // }


    // version 2
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_2d_presudo(new pcl::PointCloud<pcl::PointXYZ>);
    for (auto point: ground_filtered->points) {
        cloud_2d_presudo->points.push_back(pcl::PointXYZ(point.x, point.y, 0.0f));
    }
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(cloud_2d_presudo);

    int k_kd = 3;
    std::vector<int> point_idx_knsearch(k_kd);
    std::vector<float> point_kdn_sq_dis(k_kd);
    for (int i = 0; i < M_g; ++i) {
        for (int j = 0; j < N_g; ++j) {
            pcl::PointXYZ search_point;
            if (cn_point[i][j].z == maxDouble) {
                search_point.x = min_pt_g.x + i * grid_width - grid_width * 1; // convert to real coordiantes
                search_point.y = min_pt_g.y + j * grid_height - grid_height * 1;
                search_point.z = 0.0;

                if (kdtree.nearestKSearch(search_point, k_kd, point_idx_knsearch, point_kdn_sq_dis) > 0) {
                    double height_sum = 0;
                    for (size_t idx = 0; idx < point_idx_knsearch.size(); ++idx) {
                      height_sum += ground_filtered->points[point_idx_knsearch[idx]].z;
                    }
                    cn_point[i][j].z = height_sum / k_kd;
                    // std::cout << "complete kdtree search" << std::endl;
                }
            }
        }
    }

    auto end_1 = std::chrono::high_resolution_clock::now();
    double duration = time_inc(end_1, start_1);
    std::cout << "Bspline-fitting time: " << duration << " milliseconds" << std::endl;  

    // *************************test Bspline*******************************

    // create B-spline instance
    
    int k = 3; // 阶数
    BspSurface surface(cn_point, k);

    // create the Bspline surface
    // std::vector<Point> vertices;
    // surface.GetFittingSurface(vertices, 0.05); // 0.05为step
    

    std::vector<Point> test_points_set;
    Point test_point_;

    // double x_intr = std::stod(argv[2]);
    // double y_intr = std::stod(argv[3]);

    // double longti = std::stod(argv[2]);
    // double lati = std::stod(argv[3]);
    // convert the takeoff point to 
    SetHome(HOME_LAT_, HOME_LON_);


    double x_intr;
    double y_intr;

    // x_intr = Lon2M(longti);
    // y_intr = Lat2M(lati);

    std::vector<std::vector<double>> LonandLat = {{118.7806222, 31.8351078},
                                                    {118.7810634, 31.8347478},
                                                    {118.7814473, 31.8344676},
                                                    {118.7814108, 31.8348904},
                                                    {118.7806701, 31.8351548},
                                                    {118.7815429, 31.8352704},
                                                    {118.7811188, 31.8351964},
                                                    {118.7820834, 31.8346503}
                                                };


    for (int i = 0; i < LonandLat.size(); i++) {
        
        double longti = LonandLat[i][0];
        double lati = LonandLat[i][1];

        x_intr = Lon2M(longti);
        y_intr = Lat2M(lati);

        std::cout << x_intr << " " << y_intr << std::endl;

        if ((x_intr < min_pt.x || x_intr > max_pt.x || y_intr < min_pt.y || y_intr > max_pt.y)) {
            std::cerr << "the data point is out of the range" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        test_point_ = surface.GetFittingPoint(x_intr,  y_intr);

        std::cout << std::fixed << std::setprecision(8) << "x: " << x_intr << " " << "y: " << y_intr << std::endl;
        std::cout << std::fixed << std::setprecision(8) << "latitude: " << lati << " " << "longitude: " << longti << " " << "height: " << test_point_.z + 21.48 <<  std::endl;
        test_points_set.push_back(Point(x_intr, y_intr, test_point_.z));
    }


    // if ((x_intr < min_pt_g.x || x_intr > max_pt_g.x || y_intr < min_pt_g.y || y_intr > max_pt_g.y)) {
    //     std::cerr << "the data point is out of the range" << std::endl;
    //     std::exit(EXIT_FAILURE);
    // }

    // Eigen::Vector3d home_point_ = surface.GetFittingPoint(home_x_, home_y_);
    // std::cout << std::fixed << std::setprecision(8) << "latitude: " << home_x_ << " " << "longitude: " << home_y_ << " " << "height: " << home_point_ <<  std::endl;
    
    // test_point_ = surface.GetFittingPoint(x_intr,  y_intr);




    // double lat_ = M2Lat(x_intr);
    // double lon_ = M2Lon(y_intr);

    // std::cout << std::fixed << std::setprecision(11) << "x: " << x_intr << " " << "y: " << y_intr << std::endl;
    // std::cout << std::fixed << std::setprecision(11) << "latitude: " << lati << " " << "longitude: " << longti << " " << "height: " << test_point_.z + 21.48 <<  std::endl;


    std::cout << "The low peak is: " << min_pt_g.z << std::endl;
    std::cout << "The high peak is: " << max_pt_g.z << std::endl;
    std::cout << "The height difference is: " << max_pt_g.z - min_pt_g.z << std::endl;


    // pcl::PointCloud<pcl::PointXYZ>::Ptr bsplline_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    // for (const auto& point : vertices)
    // {
    //     bsplline_cloud->points.push_back(pcl::PointXYZ(point.x, point.y, point.z));
    // }

    std::string out_bspline_folder = "/media/taole/HHD/Doc/daily_work/work_tg/Bspline/testEncapsualtion/bspline_result";
    // std::string output_bspline_file = out_bspline_folder + "/" + "simple_bspline_hill" + ".pcd";  
    // pcl::io::savePCDFileBinary(output_bspline_file, *bsplline_cloud);


    pcl::PointCloud<pcl::PointXYZ>::Ptr control_point_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    // 提取控制点
    std::vector<Point> control_points;
    for (const auto& row : cn_point)
    {
        for (const auto& point : row)
        {
            // control_point_cloud->points.push_back(pcl::PointXYZ(point.x, point.y, point.z));
            control_points.push_back(point);
        }
    }

    // std::string output_control_file = out_bspline_folder + "/" + "control_point" + ".pcd";  
    // pcl::io::savePCDFileBinary(output_control_file, *control_point_cloud);


    // std::string output_downsample_file = out_bspline_folder + "/" + "downsample" + ".pcd";  
    // pcl::io::savePCDFileBinary(output_downsample_file, *filtered_cloud);


    // std::string output_tree_top_file = out_bspline_folder + "/" + "tree_top" + ".pcd";  
    // pcl::io::savePCDFileBinary(output_tree_top_file, *tree_top);

    // std::string output_ground_file = out_bspline_folder + "/" + "ground" + ".pcd";  
    // pcl::io::savePCDFileBinary(output_ground_file, *ground_filtered);

    // pcl::PointCloud<pcl::PointXYZ>::Ptr test_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    // for (const auto& point : test_points_set)
    // {
    //     test_cloud->points.push_back(pcl::PointXYZ(point.x, point.y, point.z));
    // }

  
    // std::string output_test_file = out_bspline_folder + "/" + "test_point" + ".pcd";  
    // pcl::io::savePCDFileBinary(output_test_file, *test_cloud);  

    // std::vector<Point> knot_points;
    // knot_points = surface.GetKnotPoints();

    // // 可视化拟合点和控制点
    VisualizePointCloudV2(filtered_cloud, control_points);

    return 0;
}
