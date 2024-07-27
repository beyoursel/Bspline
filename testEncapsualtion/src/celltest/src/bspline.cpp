#include "bspline.h"
#include <cmath>

namespace bspline {

BspSurface::BspSurface(const std::vector<std::vector<Point>>& cn_point, int k) 
{
    m_ku_ = k;
    m_kv_ = k;
    m_cn_point_ = cn_point;
    SetKnotVector(cn_point, m_knots_u_, m_knots_v_);

    // SetUniformKnotVector(cn_point, m_knots_u_, m_knots_v_);

    m_nu_ = m_cn_point_.size() - 1; 
    m_nv_ = m_cn_point_[0].size() - 1; 
    // m = n + p + 1 m+1: 节点个数; n+1：控制点个数; p: 阶数
}

// constructor
BspSurface::BspSurface(const BspSurface& surface) {
    m_cn_point_ = surface.m_cn_point_;
    m_knots_u_ = surface.m_knots_u_;
    m_knots_v_ = surface.m_knots_v_;

    m_nu_ = surface.m_nu_;
    m_nv_ = surface.m_nv_;
    m_ku_ = surface.m_kv_;
    m_kv_ = surface.m_kv_; 
}

// 函数运算符重载
BspSurface& BspSurface::operator=(const BspSurface& surface) {
    m_cn_point_ = surface.m_cn_point_;
    m_knots_u_ = surface.m_knots_u_;
    m_knots_v_ = surface.m_knots_v_;

    m_nu_ = surface.m_nu_;
    m_nv_ = surface.m_nv_;
    m_ku_ = surface.m_ku_;
    m_kv_ = surface.m_kv_;
    return *this;
}

// 根据参数u,v计算曲面上的坐标 u,v为网格坐标，不是实际点云三维坐标
Point BspSurface::CalPos(const double& u, const double& v) {

    std::vector<Point> v_constant(m_nu_ + 1); // 
    for (int i = 0; i < v_constant.size(); ++i)
    {
        v_constant[i] = CalPos(m_cn_point_[i], m_knots_v_, v); // v向Bspline插值
    }
    return CalPos(v_constant, m_knots_u_, u);
}


Point BspSurface::CalPos(const std::vector<Point>& controlpoint, const std::vector<double>& knots, const double& t)
{
    int n = controlpoint.size() - 1;
    int k = knots.size() - controlpoint.size(); // 阶数
    int L = 0;
    // 计算t所处的区间[t_L, t_(L+1)], t只在[knots[k-1], knots[n+1]]中有效 确定t在哪个区间
    if (t >= knots[n+1])
    {
        L = n;
    } else if (t <= knots[k-1])
    {
        L = k - 1;
    }
    else
    {
        for (int i = k - 1; i <= n + 1; ++i)
        {
            if (t >= knots[i] && t<knots[i+1])
            {
                L = i;
                break;
            }
        }
    }

    if (L >= n + 1) L = n;

    std::vector<Point> temp(k);
    for (int i = 0; i < k; ++i) {
        temp[i] = controlpoint[i + L - k + 1]; // 计算t处于哪个knot_span，进而选择control point
    }

    //de-BoorCox算法
    for (int r = 1; r <= k - 1; ++r)
    {
        for (int i = 0; i <= k - r - 1; ++i)
        {
            int index = L - k + 1 + r;
            double factor = 0;
            if (knots[index + i + k - r] != knots[index + i])
            {
                factor = (t - knots[index + i]) / (knots[index + i + k - r] - knots[index + i]); // 根据节点向量和坐标生成权重因子
            }
            temp[i] = factor*temp[i + 1] + (1 - factor)*temp[i]; // 根据权重因子和控制节点生成插值点

        }
    }

    return temp[0];

}


void BspSurface::GetFittingSurface(
    std::vector<Point>& vertices, double step)
{

    int m = static_cast<int>((m_knots_u_[m_nu_ + 1] - m_knots_u_[m_ku_ - 1]) / step);
    int n = static_cast<int>((m_knots_v_[m_nv_ + 1] - m_knots_v_[m_kv_ - 1]) / step);

    // std::cout << "m: " << m << " n: " << n << std::endl;

    for (int i = 0; i <= m; ++i)
    {
        for (int j = 0; j <= n; ++j)
        {
            float u = 0, v = 0;
            if (i == m)
            {
                u = m_knots_u_[m_nu_ + 1];
                v = m_knots_v_[m_kv_ - 1] + j*step;
                
            }
            else if (j == n)
            {
                u = m_knots_u_[m_ku_ - 1] + i*step;
                v = m_knots_v_[m_nv_ + 1];
                
            }
            else
            {
                u = m_knots_u_[m_ku_ - 1] + i*step;
                v = m_knots_v_[m_kv_ - 1] + j*step;
            }
            
            Point temp = CalPos(u, v);
            vertices.push_back(temp);
            // std::cout << temp(0) << " " << temp(1) << " " << temp(2) << std::endl;
        }
    }
}


void BspSurface::SetUniformKnotVector(std::vector<std::vector<Point>> cnPoint, std::vector<double>& knots_u_b, std::vector<double>& knots_v_b) {
    int num_control_points_u = cnPoint.size();
    int num_control_points_v = cnPoint[0].size();

    int num_knots_u = num_control_points_u + m_ku_;
    int num_knots_v = num_control_points_v + m_kv_;

    std::vector<double> knots_u(num_knots_u);
    std::vector<double> knots_v(num_knots_v);

    // // 均匀分布U方向节点
    for (int i = 0; i < num_knots_u; ++i) {
        knots_u[i] = static_cast<float>(i) / (num_knots_u - 1);
        // knots_u[i] = static_cast<float>(i);
    }

    // 均匀分布V方向节点
    for (int i = 0; i < num_knots_v; ++i) {
        knots_v[i] = static_cast<float>(i) / (num_knots_v - 1);
        // knots_v[i] = static_cast<float>(i);
    }

    // for (int i = 0; i <m_ku_; ++i) {
    //     knots_u[i] = 0.0;
    //     knots_u[knots_u.size() - 1 - i] = 1.0;
    // }
    // for (int i = m_ku_; i < knots_u.size() - m_ku_; ++i)
    // {
    //     knots_u[i] = i / double(num_knots_u - 1);
    // }

    // for (int i = 0; i <m_kv_; ++i) {
    //     knots_v[i] = 0.0;
    //     knots_v[knots_v.size() - 1 - i] = 1.0;
    // }
    // for (int i = m_kv_; i < knots_v.size() - m_kv_; ++i)
    // {
    //     knots_v[i] = i / double(num_knots_v - 1);
    // }


    knots_u_b = knots_u;
    knots_v_b = knots_v;

}


// 设置均匀节点向量
void BspSurface::SetKnotVector(std::vector<std::vector<Point>> cnPoint, std::vector<double>& knots_u_b, std::vector<double>& knots_v_b) {

    std::vector<double> knots_u(cnPoint.size() + m_ku_);
    std::vector<double> knots_v(cnPoint[0].size() + m_kv_);

    for (int i = 0; i < m_ku_; ++i)
    {
        knots_u[i] = 0.0f;
        knots_u[knots_u.size() - 1 - i] = static_cast<float>(cnPoint.size() - m_ku_ + 1);
    }
    for (int i = m_ku_; i < knots_u.size() - m_ku_; ++i)
    {
        knots_u[i] = static_cast<float>(i - m_ku_ + 1);
    }

    for (int i = 0; i < m_kv_; ++i)
    {
        knots_v[i] = 0.0f;
        knots_v[knots_v.size() - 1 - i] = static_cast<float>(cnPoint[0].size() - m_kv_ + 1);
    }
    for (int i = m_kv_; i < knots_v.size() - m_kv_; ++i)
    {
        knots_v[i] = static_cast<float>(i - m_kv_ + 1);
    }

    knots_u_b = knots_u;
    knots_v_b = knots_v;
}


std::vector<Point> BspSurface::GetKnotPoints() {

    std::vector<Point> ptc_est;

    for (int i = m_ku_ - 1; i <= m_nu_ + 1; i++) {
        for (int j = m_kv_ - 1; j <= m_nv_ + 1; j++) {
        float u = m_knots_u_[i];
        float v = m_knots_v_[j];
        Point temp = CalPos(u, v);
        ptc_est.push_back(temp);
        // std::cout << temp(0) << " " << temp(1) << " " << temp(2) << std::endl;
    }
    }
    return ptc_est;   

}


Point BspSurface::GetFittingPoint(double x, double y)
{   

    std::vector<double> knot_x;
    std::vector<double> knot_y;


    for (int j = m_kv_ - 1; j <= m_ku_ - 1; j++) {
        for (int i = m_ku_ - 1; i <= m_nu_ + 1; i++) 
        {
            double u = m_knots_u_[i];
            double v = m_knots_v_[j];
            Point temp = CalPos(u, v);
            knot_x.push_back(temp.x);
        }
    }

    std::cout << "knot_x[0]: " << knot_x[0] << " " << "knot_x[1]: " << knot_x[1] << std::endl;
    std::cout << "knot_x[-2]: " << knot_x[knot_x.size()-2] << " " << "knot_x[-1]: " << knot_x[knot_x.size()-1] << std::endl;


    for (int i = m_ku_ - 1; i <= m_kv_ - 1; i++) {
        for (int j = m_kv_ - 1; j <= m_nv_ + 1; j++) {
            double u = m_knots_u_[i];
            double v = m_knots_v_[j];
            Point temp = CalPos(u, v);
            knot_y.push_back(temp.y);
        }
    }

    std::cout << "knot_y[0]: " << knot_y[0] << " " << "knot_y[1]: " << knot_y[1] << std::endl;
    std::cout << "knot_y[-2]: " << knot_y[knot_y.size()-2] << " " << "knot_y[-1]: " << knot_y[knot_y.size()-1] << std::endl;

    int knot_grid_x;
    int knot_grid_y;

    for (int i = 0; i < knot_x.size() - 1; i++) {
        if (x >= knot_x[i] && x <knot_x[i+1]) {
            knot_grid_x = i;
        }
    }

    for (int i = 0; i < knot_y.size() - 1; i++) {
        if (y >= knot_y[i] && y <knot_y[i+1]) {
            knot_grid_y = i;
        }
    }


    float lu = (x - knot_x[knot_grid_x]) / (knot_x[knot_grid_x + 1] - knot_x[knot_grid_x]) * (m_knots_u_[m_ku_ - 1 + knot_grid_x + 1] - m_knots_u_[m_ku_ - 1 + knot_grid_x]); 
    float lv = (y - knot_y[knot_grid_y]) / (knot_y[knot_grid_y + 1] - knot_y[knot_grid_y]) * (m_knots_v_[m_kv_ - 1 + knot_grid_y + 1] - m_knots_v_[m_kv_ - 1 + knot_grid_y]);


    Point temp1 = CalPos(m_knots_u_[m_ku_ -1 + knot_grid_x] + lu, m_knots_v_[m_kv_ - 1 + knot_grid_y] + lv);
    // std::cout << temp1(0) << " " << temp1(1) << " " << temp1(2) << std::endl;


    return temp1;
}

}