#include "bspline.h"

namespace bspline {

BspSurface::BspSurface(const std::vector<std::vector<Eigen::Vector3d>>& cn_point, int k) 
{
    m_ku_ = k;
    m_kv_ = k;
    m_cn_point_ = cn_point;
    SetKnotVector(cn_point, m_knots_u_, m_knots_v_);

    m_nu_ = m_cn_point_.size() - 1; 
    m_nv_ = m_cn_point_[0].size() - 1; 
    // m_ku = m_knots_u.size() - 1 - m_nu; // m = n + p + 1 m+1: 节点个数; n+1：控制点个数; p: 阶数
    // m_kv = m_knots_v.size() - 1 - m_nv;
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
Eigen::Vector3d BspSurface::CalPos(const float& u, const float& v) {

    std::vector<Eigen::Vector3d> v_constant(m_nu_ + 1); // 创建v_knot_vector
    for (int i = 0; i < v_constant.size(); ++i)
    {
        v_constant[i] = CalPos(m_cn_point_[i], m_knots_v_, v); // 更新每个节点对应的knot_vector
    }
    return CalPos(v_constant, m_knots_u_, u);
}


Eigen::Vector3d BspSurface::CalPos(const std::vector<Eigen::Vector3d>& controlpoint, const std::vector<float>& knots, const float& t)
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

    std::vector<Eigen::Vector3d> temp(k);
    for (int i = 0; i < k; ++i) {
        temp[i] = controlpoint[i + L - k + 1]; // 根据阶数确定使用的controlpoint个数
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
    std::vector<Eigen::Vector3d>& vertices, float step)
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
            
            Eigen::Vector3d temp = CalPos(u, v);
            vertices.push_back(temp);
        }
    }
}


// 设置均匀节点向量
void BspSurface::SetKnotVector(std::vector<std::vector<Eigen::Vector3d>> cnPoint, std::vector<float>& knots_u_b, std::vector<float>& knots_v_b) {

    std::vector<float> knots_u(cnPoint.size() + m_ku_);
    std::vector<float> knots_v(cnPoint[0].size() + m_kv_);

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



Eigen::Vector3d BspSurface::GetFittingPoint(float x, float y, float grid_size)
{   // TODO
    float x_grid = x / grid_size;
    float y_grid = y / grid_size;
    // std::cout << x_grid << " " << y_grid << std::endl;
    Eigen::Vector3d ptc_est = CalPos(x_grid, y_grid);
    return ptc_est;
}


}