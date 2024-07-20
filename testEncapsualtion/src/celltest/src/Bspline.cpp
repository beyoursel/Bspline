#include "Bspline.h"


bspSurface::bspSurface(const vector<vector<Eigen::Vector3d>>& cnPoint, int k) 
{
    m_ku = k;
    m_kv = k;
    m_cnPoint = cnPoint;
    setknotvector(cnPoint, m_knots_u, m_knots_v);

    m_nu = m_cnPoint.size() - 1; 
    m_nv = m_cnPoint[0].size() - 1; 
    // m_ku = m_knots_u.size() - 1 - m_nu; // m = n + p + 1 m+1: 节点个数; n+1：控制点个数; p: 阶数
    // m_kv = m_knots_v.size() - 1 - m_nv;
}

// constructor
bspSurface::bspSurface(const bspSurface& surface) {
    m_cnPoint = surface.m_cnPoint;
    m_knots_u = surface.m_knots_u;
    m_knots_v = surface.m_knots_v;

    m_nu = surface.m_nu;
    m_nv = surface.m_nv;
    m_ku = surface.m_kv;
    m_kv = surface.m_kv; 
}

// 函数运算符重载
bspSurface& bspSurface::operator=(const bspSurface& surface) {
    m_cnPoint = surface.m_cnPoint;
    m_knots_u = surface.m_knots_u;
    m_knots_v = surface.m_knots_v;

    m_nu = surface.m_nu;
    m_nv = surface.m_nv;
    m_ku = surface.m_ku;
    m_kv = surface.m_kv;
    return *this;
}

// 根据参数u,v计算曲面上的坐标
Eigen::Vector3d bspSurface::calPos(const float& u, const float& v) {

    vector<Eigen::Vector3d> v_constant(m_nu + 1);
    for (int i = 0; i < v_constant.size(); ++i)
    {
        v_constant[i] = calPos(m_cnPoint[i], m_knots_v, v);
    }
    return calPos(v_constant, m_knots_u, u);
}

Eigen::Vector3d bspSurface::calPos(const vector<Eigen::Vector3d>& controlpoint, const vector<float>& knots, const float& t)
{
    int n = controlpoint.size() - 1;
    int k = knots.size() - controlpoint.size(); // 阶数
    int L = 0;
    // 计算t所处的区间[t_L, t_(L+1)], t只在[knots[k-1], knots[n+1]]中有效
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

    vector<Eigen::Vector3d> temp(k);
    for (int i = 0; i < k; ++i) {
        temp[i] = controlpoint[i + L - k + 1];
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
                factor = (t - knots[index + i]) / (knots[index + i + k - r] - knots[index + i]);
            }
            temp[i] = factor*temp[i + 1] + (1 - factor)*temp[i];

        }
    }
    return temp[0];

}


void bspSurface::getFittingSurface(
    vector<Eigen::Vector3d>& vertices, float step)
{

    int m = static_cast<int>((m_knots_u[m_nu + 1] - m_knots_u[m_ku - 1]) / step);
    int n = static_cast<int>((m_knots_v[m_nv + 1] - m_knots_v[m_kv - 1]) / step);

    for (int i = 0; i <= m; ++i)
    {
        for (int j = 0; j <= n; ++j)
        {
            float u = 0, v = 0;
            if (i == m)
            {
                u = m_knots_u[m_nu + 1];
                v = m_knots_v[m_kv - 1] + j*step;
                
            }
            else if (j == n)
            {
                u = m_knots_u[m_ku - 1] + i*step;
                v = m_knots_v[m_nv + 1];
                
            }
            else
            {
                u = m_knots_u[m_ku - 1] + i*step;
                v = m_knots_v[m_kv - 1] + j*step;
            }
            
            Eigen::Vector3d temp = calPos(u, v);
            vertices.push_back(temp);
        }
    }
}


// 设置均匀节点向量
void bspSurface::setknotvector(vector<vector<Eigen::Vector3d>> cnPoint, vector<float>& knots_u_b, vector<float>& knots_v_b) {

    vector<float> knots_u(cnPoint.size() + m_ku);
    vector<float> knots_v(cnPoint[0].size() + m_kv);

    for (int i = 0; i < m_ku; ++i)
    {
        knots_u[i] = 0.0f;
        knots_u[knots_u.size() - 1 - i] = static_cast<float>(cnPoint.size() - m_ku + 1);
    }
    for (int i = m_ku; i < knots_u.size() - m_ku; ++i)
    {
        knots_u[i] = static_cast<float>(i - m_ku + 1);
    }

    for (int i = 0; i < m_kv; ++i)
    {
        knots_v[i] = 0.0f;
        knots_v[knots_v.size() - 1 - i] = static_cast<float>(cnPoint[0].size() - m_kv + 1);
    }
    for (int i = m_kv; i < knots_v.size() - m_kv; ++i)
    {
        knots_v[i] = static_cast<float>(i - m_kv + 1);
    }

    knots_u_b = knots_u;
    knots_v_b = knots_v;
}


//生成控制顶点的顶点及边上点的索引，用于绘制
void bspSurface::getcontrolpoint()
{   // TODO
    std::cout << "It's going TODO" << std::endl;
}



