#include "bspline.h"
#include <cmath>


BspSurface::BspSurface(const std::vector<std::vector<Point>>& cn_point, int k) 
{
    m_ku_ = k;
    m_kv_ = k;
    m_cn_point_ = cn_point;

    SetKnotVector(m_cn_point_, m_knots_u_, m_knots_v_);

    m_nu_ = m_cn_point_.size() - 1; 
    m_nv_ = m_cn_point_[0].size() - 1; 

    // for vis knots 注意删除
    for (int i = m_ku_; i <= m_nu_ + 1; i++) 
    {
        for (int j = m_kv_; j <= m_nv_ + 1; j++) {
            double u_x = m_knots_u_[i];
            double v_y = m_knots_v_[j];
            Point temp_test = CalPos(u_x, v_y);
            knot_cn_point_.push_back(temp_test);
        }

    }

    
    GetActualKnotSpan(knot_x_, knot_y_);

}


int BspSurface::KnotId(const int& n, const int& k, const std::vector<double>& knots, const double& t) {

    int L = 0;
    if (t >= knots[n+1])
    {
        L = n;
    } else if (t <= knots[k]) // k需要修改，如果multiple=3改称1
    {
        L = k;
    }
    else
    {
        for (int i = k; i <= n + 1; ++i)
        {
            if (t >= knots[i] && t<knots[i+1])
            {
                L = i;
                break;
            }
        }
    }

    if (L >= n + 1) L = n;

    return L;
}

// // old version
// Point BspSurface::CalPos(const double& u, const double& v) {

//     std::vector<Point> v_constant(m_nu_ + 1);
//     for (int i = 0; i < v_constant.size(); ++i)
//     {
//         v_constant[i] = CalPos(m_cn_point_[i], m_knots_v_, v);
//     }
//     return CalPos(v_constant, m_knots_u_, u);
// }

// // old version
// Point BspSurface::CalPos(const std::vector<Point>& controlpoint, const std::vector<double>& knots, const double& t)
// {
//     int n = controlpoint.size() - 1;
//     int k = knots.size() - controlpoint.size(); // k_阶需要k_+1控制点
//     int L = 0;

//     if (t >= knots[n+1])
//     {
//         L = n;
//     } else if (t <= knots[k-1])
//     {
//         L = k - 1;
//     }
//     else
//     {
//         for (int i = k - 1; i <= n + 1; ++i)
//         {
//             if (t >= knots[i] && t<knots[i+1])
//             {
//                 L = i;
//                 break;
//             }
//         }
//     }

//     if (L >= n + 1) L = n;

//     std::vector<Point> temp(k);
//     for (int i = 0; i < k; ++i) {
//         temp[i] = controlpoint[i + L - k + 1];
//     }

//     //de-BoorCox
//     for (int r = 1; r <= k - 1; ++r)
//     {
//         for (int i = 0; i <= k - r - 1; ++i)
//         {
//             int index = L - k + 1 + r;
//             double factor = 0;
//             if (knots[index + i + k - r] != knots[index + i])
//             {
//                 factor = (t - knots[index + i]) / (knots[index + i + k - r] - knots[index + i]);
//             }
//             temp[i] = factor*temp[i + 1] + (1 - factor)*temp[i];

//         }
//     }

//     return temp[0];

// }



Point BspSurface::CalPos(const double& u, const double& v) {

    // 判断u所在的knot_span  
    int u_id = KnotId(m_nu_, m_ku_, m_knots_u_, u);
    int v_id = KnotId(m_nv_, m_kv_, m_knots_v_, v);

    // std::cout << "u_id" << u_id << std::endl;
    // std::cout << "v_id" << v_id << std::endl;

    // 根据u_id进行提前终止
    std::vector<Point> v_constant(m_nu_ + 1);

    for (int i = u_id - m_ku_ ; i <= u_id; i++) // 9 10 11 12
    {   
        v_constant[i] = CalPos(m_cn_point_[i], m_knots_v_, v, v_id); 
    }

    return CalPos(v_constant, m_knots_u_, u, u_id);
}


Point BspSurface::CalPos(const std::vector<Point>& controlpoint, const std::vector<double>& knots, const double& t, int L)
{
    int n = controlpoint.size() - 1;
    int k = knots.size() - controlpoint.size(); // k_阶需要k_+1个控制点

    std::vector<Point> temp(k);
    for (int i = 0; i < k; ++i) {
        temp[i] = controlpoint[i + L - k + 1];
    }

    //de-BoorCox
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


void BspSurface::GetFittingSurface(
    std::vector<Point>& vertices, double step)
{

    int m = static_cast<int>((m_knots_u_[m_nu_ + 1] - m_knots_u_[m_ku_ - 1]) / step);
    int n = static_cast<int>((m_knots_v_[m_nv_ + 1] - m_knots_v_[m_kv_ - 1]) / step);

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
        }
    }
}


// k-multiplicity uniform knot vector
void BspSurface::SetKnotVector(std::vector<std::vector<Point>> cnPoint, std::vector<double>& knots_u_b, std::vector<double>& knots_v_b) {

    std::vector<double> knots_u(cnPoint.size() + m_ku_ + 1);
    std::vector<double> knots_v(cnPoint[0].size() + m_kv_ + 1);

    // for (int i = 0; i < m_ku_; ++i)
    // {
    //     knots_u[i] = 0.0f;
    //     knots_u[knots_u.size() - 1 - i] = static_cast<float>(cnPoint.size() - m_ku_ + 2);
    // }
    // for (int i = m_ku_; i < knots_u.size() - m_ku_; ++i)
    // {
    //     knots_u[i] = static_cast<float>(i - m_ku_ + 1);
    // }

    // for (int i = 0; i < m_kv_; ++i)
    // {
    //     knots_v[i] = 0.0f;
    //     knots_v[knots_v.size() - 1 - i] = static_cast<float>(cnPoint[0].size() - m_kv_ + 2);
    // }
    // for (int i = m_kv_; i < knots_v.size() - m_kv_; ++i)
    // {
    //     knots_v[i] = static_cast<float>(i - m_kv_ + 1);
    // }


    // multiple k+1
    for (int i = 0; i <= m_ku_; ++i)
    {
        knots_u[i] = 0.0f;
        knots_u[knots_u.size() - 1 - i] = static_cast<float>(cnPoint.size() - m_ku_);
    }
    for (int i = m_ku_ + 1; i < knots_u.size() - m_ku_ - 1; ++i)
    {
        knots_u[i] = static_cast<float>(i - m_ku_);
    }

    for (int i = 0; i <= m_kv_; ++i)
    {
        knots_v[i] = 0.0f;
        knots_v[knots_v.size() - 1 - i] = static_cast<float>(cnPoint[0].size() - m_kv_);
    }
    for (int i = m_kv_ + 1; i < knots_v.size() - m_kv_ - 1; ++i)
    {
        knots_v[i] = static_cast<float>(i - m_kv_);
    }

    knots_u_b = knots_u;
    knots_v_b = knots_v;

    std::cout << "knots_u: ";
    for (auto& u: knots_u) {
        std::cout << u;
        std::cout << " ";
    }
    std::cout << std::endl;
    std::cout << "knots_v: " << std::endl;
    for (auto& v: knots_v) {
        std::cout << v;
        std::cout << " ";
    }
    std::cout << std::endl;
}


void BspSurface::GetActualKnotSpan(std::vector<double>& knot_x, std::vector<double>& knot_y) {

    for (int j = m_kv_ - 1; j <= m_kv_ - 1; j++) { // little bug
            for (int i = m_ku_ - 1; i <= m_nu_ + 1; i++) 
            {
                double u = m_knots_u_[i];
                double v = m_knots_v_[j];
                Point temp = CalPos(u, v);
                knot_x.push_back(temp.x);
            }
        }

        for (int i = m_ku_ - 1; i <= m_ku_ - 1; i++) {
            for (int j = m_kv_ - 1; j <= m_nv_ + 1; j++) {
                double u = m_knots_u_[i];
                double v = m_knots_v_[j];
                Point temp = CalPos(u, v);
                knot_y.push_back(temp.y);
            }
        }

}


Point BspSurface::GetFittingPoint(double x, double y)
{   
    int knot_grid_x;
    int knot_grid_y;

    for (int i = 0; i < knot_x_.size() - 1; i++) {
        if (x >= knot_x_[i] && x <knot_x_[i+1]) {
            knot_grid_x = i;
        }
    }

    for (int i = 0; i < knot_y_.size() - 1; i++) {
        if (y >= knot_y_[i] && y <knot_y_[i+1]) {
            knot_grid_y = i;
        }
    }

    float lu = (x - knot_x_[knot_grid_x]) / (knot_x_[knot_grid_x + 1] - knot_x_[knot_grid_x]) * (m_knots_u_[m_ku_ - 1 + knot_grid_x + 1] - m_knots_u_[m_ku_ - 1 + knot_grid_x]); 
    float lv = (y - knot_y_[knot_grid_y]) / (knot_y_[knot_grid_y + 1] - knot_y_[knot_grid_y]) * (m_knots_v_[m_kv_ - 1 + knot_grid_y + 1] - m_knots_v_[m_kv_ - 1 + knot_grid_y]);

    Point temp1 = CalPos(m_knots_u_[m_ku_ -1 + knot_grid_x] + lu, m_knots_v_[m_kv_ - 1 + knot_grid_y] + lv);

    return temp1;
}
