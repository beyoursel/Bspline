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

    // for vis/debug
    // for (int i = m_ku_-1; i < m_knots_u_.size() - m_ku_ + 1; i++) 
    // {
    //     for (int j = m_kv_ - 1; j < m_knots_v_.size() - m_kv_ + 1; j++) {

    //         int u_x = m_knots_u_[i];
    //         int v_y = m_knots_v_[j];

    //         Point temp_test = CalPos(u_x, v_y);
    //         knot_cn_point_.push_back(temp_test);
    //     }

    // }

    
    GetActualKnotSpan(knot_x_, knot_y_);

}


int BspSurface::KnotId(const int& n, const int& k, const std::vector<double>& knots, const double& t) {

    int L = 0;
    if (t >= knots[n+1])
    {
        L = n;
    } else if (t <= knots[k - 1])
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
//     int k = knots.size() - controlpoint.size();
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


// version 2
Point BspSurface::CalPos(const double& u, const double& v) {

    // 判断u所在的knot_span  
    int u_id = KnotId(m_nu_, m_ku_, m_knots_u_, u);
    int v_id = KnotId(m_nv_, m_kv_, m_knots_v_, v);

    std::vector<Point> v_constant(m_nu_ + 1);

    for (int i = u_id - m_ku_ + 1; i <= u_id; i++)
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

    int knots_u_num = cnPoint.size() + m_ku_; // 控制点数(n+1)+阶数k(k'+1)
    int knots_v_num = cnPoint[0].size() + m_kv_;
    knots_u_b.clear();
    knots_v_b.clear();
    knots_u_b.resize(knots_u_num);
    knots_v_b.resize(knots_v_num);

    for (int i = 0; i < m_ku_; ++i) {
    knots_u_b[i] = 0.0f;
    knots_u_b[knots_u_num - 1 - i] = static_cast<float>(cnPoint.size() - m_ku_ + 1);
    }

    for (int i = m_ku_; i < knots_u_num - m_ku_; ++i) {
    knots_u_b[i] = static_cast<float>(i - m_ku_ + 1);
    }

    for (int i = 0; i < m_kv_; ++i) {
    knots_v_b[i] = 0.0f;
    knots_v_b[knots_v_num - 1 - i] = static_cast<float>(cnPoint[0].size() - m_kv_ + 1);
    }
    for (int i = m_kv_; i < knots_v_num - m_kv_; ++i) {
    knots_v_b[i] = static_cast<float>(i - m_kv_ + 1);
    }

    std::cout << "knots_u: ";
    for (auto& u: knots_u_b) {
        std::cout << u;
        std::cout << " ";
    }
    std::cout << std::endl;
    std::cout << "knots_v: " << std::endl;
    for (auto& v: knots_v_b) {
        std::cout << v;
        std::cout << " ";
    }
    std::cout << std::endl;
}


void BspSurface::GetActualKnotSpan(std::vector<double>& knot_x, std::vector<double>& knot_y) {

    for (int j = m_kv_ - 1; j <= m_kv_ - 1; j++) {
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

void BspSurface::GetUV(double sx, double sy, double& su, double& sv) {


    for (int i = 0; i < knot_x_.size() - 1; i++) {
        if (sx >= knot_x_[i] && sx <knot_x_[i+1]) {
            su = i;
        }
    }

    for (int i = 0; i < knot_y_.size() - 1; i++) {
        if (sy >= knot_y_[i] && sy <knot_y_[i+1]) {
            sv = i;
        }
    }

}

double BspSurface::SquareError(double sx, double sy, double fx, double fy) {
    return std::sqrt((fx - sx) * (fx - sx) + (fy - sy) * (fy - sy));
}


double BspSurface::GetFittingPoint(double sx, double sy)
{   
    double su, sv;    

    std::cout << "sx: " << sx << std::endl;    
    std::cout << "sy: " << sy << std::endl;   

    GetUV(sx, sy, su, sv);

    double tolerance = 0.1;
    double epsilon = 0.001;
    double u_min = su;
    double v_min = sv;
    double u_max = su + 1;
    double v_max = sv + 1;
    int count = 1;

    // 首先根据与knot之间的距离比例计算u, v。若计算误差小于设定值，则直接返回计算的高度，否则得使用二分法逼近
    double knot_range_u = m_knots_u_[m_ku_ - 1 + su + 1] - m_knots_u_[m_ku_ - 1 + su]; 
    double knot_range_v = m_knots_v_[m_kv_ - 1 + sv + 1] - m_knots_v_[m_kv_ - 1 + sv];

    float lu = (sx - knot_x_[su]) / (knot_x_[su + 1] - knot_x_[su]) * knot_range_u;
    float lv = (sy - knot_y_[sv]) / (knot_y_[sv + 1] - knot_y_[sv]) * knot_range_v;

    Point temp0 = CalPos(m_knots_u_[m_ku_ -1 + su] + lu, m_knots_v_[m_kv_ - 1 + sv] + lv);

    if (SquareError(sx, sy, temp0.x, temp0.y) < tolerance) {
        std::cout << "x_fit: " << temp0.x << std::endl;
        std::cout << "y_fit: " << temp0.y << std::endl;
        std::cout << "u: " << m_knots_u_[m_ku_ -1 + su] + lu << std::endl;
        std::cout << "v: " << m_knots_v_[m_kv_ - 1 + sv] + lv << std::endl;
        return temp0.z;
    } 


    while ((u_max - u_min > epsilon) && (v_max - v_min > epsilon)) {
        double u_mid = (u_min + u_max) / 2.0;
        double v_mid = (v_min + v_max) / 2.0;
        
        double x_fit, y_fit;

        Point temp = CalPos(u_mid, v_mid);
        x_fit = temp.x;
        y_fit = temp.y;
        
        double error = SquareError(sx, sy, x_fit, y_fit);
        
        if (error < tolerance) {
            std::cout << "x_fit: " << x_fit << std::endl;
            std::cout << "y_fit: " << y_fit << std::endl;
            std::cout << "u: " << u_mid << std::endl;
            std::cout << "v: " << v_mid << std::endl;
            std::cout << "iterations: " << count << std::endl;
            return temp.z;
        }
        count++;
        // 更新区间
        if (x_fit < sx) {
            u_min = u_mid;
        } else {
            u_max = u_mid;
        }
        
        if (y_fit < sy) {
            v_min = v_mid;
        } else {
            v_max = v_mid;
        }
    }

    return 0;
}
