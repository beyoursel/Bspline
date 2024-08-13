#ifndef BSPLINE_H
#define BSPLINE_H

#include <vector>
#include <iostream>


struct Point {

    Point() : x(0.0), y(0.0), z(0.0) {}
    Point(double x_p, double y_p, double z_p) : x(x_p), y(y_p), z(z_p) {}

    double x;
	double y;
	double z;

    bool operator<(const Point& other) const {
        return x < other.x;
    }

    Point operator*(double scalar) const {
        return Point(x * scalar, y * scalar, z * scalar);
    }

    Point& operator*=(double scalar) {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    Point operator+(const Point& other) const {
        return Point(x + other.x, y + other.y, z + other.z);
    }

    Point& operator+=(const Point& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

};

inline Point operator*(double scalar, const Point& point) {
    return Point(point.x * scalar, point.y * scalar, point.z * scalar);
}


class BspSurface
{
public:
    BspSurface() {}

/**
 *  \brief Construct a new BspSurface object 
 * 
 *  \param cnPoint 控制点
 *  \param k 阶数
 */
	BspSurface(const std::vector<std::vector<Point>>& cnPoint, int k);

/**
 *  \brief 计算(u, v)对应的高度
 * 
 *  \param u
 *  \param v
 *  \return Point
 */
	Point CalPos(const double& u, const double& v);

/**
 *  \brief 在指定维度下对t位置进行拟合
 *  \param controlpoint 控制点
 *  \param knots 节点向量
 *  \param t 节点变量值
 *  \param l 节点向量区间id
 *  \return Point
 */
	Point CalPos(const std::vector<Point>& controlpoint, const std::vector<double>& knots, const double& t, int l);

/**
 *  \brief 设置节点向量
 *  \param cnPoint
 *  \param knots_u_b 存储u向节点向量
 *  \param knots_v_b 存储v向节点向量
 */
	void SetKnotVector(std::vector<std::vector<Point>> cnPoint, std::vector<double>& knots_u_b, std::vector<double>& knots_v_b);

/**
 *  \brief 拟合整个曲面
 *  \param vertices 存储拟合的曲面点集
 *  \param step 步长
 */
	void GetFittingSurface(std::vector<Point>& vertices, double step);

/**
 * \brief 获得knot_vector在ENU中实际位置
 * \param knot_x 保存u向knot对应的ENU位置
 * \param knot_y 保存v向knot对应的ENU位置
 */
    void GetActualKnotSpan(std::vector<double>& knot_x, std::vector<double>& knot_y);

/**
 *  \brief 输入{(x, y) ENU坐标系}获得对应位置拟合的高度
 *  \param x
 *  \param y
 *  \return Point
 */
	Point GetFittingPoint(double x, double y);

/**
 *  \brief 获得节点向量所处的knot_span区间位置
 *  \param n 控制点数-1
 *  \param k 阶数
 *  \param t 节点向量
 *  \return int 节点向量区间的id
 */
    int KnotId(const int& n, const int& k, const std::vector<double>& knots, const double& t);


    std::vector<double> knot_x_; // knot-vector在ENU对应位置
    std::vector<double> knot_y_; // 同上

private:
	int m_nu_; // 控制点的行数-1
	int m_nv_; // 控制点列数-1
	int m_ku_; // u向阶数
	int m_kv_; // v向阶数
	std::vector<std::vector<Point>> m_cn_point_; // 控制点
	std::vector<double> m_knots_u_; // u向节点向量
	std::vector<double> m_knots_v_; // v向节点向量
};

#endif