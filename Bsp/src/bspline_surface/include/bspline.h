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

    // 比较操作符，按照 x 坐标排序
    bool operator<(const Point& other) const {
        return x < other.x;
    }

    // 与 double 类型数据相乘的运算符重载
    Point operator*(double scalar) const {
        return Point(x * scalar, y * scalar, z * scalar);
    }

    // 赋值运算符，支持与 double 相乘后赋值给原对象
    Point& operator*=(double scalar) {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    // Point 之间的加法运算符重载
    Point operator+(const Point& other) const {
        return Point(x + other.x, y + other.y, z + other.z);
    }

    // 赋值运算符，支持 Point 之间的加法后赋值给原对象
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
	BspSurface(const std::vector<std::vector<Point>>& cnPoint, int k);
	// constructor
	BspSurface(const BspSurface& surface);

	// 函数运算符重载
	BspSurface& operator=(const BspSurface& surface);

	// 根据参数u,v计算曲面上的坐标
	Point CalPos(const double& u, const double& v);

	Point CalPos(const std::vector<Point>& controlpoint, const std::vector<double>& knots, const double& t);

	void SetKnotVector(std::vector<std::vector<Point>> cnPoint, std::vector<double>& knots_u_b, std::vector<double>& knots_v_b);

	void SetUniformKnotVector(std::vector<std::vector<Point>> cnPoint, std::vector<double>& knots_u_b, std::vector<double>& knots_v_b);
	void GetFittingSurface(std::vector<Point>& vertices, double step);
	
	std::vector<Point> GetKnotPoints();
	// obtain interpolate point
	Point GetFittingPoint(double x, double y);


private:
	int m_nu_; // u向 0-nu
	int m_nv_; // v向 0-nv
	int m_ku_; // u向阶
	int m_kv_; // v向阶
	std::vector<std::vector<Point>> m_cn_point_; //控制网格点坐标 （nu+1）x(nv+1)
	std::vector<double> m_knots_u_; // u向节点向量 u_0, ..., u_(nu+ku)
	std::vector<double> m_knots_v_; // v向节点向量 v_0, ..., v_(nv+kv)
};

#endif