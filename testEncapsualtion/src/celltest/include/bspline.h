#ifndef BSPLINE_H
#define BSPLINE_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>


namespace bspline {


class BspSurface
{
public:
	BspSurface() {}
	BspSurface(const std::vector<std::vector<Eigen::Vector3d>>& cnPoint, int k);
	// constructor
	BspSurface(const BspSurface& surface);

	// 函数运算符重载
	BspSurface& operator=(const BspSurface& surface);

	// 根据参数u,v计算曲面上的坐标
	Eigen::Vector3d CalPos(const float& u, const float& v);

	Eigen::Vector3d CalPos(const std::vector<Eigen::Vector3d>& controlpoint, const std::vector<float>& knots, const float& t);

	void SetKnotVector(std::vector<std::vector<Eigen::Vector3d>> cnPoint, std::vector<float>& knots_u_b, std::vector<float>& knots_v_b);

	void GetFittingSurface(std::vector<Eigen::Vector3d>& vertices, float step);
	
	// obtain interpolate point
	Eigen::Vector3d GetFittingPoint(float x, float y, float grid_size);


private:
	int m_nu_; // u向 0-nu
	int m_nv_; // v向 0-nv
	int m_ku_; // u向阶
	int m_kv_; // v向阶
	std::vector<std::vector<Eigen::Vector3d>> m_cn_point_; //控制网格点坐标 （nu+1）x(nv+1)
	std::vector<float> m_knots_u_; // u向节点向量 u_0, ..., u_(nu+ku)
	std::vector<float> m_knots_v_; // v向节点向量 v_0, ..., v_(nv+kv)
};

}

#endif