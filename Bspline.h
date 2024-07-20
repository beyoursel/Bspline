#ifndef BSPLINE_H
#define BSPLINE_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

using namespace std;


class bspSurface
{
public:
	bspSurface() {}
	bspSurface(const vector<vector<Eigen::Vector3d>>& cnPoint, 
				const vector<float>& knots_u, 
				vector<float>& knots_v);
	// constructor
	bspSurface(const bspSurface& surface);

	// 函数运算符重载
	bspSurface& operator=(const bspSurface& surface);

	// 根据参数u,v计算曲面上的坐标
	Eigen::Vector3d calPos(const float& u, const float& v);

	Eigen::Vector3d calPos(const vector<Eigen::Vector3d>& controlpoint, const vector<float>& knots, const float& t);


	void getFittingSurface(vector<Eigen::Vector3d>& vertices, float step);
	
	// obtain controlpoint TODO
	void getcontrolpoint();


private:
	int m_nu; // u向 0-nu
	int m_nv; // v向 0-nv
	int m_ku; // u向阶
	int m_kv; // v向阶
	vector<vector<Eigen::Vector3d>> m_cnPoint; //控制网格点坐标 （nu+1）x(nv+1)
	vector<float> m_knots_u; // u向节点向量 u_0, ..., u_(nu+ku)
	vector<float> m_knots_v; // v向节点向量 v_0, ..., v_(nv+kv)
};



#endif