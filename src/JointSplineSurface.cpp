#include "JointSplineSurface.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

#include "Body.h"
#include "MatrixStack.h"
#include "Program.h"
#include "SE3.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

const Matrix4d JointSplineSurface::m_B = getB();
const Matrix6d JointSplineSurface::m_E = getE();

JointSplineSurface::JointSplineSurface()
{
	
}

JointSplineSurface::JointSplineSurface(shared_ptr<Body> body, shared_ptr<Joint> parent) :
	Joint(body, 2, parent)
{

	m_cs.resize(4, 4, 6);
	m_cs.setZero();
}

void JointSplineSurface::addControlFrame(int i, int j, Vector6d C) {
	m_cs.chip(i, j) = C;
	

}

double JointSplineSurface::Cfun(Eigen::Matrix4d C, Eigen::Vector2d q) {
	double q0 = q(0);
	double q1 = q(1);

	Vector4d q0vec, q1vec;
	q0vec << 1, q0, q0*q0, q0*q0*q0;
	q1vec << 1, q1, q1*q1, q1*q1*q1;
	double f = (q1vec.transpose() * m_B.transpose()) * C * (m_B * q0vec);
	return f;

}

double JointSplineSurface::dCfun(Eigen::Matrix4d C, int i, Eigen::Vector2d q) {
	double q0 = q(0);
	double q1 = q(1);
	Vector4d q0vec, q1vec;
	if (i == 0) {
		q0vec << 0, 1, 2*q0, 3*q0*q0;
		q1vec << 1, q1, q1*q1, q1*q1*q1;
	}
	else {
		q0vec << 1, q0, q0*q0, q0*q0*q0;
		q1vec << 0, 1, 2*q1, 3*q1*q1;
	}

	double df = (q1vec.transpose() * m_B.transpose()) * C * (m_B * q0vec);
	return df;
}

double JointSplineSurface::d2Cfun(Eigen::Matrix4d C, int i, int j, Eigen::Vector2d q) {
	double q0 = q(0);
	double q1 = q(1);
	Vector4d q0vec, q1vec;

	if (i == 0 && j == 0) {
		q0vec << 0, 0, 2, 6 * q0;
		q1vec << 1, q1, q1 * q1, q1 * q1 * q1;
	}
	else if (i == 0 && j == 1) {
		q0vec << 0, 1, 2 * q0, 3 * q0 * q0;
		q1vec << 0, 1, 2 * q1, 3 * q1 * q1;
	}
	else if (i == 1 && j == 0) {
		q0vec << 0, 1, 2 * q0, 3 * q0 * q0;
		q1vec << 0, 1, 2 * q1, 3 * q1 * q1;
	}
	else if (i == 1 && j == 1) {
		q0vec << 1, q0, q0 * q0, q0 * q0 * q0;
		q1vec << 0, 0, 2, 6 * q1;
	}
	double d2f = (q1vec.transpose() * m_B.transpose()) * C * (m_B * q0vec);
	return d2f;
}


JointSplineSurface:: ~JointSplineSurface() {

}