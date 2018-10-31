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

	//m_cs.resize(4, 4, 6);
	//m_cs.setZero();
}

void JointSplineSurface::addControlFrame(int i, int j, Vector6d C) {
	//m_cs.chip(i, j) = C;
	

}

void JointSplineSurface::evalS(Eigen::Vector2d q, Matrix6x2d &S, Eigen::Tensor<double, 3> &dSdq) {
	// Evaluates spline frame derivatives
	//S.setZero();
	//dSdq.resize(6, 2, 2);
	//dSdq.setZero();
	//Vector6d e1 = m_E.col(0);
	//Eigen::array<int, 3> offsets = { 0, 0, 0 };
	//Eigen::array<int, 3> extents = { 4, 4, 1 };
	//Eigen::Tensor<double, 2> C1 = m_cs.slice(offsets, extents);
	//Eigen::Matrix4d C1_mat = Tensor_to_Matrix(C1, 4, 4);
	//for (int i = 0; i < 2; i++) {
	//	double dphi1i = JointSplineSurface::dCfun(C1_mat, i, q);
	//	S.col(i) = e1 * dphi1i;
	//	for (int j = 0; j < 2; j++) {
	//		double d2phiqij = JointSplineSurface::d2Cfun(C1_mat, i, j, q);
	//		Vector6d res = e1 * d2phiqij;
	//		for (int t = 0; t < 6; t++) {
	//			dSdq(t, i, j) = res(t);
	//		}
	//	}
	//}

	//for (int k = 1; k < 6; k++) {
	//	Vector6d ek = m_E.col(k);
	//	Eigen::array<int, 3> offsets = { 0, 0, k }; // 4x4 control point values
	//	Eigen::array<int, 3> extents = { 4, 4, 1 };
	//	Eigen::Tensor<double, 2> Ck = m_cs.slice(offsets, extents);
	//	Eigen::Matrix4d Ck_mat = Tensor_to_Matrix(Ck, 4, 4);
	//	double phik = JointSplineSurface::Cfun(Ck_mat, q);
	//	Vector6d ekphik = ek * phik;
	//	Matrix6d Ad = SE3::adjoint(SE3::inverse(SE3::exp(ekphik)));
	//	for (int i = 0; i < 2; ++i) {
	//		double dphiki = JointSplineSurface::dCfun(Ck_mat, i, q);
	//		Matrix6d ad = SE3::adjoint(S.col(i));
	//		S.col(i) = ek * dphiki + Ad * S.col(i);
	//		for (int j = 0; j < 2; ++j) {
	//			double d2phikij = JointSplineSurface::d2Cfun(Ck_mat, i, j, q);
	//			double dphikj = JointSplineSurface::dCfun(Ck_mat, j, q);
	//			
	//		}
	//	}


	//}
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