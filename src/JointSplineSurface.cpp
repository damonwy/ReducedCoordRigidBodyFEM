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
	m_cs.setZero();
}

void JointSplineSurface::addControlFrame(int i, int j, Vector6d C) {
	for (int k = 0; k < 6; k++) {
		m_cs(i, j, k) = C(k);
	}
}

Eigen::Matrix4d JointSplineSurface::evalQ(Eigen::Vector2d q)const {
	// Evaluates spline frame
	Matrix4d Q;
	Q.setIdentity();
	for (int i = 0; i < 6; ++i) {
		Matrix4d Ci;	// 4 x 4 control point values
		for (int ii = 0; ii < 4; ++ii) {
			for (int jj = 0; jj < 4; ++jj) {
				Ci(ii, jj) = m_cs(ii, jj, i); 
			}
		}

		Vector6d ei = m_E.col(i);
		double phi = JointSplineSurface::Cfun(Ci, q);
		Vector6d eiphi = ei * phi;
		Q = Q * SE3::exp(eiphi);
	}
	return Q;
}


void JointSplineSurface::evalS(Eigen::Vector2d q, MatrixXd &S, Tensor6x2x2d &dSdq) {
	// Evaluates spline frame derivatives

	S.setZero();
	//dSdq.resize(6, 2, 2);
	dSdq.setZero();
	Vector6d e1 = m_E.col(0);
	//Eigen::array<int, 3> offsets = { 0, 0, 0 };
	//Eigen::array<int, 3> extents = { 4, 4, 1 };
	//Eigen::Tensor<double, 2> C1 = m_cs.slice(offsets, extents);
	Matrix4d C0;	
	for (int ii = 0; ii < 4; ++ii) {
		for (int jj = 0; jj < 4; ++jj) {
			C0(ii, jj) = m_cs(ii, jj, 0);
		}
	}


	//Eigen::Matrix4d C1_mat = Tensor_to_Matrix(C1, 4, 4);
	for (int i = 0; i < 2; i++) {
		double dphi1i = JointSplineSurface::dCfun(C0, i, q);
		S.col(i) = e1 * dphi1i;
		for (int j = 0; j < 2; j++) {
			double d2phiqij = JointSplineSurface::d2Cfun(C0, i, j, q);
			Vector6d res = e1 * d2phiqij;
			for (int t = 0; t < 6; t++) {
				dSdq(t, i, j) = res(t);
			}
		}
	}

	for (int k = 1; k < 6; k++) {
		Vector6d ek = m_E.col(k);
	//	Eigen::array<int, 3> offsets = { 0, 0, k }; // 4x4 control point values
	//	Eigen::array<int, 3> extents = { 4, 4, 1 };
	//	Eigen::Tensor<double, 2> Ck = m_cs.slice(offsets, extents);
	//	Eigen::Matrix4d Ck_mat = Tensor_to_Matrix(Ck, 4, 4);
		Matrix4d Ck;
		for (int ii = 0; ii < 4; ++ii) {
			for (int jj = 0; jj < 4; ++jj) {
				Ck(ii, jj) = m_cs(ii, jj, k);
			}
		}

		double phik = JointSplineSurface::Cfun(Ck, q);
		Vector6d ekphik = ek * phik;
		Matrix6d Ad = SE3::adjoint(SE3::inverse(SE3::exp(ekphik)));
		for (int i = 0; i < 2; ++i) {
			double dphiki = JointSplineSurface::dCfun(Ck, i, q);
			Matrix6d ad = SE3::ad(S.col(i));
			S.col(i) = ek * dphiki + Ad * S.col(i);
			for (int j = 0; j < 2; ++j) {
				double d2phikij = JointSplineSurface::d2Cfun(Ck, i, j, q);
				double dphikj = JointSplineSurface::dCfun(Ck, j, q);
				Vector6d temp;
				for (int t = 0; t < 6; t++) {
					temp(t) = dSdq(t, i, j);
				}

				Vector6d res = ek * d2phikij + Ad * (temp + ad * ek * dphikj);

				for (int t = 0; t < 6; t++) {
					dSdq(t, i, j) = res(t);
				}
			}
		}
	}
}

void JointSplineSurface::updateSelf() {
	m_Q = evalQ(m_q);
	Tensor6x2x2d dSdq;
	evalS(m_q, m_S, dSdq);

	Matrix6x2d dSdq0;	
	for (int ii = 0; ii < 6; ++ii) {
		for (int jj = 0; jj < 2; ++jj) {
			dSdq0(ii, jj) = dSdq(ii, jj, 0);
		}
	}

	Matrix6x2d dSdq1;
	for (int ii = 0; ii < 6; ++ii) {
		for (int jj = 0; jj < 2; ++jj) {
			dSdq1(ii, jj) = dSdq(ii, jj, 1);
		}
	}

	m_Sdot = dSdq0 * m_qdot(0) + dSdq1 * m_qdot(1);

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

void JointSplineSurface::drawSelf(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {

	Matrix4d E_wp;
	if (getParent() == nullptr) {
		E_wp.setIdentity();
	}
	else {
		E_wp = getParent()->E_wj;
	}
	Matrix4d E_wj0 = E_wp * E_pj0;

	// Draw spline frames
	progSimple->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	int qmax = 1;
	Vector2d qvec;

	VectorXd q0 = VectorXd::LinSpaced(qmax * 5, 0, qmax);
	VectorXd q1 = VectorXd::LinSpaced(qmax * 5, 0, qmax);
	for (int i = 0; i < q0.size(); ++i) {
		for (int j = 0; j < q1.size(); ++j) {
			qvec(0) = q0(i);
			qvec(1) = q1(j);

			Matrix4d Q = evalQ(qvec);
			MV->pushMatrix();
			Matrix4d E = E_wj0 * Q;
			MV->multMatrix(eigen_to_glm(E));
			glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
			glLineWidth(3);
			glBegin(GL_LINES);
			// X axis
			glColor3f(1.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(1.0, 0.0, 0.0);

			// Y axis
			glColor3f(0.0, 1.0, 0.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 1.0, 0.0);

			// Z axis
			glColor3f(0.0, 0.0, 1.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, 1.0);

			glEnd();
			MV->popMatrix();
		}	
	}

	// Draw current frame
	MV->pushMatrix();
	Matrix4d E = E_wj;
	MV->multMatrix(eigen_to_glm(E));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	glLineWidth(3);
	glBegin(GL_LINES);
	// X axis
	glColor3f(1.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(3.0, 0.0, 0.0);

	// Y axis
	glColor3f(0.0, 1.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 3.0, 0.0);

	// Z axis
	glColor3f(0.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 3.0);

	glEnd();
	MV->popMatrix();

	progSimple->unbind();
}


JointSplineSurface:: ~JointSplineSurface() {

}