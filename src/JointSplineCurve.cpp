#include "JointSplineCurve.h"

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

JointSplineCurve::JointSplineCurve()
{

}

JointSplineCurve::JointSplineCurve(shared_ptr<Body> body, shared_ptr<Joint> parent):
Joint(body, 1, parent)
{


}

JointSplineCurve:: ~JointSplineCurve() {

}

void JointSplineCurve::addControlFrame(Matrix4d C) {
	m_Cs.push_back(C);
	int ncfs = m_Cs.size();
	if (ncfs >= 2) {
		Matrix4d C0 = m_Cs[ncfs - 2];
		Matrix4d C1 = C;
		VectorXd x = C0.ldlt().solve(C1);
		m_dCs.push_back(SE3::unbracket6(SE3::log(x)));
	}

	// Cyclic
	Matrix4d C0 = C;
	Matrix4d C1 = m_Cs[0];
	VectorXd x = C0.ldlt().solve(C1);
	m_dCs[0] = SE3::unbracket6(SE3::log(x));

}

void JointSplineCurve::updateSelf() {
	m_Q = evalQ(m_q(0));
	auto dSdq = evalS(m_q(0));
	m_Sdot = dSdq * m_qdot(0);

}

double JointSplineCurve::Bsum(int i, double q) {
	// Evaluates Btilde
	Vector4d qvec;
	Vector4d b;
	qvec << 1, q, q * q, q * q * q;
	if (i == 1) {
		b = m_B.row(0) + m_B.row(1) + m_B.row(2) + m_B.row(3);
	}
	else if (i == 2) {
		b = m_B.row(1) + m_B.row(2) + m_B.row(3);
	}
	else if (i == 3) {
		b = m_B.row(2) + m_B.row(3);
	}
	else {
		b = m_B.row(3);
	}
	double f = b.transpose() * qvec;
	return f;
}

double JointSplineCurve::dBsum(int i, double q) {
	// Evaluates dBtilde / dq
	Vector3d qvec;
	Vector3d b;
	qvec << 1, 2 * q, 3 * q * q;
	if (i == 1) {
		b = m_B1.row(0) + m_B1.row(1) + m_B1.row(2) + m_B1.row(3);
	}
	else if (i == 2) {
		b = m_B1.row(1) + m_B1.row(2) + m_B1.row(3);
	}
	else if (i == 3) {
		b = m_B1.row(2) + m_B1.row(3);
	}
	else {
		b = m_B1.row(3);
	}
	double df = b.transpose() * qvec;
	return df;
}

double JointSplineCurve::d2Bsum(int i, double q) {
	// Evaluates d ^ 2Btilde / dq ^ 2
	Vector2d qvec;
	qvec << 2, 6 * q;
	Vector2d b;
	if (i == 1) {
		b = m_B2.row(0) + m_B2.row(1) + m_B2.row(2) + m_B2.row(3);
	}
	else if (i == 2) {
		b = m_B2.row(1) + m_B2.row(2) + m_B2.row(3);
	}
	else if (i == 3) {
		b = m_B2.row(2) + m_B2.row(3);
	}
	else {
		b = m_B2.row(3);
	}

	double d2f = b.transpose() * qvec;
	return d2f;
}

Matrix4d JointSplineCurve::evalQ(double q) {
	Matrix4d Q;

	// Evaluate spline frame
	int ncfs = m_Cs.size();
	// Wrap around
	int qmax = ncfs;
	if (q < 0) {
		q += qmax;
	}
	else if (q >= qmax) {
		q -= qmax;
	}

	int k = floor(q); // starting control frame (0-index)
	if (k >= ncfs) {
		k -= 1; // overflow
	}

	double q_ = q - k; // local q in [0, 1]

	Q = m_Cs[k];

	for (int i = 1; i < 3; ++i) {
		int ki = k + i;
		if (ki > ncfs) {
			ki -= ncfs; // Wrap for cyclic
		}
		double Bsum = JointSplineCurve::Bsum(i, q_);
		Matrix4d dC = SE3::bracket6(m_dCs[ki]);
		Q = Q * SE3::exp(dC * Bsum);
	}

	return Q;
}