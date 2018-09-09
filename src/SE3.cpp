#include "SE3.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

Matrix4d SE3::inverse(const Matrix4d &E)
{
	Matrix4d Einv = Matrix4d::Identity();
	Matrix3d R = E.block<3, 3>(0, 0);
	Vector3d p = E.block<3, 1>(0, 3);
	Matrix3d Rt = R.transpose();
	Einv.block<3, 3>(0, 0) = Rt;
	Einv.block<3, 1>(0, 3) = -Rt * p;
	return Einv;
}

Matrix3x6d SE3::gamma(const Eigen::Vector3d &r)
{
	// Gets the 3x6 Gamma matrix, for computing the point velocity
	Matrix3x6d G = Matrix3x6d::Zero();
	G.block<3, 3>(0, 0) = bracket3(r).transpose();
	G.block<3, 3>(0, 3) = Matrix3d::Identity();
	return G;
}

Matrix6d SE3::adjoint(const Matrix4d &E)
{
	// Gets the adjoint transform
	Matrix6d Ad = Matrix6d::Zero();
	Matrix3d R = E.block<3, 3>(0, 0);
	Vector3d p = E.block<3, 1>(0, 3);
	Ad.block(0, 0, 3, 3) = R;
	Ad.block(3, 0, 3, 3) = bracket3(p) * R;
	Ad.block(3, 3, 3, 3) = R;
	return Ad;
}

Matrix6d SE3::dAddt(const Matrix4d &E, const VectorXd &phi) 
{
	// Gets the time derivative of the adjoint
	Matrix6d dA;
	dA.setZero();
	Matrix3d R = E.block<3, 3>(0, 0);
	Vector3d p = E.block<3, 1>(0, 3);
	Vector3d w = phi.segment<3>(0);
	Vector3d v = phi.segment<3>(3);
	Matrix3d wbrac = bracket3(w);
	Matrix3d vbrac = bracket3(v);
	Matrix3d pbrac = bracket3(p);
	Matrix3d Rwbrac = R * wbrac;
	dA.block<3, 3>(0, 0) = Rwbrac;
	dA.block<3, 3>(3, 3) = Rwbrac;
	dA.block<3, 3>(3, 0) = R * vbrac + pbrac * Rwbrac;

	return dA;
}

Matrix3d SE3::bracket3(const Vector3d &a)
{
	// Gets S = [x], the skew symmetric matrix
	Matrix3d A = Matrix3d::Zero();
	A(0, 1) = -a(2);
	A(0, 2) = a(1);
	A(1, 0) = a(2);
	A(1, 2) = -a(0);
	A(2, 0) = -a(1);
	A(2, 1) = a(0);
	return A;
}

Matrix4d SE3::bracket6(const Vector6d &a)
{
	Matrix4d A = Matrix4d::Zero();
	A.block<3, 3>(0, 0) = bracket3(a.segment<3>(0));
	A.block<3, 1>(0, 3) = a.segment<3>(3);
	return A;
}

Vector3d SE3::unbracket3(const Matrix3d &A)
{
	// Gets ]x[ = S, the vector corresponding to a skew symmetric matrix
	Vector3d a;
	a(0) = A(2, 1);
	a(1) = A(0, 2);
	a(2) = A(1, 0);
	return a;
}

Vector6d SE3::unbracket6(const Matrix4d &A)
{
	Vector6d a;
	a.segment<3>(0) = unbracket3(A.block<3, 3>(0, 0));
	a(3) = A(0, 3);
	a(4) = A(1, 3);
	a(5) = A(2, 3);
	return a;
}

Matrix4d SE3::integrate(const Matrix4d &E0, const VectorXd &phi, double h)
{
	Matrix3d I = Matrix3d::Identity();
	Vector3d w = phi.segment<3>(0);
	Vector3d v = phi.segment<3>(3);
	Matrix4d phib = Matrix4d::Identity();
	phib.block<3, 1>(0, 3) = h*v;
	double wlen = w.norm();
	if (wlen > 1e-10) {
		w /= wlen;
		v /= wlen;
		// Rodrigues formula %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		double wX = w(0);
		double wY = w(1);
		double wZ = w(2);
		double c = cos(wlen * h);
		double s = sin(wlen * h);
		double c1 = 1.0 - c;
		Matrix3d R;
		R << c + wX * wX * c1, -wZ * s + wX * wY * c1, wY * s + wX * wZ * c1,
			wZ * s + wX * wY * c1, c + wY * wY * c1, -wX * s + wY * wZ * c1,
			-wY * s + wX * wZ * c1, wX * s + wY * wZ * c1, c + wZ * wZ * c1;
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		Matrix3d A = I - R;
		Vector3d cc = w.cross(v);
		Vector3d d = A * cc;
		double wv = w.dot(v);
		Vector3d p = (wv * wlen * h) * w + d;
		phib.block<3, 3>(0, 0) = R;
		phib.block<3, 1>(0, 3) = p;
		//cout << phib << endl;
	}
	return E0 * phib;
}

Vector6d SE3::inertiaCuboid(Eigen::Vector3d whd, double density) {
	// Gets the diagonal inertia of a cuboid with (width, height, depth)
	Vector6d m;
	m.setZero();
	double volume = whd(0) * whd(1) * whd(2);
	double mass = density * volume;
	m(0) = (1.0 / 12.0) * mass * (whd(1) * whd(1) + whd(2) * whd(2));
	m(1) = (1.0 / 12.0) * mass * (whd(2) * whd(2) + whd(0) * whd(0));
	m(2) = (1.0 / 12.0) * mass * (whd(0) * whd(0) + whd(1) * whd(1));
	m(3) = mass;
	m(4) = mass;
	m(5) = mass;
	return m;
}
