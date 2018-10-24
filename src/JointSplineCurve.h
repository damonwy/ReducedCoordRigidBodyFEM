#pragma once
// JointSplineCurve 
//    Uses cubic B-splines.
//    Assumes a cyclic spline curve

#ifndef MUSCLEMASS_SRC_JOINTSPLINECURVE_H_
#define MUSCLEMASS_SRC_JOINTSPLINECURVE_H_

#include "Joint.h"

class Body;

class JointSplineCurve : public Joint {

public:
	JointSplineCurve();
	JointSplineCurve(std::shared_ptr<Body> body, std::shared_ptr<Joint> parent = nullptr);
	virtual ~JointSplineCurve();

	void addControlFrame(Eigen::Matrix4d C);
	void updateSelf();
	void drawSelf(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P) const;

	static double Bsum(int i, double q);
	static double dBsum(int i, double q);
	static double d2Bsum(int i, double q);

	static const Eigen::Matrix4d& getB() {
		static Eigen::Matrix4d _B(1.0 / 6.0 * (Eigen::Matrix4d() << 
			1.0, -3.0, 3.0, -1.0, 
			4.0, 0.0, -6.0, 3.0, 
			1.0, 3.0, 3.0, -3.0, 
			0.0, 0.0, 0.0, 1.0).finished());
		return _B;
	}

	static const Matrix4x3d& getB1() {
		static Matrix4x3d _B1(1.0 / 6.0 * (Matrix4x3d() <<
			-3.0, 3.0, -1.0,
			0.0, -6.0, 3.0,
			3.0, 3.0, -3.0,
			0.0, 0.0, 1.0).finished());
		return _B1;
	}

	static const Matrix4x2d& getB2() {
		static Matrix4x2d _B2(1.0 / 6.0 * (Matrix4x2d() <<
			3.0, -1.0,
			-6.0, 3.0,
			3.0, -3.0,
			0.0, 1.0).finished());
		return _B2;
	}

	static const Eigen::Matrix4d m_B; // Bspline coeffs
	static const Matrix4x3d m_B1;	 // Bspline coeffs, 1st column removed
	static const Matrix4x2d m_B2;	 // Bspline coeffs, 1,2rd colums removed

private:
	std::vector<Eigen::Matrix4d> m_Cs;
	std::vector<Vector6d> m_dCs;
	Eigen::Matrix4d evalQ(double q)const;
	void evalS(double q, Eigen::Matrix4d &S, Eigen::Matrix4d &dSdq);

};

const Eigen::Matrix4d JointSplineCurve::m_B = getB();
const Matrix4x3d JointSplineCurve::m_B1 = getB1();
const Matrix4x2d JointSplineCurve::m_B2 = getB2();


#endif //MUSCLEMASS_SRC_JOINTSPLINECURVE_H_