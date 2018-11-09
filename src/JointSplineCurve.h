#pragma once
// JointSplineCurve 
//    Uses cubic B-splines.
//    Assumes a cyclic spline curve

#ifndef REDUCEDCOORD_SRC_JOINTSPLINECURVE_H_
#define REDUCEDCOORD_SRC_JOINTSPLINECURVE_H_

#include "Joint.h"

class Body;
class Shape;

class JointSplineCurve : public Joint {

public:
	JointSplineCurve();
	JointSplineCurve(std::shared_ptr<Body> body, std::shared_ptr<Joint> parent = nullptr);
	virtual ~JointSplineCurve() {}
	void load(const std::string &RESOURCE_DIR, std::string joint_shape);
	void init(int &nm, int &nr);
	void addControlFrame(Eigen::Matrix4d C);

	static double Bsum(int i, double q);
	static double dBsum(int i, double q);
	static double d2Bsum(int i, double q);

	static const Matrix4d& getB() {
		static Matrix4d _B(1.0 / 6.0 * (Matrix4d() << 
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

	static const Matrix4d m_B; // Bspline coeffs
	static const Matrix4x3d m_B1;	 // Bspline coeffs, 1st column removed
	static const Matrix4x2d m_B2;	 // Bspline coeffs, 1,2rd colums removed

protected:
	void update_();	
	void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;

private:
	std::vector<Matrix4d> m_Cs;
	std::vector<Vector6d> m_dCs;
	Matrix4d evalQ(double q)const;
	void evalS(double q, Vector6d &S, Vector6d &dSdq);
	std::shared_ptr<Shape> m_jointSphereShape;

};

#endif // REDUCEDCOORD_SRC_JOINTSPLINECURVE_H_