#pragma once
// JointSplineSurface
//    Uses cubic B-splines.
//    

#ifndef REDUCEDCOORD_SRC_JOINTSPLINESURFACE_H_
#define REDUCEDCOORD_SRC_JOINTSPLINESURFACE_H_

#include "Joint.h"

class Body;

class JointSplineSurface : public Joint {

public:
	JointSplineSurface();
	JointSplineSurface(std::shared_ptr<Body> body, std::shared_ptr<Joint> parent = nullptr);
	virtual ~JointSplineSurface() {}

	void addControlFrame(int i, int j, Vector6d C);
	static double Cfun(Matrix4d C, Vector2d q);
	static double dCfun(Matrix4d C, int i, Vector2d q);
	static double d2Cfun(Matrix4d C, int i, int j, Vector2d q);

	static const Matrix4d& getB() {
		static Matrix4d _B(1.0 / 6.0 * (Matrix4d() <<
			1.0, -3.0, 3.0, -1.0,
			4.0, 0.0, -6.0, 3.0,
			1.0, 3.0, 3.0, -3.0,
			0.0, 0.0, 0.0, 1.0).finished());
		return _B;
	}

	static const Matrix6d& getE() {
		static Matrix6d _E((Matrix6d() <<
			0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 
			0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
			1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 1.0, 0.0, 0.0, 0.0).finished());
		return _E;
	}

	static const Matrix4d m_B;			// Bspline coeffs
	static const Matrix6d m_E;			// 6 basis twists in Eq.(25)

protected:
	void update_();
	void draw_(std::shared_ptr<MatrixStack> MV, 
		const std::shared_ptr<Program> prog, 
		const std::shared_ptr<Program> progSimple, 
		std::shared_ptr<MatrixStack> P) const;
private:
	Tensor4x4x6d m_cs;
	Matrix4d evalQ(Vector2d q)const;
	void evalS(Vector2d q, Eigen::MatrixXd &S, Tensor6x2x2d &dSdq);

};




#endif // REDUCEDCOORD_SRC_JOINTSPLINESURFACE_H_