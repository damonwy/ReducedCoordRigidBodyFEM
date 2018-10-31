#pragma once
// JointSplineSurface
//    Uses cubic B-splines.
//    

#ifndef MUSCLEMASS_SRC_JOINTSPLINESURFACE_H_
#define MUSCLEMASS_SRC_JOINTSPLINESURFACE_H_

#include "Joint.h"

class Body;

class JointSplineSurface : public Joint {

public:
	JointSplineSurface();
	JointSplineSurface(std::shared_ptr<Body> body, std::shared_ptr<Joint> parent = nullptr);
	virtual ~JointSplineSurface();

	void addControlFrame(int i, int j, Vector6d C);
	/*void updateSelf();
	void drawSelf(std::shared_ptr<MatrixStack> MV, 
		const std::shared_ptr<Program> prog, 
		const std::shared_ptr<Program> progSimple, 
		std::shared_ptr<MatrixStack> P) const;*/

	static double Cfun(Eigen::Matrix4d C, Eigen::Vector2d q);
	static double dCfun(Eigen::Matrix4d C, int i, Eigen::Vector2d q);
	static double d2Cfun(Eigen::Matrix4d C, int i, int j, Eigen::Vector2d q);

	static const Eigen::Matrix4d& getB() {
		static Eigen::Matrix4d _B(1.0 / 6.0 * (Eigen::Matrix4d() <<
			1.0, -3.0, 3.0, -1.0,
			4.0, 0.0, -6.0, 3.0,
			1.0, 3.0, 3.0, -3.0,
			0.0, 0.0, 0.0, 1.0).finished());
		return _B;
	}

	static const Matrix6d& getE() {
		static Matrix6d _E(1.0 / 6.0 * (Matrix6d() <<
			0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 
			0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
			1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 1.0, 0.0, 0.0, 0.0).finished());
		return _E;
	}

	static const Eigen::Matrix4d m_B; // Bspline coeffs
	static const Matrix6d m_E;		  // 6 basis twists in Eq.(25)

private:
	//Eigen::Tensor<double, 3> m_cs;

	std::vector<Eigen::Matrix4d> m_Cs;
	std::vector<Vector6d> m_dCs;
	Eigen::Matrix4d evalQ(Eigen::Vector2d q)const;
	void evalS(Eigen::Vector2d q, Matrix6x2d &S, Eigen::Tensor<double, 3> &dSdq);

};




#endif //MUSCLEMASS_SRC_JOINTSPLINESURFACE_H_