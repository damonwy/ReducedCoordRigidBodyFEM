#pragma once
#ifndef MUSCLEMASS_SRC_JOINT_H_
#define MUSCLEMASS_SRC_JOINT_H_
#include <vector>
#include <memory>
#define EIGEN_DONT_ALIGN_STATICALLY

#include <Eigen/Dense>
#include <iostream>
#include "MLCommon.h"

class SE3;
class Body;

class Joint {
public:

	Joint();
	virtual ~Joint();
	virtual void init();
	virtual void draw();
	virtual void update();

	double q;
	double qdot;
	double qddot;
	double tau;		// Joint torque
	Matrix4d E_pj;	// Transform of this joint wrt parent joint
	Matrix4d E_pj0;	// Transform when q is zero
	Matrix4d E_jp;	// Transform of parent joint wrt this joint
	Matrix6d Ad_jp; // Adjoint of E_jp
	Matrix4d E_wj;	// Transform of this joint wrt world
	std::shared_ptr<Body> body;	// Attached body

	void setJointTransform(Eigen::Matrix4d E);

private:

};

#endif // MUSCLEMASS_SRC_JOINT_H_