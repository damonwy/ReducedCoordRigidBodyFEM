#pragma once
#ifndef MUSCLEMASS_SRC_SOFTBODYCOROTATIONALLINEAR_H_
#define MUSCLEMASS_SRC_SOFTBODYCOROTATIONALLINEAR_H_

#include "SoftBody.h"

class SoftBodyCorotationalLinear : public SoftBody {

public:
	SoftBodyCorotationalLinear();
	SoftBodyCorotationalLinear(double density, double young, double poisson, Material material);
	virtual ~SoftBodyCorotationalLinear() {};

protected:
	void computeForce_(Vector3d grav, Eigen::VectorXd &f);
	void computeStiffnessSparse_(std::vector<T> &K_);
	void computeStiffness_(Eigen::MatrixXd &K);
private:

};

#endif // MUSCLEMASS_SRC_SOFTBODYCOROTATIONALLINEAR_H_