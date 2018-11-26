#pragma once
#ifndef MUSCLEMASS_SRC_SOFTBODYINVERTIBLEFEM_H_
#define MUSCLEMASS_SRC_SOFTBODYINVERTIBLEFEM_H_

#include "SoftBody.h"

class SoftBodyInvertibleFEM : public SoftBody {

public:
	SoftBodyInvertibleFEM();
	SoftBodyInvertibleFEM(double density, double young, double poisson, Material material);
	virtual ~SoftBodyInvertibleFEM() {};

protected:
	void computeForce_(Vector3d grav, Eigen::VectorXd &f);
	void computeStiffnessSparse_(std::vector<T> &K_);
	void computeStiffness_(Eigen::MatrixXd &K);
private:

};

#endif // MUSCLEMASS_SRC_SOFTBODYINVERTIBLEFEM_H_