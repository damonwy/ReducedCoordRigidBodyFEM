#pragma once
#ifndef MUSCLEMASS_SRC_SOFTBODYINVERTIBLEFEM_H_
#define MUSCLEMASS_SRC_SOFTBODYINVERTIBLEFEM_H_

#include "SoftBody.h"

class SoftBodyInvertibleFEM : public SoftBody {

public:
	SoftBodyInvertibleFEM();
	SoftBodyInvertibleFEM(double density, double young, double poisson, Material material);
	virtual ~SoftBodyInvertibleFEM() {};
	void computeStiffness(Eigen::MatrixXd &K);
	void computeForce(Eigen::Vector3d grav, Eigen::VectorXd &f);
	
protected:



private:

};

#endif // MUSCLEMASS_SRC_SOFTBODYINVERTIBLEFEM_H_