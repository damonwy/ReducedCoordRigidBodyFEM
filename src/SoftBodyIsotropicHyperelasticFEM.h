#pragma once
#ifndef MUSCLEMASS_SRC_SOFTBODYISOTROPICHYPERELASTICFEM_H_
#define MUSCLEMASS_SRC_SOFTBODYISOTROPICHYPERELASTICFEM_H_

#include "SoftBody.h"

class SoftBodyIsotropicHyperelasticFEM : public SoftBody {
public:
	SoftBodyIsotropicHyperelasticFEM();
	SoftBodyIsotropicHyperelasticFEM(double density, double young, double poisson, Material material);
	virtual ~SoftBodyIsotropicHyperelasticFEM() {};


protected:
	void computeForce_(Vector3d grav, Eigen::VectorXd &f);
	void computeStiffnessSparse_(std::vector<T> &K_);
	void computeStiffness_(Eigen::MatrixXd &K);



};



#endif // MUSCLEMASS_SRC_SOFTBODYISOTROPICHYPERELASTICFEM_H_