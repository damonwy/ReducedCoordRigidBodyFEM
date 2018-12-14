#pragma once

#ifndef REDUCEDCOORD_SRC_TETRAHEDRONCOROTATIONAL_H_
#define REDUCEDCOORD_SRC_TETRAHEDRONCOROTATIONAL_H_
#define EIGEN_USE_MKL_ALL

#include "Tetrahedron.h"

class TetrahedronCorotational : public Tetrahedron {

public:
	TetrahedronCorotational() {}
	TetrahedronCorotational(double young, double poisson, double density, Material material, const std::vector<std::shared_ptr<Node>> &nodes);
	virtual ~TetrahedronCorotational() {}
	void computeElasticForces(Eigen::VectorXd &f);
	void computeForceDifferentials(Eigen::MatrixXd &K);
	void computeForceDifferentialsSparse(std::vector<T> &K_);
	void precompute();

protected:
	Matrix12d Ke;		// precomputed stiffness matrix for coratational linear material

	// Corotational
	Matrix3d R;
	Matrix12d Re;
	Vector12d xx;		// undeformed position vector
	Vector12d Kexx;		// undeformed force
	Matrix12d RKR;
	Vector12d RKxx;
	Vector12d Kx;		// deformed force
	Vector12d x_vec;	// deformed positon vector
};

#endif // REDUCEDCOORD_SRC_TETRAHEDRONCOROTATIONAL_H_
