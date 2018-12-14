#pragma once
#ifndef REDUCEDCOORD_SRC_TETRAHEDRONINVERTIBLE_H_
#define REDUCEDCOORD_SRC_TETRAHEDRONINVERTIBLE_H_

#include "Tetrahedron.h"
 
class TetrahedronInvertible : public Tetrahedron {

public:
	TetrahedronInvertible() {}
	TetrahedronInvertible(double young, double poisson, double density, Material material, const std::vector<std::shared_ptr<Node>> &nodes);
	virtual ~TetrahedronInvertible() {}
	bool checkNecessityForSVD(double deltaL, double deltaU, Matrix3d F);
	//Matrix3x4d computeAreaWeightedVertexNormals();
	void computeElasticForces(Eigen::VectorXd &f);

	void computeForceDifferentials(Eigen::MatrixXd &K);
	void computeForceDifferentials(Eigen::VectorXd dx, Eigen::VectorXd &df);
	void computeForceDifferentials(Eigen::VectorXd dx, int row, int col, Eigen::MatrixXd &K);
	void computeForceDifferentialsSparse(Eigen::VectorXd dx, int row, int col, std::vector<T> &K_);

	Matrix3d clampHessian(Matrix3d &hessian, int clamped);
	void setInvertiblity(bool isInvertible) { m_isInvertible = isInvertible; }
	void precompute();
protected:
	bool m_isInvertible;

	// a given small deformation range where we do not need the diagonalization process
	double m_delta_L;
	double m_delta_U;
	bool m_isSVD;
	int clamped;

	// SVD 
	Matrix3d U;
	Matrix3d V;
	Matrix3d UT;
	Matrix3d VT;

	Matrix3d Fhat;	// diagonalized deformation gradient
	Vector3d Fhats;
	Matrix3d Phat;
	Matrix3d dPhat;

};

#endif // REDUCEDCOORD_SRC_TETRAHEDRONINVERTIBLE_H_
