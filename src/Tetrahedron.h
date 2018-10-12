#pragma once
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "MLCommon.h"

class Node;

class Tetrahedron
{
public:
	Tetrahedron();
	Tetrahedron(double young, double poisson, double density, const std::vector<std::shared_ptr<Node>> &nodes);

	virtual ~Tetrahedron();


	Eigen::Matrix3d computePKStress(Eigen::Matrix3d F, double mu, double lambda);
	Eigen::Matrix3d computePKStressDerivative(Eigen::Matrix3d F, Eigen::Matrix3d dF, double mu, double lambda);
	void computeForceDifferentials(Eigen::VectorXd dx, Eigen::VectorXd &df);
	Eigen::VectorXd computeElasticForces(Eigen::VectorXd);
	std::vector<std::shared_ptr<Node>> m_nodes;	// i, j, k, l
	bool isInverted();
	void diagDeformationGradient(Eigen::Matrix3d F);

	int i;
private:

	double m_young;
	double m_poisson;
	// Lame coefficients
	double m_mu;
	double m_lambda;
	double m_mass;
	double m_density;
	
	// precomputed
	double W;	// undeformed volume of Te
	Eigen::Matrix3d Dm;		// reference shape matrix ("material-space" shape matrix) is constant 
	Eigen::Matrix3d Bm;		// Dm.inv()

							// updated each step
	Eigen::Matrix3d Ds;		// deformed shape matrix
	Eigen::Matrix3d F;		// deformation gradient
	Eigen::Matrix3d P;		// Piola stress
	Eigen::Matrix3d H;		// forces matrix

	// SVD 
	Eigen::Matrix3d U;
	Eigen::Matrix3d V;
	Eigen::Matrix3d S;		// FTF
	Eigen::Matrix3d Fhat;	// diagonalized deformation gradient

							// for force differentials
	Eigen::Matrix3d dDs;
	Eigen::Matrix3d dF;		// the differential of the deformation gradient
	Eigen::Matrix3d dP;		// the stress differential 
	Eigen::Matrix3d dH;		// the nodal force differential of the first three vertices
	Matrix12d K;		// 12x12 stiffness matrix
	
};