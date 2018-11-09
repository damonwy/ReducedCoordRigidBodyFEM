#pragma once
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "MLCommon.h"

class Node;
class MatrixStack;
class Program;

class Tetrahedron
{
public:
	Tetrahedron();
	Tetrahedron(double young, double poisson, double density, Material material, const std::vector<std::shared_ptr<Node>> &nodes);
	virtual ~Tetrahedron();

	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;

	Eigen::Matrix3d computeDeformationGradient();

	Eigen::Matrix3d computePKStress(Eigen::Matrix3d F, double mu, double lambda);
	Eigen::Matrix3d computePKStressDerivative(Eigen::Matrix3d F, Eigen::Matrix3d dF, double mu, double lambda);
	void computeForceDifferentials(Eigen::VectorXd dx, Eigen::VectorXd &df);
	// Vector12d computeForceDifferentials(Vector12d dx, Vector12d &df);
	Eigen::VectorXd computeElasticForces(Eigen::VectorXd f);

	// Functions for Invertible FEM 
	Matrix3x4d computeAreaWeightedVertexNormals();
	Eigen::VectorXd computeInvertibleElasticForces(Eigen::VectorXd f);
	void computeInvertibleForceDifferentials(Eigen::VectorXd dx, Eigen::VectorXd &df);
	Eigen::Matrix3d computeInvertiblePKStress(Eigen::Matrix3d F, double mu, double lambda);
	Eigen::Matrix3d computeInvertiblePKStressDerivative(Eigen::Matrix3d F, Eigen::Matrix3d dF, double mu, double lambda);
	void compute_dPdF();
	void compute_dGdF();
	void compute_dFdU();


	double computeEnergy();
	std::vector<std::shared_ptr<Node>> m_nodes;	// i, j, k, l
	bool isInverted();
	void diagDeformationGradient(Eigen::Matrix3d F);
	void setInvertiblity(bool isInvertible) { m_isInvertible = isInvertible; }
	bool isInvert;
	int i;

	Matrix9d dPdF;
	Matrix9d dGdF;
	int clamped;
private:
	bool m_isInvertible;
	Material m_material;
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
							// The Dm is a 3x3 matrix where the columns are the edge vectors of a tet in rest configuration
	Eigen::Matrix3d Bm;		// Dm.inv()
	Matrix3x4d Nm;		// area-weighted vertex normals [Irving 04]	(Bm in paper)					
	
	// updated each step
	Eigen::Matrix3d Ds;		// deformed shape matrix
	Eigen::Matrix3d F;		// deformation gradient
	Eigen::Matrix3d P;		// Piola stress
	Eigen::Matrix3d H;		// forces matrix

	double psi;		// strain energy per unit undeformed volume
	double m_energy;

	// SVD 
	Eigen::Matrix3d U;
	Eigen::Matrix3d V;
	Eigen::Matrix3d Fhat;	// diagonalized deformation gradient
	Eigen::Matrix3d Phat;	
							// for force differentials
	Eigen::Matrix3d dDs;
	Eigen::Matrix3d dF;		// the differential of the deformation gradient
	Eigen::Matrix3d dP;		// the stress differential 
	Eigen::Matrix3d dPhat;
	Eigen::Matrix3d dH;		// the nodal force differential of the first three vertices
	Matrix12d K;			// 12x12 stiffness matrix
	
};