#pragma once
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "MLCommon.h"

class Node;
class MatrixStack;
class Program;

typedef Eigen::Triplet<double> T;
class Tetrahedron
{
public:
	Tetrahedron();
	Tetrahedron(double young, double poisson, double density, Material material, const std::vector<std::shared_ptr<Node>> &nodes);
	virtual ~Tetrahedron() {}

	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;
	Matrix3d computeDeformationGradient();

	Matrix3d computePKStress(Matrix3d F, double mu, double lambda);
	Matrix3d computePKStressDerivative(Matrix3d F, Matrix3d dF, double mu, double lambda);
	void computeForceDifferentials(Eigen::VectorXd dx, Eigen::VectorXd &df);
	void computeForceDifferentialsSparse(Eigen::VectorXd dx, int row, int col, std::vector<T> &K_);
	Eigen::VectorXd computeElasticForces(Eigen::VectorXd f);

	// Functions for Invertible FEM 
	Matrix3x4d computeAreaWeightedVertexNormals();
	Eigen::VectorXd computeInvertibleElasticForces(Eigen::VectorXd f);
	void computeInvertibleForceDifferentials(Eigen::VectorXd dx, Eigen::VectorXd &df);
	void computeInvertibleForceDifferentialsSparse(Eigen::VectorXd dx, int row, int col, std::vector<T> &K_);

	Matrix3d computeInvertiblePKStress(Matrix3d F, double mu, double lambda);
	Matrix3d computeInvertiblePKStressDerivative(Matrix3d F, Matrix3d dF, double mu, double lambda);

	double computeEnergy();
	std::vector<std::shared_ptr<Node>> m_nodes;	// i, j, k, l
	void diagDeformationGradient(Matrix3d F);
	void setInvertiblity(bool isInvertible) { m_isInvertible = isInvertible; }
	bool isInvert;
	int i;
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
	Matrix3d Dm;		// reference shape matrix ("material-space" shape matrix) is constant 
							// The Dm is a 3x3 matrix where the columns are the edge vectors of a tet in rest configuration
	Matrix3d Bm;		// Dm.inv()
	Matrix3x4d Nm;		// area-weighted vertex normals [Irving 04]	(Bm in paper)					
	
	// updated each step
	Matrix3d Ds;		// deformed shape matrix
	Matrix3d F;		// deformation gradient
	Matrix3d P;		// Piola stress
	Matrix3d H;		// forces matrix

	double psi;		// strain energy per unit undeformed volume
	double m_energy;

	// SVD 
	Matrix3d U;
	Matrix3d V;
	Matrix3d Fhat;	// diagonalized deformation gradient
	Vector3d Fhats;
	Matrix3d Phat;	
							// for force differentials
	Matrix3d dDs;
	Matrix3d dF;		// the differential of the deformation gradient
	Matrix3d dP;		// the stress differential 
	Matrix3d dPhat;
	Matrix3d dH;		// the nodal force differential of the first three vertices
	Matrix12d K;			// 12x12 stiffness matrix	
};