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
	virtual void precompute();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;
	Matrix3d computeDeformationGradient();
	Matrix3d computeDeformationGradientDifferential(Eigen::VectorXd dx);

	Matrix3d computePKStress(Matrix3d F, double mu, double lambda);
	Matrix3d computePKStressDerivative(Matrix3d F, Matrix3d dF, double mu, double lambda);

	virtual void computeForceDifferentials(Eigen::MatrixXd &K_global);
	virtual void computeForceDifferentials(Eigen::VectorXd dx, Eigen::VectorXd &df);
	virtual void computeForceDifferentialsSparse(Eigen::VectorXd dx, int row, int col, std::vector<T> &K_);
	virtual void computeForceDifferentialsSparse(std::vector<T> &K_) {}

	virtual void computeElasticForces(Eigen::VectorXd &f);

	bool checkSameSide(const std::shared_ptr<Node> &v0, const std::shared_ptr<Node> &v1, const std::shared_ptr<Node> &v2, const std::shared_ptr<Node> &v3, const std::shared_ptr<Node> &p);
	bool checkPointInside(const std::shared_ptr<Node>& p);
	void addEnclosedPoint(const std::shared_ptr<Node>& p) { m_enclosed_points.push_back(p); }
	Vector4d computeBarycentricWeightAndSave(const std::shared_ptr<Node>& p);
	double computeEnergy();
	double ScalarTripleProduct(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c);
	Vector3d computePositionByBarycentricWeight(Vector4d weight);
	std::vector<std::shared_ptr<Node>> m_nodes;	// i, j, k, l
	int i;
	bool m_isInverted;
	std::vector<Vector4d> m_barycentric_weights;	
	// a list of material points which are enclosed by this tetrahedral element
	std::vector<std::shared_ptr<Node> > m_enclosed_points; 

protected:
	Material m_material;
	double m_young;
	double m_poisson;
	Vector3f m_enclosed_color;

	// Lame coefficients
	double m_mu;
	double m_lambda;

	double m_mass;
	double m_density;
	double m_energy;
	double psi;		// strain energy per unit undeformed volume

	// precomputed
	double W;			// undeformed volume of Te
	Matrix3d Dm;		// reference shape matrix ("material-space" shape matrix) is constant 
						// The Dm is a 3x3 matrix where the columns are the edge vectors of a tet in rest configuration
	Matrix3d Bm;		// Dm.inv()
	Matrix3d BmT;
	//Matrix3x4d Nm;		// area-weighted vertex normals [Irving 04]	(Bm in paper)					
	
	// updated each step
	Matrix3d Ds;		// deformed shape matrix
	Matrix3d F;		// deformation gradient
	Matrix3d P;		// Piola stress
	Matrix3d H;		// forces matrix
					// for force differentials
	Matrix3d dDs;
	Matrix3d dF;		// the differential of the deformation gradient
	Matrix3d dP;		// the stress differential 
	Matrix3d dH;		// the nodal force differential of the first three vertices
	Matrix12d K;			// 12x12 stiffness matrix	
	Matrix3d FTF;
};