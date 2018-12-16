#pragma once
#define EIGEN_USE_MKL_ALL

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "MLCommon.h"

class Body;
class SoftBody;
class Program;
class MatrixStack;
typedef Eigen::Triplet<double> T;

class MeshEmbedding {
public:
	MeshEmbedding() {}
	MeshEmbedding(double density, double young, double possion, Material material, SoftBodyType type);

	virtual ~MeshEmbedding() {}

	void load(const std::string &RESOURCE_DIR, const std::string &COARSE_MESH_NAME, const std::string &DENSE_MESH_NAME);
	void init();
	void precomputeWeights();
	void updatePosNor();
	void countDofs(int &nm, int &nr);
	void transformCoarseMesh(Matrix4d E);
	void transformDenseMesh(Matrix4d E);
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;
	void setAttachmentsByYZCircle(double x, double range, Vector2d O, double r, std::shared_ptr<Body> body);
	void setAttachmentsByXZCircle(double y, double range, Vector2d O, double r, std::shared_ptr<Body> body);

	void computeMassSparse(std::vector<T> &M_);
	void computeJacobianSparse(std::vector<T> &J_);

	void computeForce(Vector3d grav, Eigen::VectorXd &f);
	void computeStiffnessSparse(std::vector<T> &K_);
	void scatterDofs(Eigen::VectorXd &y, int nr);
	void scatterDDofs(Eigen::VectorXd &ydot, int nr);
	void gatherDofs(Eigen::VectorXd &y, int nr);
	void toggleDrawingDenseMesh(bool isDenseMesh) { m_isDenseMesh = isDenseMesh; if (next != nullptr) { next->toggleDrawingDenseMesh(isDenseMesh); } }
	void toggleDrawingCoarseMesh(bool isCoarseMesh) { m_isCoarseMesh = isCoarseMesh; if (next != nullptr) { next->toggleDrawingCoarseMesh(isCoarseMesh); }
	}

	std::shared_ptr<SoftBody> getDenseMesh() { return m_dense_mesh; }
	std::shared_ptr<SoftBody> getCoarseMesh() { return m_coarse_mesh; }
	bool m_isDenseMesh;
	bool m_isCoarseMesh;


	std::shared_ptr<MeshEmbedding> next;
protected:
	std::shared_ptr<SoftBody> m_dense_mesh;
	std::shared_ptr<SoftBody> m_coarse_mesh;

};