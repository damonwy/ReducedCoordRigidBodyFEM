#pragma once

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "MLCommon.h"

class SoftBody;
class Program;
class MatrixStack;

class MeshEmbedding {
public:
	MeshEmbedding() {}
	MeshEmbedding(const std::shared_ptr<SoftBody> coarse_mesh, const std::shared_ptr<SoftBody> dense_mesh);

	virtual ~MeshEmbedding() {}

	virtual void load(const std::string &RESOURCE_DIR, const std::string &COARSE_MESH_NAME, const std::string &DENSE_MESH_NAME);
	virtual void init();
	void precomputeWeights();
	void updatePosNor();
	void countDofs(int &nm, int &nr);
	void transformCoarseMesh(Matrix4d E);
	void transformDenseMesh(Matrix4d E);
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;
	std::shared_ptr<SoftBody> getDenseMesh() { return m_dense_mesh; }
	std::shared_ptr<SoftBody> getCoarseMesh() { return m_coarse_mesh; }

protected:
	std::shared_ptr<SoftBody> m_dense_mesh;
	std::shared_ptr<SoftBody> m_coarse_mesh;

};