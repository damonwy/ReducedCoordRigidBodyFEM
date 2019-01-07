#pragma once
// MeshEmbeddingNull 

#ifndef REDUCEDCOORD_SRC_MESHEMBEDDINGNULL_H_
#define REDUCEDCOORD_SRC_MESHEMBEDDINGNULL_H_

#include "MeshEmbedding.h"
class MeshEmbeddingNull : public MeshEmbedding {

public:
	MeshEmbeddingNull() : MeshEmbedding() {
		m_dense_mesh = nullptr;
		m_coarse_mesh = nullptr;
	}
	void init() {}

	void gatherDofs(Eigen::VectorXd &y, int nr) {}
	void scatterDofs(Eigen::VectorXd &y, int nr) {}
	void scatterDDofs(Eigen::VectorXd &ydot, int nr) {}
	virtual ~MeshEmbeddingNull() {}
	void computeForce(Vector3d grav, Eigen::VectorXd &f) {}
	void computeStiffnessSparse(std::vector<T> &K_) {}
	void computeForceDamping(Eigen::VectorXd &f, Eigen::MatrixXd &D) {}
	void computeForceDampingSparse(Eigen::VectorXd &f, std::vector<T> &D_) {}
	void countDofs(int &nm, int &nr) {}
	void computeMassSparse(std::vector<T> &M_) {}
	void computeJacobianSparse(std::vector<T> &J_){}
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const {}
};

#endif // REDUCEDCOORD_SRC_MESHEMBEDDINGNULL_H_