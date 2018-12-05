#include "MeshEmbedding.h"

#include "SoftBody.h"
#include "Node.h"
#include "Tetrahedron.h"

using namespace std;
using namespace Eigen;

MeshEmbedding::MeshEmbedding(const shared_ptr<SoftBody> coarse_mesh, const shared_ptr<SoftBody> dense_mesh) {
	m_coarse_mesh = coarse_mesh;
	m_dense_mesh = dense_mesh;

}

void MeshEmbedding::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {

	m_dense_mesh->draw(MV, prog, progSimple, P);
	m_coarse_mesh->draw(MV, prog, progSimple, P);
	if (next != nullptr) {
		next->draw(MV, prog, progSimple, P);
	}
}

void MeshEmbedding::load(const string &RESOURCE_DIR, const string &COARSE_MESH_NAME, const string &DENSE_MESH_NAME) {

	m_dense_mesh->load(RESOURCE_DIR, DENSE_MESH_NAME);
	m_coarse_mesh->load(RESOURCE_DIR, COARSE_MESH_NAME);

}

void MeshEmbedding::init() {
	m_dense_mesh->init();
	m_coarse_mesh->init();
}

void MeshEmbedding::computeMassSparse(vector<T> &M_) {
	m_coarse_mesh->computeMassSparse(M_);
}

void MeshEmbedding::computeJacobianSparse(std::vector<T> &J_) {
	m_coarse_mesh->computeJacobianSparse(J_);
}

void MeshEmbedding::computeForce(Vector3d grav, VectorXd &f) {
	m_coarse_mesh->computeForce(grav, f);
}

void MeshEmbedding::computeStiffnessSparse(std::vector<T> &K_) {
	m_coarse_mesh->computeStiffnessSparse(K_);
}

void MeshEmbedding::scatterDofs(Eigen::VectorXd &y, int nr) {
	m_coarse_mesh->scatterDofs(y, nr);
	const std::vector<std::shared_ptr<Tetrahedron> > &coarse_mesh_tets = m_coarse_mesh->getTets();

	// update dense mesh using coarse mesh
	for (int i = 0; i < (int)coarse_mesh_tets.size(); i++) {
		auto tet = coarse_mesh_tets[i];
		int num_enclosed = tet->m_enclosed_points.size();
		for (int j = 0; j < num_enclosed; j++) {
			// update nodes
			auto node = tet->m_enclosed_points[j];
			node->x = tet->computePositionByBarycentricWeight(tet->m_barycentric_weights[j]);
		}
	}

	m_dense_mesh->updatePosNor();

}

void MeshEmbedding::scatterDDofs(Eigen::VectorXd &ydot, int nr) {

	m_coarse_mesh->scatterDDofs(ydot, nr);
}

void MeshEmbedding::precomputeWeights() {
	const std::vector<std::shared_ptr<Tetrahedron> > &coarse_mesh_tets = m_coarse_mesh->getTets();
	const std::vector<std::shared_ptr<Node> > &dense_mesh_nodes = m_dense_mesh->getNodes();
	for (int i = 0; i < (int)coarse_mesh_tets.size(); i++) {
		auto tet = coarse_mesh_tets[i];
		for (int j = 0; j < (int)dense_mesh_nodes.size(); j++) {
			auto node = dense_mesh_nodes[j];
			if (!node->isEnclosedByTet && tet->checkPointInside(node)) {
				// the point is inside the tet
				tet->addEnclosedPoint(node);
				tet->computeBarycentricWeightAndSave(node);
				node->isEnclosedByTet = true;
			}
		}
	}
}

void MeshEmbedding::transformCoarseMesh(Matrix4d E) {
	m_coarse_mesh->transform(E);
}

void MeshEmbedding::transformDenseMesh(Matrix4d E) {
	m_dense_mesh->transform(E);
}
void MeshEmbedding::countDofs(int &nm, int &nr) {
	m_coarse_mesh->countDofs(nm, nr);
}

void MeshEmbedding::updatePosNor() {

	m_dense_mesh->updatePosNor();
	m_coarse_mesh->updatePosNor();

}

void MeshEmbedding::setAttachmentsByYZCircle(double x, double range, Vector2d O, double r, shared_ptr<Body> body) {
	m_coarse_mesh->setAttachmentsByYZCircle(x, range, O, r, body);


}