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

}

void MeshEmbedding::load(const std::string &RESOURCE_DIR, const std::string &COARSE_MESH_NAME, const std::string &DENSE_MESH_NAME) {

	m_dense_mesh->load(RESOURCE_DIR, DENSE_MESH_NAME);
	m_coarse_mesh->load(RESOURCE_DIR, COARSE_MESH_NAME);

}

void MeshEmbedding::init() {
	m_dense_mesh->init();
	m_coarse_mesh->init();
}

void MeshEmbedding::precomputeWeights() {
	const std::vector<std::shared_ptr<Tetrahedron> > &coarse_mesh_tets = m_coarse_mesh->getTets();
	const std::vector<std::shared_ptr<Node> > &dense_mesh_nodes = m_dense_mesh->getNodes();
	for (int i = 0; i < (int)coarse_mesh_tets.size(); i++) {
		auto tet = coarse_mesh_tets[i];
		for (int j = 0; j < (int)dense_mesh_nodes.size(); j++) {
			auto node = dense_mesh_nodes[j];
			if (tet->checkPointInside(node)) {
				// the point is inside the tet
				tet->addEnclosedPoint(node);
				tet->computeBarycentricWeightAndSave(node);
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