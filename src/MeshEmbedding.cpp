#include "MeshEmbedding.h"
#include "SoftBodyInvertibleFEM.h"
#include "SoftBodyCorotationalLinear.h"
#include "Node.h"
#include "Tetrahedron.h"
#include "Body.h"
#include "FaceTriangle.h"
#include "Line.h"
#include "Surface.h"

using namespace std;
using namespace Eigen;

MeshEmbedding::MeshEmbedding(double density, double young, double possion, Material material, SoftBodyType type) {
	if (type == SOFT_INVERTIBLE) {
		auto dense_ptr = make_shared<Surface>();
		//auto dense_ptr = make_shared<SoftBodyInvertibleFEM>(density, young, possion, material);
		auto coarse_ptr = make_shared<SoftBodyInvertibleFEM>(density, young, possion, material);
		m_dense_mesh = dense_ptr;
		m_coarse_mesh = coarse_ptr;
	}
	else if (type == SOFT_COROTATED) {
		auto dense_ptr = make_shared<Surface>();
		//auto dense_ptr = make_shared<SoftBodyCorotationalLinear>(density, young, possion, material);
		m_dense_mesh = dense_ptr;
		auto coarse_ptr = make_shared<SoftBodyCorotationalLinear>(density, young, possion, material);
		m_coarse_mesh = coarse_ptr;
	}
	else {
		m_dense_mesh = make_shared<Surface>();
		//m_dense_mesh = make_shared<SoftBody>(density, young, possion, material);
		m_coarse_mesh = make_shared<SoftBody>(density, young, possion, material);
	}
	m_isDenseMesh = true;
	m_isCoarseMesh = false;
}

void MeshEmbedding::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
	if (m_isDenseMesh) {
		//m_dense_mesh->draw(MV, prog, progSimple, P);
		m_dense_mesh->draw(MV, prog, P);
	}

	if (m_isCoarseMesh) {
		m_coarse_mesh->draw(MV, prog, progSimple, P);
	}

	if (next != nullptr) {
		next->draw(MV, prog, progSimple, P);
	}
}

void MeshEmbedding::load(const string &RESOURCE_DIR, const string &COARSE_MESH_NAME, const string &TETGEN_FLAGS_0, const string &DENSE_MESH_NAME, const string &TETGEN_FLAGS_1) {
	m_coarse_mesh->load(RESOURCE_DIR, COARSE_MESH_NAME, TETGEN_FLAGS_0);
	m_dense_mesh->load(RESOURCE_DIR, DENSE_MESH_NAME, TETGEN_FLAGS_1);
}

void MeshEmbedding::init() {
	m_dense_mesh->init();
	m_coarse_mesh->init();
}

void MeshEmbedding::computeMassSparse(vector<T> &M_) {
	m_coarse_mesh->computeMassSparse(M_);

	if (next != nullptr) {
		next->computeMassSparse(M_);
	}
}

void MeshEmbedding::computeJacobianSparse(vector<T> &J_) {
	m_coarse_mesh->computeJacobianSparse(J_);

	if (next != nullptr) {
		next->computeJacobianSparse(J_);
	}
}

void MeshEmbedding::computeForce(Vector3d grav, VectorXd &f) {
	m_coarse_mesh->computeForce(grav, f);
	if (next != nullptr) {
		next->computeForce(grav, f);
	}
}

void MeshEmbedding::computeForceDamping(VectorXd &f, MatrixXd &D) {
	m_coarse_mesh->computeForceDamping(f, D);

	if (next != nullptr) {
		next->computeForceDamping(f, D);
	}
}

void MeshEmbedding::computeForceDampingSparse(VectorXd &f, vector<T> &D_) {
	m_coarse_mesh->computeForceDampingSparse(f, D_);
	if (next != nullptr) {
		next->computeForceDampingSparse(f, D_);
	}
}

void MeshEmbedding::computeStiffnessSparse(vector<T> &K_) {
	m_coarse_mesh->computeStiffnessSparse(K_);
	if (next != nullptr) {
		next->computeStiffnessSparse(K_);
	}
}

void MeshEmbedding::scatterDofs(VectorXd &y, int nr) {
	m_coarse_mesh->scatterDofs(y, nr);

	const vector<std::shared_ptr<Tetrahedron> > &coarse_mesh_tets = m_coarse_mesh->getTets();
	bool isCheckingCollision = m_dense_mesh->m_isCollisionWithFloor;

	// update dense mesh using coarse mesh
	for (int i = 0; i < (int)coarse_mesh_tets.size(); i++) {
		auto tet = coarse_mesh_tets[i];
		int num_enclosed = tet->m_enclosed_points.size();
		for (int j = 0; j < num_enclosed; j++) {
			// update nodes
			auto node = tet->m_enclosed_points[j];
			node->x = tet->computePositionByBarycentricWeight(tet->m_barycentric_weights[j]);
			if (isCheckingCollision) {

				if (node->x.y() < m_dense_mesh->m_floor_y) {
					// there is collision!

					// compute the penetration
					double deep_inside = m_dense_mesh->m_floor_y - node->x.y();

					// find the lowest node in the coarse tet
					double y_lowest = 1000.0;
					int update_id = -1;
					for (int t = 0; t < 4; t++) {
						auto pt = tet->m_nodes[t];
						double ypos = pt->x.y();
						if (ypos < y_lowest) {
							y_lowest = ypos;
							update_id = t;
						}
					}

					assert(update_id > -1);
					// push up that lowest node and kill its velocity
					tet->m_nodes[update_id]->x.y() += deep_inside;
					tet->m_nodes[update_id]->v.y() = 0.0;

					
					// kill coarse y velocity in y vector
					int idxR = tet->m_nodes[update_id]->idxR;
					y.segment<3>(idxR).y() += deep_inside;
					y.segment<3>(nr + idxR).y() = 0.0;

					// update the node in dense mesh that collides with floor and all the points before that
					//node->x = tet->computePositionByBarycentricWeight(tet->m_barycentric_weights[j]);

					for (int k = 0; k < j+1; k++) {
						auto node_update = tet->m_enclosed_points[k];
						  node_update->x = tet->computePositionByBarycentricWeight(tet->m_barycentric_weights[k]);
					}
				}
			}
		}
	}




	if (next != nullptr) {
		next->scatterDofs(y, nr);
	}
}

void MeshEmbedding::scatterDDofs(VectorXd &ydot, int nr) {

	m_coarse_mesh->scatterDDofs(ydot, nr);
	//const vector<std::shared_ptr<Tetrahedron> > &coarse_mesh_tets = m_coarse_mesh->getTets();
	//bool isCheckingCollision = m_dense_mesh->m_isCollisionWithFloor;
	//// update dense mesh using coarse mesh
	//for (int i = 0; i < (int)coarse_mesh_tets.size(); i++) {
	//	auto tet = coarse_mesh_tets[i];
	//	int num_enclosed = tet->m_enclosed_points.size();
	//	for (int j = 0; j < num_enclosed; j++) {
	//		// update nodes
	//		auto node = tet->m_enclosed_points[j];
	//		node->x = tet->computePositionByBarycentricWeight(tet->m_barycentric_weights[j]);
	//	}
	//}

	m_dense_mesh->updatePosNor();
	if (next != nullptr) {
		next->scatterDDofs(ydot, nr);
	}
}

void MeshEmbedding::gatherDofs(VectorXd &y, int nr) {
	m_coarse_mesh->gatherDofs(y, nr);
	if (next != nullptr) {
		next->gatherDofs(y, nr);
	}
}

void MeshEmbedding::precomputeWeights() {
	const vector<shared_ptr<Tetrahedron> > &coarse_mesh_tets = m_coarse_mesh->getTets();
	const vector<shared_ptr<FaceTriangle> > &dense_mesh_trifaces = m_dense_mesh->getFaces();
	//for (int i = 0; i < (int)coarse_mesh_tets.size(); i++) {
	//	auto tet = coarse_mesh_tets[i];
	//	for (int j = 0; j < (int)dense_mesh_nodes.size(); j++) {
	//		auto node = dense_mesh_nodes[j];
	//		if (!node->isEnclosedByTet && tet->checkPointInside(node)) {
	//			// the point is inside the tet
	//			tet->addEnclosedPoint(node);
	//			Vector4d weight = tet->computeBarycentricWeightAndSave(node);
	//			
	//			node->isEnclosedByTet = true;
	//		}
	//	}
	//}

	for (int i = 0; i < (int)coarse_mesh_tets.size(); i++) {
		auto tet = coarse_mesh_tets[i];
		for (int j = 0; j < (int)dense_mesh_trifaces.size(); j++) {
			for (int k = 0; k < 3; k++) {
				auto node = dense_mesh_trifaces[j]->m_nodes[k];
				if (!node->isEnclosedByTet && tet->checkPointInside(node)) {
					// the point is inside the tet
					tet->addEnclosedPoint(node);
					Vector4d weight = tet->computeBarycentricWeightAndSave(node);

					node->isEnclosedByTet = true;
				}
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

void MeshEmbedding::setAttachmentsByXZCircle(double y, double range, Vector2d O, double r, shared_ptr<Body> body) {
	m_coarse_mesh->setAttachmentsByXZCircle(y, range, O, r, body);
}

void  MeshEmbedding::setAttachmentsByLine(shared_ptr<Line> l) {
	m_coarse_mesh->setAttachmentsByLine(l);

}