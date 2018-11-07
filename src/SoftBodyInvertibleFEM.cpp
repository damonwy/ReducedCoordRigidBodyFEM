#include "SoftBodyInvertibleFEM.h"

#include <iostream>

#include "Node.h"
#include "Shape.h"
#include "Program.h"
#include "Body.h"
#include "MatrixStack.h"
#include "Tetrahedron.h"

using namespace std;
using namespace Eigen;

SoftBodyInvertibleFEM::SoftBodyInvertibleFEM() {
	m_isInvert = false;
}

SoftBodyInvertibleFEM::SoftBodyInvertibleFEM(double density, double young, double poisson, Material material):
SoftBody(density, young, poisson, material)
{
	m_isInvert = false;
	m_isGravity = true;
}

void SoftBodyInvertibleFEM::computeForce(Vector3d grav, VectorXd &f) {
	// Computes force vector

	if (m_isGravity) {
		for (int i = 0; i < (int)m_nodes.size(); i++) {
			int idxM = m_nodes[i]->idxM;
			double m = m_nodes[i]->m;
			f.segment<3>(idxM) += m * grav;
		}
	}
	//m_isInvert = false;
	// Elastic Forces
	if (m_isElasticForce) {
		for (int i = 0; i < (int)m_tets.size(); i++) {
			auto tet = m_tets[i];
			f = tet->computeInvertibleElasticForces(f);
			if (tet->isInvert) {
				//m_isInvert = true;
			}
		}
	}


	if (next != nullptr) {
		next->computeForce(grav, f);
	}
}


void SoftBodyInvertibleFEM::computeStiffness(MatrixXd &K) {
	VectorXd df(3 * m_nodes.size());
	VectorXd Dx = df;

	for (int i = 0; i < (int)m_tets.size(); i++) {
		auto tet = m_tets[i];
		for (int ii = 0; ii < 4; ii++) {
			auto node = tet->m_nodes[ii];
			int id = node->i;
			int col = node->idxM;

			for (int iii = 0; iii < 3; iii++) {
				df.setZero();
				Dx.setZero();
				Dx(3 * id + iii) = 1.0;
				tet->computeForceDifferentials(Dx, df);
				//tet->computeInvertibleForceDifferentials(Dx, df);
				//K.col(col + iii) += df;
				K.block(col - 3 * id, col + iii, 3 * m_nodes.size(), 1) += df;
			}
		}
	}
	//K.setZero();

	if (next != nullptr) {
		next->computeStiffness(K);
	}

}

