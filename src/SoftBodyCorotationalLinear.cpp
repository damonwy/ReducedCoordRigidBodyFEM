#include "SoftBodyCorotationalLinear.h"

#include <iostream>

#include "Node.h"
#include "Shape.h"
#include "Program.h"
#include "Body.h"
#include "MatrixStack.h"
#include "Tetrahedron.h"
#include "MatlabDebug.h"

using namespace std;
using namespace Eigen;

SoftBodyCorotationalLinear::SoftBodyCorotationalLinear() {
}

SoftBodyCorotationalLinear::SoftBodyCorotationalLinear(double density, double young, double poisson, Material material) :
	SoftBody(density, young, poisson, material)
{
	m_isGravity = true;
}

void SoftBodyCorotationalLinear::computeForce_(Vector3d grav, VectorXd &f) {
	// Computes force vector
	if (m_isGravity) {
		for (int i = 0; i < (int)m_nodes.size(); i++) {
			int idxM = m_nodes[i]->idxM;
			double m = m_nodes[i]->m;
			f.segment<3>(idxM) += m * grav;
		}
	}
}

void SoftBodyCorotationalLinear::computeStiffness_(MatrixXd &K) {
	for (int i = 0; i < (int)m_tets.size(); i++) {
		auto tet = m_tets[i];
		tet->computeCorotationalForceDifferentials(K);
	}

}

void SoftBodyCorotationalLinear::computeStiffnessSparse_(vector<T> &K_) {
	for (int i = 0; i < (int)m_tets.size(); i++) {
		auto tet = m_tets[i];		
		tet->computeCorotationalForceDifferentialsSparse(K_);

	}

}