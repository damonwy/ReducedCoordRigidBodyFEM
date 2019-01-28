#include "rmpch.h"
#include "SoftBodyInvertibleFEM.h"

#include "Node.h"
#include "Body.h"
#include "Tetrahedron.h"

using namespace std;
using namespace Eigen;

SoftBodyInvertibleFEM::SoftBodyInvertibleFEM() {
	m_isInverted = false;
}

SoftBodyInvertibleFEM::SoftBodyInvertibleFEM(double density, double young, double poisson, Material material):
SoftBody(density, young, poisson, material)
{
	m_isInverted = false;
	//m_isGravity = true;
	m_type = 1;
}

void SoftBodyInvertibleFEM::computeForce_(Vector3d grav, VectorXd &f) {
	m_isInverted = false;

	// Computes force vector
	if (m_isGravity) {
		for (int i = 0; i < (int)m_nodes.size(); i++) {
			int idxM = m_nodes[i]->idxM;
			double m = m_nodes[i]->m;
			f.segment<3>(idxM) += m * grav;
		}
	}

	// Elastic Forces
	if (m_isElasticForce) {
#pragma omp parallel for num_threads(getThreadsNumber((int)m_tets.size(), MIN_ITERATOR_NUM))
		for (int i = 0; i < (int)m_tets.size(); i++) {
			auto tet = m_tets[i];
			tet->computeElasticForces();			
		}

		for (int i = 0; i < (int)m_tets.size(); i++) {
			auto tet = m_tets[i];
			tet->assembleGlobalForceVector(f);
			m_isInverted = tet->checkInverted();
		}
	}
}



void SoftBodyInvertibleFEM::computeStiffness_(MatrixXd &K) {
#pragma omp parallel for num_threads(getThreadsNumber((int)m_tets.size(), MIN_ITERATOR_NUM))
	for (int i = 0; i < (int)m_tets.size(); i++) {
		auto tet = m_tets[i];
		tet->computeForceDifferentials();
	}

	for (int i = 0; i < (int)m_tets.size(); i++) {
		auto tet = m_tets[i];
		tet->assembleGlobalStiffnessMatrixDense(K);
	}
}

void SoftBodyInvertibleFEM::computeStiffnessSparse_(vector<T> &K_) {
	//VectorXd df(3 * m_nodes.size());
	//VectorXd Dx = df;
	//for (int i = 0; i < (int)m_tets.size(); i++) {
	//	auto tet = m_tets[i];
	//	for (int ii = 0; ii < 4; ii++) {
	//		auto node = tet->m_nodes[ii];
	//		int id = node->i;
	//		int col = node->idxM;

	//		for (int iii = 0; iii < 3; iii++) {
	//			df.setZero();
	//			Dx.setZero();
	//			Dx(3 * id + iii) = 1.0;
	//			
	//			//tet->computeForceDifferentials(Dx, df);
	//			int irow = col - 3 * id;
	//			int icol = col + iii;
	//			
	//			tet->computeForceDifferentialsSparse(Dx, irow, icol, K_);

	//			//tet->computeInvertibleForceDifferentials(Dx, df);
	//			//cout << df << endl;
	//		}
	//	}
	//}

#pragma omp parallel for num_threads(getThreadsNumber((int)m_tets.size(), MIN_ITERATOR_NUM))
	for (int i = 0; i < (int)m_tets.size(); i++) {
		auto tet = m_tets[i];
		tet->computeForceDifferentials();
	}

	for (int i = 0; i < (int)m_tets.size(); i++) {
		auto tet = m_tets[i];
		tet->assembleGlobalStiffnessMatrixSparse(K_);
	}
}