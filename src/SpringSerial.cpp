#include "SpringSerial.h"

#include "Body.h"
#include "Node.h"

using namespace std;
using namespace Eigen;

SpringSerial::SpringSerial() {

}


SpringSerial::SpringSerial(int n_nodes, int &countS, int &countCM):
Spring(countS, countCM),
m_K(0.0), m_mass(1.0)
{
	for (int i = 0; i < n_nodes; i++) {
		auto node = make_shared<Node>();
		node->x = Vector3d::Zero();
		node->v = Vector3d::Zero();
		node->a = Vector3d::Zero();
		m_nodes.push_back(node);

	}
}


void SpringSerial::setAttachments(shared_ptr<Body> body0, shared_ptr<Body> body1) {
	// Attaches this spring to body0 and body1
	m_body0 = body0;
	m_body1 = body1;
}

void SpringSerial::countDofs_() {


}


void SpringSerial::gatherDofs_(VectorXd &y, int nr) {
	// Gathers q and qdot into y
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		int idxR = m_nodes[i]->idxR;
		y.segment<3>(idxR) = m_nodes[i]->x;
		y.segment<3>(nr + idxR) = m_nodes[i]->v;
	}
}


void SpringSerial::gatherDDofs_(VectorXd &ydot, int nr) {
	// Gathers qdot and qddot into ydot
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		int idxR = m_nodes[i]->idxR;
		ydot.segment<3>(idxR) = m_nodes[i]->v;
		ydot.segment<3>(nr + idxR) = m_nodes[i]->a;
	}
}

void SpringSerial::scatterDofs_(VectorXd &y, int nr) {
	// Scatters q and qdot from y
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		int idxR = m_nodes[i]->idxR;
		m_nodes[i]->x = y.segment<3>(idxR);
		m_nodes[i]->v = y.segment<3>(nr + idxR);
	}
}

void SpringSerial::scatterDDofs_(VectorXd &ydot, int nr) {
	// Scatters qdot and qddot from ydot
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		int idxR = m_nodes[i]->idxR;
		m_nodes[i]->v = ydot.segment<3>(idxR);
		m_nodes[i]->a = ydot.segment<3>(nr + idxR);
	}


}

void SpringSerial::computeMassForce_(Vector3d grav, MatrixXd &M, VectorXd &f) {
	// Computes maximal mass matrix and force vector
	int n_nodes = m_nodes.size();
	double m = m_mass / n_nodes;

	Matrix3d I3 = Matrix3d::Identity();
	for (int i = 0; i < n_nodes; i++) {
		int idxM = m_nodes[i]->idxM;
		M.block<3, 3>(idxM, idxM) = m * I3;
		f.segment<3>(idxM) += m * grav;
	}

}
void SpringSerial::computeEnergies_(Vector3d grav, double &T, double &V) {
	int n_nodes = m_nodes.size();

	double m = m_mass / n_nodes;

	for (int i = 0; i < n_nodes; i++) {
		Vector3d x = m_nodes[i]->x;
		Vector3d v = m_nodes[i]->v;
		T = T + 0.5 * m * v.dot(v);
		V = V - m * grav.dot(x);
	}

	for (int i = 0; i < n_nodes - 1; i++) {
		Vector3d x0 = m_nodes[i]->x;
		Vector3d x1 = m_nodes[i+1]->x;
		Vector3d dx = x1 - x0;
		double l = sqrt(dx.dot(dx));
		double L = m_nodes[i]->L;
		double e = (l - L) / L;
		V = V + 0.5 * m_K * e * e;

	}

}
