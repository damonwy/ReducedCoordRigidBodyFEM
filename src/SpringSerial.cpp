#include "SpringSerial.h"

#include "Body.h"
#include "Node.h"
#include "MatrixStack.h"
#include "Program.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

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

void SpringSerial::init() {
	// Sets the world positions of the nodes using the attachment
	// points. r0 and r1 are in local coords.
	Matrix4d E0, E1;
	if (m_body0 == nullptr) {
		E0 = Matrix4d::Identity();
	}
	else {
		E0 = m_body0->E_wi;
	}
	
	if (m_body1 == nullptr) {
		E1 = Matrix4d::Identity();
	}
	else {
		E1 = m_body1->E_wi;
	}
	Vector4d r0, r1;
	r0.segment<3>(0) = m_r0;
	r1.segment<3>(0) = m_r1;
	r0(3) = 1.0; 
	r1(3) = 1.0;
	Vector4d x0 = E0 * r0;
	Vector4d x1 = E1 * r1;

	// Set the nodal positions
	int n_nodes = m_nodes.size();

	for (int i = 0; i < n_nodes; i++) {
		double s = i / (n_nodes - 1);
		m_nodes[i]->x = (1 - s) * x0 + s * x1;
	}

	// Compute the rest lengths
	for (int i = 0; i < n_nodes - 1; i++) {
		Vector3d x0 = m_nodes[i]->x;
		Vector3d x1 = m_nodes[i + 1]->x;
		Vector3d dx = x1 - x0;
		m_nodes[i]->L = dx.norm(); // rest length
	}
}

void SpringSerial::load(const string &RESOURCE_DIR) {


}

void SpringSerial::draw_(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
	// Draw nodes
	prog->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	MV->pushMatrix();
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	int n_nodes = m_nodes.size();
	for (int i = 0; i < n_nodes; i++) {	
		m_nodes[i]->draw(MV, prog);
	}
	MV->popMatrix();
	prog->unbind();

	// Draw line segments
	progSimple->bind();
	glLineWidth(5);
	glBegin(GL_LINES);
	glColor3f(1.0, 1.0, 1.0);
	
	for (int i = 0; i < n_nodes - 1; i++) {
		Vector3d x0 = m_nodes[i]->x;
		Vector3d x1 = m_nodes[i + 1]->x;
		glVertex3f(x0(0), x0(1), x0(2));
		glVertex3f(x1(0), x1(1), x1(2));
	}
	glEnd();
	progSimple->unbind();
	
}

void SpringSerial::setAttachments(shared_ptr<Body> body0, Vector3d r0, shared_ptr<Body> body1, Vector3d r1) {
	// Attaches this spring to body0 and body1
	m_body0 = body0;
	m_body1 = body1;
	m_r0 = r0;
	m_r1 = r1;
}

void SpringSerial::countDofs_(int &nm, int &nr) {
	// Counts maximal and reduced DOFs
	// For non-rigid DOFs, we need both maximal and reduced DOFs,
	// and the Jacobian must pass them through with the identity 
	// matrices

	for (int i = 0; i < m_nodes.size(); i++) {
		m_nodes[i]->idxM = nm;
		m_nodes[i]->idxR = nr;
		nm += 3;
		nr += 3;
	}

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

	for (int i = 0; i < n_nodes - 1; i++) {
		int row0 = m_nodes[i]->idxM;
		int row1 = m_nodes[i + 1]->idxM;
		Vector3d x0 = m_nodes[i]->x;
		Vector3d x1 = m_nodes[i + 1]->x;
		Vector3d dx = x1 - x0;
		double l = dx.norm();
		double L = m_nodes[i]->L;
		double e = (l - L) / L;
		Vector3d fs = m_K * e * (1.0 / L)* dx / l;
		f.segment<3>(row0) += fs;
		f.segment<3>(row1) -= fs;
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
		double l = dx.norm();
		double L = m_nodes[i]->L;
		double e = (l - L) / L;
		V = V + 0.5 * m_K * e * e;

	}

}

void SpringSerial::computeJacobian_(MatrixXd &J, MatrixXd &Jdot) {
	for (int i = 0; i < m_nodes.size(); i++) {
		J.block<3, 3>(m_nodes[i]->idxM, m_nodes[i]->idxR) = Matrix3d::Identity();
	}
}

SpringSerial:: ~SpringSerial() {

}