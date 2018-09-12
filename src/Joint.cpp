#include "Joint.h"
#include <iostream>

#include "Body.h"
#include "SE3.h"
#include "MatrixStack.h"
#include "Program.h"

using namespace std;
using namespace Eigen;

Joint::Joint() {

}

Joint::Joint(shared_ptr<Body> body, double ndof, shared_ptr<Joint> parent) :
m_body(body),
m_parent(parent),
m_ndof(ndof)
{
	if (parent == nullptr) {
		m_name = "NULL-" + body->getName();
	}
	else {
		m_name = parent->getBody()->getName() + "-" + body->getName();
	}

	m_q.resize(m_ndof);
	m_qdot.resize(m_ndof);
	m_qddot.resize(m_ndof);

	m_q.setZero();
	m_qdot.setZero();
	m_qddot.setZero();

	m_tau.resize(m_ndof);
	m_tau.setZero();

	m_tauCon.resize(m_ndof);
	m_tauCon.setZero();

	m_K = 0.0;

	m_S.resize(6, ndof);
	m_S.setZero();
	m_Sdot.resize(6, ndof);
	m_Sdot.setZero();
	
	

}

Joint::~Joint() {

}

void Joint::init(int &nr) {
	m_body->setJoint(getJoint());

	if (m_parent != nullptr) {
		m_parent->addChild(getJoint());
	}
	countDofs(nr);

}

void Joint::setJointTransform(Matrix4d E) {
	// Sets the transform of this joint wrt parent joint
	E_pj0 = E;
}

void Joint::update() {
	// Updates this joint and the attached body
	updateSelf();
	E_jp = SE3::inverse(E_pj);
	Ad_jp = SE3::adjoint(E_jp);

	Matrix4d E_wp;

	if (m_parent == nullptr) {
		E_wp.setIdentity();
	}
	else {
		E_wp = m_parent->E_wj;
	}

	E_wj = E_wp * E_pj;

	// Update attached body
	m_body->update();
	if (next != nullptr) {
		next->update();
	}

}

void Joint::updateSelf() {

}

void Joint::countDofs(int &nr) {
	// Counts reduced DOFs
	idxR = countR(nr, m_ndof);

}

int Joint::countR(int &nr, int data) {
	nr = nr + data;
	return (nr - data);
}

MatrixXd Joint::computeJacobian(MatrixXd J, int nm, int nr) {
	// Computes the redmax Jacobian
	




	return J;
}

MatrixXd Joint::computeJacobianDerivative(MatrixXd Jdot, int nm, int nr) {

	




	return Jdot;
}

VectorXd Joint::computerJacTransProd(VectorXd y, VectorXd x, int nr) {
	// Computes x = J'*y
	// x (nr, 1)
	VectorXd yi = y.segment(m_body->idxM, 6);
	for (int k = 0; k < (int)m_children.size(); k++) {
		yi = yi + m_children[k]->getAlpha();
	}
	m_alpha = m_body->Ad_ip.transpose() * yi;
	if (prev != nullptr) {
		x = prev->computerJacTransProd(y, x, nr);
	}
	return x;
}


Energy Joint::computeEnergies(Vector3d grav, Energy ener) {
	// Computes kinetic and potential energies
	ener = m_body->computeEnergies(grav, ener);
	ener.V += 0.5 * m_K * m_q.dot(m_q);
	if (next != nullptr) {
		ener = next->computeEnergies(grav, ener);
	}
	return ener;
}

Eigen::VectorXd Joint::gatherDofs(Eigen::VectorXd y, int nr) {
	// Gathers q and qdot into y
	y.segment(idxR, m_ndof) = m_q;
	y.segment(nr + idxR, m_ndof) = m_qdot;
	if (next != nullptr) {
		y = next->gatherDofs(y, nr);
	}
	return y;
}

Eigen::VectorXd Joint::gatherDDofs(Eigen::VectorXd ydot, int nr) {
	// Gathers qdot and qddot into ydot
	ydot.segment(idxR, m_ndof) = m_qdot;
	ydot.segment(nr + idxR, m_ndof) = m_qddot;
	if (next != nullptr) {
		ydot = next->gatherDDofs(ydot, nr);
	}
	return ydot;
}

void Joint::scatterDofs(Eigen::VectorXd y, int nr) {
	// Scatters q and qdot from y
	scatterDofsNoUpdate(y, nr);
	update();
}

void Joint::scatterDDofs(Eigen::VectorXd ydot, int nr) {
	// Scatters qdot and qddot from ydot
	m_qdot.segment(0, m_ndof) = ydot.segment(idxR, m_ndof);
	m_qddot.segment(0, m_ndof) = ydot.segment(nr + idxR, m_ndof);
	if (next != nullptr) {
		next->scatterDDofs(ydot, nr);
	}
}

void Joint::scatterDofsNoUpdate(Eigen::VectorXd y, int nr) {
	// Helper function to scatter without updating
	m_q.segment(0, m_ndof) = y.segment(idxR, m_ndof);
	m_qdot.segment(0, m_ndof) = y.segment(nr + idxR, m_ndof);
	if (next != nullptr) {
		next->scatterDofsNoUpdate(y, nr);
	}
}

void Joint::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, shared_ptr<MatrixStack> P) const {

	prog->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	MV->pushMatrix();
	MV->multMatrix(eigen_to_glm(E_wj));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	
	glLineWidth(2);
	glBegin(GL_LINES);
	
	// X axis
	glColor3f(1.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(1.0, 0.0, 0.0);

	// Y axis
	glColor3f(0.0, 1.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 1.0, 0.0);

	// Z axis
	glColor3f(0.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 1.0);

	glEnd();

	MV->popMatrix();
	prog->unbind();
}