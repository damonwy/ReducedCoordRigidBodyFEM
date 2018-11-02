#include "Joint.h"
#include <iostream>

#include "Body.h"
#include "SE3.h"
#include "MatrixStack.h"
#include "Program.h"
#include "Shape.h"

using namespace std;
using namespace Eigen;

Joint::Joint() {
	presc = false;
}

Joint::Joint(shared_ptr<Body> body, int ndof, shared_ptr<Joint> parent) :
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
	m_D = 0.0;

	
	m_S.resize(6, ndof);
	m_S.setZero();
	m_Sdot.resize(6, ndof);
	m_Sdot.setZero();
	m_Q.setIdentity();
	
	presc = false;

}

Joint::~Joint() {

}

void Joint::load(const std::string &RESOURCE_DIR, std::string joint_shape) {

	m_jointShape = make_shared<Shape>();
	m_jointShape->loadMesh(RESOURCE_DIR + joint_shape);

}

void Joint::init(int &nm, int &nr) {
	if (m_jointShape) {
		m_jointShape->init();
	}

	m_body->setJoint(getJoint());

	if (m_parent != nullptr) {
		m_parent->addChild(getJoint());
	}
	countDofs(nm, nr);

}

void Joint::setJointTransform(Matrix4d E) {
	// Sets the transform of this joint wrt parent joint
	E_pj0 = E;
}

void Joint::update() {
	// Updates this joint and the attached body
	updateSelf();
	E_pj = E_pj0 * m_Q;

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

void Joint::countDofs(int &nm, int &nr) {
	// Counts reduced DOFs
	idxR = countR(nr, m_ndof);
	m_body->countDofs(nm);
}

int Joint::countR(int &nr, int data) {
	nr = nr + data;
	return (nr - data);
}

void Joint::computeJacobian(MatrixXd &J, int nm, int nr) {
	// Computes the redmax Jacobian
	Matrix6d Ad_ij = m_body->Ad_ij;
	J.block(m_body->idxM, idxR, 6, m_ndof) = Ad_ij * m_S;
	// Loop through all ancestors
	auto jointA = m_parent;
	while (jointA != nullptr) {
		int idxM_P = m_parent->getBody()->idxM;
		Matrix6d Ad_ip = m_body->Ad_ip;
		J.block(m_body->idxM, jointA->idxR, 6, jointA->m_ndof) = Ad_ip * J.block(idxM_P, jointA->idxR, 6, jointA->m_ndof);
		jointA = jointA->getParent();
	}
	if (next != nullptr) {
		next->computeJacobian(J, nm, nr);
	}
}

void Joint::computeForceStiffness(VectorXd &fr) {
	// Computes joint stiffness force vector
	if (presc == false) {
		int row = this->idxR;
		// Add the joint torque here rather than having a separate function
		fr.segment(row, m_ndof) += m_tau - m_K * m_q;
	}

	if (next != nullptr) {
		next->computeForceStiffness(fr);
	}
}

void Joint::computeMatrixStiffness(MatrixXd &Ksr) {
	// Computes joint stiffness matrix
	if (presc == false) {
		int row = this->idxR;
		MatrixXd I(m_ndof, m_ndof);
		I.setIdentity();
		Ksr.block(row, row, m_ndof, m_ndof) -= m_K * I;
	}

	if (next != nullptr) {
		next->computeMatrixStiffness(Ksr);
	}
}

void Joint::computeForceDamping(VectorXd &fr) {
	// Computes joint damping force vector
	if (presc == false) {
		int row = this->idxR;
		fr.segment(row, m_ndof) -= m_D * m_qdot;
	}

	if (next != nullptr) {
		next->computeForceDamping(fr);
	}
}

void Joint::computeMatrixDamping(MatrixXd &Ddr) {
	// Computes joint damping matrix
	if (presc == false) {
		int row = this->idxR;
		MatrixXd I(m_ndof, m_ndof);
		I.setIdentity();
		Ddr.block(row, row, m_ndof, m_ndof) += m_D * I;
	}

	if (next != nullptr) {
		next->computeMatrixDamping(Ddr);
	}
}

void Joint::computeJacobianDerivative(MatrixXd &Jdot, MatrixXd J, int nm, int nr) {
	Matrix6d Ad_ij = m_body->Ad_ij;
	Jdot.block(m_body->idxM, idxR, 6, m_ndof) = Ad_ij * m_Sdot;
	// Loop through all ancestors
	auto jointA = m_parent;

	while (jointA != nullptr) {
		int idxM_P = m_parent->getBody()->idxM;
		Matrix6d Ad_ip = m_body->Ad_ip;
		Matrix6d Ad_iw = m_body->Ad_iw;
		Matrix6d Ad_wp = m_parent->getBody()->Ad_wi;
		Matrix6d Addot_wi = m_body->Addot_wi;
		Matrix6d Addot_wp = m_parent->getBody()->Addot_wi;

		Matrix6d Addot_ip = -Ad_iw * (Addot_wi * Ad_iw * Ad_wp - Addot_wp);
		Jdot.block(m_body->idxM, jointA->idxR, 6, jointA->m_ndof) = Ad_ip * Jdot.block(idxM_P, jointA->idxR, 6, jointA->m_ndof) + Addot_ip * J.block(idxM_P, jointA->idxR, 6, jointA->m_ndof);
		jointA = jointA->getParent();
	}
	if (next != nullptr) {
		next->computeJacobianDerivative(Jdot, J, nm, nr);
	}
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

void Joint::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {

	progSimple->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	MV->pushMatrix();
	MV->multMatrix(eigen_to_glm(E_wj));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	
	glLineWidth(5);
	glBegin(GL_LINES);
	// X axis
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(3.0f, 0.0f, 0.0f);

	// Y axis
	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 3.0f, 0.0f);

	// Z axis
	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 3.0f);

	glEnd();
	MV->popMatrix();
	progSimple->unbind();

	drawSelf(MV, prog, progSimple, P);
}

void Joint::drawSelf(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
	prog->bind();
	if (m_jointShape) {
		glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
		glUniform3f(prog->getUniform("lightPos1"), 66.0f, 25.0f, 25.0f);
		glUniform1f(prog->getUniform("intensity_1"), 0.6f);
		glUniform3f(prog->getUniform("lightPos2"), -66.0f, 25.0f, 25.0f);
		glUniform1f(prog->getUniform("intensity_2"), 0.2f);
		glUniform1f(prog->getUniform("s"), 300.0f);
		glUniform3f(prog->getUniform("ka"), 0.2f, 0.2f, 0.2f);
		glUniform3f(prog->getUniform("kd"), 0.8f, 0.7f, 0.7f);
		glUniform3f(prog->getUniform("ks"), 1.0f, 0.9f, 0.8f);

		MV->pushMatrix();
		MV->multMatrix(eigen_to_glm(E_wj));

		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		m_jointShape->draw(prog);
		MV->popMatrix();
	}
	prog->unbind();
}