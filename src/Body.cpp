#include "Body.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

#include "Wrench.h"
#include "Joint.h"
#include "SE3.h"
#include "Shape.h"
#include "MatrixStack.h"
#include "Program.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

Body::Body(double _density, Vector3d _sides):
density(_density),
sides(_sides)
{

}


Body::~Body() {
}

void Body::load(const string &RESOURCE_DIR) {
	// Inits shape
	boxShape = make_shared<Shape>();
	boxShape->loadMesh(RESOURCE_DIR + "box5.obj");

}

void Body::init(int &nm) {
	boxShape->init();
	computeInertia();
	countDofs(nm);

}

void Body::setTransform(Eigen::Matrix4d E) {
	// Sets the transform of this body wrt parent joint
	E_ji = E;
	E_ij = SE3::inverse(E_ji);
	Ad_ji = SE3::adjoint(E_ji);
	Ad_ij = SE3::adjoint(E_ij);
}

void Body::update() {
	// Updates transforms and maximal velocities
	E_wi = m_joint->E_wj * E_ji;
	E_iw = SE3::inverse(E_wi);
	Ad_wi = SE3::adjoint(E_wi);
	Ad_iw = SE3::adjoint(E_iw);

	// Joint velocity
	V = m_joint->m_S * m_joint->m_qdot;

	// Add parent velocity
	E_ip.setIdentity();
	if (m_joint->getParent() != nullptr) {
		m_parent = m_joint->getParent()->getBody();
		V = V + m_joint->Ad_jp * m_parent->V;
		E_ip = E_iw * m_parent->E_wi;
	}

	Ad_ip = SE3::adjoint(E_ip);
	phi = Ad_ij * V;
	Addot_wi = SE3::dAddt(E_wi, phi);

}

Energy Body::computeEnergies(Eigen::Vector3d grav, Energy energies) {
	// Compute kinetic and potential energies
	// Note: V is the twist at the parent joint, not at body
	// This is the same: K = 0.5 * phi.transpose()* diag(I_i) * phi;
	Energy ener;
	ener.K = energies.K + 0.5 * V.transpose() * I_j * V;
	ener.V = energies.V - I_j(5, 5)*grav.transpose() * E_wi.block<3, 1>(0, 3);

	if (next != nullptr) {
		ener = next->computeEnergies(grav, ener);
	}
	return ener;
}

void Body::computeInertiaBody() {
	I_i = SE3::inertiaCuboid(sides, density);
}

void Body::computeInertiaJoint() {
	double mass = I_i(3);
	Matrix3d R = E_ji.block<3, 3>(0, 0);
	Vector3d p = E_ji.block<3, 1>(0, 3);
	Matrix3d pBrac = SE3::bracket3(p);
	Matrix3d Ic = Matrix3d(I_i.segment<3>(0).asDiagonal());
	I_j.block<3, 3>(0, 0) = R * Ic *R.transpose() + mass *(pBrac.transpose()*pBrac);
	I_j.block<3, 3>(0, 3) = mass * pBrac;
	I_j.block<3, 3>(3, 0) = mass * pBrac.transpose();
	Matrix3d I3;
	I3.setIdentity();
	I_j.block<3, 3>(3, 3) = mass * I3;
}

void Body::computeInertia() {
	// Computes inertia at body and joint
	computeInertiaBody();
	computeInertiaJoint();
}

void Body::countDofs(int &nm) {
	// Counts maximal DOFs
	idxM = countM(nm, 6);
}

int Body::countM(int &nm, int data) {
	nm = nm + data;
	return (nm - data);
}

void Body::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, shared_ptr<MatrixStack> P) const
{
	prog->bind();
	if (boxShape) {
		glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
		MV->pushMatrix();
		glUniform3f(prog->getUniform("lightPos1"), 1.0, 1.0, 1.0);
		glUniform1f(prog->getUniform("intensity_1"), 0.8);
		glUniform3f(prog->getUniform("lightPos2"), -1.0, 1.0, 1.0);
		glUniform1f(prog->getUniform("intensity_2"), 0.2);
		glUniform1f(prog->getUniform("s"), 200);
		glUniform3f(prog->getUniform("ka"), 0.2, 0.2, 0.2);
		glUniform3f(prog->getUniform("kd"), 0.8, 0.7, 0.7);
		glUniform3f(prog->getUniform("ks"), 1.0, 0.9, 0.8);

		MV->pushMatrix();
		cout << E_wi << endl;
		MV->multMatrix(eigen_to_glm(E_wi));

		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		boxShape->draw(prog);
		MV->popMatrix();
	}
	prog->unbind();

}

MatrixXd Body::computeMass(Vector3d grav, MatrixXd M) {
	Matrix3d I;
	I.setZero();
	return I;

}

VectorXd Body::computeForce(Vector3d grav, MatrixXd f) {

	Vector3d I;
	I.setZero();
	return I;
}