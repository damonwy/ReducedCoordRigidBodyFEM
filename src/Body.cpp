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

Body::Body() {
	
}

Body::Body(double _density, Vector3d _sides):
density(_density),
sides(_sides)
{
	wext_i.setZero();
	m_attached_color << (float)(rand() % 255)/255.0,(float)(rand() % 255)/255.0,(float)(rand() % 255)/255.0;
	m_sliding_color << (float)(rand() % 255) / 255.0, (float)(rand() % 255) / 255.0, (float)(rand() % 255) / 255.0;
}


Body::~Body() {
}

void Body::load(const string &RESOURCE_DIR, string box_shape) {
	//read a JSON file
	//ifstream i(RESOURCE_DIR + "input.json");
	//json js;
	//i >> js;
	//i.close();
	//string box_shape = js[box_shape];

	// Inits shape
	boxShape = make_shared<Shape>();
	boxShape->loadMesh(RESOURCE_DIR + box_shape);

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
	I_j.block<3, 3>(3, 3) = mass * Matrix3d::Identity();
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
		glUniform3f(prog->getUniform("lightPos1"), 66.0, 25.0, 25.0);
		glUniform1f(prog->getUniform("intensity_1"), 0.6);
		glUniform3f(prog->getUniform("lightPos2"), -66.0, 25.0, 25.0);
		glUniform1f(prog->getUniform("intensity_2"), 0.2);
		glUniform1f(prog->getUniform("s"), 300);
		glUniform3f(prog->getUniform("ka"), 0.2, 0.2, 0.2);
		glUniform3f(prog->getUniform("kd"), 0.8, 0.7, 0.7);
		glUniform3f(prog->getUniform("ks"), 1.0, 0.9, 0.8);

		MV->pushMatrix();
		MV->multMatrix(eigen_to_glm(E_wi));

		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		boxShape->draw(prog);
		MV->popMatrix();

	}
	prog->unbind();

}

MatrixXd Body::computeMass(Vector3d grav, MatrixXd M) {
	// Computes maximal mass matrix 
	Matrix6d M_i = Matrix6d(I_i.asDiagonal());
	
	M.block<6, 6>(idxM, idxM) = M_i;
	if (next != nullptr) {
		M = next->computeMass(grav, M);
	}

	return M;
}

VectorXd Body::computeForce(Vector3d grav, VectorXd f) {
	// Computes maximal force vector
	Matrix6d M_i = Matrix6d(I_i.asDiagonal());
	Vector6d fcor = SE3::ad(phi).transpose() * M_i * phi;
	Matrix3d R_wi = E_wi.block<3, 3>(0, 0);
	Matrix3d R_iw = R_wi.transpose();

	Vector6d fgrav;
	fgrav.setZero();

	fgrav.segment<3>(3) = M_i(3, 3) * R_iw * grav; // wrench in body space
	f.segment<6>(idxM) = fcor + fgrav + wext_i;
	
	// Joint torque

	//if (!m_joint->presc) {
	//	
	//	Vector6d tau = m_joint->m_S * (m_joint->m_tau - m_joint->m_K * m_joint->m_q);
	//	f.segment<6>(idxM) += Ad_ji.transpose() * tau;
	//	// Also apply to parent
	//	if (m_joint->getParent() != nullptr) {
	//		m_parent = m_joint->getParent()->getBody();
	//		int idxM_P = m_parent->idxM;
	//		Matrix4d E_jp = E_ji * E_iw * m_parent->E_wi; // this joint -> parent body
	//		f.segment<6>(idxM_P) -= SE3::adjoint(E_jp).transpose() * tau;
	//	}
	//}

	if (next != nullptr) {
		f = next->computeForce(grav, f);
	}

	return f;
}