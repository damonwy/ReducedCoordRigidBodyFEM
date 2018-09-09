#include "Body.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

#include "Wrench.h"
#include "Joint.h"
#include "SE3.h"


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

void Body::setTransform(Eigen::Matrix4d E) {
	// Sets the transform of this body wrt parent joint
	E_ji = E;
	E_ij = SE3::inverse(E_ji);
	Ad_ji = SE3::adjoint(E_ji);
	Ad_ij = SE3::adjoint(E_ij);
}

void Body::update() {
	// Updates transforms and maximal velocities
	E_wi = joint->E_wj * E_ji;
	E_iw = SE3::inverse(E_wi);
	Ad_wi = SE3::adjoint(E_wi);
	Ad_iw = SE3::adjoint(E_iw);

	// Joint velocity
	V = joint->S * joint->qdot;

	// Add parent velocity
	E_ip.setIdentity();
	if (joint->parent != nullptr) {
		parent = joint->getParent->body;
		V = V + joint->Ad_jp * parent->V;
		E_ip = E_iw * parent->E_wi;
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
