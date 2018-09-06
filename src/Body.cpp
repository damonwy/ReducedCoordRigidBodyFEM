#include "Body.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

#include "Wrench.h"
#include "Joint.h"
#include "Rigid.h"

Body::Body() {


}


Body::~Body() {



}

void Body::setTransform(Eigen::Matrix4d E) {
	// Sets the transform of this body wrt parent joint
	E_ji = E;
	E_ij = Rigid::inverse(E_ji);
	Ad_ji = Rigid::adjoint(E_ji);
	Ad_ij = Rigid::adjoint(E_ij);
}

void Body::update() {
	// Updates transforms and maximal velocities
	E_wi = joint->E_wj * E_ji;
	E_iw = Rigid::inverse(E_wi);
	Ad_wi = Rigid::adjoint(E_wi);
	Ad_iw = Rigid::adjoint(E_iw);

	// Joint velocity
	V = joint->S * joint->qdot;

	// Add parent velocity
	E_ip.setIdentity();
	if (joint->parent != nullptr) {
		parent = joint->getParent->body;
		V = V + joint->Ad_jp * parent->V;
		E_ip = E_iw * parent->E_wi;
	}

	Ad_ip = Rigid::adjoint(E_ip);
	phi = Ad_ij * V;
	Addot_wi = Rigid::dAddt(E_wi, phi);

}