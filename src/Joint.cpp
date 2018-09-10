#include "Joint.h"
#include <iostream>

#include "Body.h"
#include "SE3.h"

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
	



}

void Joint::countDofs(int &nr) {
	// Counts reduced DOFs
	idxR = countR(nr, m_ndof);

}

int Joint::countR(int &nr, int data) {
	nr = nr + data;
	return (nr - data);
}