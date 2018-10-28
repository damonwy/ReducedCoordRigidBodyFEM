#include "JointRevolute.h"

#include <iostream>

#include "Body.h"
#include "SE3.h"

using namespace std;
using namespace Eigen;

JointRevolute::JointRevolute() {

}


JointRevolute::JointRevolute(std::shared_ptr<Body> body, Eigen::Vector3d axis, std::shared_ptr<Joint> parent):
Joint(body, 1, parent)
{
	m_axis = axis;
}

JointRevolute::~JointRevolute() {

}

void JointRevolute::updateSelf() {
	Matrix3d R = SE3::aaToMat(m_axis, m_q(0));
	Matrix4d Q;
	Q.setIdentity();
	Q.block<3, 3>(0, 0) = R;
	m_Q = Q;
	//E_pj = E_pj0 * Q;
	
	m_S.block<3, 1>(0, 0) = m_axis;
}