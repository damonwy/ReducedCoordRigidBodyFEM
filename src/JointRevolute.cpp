#include "JointRevolute.h"

#include <iostream>

#include "Body.h"
#include "SE3.h"

using namespace std;
using namespace Eigen;

JointRevolute::JointRevolute() {

}


JointRevolute::JointRevolute(std::shared_ptr<Body> body, Eigen::Vector3d axis, std::shared_ptr<Joint> parent = nullptr):
Joint(body, 1, parent),
m_axis(axis)
{

}

JointRevolute::~JointRevolute() {

}

void JointRevolute::draw() {

}

void JointRevolute::update() {

}