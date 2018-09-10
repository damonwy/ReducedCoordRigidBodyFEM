#include "JointFixed.h"

#include <iostream>

#include "Body.h"
#include "SE3.h"

using namespace std;
using namespace Eigen;

JointFixed::JointFixed() {

}


JointFixed::JointFixed(std::shared_ptr<Body> body, std::shared_ptr<Joint> parent = nullptr) :
	Joint(body, 0, parent)
{

}

JointFixed::~JointFixed() {

}

void JointFixed::draw() {

}

void JointFixed::update() {
	E_pj = E_pj0;
}