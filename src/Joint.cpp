#include "Joint.h"
#include <iostream>

#include "Body.h"
#include "SE3.h"

using namespace std;
using namespace Eigen;

Joint::~Joint() {

}

void Joint::setJointTransform(Matrix4d E) {
	// Sets the transform of this joint wrt parent joint
	E_pj0 = E;
}

void Joint::update() {
	// Updates this joint and the attached body
	



}