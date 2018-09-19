#include "Spring.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

#include "MatrixStack.h"
#include "Program.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

Spring::Spring() {

}

Spring::Spring(int &countS, int &countCM) {
	countS++;
	countCM++;
}

Spring::~Spring() {
}

void Spring::load(const string &RESOURCE_DIR) {

}

void Spring::init() {

}

void Spring::countDofs(int &nm, int &nr) {
	countDofs_(nm, nr);
	if (next != nullptr) {
		next->countDofs(nm, nr);
	}
}

void Spring::countDofs_(int &nm, int &nr) {

}

void Spring::gatherDofs(VectorXd &y, int nr) {
	gatherDofs_(y, nr);
	if (next != nullptr) {
		next->gatherDofs(y, nr);
	}
}

void Spring::gatherDDofs(VectorXd &ydot, int nr) {
	gatherDDofs_(ydot, nr);
	if (next != nullptr) {
		next->gatherDofs(ydot, nr);
	}
}

void Spring::computeJacobian(MatrixXd &J, MatrixXd &Jdot) {
	computeJacobian_(J, Jdot);
	if (next != nullptr) {
		next->computeJacobian(J, Jdot);
	}
}

void Spring::computeJacobian_(MatrixXd &J, MatrixXd &Jdot) {

}

void Spring::computeMassForce(Vector3d grav, MatrixXd &M, VectorXd &f) {
	computeMassForce_(grav, M, f);
	if (next != nullptr) {
		next->computeMassForce(grav, M, f);
	}
}

void Spring::computeEnergies(Vector3d grav, double &T, double &V) {
	computeEnergies_(grav, T, V);
	if (next != nullptr) {
		next->computeEnergies(grav, T, V);
	}
}

void Spring::countDofs_() {

}

void Spring::gatherDofs_(VectorXd &y, int nr) {

}

void Spring::gatherDDofs_(VectorXd &ydot, int nr) {

}

void Spring::scatterDofs_(VectorXd &y, int nr) {

}

void Spring::scatterDDofs_(VectorXd &ydot, int nr) {

}

void Spring::computeMassForce_(Vector3d grav, MatrixXd &M, VectorXd &f) {

}

void Spring::computeEnergies_(Vector3d grav, double &T, double &V) {

}


void Spring::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
	draw_(MV, prog, progSimple, P);
	if (next != nullptr)
	{
		next->draw(MV, prog, progSimple, P);
	}
}

void Spring::draw_(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {

}