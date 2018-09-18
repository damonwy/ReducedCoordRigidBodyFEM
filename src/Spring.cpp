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

void Spring::countDofs() {
	countDofs_();
	if (next != nullptr) {
		next->countDofs();
	}
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

void Spring::computeJacobian(Eigen::MatrixXd &J, Eigen::MatrixXd &Jdot) {



}


void Spring::computeMassForce(Eigen::Vector3d grav, Eigen::MatrixXd &M, Eigen::VectorXd &f) {
	computeMassForce_(grav, M, f);
	if (next != nullptr) {
		next->computeMassForce(grav, M, f);
	}

}

void Spring::computeEnergies(Eigen::Vector3d grav, double &T, double &V) {
	computeEnergies_(grav, T, V);
	if (next != nullptr) {
		next->computeEnergies(grav, T, V);
	}
}

void Spring::countDofs_() {

}

void Spring::gatherDofs_(Eigen::VectorXd &y, int nr) {


}

void Spring::gatherDDofs_(Eigen::VectorXd &ydot, int nr) {


}

void Spring::scatterDofs_(Eigen::VectorXd &y, int nr) {


}

void Spring::scatterDDofs_(Eigen::VectorXd &ydot, int nr) {


}

void Spring::computeMassForce_(Eigen::Vector3d grav, Eigen::MatrixXd &M, Eigen::VectorXd &f) {



}

void Spring::computeEnergies_(Eigen::Vector3d grav, double &T, double &V) {



}


void Spring::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P) const {



}