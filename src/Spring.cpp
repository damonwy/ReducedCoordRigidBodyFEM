#include "Spring.h"

#include "Deformable.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

#include "MatrixStack.h"
#include "Program.h"
#include "Node.h"
#include "Body.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

Spring::Spring(){

}

void Spring::init() {
	init_();
	if (next != nullptr) {
		next->init();
	}
}

void Spring::computeForceStiffnessDamping(Eigen::Vector3d grav, Eigen::Vector3d &f, Eigen::MatrixXd K, Eigen::MatrixXd D) {
	computeForceStiffnessDamping_(grav, f, K, D);

	if (next != nullptr) {
		next->computeForceStiffnessDamping(grav, f, K, D);
	}

}

void Spring::computeStiffnessProd(Eigen::VectorXd x, Eigen::VectorXd &y) {
	// Computes y=K*x
	computeStiffnessProd_(x, y);
	if (next != nullptr) {
		next->computeStiffnessProd(x, y);
	}
}

void Spring::computeDampingProd(Eigen::VectorXd x, Eigen::VectorXd &y) {
	// Computes y=D*x
	computeDampingProd_(x, y);
	if (next != nullptr) {
		next->computeDampingProd(x, y);
	}
}

void Spring::computeEnergies(Vector3d grav, Energy &ener) {
	computeEnergies_(grav, ener);
	if (next != nullptr) {
		next->computeEnergies(grav, ener);
	}
}

void Spring::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
	draw_(MV, prog, progSimple, P);
	if (next != nullptr)
	{
		next->draw(MV, prog, progSimple, P);
	}
}