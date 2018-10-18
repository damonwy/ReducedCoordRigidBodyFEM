#include "Spring.h"

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

Spring::Spring():
	m_K(0.0), m_mass(1.0) {

}

Spring::Spring(int &countS, int &countCM):
	m_K(0.0), m_mass(1.0) {
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

VectorXd Spring::gatherDofs(VectorXd y, int nr) {
	VectorXd y_ = gatherDofs_(y, nr);
	if (next != nullptr) {
		y_ = next->gatherDofs(y_, nr);
	}
	return y_;
}

VectorXd Spring::gatherDDofs(VectorXd ydot, int nr) {
	VectorXd ydot_ = gatherDDofs_(ydot, nr);
	if (next != nullptr) {
		ydot_ = next->gatherDofs(ydot_, nr);
	}
	return ydot_;
}

void Spring::scatterDofs(VectorXd &y, int nr) {
	scatterDofs_(y, nr);
	if (next != nullptr) {
		next->scatterDofs(y, nr);
	}

}

void Spring::scatterDDofs(VectorXd &ydot, int nr) {

	scatterDDofs_(ydot, nr);
	if (next != nullptr) {
		next->scatterDDofs(ydot, nr);

	}
	
}

MatrixXd Spring::computeJacobian(MatrixXd J) {
	MatrixXd J_ = computeJacobian_(J);
	if (next != nullptr) {
		J_ = next->computeJacobian(J_);
	}
	return J_;
}

MatrixXd Spring::computeJacobian_(MatrixXd J) {
	return J;
}

MatrixXd Spring::computeMass(Vector3d grav, MatrixXd M) {
	MatrixXd M_ = computeMass_(grav, M);
	if (next != nullptr) {
		M_ = next->computeMass(grav, M_);
	}
	return M_;
}

VectorXd Spring::computeForce(Vector3d grav, VectorXd f) {
	VectorXd f_ = computeForce_(grav, f);
	if (next != nullptr) {
		f_ = next->computeForce(grav, f_);
	}
	return f_;
}

Energy Spring::computeEnergies(Vector3d grav, Energy ener) {
	ener = computeEnergies_(grav, ener);
	if (next != nullptr) {
		ener = next->computeEnergies(grav, ener);
	}
	return ener;
}

VectorXd Spring::gatherDofs_(VectorXd y, int nr) {
	return y;
}

VectorXd Spring::gatherDDofs_(VectorXd ydot, int nr) {
	return ydot;
}

void Spring::scatterDofs_(VectorXd &y, int nr) {

}

void Spring::scatterDDofs_(VectorXd &ydot, int nr) {

}

MatrixXd Spring::computeMass_(Vector3d grav, MatrixXd M) {
	return M;
}

VectorXd Spring::computeForce_(Vector3d grav, VectorXd f) {
	return f;
}

Energy Spring::computeEnergies_(Vector3d grav, Energy ener) {
	return ener;
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