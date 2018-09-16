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

Spring::~Spring() {



}

void Spring::countDofs() {



}


void Spring::gatherDofs(double &y) {

}

void Spring::gatherDDofs(double &ydot) {



}

void Spring::computeJacobian(Eigen::MatrixXd &J, Eigen::MatrixXd &Jdot) {



}


void Spring::computeMassForce(Eigen::Vector3d grav, Eigen::MatrixXd M, Eigen::VectorXd f) {
	


}

void Spring::computeEnergies(Eigen::Vector3d grav, double T, double V) {



}


void Spring::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P) const {



}