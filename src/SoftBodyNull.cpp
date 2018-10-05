#include "SoftBodyNull.h"


using namespace std;
using namespace Eigen;

SoftBodyNull::SoftBodyNull():SoftBody() {



}

SoftBodyNull::~SoftBodyNull() {



}

void SoftBodyNull::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {

}

Eigen::MatrixXd SoftBodyNull::computeJacobian(Eigen::MatrixXd J) {

	return J;
}

Eigen::MatrixXd SoftBodyNull::computeMass(Eigen::Vector3d grav, Eigen::MatrixXd M) {
	return M;

}

Eigen::VectorXd SoftBodyNull::computeForce(Eigen::Vector3d grav, Eigen::VectorXd f) {

	return f;

}

Eigen::MatrixXd SoftBodyNull::computeStiffness(Eigen::MatrixXd K) {

	return K;
}

Eigen::VectorXd SoftBodyNull::gatherDofs(Eigen::VectorXd y, int nr) {

	return y;
}

Eigen::VectorXd SoftBodyNull::gatherDDofs(Eigen::VectorXd ydot, int nr) {

	return ydot;

}
void SoftBodyNull::scatterDofs(Eigen::VectorXd &y, int nr) {

}

void SoftBodyNull::scatterDDofs(Eigen::VectorXd &ydot, int nr) {


}