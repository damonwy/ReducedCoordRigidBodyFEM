#pragma once
#define EIGEN_USE_MKL_ALL
#include "SoftBody.h"


class SoftBodyNull : public SoftBody {

public:

	SoftBodyNull();
	virtual ~SoftBodyNull() {}

	virtual Eigen::MatrixXd computeJacobian(Eigen::MatrixXd J);
	virtual Eigen::MatrixXd computeMass(Eigen::Vector3d grav, Eigen::MatrixXd M);
	virtual Eigen::VectorXd computeForce(Eigen::Vector3d grav, Eigen::VectorXd f);
	virtual Eigen::MatrixXd computeStiffness(Eigen::MatrixXd K);
	virtual Eigen::VectorXd gatherDofs(Eigen::VectorXd y, int nr);
	virtual Eigen::VectorXd gatherDDofs(Eigen::VectorXd ydot, int nr);
	virtual void scatterDofs(Eigen::VectorXd &y, int nr);
	virtual void scatterDDofs(Eigen::VectorXd &ydot, int nr);

protected:
	void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const {}
};