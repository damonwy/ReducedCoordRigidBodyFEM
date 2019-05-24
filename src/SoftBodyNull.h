#pragma once
#include "SoftBody.h"

class SoftBodyNull : public SoftBody {

public:
    SoftBodyNull(){}
	virtual ~SoftBodyNull() {}

    void computeJacobian(Eigen::MatrixXd &J){}
	void computeMass(Eigen::MatrixXd &M){}
    void computeForce(Eigen::Vector3d grav, Eigen::VectorXd &f){}
    void computeStiffness(Eigen::MatrixXd &K){}
    void gatherDofs(Eigen::VectorXd &y, int nr){}
    void gatherDDofs(Eigen::VectorXd &ydot, int nr){}
    void scatterDofs(Eigen::VectorXd &y, int nr){}
    void scatterDDofs(Eigen::VectorXd &ydot, int nr){}

protected:
	void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const {}
};
