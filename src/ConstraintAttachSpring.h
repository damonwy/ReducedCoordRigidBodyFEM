#pragma once
// ConstraintAttachSpring

#include "Constraint.h"

class SpringSerial;

class ConstraintAttachSpring : public Constraint
{
public:
	ConstraintAttachSpring();
	ConstraintAttachSpring(std::shared_ptr<SpringSerial> spring);
	void computeJacEqM_(Eigen::MatrixXd &Gm, Eigen::MatrixXd &Gmdot, Eigen::VectorXd &gm);


	std::shared_ptr<SpringSerial> m_spring;

protected:


};