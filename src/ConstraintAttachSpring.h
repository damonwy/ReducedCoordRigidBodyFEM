#pragma once
// ConstraintAttachSpring

#include "Constraint.h"

class SpringSerial;
class Spring;

class ConstraintAttachSpring : public Constraint
{
public:
	ConstraintAttachSpring();
	ConstraintAttachSpring(std::shared_ptr<Spring> spring);
	void computeJacEqM_(Eigen::MatrixXd &Gm, Eigen::MatrixXd &Gmdot, Eigen::VectorXd &gm, Eigen::VectorXd &gmdot, Eigen::VectorXd &gmddot);


	std::shared_ptr<Spring> m_spring;

protected:


};