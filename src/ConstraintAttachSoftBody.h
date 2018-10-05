#pragma once
// ConstraintAttachSoftBody

#include "Constraint.h"

class SoftBody;

class ConstraintAttachSoftBody: public Constraint
{
public:
	ConstraintAttachSoftBody();
	ConstraintAttachSoftBody(std::shared_ptr<SoftBody> softbody);
	void computeJacEqM_(Eigen::MatrixXd &Gm, Eigen::MatrixXd &Gmdot, Eigen::VectorXd &gm, Eigen::VectorXd &gmdot, Eigen::VectorXd &gmddot);


	std::shared_ptr<SoftBody> m_softbody;

protected:

	int n_attachments;
};