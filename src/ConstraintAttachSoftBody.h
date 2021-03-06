#pragma once
// ConstraintAttachSoftBody
#define EIGEN_USE_MKL_ALL

#include "Constraint.h"

class SoftBody;

class ConstraintAttachSoftBody: public Constraint
{
public:
	ConstraintAttachSoftBody();
	ConstraintAttachSoftBody(std::shared_ptr<SoftBody> softbody);
	std::shared_ptr<SoftBody> m_softbody;

protected:

	int n_attachments;
	int n_sliding_nodes;
	void computeJacEqMSparse_(std::vector<T> &Gm, std::vector<T> &Gmdot, Eigen::VectorXd &gm, Eigen::VectorXd &gmdot, Eigen::VectorXd &gmddot);
	void computeJacEqM_(Eigen::MatrixXd &Gm, Eigen::MatrixXd &Gmdot, Eigen::VectorXd &gm, Eigen::VectorXd &gmdot, Eigen::VectorXd &gmddot);

};