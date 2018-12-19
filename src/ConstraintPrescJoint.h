#pragma once
#include "Constraint.h"

class Joint;

class ConstraintPrescJoint : public Constraint {
public:
	ConstraintPrescJoint() {}
	ConstraintPrescJoint(std::shared_ptr<Joint> joint, Integrator vel);
	virtual ~ConstraintPrescJoint() {}

protected:
	void computeJacEqR_(Eigen::MatrixXd &Gr, Eigen::MatrixXd &Grdot, Eigen::VectorXd &gr, Eigen::VectorXd &grdot, Eigen::VectorXd &grddot);
	void computeJacEqRSparse_(std::vector<T> &Gr, std::vector<T> &Grdot, Eigen::VectorXd &gr, Eigen::VectorXd &grdot, Eigen::VectorXd &grddot);
	void scatterForceEqR_();



private:
	std::shared_ptr<Joint> m_joint;
	Eigen::VectorXd m_q;
	Eigen::VectorXd m_qdot;
	Eigen::VectorXd m_qddot;
	Integrator m_vel;
};