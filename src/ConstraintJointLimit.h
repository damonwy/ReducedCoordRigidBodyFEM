#pragma once
// ConstraintJointLimit Inequality joint limit constraint
// Works only for revolute joints


#include "Constraint.h"

class Joint;

class ConstraintJointLimit : public Constraint 
{
public:
	ConstraintJointLimit();
	ConstraintJointLimit(std::shared_ptr<Joint> joint);
	void setLimits(double ql, double qu) { m_ql = ql; m_qu = qu; }
	void computeJacIneqR_(Eigen::MatrixXd &Cr, Eigen::MatrixXd &Crdot, Eigen::VectorXd &cr);

	std::shared_ptr<Joint> m_joint;	// Joint to apply the constriant on
	double m_ql;						// lower limit
	double m_qu;						// upper limit
protected:
	

};