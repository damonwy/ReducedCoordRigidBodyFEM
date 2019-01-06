#pragma once
#include "Constraint.h"

class Body;

class ConstraintPrescBody : public Constraint, public std::enable_shared_from_this<ConstraintPrescBody> {
public:
	ConstraintPrescBody() {}
	ConstraintPrescBody(std::shared_ptr<Body> body, Eigen::VectorXi prows, Integrator vel);	
	std::shared_ptr<ConstraintPrescBody> getSelf() { return shared_from_this(); }
	virtual ~ConstraintPrescBody() {}
	
protected:
	void computeJacEqR_(Eigen::MatrixXd &Gr, Eigen::MatrixXd &Grdot, Eigen::VectorXd &gr, Eigen::VectorXd &grdot, Eigen::VectorXd &grddot);
	void computeJacEqRSparse_(std::vector<T> &Gr, std::vector<T> &Grdot, Eigen::VectorXd &gr, Eigen::VectorXd &grdot, Eigen::VectorXd &grddot);
	void scatterForceEqM_();
	void init_();

private:
	std::shared_ptr<Body> m_body;
	Eigen::VectorXd m_q;
	Eigen::VectorXd m_qdot;
	Eigen::VectorXd m_qddot;
	Integrator m_vel;
	Eigen::VectorXi m_prows;
};