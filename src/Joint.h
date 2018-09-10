// Joint A generic joint between two rigid bodies
//    A joint is defined between a parent and the child. The DOF, q, of the joint is 
//    the relative displacement of the child wrt the parent

#pragma once
#ifndef MUSCLEMASS_SRC_JOINT_H_
#define MUSCLEMASS_SRC_JOINT_H_
#include <vector>
#include <memory>
#define EIGEN_DONT_ALIGN_STATICALLY

#include <Eigen/Dense>
#include <iostream>
#include "MLCommon.h"

class SE3;
class Body;

class Joint :public std::enable_shared_from_this<Joint> {
public:
	Joint();
	Joint(std::shared_ptr<Body> body, double ndof, std::shared_ptr<Joint> parent = nullptr);
	virtual ~Joint();

	virtual void init();
	virtual void draw();
	virtual void update();

	double m_ndof;					// Number of DOF
	Eigen::VectorXd m_q;			// Position
	Eigen::VectorXd m_qdot;			// Velocity
	Eigen::VectorXd m_qddot;		// Acceleration
	Eigen::VectorXd m_tau;			// Joint torque
	Eigen::VectorXd m_tauCon;		// Constraint torque
	double m_K;						// Joint stiffness
	Eigen::MatrixXd m_S;			// Jacobian
	Eigen::MatrixXd m_Sdot;			// dS/dt



	Eigen::Matrix4d E_pj;			// Transform of this joint wrt parent joint
	Eigen::Matrix4d E_pj0;			// Transform when q is zero
	Eigen::Matrix4d E_jp;			// Transform of parent joint wrt this joint
	Matrix6d Ad_jp;					// Adjoint of E_jp
	Eigen::Matrix4d E_wj;			// Transform of this joint wrt world

	std::shared_ptr<Joint> next;	// Forward recursive ordering
	std::shared_ptr<Joint> prev;	// Reverse recursive ordering
	int idxR;						// Reduced indices



	void setJointTransform(Eigen::Matrix4d E);

	std::string getName() const { return m_name; }
	int getUID() const { return m_uid; }
	std::shared_ptr<Body> getBody() const { return m_body; }
	std::shared_ptr<Joint> getParent() const { return m_parent; }
	void addChild(std::shared_ptr<Joint> joint) { m_children.push_back(joint); }



private:
	std::string m_name;
	int m_uid;
	std::shared_ptr<Body> m_body;						// Attached body
	std::shared_ptr<Joint> m_parent;					// Parent joint
	std::vector<std::shared_ptr<Joint> > m_children;	// Children joints

};

#endif // MUSCLEMASS_SRC_JOINT_H_