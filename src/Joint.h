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
class MatrixStack;
class Program;


class Joint : public std::enable_shared_from_this<Joint> {
public:
	Joint();
	Joint(std::shared_ptr<Body> body, double ndof, std::shared_ptr<Joint> parent = nullptr);
	virtual ~Joint();

	virtual void init(int &nm, int &nr);
	virtual void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P) const;
	virtual void drawSelf(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P) const;
	virtual void update();
	virtual void updateSelf();

	int m_ndof;					// Number of DOF
	Eigen::VectorXd m_q;			// Position
	Eigen::VectorXd m_qdot;			// Velocity
	Eigen::VectorXd m_qddot;		// Acceleration
	Eigen::VectorXd m_tau;			// Joint torque
	Eigen::VectorXd m_tauCon;		// Constraint torque
	double m_K;						// Joint stiffness
	double m_D;						// Joint damping
	Vector6d m_S;			// Jacobian
	Vector6d m_Sdot;			// dS/dt

	Eigen::Matrix4d E_pj;			// Transform of this joint wrt parent joint
	Eigen::Matrix4d E_pj0;			// Transform when q is zero
	Eigen::Matrix4d E_jp;			// Transform of parent joint wrt this joint
	Matrix6d Ad_jp;					// Adjoint of E_jp
	Eigen::Matrix4d E_wj;			// Transform of this joint wrt world
	Eigen::Vector3d m_axis;
	std::shared_ptr<Joint> next;	// Forward recursive ordering
	std::shared_ptr<Joint> prev;	// Reverse recursive ordering
	int idxR;						// Reduced indices
	bool presc;						// Use presribed motion

	void countDofs(int &nm, int &nr);
	int countR(int &nr, int data);
	void setJointTransform(Eigen::Matrix4d E);
	void setStiffness(double K) { m_K = K; } // Sets this joint's linear stiffness
	void setDamping(double D) { m_D = D; } // Sets this joint's linear velocity damping

	std::string getName() const { return m_name; }
	int getUID() const { return m_uid; }
	std::shared_ptr<Body> getBody() const { return m_body; }
	std::shared_ptr<Joint> getParent() const { return m_parent; }
	void addChild(std::shared_ptr<Joint> joint) { m_children.push_back(joint); }
	std::shared_ptr<Joint> getJoint() { return shared_from_this(); }

	Eigen::MatrixXd computeJacobian(Eigen::MatrixXd J, int nm, int nr);
	Eigen::MatrixXd computeJacobianDerivative(Eigen::MatrixXd Jdot, Eigen::MatrixXd J, int nm, int nr);
	Eigen::VectorXd computerJacTransProd(Eigen::VectorXd y, Eigen::VectorXd x, int nr);
	Eigen::VectorXd computeForceStiffness(Eigen::VectorXd fr);
	Eigen::MatrixXd computeMatrixStiffness(Eigen::MatrixXd Ksr);
	
	Eigen::VectorXd computeForceDamping(Eigen::VectorXd fr);
	Eigen::MatrixXd computeMatrixDamping(Eigen::MatrixXd Ddr);

	Energy computeEnergies(Eigen::Vector3d grav, Energy ener);
	Eigen::VectorXd gatherDofs(Eigen::VectorXd y, int nr);
	Eigen::VectorXd gatherDDofs(Eigen::VectorXd ydot, int nr);
	void scatterDofs(Eigen::VectorXd y, int nr);
	void scatterDDofs(Eigen::VectorXd ydot, int nr);
	Vector6d getAlpha() const { return m_alpha; }
protected:
	Eigen::Matrix4d m_Q;										// Transformation matrix applied about the joint
private:
	void scatterDofsNoUpdate(Eigen::VectorXd y, int nr);

	std::string m_name;
	int m_uid;
	std::shared_ptr<Body> m_body;						// Attached body
	std::shared_ptr<Joint> m_parent;					// Parent joint
	std::vector<std::shared_ptr<Joint> > m_children;	// Children joints

	Vector6d m_alpha;
	

};

#endif // MUSCLEMASS_SRC_JOINT_H_