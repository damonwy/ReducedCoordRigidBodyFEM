// Joint A generic joint between two rigid bodies
//    A joint is defined between a parent and the child. The DOF, q, of the joint is 
//    the relative displacement of the child wrt the parent

#pragma once
#ifndef REDUCEDCOORD_SRC_JOINT_H_
#define REDUCEDCOORD_SRC_JOINT_H_
#include <vector>
#include <memory>
#define EIGEN_DONT_ALIGN_STATICALLY

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "MLCommon.h"

class SE3;
class Body;
class MatrixStack;
class Program;
class Shape;
typedef Eigen::Triplet<double> T;

class Joint : public std::enable_shared_from_this<Joint> {
public:
	Joint();
	Joint(std::shared_ptr<Body> body, int ndof, std::shared_ptr<Joint> parent = nullptr);
	virtual ~Joint() {}

	virtual void load(const std::string &RESOURCE_DIR, std::string joint_shape);
	virtual void init(int &nm, int &nr);
	virtual void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const;
	virtual void update();

	int m_ndof;						// Number of DOF
	Eigen::VectorXd m_q;			// Position
	Eigen::VectorXd m_qdot;			// Velocity
	Eigen::VectorXd m_qddot;		// Acceleration
	Eigen::VectorXd m_tau;			// Joint torque
	Eigen::VectorXd m_tauCon;		// Constraint torque	
	
	bool presc;						// Presribed motion constraint
	double m_Kr;					// Joint stiffness
	double m_Dr;					// Joint damping
	Eigen::MatrixXd m_S;			// Jacobian
	Eigen::MatrixXd m_Sdot;			// dS/dt
	Matrix6d m_I_j;					// Inertia at the joint
	Vector6d V;						// Twist at parent joint
	Vector6d Vdot;					// Acceleration at parent joint
	Matrix4d E_pj;					// Transform of this joint wrt parent joint
	Matrix4d E_pj0;					// Transform when q is zero
	Matrix4d E_jp;					// Transform of parent joint wrt this joint
	Matrix4d E_wj;					// Transform of this joint wrt world
	Matrix6d Ad_jp;					// Adjoint of E_jp
	std::shared_ptr<Joint> next;	// Forward recursive ordering
	std::shared_ptr<Joint> prev;	// Reverse recursive ordering
	int idxR;						// Reduced indices
	Vector3d m_axis;

	void countDofs(int &nm, int &nr);
	int countR(int &nr, int data);
	void setJointTransform(Matrix4d E);
	void setStiffness(double K) { m_Kr = K; } // Sets this joint's linear stiffness
	void setDamping(double D) { m_Dr = D; } // Sets this joint's linear velocity damping
	void addChild(std::shared_ptr<Joint> joint) { m_children.push_back(joint); }

	std::shared_ptr<Body> getBody() const { return m_body; }
	std::shared_ptr<Joint> getParent() const { return m_parent; }
	std::shared_ptr<Joint> getJoint() { return shared_from_this(); }
	Vector6d getAlpha() const { return m_alpha; }
	std::string getName() const { return m_name; }

	void computeJacobian(Eigen::MatrixXd &J, Eigen::MatrixXd &Jdot, int nm, int nr);
	void computeJacobianSparse(std::vector<T> &J_, std::vector<T> &Jdot_, int nm, int nr);
	Eigen::VectorXd computerJacTransProd(Eigen::VectorXd y, Eigen::VectorXd x, int nr);
	void computeForceStiffness(Eigen::VectorXd &fr, Eigen::MatrixXd &Kr);
	void computeForceStiffnessSparse(Eigen::VectorXd &fr, std::vector<T> &Kr_);
	void computeForceDamping(Eigen::VectorXd &fr, Eigen::MatrixXd &Dr);
	void computeForceDampingSparse(Eigen::VectorXd &fr, std::vector<T> &Dr_);
	void computeInertia();

	void computeEnergies(Vector3d grav, Energy &ener);
	Eigen::VectorXd gatherDofs(Eigen::VectorXd y, int nr);
	Eigen::VectorXd gatherDDofs(Eigen::VectorXd ydot, int nr);
	void scatterDofs(Eigen::VectorXd y, int nr);
	void scatterDDofs(Eigen::VectorXd ydot, int nr);
	void scatterTauCon(Eigen::VectorXd tauc);

protected:
	Matrix4d m_Q;										// Transformation matrix applied about the joint
	std::shared_ptr<Shape> m_jointShape;				// Joint shape			
	std::shared_ptr<Body> m_body;						// Attached body
	std::shared_ptr<Joint> m_parent;					// Parent joint
	std::vector<std::shared_ptr<Joint>> m_children;		// Children joints
	virtual void update_() {}

private:
	void scatterDofsNoUpdate(Eigen::VectorXd y, int nr);
	std::string m_name;
	Vector6d m_alpha;									// For J'*x product
	virtual void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const;

};

#endif // REDUCEDCOORD_SRC_JOINT_H_