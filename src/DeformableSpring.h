// DeformableSpring Non-zero rest-length serial springs
//		A simple spring made up of a series of nodes. In the future, we may
//		want to make this an abstract class and derive from it. 

#pragma once
#ifndef MUSCLEMASS_SRC_DEFORMABLESPRING_H_
#define MUSCLEMASS_SRC_DEFORMABLESPRING_H_

#include "Deformable.h"

class Body;
class Node;

class DeformableSpring : public Deformable
{
public:
	DeformableSpring();
	DeformableSpring(int n_nodes, int &countS, int &countCM);
	virtual ~DeformableSpring() {}

	void setStiffness(double K) { m_K = K; }
	void setMass(double mass) { m_mass = mass; }
	void setAttachments(std::shared_ptr<Body> body0, Eigen::Vector3d r0, std::shared_ptr<Body> body1, Eigen::Vector3d r1);

protected:
	void init_();
	virtual void load(const std::string &RESOURCE_DIR);
	void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;

	void countDofs_(int &nm, int &nr);
	void gatherDofs_(Eigen::VectorXd &y, int nr);
	void gatherDDofs_(Eigen::VectorXd &ydot, int nr);
	void scatterDofs_(Eigen::VectorXd &y, int nr);
	void scatterDDofs_(Eigen::VectorXd &ydot, int nr);

	void computeMass_(Eigen::Vector3d grav, Eigen::MatrixXd &M, Eigen::VectorXd &f);
	void computeForceDamping_(Eigen::Vector3d grav, Eigen::VectorXd &f, Eigen::MatrixXd &D);
	void computeEnergies_(Eigen::Vector3d grav, Energy &ener);
	void computeJacobian_(Eigen::MatrixXd &J, Eigen::MatrixXd &Jdot);

};


#endif // MUSCLEMASS_SRC_DEFORMABLESPRING_H_