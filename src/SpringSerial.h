// SpringSerial Non-zero rest-length serial springs
//		A simple spring made up of a series of nodes. In the future, we may
//		want to make this an abstract class and derive from it. 

#pragma once
#ifndef MUSCLEMASS_SRC_SPRINGSERIAL_H_
#define MUSCLEMASS_SRC_SPRINGSERIAL_H_

#include "Spring.h"

class Body;
class Node;

class SpringSerial : public Spring 
{
public:
	SpringSerial();
	SpringSerial(int n_nodes, int &countS, int &countCM);
	virtual ~SpringSerial();
	void setStiffness(double K) { m_K = K; }
	void setMass(double mass) { m_mass = mass; }
	void setAttachments(std::shared_ptr<Body> body0, Eigen::Vector3d r0, std::shared_ptr<Body> body1, Eigen::Vector3d r1);

	virtual void init();
	virtual void load(const std::string &RESOURCE_DIR);
	virtual void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;

	virtual void countDofs_(int &nm, int &nr);
	virtual Eigen::VectorXd gatherDofs_(Eigen::VectorXd y, int nr);
	virtual Eigen::VectorXd gatherDDofs_(Eigen::VectorXd ydot, int nr);
	virtual void scatterDofs_(Eigen::VectorXd &y, int nr);
	virtual void scatterDDofs_(Eigen::VectorXd &ydot, int nr);

	virtual Eigen::MatrixXd computeMass_(Eigen::Vector3d grav, Eigen::MatrixXd M);
	virtual Eigen::VectorXd computeForce_(Eigen::Vector3d grav, Eigen::VectorXd f);
	virtual void computeEnergies_(Eigen::Vector3d grav, double &T, double &V);
	virtual Eigen::MatrixXd computeJacobian_(Eigen::MatrixXd J);


	

};


#endif // MUSCLEMASS_SRC_SPRINGSERIAL_H_