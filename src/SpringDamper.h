// SpringDamper 
//		f = k * (l-L)/L - d * v
//		
#pragma once
#ifndef MUSCLEMASS_SRC_SPRINGDAMPER_H_
#define MUSCLEMASS_SRC_SPRINGDAMPER_H_

#include "Spring.h"

class Body;
class Node;

class SpringDamper : public Spring
{
public:
	SpringDamper();
	SpringDamper(int &countS, int &countCM);
	virtual ~SpringDamper() {}

	void setStiffness(double K) { m_K = K; }
	void setMass(double mass) { m_mass = mass; }
	void setDamping(double damping) { m_d = damping; }
	void setRestLength(double L) { m_L = L; }
	void setAttachments(std::shared_ptr<Body> body0, Eigen::Vector3d r0, std::shared_ptr<Body> body1, Eigen::Vector3d r1);

	void init();
	void load(const std::string &RESOURCE_DIR);
	void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;

	void countDofs_(int &nm, int &nr);
	Eigen::VectorXd gatherDofs_(Eigen::VectorXd y, int nr);
	Eigen::VectorXd gatherDDofs_(Eigen::VectorXd ydot, int nr);
	void scatterDofs_(Eigen::VectorXd &y, int nr);
	void scatterDDofs_(Eigen::VectorXd &ydot, int nr);

	Eigen::MatrixXd computeMass_(Eigen::Vector3d grav, Eigen::MatrixXd M);
	Eigen::VectorXd computeForce_(Eigen::Vector3d grav, Eigen::VectorXd f);
	Energy computeEnergies_(Eigen::Vector3d grav, Energy ener);
	Eigen::MatrixXd computeJacobian_(Eigen::MatrixXd J);

protected:
	double m_L;		// Rest Length
	double m_d;		// Damping parameter
	double m_l;		// Current Length



};


#endif // MUSCLEMASS_SRC_SPRINGSERIAL_H_