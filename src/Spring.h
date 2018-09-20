#pragma once
#ifndef MUSCLEMASS_SRC_SPRING_H_
#define MUSCLEMASS_SRC_SPRING_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "MLCommon.h"

class MatrixStack;
class Program;
class Node;
class Body;

class Spring {

public:
	Spring();
	Spring(int &countS, int &countCM);
	virtual ~Spring();

	void countDofs(int &nm, int &nr);
	Eigen::VectorXd gatherDofs(Eigen::VectorXd y, int nr);
	Eigen::VectorXd gatherDDofs(Eigen::VectorXd ydot, int nr);
	void scatterDofs(Eigen::VectorXd &y, int nr);
	void scatterDDofs(Eigen::VectorXd &ydot, int nr);

	Eigen::MatrixXd computeJacobian(Eigen::MatrixXd J);
	Eigen::MatrixXd computeMass(Eigen::Vector3d grav, Eigen::MatrixXd M);
	Eigen::VectorXd computeForce(Eigen::Vector3d grav, Eigen::VectorXd f);
	void computeEnergies(Eigen::Vector3d grav, double &T, double &V);

	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;
	virtual void load(const std::string &RESOURCE_DIR);
	virtual void init();

	virtual void countDofs_(int &nm, int &nr);
	virtual Eigen::VectorXd gatherDofs_(Eigen::VectorXd y, int nr);
	virtual Eigen::VectorXd gatherDDofs_(Eigen::VectorXd ydot, int nr);
	virtual void scatterDofs_(Eigen::VectorXd &y, int nr);
	virtual void scatterDDofs_(Eigen::VectorXd &ydot, int nr);
	virtual Eigen::MatrixXd computeMass_(Eigen::Vector3d grav, Eigen::MatrixXd M);
	virtual Eigen::VectorXd computeForce_(Eigen::Vector3d grav, Eigen::VectorXd f);
	virtual void computeEnergies_(Eigen::Vector3d grav, double &T, double &V);
	virtual Eigen::MatrixXd computeJacobian_(Eigen::MatrixXd J);
	virtual void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;
	std::shared_ptr<Spring> next;
	std::string m_name;
	int m_uid;

	std::vector<std::shared_ptr<Node>> m_nodes;
	double m_K;
	std::shared_ptr<Body> m_body0;
	std::shared_ptr<Body> m_body1;
	Eigen::Vector3d m_r0;
	Eigen::Vector3d m_r1;

	double m_mass;
};

#endif // MUSCLEMASS_SRC_SPRING_H_