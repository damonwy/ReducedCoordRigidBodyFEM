#pragma once
#ifndef REDUCEDCOORD_SRC_DEFORMABLE_H_
#define REDUCEDCOORD_SRC_DEFORMABLE_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "MLCommon.h"

class MatrixStack;
class Program;
class Node;
class Body;

class Deformable {

public:
	Deformable();
	Deformable(int &countS, int &countCM);
	virtual ~Deformable() {}

	void setDamping(double damping) { m_damping = damping; }

	void countDofs(int &nm, int &nr);
	void gatherDofs(Eigen::VectorXd &y, int nr);
	void gatherDDofs(Eigen::VectorXd &ydot, int nr);
	void scatterDofs(Eigen::VectorXd &y, int nr);
	void scatterDDofs(Eigen::VectorXd &ydot, int nr);

	void computeJacobian(Eigen::MatrixXd &J, Eigen::MatrixXd &Jdot);
	void computeMass(Eigen::Vector3d grav, Eigen::MatrixXd &M, Eigen::VectorXd &f);
	void computeForceDamping(Eigen::Vector3d grav, Eigen::VectorXd &f, Eigen::MatrixXd D);
	void computeEnergies(Eigen::Vector3d grav, Energy &ener);

	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;
	void init();

	virtual void load(const std::string &RESOURCE_DIR) {}
	virtual void init_() {}
	virtual void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const {}

	virtual void countDofs_(int &nm, int &nr) {}
	virtual void gatherDofs_(Eigen::VectorXd &y, int nr) {}
	virtual void gatherDDofs_(Eigen::VectorXd &ydot, int nr) {}
	virtual void scatterDofs_(Eigen::VectorXd &y, int nr) {}
	virtual void scatterDDofs_(Eigen::VectorXd &ydot, int nr) {}

	virtual void computeMass_(Eigen::Vector3d grav, Eigen::MatrixXd &M, Eigen::VectorXd &f) {}
	virtual void computeForceDamping_(Eigen::Vector3d grav, Eigen::VectorXd &f, Eigen::MatrixXd D) {}
	virtual void computeEnergies_(Eigen::Vector3d grav, Energy &ener) {}
	virtual void computeJacobian_(Eigen::MatrixXd &J, Eigen::MatrixXd &Jdot) {}
	
	std::shared_ptr<Deformable> next;
	std::string m_name;
	int m_uid;

	std::vector<std::shared_ptr<Node>> m_nodes;
	double m_K;
	double m_damping;	
	std::shared_ptr<Body> m_body0;
	std::shared_ptr<Body> m_body1;
	Eigen::Vector3d m_r0;
	Eigen::Vector3d m_r1;

	double m_mass;
};

#endif // REDUCEDCOORD_SRC_DEFORMABLE_H_