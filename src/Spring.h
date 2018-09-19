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


class Spring {

public:
	Spring();
	Spring(int &countS, int &countCM);
	virtual ~Spring();

	void countDofs();
	void gatherDofs(Eigen::VectorXd &y, int nr);
	void gatherDDofs(Eigen::VectorXd &ydot, int nr);

	void computeJacobian(Eigen::MatrixXd &J, Eigen::MatrixXd &Jdot);
	void computeMassForce(Eigen::Vector3d grav, Eigen::MatrixXd &M, Eigen::VectorXd &f);

	void computeEnergies(Eigen::Vector3d grav, double &T, double &V);

	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;
	virtual void load(const std::string &RESOURCE_DIR);
	virtual void init();

	virtual void countDofs_();
	virtual void gatherDofs_(Eigen::VectorXd &y, int nr);
	virtual void gatherDDofs_(Eigen::VectorXd &ydot, int nr);
	virtual void scatterDofs_(Eigen::VectorXd &y, int nr);
	virtual void scatterDDofs_(Eigen::VectorXd &ydot, int nr);
	virtual void computeMassForce_(Eigen::Vector3d grav, Eigen::MatrixXd &M, Eigen::VectorXd &f);
	virtual void computeEnergies_(Eigen::Vector3d grav, double &T, double &V);
	virtual void computeJacobian_(Eigen::MatrixXd &J, Eigen::MatrixXd &Jdot);
	virtual void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;
	std::shared_ptr<Spring> next;
	std::string m_name;
	int m_uid;
};

#endif // MUSCLEMASS_SRC_SPRING_H_