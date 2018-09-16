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
	virtual ~Spring();

	void countDofs();
	void gatherDofs(double &y);
	void gatherDDofs(double &ydot);

	void computeJacobian(Eigen::MatrixXd &J, Eigen::MatrixXd &Jdot);
	void computeMassForce(Eigen::Vector3d grav, Eigen::MatrixXd M, Eigen::VectorXd f);

	void computeEnergies(Eigen::Vector3d grav, double T, double V);


	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P) const;

	std::shared_ptr<Spring> next;
	std::string m_name;
	int m_uid;
};



#endif // MUSCLEMASS_SRC_SPRING_H_