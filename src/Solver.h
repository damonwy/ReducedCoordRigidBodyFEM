#pragma once
#ifndef MUSCLEMASS_SRC_SOLVER_H_
#define MUSCLEMASS_SRC_SOLVER_H_

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <json.hpp>
#include "MLCommon.h"

class World;

struct Solution {
	VectorXd t;
	Eigen::MatrixXd y;
};

class Solver 
{
public:
	Solver();
	Solver(std::shared_ptr<World> world, Integrator integrator);
	virtual ~Solver();
	std::shared_ptr<Solution> solve();


private:
	std::shared_ptr<World> m_world;
	Integrator m_integrator;
	std::shared_ptr<Solution> m_solutions;

	Eigen::MatrixXd M;
	Eigen::VectorXd f;
	Eigen::MatrixXd J;
	Eigen::MatrixXd Jdot;
	Eigen::VectorXd q0;
	Eigen::VectorXd qdot0;
	Eigen::MatrixXd Mtilde;
	Eigen::VectorXd ftilde;


};

#endif // MUSCLEMASS_SRC_SOLVER_H_