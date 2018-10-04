#pragma once
#ifndef MUSCLEMASS_SRC_SOLVER_H_
#define MUSCLEMASS_SRC_SOLVER_H_

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen\src\Core\util\IndexedViewHelper.h>

#include <json.hpp>
#include "MLCommon.h"

class World;

struct Solution {
	Eigen::VectorXd t;
	Eigen::MatrixXd y;

	Solution(){}

	Eigen::VectorXd step(int time_step) {
		return y.row(time_step);
	}

	int getNsteps() {
		return y.rows();
	}

	void searchTime(double ti, int search_index, int &output_index, double &s) {
		// Finds the index of the time interval around ti
		
		while (ti > t(search_index))
		{
			search_index++;
		}
		output_index = search_index - 1;
		if (output_index < 0) {
			// Beginning of time
			output_index = 0;
			s = 0.0;
		}
		else {
			if (output_index == y.rows() - 1) {
				// End of time
				output_index -= 1;
				s = 1.0;
			}
			else {
				// Somewhere between
				double t0 = t(output_index);
				double t1 = t(output_index + 1);
				s = (ti - t0) / (t1 - t0);
			}
		}
	}

};

class Solver 
{
public:
	Solver();
	Solver(std::shared_ptr<World> world, Integrator integrator);
	virtual ~Solver();
	std::shared_ptr<Solution> solve();
	void init();
	void load(const std::string &RESOURCE_DIR);
	
private:
	int nr;
	int nm;

	std::shared_ptr<World> m_world;
	Integrator m_integrator;
	std::shared_ptr<Solution> m_solutions;

	Eigen::MatrixXd M;
	Eigen::MatrixXd K;
	Eigen::VectorXd f;
	Eigen::MatrixXd J;
	Eigen::MatrixXd Jdot;
	Eigen::VectorXd q0;
	Eigen::VectorXd q1;
	Eigen::VectorXd qdot0;
	Eigen::VectorXd qdot1;
	Eigen::VectorXd qddot;

	Eigen::MatrixXd Mtilde;
	Eigen::VectorXd ftilde;

	Eigen::MatrixXd Gm;
	Eigen::MatrixXd Gmdot;
	Eigen::VectorXd gm;
	Eigen::VectorXd gmdot;
	Eigen::VectorXd gmddot;

	Eigen::MatrixXd Gr;
	Eigen::MatrixXd Grdot;
	Eigen::VectorXd gr;

	Eigen::MatrixXd G;
	Eigen::VectorXd g;
	Eigen::VectorXd gdot;
	Eigen::VectorXd rhsG;

	Eigen::MatrixXd Cm;
	Eigen::MatrixXd Cmdot;
	Eigen::VectorXd cm;

	Eigen::MatrixXd Cr;
	Eigen::MatrixXd Crdot;
	Eigen::VectorXd cr;

	Eigen::MatrixXd C;
	Eigen::VectorXd c;

	std::vector<int> rowsM;
	std::vector<int> rowsR;

};

#endif // MUSCLEMASS_SRC_SOLVER_H_