#pragma once
#ifndef MUSCLEMASS_SRC_SOLVER_H_
#define MUSCLEMASS_SRC_SOLVER_H_

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen\src\Core\util\IndexedViewHelper.h>

#include <json.hpp>
#include "MLCommon.h"

class World;

typedef Eigen::Triplet<double> T;

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
	virtual ~Solver() {}
	std::shared_ptr<Solution> solve();
	Eigen::VectorXd dynamics(Eigen::VectorXd y);
	Eigen::VectorXd dynamics_sparse(Eigen::VectorXd y);
	void initMatrix(int nm, int nr, int nem, int ner, int nim, int nir);
	void initMatrixSparse(int nm, int nr, int nem, int ner, int nim, int nir);
	void init();
	void reset();
	void load(const std::string &RESOURCE_DIR);
	bool isSparse;
private:
	int nr;
	int nm;

	std::shared_ptr<World> m_world;
	Integrator m_integrator;
	std::shared_ptr<Solution> m_solutions;

	Eigen::MatrixXd Mm;
	Eigen::SparseMatrix<double> Mm_sp;
	std::vector<T> Mm_;

	Eigen::MatrixXd MDKr_;
	Eigen::SparseMatrix<double> MDKr_sp;
	Eigen::MatrixXd K;
	Eigen::SparseMatrix<double> K_sp;
	std::vector<T> K_;

	Eigen::MatrixXd Km;
	Eigen::SparseMatrix<double> Km_sp;
	std::vector<T> Km_;

	Eigen::VectorXd fm;
	Eigen::MatrixXd J;
	Eigen::MatrixXd J_dense;	// dense_nm x dense_nr
	Eigen::MatrixXd Jdot_dense;

	Eigen::SparseMatrix<double> J_sp;
	std::vector<T> J_;
	Eigen::MatrixXd Jdot;
	Eigen::SparseMatrix<double> Jdot_sp;
	std::vector<T> Jdot_;
	Eigen::VectorXd q0;
	Eigen::VectorXd q1;
	Eigen::VectorXd qdot0;
	Eigen::VectorXd qdot1;
	Eigen::VectorXd qddot;

	Eigen::MatrixXd Mr;
	Eigen::SparseMatrix<double> Mr_sp;
	Eigen::MatrixXd Mr_temp;
	Eigen::SparseMatrix<double> Mr_sp_temp;

	Eigen::MatrixXd Dm; // nm x nm
	Eigen::SparseMatrix<double> Dm_sp;
	std::vector<T> Dm_;
	Eigen::VectorXd tmp; // nm x 1
	Eigen::MatrixXd Dr;
	Eigen::SparseMatrix<double> Dr_sp;
	std::vector<T> Dr_;

	Eigen::MatrixXd Kr;
	Eigen::SparseMatrix<double> Kr_sp;
	std::vector<T> Kr_;
	Eigen::VectorXd fr;
	Eigen::VectorXd fr_;

	Eigen::MatrixXd Gm;
	Eigen::SparseMatrix<double> Gm_sp;
	std::vector<T> Gm_;

	Eigen::MatrixXd Gmdot;
	Eigen::SparseMatrix<double> Gmdot_sp;
	std::vector<T> Gmdot_;

	Eigen::VectorXd gm;
	Eigen::VectorXd gmdot;
	Eigen::VectorXd gmddot;

	Eigen::MatrixXd Gr;
	Eigen::SparseMatrix<double> Gr_sp;
	std::vector<T> Gr_;

	Eigen::MatrixXd Grdot;
	Eigen::SparseMatrix<double> Grdot_sp;
	std::vector<T> Grdot_;

	Eigen::VectorXd gr;
	Eigen::VectorXd grdot;
	Eigen::VectorXd grddot;

	Eigen::MatrixXd G;
	//Eigen::SparseMatrix<double> G_sp;
	
	Eigen::VectorXd g;
	Eigen::VectorXd gdot;
	Eigen::VectorXd rhsG;

	Eigen::MatrixXd Cm;
	Eigen::SparseMatrix<double> Cm_sp;
	Eigen::MatrixXd Cmdot;
	Eigen::SparseMatrix<double> Cmdot_sp;
	Eigen::VectorXd cm;
	Eigen::VectorXd cmdot;
	Eigen::VectorXd cmddot;

	Eigen::MatrixXd Cr;
	Eigen::SparseMatrix<double> Cr_sp;
	Eigen::MatrixXd Crdot;
	Eigen::SparseMatrix<double> Crdot_sp;
	Eigen::VectorXd cr;
	Eigen::VectorXd crdot;
	Eigen::VectorXd crddot;

	Eigen::VectorXd rhsC;

	Eigen::MatrixXd C;
	Eigen::SparseMatrix<double> C_sp;
	Eigen::VectorXd c;

	std::vector<int> rowsM;
	std::vector<int> rowsR;

};

#endif // MUSCLEMASS_SRC_SOLVER_H_