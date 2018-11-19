#pragma once
#include "Solver.h"

class SolverSparse : public Solver {
public:
	SolverSparse() {}
	SolverSparse(std::shared_ptr<World> world, Integrator integrator) : Solver(world, integrator) {}
	Eigen::VectorXd dynamics(Eigen::VectorXd y);
	void initMatrix(int nm, int nr, int nem, int ner, int nim, int nir);

private:
	Eigen::SparseMatrix<double> Mm_sp;
	std::vector<T> Mm_;

	Eigen::SparseMatrix<double> MDKr_sp;
	Eigen::SparseMatrix<double> K_sp;
	std::vector<T> K_;

	Eigen::SparseMatrix<double> Km_sp;
	std::vector<T> Km_;

	Eigen::VectorXd fm;
	Eigen::MatrixXd J_dense;	// dense_nm x dense_nr
	Eigen::MatrixXd Jdot_dense;

	Eigen::SparseMatrix<double> J_sp;
	std::vector<T> J_;
	Eigen::SparseMatrix<double> Jdot_sp;
	std::vector<T> Jdot_;
	Eigen::VectorXd q0;
	Eigen::VectorXd q1;
	Eigen::VectorXd qdot0;
	Eigen::VectorXd qdot1;
	Eigen::VectorXd qddot;

	Eigen::SparseMatrix<double> Mr_sp;
	Eigen::SparseMatrix<double> Mr_sp_temp;

	Eigen::SparseMatrix<double> Dm_sp;
	std::vector<T> Dm_;
	Eigen::VectorXd tmp; // nm x 1
	Eigen::SparseMatrix<double> Dr_sp;
	std::vector<T> Dr_;

	Eigen::SparseMatrix<double> Kr_sp;
	std::vector<T> Kr_;
	Eigen::VectorXd fr;
	Eigen::VectorXd fr_;

	Eigen::SparseMatrix<double> Gm_sp;
	std::vector<T> Gm_;

	Eigen::SparseMatrix<double> Gmdot_sp;
	std::vector<T> Gmdot_;

	Eigen::VectorXd gm;
	Eigen::VectorXd gmdot;
	Eigen::VectorXd gmddot;

	Eigen::SparseMatrix<double> Gr_sp;
	std::vector<T> Gr_;

	Eigen::SparseMatrix<double> Grdot_sp;
	std::vector<T> Grdot_;

	Eigen::VectorXd gr;
	Eigen::VectorXd grdot;
	Eigen::VectorXd grddot;

	Eigen::VectorXd g;
	Eigen::VectorXd gdot;
	Eigen::VectorXd rhsG;

	Eigen::SparseMatrix<double> Cm_sp;
	Eigen::SparseMatrix<double> Cmdot_sp;
	Eigen::VectorXd cm;
	Eigen::VectorXd cmdot;
	Eigen::VectorXd cmddot;

	Eigen::SparseMatrix<double> Cr_sp;
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