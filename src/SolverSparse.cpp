#include "SolverSparse.h"
#include "World.h"
#include "Body.h"
#include "SoftBody.h"
#include "Joint.h"
#include "Spring.h"
#include "SpringDamper.h"
#include "Deformable.h"
#include "DeformableSpring.h"
#include "ConstraintJointLimit.h"
#include "ConstraintLoop.h"
#include "ConstraintAttachSpring.h"
#include "QuadProgMosek.h"
#include "MatlabDebug.h"
#include "ChronoTimer.h"
#include <iostream>
#include <fstream>
#include <json.hpp>

using namespace std;
using namespace Eigen;

void SolverSparse::initMatrix(int nm, int nr, int nem, int ner, int nim, int nir) {
	int ne = ner + nem;
	int ni = nim + nir;
	Mm_sp.resize(nm, nm);
	Mm_sp.setZero();
	Mr_sp.resize(nr, nr);
	Mr_sp.setZero();
	Mm_.clear();

	Mr_sp_temp.resize(nr, nr);
	Mr_sp_temp.setZero();
	MDKr_sp.resize(nr, nr);
	MDKr_sp.setZero();

	fm.resize(nm);
	fr.resize(nr);
	fr_.resize(nr);
	fm.setZero();
	fr.setZero();
	fr_.setZero();
	tmp.resize(nm);
	tmp.setZero();

	Kr_sp.resize(nr, nr);
	Kr_sp.setZero();
	Kr_.clear();
	Dr_sp.resize(nr, nr);
	Dr_sp.setZero();
	Dr_.clear();
	Dm_sp.resize(nm, nm);
	Dm_sp.setZero();
	Dm_.clear();

	K_sp.resize(nm, nm);
	K_sp.data().squeeze();
	K_sp.setZero();
	K_.clear();

	Km_sp.resize(nm, nm);
	Km_sp.setZero();
	Km_.clear();

	J_dense.resize(m_world->m_dense_nm, m_world->m_dense_nr);
	J_dense.setZero();
	Jdot_dense.resize(m_world->m_dense_nm, m_world->m_dense_nr);
	Jdot_dense.setZero();

	J_sp.resize(nm, nr);
	J_sp.setZero();
	J_.clear();
	Jdot_sp.resize(nm, nr);
	Jdot_sp.setZero();
	Jdot_.clear();

	Gm_sp.resize(nem, nm);
	Gm_sp.data().squeeze();
	Gm_sp.setZero();
	Gmdot_sp.resize(nem, nm);
	Gmdot_sp.setZero();
	Gm_.clear();
	Gmdot_.clear();

	gm.resize(nem);
	gm.setZero();
	gmdot.resize(nem);
	gmdot.setZero();
	gmddot.resize(nem);
	gmddot.setZero();

	Gr_sp.resize(ner, nr);
	Gr_sp.data().squeeze();
	Gr_sp.setZero();
	Grdot_sp.resize(ner, nr);
	Grdot_sp.setZero();
	Gr_.clear();
	Grdot_.clear();

	gr.resize(ner);
	gr.setZero();
	grdot.resize(ner);
	grdot.setZero();
	grddot.resize(ner);
	grddot.setZero();

	g.resize(ne);
	rhsG.resize(ne);
	gdot.resize(ne);
	g.setZero();
	gdot.setZero();
	rhsG.setZero();

	Cm_sp.resize(nim, nm);
	Cm_sp.setZero();
	Cmdot_sp.resize(nim, nm);
	Cmdot_sp.setZero();
	cm.resize(nim);
	cm.setZero();
	cmdot.resize(nim);
	cmdot.setZero();
	cmddot.resize(nim);
	cmddot.setZero();

	Cr_sp.resize(nir, nr);
	Cr_sp.setZero();
	Crdot_sp.resize(nir, nr);
	Crdot_sp.setZero();
	cr.resize(nir);
	cr.setZero();
	crdot.resize(nir);
	crdot.setZero();
	crddot.resize(nir);
	crddot.setZero();
}

VectorXd SolverSparse::dynamics(VectorXd y)
{
	ChronoTimer timer("Test");
	SparseMatrix<double, RowMajor> G_sp;

	switch (m_integrator)
	{
	case REDMAX_EULER:
	{
		int nr = m_world->nr;
		int nm = m_world->nm;

		int nem = m_world->nem;
		int ner = m_world->ner;
		int ne = nem + ner;

		int nim = m_world->nim;
		int nir = m_world->nir;
		int ni = nim + nir;
		G_sp.resize(ne, nr);

		initMatrix(nm, nr, nem, ner, nim, nir);

		auto body0 = m_world->getBody0();
		auto joint0 = m_world->getJoint0();
		auto deformable0 = m_world->getDeformable0();
		auto softbody0 = m_world->getSoftBody0();
		auto constraint0 = m_world->getConstraint0();
		auto spring0 = m_world->getSpring0();

		double t = m_world->getTspan()(0);
		double h = m_world->getH();
		Vector3d grav = m_world->getGrav();

		VectorXd yk(2 * nr);
		VectorXd ydotk(2 * nr);

		body0->computeMassGravSparse(grav, Mm_, fm);
		body0->computeForceDampingSparse(tmp, Dm_);

		deformable0->computeMassSparse(grav, Mm_, fm);
		deformable0->computeForceDampingSparse(grav, tmp, Dm_);

		softbody0->computeMassSparse(grav, Mm_);
		softbody0->computeForce(grav, fm);

		softbody0->computeStiffnessSparse(K_);

		joint0->computeForceStiffnessSparse(fr, Kr_);
		joint0->computeForceDampingSparse(tmp, Dr_);

		//// First get dense jacobian (only a small part of the matrix)
		joint0->computeJacobian(J_dense, Jdot_dense);

		//// Push back the dense part
		for (int i = 0; i < J_dense.rows(); ++i) {
			for (int j = 0; j < J_dense.cols(); ++j) {

				J_.push_back(T(i, j, J_dense(i, j)));
				Jdot_.push_back(T(i, j, Jdot_dense(i, j)));
			}
		}

		deformable0->computeJacobianSparse(J_, Jdot_);
		softbody0->computeJacobianSparse(J_);

		spring0->computeForceStiffnessDampingSparse(fm, Km_, Dm_);

		Mm_sp.setFromTriplets(Mm_.begin(), Mm_.end());
		Km_sp.setFromTriplets(Km_.begin(), Km_.end());
		Dm_sp.setFromTriplets(Dm_.begin(), Dm_.end());
		Dr_sp.setFromTriplets(Dr_.begin(), Dr_.end());
		K_sp.setFromTriplets(K_.begin(), K_.end()); // check
		Kr_sp.setFromTriplets(Kr_.begin(), Kr_.end());
		J_sp.setFromTriplets(J_.begin(), J_.end()); // check
		Jdot_sp.setFromTriplets(Jdot_.begin(), Jdot_.end());

		q0 = y.segment(0, nr);
		qdot0 = y.segment(nr, nr);

		Mr_sp = J_sp.transpose() * (Mm_sp - h * h * K_sp) * J_sp;
		//sparse_to_file_as_dense(K_sp, "K_sp");

		//Mr_sp_temp = Mr_sp.transpose();
		//Mr_sp += Mr_sp_temp;
		//Mr_sp *= 0.5;

		fr_ = Mr_sp * qdot0 + h * (J_sp.transpose() * (fm - Mm_sp * Jdot_sp * qdot0) + fr); // check
		MDKr_sp = Mr_sp + J_sp.transpose() * (h * Dm_sp - h * h * Km_sp) * J_sp + h * Dr_sp - h * h * Kr_sp;
		/*sparse_to_file_as_dense(Mr_sp, "Mr_sp");
		sparse_to_file_as_dense(J_sp, "J_sp");
		sparse_to_file_as_dense(Jdot_sp, "Jdot_sp");
		sparse_to_file_as_dense(Km_sp, "Km_sp");
		sparse_to_file_as_dense(Dr_sp, "Dr_sp");
		sparse_to_file_as_dense(Kr_sp, "Kr_sp");
		sparse_to_file_as_dense(K_sp, "K_sp");
		sparse_to_file_as_dense(Mm_sp, "Mm_sp");
*/
		
		if (ne > 0) {
			constraint0->computeJacEqMSparse(Gm_, Gmdot_, gm, gmdot, gmddot);
			constraint0->computeJacEqRSparse(Gr_, Grdot_, gr, grdot, grddot);

			Gm_sp.setFromTriplets(Gm_.begin(), Gm_.end());
			Gmdot_sp.setFromTriplets(Gmdot_.begin(), Gmdot_.end());
			Gr_sp.setFromTriplets(Gr_.begin(), Gr_.end());
			Grdot_sp.setFromTriplets(Grdot_.begin(), Grdot_.end());

			//sparse_to_file_as_dense(Gm_sp * J_sp, "Gm_sp * J_sp");
			G_sp.topRows(nem) = Gm_sp * J_sp;
			G_sp.bottomRows(ner) = Gr_sp;
			
			//sparse_to_file_as_dense(G_sp, "G_sp");

			g.segment(0, nem) = gm;
			g.segment(nem, ner) = gr;
			rhsG = -gdot - 100.0 * g;// todo!!!!!
		}

		//if (ni > 0) {
		//	// Check for active inequality constraint
		//	constraint0->computeJacIneqMSparse(Cm, Cmdot, cm, cmdot, cmddot);
		//	constraint0->computeJacIneqRSparse(Cr, Crdot, cr, crdot, crddot);
		//	rowsR.clear();
		//	rowsM.clear();

		//	constraint0->getActiveList(rowsM, rowsR);
		//	nim = rowsM.size();
		//	nir = rowsR.size();
		//	ni = nim + nir;

		//	if (ni > 0) {
		//		Eigen::VectorXi m_rowsM = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(rowsM.data(), rowsM.size());
		//		Eigen::VectorXi m_rowsR = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(rowsR.data(), rowsR.size());

		//		MatrixXd m_Cm = Cm(m_rowsM, Eigen::placeholders::all);
		//		MatrixXd m_Cr = Cr(m_rowsR, Eigen::placeholders::all);
		//		VectorXd m_cm = cm(m_rowsM);
		//		VectorXd m_cr = cr(m_rowsR);
		//		VectorXd m_cmdot = cmdot(m_rowsM);
		//		VectorXd m_crdot = crdot(m_rowsR);

		//		MatrixXd CmJ = m_Cm * J;
		//		C.resize(CmJ.rows() + m_Cr.rows(), m_Cr.cols());
		//		C << CmJ, m_Cr;
		//		rhsC.resize(C.rows());
		//		VectorXd c(C.rows());
		//		c << m_cm, m_cr;
		//		VectorXd cdot(C.rows());
		//		cdot << m_cmdot, m_crdot;
		//		rhsC = -cdot - 105.0 * c;
		//	}
		//}

		if (ne == 0 && ni == 0) {	// No constraints
			ConjugateGradient< SparseMatrix<double> > cg;
			cg.setMaxIterations(25);
			cg.setTolerance(1e-3);
			cg.compute(MDKr_sp);
			qdot1 = cg.solveWithGuess(fr_, qdot0);
		}
		else if (ne > 0 && ni == 0) {  // Just equality
			int rows = MDKr_sp.rows() + G_sp.rows();
			int cols = MDKr_sp.cols() + G_sp.rows();
			MatrixXd LHS(rows, cols);
			VectorXd rhs(rows);
			LHS.setZero();
			rhs.setZero();

			MatrixXd MDKr_ = MatrixXd(MDKr_sp);
			MatrixXd G = MatrixXd(G_sp);

			LHS.block(0, 0, MDKr_.rows(), MDKr_.cols()) = MDKr_;
			LHS.block(0, MDKr_.cols(), MDKr_.rows(), G.rows()) = G.transpose();
			LHS.block(MDKr_.rows(), 0, G.rows(), G.cols()) = G;
			rhs.segment(0, fr_.rows()) = fr_;
			rhs.segment(fr_.rows(), g.rows()) = rhsG;

			VectorXd sol = LHS.ldlt().solve(rhs);
			qdot1 = sol.segment(0, nr);
			VectorXd l = sol.segment(nr, sol.rows() - nr);
			//shared_ptr<QuadProgMosek> program_ = make_shared <QuadProgMosek>();
			//program_->setParamInt(MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT);
			//program_->setParamInt(MSK_IPAR_LOG, 10);
			//program_->setParamInt(MSK_IPAR_LOG_FILE, 1);
			//program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_DFEAS, 1e-8);
			//program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_INFEAS, 1e-10);
			//program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_MU_RED, 1e-8);
			//program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_NEAR_REL, 1e3);
			//program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_PFEAS, 1e-8);
			//program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_REL_GAP, 1e-8);
			//program_->setNumberOfVariables(nr);
			//program_->setObjectiveMatrix(MDKr_sp);
			////sparse_to_file_as_dense(MDKr_sp, "MDKr_s");

			//program_->setObjectiveVector(-fr_);
			////vec_to_file(fr_, "fr_s");

			//program_->setNumberOfEqualities(ne);
			//program_->setEqualityMatrix(G_sp);
			////sparse_to_file_as_dense(G_sp, "G_s");

			//program_->setEqualityVector(rhsG);

			////vec_to_file(rhsG, "rhsG_s");

			//bool success = program_->solve();
			//VectorXd sol = program_->getPrimalSolution();
			//qdot1 = sol.segment(0, nr);
			////vec_to_file(qdot1, "qdot1_s");
			//VectorXd l = program_->getDualEquality();
			//vec_to_file(l, "l_s");
			constraint0->scatterForceEqM(MatrixXd(Gm_sp.transpose()), l.segment(0, nem) / h);
			constraint0->scatterForceEqR(MatrixXd(Gr_sp.transpose()), l.segment(nem, l.rows() - nem) / h);

		}
		else if (ne == 0 && ni > 0) {  // Just inequality
			shared_ptr<QuadProgMosek> program_ = make_shared <QuadProgMosek>();
			program_->setParamInt(MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT);
			program_->setParamInt(MSK_IPAR_LOG, 10);
			program_->setParamInt(MSK_IPAR_LOG_FILE, 1);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_DFEAS, 1e-8);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_INFEAS, 1e-10);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_MU_RED, 1e-8);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_NEAR_REL, 1e3);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_PFEAS, 1e-8);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_REL_GAP, 1e-8);

			program_->setNumberOfVariables(nr);
			program_->setObjectiveMatrix(MDKr_sp);
			program_->setObjectiveVector(-fr_);
			program_->setNumberOfInequalities(ni);
			program_->setInequalityMatrix(C.sparseView());

			VectorXd cvec(ni);
			cvec.setZero();
			program_->setInequalityVector(cvec);

			bool success = program_->solve();
			VectorXd sol = program_->getPrimalSolution();
			qdot1 = sol.segment(0, nr);
		}
		//else {  // Both equality and inequality
		//	shared_ptr<QuadProgMosek> program_ = make_shared <QuadProgMosek>();
		//	program_->setParamInt(MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT);
		//	program_->setParamInt(MSK_IPAR_LOG, 10);
		//	program_->setParamInt(MSK_IPAR_LOG_FILE, 1);
		//	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_DFEAS, 1e-8);
		//	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_INFEAS, 1e-10);
		//	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_MU_RED, 1e-8);
		//	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_NEAR_REL, 1e3);
		//	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_PFEAS, 1e-8);
		//	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_REL_GAP, 1e-8);
		//	program_->setNumberOfVariables(nr);

		//	program_->setObjectiveMatrix(MDKr_sp);
		//	program_->setObjectiveVector(-fr_);
		//	program_->setNumberOfInequalities(ni);
		//	program_->setInequalityMatrix(C.sparseView());
		//	program_->setNumberOfEqualities(ne);
		//	VectorXd cvec(ni);
		//	cvec.setZero();

		//	program_->setInequalityVector(cvec);
		//	program_->setEqualityMatrix(G_sp);

		//	VectorXd gvec(ne);
		//	gvec.setZero();
		//	program_->setEqualityVector(rhsG);

		//	bool success = program_->solve();
		//	VectorXd sol = program_->getPrimalSolution();
		//	qdot1 = sol.segment(0, nr);
		//}

		qddot = (qdot1 - qdot0) / h;
		q1 = q0 + h * qdot1;
		//cout << "ddot" << qddot << endl;
		//cout << "qdot1" << qdot1 << endl;

		yk.segment(0, nr) = q1;
		yk.segment(nr, nr) = qdot1;

		ydotk.segment(0, nr) = qdot1;
		ydotk.segment(nr, nr) = qddot;

		joint0->scatterDofs(yk, nr);
		joint0->scatterDDofs(ydotk, nr);

		deformable0->scatterDofs(yk, nr);
		deformable0->scatterDDofs(ydotk, nr);

		softbody0->scatterDofs(yk, nr);
		softbody0->scatterDDofs(ydotk, nr);

		//Energy ener = m_world->computeEnergy();
		/*cout << "V" << ener.V << endl;
		cout << "K" << ener.K << endl;
		cout << " sum " << ener.V + ener.K << endl;*/
		return yk;
	}
	break;

	case REDUCED_ODE45:
		break;
	case REDMAX_ODE45:
		break;
	default:
		break;
	}
}

