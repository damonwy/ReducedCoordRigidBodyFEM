#include "Solver.h"
#include "MatlabDebug.h"
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
#include "ChronoTimer.h"
#include <iostream>
#include <fstream>
#include <json.hpp>

using namespace std;
using namespace Eigen;


Solver::Solver()
{
	m_solutions = make_shared<Solution>();
}

Solver::Solver(shared_ptr<World> world, Integrator integrator) :
	m_world(world),
	m_integrator(integrator)
{
	m_solutions = make_shared<Solution>();
}

void Solver::init() {

}

void Solver::load(const string &RESOURCE_DIR) {

}

void Solver::reset() {
	int nr = m_world->nr;
	int nm = m_world->nm;
	// constraints
	int nem = m_world->nem;
	int ner = m_world->ner;
	int ne = nem + ner;

}

void Solver::initMatrix(int nm, int nr, int nem, int ner, int nim, int nir) {
	int ne = ner + nem;
	int ni = nim + nir;

	Mm.resize(nm, nm);
	Mm.setZero();
	Mr.resize(nr, nr);
	Mr.setZero();
	MDKr_.resize(nr, nr);
	MDKr_.setZero();

	fm.resize(nm);
	fr.resize(nr);
	fr_.resize(nr);
	fm.setZero();
	fr.setZero();
	fr_.setZero();
	tmp.resize(nm);
	tmp.setZero();

	Kr.resize(nr, nr);
	Dr.resize(nr, nr);
	Kr.setZero();
	Dr.setZero();
	Dm.resize(nm, nm);
	Dm.setZero();

	K.resize(nm, nm);
	K.setZero();
	Km.resize(nm, nm);
	Km.setZero();

	J.resize(nm, nr);
	Jdot.resize(nm, nr);
	J.setZero();
	Jdot.setZero();

	Gm.resize(nem, nm);
	Gm.setZero();
	Gmdot.resize(nem, nm);
	Gmdot.setZero();

	gm.resize(nem);
	gm.setZero();
	gmdot.resize(nem);
	gmdot.setZero();
	gmddot.resize(nem);
	gmddot.setZero();

	Gr.resize(ner, nr);
	Gr.setZero();
	Grdot.resize(ner, nr);
	Grdot.setZero();

	gr.resize(ner);
	gr.setZero();
	grdot.resize(ner);
	grdot.setZero();
	grddot.resize(ner);
	grddot.setZero();

	G.resize(ne, nr);
	g.resize(ne);
	rhsG.resize(ne);
	gdot.resize(ne);
	G.setZero();
	g.setZero();
	gdot.setZero();
	rhsG.setZero();

	Cm.resize(nim, nm);
	Cm.setZero();
	Cmdot.resize(nim, nm);
	Cmdot.setZero();
	cm.resize(nim);
	cm.setZero();
	cmdot.resize(nim);
	cmdot.setZero();
	cmddot.resize(nim);
	cmddot.setZero();

	Cr.resize(nir, nr);
	Cr.setZero();
	Crdot.resize(nir, nr);
	Crdot.setZero();
	cr.resize(nir);
	cr.setZero();
	crdot.resize(nir);
	crdot.setZero();
	crddot.resize(nir);
	crddot.setZero();
}

void Solver::initMatrixSparse(int nm, int nr, int nem, int ner, int nim, int nir) {
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
	K_sp.setZero();
	Km_sp.resize(nm, nm);
	Km_sp.setZero();

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

	//G_sp.resize(ne, nr);
	//G_sp.setZero();
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

Eigen::VectorXd Solver::dynamics(Eigen::VectorXd y)
{
	ChronoTimer timer("Test");

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

		// sceneFcn()
		body0->computeMassGrav(grav, Mm, fm);
		body0->computeForceDamping(tmp, Dm);

		timer.toc();
		timer.print();

		deformable0->computeMass(grav, Mm, fm);
		deformable0->computeForceDamping(grav, tmp, Dm);
		timer.toc();
		timer.print();

		softbody0->computeMass(grav, Mm);
		softbody0->computeForce(grav, fm);
		softbody0->computeStiffness(K);
		timer.toc();
		timer.print();

		joint0->computeForceStiffness(fr, Kr);
		joint0->computeForceDamping(tmp, Dr);
		joint0->computeJacobian(J, Jdot);

		timer.toc();
		timer.print();

		deformable0->computeJacobian(J, Jdot);
		softbody0->computeJacobian(J);
		
		timer.toc();
		timer.print();
		
		spring0->computeForceStiffnessDamping(fm, Km, Dm);

		timer.toc();
		timer.print();

		q0 = y.segment(0, nr);
		qdot0 = y.segment(nr, nr);

		Mr = J.transpose() * (Mm - h * h * K) * J;
		Mr = 0.5 * (Mr + Mr.transpose());

		timer.toc();
		timer.print();

		fr_ = Mr * qdot0 + h * (J.transpose() * (fm - Mm * Jdot * qdot0) + fr);
		MDKr_ = Mr + J.transpose() * (h * Dm - h * h * Km)*J + h * Dr - h * h * Kr;
		cout << "MDKr_" << MDKr_ << endl;
		cout << "fr_" << fr_ << endl;

		if (ne > 0) {		
			constraint0->computeJacEqM(Gm, Gmdot, gm, gmdot, gmddot);
			constraint0->computeJacEqR(Gr, Grdot, gr, grdot, grddot);
			G.block(0, 0, nem, nr) = Gm * J;
			G.block(nem, 0, ner, nr) = Gr;
			g.segment(0, nem) = gm;
			g.segment(nem, ner) = gr;
			rhsG = -gdot - 100.0 * g;// todo!!!!!
		}

		if (ni > 0) {
			// Check for active inequality constraint
			constraint0->computeJacIneqM(Cm, Cmdot, cm, cmdot, cmddot);
			constraint0->computeJacIneqR(Cr, Crdot, cr, crdot, crddot);
			rowsR.clear();
			rowsM.clear();

			constraint0->getActiveList(rowsM, rowsR);
			nim = rowsM.size();
			nir = rowsR.size();
			ni = nim + nir;

			if (ni > 0) {
				Eigen::VectorXi m_rowsM = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(rowsM.data(), rowsM.size());
				Eigen::VectorXi m_rowsR = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(rowsR.data(), rowsR.size());

				MatrixXd m_Cm = Cm(m_rowsM, Eigen::placeholders::all);
				MatrixXd m_Cr = Cr(m_rowsR, Eigen::placeholders::all);
				VectorXd m_cm = cm(m_rowsM);
				VectorXd m_cr = cr(m_rowsR);
				VectorXd m_cmdot = cmdot(m_rowsM);
				VectorXd m_crdot = crdot(m_rowsR);

				MatrixXd CmJ = m_Cm * J;
				C.resize(CmJ.rows() + m_Cr.rows(), m_Cr.cols());
				C << CmJ, m_Cr;
				rhsC.resize(C.rows());
				VectorXd c(C.rows());
				c << m_cm, m_cr;
				VectorXd cdot(C.rows());
				cdot << m_cmdot, m_crdot;
				rhsC = -cdot - 105.0 * c;		
			}
		}

		if (ne == 0 && ni == 0) {	// No constraints	
			qdot1 = MDKr_.ldlt().solve(fr_);
		}
		else if (ne > 0 && ni == 0) {  // Just equality
			int rows = MDKr_.rows() + G.rows();
			int cols = MDKr_.cols() + G.rows();
			MatrixXd LHS(rows, cols);
			VectorXd rhs(rows);
			LHS.setZero();
			rhs.setZero();
			LHS.block(0, 0, MDKr_.rows(), MDKr_.cols()) = MDKr_;
			LHS.block(0, MDKr_.cols(), MDKr_.rows(), G.rows()) = G.transpose();
			LHS.block(MDKr_.rows(), 0, G.rows(), G.cols()) = G;
			rhs.segment(0, fr_.rows()) = fr_;	
			rhs.segment(fr_.rows(), g.rows()) = rhsG;

			VectorXd sol = LHS.ldlt().solve(rhs);
			qdot1 = sol.segment(0, nr);

			VectorXd l = sol.segment(nr, sol.rows() - nr);

			constraint0->scatterForceEqM(Gm.transpose(), l.segment(0, nem) / h);
			constraint0->scatterForceEqR(Gr.transpose(), l.segment(nem, l.rows() - nem) / h);

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
			program_->setObjectiveMatrix(MDKr_.sparseView());
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
		else {  // Both equality and inequality
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

			program_->setObjectiveMatrix(MDKr_.sparseView());
			program_->setObjectiveVector(-fr_);
			program_->setNumberOfInequalities(ni);
			program_->setInequalityMatrix(C.sparseView());
			program_->setNumberOfEqualities(ne);
			VectorXd cvec(ni);
			cvec.setZero();

			program_->setInequalityVector(cvec);
			program_->setEqualityMatrix(G.sparseView());

			VectorXd gvec(ne);
			gvec.setZero();
			program_->setEqualityVector(rhsG);

			bool success = program_->solve();
			VectorXd sol = program_->getPrimalSolution();
			qdot1 = sol.segment(0, nr);
		}

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

		Energy ener = m_world->computeEnergy();
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

Eigen::VectorXd Solver::dynamics_sparse(Eigen::VectorXd y)
{
	ChronoTimer timer("Test");
	Eigen::SparseMatrix<double, Eigen::RowMajor> G_sp;
	
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

		initMatrixSparse(nm, nr, nem, ner, nim, nir);

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

		// sceneFcn()
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

				J_.push_back(T(i, j, J_dense(i,j)));
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
		Kr_sp.setFromTriplets(Kr_.begin(), Kr_.end());
		J_sp.setFromTriplets(J_.begin(), J_.end());
		Jdot_sp.setFromTriplets(Jdot_.begin(), Jdot_.end());

		q0 = y.segment(0, nr);
		qdot0 = y.segment(nr, nr);

		Mr_sp = J_sp.transpose() * (Mm_sp - h * h * K_sp) * J_sp;
		//Mr_sp = 0.5 * (Mr_sp_temp + Mr_sp_temp.transpose());

		//Mr_sp = 0.5 * (Mr_sp + Mr_sp.transpose());

		fr_ = Mr_sp * qdot0 + h * (J_sp.transpose() * (fm - Mm_sp * Jdot_sp * qdot0) + fr);
		MDKr_sp = Mr_sp + J_sp.transpose() * (h * Dm_sp - h * h * Km_sp)*J_sp + h * Dr_sp - h * h * Kr_sp;

		if (ne > 0) {
			constraint0->computeJacEqMSparse(Gm_, Gmdot_, gm, gmdot, gmddot);
			constraint0->computeJacEqRSparse(Gr_, Grdot_, gr, grdot, grddot);

			Gm_sp.setFromTriplets(Gm_.begin(), Gm_.end());
			Gmdot_sp.setFromTriplets(Gmdot_.begin(), Gmdot_.end());
			Gr_sp.setFromTriplets(Gr_.begin(), Gr_.end());
			Grdot_sp.setFromTriplets(Grdot_.begin(), Grdot_.end());

			sparse_to_file_as_dense(Gm_sp * J_sp, "Gm_sp * J_sp");
			G_sp.topRows(nem) = Gm_sp * J_sp;
			G_sp.bottomRows(ner) = Gr_sp;
			sparse_to_file_as_dense(G_sp, "G_sp");
			
			g.segment(0, nem) = gm;
			g.segment(nem, ner) = gr;
			rhsG = -gdot - 100.0 * g;// todo!!!!!
		}

		//if (ni > 0) {
		//	// Check for active inequality constraint
		//	constraint0->computeJacIneqM(Cm, Cmdot, cm, cmdot, cmddot);
		//	constraint0->computeJacIneqR(Cr, Crdot, cr, crdot, crddot);
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
			sparse_to_file_as_dense(MDKr_sp, "MDKr_sp");
			cout << "fr_" << fr_ << endl;
			qdot1 = cg.solveWithGuess(fr_, qdot0);
			cout << "qdot1" << qdot1 << endl;
		}
		else if (ne > 0 && ni == 0) {  // Just equality
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
			sparse_to_file_as_dense(MDKr_sp, "MDKr_sp");

			program_->setObjectiveVector(-fr_);
			program_->setNumberOfEqualities(ne);
sparse_to_file_as_dense(G_sp, "G_sp");
			program_->setEqualityMatrix(G_sp);
			
			VectorXd gvec(ne);
			gvec.setZero();
			program_->setEqualityVector(rhsG);

			bool success = program_->solve();
			VectorXd sol = program_->getPrimalSolution();
			qdot1 = sol.segment(0, nr);
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
		cout << "ddot" << qddot << endl;
		////cout << "qdot1" << qdot1 << endl;

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


shared_ptr<Solution> Solver::solve() {
	switch (m_integrator)
	{
	case REDMAX_EULER:
	{
		int nr = m_world->nr;
		int nm = m_world->nm;

		int nem = m_world->nem;
		int ner = m_world->ner;
		int ne = nem + ner;

		auto body0 = m_world->getBody0();
		auto joint0 = m_world->getJoint0();
		auto deformable0 = m_world->getDeformable0();
		auto softbody0 = m_world->getSoftBody0();
		auto constraint0 = m_world->getConstraint0();

		int nsteps = m_world->getNsteps();
		m_solutions->t.resize(nsteps);
		m_solutions->y.resize(nsteps, 2 * nr);
		m_solutions->y.setZero();

		// initial state
		m_solutions->t(0) = m_world->getTspan()(0);
		m_solutions->y.row(0) = joint0->gatherDofs(m_solutions->y.row(0), nr);
		VectorXd sol_y = m_solutions->y.row(0);
		deformable0->gatherDofs(sol_y, nr);
		m_solutions->y.row(0) = sol_y;
		m_solutions->y.row(0) = softbody0->gatherDofs(m_solutions->y.row(0), nr);

		double t = m_world->getTspan()(0);
		double h = m_world->getH();
		Vector3d grav = m_world->getGrav();

		VectorXd yk(2 * nr);
		VectorXd ydotk(2 * nr);

		for (int k = 1; k < nsteps; k++) {
			int nim = m_world->nim;
			int nir = m_world->nir;
			int ni = nim + nir;
			initMatrix(nm, nr, nem, ner, nim, nir);
			// sceneFcn()
			body0->computeMassGrav(grav, Mm, fm);
			deformable0->computeMass(grav, Mm, fm);

			softbody0->computeMass(grav, Mm);
			softbody0->computeForce(grav, fm);
			softbody0->computeStiffness(K);

			joint0->computeForceStiffness(fr, Kr);
			joint0->computeForceDamping(tmp, Dr);

			joint0->computeJacobian(J, Jdot);
			//Jdot = joint0->computeJacobianDerivative(Jdot, J, nm, nr);
			// spring jacobian todo
			deformable0->computeJacobian(J, Jdot);

			softbody0->computeJacobian(J);

			q0 = m_solutions->y.row(k - 1).segment(0, nr);
			qdot0 = m_solutions->y.row(k - 1).segment(nr, nr);

			Mr = J.transpose() * (Mm - h * h * K) * J;
			Mr = 0.5 * (Mr + Mr.transpose());

			fr_ = Mr * qdot0 + h * (J.transpose() * (fm - Mm * Jdot * qdot0) + fr);
			MDKr_ = Mr + J.transpose() * (h * Dm - h * h * Km)*J + h * Dr - h * h * Kr;

			if (ne > 0) {
				constraint0->computeJacEqM(Gm, Gmdot, gm, gmdot, gmddot);
				constraint0->computeJacEqR(Gr, Grdot, gr, grdot, grddot);
				G.block(0, 0, nem, nr) = Gm * J;
				G.block(nem, 0, ner, nr) = Gr;
				g.segment(0, nem) = gm;
				g.segment(nem, ner) = gr;
				rhsG = -gdot - 5.0 * g;// todo!!!!!

			}

			if (ni > 0) {
				// Check for active inequality constraint
				constraint0->computeJacIneqM(Cm, Cmdot, cm, cmdot, cmddot);
				constraint0->computeJacIneqR(Cr, Crdot, cr, crdot, crddot);
				rowsR.clear();
				rowsM.clear();

				constraint0->getActiveList(rowsM, rowsR);
				nim = rowsM.size();
				nir = rowsR.size();
				ni = nim + nir;

				if (ni > 0) {
					Eigen::VectorXi m_rowsM = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(rowsM.data(), rowsM.size());
					Eigen::VectorXi m_rowsR = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(rowsR.data(), rowsR.size());
					MatrixXd m_Cm = Cm(m_rowsM, Eigen::placeholders::all);
					MatrixXd m_Cr = Cr(m_rowsR, Eigen::placeholders::all);
					MatrixXd CmJ = m_Cm * J;
					C.resize(CmJ.rows() + m_Cr.rows(), m_Cr.cols());
					C << CmJ, m_Cr;
				}
			}

			if (ne == 0 && ni == 0) {	// No constraints	
				qdot1 = MDKr_.ldlt().solve(fr_);
			}
			else if (ne > 0 && ni == 0) {  // Just equality
				int rows = MDKr_.rows() + G.rows();
				int cols = MDKr_.cols() + G.rows();
				MatrixXd LHS(rows, cols);
				VectorXd rhs(rows);
				LHS.setZero();
				rhs.setZero();
				LHS.block(0, 0, MDKr_.rows(), MDKr_.cols()) = MDKr_;
				LHS.block(0, MDKr_.cols(), MDKr_.rows(), G.rows()) = G.transpose();
				LHS.block(MDKr_.rows(), 0, G.rows(), G.cols()) = G;
				rhs.segment(0, fr_.rows()) = fr_;
				//rhs.segment(ftilde.rows(), g.rows()) = g;
				//cout << "g" << endl << g << endl;
				rhs.segment(fr_.rows(), g.rows()) = rhsG;

				VectorXd sol = LHS.ldlt().solve(rhs);
				//cout << LHS << endl;
				//cout << rhs << endl;
				qdot1 = sol.segment(0, nr);

				VectorXd l = sol.segment(nr, sol.rows() - nr);

				constraint0->scatterForceEqM(Gm.transpose(), l.segment(0, nem) / h);
				constraint0->scatterForceEqR(Gr.transpose(), l.segment(nem, l.rows() - nem) / h);

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
				program_->setObjectiveMatrix(MDKr_.sparseView());
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
			else {  // Both equality and inequality
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

				program_->setObjectiveMatrix(MDKr_.sparseView());
				program_->setObjectiveVector(-fr_);
				program_->setNumberOfInequalities(ni);
				program_->setInequalityMatrix(C.sparseView());

				VectorXd cvec(ni);
				cvec.setZero();

				program_->setInequalityVector(cvec);
				program_->setEqualityMatrix(G.sparseView());

				VectorXd gvec(ne);
				gvec.setZero();
				program_->setEqualityVector(gvec);

				bool success = program_->solve();
				VectorXd sol = program_->getPrimalSolution();
				qdot1 = sol.segment(0, nr);

			}
			qddot = (qdot1 - qdot0) / h;
			cout << "ddot" << qddot << endl;
			cout <<"qdot1"<< qdot1 << endl;
			q1 = q0 + h * qdot1;
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

			t += h;
			m_solutions->y.row(k) = yk;
			m_solutions->t(k) = t;

		}
		return m_solutions;
		break;
	}

	case REDUCED_ODE45:
		break;
	case REDMAX_ODE45:
		break;
	default:
		break;
	}
}