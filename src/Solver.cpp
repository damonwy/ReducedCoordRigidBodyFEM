#include "Solver.h"

#include "World.h"
#include "Body.h"
#include "SoftBody.h"
#include "Joint.h"
#include "Deformable.h"
#include "DeformableSpring.h"
#include "ConstraintJointLimit.h"
#include "ConstraintLoop.h"
#include "ConstraintAttachSpring.h"
#include "QuadProgMosek.h"

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


Eigen::VectorXd Solver::dynamics(Eigen::VectorXd y)
{
	switch (m_integrator)
	{
	case REDMAX_EULER:
	{
		int nr = m_world->nr;
		int nm = m_world->nm;

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

		K.resize(nm, nm);
		K.setZero();
		Km.resize(nm, nm);
		Km.setZero();

		J.resize(nm, nr);
		Jdot.resize(nm, nr);
		J.setZero();
		Jdot.setZero();

		// constraints
		int nem = m_world->nem;
		int ner = m_world->ner;
		int ne = nem + ner;

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

		auto body0 = m_world->getBody0();
		auto joint0 = m_world->getJoint0();
		auto deformable0 = m_world->getDeformable0();
		auto softbody0 = m_world->getSoftBody0();
		auto constraint0 = m_world->getConstraint0();

		double t = m_world->getTspan()(0);
		double h = m_world->getH();
		Vector3d grav = m_world->getGrav();

		VectorXd yk(2 * nr);
		VectorXd ydotk(2 * nr);

		int nim = m_world->nim;
		int nir = m_world->nir;
		int ni = nim + nir;

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

		Mm.setZero();
		K.setZero();
		fm.setZero();
		J.setZero();
		Jdot.setZero();
		g.setZero();
		gdot.setZero();
		gm.setZero();
		gmdot.setZero();
		gmddot.setZero();
		gr.setZero();
		grdot.setZero();
		grddot.setZero();
		Gr.setZero();
		Grdot.setZero();
		Gm.setZero();
		Gmdot.setZero();
		G.setZero();
		rhsG.setZero();
		tmp.resize(nm);
		Dm.resize(nm, nm);
		tmp.setZero();
		Dm.setZero();

		// sceneFcn()
		body0->computeMassGrav(grav, Mm, fm);
		body0->computeForceDamping(tmp, Dm);

		deformable0->computeMass(grav, Mm, fm);
		deformable0->computeForceDamping(grav, tmp, Dm);

		softbody0->computeMass(grav, Mm);
		softbody0->computeForce(grav, fm);
		softbody0->computeStiffness(K);

		joint0->computeForceStiffness(fr, Kr);
		joint0->computeForceDamping(tmp, Dr);
		joint0->computeJacobian(J, Jdot, nm, nr);

		deformable0->computeJacobian(J, Jdot);
		softbody0->computeJacobian(J);

		q0 = y.segment(0, nr);
		qdot0 = y.segment(nr, nr);

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
			//cout << qdot1 << endl;

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

shared_ptr<Solution> Solver::solve() {
	switch (m_integrator)
	{
	case REDMAX_EULER:
	{
		int nr = m_world->nr;
		int nm = m_world->nm;

		Mm.resize(nm, nm);
		Mm.setZero();
		MDKr_.resize(nr, nr);
		MDKr_.setZero();
		fr_.resize(nr, 1);
		fr_.setZero();
		fr.resize(nr, 1);
		fr.setZero();
		Kr.resize(nr, nr);
		Kr.setZero();
		Dr.resize(nr, nr);
		Dr.setZero();

		K.resize(nm, nm);
		K.setZero();
		fm.resize(nm);
		fm.setZero();
		J.resize(nm, nr);
		Jdot.resize(nm, nr);
		J.setZero();
		Jdot.setZero();

		// constraints
		int nem = m_world->nem;
		int ner = m_world->ner;
		int ne = nem + ner;

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

		G.resize(ne, nr);
		g.resize(ne);
		rhsG.resize(ne);
		gdot.resize(ne);
		G.setZero();
		g.setZero();
		gdot.setZero();
		rhsG.setZero();

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

			Cm.resize(nim, nm);
			Cm.setZero();
			Cmdot.resize(nim, nm);
			Cmdot.setZero();
			cm.resize(nim);
			cm.setZero();

			Cr.resize(nir, nr);
			Cr.setZero();
			Crdot.resize(nir, nr);
			Crdot.setZero();
			cr.resize(nir);
			cr.setZero();

			Mm.setZero();
			K.setZero();
			fm.setZero();
			J.setZero();
			Jdot.setZero();
			g.setZero();
			gdot.setZero();
			gm.setZero();
			gmdot.setZero();
			gmddot.setZero();
			gr.setZero();
			Gr.setZero();
			Grdot.setZero();
			Gm.setZero();
			Gmdot.setZero();
			G.setZero();
			rhsG.setZero();
			// sceneFcn()
			body0->computeMassGrav(grav, Mm, fm);
			deformable0->computeMass(grav, Mm, fm);

			softbody0->computeMass(grav, Mm);
			softbody0->computeForce(grav, fm);

			softbody0->computeStiffness(K);

			joint0->computeForceStiffness(fr, Kr);
			joint0->computeForceDamping(tmp, Dr);

			joint0->computeJacobian(J, Jdot, nm, nr);
			//Jdot = joint0->computeJacobianDerivative(Jdot, J, nm, nr);
			// spring jacobian todo
			deformable0->computeJacobian(J, Jdot);

			softbody0->computeJacobian(J);

			q0 = m_solutions->y.row(k - 1).segment(0, nr);
			//cout << "q0"<<q0 << endl;
			qdot0 = m_solutions->y.row(k - 1).segment(nr, nr);
			//cout << "q0" << qdot0 << endl;
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

				//cout << Mtilde << endl;
				//cout << ftilde << endl;

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
				//cout << qdot1 << endl;

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
			//cout << "q1" << q1 << endl;
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