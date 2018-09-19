#include "Solver.h"

#include "World.h"
#include "Body.h"
#include "Joint.h"
#include "Spring.h"
#include "SpringSerial.h"
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

shared_ptr<Solution> Solver::solve() {
	switch (m_integrator)
	{
	case REDMAX_EULER:
		{
			int nr = m_world->nr;
			int nm = m_world->nm;

			M.resize(nm, nm);
			M.setZero();
			f.resize(nm);
			f.setZero();
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

			Gr.resize(ner, nr);
			Gr.setZero();
			Grdot.resize(ner, nr);
			Grdot.setZero();

			gr.resize(ner);
			gr.setZero();

			G.resize(ne, nr);
			g.resize(ne);
			G.setZero();
			g.setZero();

			auto body0 = m_world->getBody0();
			auto joint0 = m_world->getJoint0();
			auto spring0 = m_world->getSpring0();
			auto constraint0 = m_world->getConstraint0();

			int nsteps = m_world->getNsteps();
			m_solutions->t.resize(nsteps);
			m_solutions->y.resize(nsteps, 2 * nr);
			m_solutions->y.setZero();

			// initial state
			m_solutions->t(0) = m_world->getTspan()(0);
			m_solutions->y.row(0) = joint0->gatherDofs(m_solutions->y.row(0), nr);
			spring0->gatherDofs(m_solutions->y.row(0), nr);

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

				// sceneFcn()
				M = body0->computeMass(grav, M);
				f = body0->computeForce(grav, f);

				// spring..
				J = joint0->computeJacobian(J, nm, nr);	
				Jdot = joint0->computeJacobianDerivative(Jdot, J, nm, nr);
				
				// spring jacobian todo
				
				q0 = m_solutions->y.row(k - 1).segment(0, nr);
				qdot0 = m_solutions->y.row(k - 1).segment(nr, nr);
				
				Mtilde = J.transpose() * M * J;
				Mtilde = 0.5 * (Mtilde + Mtilde.transpose());
				ftilde = Mtilde * qdot0 + h * J.transpose() * (f - M * Jdot * qdot0);
				
				if (ne > 0) {
					constraint0->computeJacEqM(Gm, Gmdot, gm);
					constraint0->computeJacEqR(Gr, Grdot, gr);
					G.block(0, 0, nem, nr) = Gm * J;
					G.block(nem, 0, ner, nr) = Gr;
					g.segment(0, nem) = gm;
					g.segment(nem, ner) = gr;
				}

				if (ni > 0) {
					// Check for active inequality constraint
					constraint0->computeJacIneqM(Cm, Cmdot, cm);
					constraint0->computeJacIneqR(Cr, Crdot, cr);
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
					qdot1 = Mtilde.ldlt().solve(ftilde);

				}
				else if (ne > 0 && ni == 0) {  // Just equality
					int rows = Mtilde.rows() + G.rows();
					int cols = Mtilde.cols() + G.rows();
					MatrixXd LHS(rows, cols);
					VectorXd rhs(rows);
					LHS.setZero();
					rhs.setZero();
					LHS.block(0, 0, Mtilde.rows(), Mtilde.cols()) = Mtilde;
					LHS.block(0, Mtilde.cols(), Mtilde.rows(), G.rows()) = G.transpose();
					LHS.block(Mtilde.rows(), 0, G.rows(), G.cols()) = G;
					rhs.segment(0, ftilde.rows()) = ftilde;
					rhs.segment(ftilde.rows(), g.rows()) = g;

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
					program_->setObjectiveMatrix(Mtilde.sparseView());
					program_->setObjectiveVector(-ftilde);
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

					program_->setObjectiveMatrix(Mtilde.sparseView());
					program_->setObjectiveVector(-ftilde);
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
				q1 = q0 + h * qdot1;

				yk.segment(0, nr) = q1;
				yk.segment(nr, nr) = qdot1;

				ydotk.segment(0, nr) = qdot1;
				ydotk.segment(nr, nr) = qddot;

				joint0->scatterDofs(yk, nr);
				joint0->scatterDDofs(ydotk, nr);

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

Solver::~Solver() {
}