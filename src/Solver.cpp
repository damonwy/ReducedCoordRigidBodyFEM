#include "Solver.h"

#include "World.h"
#include "Body.h"
#include "Joint.h"
#include "ConstraintJointLimit.h"


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

			int ni = 0;

			auto body0 = m_world->getBody0();
			auto joint0 = m_world->getJoint0();
			// auto spring0 = m_world->getSpring0();
			auto constraint0 = m_world->getConstraint0();

			int nsteps = m_world->getNsteps();
			m_solutions->t.resize(nsteps);
			m_solutions->y.resize(nsteps, 2 * nr);
			m_solutions->y.setZero();

			// initial state
			m_solutions->t(0) = m_world->getTspan()(0);
			m_solutions->y.row(0) = joint0->gatherDofs(m_solutions->y.row(0), nr);
			// m_solutions->y.row(0) = spring0->gatherDofs(m_solutions->y.row(0));

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
				//cout << "k" << k << endl << endl;
				M = body0->computeMass(grav, M);

				//cout << "M" << M << endl << endl; 

				f = body0->computeForce(grav, f);
				//cout << "f" << f << endl << endl;

				// spring..
				J = joint0->computeJacobian(J, nm, nr);
				//cout << "J" << J << endl << endl;
				
				Jdot = joint0->computeJacobianDerivative(Jdot, J, nm, nr);
				//cout << "JD" << Jdot << endl << endl;
				
				// spring jacobian todo
				
				q0 = m_solutions->y.row(k - 1).segment(0, nr);
				qdot0 = m_solutions->y.row(k - 1).segment(nr, nr);
				//cout << "q" << q0 << endl << endl;
				//cout << "qdot0" << qdot0 << endl << endl;

				Mtilde = J.transpose() * M * J;
				
				Mtilde = 0.5 * (Mtilde + Mtilde.transpose());
				//cout << "Mtilde" << Mtilde << endl << endl;

				ftilde = Mtilde * qdot0 + h * J.transpose() * (f - M * Jdot * qdot0);
				//cout << "ftilde" << ftilde << endl << endl;
				
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


				}
				// Solve 
				if (ne == 0 && ni == 0) {	// No constraints	
					qdot1 = Mtilde.ldlt().solve(ftilde);
					//cout << "qdot1" << qdot1 << endl << endl;  

				}
				else if (ne > 0 && ni == 0) {  // Just equality
					// todo
					// kkt
				}
				else if (ne == 0 && ni > 0) {  // Just inequality
					// todo
					// qp
				}
				else {  // Both equality and inequality
					// todo
					// qp
				}

				qddot = (qdot1 - qdot0) / h;
				//cout << "qddot" << qddot << endl << endl;

				q1 = q0 + h * qdot1;
				//cout << "q1" << q1 << endl << endl;

				yk.segment(0, nr) = q1;
				yk.segment(nr, nr) = qdot1;

				//cout << "yk" << yk << endl << endl;

				ydotk.segment(0, nr) = qdot1;
				ydotk.segment(nr, nr) = qddot;
				//cout << "ydotk" << ydotk << endl << endl;

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