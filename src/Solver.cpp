#include "Solver.h"

#include "World.h"
#include "Body.h"
#include "Joint.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

using namespace std;
using namespace Eigen;


Solver::Solver() {

}

Solver::Solver(shared_ptr<World> world, Integrator integrator) :
m_world(world),
m_integrator(integrator)
{

}

shared_ptr<Solution> Solver::solve() {
	switch (m_integrator)
	{
	case REDMAX_EULER:
		{
			int nr = m_world->nr;
			int nm = m_world->nm;
			// constraints
			// int nem = m_world->nem;
			// int ner = m_world->ner;
			// int ne = nem + ner;
			int nem = 0;
			int ner = 0;
			int ne = nem + ner;
			int ni = 0;

			auto body0 = m_world->getBody0();
			auto joint0 = m_world->getJoint0();
			// auto spring0 = m_world->getSpring0();
			// auto constraint0 = m_world->getConstraint0();

			int nsteps = m_world->getNsteps();
			m_solutions->t.resize(nsteps);
			m_solutions->y.resize(nsteps, 2 * nr);
			m_solutions->y.setZero();

			// initial state

			m_solutions->y.row(0) = joint0->gatherDofs(m_solutions->y.row(0), nr);
			// m_solutions->y.row(0) = spring0->gatherDofs(m_solutions->y.row(0));

			double t = m_world->getTspan()(0);
			double h = m_world->getH();
			Vector3d grav = m_world->getGrav();

			for (int k = 1; k < nsteps; k++) {
				// sceneFcn()
				M = body0->computeMass(grav, M);
				f = body0->computeForce(grav, f);
				// spring..
				J = joint0->computeJacobian(J, nm, nr);
				Jdot = joint0->computeJacobianDerivative(Jdot, nm, nr);
				// spring jacobian todo
				
				q0 = m_solutions->y.row(k - 1).segment(0, nr);
				qdot0 = m_solutions->y.row(k - 1).segment(nr, nr);


				Mtilde = J.transpose() * M * J;
				Mtilde = 0.5 * (Mtilde + Mtilde.transpose());

				ftilde = Mtilde * qdot0 + h * J.transpose() * (f - M * Jdot * qdot0);

				// Solve 
				if (ne == 0 && ni == 0) {
					// No constraints
					qdot1 = Mtilde.ldlt().solve(ftilde);
				}
				else if (ne > 0 && ni == 0) {
					// Just equality
					// todo
					// kkt
				}
				else if (ne == 0 && ni > 0) {
					// Just inequality
					// todo
					// qp
				}
				else {
					// Both equality and inequality
					// todo
					// qp
				}

				qddot = (qdot1 - qdot0) / h;
				q1 = q0 + h * qdot1;
				
				VectorXd yk(2 * nr);
				VectorXd ydotk(2 * nr);
				yk.segment(0, nr) = q1;
				yk.segment(nr, nr) = qdot1;
				ydotk.segment(0, nr) = qdot1;
				ydotk.segment(nr, nr) = qddot;

				joint0->scatterDofs(yk, nr);
				joint0->scatterDDofs(ydotk, nr);

				t += h;
				m_solutions->y.row(k) = yk;
				m_solutions->t(k) = t;

				return m_solutions;
			}
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