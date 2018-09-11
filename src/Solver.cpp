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
			auto body0 = m_world->getBody0();
			auto joint0 = m_world->getJoint0();
			// auto spring0 = m_world->getSpring0();
			// auto constraint0 = m_world->getConstraint0();

			int nsteps = m_world->getNsteps();
			m_solutions->t.resize(nsteps);
			m_solutions->y.resize(nsteps, 2 * nr);

			// initial state
			m_solutions->y(0) = 


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


				t += h;
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