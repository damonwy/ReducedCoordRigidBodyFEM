#include "Solver.h"
#include "World.h"

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

void Solver::reset() {
	int nr = m_world->nr;
	int nm = m_world->nm;
	// constraints
	int nem = m_world->nem;
	int ner = m_world->ner;
	int ne = nem + ner;
}
