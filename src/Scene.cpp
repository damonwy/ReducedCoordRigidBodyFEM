#include "Scene.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

#include "Node.h"
#include "Shape.h"
#include "Program.h"
#include "Joint.h"
#include "MatlabDebug.h"
#include "Vector.h"
#include "JsonEigen.h"
#include "Stepper.h"
#include "World.h"
#include "Solver.h"
#include "Spring.h"
#include "SpringSerial.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

#include <unsupported/Eigen/MatrixFunctions> // TODO: avoid using this later, write a func instead

Scene::Scene() :
	t(0.0),
	h(1e-2),
	grav(0.0, 0.0, 0.0),
	time_step(0),
	drawHz(10)
{

}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{	

	//read a JSON file
	ifstream i(RESOURCE_DIR + "input.json");
	i >> js;
	i.close();

	// Units: meters, kilograms, seconds
	h = js["h"];
	Eigen::from_json(js["grav"], grav);
	drawHz = js["drawHz"];

	m_world = make_shared<World>(SOFT_BODIES);
	m_world->load(RESOURCE_DIR);

	m_solver = make_shared<Solver>(m_world, REDMAX_EULER);
	m_solver->load(RESOURCE_DIR);

}


void Scene::init()
{
	m_world->init();
	m_solver->init();
	m_solution = m_solver->solve();

	tk = m_solution->t(0);
	drawH = 1.0 / drawHz;
	search_idx = 0;
}

void Scene::reset()
{
	
}

void Scene::solve() {
	

}

void Scene::step()
{	
	int n_steps = m_solution->getNsteps();
	
	int output_idx;
	double s;
	VectorXd ys;

	if(tk < m_solution->t(n_steps-1)) {
		m_solution->searchTime(tk, search_idx, output_idx, s);
		search_idx = output_idx;
		ys = (1 - s)* m_solution->y.row(output_idx) + s * m_solution->y.row(output_idx + 1);

		m_world->getJoint0()->scatterDofs(ys, m_world->nr);
		m_world->getSpring0()->scatterDofs(ys, m_world->nr);

		tk = tk + drawH;
	}
	
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const
{
	m_world->draw(MV, prog, progSimple, P);
	
}
