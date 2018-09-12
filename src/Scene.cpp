#include "Scene.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

#include "Particle.h"
#include "Shape.h"
#include "Program.h"
#include "Joint.h"
#include "MatlabDebug.h"
#include "Vector.h"
#include "JsonEigen.h"
#include "Stepper.h"
#include "World.h"
#include "Solver.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

#include <unsupported/Eigen/MatrixFunctions> // TODO: avoid using this later, write a func instead

Scene::Scene() :
	t(0.0),
	h(1e-2),
	grav(0.0, 0.0, 0.0),
	time_step(0)
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

	m_world = make_shared<World>(SERIAL_CHAIN);
	m_world->load(RESOURCE_DIR);

	m_solver = make_shared<Solver>(m_world, REDMAX_EULER);
	m_solver->load(RESOURCE_DIR);

}


void Scene::init()
{
	m_world->init();
	m_solver->init();
	m_solution = m_solver->solve();
}

void Scene::reset()
{
	
}

void Scene::solve() {
	

}

void Scene::step()
{	
	//world->update();
	if (time_step < m_solution->getNsteps()) {
		m_solution->step(time_step);
		time_step++;
		t += h;
	}	
}


void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const
{
	m_world->draw(MV, prog, progSimple, P);
	
}
