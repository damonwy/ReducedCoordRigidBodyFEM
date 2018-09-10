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

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

#include <unsupported/Eigen/MatrixFunctions> // TODO: avoid using this later, write a func instead

Scene::Scene() :
	t(0.0),
	h(1e-2),
	grav(0.0, 0.0, 0.0)
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

	world = make_shared<World>(SERIAL_CHAIN);
	world->load(RESOURCE_DIR);

}


void Scene::init()
{
	world->init();

}

void Scene::reset()
{
	
}

void Scene::step()
{	
	world->update();
	t += h;
}


void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> prog2, shared_ptr<MatrixStack> P) const
{
	world->draw(MV, prog, P);
	
}
