

#include "Scene.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

#include "Particle.h"
#include "Shape.h"
#include "Program.h"
#include "Rigid.h"
#include "Joint.h"
#include "MatlabDebug.h"
#include "Vector.h"
#include "JsonEigen.h"
#include "Stepper.h"
#include "SymplecticIntegrator.h"

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
	// Init shapes
	boxShape = make_shared<Shape>();
	boxShape->loadMesh(RESOURCE_DIR + "box5.obj");

	//read a JSON file
	ifstream i(RESOURCE_DIR + "input.json");
	i >> js;
	i.close();

	// Units: meters, kilograms, seconds
	h = js["h"];
	Eigen::from_json(js["grav"], grav);

	time_integrator = SYMPLECTIC;

	// Init boxes	
	auto box0 = addBox(js["Rz"], js["p0"], js["dimension"], js["scale"], js["mass"], boxShape, js["isReduced"], 0);	// No parent
	auto box1 = addBox(js["Rz"], js["p1"], js["dimension"], js["scale"], js["mass"], boxShape, js["isReduced"], 1, box0);
	auto box2 = addBox(js["Rz"], js["p2"], js["dimension"], js["scale"], js["mass"], boxShape, js["isReduced"], 2, box1);

	// Init joints
	auto joint1 = addJoint(box0, box1, js["E_P_J_1"], js["min_theta_1"], js["max_theta_1"]);
	auto joint2 = addJoint(box1, box2, js["E_P_J_2"], js["min_theta_2"], js["max_theta_2"]);

	if (time_integrator == SYMPLECTIC) {
		symplectic_solver = make_shared<SymplecticIntegrator>(boxes, joints, js["isReduced"], js["isMuscle"], js["num_samples_on_muscle"], js["grav"], js["epsilon"]);
	}

}

shared_ptr<Joint> Scene::addJoint(shared_ptr<Rigid> parent, shared_ptr<Rigid> child, json jE_P_J, json jmin_theta, json jmax_theta) {
	Matrix4d E_P_J;
	Eigen::from_json(jE_P_J, E_P_J);
	Matrix4d E_C_J = child->getE().inverse() * parent->getE() * E_P_J;
	auto joint = make_shared<Joint>(E_P_J, E_C_J, jmin_theta.get<double>() / 180.0 * M_PI, jmax_theta.get<double>() / 180.0 * M_PI);
	child->setJoint(joint);
	joint->setChild(child);
	joint->setParent(parent);
	joints.push_back(joint);
	return joint;
}

shared_ptr<Rigid> Scene::addBox(json _R, json _p, json _dimension, json _scale, json _mass, const shared_ptr<Shape> shape, json _isReduced, int id, shared_ptr<Rigid> parent) {
	Matrix3d R;	
	Vector3d p, dimension;
	Eigen::from_json(_R, R);
	Eigen::from_json(_p, p);
	Eigen::from_json(_dimension, dimension);
	double scale = _scale;
	double mass = _mass;
	auto box = make_shared<Rigid>(shape, R, p, dimension, scale, mass, _isReduced, js["grav"]);
	box->setIndex(id);
	box->setParent(parent);
	boxes.push_back(box);
	return box;
}

void Scene::init()
{
	boxShape->init();

}

void Scene::tare()
{
	for (int i = 0; i < (int)boxes.size(); ++i) {
		boxes[i]->tare();
	}
}

void Scene::reset()
{
	t = 0.0;
	for (int i = 0; i < (int)boxes.size(); ++i) {
		boxes[i]->reset();
	}
}

void Scene::step()
{	

	if (time_integrator == SYMPLECTIC) {
		symplectic_solver->step(h);
		for (int i = 0; i < (int)boxes.size(); ++i) {
			boxes[i]->step(h);			
			if (i != 0) {
				//this->phi.segment<6>(6 * (i - 1)) = boxes[i]->getTwist();
			}
		}
	}

	t += h;
}


void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> prog2, shared_ptr<MatrixStack> P) const
{
	for (int i = 0; i < (int)boxes.size(); ++i) {
		boxes[i]->draw(MV, prog, prog2, P);
	}

	symplectic_solver->draw(MV, prog2, P);
}
