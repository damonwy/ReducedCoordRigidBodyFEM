#include "World.h"


#include <iostream>
#include <fstream>
#include <json.hpp>
#include "Joint.h"
#include "JointRevolute.h"
#include "Body.h"
#include "MatrixStack.h"
#include "Program.h"
#include "SE3.h"
#include "JsonEigen.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

World::World():
nr(0),
nm(0)
{

}

World::World(WorldType type):
m_type(type),
nr(0),
nm(0)
{
}

World::~World() {
}

void World::load(const std::string &RESOURCE_DIR) {

	//read a JSON file
	ifstream i(RESOURCE_DIR + "input.json");
	json js;
	i >> js;
	i.close();

	double density;
	Eigen::Vector3d sides;
	Matrix4d E;
	Vector3d p;


	switch (m_type)
	{
	case SERIAL_CHAIN: 
		{	
			m_h = 1.0e-3;
			density = 1.0;
			m_grav << 0.0, -9.81, 0.0;
			Eigen::from_json(js["sides"], sides);
			m_nbodies = 5;
			m_njoints = 5;
			m_Hexpected = 10000; // todo
			m_tspan << 0.0, 5.0;
			
			// Inits rigid bodies
			for (int i = 0; i < m_nbodies; i++) {

				auto body = make_shared<Body>(density, sides);

				// Inits joints
				if (i == 0) {
					auto joint = make_shared<JointRevolute>(body, Vector3d::UnitZ());
					joint->setJointTransform(Matrix4d::Identity());
					m_joints.push_back(joint);
				}
				else {
					auto joint = make_shared<JointRevolute>(body, Vector3d::UnitZ(), m_joints[i-1]);
					p << 10.0, 0.0, 0.0; // todo
					E = SE3::RpToE(Matrix3d::Identity(), p);

					joint->setJointTransform(E);
					m_joints.push_back(joint);
				}

				p << 5.0, 0.0, 0.0;
				E = SE3::RpToE(Matrix3d::Identity(), p);

				body->setTransform(E);
				body->load(RESOURCE_DIR);

				m_joints[i]->m_q(0) = 0.0;

				m_bodies.push_back(body);
			}
			break;
		}
	case DIFF_REVOLUTE_AXES:
		break;
	case BRANCHING:
		break;
	case SHPERICAL_JOINT:
		break;
	case LOOP:
		break;
	case JOINT_TORQUE:
		break;
	case JOINT_LIMITS:
		break;
	case EQUALITY_CONSTRAINED_ANGLES:
		break;
	case EQUALITY_AND_LOOP:
		break;
	case HYBRID_DYNAMICS:
		break;
	case EXTERNAL_WORLD_FORCE:
		break;
	case JOINT_STIFFNESS:
		break;
	case SPRINGS:
		break;
	default:
		break;
	}
	
}

void World::init() {
	for (int i = 0; i < m_nbodies; i++) {
		m_bodies[i]->init(nm);
		if (i < m_nbodies-1) {
			m_bodies[i]->next = m_bodies[i + 1];
		}
	}

	for (int i = 0; i < m_njoints; i++) {
		m_joints[i]->init(nr);
	}

	//joint ordering
	// todo

	for (int i = 0; i < m_njoints; i++) {
		if (i < m_njoints - 1) {
			m_joints[i]->next = m_joints[i + 1];
		}
		if (i > 0) {
			m_joints[i]->prev = m_joints[i - 1];
		}
	}

	m_joints[0]->update();


	// init constraints

	// init collisions

}

void World::update() {
	for (int i = 0; i < m_nbodies; i++) {
		m_bodies[i]->update();
	}
}

void World::updateQ() {
	for (int i = 0; i < m_nbodies; i++) {
		auto body = m_bodies[i];
		//body->updateQ();

	}
}

void World::updateQDot() {

}

int World::getNsteps() {
	// Computes the number of results
	int nsteps = (m_tspan(1) - m_tspan(0)) / m_h;
	return nsteps;
}

void World::addBody(shared_ptr<Body> body) {

}

void World::addJoint(shared_ptr<Joint> joint) {

}

std::shared_ptr<Body> World::getBody(int uid) {
	MapBodyUID::const_iterator it = m_bodyUID.find(uid);
	return (it == m_bodyUID.end() ? NULL : it->second);
}

std::shared_ptr<Body> World::getBody(const std::string &name) {
	MapBodyName::const_iterator it = m_bodyName.find(name);
	return (it == m_bodyName.end() ? NULL : it->second);
}

std::shared_ptr<Joint> World::getJoint(int uid) {
	MapJointUID::const_iterator it = m_jointUID.find(uid);
	return (it == m_jointUID.end() ? NULL : it->second);
}

std::shared_ptr<Joint> World::getJoint(const std::string &name) {
	MapJointName::const_iterator it = m_jointName.find(name);
	return (it == m_jointName.end() ? NULL : it->second);
}

void World::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) {
	// Draw rigid bodies
	for (int i = 0; i < m_nbodies; i++) {
		m_bodies[i]->draw(MV, prog, P);
	}

	// Draw joints
	for (int i = 0; i < m_njoints; i++) {
		m_joints[i]->draw(MV, progSimple, P);
	}

	// Draw constraints
	for (int i = 0; i < 1; i++) {

	}

	// Draw collisions


}