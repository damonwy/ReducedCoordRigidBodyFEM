#include "World.h"


#include <iostream>
#include <fstream>
#include <json.hpp>
#include "Rigid.h"
#include "Joint.h"
#include "Body.h"
#include "MatrixStack.h"
#include "Program.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

World::World() {

}

World::~World() {
}

void World::load(const std::string &RESOURCE_DIR) {
	m_nlinks = 2;
	m_grav << 0.0, -9.81, 0.0;
	double density = 1.0;
	Eigen::Vector3d sides;
	sides << 1.0, 4.0, 1.0;

	for (int i = 0; i < m_nlinks; i++) {

		auto body = make_shared<Body>(density, sides);


		m_bodies.push_back(body);
	}
}

void World::init() {
	for (int i = 0; i < (int)m_bodies.size(); i++) {
		m_bodies[i]->init();
	}

	// init constraints

	// init collisions

}

void World::update() {
	for (int i = 0; i < (int)m_bodies.size(); i++) {
		m_bodies[i]->update();
	}
}

void World::updateQ() {
	for (int i = 0; i < (int)m_bodies.size(); i++) {
		auto body = m_bodies[i];
		body->updateQ();

	}
}

void World::updateQDot() {

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


}

void World::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, shared_ptr<MatrixStack> P) {
	// Draw rigid bodies
	for (int i = 0; i < (int)m_bodies.size(); i++) {
		m_bodies[i]->draw(MV, prog, P);
	}

	// Draw constraints
	for (int i = 0; i < 100; i++) {

	}

	// Draw collisions


}