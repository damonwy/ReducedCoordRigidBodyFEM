#include "World.h"


#include <iostream>
#include <fstream>
#include <json.hpp>
#include "Rigid.h"
#include "Joint.h"
#include "Body.h"

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
}


void World::addBody(shared_ptr<Body> body) {

}

void World::addJoint(shared_ptr<Joint> joint) {

}

std::shared_ptr<Body> World::getBody(int uid) {

}

std::shared_ptr<Joint> World::getJoint(int uid) {


}