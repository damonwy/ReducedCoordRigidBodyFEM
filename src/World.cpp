#include "World.h"


#include <iostream>
#include <fstream>
#include <json.hpp>
#include "Rigid.h"
#include "Joint.h"
#include "Body.h"

World::World() {

}

World::~World() {



}

void World::load(const std::string &RESOURCE_DIR) {
	m_nlinks = 2;
	m_grav << 0.0, -9.81, 0.0;

	for (int i = 0; i < m_nlinks; i++) {

		Body *body = new Body();

		// Set dimensions



	}



}

void World::addBody(Rigid *body) {

}

void World::addJoint(Joint *joint) {

}