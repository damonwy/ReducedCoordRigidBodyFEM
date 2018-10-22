#include "CompSphere.h"

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <iostream>
#include <fstream>

#include "Body.h"
#include "Shape.h"
#include "SE3.h"
#include "MatrixStack.h"
#include "Program.h"
#include "Node.h"

#include <json.hpp>

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

CompSphere::CompSphere() {
	m_O = make_shared<Node>();
	m_O->x0.setZero();

}

CompSphere::CompSphere(std::shared_ptr<Body> parent, double r) :m_parent(parent), m_r(r){
	m_O = make_shared<Node>();
	m_O->x0.setZero();

}

CompSphere::~CompSphere() {
}

void CompSphere::init() {
	m_shape->init();
	m_O->init();
}

void CompSphere::load(const std::string &RESOURCE_DIR) {
	m_shape = make_shared<Shape>();
	m_shape->loadMesh(RESOURCE_DIR + "sphere2.obj");
	m_O->load(RESOURCE_DIR);
}


void CompSphere::update() {
	E_wi = m_parent->E_wi * E_ji;
	m_O->update(E_wi);

	if (next != nullptr) {
		this->next->update();
	}
}

void CompSphere::setTransform(Eigen::Matrix4d E) {

	E_ji = E;

}

void CompSphere::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P)const {

	prog->bind();
	if (m_shape) {
		glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
		glUniform3f(prog->getUniform("lightPos1"), 66.0, 25.0, 25.0);
		glUniform1f(prog->getUniform("intensity_1"), 0.6);
		glUniform3f(prog->getUniform("lightPos2"), -66.0, 25.0, 25.0);
		glUniform1f(prog->getUniform("intensity_2"), 0.2);
		glUniform1f(prog->getUniform("s"), 300);
		glUniform3f(prog->getUniform("ka"), 0.2, 0.2, 0.2);
		glUniform3f(prog->getUniform("kd"), 0.8, 0.7, 0.7);
		glUniform3f(prog->getUniform("ks"), 1.0, 0.9, 0.8);
		m_O->draw(MV, prog);
		
		MV->pushMatrix();
		MV->multMatrix(eigen_to_glm(E_wi));
		MV->scale(m_r);
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		m_shape->draw(prog);
		MV->popMatrix();


	}
	prog->unbind();

}

