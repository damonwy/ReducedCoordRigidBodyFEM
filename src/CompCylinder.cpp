#include "CompCylinder.h"

#include "Body.h"
#include "Shape.h"
#include "SE3.h"
#include "MatrixStack.h"
#include "Program.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

CompCylinder::CompCylinder() {
}

CompCylinder::CompCylinder(shared_ptr<Body> parent, double r) : m_parent(parent), m_r(r){
}


CompCylinder::~CompCylinder() {

}

void CompCylinder::init() {
	m_shape->init();
}

void CompCylinder::update() {
	E_wi = m_parent->E_wi * E_ji;
	if (next != nullptr) {
		this->next->update();
	}

}

void CompCylinder::setTransform(Matrix4d E) {
	E_ji = E;
	update();
}

void CompCylinder::load(const string &RESOURCE_DIR, string shape) {
	m_shape = make_shared<Shape>();
	m_shape->loadMesh(RESOURCE_DIR + shape);
}

void CompCylinder::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, shared_ptr<MatrixStack> P)const {
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

		MV->pushMatrix();
		MV->multMatrix(eigen_to_glm(E_wi));

		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		m_shape->draw(prog);
		MV->popMatrix();

	}
	prog->unbind();

}
