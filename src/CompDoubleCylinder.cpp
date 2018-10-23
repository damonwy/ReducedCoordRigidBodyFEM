#include "CompDoubleCylinder.h"

#include "Body.h"
#include "Shape.h"
#include "SE3.h"
#include "MatrixStack.h"
#include "Program.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

CompDoubleCylinder::CompDoubleCylinder(): Comp() {
}

CompDoubleCylinder::CompDoubleCylinder(shared_ptr<Body> parentA, double rA, shared_ptr<Body> parentB, double rB) : 
Comp(), m_rA(rA), m_rB(rB), m_parentA(parentA), m_parentB(parentB)
{

}

CompDoubleCylinder::~CompDoubleCylinder() {
}

void CompDoubleCylinder::load(const std::string &RESOURCE_DIR, std::string shapeA, std::string shapeB) {
	m_shapeA = make_shared<Shape>();
	m_shapeA->loadMesh(RESOURCE_DIR + shapeA);

	m_shapeB = make_shared<Shape>();
	m_shapeB->loadMesh(RESOURCE_DIR + shapeB);
}

void CompDoubleCylinder::init() {
	m_shapeA->init();
	m_shapeB->init();
}

void CompDoubleCylinder::update() {
	E_wiA = m_parentA->E_wi * E_jiA;
	E_wiB = m_parentB->E_wi * E_jiB;
	
	m_OA->update(E_wiA);
	m_OB->update(E_wiB);
	m_ZA->update(E_wiA);
	m_ZB->update(E_wiB);

	if (next != nullptr) {
		this->next->update();
	}
}

void CompDoubleCylinder::setTransformA(Matrix4d E) {
	E_jiA = E;
}

void CompDoubleCylinder::setTransformB(Matrix4d E) {

	E_jiB = E;
}

void CompDoubleCylinder::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P)const {
	prog->bind();
	if (m_shapeA && m_shapeB) {
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
		MV->multMatrix(eigen_to_glm(E_wiA));

		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		m_shapeA->draw(prog);
		MV->popMatrix();

		MV->pushMatrix();
		MV->multMatrix(eigen_to_glm(E_wiB));
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		m_shapeB->draw(prog);
		MV->popMatrix();

	}
	prog->unbind();
}
