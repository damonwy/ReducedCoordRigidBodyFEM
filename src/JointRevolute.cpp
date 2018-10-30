#include "JointRevolute.h"

#include <iostream>

#include "Body.h"
#include "SE3.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"

using namespace std;
using namespace Eigen;

JointRevolute::JointRevolute() {

}


JointRevolute::JointRevolute(std::shared_ptr<Body> body, Eigen::Vector3d axis, std::shared_ptr<Joint> parent):
Joint(body, 1, parent)
{
	m_axis = axis;
}

void JointRevolute::load(const std::string &RESOURCE_DIR, std::string joint_shape) {

	m_jointShape = make_shared<Shape>();
	m_jointShape->loadMesh(RESOURCE_DIR + "sphere2.obj");

}


JointRevolute::~JointRevolute() {

}

void JointRevolute::updateSelf() {
	Matrix3d R = SE3::aaToMat(m_axis, m_q(0));
	Matrix4d Q;
	Q.setIdentity();
	Q.block<3, 3>(0, 0) = R;
	m_Q = Q;
	//E_pj = E_pj0 * Q;
	
	m_S.block<3, 1>(0, 0) = m_axis;
}

void JointRevolute::drawSelf(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
	prog->bind();

	double r = 0.5;
	if (m_jointShape) {
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
		MV->multMatrix(eigen_to_glm(E_wj));
		MV->scale(r);
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		m_jointShape->draw(prog);
		MV->popMatrix();
	}
	prog->unbind();
}

