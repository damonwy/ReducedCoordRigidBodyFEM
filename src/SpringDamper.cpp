#include "SpringDamper.h"

#include <iostream>
#include "Body.h"
#include "Node.h"
#include "MatrixStack.h"
#include "Program.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

SpringDamper::SpringDamper() {



}


SpringDamper::SpringDamper(int &countS, int &countCM):
Spring(countS, countCM)
{

	m_K = 1.0;
	m_L = 0.0;	// 
	m_d = 1.0;	

}

void SpringDamper::setAttachments(shared_ptr<Body> body0, Vector3d r0, shared_ptr<Body> body1, Vector3d r1) {
	// Attaches this spring to body0 and body1
	m_body0 = body0;
	m_body1 = body1;
	m_r0 = r0;
	m_r1 = r1;
}

void SpringDamper::init() {

	for (int i = 0; i < (int)m_nodes.size(); i++) {
		m_nodes[i]->init();
	}



}

void SpringDamper::load(const std::string &RESOURCE_DIR) {

	for (int i = 0; i < (int)m_nodes.size(); i++) {
		m_nodes[i]->load(RESOURCE_DIR);
	}

}

void SpringDamper::countDofs_(int &nm, int &nr) {




}

Eigen::VectorXd SpringDamper::gatherDofs_(Eigen::VectorXd y, int nr) {



}

Eigen::VectorXd SpringDamper::gatherDDofs_(Eigen::VectorXd ydot, int nr) {



}

void SpringDamper::scatterDofs_(Eigen::VectorXd &y, int nr) {



}

void SpringDamper::scatterDDofs_(Eigen::VectorXd &ydot, int nr) {



}

Eigen::MatrixXd SpringDamper::computeMass_(Eigen::Vector3d grav, Eigen::MatrixXd M) {



}

Eigen::VectorXd SpringDamper::computeForce_(Eigen::Vector3d grav, Eigen::VectorXd f) {



}

Energy SpringDamper::computeEnergies_(Eigen::Vector3d grav, Energy ener) {
	Matrix4d E0, E1;
	E0.setIdentity();
	E1.setIdentity();

	if (m_body0 != nullptr) {
		E0 = m_body0->E_wi;
	}

	if (m_body1 != nullptr) {
		E1 = m_body1->E_wi;
	}

	Vector4d temp0, temp1;
	temp0 << m_r0, 1.0;
	temp1 << m_r1, 1.0;
	Vector3d x0_w = E0.block<3, 4>(0, 0)*temp0;
	Vector3d x1_w = E1.block<3, 4>(0, 0)*temp1;

	m_l = (x1_w - x0_w).norm();
	if (m_L == 0.0) {
		m_L = m_l;
	}

	double e = (m_l - m_L) / m_L;
	ener.V = ener.V + 0.5 * m_K * e * e;
	
	return ener;
}

Eigen::MatrixXd SpringDamper::computeJacobian_(Eigen::MatrixXd J) {




}

void SpringDamper::draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const {
	// Draw nodes
	prog->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	MV->pushMatrix();
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	int n_nodes = (int)m_nodes.size();
	for (int i = 0; i < n_nodes; i++) {
		m_nodes[i]->draw(MV, prog);
	}

	MV->popMatrix();
	prog->unbind();

	// Draw line segments

	progSimple->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	glLineWidth(5);
	glBegin(GL_LINES);
	glColor3f(0.0f, 0.0f, 0.0f);

	for (int i = 0; i < n_nodes - 1; i++) {

		Vector3f x0 = m_nodes[i]->x.cast<float>();
		Vector3f x1 = m_nodes[i + 1]->x.cast<float>();

		glVertex3f(x0(0), x0(1), x0(2));
		glVertex3f(x1(0), x1(1), x1(2));
	}

	glEnd();
	progSimple->unbind();
	
}


