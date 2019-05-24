#include "rmpch.h"
#include "JointSplineCurve.h"
#include "Body.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

const Matrix4d JointSplineCurve::m_B = getB();
const Matrix4x3d JointSplineCurve::m_B1 = getB1();
const Matrix4x2d JointSplineCurve::m_B2 = getB2();

JointSplineCurve::JointSplineCurve()
{
}

JointSplineCurve::JointSplineCurve(shared_ptr<Body> body, shared_ptr<Joint> parent):
Joint(body, 1, parent)
{
}

void JointSplineCurve::load(const string &RESOURCE_DIR, string joint_shape) {

	m_jointShape = make_shared<Shape>();
	m_jointShape->loadMesh(RESOURCE_DIR + joint_shape);
	m_jointSphereShape = make_shared<Shape>();
	m_jointSphereShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
}

void JointSplineCurve::init(int &nm, int &nr) {
	if (m_jointShape) {
		m_jointShape->init();
	}
	if (m_jointSphereShape) {
		m_jointSphereShape->init();
	}

	m_body->setJoint(getJoint());

	if (m_parent != nullptr) {
		m_parent->addChild(getJoint());
	}
	countDofs(nm, nr);
}

void JointSplineCurve::addControlFrame(Matrix4d C) {
	m_Cs.push_back(C);
	int ncfs = m_Cs.size();
	if (ncfs >= 2) {
		Matrix4d C0 = m_Cs[ncfs - 2];
		Matrix4d C1 = C;
		Matrix4d x = C0.colPivHouseholderQr().solve(C1);
		m_dCs.push_back(SE3::log(x));
	}

	// Cyclic
	Matrix4d C0 = C;
	Matrix4d C1 = m_Cs[0];
	Matrix4d x = C0.colPivHouseholderQr().solve(C1);

	if (m_dCs.size() < 1) {
		m_dCs.push_back(SE3::log(x));
	}
	else {
		m_dCs[0] = SE3::log(x);
	}
}

void JointSplineCurve::update_() {
	m_Q = evalQ(m_q(0));
	Vector6d S, dSdq;
	evalS(m_q(0), S, dSdq);
	m_S = S;
	m_Sdot = dSdq * m_qdot(0);
}

double JointSplineCurve::Bsum(int i, double q) {
	// Evaluates Btilde
	Vector4d qvec;
	Vector4d b;
	qvec << 1, q, q * q, q * q * q;
	if (i == 0) {
		b = m_B.row(0) + m_B.row(1) + m_B.row(2) + m_B.row(3);
	}
	else if (i == 1) {
		b = m_B.row(1) + m_B.row(2) + m_B.row(3);
	}
	else if (i == 2) {
		b = m_B.row(2) + m_B.row(3);
	}
	else {
		b = m_B.row(3);
	}
	double f = b.transpose() * qvec;
	return f;
}

double JointSplineCurve::dBsum(int i, double q) {
	// Evaluates dBtilde / dq
	Vector3d qvec;
	Vector3d b;
	qvec << 1, 2 * q, 3 * q * q;
	if (i == 0) {
		b = m_B1.row(0) + m_B1.row(1) + m_B1.row(2) + m_B1.row(3);
	}
	else if (i == 1) {
		b = m_B1.row(1) + m_B1.row(2) + m_B1.row(3);
	}
	else if (i == 2) {
		b = m_B1.row(2) + m_B1.row(3);
	}
	else {
		b = m_B1.row(3);
	}
	double df = b.transpose() * qvec;
	return df;
}

double JointSplineCurve::d2Bsum(int i, double q) {
	// Evaluates d ^ 2Btilde / dq ^ 2
	Vector2d qvec;
	qvec << 2, 6 * q;
	Vector2d b;
	if (i == 0) {
		b = m_B2.row(0) + m_B2.row(1) + m_B2.row(2) + m_B2.row(3);
	}
	else if (i == 1) {
		b = m_B2.row(1) + m_B2.row(2) + m_B2.row(3);
	}
	else if (i == 2) {
		b = m_B2.row(2) + m_B2.row(3);
	}
	else {
		b = m_B2.row(3);
	}

	double d2f = b.transpose() * qvec;
	return d2f;
}

Matrix4d JointSplineCurve::evalQ(double q) const {
	Matrix4d Q;

	// Evaluate spline frame
	int ncfs = m_Cs.size();
	// Wrap around
	int qmax = ncfs;
	if (q < 0) {
		q += qmax;
	}
	else if (q >= qmax) {
		q -= qmax;
	}

	int k = int(floor(q)); // starting control frame (0-index)
	if (k >= ncfs) {
		k -= 1; // overflow
	}

	double q_ = q - k; // local q in [0, 1]
	
	Q = m_Cs[k];
	for (int i = 1; i < 4; ++i) {
		int ki = k + i;
		if (ki > ncfs - 1) {
			ki = ki - ncfs; // Wrap for cyclic
		}
		double Bsum = JointSplineCurve::Bsum(i, q_);
		Matrix4d dC = SE3::bracket6(m_dCs[ki]);

		Matrix4d dCBsum = dC * Bsum;
		Q = Q * SE3::exp(dCBsum);
	}

	return Q;
}

void JointSplineCurve::evalS(double q, Vector6d &S, Vector6d &dSdq) {

	// Evaluates spline frame derivatives
	int ncfs = m_Cs.size();
	// Wrap around
	int qmax = ncfs;
	if (q < 0) {
		q += qmax;
	}
	else if (q >= qmax) {
		q -= qmax;
	}

	int k = int(floor(q)); // starting control frame 
	if (k >= ncfs) {
		k -= 1; // overflow
	}

	double q_ = q - k; // local q in [0, 1]
	for (int i = 1; i < 4; ++i) {
		int ki = k + i;
		if (ki > ncfs - 1) {
			ki = ki - ncfs; // Wrap for cyclic
		}
		Vector6d dC = m_dCs[ki];
		double dBsum = JointSplineCurve::dBsum(i, q_);
		double d2Bsum = JointSplineCurve::d2Bsum(i, q_);
		if (i == 1) {
			S = dC *dBsum;
			dSdq = dC * d2Bsum;
		}
		else {
			double Bsum = JointSplineCurve::Bsum(i, q_);
			Vector6d dCBsum = dC * Bsum;
			Matrix6d Ad = SE3::adjoint(SE3::inverse(SE3::exp(dCBsum)));
			Matrix6d ad = SE3::ad(S);
			Vector6d AdS = Ad * S;
			S = dC * dBsum + AdS;
			Vector6d addCdBsum = ad * dC * dBsum;
			Vector6d AdSum = Ad *(dSdq + addCdBsum);
			dSdq = dC * d2Bsum + AdSum;
		}
	}
}

void JointSplineCurve::draw_(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
	int ncfs = m_Cs.size();
	Matrix4d E_wp;

	if (getParent() == nullptr) {
		E_wp.setIdentity();
	}
	else {
		E_wp = getParent()->E_wj;
	}

	Matrix4d E_wj0 = E_wp * E_pj0;

	progSimple->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	// Draw control frames
	for (int i = 0; i < ncfs; ++i) {
		MV->pushMatrix();
		Matrix4d E = E_wj0 * m_Cs[i];
		MV->multMatrix(eigen_to_glm(E));
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		glLineWidth(5);
		glBegin(GL_LINES);
		// X axis
		glColor3f(1.0, 0.0, 0.0);
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(2.0, 0.0, 0.0);

		// Y axis
		glColor3f(0.0, 1.0, 0.0);
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, 2.0, 0.0);

		// Z axis
		glColor3f(0.0, 0.0, 1.0);
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, 0.0, 2.0);

		glEnd();
		MV->popMatrix();
	}

	int qmax = ncfs;
	VectorXd q = VectorXd::LinSpaced(qmax*8, 0, qmax);
	
	// Draw spline frames
	for (int i = 0; i < q.size(); ++i) {
		double qi = q(i);
		Matrix4d Q = evalQ(qi);
		MV->pushMatrix();
		Matrix4d E = E_wj0 * Q;
		MV->multMatrix(eigen_to_glm(E));
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		glLineWidth(2);
		glBegin(GL_LINES);
		// X axis
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(1.0f, 0.0f, 0.0f);

		// Y axis
		glColor3f(0.0f, 1.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 1.0f, 0.0f);

		// Z axis
		glColor3f(0.0f, 0.0f, 1.0f);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 1.0f);

		glEnd();
		MV->popMatrix();
	}
	
	// Draw current frame
	MV->pushMatrix();
	Matrix4d E = E_wj;
	MV->multMatrix(eigen_to_glm(E));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	glLineWidth(4);
	glBegin(GL_LINES);
	// X axis
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(3.0f, 0.0f, 0.0f);

	// Y axis
	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 3.0f, 0.0f);

	// Z axis
	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 3.0f);

	glEnd();
	MV->popMatrix();
	progSimple->unbind();

	float r = 0.5f;
	prog->bind();
	if (m_jointShape) {
		glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
		glUniform3f(prog->getUniform("lightPos1"), 66.0f, 25.0f, 25.0f);
		glUniform1f(prog->getUniform("intensity_1"), 0.6f);
		glUniform3f(prog->getUniform("lightPos2"), -66.0f, 25.0f, 25.0f);
		glUniform1f(prog->getUniform("intensity_2"), 0.2f);
		glUniform1f(prog->getUniform("s"), 300.0f);
		glUniform3f(prog->getUniform("ka"), 0.2f, 0.2f, 0.2f);
		glUniform3f(prog->getUniform("kd"), 0.8f, 0.7f, 0.7f);
		glUniform3f(prog->getUniform("ks"), 1.0f, 0.9f, 0.8f);

		MV->pushMatrix();
		Matrix4d E_wj0 = E_wp * E_pj0;
		MV->multMatrix(eigen_to_glm(E_wj0));

		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		m_jointShape->draw(prog);
		MV->popMatrix();

	}
	if (m_jointSphereShape) {
		MV->pushMatrix();
		MV->multMatrix(eigen_to_glm(E_wj));
		MV->scale(r);
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		m_jointSphereShape->draw(prog);
		MV->popMatrix();
	}
	prog->unbind();
}
