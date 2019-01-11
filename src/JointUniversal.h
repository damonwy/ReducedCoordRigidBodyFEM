#pragma once
#include "Joint.h"
#ifndef REDUCEDCOORD_SRC_JOINTUNIVERSAL_H_
#define REDUCEDCOORD_SRC_JOINTUNIVERSAL_H_

#define EIGEN_USE_MKL_ALL


#include "Body.h"
#include "SE3.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"
#include "ConstraintPrescJoint.h"


class JointUniversal : public Joint {

public:
	JointUniversal() {}
	JointUniversal(std::shared_ptr<Body> body, std::shared_ptr<Joint> parent = nullptr):
	Joint(body, 2, parent)
	{


	}

	void load(const std::string &RESOURCE_DIR, std::string joint_shape) {


	}

	virtual ~JointUniversal() {}
	void update_() {
		double q0 = m_q(0);
		double q1 = m_q(1);
		double dq0 = m_qdot(0);
		double dq1 = m_qdot(1);

		double c0 = cos(q0);
		double c1 = cos(q1);
		double s0 = sin(q0);
		double s1 = sin(q1);

		m_Q = Matrix4d::Identity();
		Vector3d temp;
		temp << c1, s0*s1, -c0 * s1;
		m_Q.block<3, 1>(0, 0) = temp;
		temp << 0, c0, s0;
		m_Q.block<3, 1>(0, 1) = temp;
		temp << s1, -s0 * c1, c0 * c1;
		m_Q.block<3, 1>(0, 2) = temp;

		m_S(0, 0) = c1;
		m_S(2, 0) = s1;
		m_S(1, 1) = 1.0;
		m_Sdot(0, 0) = -s1 * dq1;
		m_Sdot(2, 0) = c1 * dq1;
	}

protected:

	void draw_(std::shared_ptr<MatrixStack> MV,
		const std::shared_ptr<Program> prog,
		const std::shared_ptr<Program> prog2,
		std::shared_ptr<MatrixStack> P) const {

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
			MV->multMatrix(eigen_to_glm(E_wj));
			double alpha = 1.0;
			if (presc->activeER) {
				alpha = 2.0;
			}

			MV->scale(alpha * m_draw_radius);
			glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
			m_jointShape->draw(prog);
			MV->popMatrix();
		}

		prog->unbind();
	}
};

#endif // REDUCEDCOORD_SRC_JOINTUNIVERAL_H_
