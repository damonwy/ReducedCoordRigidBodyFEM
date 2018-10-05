#include "ConstraintAttachSoftBody.h"

#include <iostream>
#include <fstream>
#include <json.hpp>
#include "SoftBody.h"
#include "Node.h"
#include "Body.h"
#include "SE3.h"

using namespace std;
using namespace Eigen;

ConstraintAttachSoftBody::ConstraintAttachSoftBody() {


}

ConstraintAttachSoftBody::ConstraintAttachSoftBody(shared_ptr<SoftBody> softbody) :
	m_softbody(softbody), 
	n_attachments (softbody->m_attach_bodies.size()), 
	Constraint(3 * softbody->m_attach_bodies.size(), 0, 0, 0)
{
	
	
}


void ConstraintAttachSoftBody::computeJacEqM_(MatrixXd &Gm, MatrixXd &Gmdot, VectorXd &gm, VectorXd &gmdot, VectorXd &gmddot) {

	int rowi = idxEM;
	int colSi, colBi;
	Matrix4d E;
	Matrix3d R, W;
	Matrix3x6d G;

	for (int i = 0; i < n_attachments; i++) {

		colSi = m_softbody->m_attach_nodes[i]->idxM;
		auto body = m_softbody->m_attach_bodies[i];

		if (body == nullptr) {
			E = Matrix4d::Identity();
		}
		else {
			E = body->E_wi;
			colBi = body->idxM;
		}

		R = E.block<3, 3>(0, 0);
		G = SE3::gamma(m_softbody->m_r[i]);

		if (body != nullptr) {
			W = SE3::bracket3(body->phi.segment<3>(0));
			Gm.block<3, 6>(rowi, colBi) = R * G;
			Gmdot.block<3, 6>(rowi, colBi) = R * W * G;

		}

		Gm.block<3, 3>(rowi, colSi) = -Matrix3d::Identity();

		Vector4d tem0;
		tem0.segment<3>(0) = m_softbody->m_r[i];
		tem0(3) = 1.0;
		Vector4d tem1;
		tem1.segment<3>(0) = m_softbody->m_attach_nodes[i]->x;
		tem1(3) = 1.0;

		Vector4d gmi = E * tem0 - tem1;
		gm.segment<3>(rowi) = gmi.segment<3>(0);
		rowi += 3;
	}
}
