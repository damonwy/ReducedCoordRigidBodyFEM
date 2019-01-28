#include "rmpch.h"
#include "ConstraintAttachSoftBody.h"

#include "SoftBody.h"
#include "Node.h"
#include "Body.h"
#include "Vector.h"

using namespace std;
using namespace Eigen;

ConstraintAttachSoftBody::ConstraintAttachSoftBody() {


}

ConstraintAttachSoftBody::ConstraintAttachSoftBody(shared_ptr<SoftBody> softbody) :
	m_softbody(softbody), 
	n_attachments (softbody->m_attach_bodies.size()), 
	n_sliding_nodes(softbody->m_sliding_nodes.size()),
	Constraint(3 * softbody->m_attach_bodies.size()+ softbody->m_sliding_nodes.size(), 0, 0, 0)
{
	
}


void ConstraintAttachSoftBody::computeJacEqM_(MatrixXd &Gm, MatrixXd &Gmdot, VectorXd &gm, VectorXd &gmdot, VectorXd &gmddot) {

	int rowi = idxEM;
	int colSi, colBi;
	Matrix4d E;
	Matrix3d R, W;
	Matrix3x6d G;

	/*for (int i = 0; i < n_attachments; i++) {

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
	}*/

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

	for (int i = 0; i < n_sliding_nodes; i++) {

		colSi = m_softbody->m_sliding_nodes[i]->idxM;
		auto body = m_softbody->m_sliding_bodies[i];

		if (body == nullptr) {
			E = Matrix4d::Identity();
		}
		else {
			E = body->E_wi;
			colBi = body->idxM;
		}

		R = E.block<3, 3>(0, 0);
		Vector3d xi = m_softbody->m_sliding_nodes[i]->x - body->E_iw.block<3, 1>(0, 3);
		//G = SE3::gamma(xi);
		G = SE3::gamma(m_softbody->m_r_sliding[i]);
		auto normal = m_softbody->m_normals_sliding[i];
		Vector3d nor = normal->dir;

		// No velocity in the normal direction
		if (body != nullptr) {
			W = SE3::bracket3(body->phi.segment<3>(0));
			Gm.block<1, 6>(rowi, colBi) = nor.transpose() * R * G;
			//cout << nor.transpose() * R * G;
			Gmdot.block<1, 6>(rowi, colBi) = nor.transpose() * R * W * G;

		}

		Gm.block<1, 3>(rowi, colSi) = -nor;
		//Gm.block<3, 3>(rowi, colSi) = -Matrix3d::Identity();

		Vector4d tem0;
		tem0.segment<3>(0) = m_softbody->m_r_sliding[i];
		tem0(3) = 1.0;
		Vector4d tem1;
		tem1.segment<3>(0) = m_softbody->m_sliding_nodes[i]->x;
		tem1(3) = 1.0;

		Vector4d gmi = E * tem0 - tem1;
		gm.segment<1>(rowi) = nor.transpose() * gmi.segment<3>(0);

		rowi += 1;

	}

}


void ConstraintAttachSoftBody::computeJacEqMSparse_(vector<T> &Gm, vector<T> &Gmdot, VectorXd &gm, VectorXd &gmdot, VectorXd &gmddot) {
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

			Matrix3x6d Gm_temp = R * G;
			Matrix3x6d Gmdot_temp = R * W * G;

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 6; j++) {
					Gm.push_back(T(rowi + i, colBi + j, Gm_temp(i, j)));
					Gmdot.push_back(T(rowi + i, colBi + j, Gmdot_temp(i, j)));
				}
			}
		}

		for (int i = 0; i < 3; i++) {
			Gm.push_back(T(rowi + i, colSi + i, -1.0));		
		}

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

	for (int i = 0; i < n_sliding_nodes; i++) {

		colSi = m_softbody->m_sliding_nodes[i]->idxM;
		auto body = m_softbody->m_sliding_bodies[i];

		if (body == nullptr) {
			E = Matrix4d::Identity();
		}
		else {
			E = body->E_wi;
			colBi = body->idxM;
		}

		R = E.block<3, 3>(0, 0);
		Vector3d xi = m_softbody->m_sliding_nodes[i]->x - body->E_iw.block<3, 1>(0, 3);
		G = SE3::gamma(m_softbody->m_r_sliding[i]);
		auto normal = m_softbody->m_normals_sliding[i];
		Vector3d nor = normal->dir;

		// No velocity in the normal direction
		if (body != nullptr) {
			W = SE3::bracket3(body->phi.segment<3>(0));
			Vector6d temp = nor.transpose() * R * G;
			Vector6d dottemp = nor.transpose() * R * W * G;
			for (int i = 0; i < 6; ++i) {
				Gm.push_back(T(rowi, colBi + i, temp(i)));
				Gmdot.push_back(T(rowi, colBi + i, dottemp(i)));
			}
		}

		for (int i = 0; i < 3; i++) {
			Gm.push_back(T(rowi, colSi + i, -nor(i)));
		}

		Vector4d tem0;
		tem0.segment<3>(0) = m_softbody->m_r_sliding[i];
		tem0(3) = 1.0;
		Vector4d tem1;
		tem1.segment<3>(0) = m_softbody->m_sliding_nodes[i]->x;
		tem1(3) = 1.0;

		Vector4d gmi = E * tem0 - tem1;
		gm.segment<1>(rowi) = nor.transpose() * gmi.segment<3>(0);

		rowi += 1;
	}
}
