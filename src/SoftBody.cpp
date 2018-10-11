#define TETLIBRARY
#include <tetgen.h>

#include "SoftBody.h"

#include <iostream>
#include <fstream>
#include <cmath>        // std::abs

#include <json.hpp>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "GLSL.h"
#include "MatrixStack.h"
#include "Program.h"

#include "Node.h"
#include "FaceTriangle.h"
#include "Tetrahedron.h"
#include "Body.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

SoftBody::SoftBody() {
	m_color << 1.0f, 1.0f, 0.0f;
}

SoftBody::SoftBody(double density, double young, double poisson):
m_density(density), m_young(young), m_poisson(poisson)
{
	m_color << 1.0f, 1.0f, 0.0f;
}

void SoftBody::load(const string &RESOURCE_DIR, const string &MESH_NAME) {

	// Tetrahedralize 3D mesh
	tetgenio input_mesh, output_mesh;
	input_mesh.load_ply((char *)(RESOURCE_DIR + MESH_NAME).c_str());
	tetrahedralize("pq1.5z", &input_mesh, &output_mesh);

	double r = 0.1;

	// Create Nodes
	for (int i = 0; i < output_mesh.numberofpoints; i++) {
		auto node = make_shared<Node>();	
		node->r = r;


		node->x0 << output_mesh.pointlist[3 * i + 0],
					output_mesh.pointlist[3 * i + 1],
					output_mesh.pointlist[3 * i + 2];

		node->x = node->x0;
		node->v0.setZero();
		node->v = node->v0;
		node->m = 0.0;
		node->i = i;
		node->load(RESOURCE_DIR);
		/*if (node->x(1) > 0.5) {
			node->fixed = true;
		}
		else {
			node->fixed = false;
		}*/

		m_nodes.push_back(node);
	}

	// Create Faces
	for (int i = 0; i < output_mesh.numberoftrifaces; i++) {
		auto triface = make_shared<FaceTriangle>();
		
		for (int ii = 0; ii < 3; ii++) {
			auto node = m_nodes[output_mesh.trifacelist[3 * i + ii]];
			node->m_nfaces ++;
			triface->m_nodes.push_back(node);
		}
		triface->update();
		m_trifaces.push_back(triface);
	}

	// Create Tets
	vector<shared_ptr<Node>> tet_nodes;
	for (int i = 0; i < output_mesh.numberoftetrahedra; i++) {
		tet_nodes.clear();
		for (int ii = 0; ii < 4; ii++) {
			tet_nodes.push_back(m_nodes[output_mesh.tetrahedronlist[4 * i + ii]]);
		}
		auto tet = make_shared<Tetrahedron>(m_young, m_poisson, m_density, tet_nodes);
		m_tets.push_back(tet);
	}

	// Fix the normal of top and bottom surface
	for (int i = 0; i < m_trifaces.size(); i++) {
		auto triface = m_trifaces[i];

		Vector3d p0 = triface->m_nodes[0]->x;
		Vector3d p1 = triface->m_nodes[1]->x;
		Vector3d p2 = triface->m_nodes[2]->x;

		Vector3d normal = triface->computeNormal();

		if ((normal - m_trifaces[0]->m_normal).norm() < 0.1) {
			triface->isFlat = true;
		}
		if ((normal + m_trifaces[0]->m_normal).norm() < 0.1) {
			triface->isFlat = true;
		}
	}
}

void SoftBody::init() {

	for (int i = 0; i < m_nodes.size(); i++) {
		m_nodes[i]->init();
	}
	// Init Buffers
	posBuf.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();

	posBuf.resize(m_trifaces.size() * 9);
	norBuf.resize(m_trifaces.size() * 9);
	eleBuf.resize(m_trifaces.size() * 3);
	updatePosNor();

	for (int i = 0; i < m_trifaces.size(); i++) {
		eleBuf[3 * i + 0] = 3 * i;
		eleBuf[3 * i + 1] = 3 * i + 1;
		eleBuf[3 * i + 2] = 3 * i + 2;
	}

	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);

	glGenBuffers(1, &norBufID);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size() * sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);

	/*glGenBuffers(1, &texBufID);
	glBindBuffer(GL_ARRAY_BUFFER, texBufID);
	glBufferData(GL_ARRAY_BUFFER, texBuf.size() * sizeof(float), &texBuf[0], GL_STATIC_DRAW);*/

	glGenBuffers(1, &eleBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size() * sizeof(unsigned int), &eleBuf[0], GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	assert(glGetError() == GL_NO_ERROR);
}


void SoftBody::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
	// Draw mesh
	
	prog->bind();
	//glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	//glUniform3fv(prog->getUniform("kdFront"), 1, m_color.data());
	//glUniform3fv(prog->getUniform("kdBack"), 1, Vector3f(1.0, 1.0, 0.0).data());

	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	MV->pushMatrix();
	glUniform3f(prog->getUniform("lightPos1"), 66.0, 50.0, 50.0);
	glUniform1f(prog->getUniform("intensity_1"), 0.6);
	glUniform3f(prog->getUniform("lightPos2"), -66.0, 50.0, 50.0);
	glUniform1f(prog->getUniform("intensity_2"), 0.2);
	glUniform1f(prog->getUniform("s"), 200);
	glUniform3f(prog->getUniform("ka"), 0.2, 0.2, 0.2);
	//glUniform3f(prog->getUniform("kd"), 0.8, 0.7, 0.7);
	glUniform3f(prog->getUniform("ks"), 1.0, 0.9, 0.8);
	glUniform3fv(prog->getUniform("kd"), 1, this->m_color.data());
	MV->pushMatrix();

	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	int h_pos = prog->getAttribute("aPos");
	glEnableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);

	int h_nor = prog->getAttribute("aNor");
	glEnableVertexAttribArray(h_nor);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size() * sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);

	glDrawElements(GL_TRIANGLES, 3 * m_trifaces.size(), GL_UNSIGNED_INT, (const void *)(0 * sizeof(unsigned int)));

	glDisableVertexAttribArray(h_nor);
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	for (int i = 0; i < m_attach_nodes.size(); i++) {
		auto node = m_attach_nodes[i];
		glUniform3fv(prog->getUniform("kd"), 1, node->m_color.data());
		//glUniform3fv(prog->getUniform("kdFront"), 1, node->m_color.data());
		node->draw(MV, prog);
	}

	MV->popMatrix();
	prog->unbind();

	progSimple->bind();
	for (int i = 0; i < 65; i++) {
		auto node = m_nodes[i];
		node->drawNormal(MV, P, progSimple);
	}
	progSimple->unbind();

}

void SoftBody::countDofs(int &nm, int &nr) {
	// Counts maximal and reduced DOFs
	// For non-rigid DOFs, we need both maximal and reduced DOFs,
	// and the Jacobian must pass them through with the identity 
	// matrices

	for (int i = 0; i < m_nodes.size(); i++) {
		m_nodes[i]->idxM = nm;
		m_nodes[i]->idxR = nr;
		nm += 3;
		nr += 3;
	}
}

void SoftBody::transform(Eigen::Vector3d dx) {

	for (int i = 0; i < m_nodes.size(); i++) {
		auto node = m_nodes[i];		
		node->x = node->x + dx;	
	}

}


void SoftBody::updatePosNor() {
	for (int i = 0; i < m_nodes.size(); i++) {
		auto node = m_nodes[i];
		node->clearNormals();	
	}

	for (int i = 0; i < m_trifaces.size(); i++) {
		auto triface = m_trifaces[i];

		Vector3d p0 = triface->m_nodes[0]->x;
		Vector3d p1 = triface->m_nodes[1]->x;
		Vector3d p2 = triface->m_nodes[2]->x;

		Vector3d normal = triface->computeNormal();
		
		for (int ii = 0; ii < 3; ii++) {
			posBuf[9 * i + 0 + ii] = p0(ii);
			posBuf[9 * i + 3 + ii] = p1(ii);
			posBuf[9 * i + 6 + ii] = p2(ii);
		
			if (!triface->isFlat) {
				// If the node is on the curved surface, average the normals after this loop
				auto node = triface->m_nodes[ii];
				node->addNormal(normal);
			}
			else {
				// Don't average normals if it's a flat surface
				norBuf[9 * i + 0 + ii] = normal(ii);
				norBuf[9 * i + 3 + ii] = normal(ii);
				norBuf[9 * i + 6 + ii] = normal(ii);
			}
		}
	}

	for (int i = 0; i < m_trifaces.size(); i++) {
		auto triface = m_trifaces[i];
		if (!triface->isFlat) {
			// Use the average normals if it's a curved surface
			for (int ii = 0; ii < 3; ii++) {
				// Average normals here
				Vector3d normal = triface->m_nodes[ii]->computeNormal();
				for (int iii = 0; iii < 3; iii++) {
				norBuf[9 * i + 3 * ii + iii] = normal(iii);
				}
			}
		}	
	}

}

void SoftBody::setAttachments(int id, shared_ptr<Body> body) {
	auto node = m_nodes[id];
	node->setParent(body);
	node->setColor(body->m_attached_color);

	m_attach_bodies.push_back(body);
	m_attach_nodes.push_back(node);

	Matrix4d E_ws = Matrix4d::Identity();
	E_ws.block<3, 1>(0, 3) = node->x;

	Matrix4d E_is = body->E_iw * E_ws;
	Vector3d r = E_is.block<3, 1>(0, 3);
	m_r.push_back(r);

}


void SoftBody::setAttachmentsByLine(Vector3d direction, Vector3d orig, shared_ptr<Body> body) {	
	double t, u, v;
	Vector3d xa, xb, xc, xd;
	int numIntersects = 0;

	for (int i = 0; i < m_tets.size(); i++) {
		auto tet = m_tets[i];
		xa = tet->m_nodes[0]->x;
		xb = tet->m_nodes[1]->x;
		xc = tet->m_nodes[2]->x;
		xd = tet->m_nodes[3]->x;

		if (rayTriangleIntersects(xa, xc, xb, direction, orig, t, u, v)) {
			numIntersects += 1;
			setAttachments(tet->m_nodes[0]->i, body);
			setAttachments(tet->m_nodes[2]->i, body);
			setAttachments(tet->m_nodes[1]->i, body);
		}

		if (rayTriangleIntersects(xa, xb, xd, direction, orig, t, u, v)) {
			numIntersects += 1;
			setAttachments(tet->m_nodes[0]->i, body);
			setAttachments(tet->m_nodes[3]->i, body);
			setAttachments(tet->m_nodes[1]->i, body);
		}

		if (rayTriangleIntersects(xb, xd, xc, direction, orig, t, u, v)) {
			numIntersects += 1;
			setAttachments(tet->m_nodes[3]->i, body);
			setAttachments(tet->m_nodes[2]->i, body);
			setAttachments(tet->m_nodes[1]->i, body);
		}

		if (rayTriangleIntersects(xa, xc, xd, direction, orig, t, u, v)) {
			numIntersects += 1;
			setAttachments(tet->m_nodes[0]->i, body);
			setAttachments(tet->m_nodes[2]->i, body);
			setAttachments(tet->m_nodes[3]->i, body);
		}

	}



}


VectorXd SoftBody::gatherDofs(VectorXd y, int nr) {
	// Gathers qdot and qddot into y
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		int idxR = m_nodes[i]->idxR;
		y.segment<3>(idxR) = m_nodes[i]->x;
		y.segment<3>(nr + idxR) = m_nodes[i]->v;
	}

	if (next != nullptr) {
		y = next->gatherDofs(y, nr);
	}

	return y;
}

VectorXd SoftBody::gatherDDofs(VectorXd ydot, int nr) {
	// Gathers qdot and qddot into ydot
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		int idxR = m_nodes[i]->idxR;
		ydot.segment<3>(idxR) = m_nodes[i]->v;
		ydot.segment<3>(nr + idxR) = m_nodes[i]->a;
	}

	if (next != nullptr) {
		ydot = next->gatherDDofs(ydot, nr);
	}

	return ydot;
}

void SoftBody::scatterDofs(VectorXd &y, int nr) {
	// Scatters q and qdot from y

	for (int i = 0; i < (int)m_nodes.size(); i++) {
		int idxR = m_nodes[i]->idxR;
		if (!m_nodes[i]->fixed) {
			m_nodes[i]->x = y.segment<3>(idxR);
			m_nodes[i]->v = y.segment<3>(nr + idxR);
		}

	}
	updatePosNor();

	if (next != nullptr) {
		next->scatterDofs(y, nr);
	}
}

void SoftBody::scatterDDofs(VectorXd &ydot, int nr) {
	// Scatters qdot and qddot from ydot
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		int idxR = m_nodes[i]->idxR;
		if (!m_nodes[i]->fixed) {
			m_nodes[i]->v = ydot.segment<3>(idxR);
			m_nodes[i]->a = ydot.segment<3>(nr + idxR);
		}
	}

	if (next != nullptr) {
		next->scatterDDofs(ydot, nr);

	}

}

void SoftBody::computeEnergies(Vector3d grav, double &T, double &V) {



}

MatrixXd SoftBody::computeMass(Vector3d grav, MatrixXd M) {
	// Computes maximal mass matrix
	
	Matrix3d I3 = Matrix3d::Identity();
	for (int i = 0; i < m_nodes.size(); i++) {
		int idxM = m_nodes[i]->idxM;
		double m = m_nodes[i]->m;
		M.block<3, 3>(idxM, idxM) = m * I3;
	}

	if (next != nullptr) {
		M = next->computeMass(grav, M);
	}

	return M;
}

VectorXd SoftBody::computeForce(Vector3d grav, VectorXd f) {
	// Computes force vector
	
	// Gravity
	Matrix3d I3 = Matrix3d::Identity();
	for (int i = 0; i < m_nodes.size(); i++) {
		int idxM = m_nodes[i]->idxM;
		double m = m_nodes[i]->m;
		f.segment<3>(idxM) += m * grav;
	}

	// Elastic Forces
	for (int i = 0; i < m_tets.size(); i++) {
		auto tet = m_tets[i];
		f = tet->computeElasticForces(f);


	}

	if (next != nullptr) {
		f = next->computeForce(grav, f);
	}

	return f;
}

MatrixXd SoftBody::computeStiffness(MatrixXd K) {
	VectorXd df(3*m_nodes.size());
	VectorXd Dx = df;

	for (int i = 0; i < m_tets.size(); i++) {
		auto tet = m_tets[i];
		for (int ii = 0; ii < 4; ii++) {
			auto node = tet->m_nodes[ii];
			int id = node->i;
			int col = node->idxM;
			
			for (int iii = 0; iii < 3; iii++) {
				df.setZero();
				Dx.setZero();
				Dx(3 * id + iii) = 1.0;
				tet->computeForceDifferentials(Dx, df);
				//K.col(col + iii) += df;
				K.block(col - 3 * id, col + iii, 3 * m_nodes.size(), 1) += df;
			}
		}
	}


	if (next != nullptr) {
		K = next->computeStiffness(K);
	}
	return K;
}

MatrixXd SoftBody::computeJacobian(MatrixXd J) {
	for (int i = 0; i < m_nodes.size(); i++) {
		J.block<3, 3>(m_nodes[i]->idxM, m_nodes[i]->idxR) = Matrix3d::Identity();
	}

	if (next != nullptr) {
		J = next->computeJacobian(J);
	}

	return J;
}

SoftBody:: ~SoftBody() {



}