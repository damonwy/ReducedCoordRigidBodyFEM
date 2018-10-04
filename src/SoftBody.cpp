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

}

SoftBody::SoftBody(double density, double young, double poisson):
m_density(density), m_young(young), m_poisson(poisson)
{

}

void SoftBody::load(const string &RESOURCE_DIR, const string &MESH_NAME) {

	// Tetrahedralize 3D mesh
	tetgenio input_mesh, output_mesh;
	input_mesh.load_ply((char *)(RESOURCE_DIR + MESH_NAME).c_str());
	tetrahedralize("pqz", &input_mesh, &output_mesh);

	double r = 0.01;

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

		/*if (node->x(1) > 1.3) {
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
			triface->m_nodes.push_back(m_nodes[output_mesh.trifacelist[3 * i + ii]]);
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

}

void SoftBody::init() {
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
	//glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 0.0, 0.0).data());
	//glUniform3fv(prog->getUniform("kdBack"), 1, Vector3f(1.0, 1.0, 0.0).data());
	prog->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	MV->pushMatrix();
	glUniform3f(prog->getUniform("lightPos1"), 66.0, 25.0, 25.0);
	glUniform1f(prog->getUniform("intensity_1"), 0.6);
	glUniform3f(prog->getUniform("lightPos2"), -66.0, 25.0, 25.0);
	glUniform1f(prog->getUniform("intensity_2"), 0.2);
	glUniform1f(prog->getUniform("s"), 300);
	glUniform3f(prog->getUniform("ka"), 0.2, 0.2, 0.2);
	glUniform3f(prog->getUniform("kd"), 0.8, 0.7, 0.7);
	glUniform3f(prog->getUniform("ks"), 1.0, 0.9, 0.8);

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
	MV->popMatrix();
	prog->unbind();

}

void SoftBody::updatePosNor() {
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

			norBuf[9 * i + 0 + ii] = normal(ii);
			norBuf[9 * i + 3 + ii] = normal(ii);
			norBuf[9 * i + 6 + ii] = normal(ii);
		}
	}
}

void SoftBody::computeEnergies(Vector3d grav, double &T, double &V) {





}

SoftBody:: ~SoftBody() {



}