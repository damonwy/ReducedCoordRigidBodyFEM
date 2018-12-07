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
#include "Vector.h"
#include "TetrahedronCorotational.h"
#include "TetrahedronInvertible.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

SoftBody::SoftBody(): m_isInvertible(true), m_isGravity(false), m_isElasticForce(true){
	m_color << 1.0f, 1.0f, 0.0f;
	m_isInverted = false;
}

SoftBody::SoftBody(double density, double young, double poisson, Material material) :
	m_density(density), m_young(young), m_poisson(poisson), m_material(material), 
	m_isInvertible(true), m_isGravity(false), m_isElasticForce(true)
{
	m_color << 1.0f, 1.0f, 0.0f;
	m_isInverted = false;
	//m_isGravity = true;
	m_type = 0;
}

void SoftBody::load(const string &RESOURCE_DIR, const string &MESH_NAME) {

	// Tetrahedralize 3D mesh
	tetgenio input_mesh, output_mesh;
	input_mesh.load_ply((char *)(RESOURCE_DIR + MESH_NAME).c_str());
	tetrahedralize("pqzRa20.0", &input_mesh, &output_mesh);

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
		node->load(RESOURCE_DIR);
		m_nodes.push_back(node);
	}

	// Create Faces
	for (int i = 0; i < output_mesh.numberoftrifaces; i++) {
		auto triface = make_shared<FaceTriangle>();

		for (int ii = 0; ii < 3; ii++) {
			auto node = m_nodes[output_mesh.trifacelist[3 * i + ii]];
			node->m_nfaces++;
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

		shared_ptr<Tetrahedron> tet;
		if (m_type == 0) {
			tet = make_shared<Tetrahedron>(m_young, m_poisson, m_density, m_material, tet_nodes);
		}
		else if (m_type == 1) {
			tet = make_shared<TetrahedronInvertible>(m_young, m_poisson, m_density, m_material, tet_nodes);
		}
		else
		{
			tet = make_shared<TetrahedronCorotational>(m_young, m_poisson, m_density, m_material, tet_nodes);
		}

		tet->i = i;
		tet->precompute();
		//tet->setInvertiblity(m_isInvertible);
		m_tets.push_back(tet);
	}

	// Fix the normal of top and bottom surface
	for (int i = 0; i < (int)m_trifaces.size(); i++) {
		auto triface = m_trifaces[i];

		Vector3d p0 = triface->m_nodes[0]->x;
		Vector3d p1 = triface->m_nodes[1]->x;
		Vector3d p2 = triface->m_nodes[2]->x;

		Vector3d normal = triface->computeNormal();

		if ((normal - m_trifaces[0]->m_normal).norm() < 0.1) {
			//triface->isFlat = true;
		}
		if ((normal + m_trifaces[0]->m_normal).norm() < 0.1) {
			//triface->isFlat = true;
		}
	}
}

void SoftBody::init() {

	for (int i = 0; i < (int)m_nodes.size(); i++) {
		m_nodes[i]->init();
	}

	for (int i = 0; i < (int)m_compared_nodes.size(); ++i) {
		//m_compared_nodes[i]->init();
	}

	for (int i = 0; i < (int)m_tets.size(); ++i) {
		auto tet = m_tets[i];
		for (int j = 0; j < (int)tet->m_enclosed_points.size(); ++j) {
			//tet->m_enclosed_points[j]->init();

		}
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

	for (int i = 0; i < (int)m_trifaces.size(); i++) {
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
	draw_(MV, prog, progSimple, P);
	if (next != nullptr) {
		next->draw(MV, prog, progSimple, P);
	}
}

void SoftBody::draw_(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
	// Draw mesh
	prog->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	glUniform3f(prog->getUniform("lightPos1"), 66.0f, 50.0f, 50.0f);
	glUniform1f(prog->getUniform("intensity_1"), 0.6f);
	glUniform3f(prog->getUniform("lightPos2"), -66.0f, 50.0f, 50.0f);
	glUniform1f(prog->getUniform("intensity_2"), 0.2f);
	glUniform1f(prog->getUniform("s"), 200.0f);
	glUniform3f(prog->getUniform("ka"), 0.2f, 0.2f, 0.2f);
	glUniform3f(prog->getUniform("ks"), 1.0f, 0.9f, 0.8f);
	glUniform3fv(prog->getUniform("kd"), 1, this->m_color.data());

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

	for (int i = 0; i < (int)m_attach_nodes.size(); i++) {
		//auto node = m_attach_nodes[i];
		//glUniform3fv(prog->getUniform("kd"), 1, node->m_color.data());
		//node->draw(MV, prog);
	}

	for (int i = 0; i < (int)m_sliding_nodes.size(); i++) {
		//auto node = m_sliding_nodes[i];
		//glUniform3fv(prog->getUniform("kd"), 1, node->m_color.data());
		//node->draw(MV, prog);
	}

	for (int i = 0; i < (int)m_compared_nodes.size(); i++) {
		//auto node = m_compared_nodes[i];
		//glUniform3fv(prog->getUniform("kd"), 1, node->m_color.data());
		//node->draw(MV, prog);
	}

	for (int i = 0; i < (int)m_tets.size(); ++i) {
		auto tet = m_tets[i];
		tet->draw(MV, prog, progSimple, P);
	}
	prog->unbind();

	// Draw the normals of the sliding nodes
	for (int i = 0; i < (int)m_normals_sliding.size(); ++i) {
		//auto vec = m_normals_sliding[i];
		//vec->draw(MV, P, progSimple);
	}

}

void SoftBody::countDofs(int &nm, int &nr) {
	// Counts maximal and reduced DOFs
	// For non-rigid DOFs, we need both maximal and reduced DOFs,
	// and the Jacobian must pass them through with the identity 
	// matrices

	for (int i = 0; i < (int)m_nodes.size(); i++) {
		m_nodes[i]->idxM = nm;
		m_nodes[i]->idxR = nr;
		nm += 3;
		nr += 3;
	}
}

void SoftBody::transform(Vector3d dx) {
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		auto node = m_nodes[i];
		node->x = node->x + dx;
	}
}

void SoftBody::transform(Matrix4d E) {
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		auto node = m_nodes[i];
		node->update(E);
	}
}

void SoftBody::updatePosNor() {
	// update normals
	for (int i = 0; i < (int)m_normals_sliding.size(); ++i) {
		auto vec = m_normals_sliding[i];
		vec->update();
	}

	for (int i = 0; i < (int)m_nodes.size(); i++) {
		auto node = m_nodes[i];
		node->clearNormals();
	}

	for (int i = 0; i < (int)m_trifaces.size(); i++) {
		auto triface = m_trifaces[i];

		Vector3d p0 = triface->m_nodes[0]->x;
		Vector3d p1 = triface->m_nodes[1]->x;
		Vector3d p2 = triface->m_nodes[2]->x;

		Vector3d normal = triface->computeNormal();

		for (int ii = 0; ii < 3; ii++) {
			posBuf[9 * i + 0 + ii] = float(p0(ii));
			posBuf[9 * i + 3 + ii] = float(p1(ii));
			posBuf[9 * i + 6 + ii] = float(p2(ii));

			if (!triface->isFlat) {
				// If the node is on the curved surface, average the normals after this loop
				auto node = triface->m_nodes[ii];
				node->addNormal(normal);
			}
			else {
				// Don't average normals if it's a flat surface
				norBuf[9 * i + 0 + ii] = float(normal(ii));
				norBuf[9 * i + 3 + ii] = float(normal(ii));
				norBuf[9 * i + 6 + ii] = float(normal(ii));
			}
		}
	}

	for (int i = 0; i < (int)m_trifaces.size(); i++) {
		auto triface = m_trifaces[i];
		if (!triface->isFlat) {
			// Use the average normals if it's a curved surface
			for (int ii = 0; ii < 3; ii++) {
				// Average normals here
				Vector3d normal = triface->m_nodes[ii]->computeNormal();
				for (int iii = 0; iii < 3; iii++) {
					norBuf[9 * i + 3 * ii + iii] = float(normal(iii));
				}
			}
		}
	}
}

void SoftBody::setAttachments(int id, shared_ptr<Body> body) {
	auto node = m_nodes[id];
	node->setParent(body);
	node->setColor(body->m_attached_color);
	node->r = 0.1;
	m_attach_bodies.push_back(body);
	m_attach_nodes.push_back(node);

	Matrix4d E_ws = Matrix4d::Identity();
	E_ws.block<3, 1>(0, 3) = node->x;

	Matrix4d E_is = body->E_iw * E_ws;
	Vector3d r = E_is.block<3, 1>(0, 3);
	m_r.push_back(r);
}

void SoftBody::setSlidingNodes(int id, std::shared_ptr<Body> body, Eigen::Vector3d init_dir) {
	auto node = m_nodes[id];
	node->setParent(body);
	node->setColor(body->m_sliding_color);
	node->r = 0.1;
	m_sliding_bodies.push_back(body);
	m_sliding_nodes.push_back(node);

	Matrix4d E_ws = Matrix4d::Identity();
	E_ws.block<3, 1>(0, 3) = node->x;

	Matrix4d E_is = body->E_iw * E_ws;
	Vector3d r = E_is.block<3, 1>(0, 3);
	m_r_sliding.push_back(r);

	// add normals
	auto vec = make_shared<Vector>(node, body, init_dir);
	m_normals_sliding.push_back(vec);
}

void SoftBody::setAttachmentsByLine(Vector3d direction, Vector3d orig, shared_ptr<Body> body) {
	double t, u, v;
	Vector3d xa, xb, xc, xd;
	int numIntersects = 0;

	for (int i = 0; i < (int)m_tets.size(); i++) {
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

void SoftBody::setAttachmentsByXYSurface(double z, double range, Vector2d xrange, Vector2d yrange, shared_ptr<Body> body) {
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		auto node = m_nodes[i];
		Vector3d xi = node->x;

		if (abs(xi(2) - z) < range && xi(0) <= xrange(1) && xi(0) >= xrange(0) && xi(1) <= yrange(1) && xi(1) >= yrange(0)) {
			setAttachments(i, body);
		}
	}
}

void SoftBody::setAttachmentsByYZSurface(double x, double range, Vector2d yrange, Vector2d zrange, shared_ptr<Body> body) {
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		auto node = m_nodes[i];
		Vector3d xi = node->x;

		if (abs(xi(0) - x) < range && xi(2) <= zrange(1) && xi(2) >= zrange(0) && xi(1) <= yrange(1) && xi(1) >= yrange(0)) {
			setAttachments(i, body);
		}
	}

}

void SoftBody::setAttachmentsByYZCircle(double x, double range, Vector2d O, double r, shared_ptr<Body> body) {
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		auto node = m_nodes[i];
		Vector3d xi = node->x;
		double diff = pow((xi(1) - O(0)), 2) + pow((xi(2) - O(1)), 2) - r * r;

		if (abs(xi(0) - x) < range && diff < 0.0001) {
			setAttachments(i, body);
		}
	}
}

void SoftBody::setAttachmentsByXZCircle(double y, double range, Vector2d O, double r, shared_ptr<Body> body) {
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		auto node = m_nodes[i];
		Vector3d xi = node->x;
		double diff = pow((xi(0) - O(0)), 2) + pow((xi(2) - O(1)), 2) - r * r;

		if (abs(xi(1 ) - y) < range && diff < 0.0001) {
			setAttachments(i, body);
		}
	}
}

void SoftBody::setAttachmentsByXZSurface(double y, double range, Eigen::Vector2d xrange, Eigen::Vector2d zrange, shared_ptr<Body> body) {
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		auto node = m_nodes[i];
		Vector3d xi = node->x;

		if (abs(xi(1) - y) < range && xi(0) <= xrange(1) && xi(0) >= xrange(0) && xi(2) <= zrange(1) && xi(2) >= zrange(0)) {
			setAttachments(i, body);
		}
	}

}

void SoftBody::setSlidingNodesByXYSurface(double z, Eigen::Vector2d xrange, Eigen::Vector2d yrange, double dir, std::shared_ptr<Body> body) {
	Vector3d z_axis;
	z_axis << 0.0, 0.0, 1.0;
	z_axis *= dir;

	for (int i = 0; i < (int)m_nodes.size(); i++) {
		auto node = m_nodes[i];
		Vector3d xi = node->x;

		if (abs(xi(2) - z) < 0.0001 && xi(0) <= xrange(1) && xi(0) >= xrange(0) && xi(1) <= yrange(1) && xi(1) >= yrange(0)) {
			setSlidingNodes(i, body, z_axis);
		}
	}
}

void SoftBody::setSlidingNodesByYZSurface(double x, Eigen::Vector2d yrange, Eigen::Vector2d zrange, double dir, std::shared_ptr<Body> body) {
	Vector3d x_axis;
	x_axis << 1.0, 0.0, 0.0;
	x_axis *= dir;

	for (int i = 0; i < (int)m_nodes.size(); i++) {
		auto node = m_nodes[i];
		Vector3d xi = node->x;

		if (abs(xi(0) - x) < 0.0001 && xi(2) <= zrange(1) && xi(2) >= zrange(0) && xi(1) <= yrange(1) && xi(1) >= yrange(0)) {
			setSlidingNodes(i, body, x_axis);
		}
	}
}

void SoftBody::setSlidingNodesByYZCircle(double x, double range_x, Eigen::Vector2d O, double r, std::shared_ptr<Body> body) {
	Vector3d x_axis;

	for (int i = 0; i < (int)m_nodes.size(); i++) {
		auto node = m_nodes[i];
		Vector3d xi = node->x;
		x_axis << 0.0, O(0) - xi(1), O(1) - xi(2);
		x_axis.normalized();

		double diff = pow((xi(1) - O(0)), 2) + pow((xi(2) - O(1)), 2) - r * r;
		if (abs(xi(0) - x) < range_x && diff < 0.01) {
			setSlidingNodes(i, body, x_axis);//1.7

			// Create Attached node to compare if they are really sliding
			auto comp_node = make_shared<Node>();
			comp_node->x0 = m_nodes[i]->x0;
			comp_node->x0(0) -= 5.0;
			//comp_node->x = m_nodes[i]->x;
			comp_node->r = m_nodes[i]->r;
			comp_node->sphere = m_nodes[i]->sphere;
			comp_node->setParent(body);
			m_compared_nodes.push_back(comp_node);
		}
	}
}

void SoftBody::setSlidingNodesByXZSurface(double y, Eigen::Vector2d xrange, Eigen::Vector2d zrange, double dir, std::shared_ptr<Body> body) {
	Vector3d y_axis;
	y_axis << 0.0, 1.0, 0.0;
	y_axis *= dir;
	
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		auto node = m_nodes[i];
		Vector3d xi = node->x;

		if (abs(xi(1) - y) < 0.0001 && xi(0) <= xrange(1) && xi(0) >= xrange(0) && xi(2) <= zrange(1) && xi(2) >= zrange(0)) {
			setSlidingNodes(i, body, y_axis);
		}
	}
}

void SoftBody::gatherDofs(VectorXd &y, int nr) {
	// Gathers qdot and qddot into y
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		int idxR = m_nodes[i]->idxR;
		y.segment<3>(idxR) = m_nodes[i]->x;
		y.segment<3>(nr + idxR) = m_nodes[i]->v;
	}

	if (next != nullptr) {
		next->gatherDofs(y, nr);
	}
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

	// Update points
	for (int i = 0; i < (int)m_compared_nodes.size(); ++i) {
		m_compared_nodes[i]->update();

	}

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

void SoftBody::computeMass(MatrixXd &M) {
	// Computes maximal mass matrix

	Matrix3d I3 = Matrix3d::Identity();
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		int idxM = m_nodes[i]->idxM;
		double m = m_nodes[i]->m;

		M.block<3, 3>(idxM, idxM) = m * I3;
	}

	if (next != nullptr) {
		next->computeMass(M);
	}
}

void SoftBody::computeMassSparse(vector<T> &M_) {
	Matrix3d I3 = Matrix3d::Identity();
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		int idxM = m_nodes[i]->idxM;
		double m = m_nodes[i]->m;		 
		for (int j = 0; j < 3; ++j) {
			M_.push_back(T(idxM + j, idxM + j, m));
		}
	}

	if (next != nullptr) {
		next->computeMassSparse(M_);
	}
}

void SoftBody::computeForce(Vector3d grav, VectorXd &f) {

	computeForce_(grav, f);
	if (next != nullptr) {
		next->computeForce(grav, f);
	}
}

void SoftBody::computeForce_(Vector3d grav, VectorXd &f) {
	// Computes force vector

	if (m_isGravity) {
		for (int i = 0; i < (int)m_nodes.size(); i++) {
			int idxM = m_nodes[i]->idxM;
			double m = m_nodes[i]->m;

			f.segment<3>(idxM) += m * grav;
		}
	}

	// Elastic Forces
	if (m_isElasticForce) {
		for (int i = 0; i < (int)m_tets.size(); i++) {
			auto tet = m_tets[i];
			f = tet->computeElasticForces(f);
		}
	}
}

void SoftBody::computeStiffness(MatrixXd &K) {
	computeStiffness_(K);

	if (next != nullptr) {
		next->computeStiffness(K);
	}
}

void SoftBody::computeStiffness_(MatrixXd &K) {
	for (int i = 0; i < (int)m_tets.size(); i++) {
		auto tet = m_tets[i];
		tet->computeForceDifferentials(K);
	}
}

void SoftBody::computeStiffnessSparse(vector<T> &K_) {
	computeStiffnessSparse_(K_);

	if (next != nullptr) {
		next->computeStiffnessSparse(K_);
	}
}

void SoftBody::computeStiffnessSparse_(vector<T> &K_) {
	VectorXd df(3 * m_nodes.size());
	VectorXd Dx = df;

	for (int i = 0; i < (int)m_tets.size(); i++) {
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
				int irow = col - 3 * id;
				int icol = col + iii;
				tet->computeForceDifferentialsSparse(Dx, irow, icol, K_);
			}
		}
	}
}

void SoftBody::computeJacobian(MatrixXd &J) {
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		J.block<3, 3>(m_nodes[i]->idxM, m_nodes[i]->idxR) = Matrix3d::Identity();
	}

	if (next != nullptr) {
		next->computeJacobian(J);
	}
}

void SoftBody::computeJacobianSparse(vector<T> &J_) {
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		for (int j = 0; j < 3; ++j) {
			J_.push_back(T(m_nodes[i]->idxM + j, m_nodes[i]->idxR + j, 1.0));
		}
	}

	if (next != nullptr) {
		next->computeJacobianSparse(J_);
	}
}

Energy SoftBody::computeEnergies(Eigen::Vector3d grav, Energy ener) {
	int n_nodes = (int)m_nodes.size();

	for (int i = 0; i < n_nodes; i++) {
		Vector3d x = m_nodes[i]->x;
		Vector3d v = m_nodes[i]->v;
		double m = m_nodes[i]->m;
		ener.K = ener.K + 0.5 * m * v.dot(v);
		ener.V = ener.V - m * grav.dot(x);
	}

	for (int i = 0; i < (int)m_tets.size(); i++) {
		double vi = m_tets[i]->computeEnergy();
		ener.V = ener.V + vi;
	}

	if (next != nullptr) {
		ener = next->computeEnergies(grav, ener);
	}

	return ener;
}
