#define TETLIBRARY
#include <tetgen.h>

#include "Surface.h"
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "GLSL.h"
#include "MatrixStack.h"
#include "Program.h"
#include "Node.h"
#include "FaceTriangle.h"
#include <omp.h>
using namespace std;
using namespace Eigen;
using json = nlohmann::json;

void Surface::load(const std::string &RESOURCE_DIR, const std::string &MESH_NAME) {
	// Tetrahedralize 3D mesh
	tetgenio input_mesh, output_mesh;
	input_mesh.load_ply((char *)(RESOURCE_DIR + MESH_NAME).c_str());
	tetrahedralize("pqz", &input_mesh, &output_mesh);//

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
}


void Surface::init() {
	for (int i = 0; i < (int)m_nodes.size(); i++) {
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

void Surface::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P) const {
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
	prog->unbind();

}

void Surface::updatePosNor() {

#pragma omp parallel for num_threads(getThreadsNumber((int)m_nodes.size(), MIN_ITERATOR_NUM))
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		auto node = m_nodes[i];
		node->clearNormals();
	}


	for (int i = 0; i < (int)m_trifaces.size(); i++) {
		auto triface = m_trifaces[i];

		Vector3d p0 = triface->m_nodes[0]->x;
		Vector3d p1 = triface->m_nodes[1]->x;
		Vector3d p2 = triface->m_nodes[2]->x;
		{
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
	}

	for (int i = 0; i < (int)m_trifaces.size(); i++) {
		auto triface = m_trifaces[i];
		if (!triface->isFlat) {
			// Use the average normals if it's a curved surface
			for (int ii = 0; ii < 3; ii++) {
				// Average normals here
				{
					Vector3d normal = triface->m_nodes[ii]->computeNormal();
					for (int iii = 0; iii < 3; iii++) {
						norBuf[9 * i + 3 * ii + iii] = float(normal(iii));
					}
				}
			}
		}
	}


}

void Surface::transform(Matrix4d E) {
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		auto node = m_nodes[i];
		node->update(E);
	}
}