#include "FaceTriangle.h"

#include <iostream>
#include <fstream>
#include "MatrixStack.h"
#include "Program.h"
#include "Node.h"



using namespace std;
using namespace Eigen;

FaceTriangle::FaceTriangle():Face() {



}


FaceTriangle::FaceTriangle(vector<shared_ptr<Node>> nodes):Face(nodes) {



}

void FaceTriangle::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {



}

void FaceTriangle::update() {
	computeNormal();
	computeArea();
}

Eigen::Vector3d FaceTriangle::computeNormal() {
	Vector3d p0 = m_nodes[0]->x;
	Vector3d p1 = m_nodes[1]->x;
	Vector3d p2 = m_nodes[2]->x;

	Vector3d e0 = p1 - p0;
	Vector3d e1 = p2 - p0;
	Vector3d normal = e0.cross(e1);
	normal.normalize();
	this->m_normal = normal;

	return m_normal;
}


double FaceTriangle::computeArea() {

	return m_area;

}