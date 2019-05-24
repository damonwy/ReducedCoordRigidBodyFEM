#include "Face.h"

#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

#include "Program.h"
#include "MatrixStack.h"
#include "Node.h"

using namespace std;
using namespace Eigen;

Face::Face() {

}


Face::Face(std::vector<std::shared_ptr<Node>> nodes):
m_nodes(nodes)
{

}

void Face::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {

}

Vector3d Face::computeNormal() {

	return m_normal;
}


double Face::computeArea() {

	return m_area;
}

Face:: ~Face() {

}
