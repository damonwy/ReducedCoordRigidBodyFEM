#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Node.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"

using namespace std;
using namespace Eigen;

Node::Node() :
	r(1.0),
	m(1.0),
	i(-1),
	x(0.0, 0.0, 0.0),
	v(0.0, 0.0, 0.0),
	fixed(false)
{
	
}

Node::Node(const shared_ptr<Shape> s) :
	r(1.0),
	m(1.0),
	i(-1),
	x(0.0, 0.0, 0.0),
	v(0.0, 0.0, 0.0),
	x_temp(0.0, 0.0, 0.0),
	fixed(false),
	normal(0.0,0.0,0.0),
	sphere(s)
{
	
}

Node::~Node()
{
}

void Node::tare()
{
	x0 = x;
	v0 = v;
}

void Node::reset()
{
	x = x0;
	v = v0;
}

void Node::update(Matrix4d E) {
	Vector4d pos;
	pos.segment<3>(0) = this->x0;
	pos(3) = 1.0;

	pos = E * pos;
	this->x = pos.segment<3>(0);
}

void Node::updateTemp(Matrix4d E) {
	Vector4d pos;
	pos.segment<3>(0) = this->x0;
	pos(3) = 1.0;

	pos = E * pos;
	this->x_temp = pos.segment<3>(0);
}

double Node::computePotentialEnergy(Vector3d grav) {
	this->V = this->m * grav.transpose() * this->x;
	return this->V;
}

double Node::computeKineticEnergy(VectorXd phi) {
	//assert(this->J.cols == phi.size());
	this->v = this->J * phi;
	this->K =0.5 * this->m * this->v.transpose() * this->v;
	return this->K;
}

void Node::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	if(sphere) {
		MV->pushMatrix();
		MV->translate(x(0), x(1), x(2));
		MV->scale(r);
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		sphere->draw(prog);
		MV->popMatrix();
	}
}
