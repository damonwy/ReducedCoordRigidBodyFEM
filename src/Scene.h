#pragma once
#ifndef MUSCLEMASS_SRC_SCENE_H_
#define MUSCLEMASS_SRC_SCENE_H_

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <json.hpp>
#include "MLCommon.h"

class Particle;
class MatrixStack;
class Program;
class Shape;
class Rigid;
class Solver;
class Spring;
class Joint;
class Vector;
class WrapSphere;
class WrapCylinder;
class WrapDoubleCylinder;
class SymplecticIntegrator;
class Stepper;

class Scene
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW	
	Scene();
	virtual ~Scene();
	
	void load(const std::string &RESOURCE_DIR);
	void init();
	void tare();
	void reset();
	void step();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const;

	double getTime() const { return t; }

private:
	double t;
	double h;
	Eigen::Vector3d grav;

	std::shared_ptr<Shape> boxShape;

	std::vector< std::shared_ptr<Rigid> > boxes;
	std::vector< std::shared_ptr<Joint> > joints;

	std::shared_ptr<SymplecticIntegrator> symplectic_solver;

	nlohmann::json js;
	Integrator time_integrator;
};

#endif // MUSCLEMASS_SRC_SCENE_H_