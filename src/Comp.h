#pragma once

#ifndef MUSCLEMASS_SRC_COMP_H_
#define MUSCLEMASS_SRC_COMP_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "MLCommon.h"

class Shape;
class Body;
class Joint;
class Program;
class MatrixStack;

class Comp
{
public:
	Comp();
	virtual ~Comp();

	virtual void load(const std::string &RESOURCE_DIR, std::string shape);
	virtual void init();
	virtual void update();
	virtual void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P)const;
	
	std::shared_ptr<Comp> next;				// Next component in traversal order

private:
	
};

#endif // MUSCLEMASS_SRC_COMP_H_
