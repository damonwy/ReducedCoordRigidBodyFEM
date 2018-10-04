#pragma once
#ifndef MUSCLEMASS_SRC_SOFTBODY_H_
#define MUSCLEMASS_SRC_SOFTBODY_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "MLCommon.h"

class MatrixStack;
class Program;
class Node;
class Body;

class SoftBody {

public:
	SoftBody();
	//SoftBody();
	virtual ~SoftBody();

	virtual void load(const std::string &RESOURCE_DIR);
	virtual void init();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;

	void computeEnergies(Eigen::Vector3d grav, double &T, double &V);

private:
	std::vector<std::shared_ptr<Node> > m_nodes;
	//std::vector<std::shared_ptr<

};


#endif // MUSCLEMASS_SRC_SOFTBODY_H_