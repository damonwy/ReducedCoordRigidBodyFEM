#pragma once
#ifndef MUSCLEMASS_SRC_FACE_H_
#define MUSCLEMASS_SRC_FACE_H_
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "MLCommon.h"

class MatrixStack;
class Program;
class Node;

class Face
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Face();
	Face(std::vector<std::shared_ptr<Node>> nodes);
	virtual ~Face();

	virtual void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;

	virtual Eigen::Vector3d computeNormal();
	virtual double computeArea();

	std::vector<std::shared_ptr<Node>> m_nodes;
	Eigen::Vector3d m_normal;
	double m_area;

};

#endif // MUSCLEMASS_SRC_FACE_H_