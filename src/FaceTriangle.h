#pragma once

#include "Face.h"

class FaceTriangle : public Face
{
public:
	FaceTriangle();
	FaceTriangle(std::vector<std::shared_ptr<Node>> nodes);

	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;

	void update();

	Eigen::Vector3d computeNormal();
	double computeArea();

	bool isFlat;
};