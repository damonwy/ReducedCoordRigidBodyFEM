#pragma once

#define EIGEN_USE_MKL_ALL
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "MLCommon.h"
class Node;

class TetgenHelper {

public:
	static void createNodeFile(std::vector<std::shared_ptr<Node>> &node, char *filename);
	



};