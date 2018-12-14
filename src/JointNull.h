#pragma once

#ifndef REDUCEDCOORD_SRC_JOINTNULL_H_
#define REDUCEDCOORD_SRC_JOINTNULL_H_
#define EIGEN_USE_MKL_ALL

#include "Joint.h"

class JointNull :public Joint {
public:
	JointNull() {}
	virtual ~JointNull() {}
	void update() {}
};

#endif // REDUCEDCOORD_SRC_JOINTNULL_H_