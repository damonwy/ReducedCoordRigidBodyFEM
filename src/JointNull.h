#pragma once

#ifndef MUSCLEMASS_SRC_JOINTNULL_H_
#define MUSCLEMASS_SRC_JOINTNULL_H_

#include "Joint.h"

class JointNull :public Joint {
public:
	JointNull();
	virtual ~JointNull();
	void update();
};

#endif // MUSCLEMASS_SRC_JOINTNULL_H_