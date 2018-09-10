// Constraint Generic constraint
// A constraint can be applied in reduced or maximal coordinates.
// Reduced: keeping a quaternion to be of unit length.
// Maximal: holding two bodies together in world space.

#pragma once
#ifndef MUSCLEMASS_SRC_CONSTRAINT_H_
#define MUSCLEMASS_SRC_CONSTRAINT_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

#include <iostream>

#include "MLCommon.h"

class SE3;
class Body;


class Constraint
{

public:
	Constraint();

	virtual ~Constraint();

	virtual void update();
	virtual void draw();

	
private:
	std::string m_name;
	int m_uid;
};



#endif MUSCLEMASS_SRC_CONSTRAINT_H_