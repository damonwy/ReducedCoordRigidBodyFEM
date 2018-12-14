#pragma once
#define EIGEN_USE_MKL_ALL


#include "Constraint.h"

class ConstraintNull : public Constraint
{
public:
	ConstraintNull() : Constraint(0, 0, 0, 0) {

		m_name = "null";
	}

	virtual ~ConstraintNull() {}
};