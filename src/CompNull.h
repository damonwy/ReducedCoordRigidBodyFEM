#pragma once
#define EIGEN_USE_MKL_ALL

#include "Comp.h"


class CompNull : public Comp
{
public:
	CompNull() {}
	virtual ~CompNull() {}
};