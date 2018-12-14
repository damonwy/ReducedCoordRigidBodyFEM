#pragma once
// SpringNull 

#ifndef REDUCEDCOORD_SRC_SPRINGNULL_H_
#define REDUCEDCOORD_SRC_SPRINGNULL_H_
#define EIGEN_USE_MKL_ALL

#include "Spring.h"
class SpringNull : public Spring {

public:
	SpringNull() : Spring() {}
	virtual ~SpringNull() {}

};

#endif // REDUCEDCOORD_SRC_SPRINGNULL_H_