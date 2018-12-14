#pragma once
// DeformableNull 

#ifndef REDUCEDCOORD_SRC_DEFORMABLENULL_H_
#define REDUCEDCOORD_SRC_DEFORMABLENULL_H_
#define EIGEN_USE_MKL_ALL

#include "Deformable.h"
class DeformableNull : public Deformable {

public:
	DeformableNull() : Deformable() {}
	virtual ~DeformableNull() {}

};

#endif // REDUCEDCOORD_SRC_SPRINGNULL_H_