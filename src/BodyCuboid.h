#pragma once
#define EIGEN_USE_MKL_ALL

#ifndef REDUCEDCOORD_SRC_BODYCUBOID_H_
#define REDUCEDCOORD_SRC_BODYCUBOID_H_

#include "Body.h"

class BodyCuboid : public Body 
{
public:
	BodyCuboid() {}
	BodyCuboid(double density, Vector3d sides);
	void computeInertia_();

protected:
	Vector3d m_sides;	// Side lengths

};

#endif // REDUCEDCOORD_SRC_BODYCUBOID_H_
