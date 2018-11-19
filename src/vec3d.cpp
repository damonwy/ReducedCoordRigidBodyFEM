#include "vec3d.h"

/*
Given an input vector v, find a unit vector that is orthogonal to it
*/
Vec3d Vec3d::findOrthonormalVector() const
{
	// find smallest abs component of v
	int smallestIndex = 0;
	for (int dim = 1; dim<3; dim++)
		if (fabs(elt[dim]) < fabs(elt[smallestIndex]))
			smallestIndex = dim;

	Vec3d axis(0.0, 0.0, 0.0);
	axis[smallestIndex] = 1.0;

	// this cross-product will be non-zero (as long as v is not zero)
	Vec3d result = norm(cross(elt, axis));
	return result;
}
