#pragma once
#define EIGEN_USE_MKL_ALL

#include <vector>
#include <memory>
#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "MLCommon.h"

class Skeleton {
public:
	Skeleton() {}
	Skeleton(Matrix4d init_frame, double length_skeleton, int nsegments);
	~Skeleton() {}

protected:
	Matrix4d m_E0;
	double m_length_skeleton;
	int m_nsegments;
	double m_length_segment;


};