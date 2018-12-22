#include "Skeleton.h"

using namespace std;
using namespace Eigen;

Skeleton::Skeleton(Matrix4d init_frame, double length_skeleton, int nsegments):
m_E0(init_frame), m_length_skeleton(length_skeleton), m_nsegments(nsegments)
{

	m_length_segment = length_skeleton / nsegments;



}