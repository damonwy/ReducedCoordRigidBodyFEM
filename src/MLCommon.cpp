#include "MLCommon.h"

glm::mat3 eigen_to_glm(const Eigen::Matrix3d &m) {
	return glm::make_mat3x3((const double *)m.data());
}

glm::mat4 eigen_to_glm(const Eigen::Matrix4d &m) {
	return glm::make_mat4x4((const double *)m.data());
}

glm::vec3 eigen_to_glm(const Eigen::Vector3d &v) {
	return glm::vec3(v[0], v[1], v[2]);
}

Eigen::Vector3d glm_to_eigen(const glm::vec3 &v) {
	return Eigen::Vector3d(v[0], v[1], v[2]);
}

Eigen::Matrix3d glm_to_eigen(const glm::mat3 &m) {
	const float *_m = glm::value_ptr(m);
	float m_cp[3 * 3];
	memcpy(m_cp, _m, sizeof(m_cp));
	Eigen::Map<Eigen::Matrix3f> result(m_cp);
	return result.cast<double>();
}

Eigen::Matrix4d glm_to_eigen(const glm::mat4 &m) {
	const float *_m = glm::value_ptr(m);
	float m_cp[4 * 4];
	memcpy(m_cp, _m, sizeof(m_cp));
	Eigen::Map<Eigen::Matrix4f> result(m_cp);
	return result.cast<double>();
}

bool rayTriangleIntersects(Eigen::Vector3d v1, Eigen::Vector3d v2, Eigen::Vector3d v3, Eigen::Vector3d dir, Eigen::Vector3d pos, double &t, double &u, double &v) {

	Eigen::Vector3d e1 = v2 - v1;
	Eigen::Vector3d e2 = v3 - v1;

	// Calculate planes normal vector
	//cross product
	Eigen::Vector3d pvec = dir.cross(e2);

	//dot product
	double det = e1.dot(pvec);

	// Ray is parallel to plane
	if (det <1e-8 && det > -1e-8) {
		return false;
	}

	double inv_det = 1 / det;

	// Distance from v1 to ray pos
	Eigen::Vector3d tvec = pos - v1;
	u = (tvec.dot(pvec))*inv_det;
	if (u < 0 || u > 1) {
		return false;
	}

	Eigen::Vector3d qvec = tvec.cross(e1);
	v = dir.dot(qvec) * inv_det;
	if (v<0 || u + v>1) {
		return false;
	}

	t = e2.dot(qvec) * inv_det;
	if (t > 1e-8) { return true; }
	return false;
}