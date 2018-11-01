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

void eigen_sym(Eigen::Matrix3d &a, Eigen::Vector3d &eig_val, Eigen::Matrix3d &eig_vec) {
	Eigen::EigenSolver<Eigen::Matrix3d> es(a);

	for (int i = 0; i < 3; i++) {
		std::complex<double> ev = es.eigenvalues()[i];
		eig_val(i) = ev.real();
		Eigen::Vector3cd v = es.eigenvectors().col(i);
		eig_vec(0, i) = v(0).real();
		eig_vec(1, i) = v(1).real();
		eig_vec(2, i) = v(2).real();
	}


}

int SVD(Eigen::Matrix3d &F,
	Eigen::Matrix3d &U,
	Eigen::Vector3d &Sigma,
	Eigen::Matrix3d &V,
	double sv_eps,
	int modifiedSVD) {

	// Adapted from Jernej Barbic's code
	// https://github.com/starseeker/VegaFEM/blob/master/libraries/minivector/mat3d.cpp

	// The code handles the following special situations:

	//---------------------------------------------------------
	// 1. det(V) == -1
	//    - multiply the first column of V by -1
	//---------------------------------------------------------
	// 2. An entry of Sigma is near zero
	//---------------------------------------------------------
	// (if modifiedSVD == 1) :
	// 3. negative determinant (Tet is inverted in solid mechanics).
	//    - check if det(U) == -1
	//    - If yes, then negate the minimal element of Sigma
	//      and the corresponding column of U
	//---------------------------------------------------------

	// form F^T F and do eigendecomposition

	Eigen::Matrix3d normalEq = F.transpose() * F;
	Eigen::Vector3d eigenValues;
	Eigen::Matrix3d eigenVectors;

	eigen_sym(normalEq, eigenValues, eigenVectors);
	V = eigenVectors;

	// Handle situation:
	// 1. det(V) == -1
	//    - multiply the first column of V by -1
	if (V.determinant() < -0.0000001) {
		// convert V into a rotation (multiply column 1 by -1)
		V.col(0) *= -1.0;
	}

	Sigma(0) = (eigenValues(0) > 0.0) ? sqrt(eigenValues(0)) : 0.0;
	Sigma(1) = (eigenValues(1) > 0.0) ? sqrt(eigenValues(1)) : 0.0;
	Sigma(2) = (eigenValues(2) > 0.0) ? sqrt(eigenValues(2)) : 0.0;

	// compute inverse of singular values
	// also check if singular values are close to zero
	Eigen::Vector3d SigmaInverse;
	SigmaInverse(0) = (Sigma(0) > sv_eps) ? (1.0 / Sigma(0)) : 0.0;
	SigmaInverse(1) = (Sigma(1) > sv_eps) ? (1.0 / Sigma(1)) : 0.0;
	SigmaInverse(2) = (Sigma(2) > sv_eps) ? (1.0 / Sigma(2)) : 0.0;

	// compute U using the formula:
	// U = F * V * diag(SigmaInverse)
	U = F * V;
	U = U * SigmaInverse.asDiagonal();

	// In theory, U is now orthonormal, U^T U = U U^T = I .. it may be a rotation or a reflection, depending on F.
	// But in practice, if singular values are small or zero, it may not be orthonormal, so we need to fix it.
	// Handle situation:
	// 2. An entry of Sigma is near zero
	// ---------------------------------------------------------
	if ((Sigma(0) < sv_eps) && (Sigma(1) < sv_eps) && (Sigma(2) < sv_eps))
	{
		// extreme case, all singular values are small, material has collapsed almost to a point
		// see [Irving 04], p. 4
		U.setIdentity();
	}
	else
	{
		// handle the case where two singular values are small, but the third one is not
		// handle it by computing two (arbitrary) vectors orthogonal to the eigenvector for the large singular value
		int done = 0;
		for (int dim = 0; dim<3; dim++)
		{
			int dimA = dim;
			int dimB = (dim + 1) % 3;
			int dimC = (dim + 2) % 3;
			if ((Sigma(dimB) < sv_eps) && (Sigma(dimC) < sv_eps))
			{
				// only the column dimA can be trusted, columns dimB and dimC correspond to tiny singular values
				Eigen::Vector3d tmpVec1 = U.col(dimA); // column dimA
				Eigen::Vector3d tmpVec2;
				tmpVec2 = findOrthonormalVector(tmpVec1);
				Eigen::Vector3d tmpVec3 = tmpVec1.cross(tmpVec2).normalized();
				U.col(dimB) = tmpVec2;
				U.col(dimC) = tmpVec3;
				
				if (U.determinant() < -0.0000001)
				{	
					U.col(dimB) *= -1.0;
				}
				done = 1;
				break; // out of for
			}
		}


		// handle the case where one singular value is small, but the other two are not
		// handle it by computing the cross product of the two eigenvectors for the two large singular values
		if (!done)
		{
			for (int dim = 0; dim<3; dim++)
			{
				int dimA = dim;
				int dimB = (dim + 1) % 3;
				int dimC = (dim + 2) % 3;

				if (Sigma(dimA) < sv_eps)
				{
					// columns dimB and dimC are both good, but column dimA corresponds to a tiny singular value
					Eigen::Vector3d tmpVec1 = U.col(dimB); // column dimB
					Eigen::Vector3d tmpVec2 = U.col(dimC); // column dimC
					Eigen::Vector3d tmpVec3 = tmpVec1.cross(tmpVec2).normalized();
					U.col(dimA) = tmpVec3;
				
					if (U.determinant() < -0.0000001)
					{	
						U.col(dimA) *= -1.0;
					}

					done = 1;
					break; // out of for
				}
			}
		}

		if ((!done) && (modifiedSVD == 1))
		{
			// Handle situation:
			// 3. negative determinant (Tet is inverted in solid mechanics)
			//    - check if det(U) == -1
			//    - If yes, then negate the minimal element of Sigma
			//      and the corresponding column of U

			double detU = U.determinant();
			if (detU < -0.0000001)
			{
				// negative determinant
				// find the smallest singular value (they are all non-negative)
				int smallestSingularValueIndex = 0;
				for (int dim = 1; dim<3; dim++)
					if (Sigma(dim) < Sigma(smallestSingularValueIndex))
						smallestSingularValueIndex = dim;

				// negate the smallest singular value
				Sigma(smallestSingularValueIndex) *= -1.0;
				U.col(smallestSingularValueIndex) *= -1.0;
				
			}
		}
	}

	return 0;
}

Eigen::Vector3d findOrthonormalVector(Eigen::Vector3d input) {
	// Find a unit vector that is orthogonal to an input vector v

	// find smallest abs component of v
	int smallestIndex = 0;
	for (int dim = 1; dim<3; dim++)
		if (fabs(input[dim]) < fabs(input[smallestIndex]))
			smallestIndex = dim;

	Eigen::Vector3d axis;
	axis.setZero();
	axis(smallestIndex) = 1.0;

	// this cross-product will be non-zero (as long as v is not zero)
	Eigen::Vector3d result = input.cross(axis).normalized();
	return result;

}