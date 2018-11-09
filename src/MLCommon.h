#pragma once
#ifndef MUSCLEMASS_SRC_MLCOMMON_H_
#define MUSCLEMASS_SRC_MLCOMMON_H_

#include <json.hpp>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <complex>
#include <unsupported/Eigen/CXX11/Tensor>

typedef Eigen::Matrix<int, 6, 1> Vector6i;
typedef Eigen::Matrix<float, 3, 1> Vector3f;
typedef Eigen::Matrix<double, 2, 1> Vector2d;

typedef Eigen::Matrix<double, 3, 1> Vector3d;
typedef Eigen::Matrix<double, 4, 1> Vector4d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 12, 1> Vector12d;

typedef Eigen::Matrix<double, 3, 3> Matrix3d;
typedef Eigen::Matrix<double, 4, 4> Matrix4d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 9, 9> Matrix9d;
typedef Eigen::Matrix<double, 12, 12> Matrix12d;

typedef Eigen::Matrix<double, 3, 4> Matrix3x4d;
typedef Eigen::Matrix<double, 3, 6> Matrix3x6d;
typedef Eigen::Matrix<double, 3, 12> Matrix3x12d;

typedef Eigen::Matrix<double, 4, 2> Matrix4x2d;
typedef Eigen::Matrix<double, 4, 3> Matrix4x3d;

typedef Eigen::Matrix<double, 5, 6> Matrix5x6d;

typedef Eigen::Matrix<double, 6, 2> Matrix6x2d;
typedef Eigen::Matrix<double, 6, 3> Matrix6x3d;

typedef Eigen::TensorFixedSize<double, Eigen::Sizes<4, 4, 6>> Tensor4x4x6d;
typedef Eigen::TensorFixedSize<double, Eigen::Sizes<6, 2, 2>> Tensor6x2x2d;

enum Integrator { REDMAX_EULER, REDUCED_ODE45, REDMAX_ODE45 };
enum Material {LINEAR, CO_ROTATED, STVK, NEO_HOOKEAN, MOONEY_RIVLIN};
enum Axis {X_AXIS, Y_AXIS, Z_AXIS};

struct Energy {
	double K;
	double V;
};


// Eigen types to/from GLM types
glm::mat3 eigen_to_glm(const Eigen::Matrix3d &m);
glm::mat4 eigen_to_glm(const Eigen::Matrix4d &m);
glm::vec3 eigen_to_glm(const Eigen::Vector3d &v);
Eigen::Vector3d glm_to_eigen(const glm::vec3 &v);
Eigen::Matrix3d glm_to_eigen(const glm::mat3 &m);
Eigen::Matrix4d glm_to_eigen(const glm::mat4 &m);

bool rayTriangleIntersects(Eigen::Vector3d v1, Eigen::Vector3d v2, Eigen::Vector3d v3, Eigen::Vector3d dir, Eigen::Vector3d pos, double &t, double &u, double &v);

template<typename T>
using  MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar, int rank, typename sizeType>
auto Tensor_to_Matrix(const Eigen::Tensor<Scalar, rank> &tensor, const sizeType rows, const sizeType cols)
{
	return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), rows, cols);
}


template<typename Scalar, typename... Dims>
auto Matrix_to_Tensor(const MatrixType<Scalar> &matrix, Dims... dims)
{
	constexpr int rank = sizeof... (Dims);
	return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.data(), { dims... });
}

int SVD(Eigen::Matrix3d &F,
	Eigen::Matrix3d &U,
	Eigen::Vector3d &Sigma,
	Eigen::Matrix3d &V,
	double sv_eps,
	int modifiedSVD);

void eigen_sym(Eigen::Matrix3d &a, Eigen::Vector3d &eig_val, Eigen::Matrix3d &eig_vec);
Eigen::Vector3d findOrthonormalVector(Eigen::Vector3d input);
#endif // MUSCLEMASS_SRC_MLCOMMON_H_