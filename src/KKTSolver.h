#pragma once
#define EIGEN_DONT_ALIGN_STATICALLY

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include <cstddef>
#include <memory>
#include <iostream>
template<
	typename _Scalar,
	typename ASolver,
	typename AMatType = Eigen::SparseMatrix<_Scalar>,
	typename GMatType = Eigen::SparseMatrix<_Scalar>
>
class KKTMatrix : public Eigen::EigenBase< KKTMatrix<_Scalar, ASolver, AMatType, GMatType> > {

public:
	typedef _Scalar Scalar;
	typedef Scalar RealScalar;
	typedef int StorageIndex;
	
	enum {
		ColsAtCompileTime = Eigen::Dynamic,
		MaxColsAtCompileTime = Eigen::Dynamic,
		IsRowMajor = false
	};

	KKTMatrix() : m_A_Solver(nullptr), m_A(nullptr), m_G(nullptr) {}
	
	template<typename Rhs>
	Eigen::Product<KKTMatrix<Scalar, ASolver>, Rhs, Eigen::AliasFreeProduct>
		operator*(const Eigen::MatrixBase<Rhs> &x) const {
		return Eigen::Product<
			KKTMatrix<Scalar, ASolver>,
			Rhs,
			Eigen::AliasFreeProduct>(*this, x.derived());
	}

	KKTMatrix & setGMatrix(const GMatType &G) {
		m_G = &G;
		return *this;
	}

	KKTMatrix & setAMatrix(const AMatType &A, const ASolver & A_Solver) {
		m_A = &A;
		m_A_Solver = &A_Solver;
		return *this;
	}

	const GMatType &  getGMatrix() const { return *m_G; }
	const AMatType &  getAMatrix() const { return *m_A; }
	const ASolver  &  getASolver() const { return *m_A_Solver; }

	bool isInitialized() const
	{
		return m_A_Solver != nullptr && m_G != nullptr && m_A != nullptr;
	}

	Eigen::Index rows() const { int i = m_A->rows(); return i; }
	Eigen::Index cols() const { return (m_A->rows()); }

private:
	const AMatType *m_A;
	const GMatType *m_G;
	const ASolver *m_A_Solver;

};

namespace Eigen {
	namespace internal {
		// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
		template<
			typename _Scalar,
			typename ASolver
			>
			struct traits<KKTMatrix<_Scalar, ASolver>> : public Eigen::internal::traits<Eigen::SparseMatrix<double> >
		{};
	}
}

// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
	namespace internal {
		template<
			typename Rhs,
			typename _Scalar,
			typename ASolver
			//typename AMatType = Eigen::SparseMatrix<_Scalar>,
			//typename GMatType = Eigen::SparseMatrix<_Scalar>
		>
		struct generic_product_impl<KKTMatrix<_Scalar, ASolver>, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
			: generic_product_impl_base<KKTMatrix<_Scalar, ASolver>, Rhs, generic_product_impl<KKTMatrix<_Scalar, ASolver>, Rhs> >
		{
			typedef typename Product<KKTMatrix<_Scalar, ASolver>, Rhs>::Scalar Scalar;
			template<typename Dest>
			static void scaleAndAddTo(Dest& dst, const KKTMatrix<_Scalar, ASolver>& lhs, const Rhs& rhs, const Scalar& alpha)
			{
				// This method should implement "dst += alpha * lhs * rhs" inplace,
				// however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
				assert(alpha == Scalar(1) && "scaling is not implemented");
				EIGEN_ONLY_USED_FOR_DEBUG(alpha);
				// Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
				// but let's do something fancier (and less efficient):
				for (Index i = 0; i < lhs.cols(); ++i) {
					VectorXd acol = lhs.getAMatrix().col(i);
					//VectorXd gcol = lhs.getGMatrix().col(i);
					//VectorXd lhscol(acol.size() + gcol.size());
					//lhscol.topRows(acol.size()) = acol;
					//lhscol.bottomRows(gcol.size()) = gcol;
					dst += rhs(i) * acol;

				}
					//dst += rhs(i) * lhs.my_matrix().col(i);
			}
		};
	}
}

template <typename _Scalar>
class SchurComplementPreconditioner
{
	typedef _Scalar Scalar;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
	typedef int StorageIndex;
public:
	enum {
		ColsAtCompileTime = Eigen::Dynamic,
		MaxColsAtCompileTime = Eigen::Dynamic
	};

	SchurComplementPreconditioner():m_isInitialized(true) {}

	template<typename MatType>
	SchurComplementPreconditioner& analyzePattern(const MatType&) { return *this; }


	template<typename MatType>
	SchurComplementPreconditioner& compute(const MatType&) { return *this; }

	template<
		typename ASolver,
		typename AMatType,
		typename GMatType>
	SchurComplementPreconditioner& compute(
		const KKTMatrix<
			Scalar,
			ASolver,
			AMatType,
			GMatType> &kkt_mat) 
	{ 
		// get diag(A)
		m_diag_precon.compute(kkt_mat.getAMatrix());
		m_nA = kkt_mat.rows();
		// init KKT mat
		m_mat.setAMatrix(kkt_mat.getAMatrix(), kkt_mat.getASolver());
		m_mat.setGMatrix(kkt_mat.getGMatrix());

		m_outer_solver.compute(m_mat);

		return *this; }

	template<typename MatType>
	SchurComplementPreconditioner& factorize(const MatType&mat) { 
		return *this; }

	Eigen::Index rows() const { return m_n; }
	Eigen::Index cols() const { return m_n; }
	
	inline const Vector solve(const Vector& b) const
	{
		Vector top = m_outer_solver.solve(b.topRows(m_nA));
		std::cout << b.topRows(m_nA) << std::endl;
		std::cout << top << std::endl;
		Vector bottom = m_D * b.bottomRows(m_nG);
		Vector all(m_n, 1);
		all.topRows(m_nA) = top;
		all.bottomRows(m_nG) = bottom;
		std::cout << all << std::endl;
		return all;
	}
	
	template <typename MatrixType>
	void setDMatrix(const MatrixType &D)
	{
		m_nG = D.rows();
		m_n = m_nA + m_nG;
		m_D = D;
	}

	Eigen::ComputationInfo info() { return Eigen::Success; }

private:
	Eigen::SparseMatrix<Scalar> m_D;
	bool m_isInitialized;
	Eigen::DiagonalPreconditioner<Scalar> m_diag_precon;
	KKTMatrix<Scalar, Eigen::DiagonalPreconditioner<Scalar>> m_mat;
	Eigen::ConjugateGradient<KKTMatrix<Scalar, Eigen::DiagonalPreconditioner<Scalar>>,
		Eigen::Lower | Eigen::Upper,
		Eigen::IdentityPreconditioner> m_outer_solver;
	StorageIndex m_nA;
	StorageIndex m_nG;
	StorageIndex m_n;
};