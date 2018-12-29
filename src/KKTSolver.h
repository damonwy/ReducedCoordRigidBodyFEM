#pragma once
#define EIGEN_DONT_ALIGN_STATICALLY

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <cstddef>

#include <memory>

template<
	typename _Scalar,
	typename ASolver,
	typename AMatType = Eigen::SparseMatrix<_Scalar>,
	typename GMatType = Eigen::SparseMatrix<_Scalar>
>
class KKTMatrix : public Eigen::EigenBase< KKTMatrix<_Scalar, ASolver, GMatType> > {

public:
	typedef _Scalar Scalar;
	typedef Scalar RealScalar;
	typedef int StorageIndex;
	
	enum {
		ColsAtCompileTime = Eigen::Dynamic,
		MaxColsAtCompileTime = Eigen::Dynamic,
		IsRowMajor = false
	};

	KKTMatrix() : m_A_Solver(null_ptr), m_A(null_ptr), m_G(null_ptr) {}
	
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

	Eigen::Index rows() const { return (m_G->rows() + m_A->rows()); }
	Eigen::Index cols() const { return (m_G->rows() + m_A->rows()); }

private:
	const AMatType *m_A;
	const GMatType *m_G;
	const ASolver *m_A_Solver;

};

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
		m_n = kkt_mat.rows();
		// init KKT mat
		m_mat.setAMatrix(kkt_mat.getAMatrix(), kkt_mat.getASolver());
		m_mat.setGMatrix(kkt_mat.getGMatrix());

		m_out_solver.compute(m_mat);

		return *this; }

	template<typename MatType>
	SchurComplementPreconditioner& factorize(const MatType&mat) { 
		return *this; }

	
	inline const Vector solve(const Vector& b) const
	{
		Vector top = m_outer_solver.solve(b.topRows(m_nA));
		Vector bottom = m_D * b.bottomRows(m_nG);
		Vector all(m_n, 1);
		all.topRows(m_nA) = top;
		all.bottomRows(m_nG) = bottom;
		return all;
	}
	
	template <typename MatrixType>
	void setDMatrix(const MatrixType &D)
	{
		m_nG = D.rows();
		m_nA = m_n - m_nG;
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
