#ifndef S4R_LEAST_NORM_SOLVE_HPP_INCLUDED
#define S4R_LEAST_NORM_SOLVE_HPP_INCLUDED

#include <S4r/Types.hpp>

namespace S4r{

int LeastNormSolve(
	size_t m, size_t n, size_t nrhs,
	doublecomplex *a, size_t lda,
	doublecomplex *b, size_t ldb
);

template<typename TA, typename TB>
int LeastNormSolve(Eigen::MatrixBase<TA> &A, Eigen::MatrixBase<TB> &b){
	int ret = LeastNormSolve(
		A.rows(), A.cols(), 1, A.data(), A.outerStride(), b.data(), b.size()
	);
	b.conservativeResize(A.cols());
	return ret;
}

int Pseudoinverse(
	size_t m, size_t n,
	const doublecomplex *A, size_t lda,
	doublecomplex *P, size_t ldp
);

} // namespace S4r

#endif // S4R_LEAST_NORM_SOLVE_HPP_INCLUDED
