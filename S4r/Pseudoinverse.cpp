#include <S4r/Eigensystems.hpp>
#include <S4r/IRA.hpp>
#include <Eigen/Dense>

using S4r::doublecomplex;

typedef int integer;
extern "C" void zgelsy_(
	const integer &m, const integer &n, const integer &nrhs,
	doublecomplex *a, const integer &lda,
	doublecomplex *b, const integer &ldb,
	integer *jpvt, const double &rcond, integer *rank,
	doublecomplex *work, const integer &lwork, double *rwork,
	integer *info
);
extern "C" int zgesvd_(
	const char *jobu, const char *jobvt, const integer &m,
	const integer &n, doublecomplex *a, const integer &lda, double *s,
	doublecomplex *u, const integer &ldu, doublecomplex *vt, const integer &ldvt,
	doublecomplex *work, const integer &lwork, double *rwork, integer *info
);

namespace S4r{

int LeastNormSolve(
	size_t m, size_t n, size_t nrhs,
	doublecomplex *a, size_t lda,
	doublecomplex *b, size_t ldb
){
	static const double rcond = std::numeric_limits<double>::epsilon();
	doublecomplex dum;
	doublecomplex *work;
	integer *jpvt = (integer*)malloc(sizeof(integer) * n);
	memset(jpvt, 0, sizeof(integer) * n);
	integer rank, info;
	integer lwork = -1;
	double *rwork = (double*)malloc(sizeof(double) * 2*n);
	zgelsy_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, &rank, &dum, lwork, rwork, &info);
	lwork = (int)dum.real();
	work = (doublecomplex*)malloc(sizeof(doublecomplex) * lwork);
	zgelsy_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, &rank, work, lwork, rwork, &info);
	free(work);
	free(rwork);
	free(jpvt);
	return info;
}


int Pseudoinverse(
	size_t m, size_t n,
	const doublecomplex *A, size_t lda,
	doublecomplex *P, size_t ldp
){
	integer info;
	CMat Acopy(Eigen::Map<const CMat,Eigen::Unaligned,Eigen::OuterStride<> >(A, m, n, Eigen::OuterStride<>(lda)));
	Eigen::Map<CMat,Eigen::Unaligned,Eigen::OuterStride<> > mP(P, n, m, Eigen::OuterStride<>(ldp));
	if(m >= n){ // tall case
		RVec S(n);
		CMat VH(n,n);
		doublecomplex dum;
		integer lwork = -1;
		RVec rwork(5*n);
		zgesvd_(
			"O","A", m,n, Acopy.data(), Acopy.outerStride(),
			S.data(), NULL, m, VH.data(), VH.outerStride(),
			&dum, lwork, rwork.data(), &info
		);
		lwork = (integer)dum.real();
		CVec work(lwork);
		zgesvd_(
			"O","A", m,n, Acopy.data(), Acopy.outerStride(),
			S.data(), NULL, m, VH.data(), VH.outerStride(),
			work.data(), lwork, rwork.data(), &info
		);
		mP = Acopy.adjoint();
		
		{
			double threshold = 2 * std::numeric_limits<double>::epsilon() * S[0];
			for(size_t i = 0; i < n; ++i){
				if(S[i] < threshold){
					break;
				}
				S[i] = 1./S[i];
			}
		}
		mP = VH.adjoint() * S.asDiagonal() * mP;
	}else{ // wide case
		RVec S(m);
		CMat U(m,m);
		doublecomplex dum;
		integer lwork = -1;
		RVec rwork(5*m);
		zgesvd_(
			"A","O", m,n, Acopy.data(), Acopy.outerStride(),
			S.data(), U.data(), U.outerStride(), NULL, m,
			&dum, lwork, rwork.data(), &info
		);
		lwork = (integer)dum.real();
		CVec work(lwork);
		zgesvd_(
			"A","O", m,n, Acopy.data(), Acopy.outerStride(),
			S.data(), U.data(), U.outerStride(), NULL, m,
			work.data(), lwork, rwork.data(), &info
		);
		mP = Acopy.adjoint();
		
		{
			double threshold = 2 * std::numeric_limits<double>::epsilon() * S[0];
			for(size_t i = 0; i < m; ++i){
				if(S[i] < threshold){
					break;
				}
				S[i] = 1./S[i];
			}
		}
		mP = mP * S.asDiagonal() * U.adjoint();
	}
	return info;
}

} // namespace S4r
