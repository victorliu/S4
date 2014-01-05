#include <S4r/Eigensystems.hpp>
#include <S4r/IRA.hpp>
#include <iostream>

using S4r::doublecomplex;

typedef int integer;
extern "C" void zgeev_(
	const char *jobvl, const char *jobvr,
	const integer &n, doublecomplex *a, const integer &lda,
	doublecomplex *w,
	doublecomplex *vl, const integer &ldvl,
	doublecomplex *vr, const integer &ldvr,
	doublecomplex *work, const integer &lwork, double *rwork,
	integer *info
);
extern "C" void zheev_(
	const char *jobz, const char *uplo,
	const integer &n, doublecomplex *a, const integer &lda,
	double *w,
	doublecomplex *work, const integer &lwork, double *rwork,
	integer *info
);

namespace S4r{

void SortEigenpairs(CVec &vals, CMat &vecs);

static inline int LapackEigensystem(
	size_t n,
	doublecomplex *a, size_t lda,
	doublecomplex *eval,
	doublecomplex *vl, size_t ldvl, doublecomplex *vr, size_t ldvr,
	doublecomplex *work_ = NULL, double *rwork_ = NULL, size_t lwork_ = 0
){
	if(n == 0) {
		return 0;
	}
	char jobvl[2] = "N";
	char jobvr[2] = "N";
	if(vl != NULL){ jobvl[0] = 'V'; }
	if(vr != NULL){ jobvr[0] = 'V'; }
	integer info;

	if((size_t)-1 == lwork_){
		zgeev_(jobvl, jobvr, n, a, lda, eval, vl, ldvl, vr, ldvr, work_, lwork_, rwork_, &info);
		return 0;
	}
	
	doublecomplex *work = work_;
	double *rwork = rwork_;
	if(NULL == rwork_){
		rwork = new double[2*n];
	}
	if(0 == lwork_){ lwork_ = 2*n; }
	integer lwork = lwork_;
	if(NULL == work_ || lwork < (integer)(2*n)){
		lwork = -1;
		doublecomplex zlen;
		zgeev_(jobvl, jobvr, n, a, lda, eval, vl, ldvl, vr, ldvr, &zlen, lwork, rwork, &info);
		lwork = (integer)(zlen.real());
		work = new doublecomplex[lwork];
	}
	zgeev_(jobvl, jobvr, n, a, lda, eval, vl, ldvl, vr, ldvr, work, lwork, rwork, &info);
	//std::cout << "info = " << info << "\n";
	if(NULL == work_ || lwork < (integer)(2*n)){
		delete [] work;
	}
	if(NULL == rwork_){
		delete [] rwork;
	}
	
	return info;
}
static inline int LapackHermitianEigensystem(
	size_t n,
	doublecomplex *a, size_t lda,
	doublecomplex *eval,
	doublecomplex *work_ = NULL, double *rwork_ = NULL, size_t lwork_ = 0
){
	if(n == 0) {
		return 0;
	}
	integer info;

	if((size_t)-1 == lwork_){
		zheev_("V", "L", n, a, lda, rwork_, work_, lwork_, rwork_, &info);
		return 0;
	}
	
	doublecomplex *work = work_;
	double *rwork = rwork_;
	if(NULL == rwork_){
		rwork = new double[4*n];
	}
	if(0 == lwork_){ lwork_ = 2*n; }
	integer lwork = lwork_;
	if(NULL == work_ || lwork < (integer)(2*n)){
		lwork = -1;
		doublecomplex zlen;
		zheev_("V", "L", n, a, lda, rwork, &zlen, lwork, rwork+n, &info);
		lwork = (integer)(zlen.real());
		work = new doublecomplex[lwork];
	}
	zheev_("V", "L", n, a, lda, rwork, work, lwork, rwork+n, &info);
	for(size_t i = 0; i < n; ++i){
		eval[i] = rwork[i];
	}
	//std::cout << "info = " << info << "\n";
	if(NULL == work_ || lwork < (integer)(2*n)){
		delete [] work;
	}
	if(NULL == rwork_){
		delete [] rwork;
	}
	
	return info;
}


class ComplexEigensolver : public IRA::ComplexEigensystem{
	const SpMat &A;
public:
	ComplexEigensolver(
		const SpMat &A, size_t nwanted, size_t narnoldi,
		const IRA::ComplexEigensystem::Params &params
	):	ComplexEigensystem(
			A.rows(), nwanted, narnoldi, false, 0., params, NULL
		),
		A(A)
	{
	}
	
	// Common interface routines
	bool IsOpInPlace() const{ return false; }
	// If IsAInPlace returns true, then y contains the input
	// vector rather than x, and x is NULL.
	void ApplyOp(size_t n, const doublecomplex *x, doublecomplex*y){
		Eigen::Map<CVec>(y, n) = A * Eigen::Map<const CVec>(x, n);
	}
	
	// Returns whether Bop is the identity.
	bool IsBIdentity() const{ return true; }
	// True if the B operator can be applied in place: y <- Bop*y
	// instead of the usual y <- Bop*x
	bool IsBInPlace() const{ return true; }
	// If IsBIdentity returns true, then ApplyB will never be called.
	// If IsBInPlace returns true, then y contains the input
	// vector rather than x, and x is NULL.
	void ApplyB(size_t n, const doublecomplex *x, doublecomplex *y){ return; }
	
	// Return true to use user specified shifts (typically this is not used).
	bool CanSupplyShifts() const{ return false; }
	// Requesting n shifts returned in s. The array s0 contains the current
	// Ritz eigenvalue estimates.
	void GetShifts(size_t n, const doublecomplex *s0, doublecomplex *s) const{}
	
	// If a is more desireable than b, then Compare(a,b) should return false.
	bool EigenvalueCompare(const doublecomplex &a, const doublecomplex &b) const{
		return a.real() < b.real();
	}
};

class HermitianEigensolver : public IRA::ComplexEigensystem{
	const SpMat &A;
public:
	HermitianEigensolver(
		const SpMat &A, size_t nwanted, size_t narnoldi,
		const IRA::ComplexEigensystem::Params &params
	):	ComplexEigensystem(
			A.rows(), nwanted, narnoldi, false, 0., params, NULL
		),
		A(A)
	{
	}
	
	// Common interface routines
	bool IsOpInPlace() const{ return false; }
	// If IsAInPlace returns true, then y contains the input
	// vector rather than x, and x is NULL.
	void ApplyOp(size_t n, const doublecomplex *x, doublecomplex*y){
		Eigen::Map<CVec>(y, n) = A.selfadjointView<Eigen::Lower>() * Eigen::Map<const CVec>(x, n);
	}
	
	// Returns whether Bop is the identity.
	bool IsBIdentity() const{ return true; }
	// True if the B operator can be applied in place: y <- Bop*y
	// instead of the usual y <- Bop*x
	bool IsBInPlace() const{ return true; }
	// If IsBIdentity returns true, then ApplyB will never be called.
	// If IsBInPlace returns true, then y contains the input
	// vector rather than x, and x is NULL.
	void ApplyB(size_t n, const doublecomplex *x, doublecomplex *y){ return; }
	
	// Return true to use user specified shifts (typically this is not used).
	bool CanSupplyShifts() const{ return false; }
	// Requesting n shifts returned in s. The array s0 contains the current
	// Ritz eigenvalue estimates.
	void GetShifts(size_t n, const doublecomplex *s0, doublecomplex *s) const{
		for(size_t i = 0; i < n; ++i){
			s[i] = s0[i].real();
		}
	}
	
	// If a is more desireable than b, then Compare(a,b) should return false.
	bool EigenvalueCompare(const doublecomplex &a, const doublecomplex &b) const{
		return a.real() < b.real();
	}
};

int Eigensystem(
	SpMat::Index cols, const SpMat &A,
	CVec &vals, CMat &vecs
){
	static const SpMat::Index crossover_dense = 512;
	static const SpMat::Index crossover_dense_ira = 64;
	
	int ret;
	
	if(
		A.rows() < crossover_dense && // small enough to solve densely ...
		cols > crossover_dense_ira    // and want more columns than IRA can handle for a dense problem
	){
		// Solve the full dense problem
		CMat Ad(A);
		if(cols == A.rows()){
			vals.resize(Ad.rows());
			vecs.resize(Ad.rows(), Ad.rows());
			ret = LapackEigensystem(
				Ad.rows(), Ad.data(), Ad.outerStride(),
				vals.data(), NULL, Ad.rows(), vecs.data(), vecs.outerStride()
			);
		}else{
			vals.resize(Ad.rows());
			vecs.resize(Ad.rows(), Ad.cols());
			ret = LapackEigensystem(
				Ad.rows(), Ad.data(), Ad.outerStride(),
				vals.data(), NULL, Ad.rows(), vecs.data(), vecs.outerStride()
			);
			SortEigenpairs(vals, vecs);
			vals.conservativeResize(cols);
			vecs.conservativeResize(Ad.rows(), cols);
		}
	}else{
		// The number of Arnoldi vectors should be larger than the number of
		// desired eigenpairs by an amount that is greater than the maximum
		// degeneracy of the eigenvalues. We know that the maximally symmetric
		// lattices have a 6-fold degeneracy, times 2 for polarization. We round
		// up to 16 to be safe.
		SpMat::Index n_arnoldi = (cols > 16 ? 2*cols+16 : cols + 16);
		if(n_arnoldi > A.rows()){ n_arnoldi = A.rows(); }
		
		IRA::ComplexEigensystem::Params params;
		params.max_iterations = 4*A.rows();
		if((size_t)(4*cols) < params.max_iterations){
			params.max_iterations = 4*cols;
		}
		params.tol = 2. * std::numeric_limits<double>::epsilon();
		
		ComplexEigensolver solver(A, cols, n_arnoldi, params);
		SpMat::Index nconv = solver.GetConvergedCount();
		if(nconv >= cols){
			vals = Eigen::Map<const CVec>(solver.GetEigenvalues(), cols);
			vecs = Eigen::Map<const CMat>(solver.GetEigenvectors(), A.rows(), cols);
			SortEigenpairs(vals, vecs);
			ret = 0;
		}else{
			vals.resize(cols); vals.setZero();
			vecs.resize(A.rows(), cols); vecs.setZero();
			ret = (int)nconv - (int)cols;
		}
	}
	return ret;
}

int HermitianEigensystem(
	SpMat::Index cols, const SpMat &A,
	CVec &vals, CMat &vecs
){
	static const SpMat::Index crossover_dense = 512;
	static const SpMat::Index crossover_dense_ira = 64;
	
	int ret;
	
	if(
		A.rows() < crossover_dense && // small enough to solve densely ...
		cols > crossover_dense_ira    // and want more columns than IRA can handle for a dense problem
	){
		// Solve the full dense problem
		vals.resize(A.rows());
		vecs = A;
		ret = LapackHermitianEigensystem(
			A.rows(), vecs.data(), vecs.outerStride(),
			vals.data()
		);
		if(cols < A.rows()){
			SortEigenpairs(vals, vecs);
			vals.conservativeResize(cols);
			vecs.conservativeResize(A.rows(), cols);
		}
	}else{
		// The number of Arnoldi vectors should be larger than the number of
		// desired eigenpairs by an amount that is greater than the maximum
		// degeneracy of the eigenvalues. We know that the maximally symmetric
		// lattices have a 6-fold degeneracy, times 2 for polarization. We round
		// up to 16 to be safe.
		SpMat::Index n_arnoldi = (cols > 16 ? 2*cols+16 : cols + 16);
		if(n_arnoldi > A.rows()){ n_arnoldi = A.rows(); }
		
		IRA::ComplexEigensystem::Params params;
		params.max_iterations = 4*A.rows();
		if((size_t)(4*cols) < params.max_iterations){
			params.max_iterations = 4*cols;
		}
		params.tol = 2. * std::numeric_limits<double>::epsilon();
		
		HermitianEigensolver solver(A, cols, n_arnoldi, params);
		SpMat::Index nconv = solver.GetConvergedCount();
		if(nconv >= cols){
			vals = Eigen::Map<const CVec>(solver.GetEigenvalues(), cols);
			vecs = Eigen::Map<const CMat>(solver.GetEigenvectors(), A.rows(), cols);
			SortEigenpairs(vals, vecs);
			ret = 0;
		}else{
			vals.resize(cols); vals.setZero();
			vecs.resize(A.rows(), cols); vecs.setZero();
			ret = (int)nconv - (int)cols;
		}
	}
	return ret;
}

bool LeftIsBetter(
	const doublecomplex &left, const doublecomplex &right
){
	return left.real() > right.real();
}

void SortEigenpairs(CVec &vals, CMat &vecs){
	const size_t n = vals.size();
	// Use selection sort (cycle sort may be better)
	for(size_t j = 0; j+1 < n; ++j){
		size_t imax = j;
		for(size_t i = j+1; i < n; ++i){
			if(LeftIsBetter(vals[i], vals[imax])){
				imax = i;
			}
		}
		if(imax != j){
			std::swap(vals[j], vals[imax]);
			vecs.col(j).swap(vecs.col(imax));
		}
	}
}

} // namespace S4r
