#include "IRA.hpp"
#include <limits>
#include <algorithm>
#include <cmath>

#include <iostream>
#include <cassert>

// References:
//  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
//     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
//     pp 357-385.
//  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
//     Restarted Arnoldi Iteration", Rice University Technical Report
//     TR95-13, Department of Computational and Applied Mathematics.
//  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
//     real_type precision Matrices", Linear Algebra and its Applications,
//     vol 88/89, pp 575-595, (1987).
//  4. B. Nour-Omid, B. N. Parlett, T. Ericsson and P. S. Jensen,
//     "How to Implement the Spectral Transformation", Math Comp.,
//     Vol. 48, No. 178, April, 1987 pp. 664-673.
//
// Authors
//     Danny Sorensen               Phuong Vu
//     Richard Lehoucq              CRPC / Rice University
//     Chao Yang                    Houston, Texas
//     Dept. of Computational &
//     Applied Mathematics
//     Rice University
//     Houston, Texas
//
//     Victor Liu
//     Stanford University
//     Stanford, CA

namespace IRA{

static inline real_type _abs1(const complex_type &z){
	return fabs(z.real()) + fabs(z.imag());
}

static real_type frand(){
	static const real_type nrm(real_type(1) / real_type(RAND_MAX));
	real_type a = rand();
	real_type b = rand();
	return ((a * nrm) + b) * real_type(2) - real_type(1);
}

static real_type HessenbergNorm1(size_t n, const complex_type *a, size_t lda){
	real_type result = 0;
	for(size_t j = 0; j < n; ++j){
		size_t ilimit = (n < j+2 ? n : j+2);
		real_type sum = 0;
		for(size_t i = 0; i < ilimit; ++i){
			sum += std::abs(a[i+j*lda]);
		}
		if(sum > result){ result = sum; }
	}
	return result;
}

template <typename TS, typename T>
static void RescaleMatrix(const TS &cfrom, const TS &cto, size_t m, size_t n, T *a, size_t lda){
	if(n == 0 || m == 0){ return; }

	const TS smlnum = std::numeric_limits<TS>::min();
	const TS bignum = TS(1) / smlnum;

	TS cfromc = cfrom;
	TS ctoc = cto;

	bool done = true;
	do{
		const TS cfrom1 = cfromc * smlnum;
		TS mul;
		if(cfrom1 == cfromc){
			// CFROMC is an inf.  Multiply by a correctly signed zero for
			// finite CTOC, or a NaN if CTOC is infinite.
			mul = ctoc / cfromc;
			done = true;
			//cto1 = ctoc;
		}else{
			const TS cto1 = ctoc / bignum;
			if(cto1 == ctoc){
				// CTOC is either 0 or an inf.  In both cases, CTOC itself
				// serves as the correct multiplication factor.
				mul = ctoc;
				done = true;
				//cfromc = 1.;
			}else if(std::abs(cfrom1) > std::abs(ctoc) && ctoc != TS(0)){
				mul = smlnum;
				done = false;
				cfromc = cfrom1;
			}else if(std::abs(cto1) > std::abs(cfromc)){
				mul = bignum;
				done = false;
				ctoc = cto1;
			}else{
				mul = ctoc / cfromc;
				done = true;
			}
		}

		for(size_t j = 0; j < n; ++j){
			for(size_t i = 0; i < m; ++i){
				a[i+j*lda] *= mul;
			}
		}
	}while(!done);
}

#ifdef IRA_STANDALONE
# include "IRA_deps.cpp"
#else
// External Lapack routines
namespace Lapack{

typedef int integer;
typedef int logical;

extern "C" void zlartg_(
	const complex_type &f, const complex_type &g,
	real_type *cs, complex_type *sn, complex_type *r
);

void RotationGenerate(
	const complex_type &f, const complex_type &g,
	real_type *cs, complex_type *sn, complex_type *r
){
	zlartg_(f, g, cs, sn, r);
}


extern "C" void ztrexc_(
	const char *compq, const integer &n,
	complex_type *t, const integer &ldt,
	complex_type *q, const integer &ldq,
	const integer &ifst, const integer &ilst,
	integer *info
);

// Reorders the Schur factorization of a complex matrix
// A = Q*T*Q**H, so that the diagonal element of T with row index ifst
// is moved to row ilst.
void SchurReorder(
	const char *compq, size_t n, complex_type *t, size_t ldt,
	complex_type *q, size_t ldq, size_t ifst, size_t ilst
){
	integer info;
	ztrexc_(compq, n, t, ldt, q, ldq, ifst+1, ilst+1, &info);
}

extern "C" void zlahqr_(
	const logical &wantt, const logical &wantz,
	const integer &n, const integer &ilo, const integer &ihi,
	complex_type *h, const integer &ldh,
	complex_type *w,
	const integer &iloz, const integer &ihiz,
	complex_type *z, const integer &ldz, integer *info
);

// We assume that H is already upper triangular in columns indices
// ihi to n-1, inclusive, and that h[ilo+(ilo-1)*ldh] == 0.
// iloz and ihiz specify the rows of Z to which transformations
// should be applied if wantz == true. iloz is the lower index (inclusive),
// and ihi is the upper index (exclusive).
int HessenbergQRIteration(
	bool wantt, bool wantz, size_t n, size_t ilo, size_t ihi,
	complex_type *h, int ldh,
	complex_type *w, size_t iloz, size_t ihiz,
	complex_type *z, size_t ldz
){
	integer info;
	zlahqr_(
		wantt ? 1 : 0, wantz ? 1 : 0, n, ilo+1, ihi,
		h, ldh, w, iloz+1, ihiz, z, ldz, &info
	);
	return info;
}

extern "C" void ztrevc_(
	const char *side, const char *howmny, const logical *select,
	const integer &n, complex_type *t, const integer &ldt,
	complex_type *vl, const integer &ldvl,
	complex_type *vr, const integer &ldvr,
	const integer &mm, integer *m, complex_type *work,
	real_type *rwork, integer *info
);

int TriangularEigenvectors(
	const char *howmny, const bool *select,
	size_t n, complex_type *t, size_t ldt,
	complex_type *vr, size_t ldvr,
	size_t mm, complex_type *work, real_type *rwork
){
	integer im, info;
	logical *bsel = NULL;
	if(NULL != select){
		bsel = new logical[n];
		for(size_t i = 0; i < n; ++i){
			bsel[i] = (select[i] ? 1 : 0);
		}
	}
	ztrevc_(
		"R", howmny, bsel, n, t, ldt, NULL, n, vr, ldvr,
		mm, &im, work, rwork, &info
	);
	if(NULL != select){
		delete [] bsel;
	}
	return info;
}

extern "C" void zgeqr2_(
	const integer &m, const integer &n,
	complex_type *a, const integer &lda, complex_type *tau,
	complex_type *work, integer *info
);

void QRFactorization(
	size_t m, size_t n, complex_type *a, size_t lda,
	complex_type *tau, complex_type *work
){
	integer info;
	zgeqr2_(m, n, a, lda, tau, work, &info);
}

extern "C" void zunm2r_(
	const char *side, const char *trans,
	const integer &m, const integer &n, const integer &k,
	complex_type *a, const integer &lda, const complex_type *tau,
	complex_type *c, const integer &ldc, complex_type *work, integer *info
);

void QRMultQ_RN(
	size_t m, size_t n, size_t k, complex_type *a, size_t lda,
	const complex_type *tau, complex_type *c, size_t ldc,
	complex_type *work
){
	integer info;
	zunm2r_("R", "N", m, n, k, a, lda, tau, c, ldc, work, &info);
}

} // namespace Lapack

namespace BLAS{

typedef int integer;

extern "C" void zgemv_(
	const char *trans,
	const integer &m, const integer &n,
	const complex_type &alpha, const complex_type *a, const integer &lda,
	const complex_type *x, const integer &incx,
	const complex_type &beta, complex_type *y, const integer &incy
);

void MultMV(const char *trans, size_t m, size_t n,
	const complex_type &alpha, const complex_type *a, size_t lda,
	const complex_type *x, size_t incx,
	const complex_type &beta, complex_type *y, size_t incy
){
	zgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

extern "C" void ztrmm_(
	const char *side, const char *uplo, const char *transa, const char *diag,
	const integer &m, const integer &n,
	const complex_type &alpha, const complex_type *a, const integer &lda,
	complex_type *b, const integer &ldb
);

void MultTrM_RUNN1(
	size_t m, size_t n,
	const complex_type *a, size_t lda,
	complex_type *b, size_t ldb
){
	ztrmm_("R","U","N","N", m, n, complex_type(1), a, lda, b, ldb);
}

} // namespace BLAS

#endif // IRA_STANDALONE

class ComplexEigensystemImpl{
	friend class ComplexEigensystem;
protected:
	ComplexEigensystem *parent;
	
	size_t nn, n_wanted, n_arnoldi;
	complex_type *workd;           // length 2*n
	complex_type *workl;           // length (3*n_arnoldi + 5)*n_arnoldi.
	complex_type *workev;          // length 2*n_arnoldi;
	complex_type *resid;           // length n
	complex_type *v;               // length n*n_arnoldi
	complex_type *w;               // length n_arnoldi
	real_type *rwork;   // length n_arnoldi
	bool *bwork;        // length n_arnoldi
	
	// Pointers into workspaces
	size_t ldh, ldq;
	complex_type *h;
	complex_type *ritz, *bounds;
	complex_type *q;
	complex_type *subwork;
	
	bool shift_invert;
	complex_type shift;
	size_t maxiter;
	real_type tol;
	bool bsubspace;
	bool resid0_provided;
	
	bool solved;
	
	size_t nconv; // number of converged pairs
	real_type rnorm;
	
	static const size_t nb = 1; // block size
public:
	ComplexEigensystemImpl(ComplexEigensystem *sup,
		size_t n, size_t nwanted,
		size_t narnoldi, // n_wanted+1 <= n_arnoldi <= n
		bool shift_invert_mode, const complex_type &shift_value,
		const ComplexEigensystem::Params &params,
		complex_type *resid0
	):parent(sup),
		nn(n), n_wanted(nwanted), n_arnoldi(narnoldi),
		shift_invert(shift_invert_mode), shift(shift_value),
		maxiter(params.max_iterations), tol(params.tol),
		bsubspace(params.invariant_subspace),
		solved(false), nconv(0)
	{
		if(n_arnoldi > n){ n_arnoldi = n; }
		const size_t ncv = n_arnoldi;
		ldh = ncv;
		ldq = ncv;
		
		w      = new complex_type[n_arnoldi];
		v      = new complex_type[n*n_arnoldi];
		workd  = new complex_type[2*n];
		workl  = new complex_type[(3*n_arnoldi+5)*n_arnoldi];
		workev = new complex_type[2*n_arnoldi];
		rwork  = new real_type[n_arnoldi];
		bwork  = new bool[n_arnoldi];
		if(NULL != resid0){
			resid0_provided = true;
			resid = resid0;
		}else{
			resid0_provided = false;
			resid  = new complex_type[n];
		}
		if(tol <= real_type(0)){
			tol = std::numeric_limits<real_type>::epsilon();
		}
		
		h = workl;
		ritz = h + ldh*ncv;
		bounds = ritz + ncv;
		q = bounds + ncv;
		subwork = q + ldq*ncv;
	}
	~ComplexEigensystemImpl(){
		delete [] v;
		delete [] w;
		delete [] bwork;
		delete [] rwork;
		delete [] workev;
		delete [] workl;
		delete [] workd;
		if(!resid0_provided){
			delete [] resid;
		}
	}
	
protected:
	int RunArnoldi(){ // the old znaupd
		const size_t ncv = n_arnoldi;
		const size_t nev = n_wanted;
		
		memset(workl, 0, sizeof(complex_type) * (ncv * ncv * 3 + ncv * 5));

		nconv = nev;
		int info = Phase1(
			ncv - nev, &maxiter, subwork
		);
		info = ComputeEigenvectors(
			!bsubspace, 'A', bwork, workev, n_wanted, n_arnoldi, workd, workl
		);
		return info;
	}
	
	int Phase1( // the old znaup2
		size_t np, size_t *mxiter, complex_type *workl
	){
		const size_t n = nn;
		int info;
		static const real_type eps23 = pow(std::numeric_limits<real_type>::epsilon(), (real_type)2/(real_type)3);

		size_t nev = n_wanted;
		// Save initial values
		const size_t nev0 = n_wanted;
		const size_t np0 = np;

		// kplusp is the bound on the largest Lanczos factorization built.
		// nconv is the current number of "converged" eigenvalues.
		const size_t kplusp = nev + np;
		nconv = 0;
		
		// Get a possibly random starting vector and
		// force it into the range of the operator OP.
		info = GetInitialVector(resid0_provided, 0, &rnorm);
		if(real_type(0) == rnorm){ // The initial vector is zero. Error exit.
			*mxiter = 0;
			return -9;
		}

		// Compute the first nev steps of the Arnoldi factorization
		info = Iterate(0, nev, &rnorm, workd);
		if(info > 0){
			nconv = info;
			*mxiter = 0;
			return -9999;
		}

		// MAIN ARNOLDI ITERATION LOOP
		// Each iteration implicitly restarts the Arnoldi factorization in place.
		size_t iter = 0;
		do{
			++iter;

			// Compute NP additional steps of the Arnoldi factorization.
			// Adjust NP since NEV might have been updated by last call
			// to the shift application routine znapps .
			np = kplusp - nev;

			// Compute np additional steps of the Arnoldi factorization.
			info = Iterate(nev, np, &rnorm, workd);
			if(info > 0){
				nconv = info;
				*mxiter = iter;
				return -9999;
			}

			// Compute the eigenvalues and corresponding error bounds
			// of the current upper Hessenberg matrix.
			if(0 != ComputeHessenbergEigenvalues(rnorm, kplusp, workl, rwork)){
				nconv = np;
				return -8;
			}
			// Select the wanted Ritz values and their bounds
			// to be used in the convergence test.
			// The wanted part of the spectrum and corresponding
			// error bounds are in the last NEV loc. of RITZ,
			// and BOUNDS respectively.
			nev = nev0;
			np = np0;

			// Make a copy of Ritz values and the corresponding
			// Ritz estimates obtained from zneigh .
			BLAS::Copy(kplusp, ritz  , 1, &workl[kplusp * kplusp], 1);
			BLAS::Copy(kplusp, bounds, 1, &workl[kplusp * kplusp + kplusp], 1);

			// Select the wanted Ritz values and their bounds
			// to be used in the convergence test.
			// The wanted part of the spectrum and corresponding
			// bounds are in the last NEV loc. of RITZ
			// BOUNDS respectively.
			SortRitz(!parent->CanSupplyShifts(), nev, np, ritz, bounds);

			// Convergence test: currently we use the following criteria.
			// The relative accuracy of a Ritz value is considered
			// acceptable if:
			//
			// error_bounds(i) < tol*max(eps23, magnitude_of_ritz(i)).
			nconv = 0;
			for(size_t i = 0; i < nev; ++i){
				real_type rtemp = std::abs(ritz[np + i]);
				if(eps23 > rtemp){ rtemp = eps23; }
				if(std::abs(bounds[np + i]) <= tol * rtemp){
					++nconv;
				}else{
				}
			}

			// Count the number of unwanted Ritz values that have zero
			// Ritz estimates. If any Ritz estimates are equal to zero
			// then a leading block of H of order equal to at least
			// the number of Ritz values with zero Ritz estimates has
			// split off. None of these Ritz values may be removed by
			// shifting. Decrease NP the number of shifts to apply. If
			// no shifts may be applied, then prepare to exit
			{ size_t nptemp = np;
				for(size_t j = 0; j < nptemp; ++j){
					if(real_type(0) == bounds[j]){
						--np;
						++nev;
					}
				}
			}

			if(nconv >= nev0 || iter > *mxiter || np == 0){
				// Prepare to exit. Put the converged Ritz values
				// and corresponding bounds in RITZ(1:NCONV) and
				// BOUNDS(1:NCONV) respectively. Then sort.
				// Be careful when NCONV > NP

				//  Use h( 3,1 ) as storage to communicate
				//  rnorm to zneupd if needed
				h[2] = rnorm;

				// Sort Ritz values so that converged Ritz
				// values appear within the first NEV locations
				// of ritz and bounds, and the most desired one
				// appears at the front.
				SortUtil(kplusp, ritz, bounds, true, false);

				// Scale the Ritz estimate of each Ritz value
				// by 1 / max(eps23, magnitude of the Ritz value).
				for(size_t j = 0; j < nev0; ++j){
					real_type rtemp = std::abs(ritz[j]);
					if(eps23 > rtemp){ rtemp = eps23; }
					bounds[j] /= rtemp;
				}

				// Sort the Ritz values according to the scaled Ritz
				// estimates.  This will push all the converged ones
				// towards the front of ritz, bounds (in the case
				// when NCONV < NEV.)
				SortUtil(nev0, bounds, ritz, true, true);

				// Scale the Ritz estimate back to its original
				// value.
				for(size_t j = 0; j < nev0; ++j){
					real_type rtemp = std::abs(ritz[j]);
					if(eps23 > rtemp){ rtemp = eps23; }
					bounds[j] *= rtemp;
				}

				// Sort the converged Ritz values again so that
				// the "threshold" value appears at the front of
				// ritz and bound.
				SortUtil(nconv, ritz, bounds, false, false);

				// Max iterations have been exceeded.
				if(iter > *mxiter && nconv < nev0){
					info = 1;
				}

				// No shifts to apply.
				if(0 == np && nconv < nev0){
					info = 2;
				}

				np = nconv;
				break; // Exit main loop
			}else if(nconv < nev0 && !parent->CanSupplyShifts()){
				// Do not have all the requested eigenvalues yet.
				// To prevent possible stagnation, adjust the size
				// of NEV.
				size_t nevbef = nev;
				nev += (nconv < np/2 ? nconv : np/2);
				if(nev == 1 && kplusp >= 6){
					nev = kplusp / 2;
				} else if(nev == 1 && kplusp > 3){
					nev = 2;
				}
				np = kplusp - nev;

				// If the size of NEV was just increased,
				// resort the eigenvalues.
				if(nevbef < nev){
					SortRitz(!parent->CanSupplyShifts(), nev, np, ritz, bounds);
				}
			}

			if(parent->CanSupplyShifts()){
				// User specified shifts: pop back out to get the shifts
				// and return them in the first 2*NP locations of WORKL.
				parent->GetShifts(np, ritz, subwork);
				// Move the NP shifts from WORKL to RITZ, to free up WORKL
				BLAS::Copy(np, subwork, 1, ritz, 1);
			}

			// Apply the NP implicit shifts by QR bulge chasing.
			// Each shift is applied to the whole upper Hessenberg matrix H.
			// The first 2*N locations of WORKD are used as workspace.
			ApplyShifts(&nev, np, ritz, workl);

			// Compute the B-norm of the updated residual.
			// Keep B*RESID in WORKD(1:N) to be used in
			// the first step of the next call to znaitr .
			// WORKD(1:N) := B*RESID
			if(!parent->IsBIdentity()){
				if(parent->IsBInPlace()){
					BLAS::Copy(n, resid, 1, workd, 1);
					parent->ApplyB(n, NULL, workd);
				}else{
					parent->ApplyB(n, resid, workd);
				}
				complex_type cmpnorm = BLAS::ConjugateDot(n, resid, 1, workd, 1);
				rnorm = sqrt(std::abs(cmpnorm));
			}else{
				BLAS::Copy(n, resid, 1, workd, 1);
				rnorm = BLAS::Norm2(n, resid, 1);
			}
		}while(1);

		*mxiter = iter;
		return info;
	}
	// (the old zsortc)
	// Sort the arrays x and y of length n according to the comparison function.
	// The comparison function indicates when two elements needs to be swapped.
	// If reverse is true, then the sort order is reversed.
	// If sm is true, then the sort order is to place smallest magnitude values
	// at the beginning of the arrays and the comparison function is not called.
	void SortUtil(size_t n, complex_type *x, complex_type *y, bool reverse, bool sm){
		size_t igap = n / 2;

		while(igap != 0){
			for(size_t i = igap; i < n; ++i){
				int j = i - igap;
				while(j >= 0){
					if(j < 0){ break; }
					bool cmpres;
					if(sm){
						cmpres = (std::abs(x[j]) < std::abs(x[j+igap]));
					}else{
						cmpres = parent->EigenvalueCompare(x[j], x[j+igap]);
					}
					if(reverse == cmpres){
						std::swap(x[j], x[j+igap]);
						if(NULL != y){
							std::swap(y[j], y[j+igap]);
						}
					}else{
						break;
					}
					j -= igap;
				}
			}
			igap /= 2;
		}
	}
	// (the old zngets)
	// Given the eigenvalues of the upper Hessenberg matrix H,
	// computes the NP shifts AMU that are zeros of the polynomial of
	// degree NP which filters out components of the unwanted eigenvectors
	// corresponding to the AMU's based on some given criteria.
	//
	// NOTE: call this even in the case of user specified shifts in order
	// to sort the eigenvalues, and error bounds of H for later use.
	void SortRitz(bool srtunwanted, size_t kev, size_t np, complex_type *r, complex_type *b){
		SortUtil(kev + np, r, b, false, false);

		if(srtunwanted){
			// Sort the unwanted Ritz values used as shifts so that
			// the ones with largest Ritz estimates are first
			// This will tend to minimize the effects of the
			// forward instability of the iteration when the shifts
			// are applied in subroutine znapps.
			// Be careful and use 'SM' since we want to sort BOUNDS!
			SortUtil(np, b, r, false, true);
		}
	}

	// (the old zneigh)
	// Compute the eigenvalues of the current upper Hessenberg matrix
	// and the corresponding Ritz estimates given the current residual norm.
	int ComputeHessenbergEigenvalues(
		const real_type &rnorm, size_t n, complex_type *workl, real_type *rwork
	){
		int ierr;

		// 1. Compute the eigenvalues, the last components of the
		//    corresponding Schur vectors and the full Schur form T
		//    of the current upper Hessenberg matrix H.
		//    zlahqr returns the full Schur form of H
		//    in WORKL(1:N**2), and the Schur vectors in q.
		BLAS::Copy(n, n, h, ldh, workl, n);
		BLAS::Set(n, n, 0.0, 1.0, q, ldq);
		ierr = Lapack::HessenbergQRIteration(true, true, n, 0, n, workl, ldh, ritz, 0, n, q, ldq);
		if(ierr != 0){
			return ierr;
		}
		BLAS::Copy(n, &q[n-2+0*ldq], ldq, bounds, 1);

		// 2. Compute the eigenvectors of the full Schur form T and
		//    apply the Schur vectors to get the corresponding
		//    eigenvectors.
		ierr = Lapack::TriangularEigenvectors("B", NULL, n, workl, n, q, ldq, n, &workl[n*n], rwork);

		if(0 == ierr){
			// Scale the returning eigenvectors so that their
			// Euclidean norms are all one. LAPACK subroutine
			// ztrevc returns each eigenvector normalized so
			// that the element of largest magnitude has
			// magnitude 1; here the magnitude of a complex
			// number (x,y) is taken to be |x| + |y|.
			for(size_t j = 0; j < n; ++j){
				real_type temp = BLAS::Norm2(n, &q[0+j*ldq], 1);
				BLAS::Scale(n, real_type(1) / temp, &q[0+j*ldq], 1);
			}

			// Compute the Ritz estimates
			BLAS::Copy(n, &q[n-1+0*ldq], n, bounds, 1);
			BLAS::Scale(n, rnorm, bounds, 1);
		}
		return ierr;
	}

	// (the old znapps)
	// Given the Arnoldi factorization
	//    A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
	// apply NP implicit shifts resulting in
	//    A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
	// where Q is an orthogonal matrix which is the product of rotations
	// and reflections resulting from the NP bulge change sweeps.
	// The updated Arnoldi factorization becomes:
	//    A*Vnew_{k} - Vnew_{k}*Hnew_{k} = rnew_{k}*e_{k}^T.
	void ApplyShifts(size_t *kev, size_t np, complex_type *shift, complex_type *workl){
		const size_t n = nn;

		real_type c;
		complex_type f, g;
		complex_type r, s, t;
		complex_type h11, h21;
		size_t iend;
		complex_type sigma;
		size_t istart, kplusp;

		// Set machine-dependent constants for the
		// stopping criterion. If norm(H) <= sqrt(OVFL),
		// overflow should not occur.
		// REFERENCE: LAPACK subroutine zlahqr
		const real_type unfl = std::numeric_limits<real_type>::min();
		//const real_type ovfl = 1. / unfl;
		const real_type ulp = 2.0*std::numeric_limits<real_type>::epsilon();
		const real_type smlnum = unfl * (n / ulp);

		kplusp = *kev + np;

		// Initialize Q to the identity to accumulate
		// the rotations and reflections
		BLAS::Set(kplusp, kplusp, 0.0, 1.0, q, ldq);

		// Quick return if there are no shifts to apply
		if(np == 0){
			return;
		}

		// Chase the bulge with the application of each
		// implicit shift. Each shift is applied to the
		// whole matrix including each block.
		for(size_t jj = 0; jj < np; ++jj){
			sigma = shift[jj];

			istart = 0;
			do{
				iend = kplusp-1;
				for(size_t i = istart; i < iend; ++i){
					// Check for splitting and deflation. Use
					// a standard test as in the QR algorithm
					// REFERENCE: LAPACK subroutine zlahqr
					real_type tst1 = _abs1(h[i+i*ldh]) + _abs1(h[(i+1)+(i+1)*ldh]);
					if(tst1 == 0.){
						tst1 = HessenbergNorm1(kplusp - jj, h, ldh);
					}
					if(std::abs(h[(i+1)+i*ldh].real()) <= (ulp*tst1 > smlnum ? ulp*tst1 : smlnum)){
						iend = i;
						h[(i+1)+i*ldh] = 0;
						break;
					}
				}

				// No reason to apply a shift to block of order 1
				// or if the current block starts after the point
				// of compression since we'll discard this stuff
				if(!(istart == iend || istart >= *kev)){
						
					h11 = h[istart+istart*ldh];
					h21 = h[(istart+1)+istart*ldh];
					f = h11 - sigma;
					g = h21;

					for(size_t i = istart; i < iend; ++i){
						// Construct the plane rotation G to zero out the bulge
						Lapack::RotationGenerate(f, g, &c, &s, &r);
						if(i > istart){
							h[i+(i-1)*ldh] = r;
							h[(i+1)+(i-1)*ldh] = 0;
						}

						// Apply rotation to the left of H;  H <- G'*H
						for(size_t j = i; j < kplusp; ++j){
							t = c*h[i+j*ldh] + s*h[(i+1)+j*ldh];
							h[(i+1)+j*ldh] = -std::conj(s)*h[i+j*ldh] + c * h[(i+1)+j*ldh];
							h[i+j*ldh] = t;
						}

						// Apply rotation to the right of H;  H <- H*G
						{ size_t jend = (i+2 < iend ? i+2 : iend);
						for(size_t j = 0; j <= jend; ++j){
							t = c*h[j+i*ldh] + std::conj(s)*h[j+(i+1)*ldh];
							h[j+(i+1)*ldh] = -s*h[j+i*ldh] + c*h[j+(i+1)*ldh];
							h[j+i*ldh] = t;
						} }

						// Accumulate the rotation in the matrix Q;  Q <- Q*G'
						{ size_t jend = (i+jj+2 < kplusp ? i+jj+2 : kplusp);
						for(size_t j = 0; j < jend; ++j){
							t = c*q[j+i*ldq] + std::conj(s)*q[j+(i+1)*ldq];
							q[j+(i+1)*ldq] = -s*q[j+i*ldq] + c*q[j+(i+1)*ldq];
							q[j+i*ldq] = t;
						} }

						// Prepare for next rotation
						if(i+1 < iend){
							f = h[(i+1)+i*ldh];
							g = h[(i+2)+i*ldh];
						}
					}
				} // Finished applying the shift.

				// Apply the same shift to the next block if there is any.
				istart = iend + 1;
				// Loop back to the top to get the next shift.
			}while(iend+1 < kplusp);
		}

		// Perform a similarity transformation that makes
		// sure that the compressed H will have non-negative
		// real subdiagonal elements.
		for(size_t j = 0; j < *kev; ++j){
			if(h[(j+1)+j*ldh].real() < 0. || h[(j+1)+j*ldh].imag() != 0.){
				t = h[(j+1)+j*ldh] / std::abs(h[(j+1)+j*ldh]);
				BLAS::Scale(kplusp - j, std::conj(t), &h[(j+1)+j*ldh], ldh);
				/* Computing MIN */
				size_t sublen = (j+3 < kplusp ? j+3 : kplusp);
				BLAS::Scale(sublen, t, &h[0+(j+1)*ldh], 1);
				sublen = (j+1+np+1 < kplusp ? j+1+np+1 : kplusp);
				BLAS::Scale(sublen, t, &q[0+(j+1)*ldq], 1);
				h[(j+1)+j*ldh] = h[(j+1)+j*ldh].real();
			}
		}

		for(size_t i = 0; i < *kev; ++i){
			// Final check for splitting and deflation.
			// Use a standard test as in the QR algorithm
			// REFERENCE: LAPACK subroutine zlahqr.
			// Note: Since the subdiagonals of the
			// compressed H are nonnegative real numbers,
			// we take advantage of this.
			real_type tst1 = _abs1(h[i+i*ldh]) + _abs1(h[(i+1)+(i+1)*ldh]);
			if(tst1 == 0.){
				tst1 = HessenbergNorm1(*kev, h, ldh);
			}
			if(h[(i+1)+i*ldh].real() <= (ulp*tst1 > smlnum ? ulp*tst1 : smlnum)){
				h[(i+1)+i*ldh] = 0;
			}
		}

		// Compute the (kev+1)-st column of (V*Q) and
		// temporarily store the result in WORKD(N+1:2*N).
		// This is needed in the residual update since we
		// cannot GUARANTEE that the corresponding entry
		// of H would be zero as in exact arithmetic.
		if(h[*kev+(*kev-1)*ldh].real() > 0.){
			BLAS::MultMV("N", n, kplusp, 1.0, v, n, &q[0+(*kev)*ldq], 1, 0.0, &workd[n], 1);
		}

		// Compute column 1 to kev of (V*Q) in backward order
		// taking advantage of the upper Hessenberg structure of Q.
		for(size_t i = 0; i < *kev; ++i){
			BLAS::MultMV("N", n, kplusp - i, 1.0, v, n, &q[0+(*kev - i-1)*ldq], 1, 0.0, workd, 1);
			BLAS::Copy(n, workd, 1, &v[0+(kplusp - i-1)*n], 1);
		}

		//  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev).
		for(size_t i = 0; i < *kev; ++i){
			BLAS::Copy(n, &v[0+(kplusp - *kev+i)*n], 1, &v[i*n], 1);
		}
		// Can't use this due to potential column overlap
		//BLAS::Copy(n, *kev, &v[0+(kplusp - *kev)*n], n, v, n);

		// Copy the (kev+1)-st column of (V*Q) in the appropriate place
		if(h[*kev+(*kev-1)*ldh].real() > 0.){
			BLAS::Copy(n, &workd[n], 1, &v[0+(*kev)*n], 1);
		}

		// Update the residual vector:
		//    r <- sigmak*r + betak*v(:,kev+1)
		// where
		//    sigmak = (e_{kev+p}'*Q)*e_{kev}
		//    betak = e_{kev+1}'*H*e_{kev}
		BLAS::Scale(n, q[(kplusp-1)+(*kev-1)*ldq], resid, 1);
		if(h[*kev+(*kev-1)*ldh].real() > 0.){
			for(size_t i = 0; i < n; ++i){
				resid[i] += h[*kev+(*kev-1)*ldh] * v[i+(*kev)*n];
			}
		}
	}

	// the old zgetv0
	int GetInitialVector(bool initv, size_t j, real_type *rnorm){
		const size_t n = nn;
		// State data
		size_t iter;
		real_type rnorm0;
		
		iter = 0;
		
		// Possibly generate a random starting vector in RESID
		// Use a LAPACK random number generator used by the
		// matrix generation routines.
		//    idist = 1: uniform (0,1)  distribution;
		//    idist = 2: uniform (-1,1) distribution;
		//    idist = 3: normal  (0,1)  distribution;
		if(!initv){
			for(size_t i = 0; i < n; ++i){
				resid[i] = complex_type(frand(), frand());
			}
		}

		if(!parent->IsBIdentity()){
			// Force the starting vector into the range of Op to handle
			// the generalized problem when B is possibly (singular).
			// So we perform:
			//    resid <- Op*B*resid (shift-invert)
			//    resid <- Op*resid   (non-shift-invert)
			if(shift_invert){
				if(parent->IsOpInPlace()){
					if(parent->IsBInPlace()){
						parent->ApplyB(n, NULL, resid);
						parent->ApplyOp(n, NULL, resid);
					}else{
						parent->ApplyB(n, resid, workd);
						BLAS::Copy(n, workd, 1, resid, 1);
						parent->ApplyOp(n, NULL, resid);
					}
				}else{
					if(parent->IsBInPlace()){
						parent->ApplyB(n, NULL, resid);
						parent->ApplyOp(n, resid, workd);
						BLAS::Copy(n, workd, 1, resid, 1);
					}else{
						parent->ApplyB(n, resid, workd);
						parent->ApplyOp(n, workd, resid);
					}
				}
			}else{
				if(parent->IsOpInPlace()){
					parent->ApplyOp(n, NULL, resid);
				}else{
					parent->ApplyOp(n, resid, workd);
					BLAS::Copy(n, workd, 1, resid, 1);
				}
			} // we need to end up with the result in resid
			// Starting vector is now in the range of OP; r = OP*r;
			// Compute B-norm of starting vector.
			parent->ApplyB(n, resid, workd);
			complex_type cnorm = BLAS::ConjugateDot(n, resid, 1, workd, 1);
			rnorm0 = sqrt(std::abs(cnorm));
		}else{
			BLAS::Copy(n, resid, 1, workd, 1);
			rnorm0 = BLAS::Norm2(n, resid, 1);
		}
		*rnorm = rnorm0;
		// At this point, we need resid to be valid, and
		// workd[0..n] contains B*resid

		// Exit if this is the very first Arnoldi step
		if(j > 0){
			// Otherwise need to B-orthogonalize the starting vector against
			// the current Arnoldi basis using Gram-Schmidt with iter. ref.
			// This is the case where an invariant subspace is encountered
			// in the middle of the Arnoldi factorization.
			//
			//       s = V^{T}*B*r;   r = r - V*s;
			//
			// Stopping criteria used for iter. ref. is discussed in
			// Parlett's book, page 107 and in Gragg & Reichel TOMS paper.
			do{
				assert(j <= n);
				BLAS::MultMV("C", n, j, 1.0, v, n, workd, 1, 0.0, &workd[n], 1);
				BLAS::MultMV("N", n, j, -1.0, v, n, &workd[n], 1, 1.0, resid, 1);

				// Compute the B-norm of the orthogonalized starting vector
				if(!parent->IsBIdentity()){
					if(parent->IsBInPlace()){
						BLAS::Copy(n, resid, 1, workd, 1);
						parent->ApplyB(n, NULL, workd);
					}else{
						parent->ApplyB(n, resid, workd);
					}
					complex_type cnorm = BLAS::ConjugateDot(n, resid, 1, workd, 1);
					*rnorm = sqrt(std::abs(cnorm));
				}else{
					BLAS::Copy(n, resid, 1, workd, 1);
					*rnorm = BLAS::Norm2(n, resid, 1);
				}

				// Check for further orthogonalization.
				if(*rnorm > rnorm0 * .717f){
					break;
				}

				++iter;
				if(iter > 1){
					// Iterative refinement step "failed"
					for(size_t jj = 0; jj < n; ++jj){
						resid[jj] = 0;
					}
					*rnorm = 0.;
					return -1;
				}
				// Perform iterative refinement step
				rnorm0 = *rnorm;
			}while(1);
		}
		return 0;
	}

	// (the old znaitr)
	int Iterate(size_t k, size_t np, real_type *rnorm, complex_type *workd){
		const size_t n = nn;
		const size_t ldv = nn;

		using namespace std;

		// Set machine-dependent constants for the
		// the splitting and deflation criterion.
		// If norm(H) <= sqrt(OVFL),
		// overflow should not occur.
		// REFERENCE: LAPACK subroutine zlahqr
		static const real_type unfl = std::numeric_limits<real_type>::min();
		//const real_type ovfl = 1. / unfl;
		static const real_type ulp = 2.*std::numeric_limits<real_type>::epsilon();
		const real_type smlnum = unfl * (n / ulp);

		// Initial call to this routine
		size_t j = k;
		
		complex_type *work_pj = &workd[0];
		complex_type *work_rj = &workd[n];

		// Note:  B*r_{j-1} is already in workd[0..n] = work_pj
		do{
			// STEP 1: Check if the B norm of j-th residual
			// vector is zero. Equivalent to determine whether
			// an exact j-step Arnoldi factorization is present.
			real_type betaj = *rnorm;
			if(!(*rnorm > 0.)){
				// Invariant subspace found, generate a new starting
				// vector which is orthogonal to the current Arnoldi
				// basis and continue the iteration.
				betaj = 0.;
				
				// ITRY is the loop variable that controls the
				// maximum amount of times that a restart is
				// attempted. NRSTRT is used by stat.h
				size_t itry = 0;

				int ierr;
				do{
					ierr = GetInitialVector(false, j, rnorm);
					++itry;
				}while(ierr < 0 && itry <= 3);
				if(ierr < 0){
					// Give up after several restart attempts.
					// Set INFO to the size of the invariant subspace
					// which spans OP and exit.
					return j;
				}
			}

			// STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm
			// Note that p_{j} = B*r_{j-1}. In order to avoid overflow
			// when reciprocating a small RNORM, test against lower
			// machine bound.

			BLAS::Copy(n, resid, 1, &v[0+j*ldv], 1);
			if(*rnorm >= unfl){
				real_type temp1 = 1. / *rnorm;
				BLAS::Scale(n, temp1, &v[0+j*ldv], 1);
				BLAS::Scale(n, temp1, work_pj, 1);
			} else {
				// To scale both v_{j} and p_{j} carefully
				// use LAPACK routine zlascl
				RescaleMatrix(*rnorm, 1.0, n, 1, &v[0+j*ldv], n);
				RescaleMatrix(*rnorm, 1.0, n, 1, work_pj, n);
			}

			// STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j}
			// Note that this is not quite yet r_{j}. See STEP 4
			if(shift_invert){
				// Compute inv(A - s B) * B * x
				// workd[ipj] already contains B*x
				if(parent->IsOpInPlace()){
					BLAS::Copy(n, work_pj, 1, work_rj, 1);
					parent->ApplyOp(n, NULL, work_rj);
				}else{
					parent->ApplyOp(n, work_pj, work_rj);
				}
			}else{
				// Compute (inv(B)*A) * x
				if(parent->IsOpInPlace()){
					BLAS::Copy(n, &v[0+j*ldv], 1, work_rj, 1);
					parent->ApplyOp(n, NULL, work_rj);
				}else{
					parent->ApplyOp(n, &v[0+j*ldv], work_rj);
				}
			}

			// Put another copy of OP*v_{j} into RESID.
			BLAS::Copy(n, work_rj, 1, resid, 1);

			// STEP 4:  Finish extending the Arnoldi
			//          factorization to length j.
			// work_pj contains B*Op*v_{j}
			// The following is needed for STEP 5.
			// Compute the B-norm of OP*v_{j}.
			real_type wnorm;
			if(!parent->IsBIdentity()){
				if(parent->IsBInPlace()){
					BLAS::Copy(n, work_rj, 1, work_pj, 1);
					parent->ApplyB(n, NULL, work_pj);
				}else{
					parent->ApplyB(n, work_rj, work_pj);
				}
				complex_type cnorm = BLAS::ConjugateDot(n, resid, 1, work_pj, 1);
				wnorm = sqrt(std::abs(cnorm));
			}else{
				BLAS::Copy(n, resid, 1, work_pj, 1);
				wnorm = BLAS::Norm2(n, resid, 1);
			}

			// Compute the j-th residual corresponding
			// to the j step factorization.
			// Use Classical Gram Schmidt and compute:
			// w_{j} <-  V_{j}^T * B * OP * v_{j}
			// r_{j} <-  OP*v_{j} - V_{j} * w_{j}

			// Compute the j Fourier coefficients w_{j}
			// work_pj contains B*Op*v_{j}.
			BLAS::MultMV("C", n, j+1, 1.0, v, ldv, work_pj, 1, 0.0, &h[0+j*ldh], 1);

			// Orthogonalize r_{j} against V_{j}.
			// resid contains Op*v_{j}. See STEP 3.
			BLAS::MultMV("N", n, j+1, -1.0, v, ldv, &h[0+j*ldh], 1, 1.0, resid, 1);

			if(j > 0){
				h[j+(j-1)*ldh] = betaj;
			}

			// work_pj contains B*r_{j}.
			// Compute the B-norm of r_{j}.
			if(!parent->IsBIdentity()){
				if(parent->IsBInPlace()){
					BLAS::Copy(n, resid, 1, work_pj, 1);
					parent->ApplyB(n, NULL, work_pj);
				}else{
					parent->ApplyB(n, resid, work_pj);
				}
				complex_type cnorm = BLAS::ConjugateDot(n, resid, 1, work_pj, 1);
				*rnorm = sqrt(std::abs(cnorm));
			}else{
				BLAS::Copy(n, resid, 1, work_pj, 1);
				*rnorm = BLAS::Norm2(n, resid, 1);
			}

			// STEP 5: Re-orthogonalization / Iterative refinement phase
			// Maximum NITER_ITREF tries.
			//
			//          s      = V_{j}^T * B * r_{j}
			//          r_{j}  = r_{j} - V_{j}*s
			//          alphaj = alphaj + s_{j}
			//
			// The stopping criteria used for iterative refinement is
			// discussed in Parlett's book SEP, page 107 and in Gragg &
			// Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.
			// Determine if we need to correct the residual. The goal is
			// to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||
			// The following test determines whether the sine of the
			// angle between  OP*x and the computed residual is less
			// than or equal to 0.717.
			if(*rnorm <= wnorm * .717f){
				size_t iter = 0;

				// Enter the Iterative refinement phase. If further
				// refinement is necessary, loop back here. The loop
				// variable is ITER. Perform a step of Classical
				// Gram-Schmidt using all the Arnoldi vectors V_{j}
				do{
					// Compute V_{j}^T * B * r_{j}.
					// WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1).

					BLAS::MultMV("C", n, j+1, 1.0, v, ldv, work_pj, 1, 0.0, work_rj, 1);

					// Compute the correction to the residual:
					// r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1).
					// The correction to H is v(:,1:J)*H(1:J,1:J)
					// + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.

					BLAS::MultMV("N", n, j+1, -1.0, v, ldv, work_rj, 1, 1.0, resid, 1);
					for(size_t i = 0; i < j; ++i){
						h[i+j*ldh] += work_rj[i];
					}

					// Compute the B-norm of the corrected residual r_{j}.
					real_type rnorm1;
					if(!parent->IsBIdentity()){
						if(parent->IsBInPlace()){
							BLAS::Copy(n, resid, 1, work_pj, 1);
							parent->ApplyB(n, NULL, work_pj);
						}else{
							parent->ApplyB(n, resid, work_pj);
						}
						complex_type cnorm = BLAS::ConjugateDot(n, resid, 1, work_pj, 1);
						rnorm1 = sqrt(std::abs(cnorm));
					}else{
						BLAS::Copy(n, resid, 1, work_pj, 1);
						rnorm1 = BLAS::Norm2(n, resid, 1);
					}

					// Determine if we need to perform another
					// step of re-orthogonalization.
					if(rnorm1 > *rnorm * .717f){
						// No need for further refinement.
						// The cosine of the angle between the
						// corrected residual vector and the old
						// residual vector is greater than 0.717
						// In other words the corrected residual
						// and the old residual vector share an
						// angle of less than arcCOS(0.717)
						*rnorm = rnorm1;
					} else {
						// Another step of iterative refinement step
						// is required. NITREF is used by stat.h
						*rnorm = rnorm1;
						++iter;
						if(iter <= 1){
							continue;
						}

						// Otherwise RESID is numerically in the span of V
						for(size_t jj = 0; jj < n; ++jj){
							resid[jj] = 0;
						}
						*rnorm = 0.;
					}
				}while(0);
			}

			// STEP 6: Update  j = j+1;  Continue
			++j;
			if(j >= k + np){
				size_t istart = (k > 1 ? k : 1);
				for(size_t i = istart; i < k+np; ++i){
					// Check for splitting and deflation.
					// Use a standard test as in the QR algorithm
					// REFERENCE: LAPACK subroutine zlahqr
					real_type tst1 = std::abs(h[(i-1)+(i-1)*ldh]) + std::abs(h[i+i*ldh]);
					if(tst1 == 0.){
						tst1 = HessenbergNorm1(k + np, h, ldh);
					}
					/* Computing MAX */
					if(std::abs(h[i+i*ldh]) <= max(ulp * tst1,smlnum)){
						h[i+(i-1)*ldh] = 0;
					}
				}

				return 0;
			}
		}while(1);

		return 0;
	}

	// (the old zneupd)
	// We assume that rnorm and nconv are set from a previous call to
	// Phase1. Furthermore, we assume that the work arrays have not
	// been modified.
	// 
	//   rvec:   false if only Ritz values are desired, true if Ritz
	//           vectors or Schur vectors are also desired.
	//   howmny: 'A' for NEV Ritz vectors
	//           'P' for NEV Schur vectors
	//           'S' for some Ritz vectors, selected by select
	//   select: bool array of length ncv
	//
	// Ritz values returned in w, vectors in v.
	//
	// Notes:
	// 1. Currently only HOWMNY = 'A' and 'P' are implemented. 
	//
	// 2. Schur vectors are an orthogonal representation for the basis of
	//    Ritz vectors. Thus, their numerical properties are often superior.
	//    If RVEC = .true. then the relationship
	//            A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
	//      transpose( V(:,1:IPARAM(5)) ) * V(:,1:IPARAM(5)) = I
	//    are approximately satisfied.
	//    Here T is the leading submatrix of order IPARAM(5) of the 
	//    upper triangular matrix stored workl(ipntr(12)). 
	//
	int ComputeEigenvectors(
		bool rvec, char howmny, bool *select, complex_type *workev,
		size_t nev, size_t ncv, complex_type *workd, complex_type *workl
	){
		const size_t n = nn;
		const size_t ldv = nn;
		complex_type *z = v;
		const size_t ldz = nn;
		using namespace std;
		complex_type temp;

		bool reord;
		real_type rtemp;
		size_t numcnv;

		// Get machine dependent constant.
		const real_type eps23 = pow(std::numeric_limits<real_type>::epsilon(), (real_type)2/(real_type)3);
		
		//const size_t bounds = ncv*ncv+3*ncv;
		const size_t ldh = ncv;
		const size_t ldq = ncv;

		// Input pointers from Phase 1
		complex_type *work_h = h;
		complex_type *work_ritz = ritz;
		complex_type *work_bounds = bounds;
		complex_type *work_heig = q;
		
		// Output pointers
		complex_type *work_hbds = q + ldh;
		complex_type *work_uptri = work_hbds + ldh;
		complex_type *work_invsub = work_uptri + ldh*ncv;
		
		// rz points to ritz values      computed by zneigh before exiting znaup2.
		// bd points to bounds estimates computed by zneigh before exiting znaup2.
		complex_type *work_rz = &workl[2*ncv*ncv+2*ncv + ncv * ncv];
		complex_type *work_bd = &workl[2*ncv*ncv+2*ncv + ncv * ncv + ncv];
		if(rvec){
			reord = false;

			// Use the temporary bounds array to store indices
			// These will be used to mark the select array later
			for(size_t j = 0; j < ncv; ++j){
				work_bounds[j] = (real_type)j;
				select[j] = false;
			}

			// Select the wanted Ritz values.
			// Sort the Ritz values so that the
			// wanted ones appear at the tailing
			// NEV positions of workl(irr) and
			// workl(iri).  Move the corresponding
			// error estimates in workl(ibd)
			// accordingly.
			SortRitz(false, nev, ncv - nev, work_rz, work_bounds);

			// Record indices of the converged wanted Ritz values
			// Mark the select array for possible reordering
			numcnv = 0;
			for(size_t j = 0; j < ncv; ++j){
				/* Computing MAX */
				rtemp = std::abs(work_rz[ncv - j-1]);
				if(eps23 > rtemp){ rtemp = eps23; }
				size_t jj = (int)work_bounds[ncv - j-1].real();
				if(numcnv < nconv && std::abs(work_bd[jj]) <= tol * rtemp){
					select[jj] = true;
					++numcnv;
					if(jj >= nev){
						reord = true;
					}
				}
			}

			// Check the count (numcnv) of converged Ritz values with
			// the number (nconv) reported by dnaupd.  If these two
			// are different then there has probably been an error
			// caused by incorrect passing of the dnaupd data.
			if(numcnv != nconv){
				return -15;
			}

			// Call LAPACK routine zlahqr to compute the Schur form
			// of the upper Hessenberg matrix returned by ZNAUPD.
			// Make a copy of the upper Hessenberg matrix.
			// Initialize the Schur vector matrix Q to the identity.
			BLAS::Copy(ldh * ncv, work_h, 1, work_uptri, 1);
			BLAS::Set(ncv, ncv, 0.0, 1.0, work_invsub, ldq);
			int ierr = Lapack::HessenbergQRIteration(true, true, ncv, 0, ncv, work_uptri, ldh, work_heig, 0, ncv, work_invsub, ldq);
			BLAS::Copy(ncv, &work_invsub[ncv - 1], ldq, work_hbds, 1);

			if(ierr != 0){
				return -8;
			}

			if(reord){
				// Reorder the computed upper triangular matrix.
				//ztrsen_(&select[1], ncv, &workl[iuptri], ldh, &workl[invsub], ldq, &workl[iheig], &nconv, NULL, NULL);
				{
					nconv = 0;
					for(size_t k = 0; k < ncv; ++k){
						if(select[k]){
							++nconv;
						}
					}

					if(!(nconv == ncv || nconv == 0)){
						// Collect the selected eigenvalues at the top left corner of T.
						size_t ks = 0;
						for(size_t k = 0; k < ncv; ++k){
							if(select[k]){
								// Swap the K-th eigenvalue to position KS.
								if(k != ks){
									Lapack::SchurReorder("V", ncv, work_uptri, ldh, work_invsub, ldq, k, ks);
								}
								++ks;
							}
						}
					}
					// Copy reordered eigenvalues to W.
					for(size_t k = 0; k < ncv; ++k){
						work_heig[k] = work_uptri[k+k*ldh];
					}
				}
				if(ierr == 1){
					return 1;
				}
			}

			// Copy the last row of the Schur basis matrix
			// to workl(ihbds).  This vector will be used
			// to compute the Ritz estimates of converged
			// Ritz values.
			BLAS::Copy(ncv, &work_invsub[ncv-1], ldq, work_hbds, 1);

			// Place the computed eigenvalues of H into D
			// if a spectral transformation was not used.
			if(!shift_invert){
				BLAS::Copy(nconv, work_heig, 1, w, 1);
			}

			// Compute the QR factorization of the matrix representing
			// the wanted invariant subspace located in the first NCONV
			// columns of workl(invsub,ldq).
			Lapack::QRFactorization(ncv, nconv, work_invsub, ldq, workev, &workev[ncv]);

			// * Postmultiply V by Q using zunm2r.
			// * Copy the first NCONV columns of VQ into Z.
			// * Postmultiply Z by R.
			// The N by NCONV matrix Z is now a matrix representation
			// of the approximate invariant subspace associated with
			// the Ritz values in workl(iheig). The first NCONV
			// columns of V are now approximate Schur vectors
			// associated with the upper triangular matrix of order
			// NCONV in workl(iuptri).
			Lapack::QRMultQ_RN(n, ncv, nconv, work_invsub, ldq, workev, v, ldv, &workd[n]);
			//BLAS::Copy(n, nconv, v, ldv, z, ldz); // this is a Nop

			for(size_t j = 0; j < nconv; ++j){
				// Perform both a column and row scaling if the
				// diagonal element of workl(invsub,ldq) is negative
				// I'm lazy and don't take advantage of the upper
				// triangular form of workl(iuptri,ldq).
				// Note that since Q is orthogonal, R is a diagonal
				// matrix consisting of plus or minus ones.
				if(work_invsub[j+j*ldq].real() < 0.){
					BLAS::Scale(nconv, -1.0, &work_uptri[j], ldq);
					BLAS::Scale(nconv, -1.0, &work_uptri[j*ldq], 1);
				}
			}

			if(howmny == 'A'){
				// Compute the NCONV wanted eigenvectors of T
				// located in workl(iuptri,ldq).
				for(size_t j = 0; j < ncv; ++j){
					if(j < nconv){
						select[j] = true;
					} else {
						select[j] = false;
					}
				}

				ierr = Lapack::TriangularEigenvectors("S", select, ncv, work_uptri, ldq, work_invsub, ldq, ncv, workev, rwork);
				if(ierr != 0){
					return -9;
				}

				// Scale the returning eigenvectors so that their
				// Euclidean norms are all one. LAPACK subroutine
				// ztrevc returns each eigenvector normalized so
				// that the element of largest magnitude has
				// magnitude 1.

				for(size_t j = 0; j < nconv; ++j){
					rtemp = BLAS::Norm2(ncv, &work_invsub[j*ldq], 1);
					BLAS::Scale(ncv, 1. / rtemp, &work_invsub[j*ldq], 1);

					// Ritz estimates can be obtained by taking
					// the inner product of the last row of the
					// Schur basis of H with eigenvectors of T.
					// Note that the eigenvector matrix of T is
					// upper triangular, thus the length of the
					// inner product can be set to j.
					workev[j] = BLAS::ConjugateDot(j+1, work_hbds, 1, &work_invsub[j*ldq], 1);
				}


				// Copy Ritz estimates into workl(ihbds)
				BLAS::Copy(nconv, workev, 1, work_hbds, 1);

				// The eigenvector matrix Q of T is triangular.
				// Form Z*Q.
				BLAS::MultTrM_RUNN1(n, nconv, work_invsub, ldq, z, ldz);
			}
		}else{
			// An approximate invariant subspace is not needed.
			// Place the Ritz values computed ZNAUPD into D.
			BLAS::Copy(nconv, work_ritz, 1, w, 1);
			BLAS::Copy(nconv, work_ritz, 1, work_heig, 1);
			BLAS::Copy(nconv, work_bounds, 1, work_hbds, 1);
		}

		// Transform the Ritz values and possibly vectors
		// and corresponding error bounds of OP to those
		// of A*x = lambda*B*x.
		if(!shift_invert){
			if(rvec){
				BLAS::Scale(ncv, rnorm, work_hbds, 1);
			}
		}else{
			// A spectral transformation was used.
			// Determine the Ritz estimates of the
			// Ritz values in the original system.
			if(rvec){
				BLAS::Scale(ncv, rnorm, work_hbds, 1);
			}
			for(size_t k = 0; k < ncv; ++k){
				complex_type temp = work_heig[k];
				work_hbds[k] /= temp;
				work_hbds[k] /= temp;
			}
		}

		// Transform the Ritz values back to the original system.
		// For shift-invert the transformation is
		//          lambda = 1/theta + sigma
		// Note the Ritz vectors are not affected by the transformation.
		if(shift_invert){
			for(size_t k = 0; k < nconv; ++k){
				w[k] = 1.0/work_heig[k] + shift;
			}
		}

		// Eigenvector Purification step. Formally perform
		// one of inverse subspace iteration. Only used
		// for shift-invert. See references.
		if(rvec && howmny == 'A' && shift_invert){
			// Purify the computed Ritz vectors by adding a
			// little bit of the residual vector:
			//                      T
			//          resid(:)*( e    s ) / theta
			//                      NCV
			// where H s = s theta.
			for(size_t j = 0; j < nconv; ++j){
				if(work_heig[j] != 0.){
					workev[j] = work_invsub[j*ldq + ncv-1] / work_heig[j];
				}
			}
			// Perform a rank one update to Z and
			// purify all the Ritz vectors together.
			for(size_t j = 0; j < nconv; ++j){
				if(complex_type(0) != workev[j]){
					for(size_t i = 0; i < n; ++i){
						z[i+j*ldz] += resid[i]*workev[j];
					}
				}
			}
		}

		return 0;
	}
};

ComplexEigensystem::ComplexEigensystem(
	size_t n, size_t n_wanted,
	size_t n_arnoldi, // n_wanted+1 <= n_arnoldi <= n
	bool shift_invert, const complex_type &shift,
	const Params &params,
	complex_type *resid0
){
	impl = new ComplexEigensystemImpl(this,
		n, n_wanted, n_arnoldi, shift_invert, shift, params, resid0
	);
}
ComplexEigensystem::~ComplexEigensystem(){ delete impl; }

size_t ComplexEigensystem::GetConvergedCount(){
	if(!impl->solved){
		impl->RunArnoldi();
		impl->solved = true;
	}
	return impl->nconv;
}
const complex_type *ComplexEigensystem::GetEigenvalues(){
	size_t n = GetConvergedCount();
	return (n > 0 ? impl->w : NULL);
}
const complex_type *ComplexEigensystem::GetEigenvectors(){
	size_t n = GetConvergedCount();
	return (n > 0 ? impl->v : NULL);
}

} // namespace IRA

