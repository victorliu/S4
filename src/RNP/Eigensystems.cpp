#include <cstddef>
#include <cmath>
#include <complex>
#include <limits>
#include "Eigensystems.h"

// A non-parallel implementation of:
//   A parallel algorithm for the eigenvalues
//   and eigenvectors of a general complex matrix
//   by Guatam M. Shroff
//   Numerische Mathematik 58, pages 779-805 (1991)
// The matrix is only reduced to upper triangular form,
// but in practice it is always diagonal.

// Returns the unitary transformation parameters
// Returns the reduction in lower triangular norm
static inline double
pnrj_unitary(size_t n, size_t lda, std::complex<double> *A, size_t p, size_t q, std::complex<double> cs[3]){
	std::complex<double> Apq = A[p+q*lda];
	double aApq2 = std::norm(Apq);
	if(0 != aApq2){
		std::complex<double> Aqp = A[q+p*lda];
		std::complex<double> dpq = 0.5*(A[p+p*lda] - A[q+q*lda]);
		std::complex<double> root = sqrt(dpq*dpq + Apq*Aqp);
		
		std::complex<double> dmax_plus  = dpq + root;
		std::complex<double> dmax_minus = dpq - root;
		std::complex<double> dmax;
		if(std::norm(dmax_plus) > std::norm(dmax_minus)){
			dmax = dmax_plus;
		}else{
			dmax = dmax_minus;
		}
		// At this point we should choose
		// tan x = exp(i theta) Apq/dmax
		// such that theta makes tan x real
		// The resulting parameters are then
		// [ cs[0]  cs[2] ]
		// [ cs[1]  cs[3] ]
		// where
		// cs[0] = cs[3] = cos x, cs[1] = exp(i theta) sin x
		// cs[2] = -exp(-i theta) sin x
		//
		// Thus, tan x = std::abs(Apq/dmax)
		dmax = Apq/dmax; // we use dmax as Aqp/dmax
		double tanx = std::abs(dmax);
		std::complex<double> exp_i_theta = tanx/dmax; // this is actually exp(-i theta)
		if(tanx > 1){ tanx = 1; } // limit the maximum rotation angle
		// Now cos^2 = 1/(1+tan^2), sin^2 = tan^2 cos^2
		cs[0] = double(1)/sqrt(tanx*tanx+1);
		//cs[3] = cs[0];
		cs[1] = tanx*cs[0]; // cs[1] = sin x
		//cs[2] = cs[1];
		cs[1] *= exp_i_theta;
		//cs[2] /= -exp_i_theta;
		cs[2] = -std::conj(cs[1]);
	}else{
		cs[0] = 1; // cs[3] = cs[0];
		cs[1] = cs[2] = 0;
	}

	return aApq2;
}

static inline double
pnrj_shear(size_t n, size_t lda, std::complex<double> *A, size_t p, size_t q, std::complex<double> cs[3]){
	double Gpq = 0;
	std::complex<double> cpq(0);
	for(size_t j = 0; j < n; ++j){
		cpq += (A[p+j*lda]*std::conj(A[q+j*lda]) - std::conj(A[j+p*lda])*A[j+q*lda]);
		
		if(j == p || j == q){ continue; }
		double Gterm = 0;
		Gterm += std::norm(A[p+j*lda]);
		Gterm += std::norm(A[q+j*lda]);
		Gterm += std::norm(A[j+p*lda]);
		Gterm += std::norm(A[j+q*lda]);
		Gpq += Gterm;
	}
	std::complex<double> dpq = A[q+q*lda] - A[p+p*lda];
	// xi_pq = exp(i alpha) Aqp + exp(-i alpha) Apq
	// alpha = arg(cpq) - pi/2
	// Thus, xi_pq = -i exp(i arg(cpq)) Aqp + i exp(-i arg(cpq)) Apq
	// But exp(i arg(cpq)) is simply cpq/|cpq|
	double acpq = std::abs(cpq);
	if(0 == acpq){
		cs[0] = 1;
		cs[1] = cs[2] = 0;
		return 0;
	}
	std::complex<double> eialpha = std::complex<double>(0,-1)*(cpq/acpq);
	std::complex<double> xipq = eialpha*A[q+p*lda] + A[p+q*lda]/eialpha;
	// Now, we will generate the transformation
	// [ cs[0]  cs[2] ]
	// [ cs[1]  cs[3] ]
	// where
	// cs[0] = cs[3] = cosh y,
	// cs[1] =  i exp(-i alpha) sinh y
	// cs[2] = -i exp( i alpha) sinh y
	// and
	// tanh y = -|cpq| / (2*(|dpq|^2 + |xipq|^2) + Gpq)
	double tanhy = -acpq / (2*(std::norm(dpq) + std::norm(xipq)) + Gpq);
	// cosh^2 - sinh^2 = 1, tanh = sinh/cosh
	double coshy = double(1)/sqrt(double(1) - tanhy*tanhy);
	cs[0] = coshy; // cs[3] = cs[0]
	double sinhy = coshy*tanhy;
	cs[1] = std::complex<double>(0, sinhy)/eialpha;
	//cs[2] = std::complex<double>(0,-sinhy)*eialpha;
	cs[2] = std::conj(cs[1]);

	return 0;
}

static inline double
pnrj_diagonal(size_t n, size_t lda, std::complex<double> *A, size_t j, double *t){
	double g = 0;
	double h = 0;
	for(size_t i = 0; i < n; ++i){
		if(i == j){ continue; }
		g += std::norm(A[i+j*lda]);
		h += std::norm(A[j+i*lda]);
	}
	g = sqrt(g); h = sqrt(h);
	*t = sqrt(h/g);
	
	static const double tlimit = 1e8;
	
	if(*t > tlimit){
		*t = tlimit;
		h = tlimit*tlimit*g;
	}else if(*t < double(1)/tlimit){
		*t = double(1)/tlimit;
		g = h/(tlimit*tlimit);
	}
	
	// Compute the estimate of norm reduction
	g -= h;
	return g*g;
}

inline static void
pnrj_apply_rotation_L(size_t n, size_t lda, std::complex<double>* A, size_t p, size_t q, std::complex<double> cs[3]){
	// Apply rotation to matrix A,  A' = J^{-1} A
	for (size_t j = 0; j < n; j++){
		std::complex<double> Apj = A[p+j*lda];
		std::complex<double> Aqj = A[q+j*lda];
		A[p+j*lda] =  Apj * cs[0] - Aqj * cs[2];
		A[q+j*lda] = -Apj * cs[1] + Aqj * cs[0]; // cs[3]
	}
}

inline static void
pnrj_apply_rotation_R(size_t n, size_t lda, std::complex<double>* A, size_t p, size_t q, std::complex<double> cs[3]){
	// Apply rotation to matrix A,  A' = A J
	for(size_t i = 0; i < n; i++){
		std::complex<double> Aip = A[i+p*lda];
		std::complex<double> Aiq = A[i+q*lda];
		A[i+p*lda] = Aip * cs[0] + Aiq * cs[1];
		A[i+q*lda] = Aip * cs[2] + Aiq * cs[0]; // cs[3]
	}
}

inline static void
pnrj_diagonal_L(size_t n, size_t lda, std::complex<double>* A, size_t j, double t){
	t = (double)1/t;
	// Apply diagonal to matrix A,  A' = A J
	for(size_t i = 0; i < n; i++){
		A[j+i*lda] *= t;
	}
}

inline static void
pnrj_diagonal_R(size_t n, size_t lda, std::complex<double>* A, size_t j, double t){
	// Apply diagonal to matrix A,  A' = A J
	for(size_t i = 0; i < n; i++){
		A[i+j*lda] *= t;
	}
}

int RNP::Eigensystem_jacobi(size_t n, 
	std::complex<double> *A, size_t lda,
	std::complex<double> *eval,
	std::complex<double> *vl, size_t ldvl, std::complex<double> *vr, size_t ldvr,
	std::complex<double> *work, double *rwork)
{
	if(NULL != vr){
		for(size_t p = 0; p < n; ++p){
			for(size_t q = 0; q < n; ++q){
				if(p == q){
					vr[q+p*ldvr] = 1;
				}else{
					vr[q+p*ldvr] = 0;
				}
			}
			eval[p] = 0;
		}
	}
	if(NULL != vl){
		for(size_t p = 0; p < n; ++p){
			for(size_t q = 0; q < n; ++q){
				if(p == q){
					vl[q+p*ldvl] = 1;
				}else{
					vl[q+p*ldvl] = 0;
				}
			}
		}
	}
	// This is a heuristic that seems to work well
	size_t max_iter = (int)(2+4.04*log((double)n));
	if(max_iter < 8){ max_iter = 8; }

	double normL = std::numeric_limits<double>::max();
	double normL_prev;
	const double normL_threshold = n*n/2 * std::numeric_limits<double>::epsilon();

	size_t iter = 0;
	bool converged = false;
	do{
		normL_prev = normL;
		// Compute the Frobenius norm of the lower triangle
		normL = 0;
		for(size_t q = 0; q < n-1; ++q){ // column
			for(size_t p = q+1; p < n; ++p){
				double absLpq2 = std::norm(A[p+q*lda]);
				normL += absLpq2;
			}
		}
		normL = sqrt(normL);
//		std::cout << "norm = " << normL << std::endl;
		
		if(normL < normL_threshold){
			converged = true;
			break;
		}
		
		// Perform a sweep
		//   A sweep is a set of rotations followed by a set of diagonal
		//   transformations. A rotation is a shear followed by a unitary
		//   transformation. Rotations are performed on all subdiagonal
		//   elements, while diagonal transformations are applied to each
		//   element of the diagonal.
		
		// Apply all rotations
		for(size_t q = 0; q < n-1; ++q){
			for(size_t p = q+1; p < n; ++p){
				std::complex<double> cs[3];
				
				pnrj_shear(n, lda, A, p, q, cs);
				pnrj_apply_rotation_L(n, lda, A, p, q, cs);
				pnrj_apply_rotation_R(n, lda, A, p, q, cs);
				if(NULL != vr){
					pnrj_apply_rotation_R(n, ldvr, vr, p, q, cs);
				}
				if(NULL != vl){
					cs[1] = -cs[1];
					cs[2] = -cs[2];
					pnrj_apply_rotation_R(n, ldvl, vl, p, q, cs);
				}
				
				pnrj_unitary(n, lda, A, p, q, cs);
				pnrj_apply_rotation_L(n, lda, A, p, q, cs);
				pnrj_apply_rotation_R(n, lda, A, p, q, cs);
				if(NULL != vr){
					pnrj_apply_rotation_R(n, ldvr, vr, p, q, cs);
				}
				if(NULL != vl){
					pnrj_apply_rotation_R(n, ldvl, vl, p, q, cs);
				}
			}
		}
		// Apply all diagonal transformations
		for(size_t j = 0; j < n; ++j){
			double t;
			pnrj_diagonal(n, lda, A, j, &t);
			pnrj_diagonal_L(n, lda, A, j, t);
			pnrj_diagonal_R(n, lda, A, j, t);
			if(NULL != vr){
				pnrj_diagonal_R(n, ldvr, vr, j, t);
			}
			if(NULL != vl){
				pnrj_diagonal_R(n, ldvl, vl, j, 1./t);
			}
		}
	}while(iter++ <= max_iter);
	
	int info = 0;
	if(!converged){
		info = iter;
	}

	for(size_t p = 0; p < n; p++){
		eval[p] = A[p+p*lda];
	}

	return info;
}



////////////// BEGIN ZGEEV

#include "TBLAS.h"
#include "TLASupport.h"
#include "Eigensystems.h"

static inline double abs1(const std::complex<double> &z){
	return std::abs(z.real()) + std::abs(z.imag());
}

/*
Copyright (c) 1992-2010 The University of Tennessee.  All rights reserved.

Additional copyrights may follow

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer. 
  
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.
  
- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
  
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

*/

static void zgebal_(char job, size_t n, std::complex<double> *a, size_t lda, size_t *ilo, size_t *ihi, double *scale)
{
	using namespace std;
	// Purpose
	// =======

	// ZGEBAL balances a general complex matrix A.  This involves, first,
	// permuting A by a similarity transformation to isolate eigenvalues
	// in the first 1 to ILO-1 and last IHI+1 to N elements on the
	// diagonal; and second, applying a diagonal similarity transformation
	// to rows and columns ILO to IHI to make the rows and columns as
	// close in norm as possible.  Both steps are optional.

	// Balancing may reduce the 1-norm of the matrix, and improve the
	// accuracy of the computed eigenvalues and/or eigenvectors.

	// Arguments
	// =========

	// JOB     Specifies the operations to be performed on A:
	//         = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0
	//                 for i = 1,...,N;
	//         = 'P':  permute only;
	//         = 'S':  scale only;
	//         = 'B':  both permute and scale.

	// N       The order of the matrix A.  N >= 0.

	// A       (input/output) COMPLEX*16 array, dimension (LDA,N)
	//         On entry, the input matrix A.
	//         On exit,  A is overwritten by the balanced matrix.
	//         If JOB = 'N', A is not referenced.
	//         See Further Details.

	// LDA     The leading dimension of the array A.  LDA >= max(1,N).

	// ILO     (output) INTEGER
	// IHI     (output) INTEGER
	//         ILO and IHI are set to integers such that on exit
	//         A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.
	//         If JOB = 'N' or 'S', ILO = 1 and IHI = N.

	// SCALE   (output) DOUBLE PRECISION array, dimension (N)
	//         Details of the permutations and scaling factors applied to
	//         A.  If P(j) is the index of the row and column interchanged
	//         with row and column j and D(j) is the scaling factor
	//         applied to row and column j, then
	//         SCALE(j) = P(j)    for j = 1,...,ILO-1
	//                  = D(j)    for j = ILO,...,IHI
	//                  = P(j)    for j = IHI+1,...,N.
	//         The order in which the interchanges are made is N to IHI+1,
	//         then 1 to ILO-1.

	// Further Details
	// ===============

	// The permutations consist of row and column interchanges which put
	// the matrix in the form

	//            ( T1   X   Y  )
	//    P A P = (  0   B   Z  )
	//            (  0   0   T2 )

	// where T1 and T2 are upper triangular matrices whose eigenvalues lie
	// along the diagonal.  The column indices ILO and IHI mark the starting
	// and ending columns of the submatrix B. Balancing consists of applying
	// a diagonal similarity transformation inv(D) * B * D to make the
	// 1-norms of each row of B and its corresponding column nearly equal.
	// The output matrix is

	//    ( T1     X*D          Y    )
	//    (  0  inv(D)*B*D  inv(D)*Z ).
	//    (  0      0           T2   )

	// Information about the permutations P and the diagonal matrix D is
	// returned in the vector SCALE.

	// This subroutine is based on the EISPACK routine CBAL.

	// Modified by Tzu-Yi Chen, Computer Science Division, University of
	//   California at Berkeley, USA

	// Parameter adjustments
	const size_t a_offset = 1 + lda;
	a -= a_offset;
	--scale;

	size_t k, l;
	*ilo = k = 1;
	*ihi = l = n;

	if(n == 0){
		return;
	}

	if(job == 'N') {
		for(size_t i = 1; i <= n; ++i){
			scale[i] = 1.;
		}
		return;
	}

	if(job != 'S') {
		// Permutation to isolate eigenvalues if possible
		// Search for rows isolating an eigenvalue and push them down.

		for(size_t j = l; j >= 1; --j){
			bool found_nonzero = false;
			for(size_t i = 1; i <= l; ++i){
				if(i == j){
					continue;
				}
				if(a[j+i*lda] != 0.) {
					found_nonzero = true;
					break;
				}
			}
			if(!found_nonzero){
				scale[l] = (double)j;
				if(j != l){
					RNP::TBLAS::Swap(l, &a[1+j*lda], 1, &a[1+l*lda], 1);
					RNP::TBLAS::Swap(n-k+1, &a[j+k*lda], lda, &a[l+k*lda], lda);
				}
				if(l == 1){
					*ilo = k;
					*ihi = l;
					return;
				}
				--l;
				j = l;
			}
		}

		// Search for columns isolating an eigenvalue and push them left.
		for(size_t j = k; j <= l; ++j){
			bool found_nonzero = false;
			for(size_t i = k; i <= l; ++i){
				if(i == j){
					continue;
				}
				if(a[i+j*lda] != 0.) {
					found_nonzero = true;
					break;
				}
			}
			if(!found_nonzero){
				scale[k] = (double)j;
				if(j != k){
					RNP::TBLAS::Swap(l, &a[1+j*lda], 1, &a[1+k*lda], 1);
					RNP::TBLAS::Swap(n-k+1, &a[j+k*lda], lda, &a[k+k*lda], lda);
				}
				++k;
				j = k;
			}
		}
	}

	for(size_t i = k; i <= l; ++i){
		scale[i] = 1.;
	}

	if(job == 'P') {
		*ilo = k;
		*ihi = l;
		return;
	}

	// Balance the submatrix in rows K to L.
	// Iterative loop for norm reduction
	
	const double sfmin1 = 0.5*std::numeric_limits<double>::min() / std::numeric_limits<double>::epsilon();
	const double sfmax1 = 1. / sfmin1;
	const double sfmin2 = sfmin1 * 2.;
	const double sfmax2 = 1. / sfmin2;

	bool noconv;
	do{
		noconv = false;

		for(size_t i = k; i <= l; ++i){
			double c = 0.;
			double r = 0.;

			for(size_t j = k; j <= l; ++j){
				if(j != i){
					c += abs1(a[j+i*lda]);
					r += abs1(a[i+j*lda]);
				}
			}
			size_t ica = 1+RNP::TBLAS::MaximumIndex(l, &a[1+i*lda], 1);
			double ca = std::abs(a[ica+i*lda]);
			size_t ira = 1+RNP::TBLAS::MaximumIndex(n-k+1, &a[i+k*lda], lda);
			double ra = std::abs(a[i+(ira+k-1)*lda]);

			// Guard against zero C or R due to underflow.
			if(c == 0. || r == 0.) {
				continue;
			}
			double g = r / 2.;
			double f = 1.;
			double s = c + r;
			while(!(c >= g || max(max(f,c),ca) >= sfmax2 || min(min(r,g),ra) <= sfmin2)){
				f *= 2.;
				c *= 2.;
				ca *= 2.;
				r /= 2.;
				g /= 2.;
				ra /= 2.;
			}

			g = c / 2.;
			while(!(g < r || max(r,ra) >= sfmax2 || min(min(min(f,c),g),ca) <= sfmin2)){
				f /= 2.;
				c /= 2.;
				g /= 2.;
				ca /= 2.;
				r *= 2.;
				ra *= 2.;
			}

			// Now balance.
			if(c + r >= s * 0.95){
				continue;
			}
			if(f < 1. && scale[i] < 1.) {
				if(f * scale[i] <= sfmin1){
					continue;
				}
			}
			if(f > 1. && scale[i] > 1.) {
				if(scale[i] >= sfmax1 / f){
					continue;
				}
			}
			g = 1. / f;
			scale[i] *= f;
			noconv = true;

			RNP::TBLAS::Scale(n-k+1, g, &a[i+k*lda], lda);
			RNP::TBLAS::Scale(l, f, &a[1+i*lda], 1);
		}

	}while(noconv);

	*ilo = k;
	*ihi = l;
}


static void zgebak_(char job, char side, size_t n, size_t ilo, size_t ihi, double *scale, size_t m, std::complex<double> *v, size_t ldv)
{

	// Purpose
	// =======

	// ZGEBAK forms the right or left eigenvectors of a complex general
	// matrix by backward transformation on the computed eigenvectors of the
	// balanced matrix output by ZGEBAL.

	// Arguments
	// =========

	// JOB     (input) CHARACTER*1
	//         Specifies the type of backward transformation required:
	//         = 'N', do nothing, return immediately;
	//         = 'P', do backward transformation for permutation only;
	//         = 'S', do backward transformation for scaling only;
	//         = 'B', do backward transformations for both permutation and
	//                scaling.
	//         JOB must be the same as the argument JOB supplied to ZGEBAL.

	// SIDE    (input) CHARACTER*1
	//         = 'R':  V contains right eigenvectors;
	//         = 'L':  V contains left eigenvectors.

	// N       The number of rows of the matrix V.  N >= 0.

	// ILO     (input) INTEGER
	// IHI     (input) INTEGER
	//         The integers ILO and IHI determined by ZGEBAL.
	//         1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.

	// SCALE   (input) DOUBLE PRECISION array, dimension (N)
	//         Details of the permutation and scaling factors, as returned
	//         by ZGEBAL.

	// M       The number of columns of the matrix V.  M >= 0.

	// V       (input/output) COMPLEX*16 array, dimension (LDV,M)
	//         On entry, the matrix of right or left eigenvectors to be
	//         transformed, as returned by ZHSEIN or ZTREVC.
	//         On exit, V is overwritten by the transformed eigenvectors.

	// LDV     The leading dimension of the array V. LDV >= max(1,N).

	const bool rightv = (side == 'R');
	const bool leftv = (side == 'L');

	if(n == 0 || m == 0 || job == 'N'){
		return;
	}

	if(ilo != ihi){ // Backward balance
		if(job == 'S' || job == 'B'){
			if(rightv){
				for(size_t i = ilo; i <= ihi; ++i){
					RNP::TBLAS::Scale(m, scale[i], &v[i+0*ldv], ldv);
				}
			}
			if(leftv){
				for(size_t i = ilo; i <= ihi; ++i){
					RNP::TBLAS::Scale(m, 1. / scale[i], &v[i+0*ldv], ldv);
				}
			}

		}
	}
	
	if(job == 'P' || job == 'B'){
		// Backward permutation
		// For  I = ILO-1 step -1 until 1,
		//        IHI+1 step 1 until N do
		if(rightv){
			for(size_t ii = 0; ii < n; ++ii){
				size_t i = ii;
				if(i >= ilo && i <= ihi){
					continue;
				}
				if(i < ilo){
					i = ilo - ii;
				}
				size_t k = (int)scale[i];
				if(k != i){
					RNP::TBLAS::Swap(m, &v[i+0*ldv], ldv, &v[k+0*ldv], ldv);
				}
			}
		}

		if(leftv){
			for(size_t ii = 0; ii < n; ++ii){
				size_t i = ii;
				if(i >= ilo && i <= ihi){
					continue;
				}
				if(i < ilo){
					i = ilo - ii;
				}
				size_t k = (int)scale[i];
				if(k != i){
					RNP::TBLAS::Swap(m, &v[i+0*ldv], ldv, &v[k+0*ldv], ldv);
				}
			}
		}
	}
}

static inline int iparmq_(int ispec, size_t n, size_t ilo, size_t ihi){
	// Purpose
	// =======

	//      This program sets problem and machine dependent parameters
	//      useful for xHSEQR and its subroutines. It is called whenever
	//      ILAENV is called with 12 <= ISPEC <= 16

	// Arguments
	// =========

	//      ISPEC  (input) integer scalar
	//             ISPEC specifies which tunable parameter IPARMQ should
	//             return.

	//             ISPEC=12: (INMIN)  Matrices of order nmin or less
	//                       are sent directly to xLAHQR, the implicit
	//                       double shift QR algorithm.  NMIN must be
	//                       at least 11.

	//             ISPEC=13: (INWIN)  Size of the deflation window.
	//                       This is best set greater than or equal to
	//                       the number of simultaneous shifts NS.
	//                       Larger matrices benefit from larger deflation
	//                       windows.

	//             ISPEC=14: (INIBL) Determines when to stop nibbling and
	//                       invest in an (expensive) multi-shift QR sweep.
	//                       If the aggressive early deflation subroutine
	//                       finds LD converged eigenvalues from an order
	//                       NW deflation window and LD.GT.(NW*NIBBLE)/100,
	//                       then the next QR sweep is skipped and early
	//                       deflation is applied immediately to the
	//                       remaining active diagonal block.  Setting
	//                       IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
	//                       multi-shift QR sweep whenever early deflation
	//                       finds a converged eigenvalue.  Setting
	//                       IPARMQ(ISPEC=14) greater than or equal to 100
	//                       prevents TTQRE from skipping a multi-shift
	//                       QR sweep.

	//             ISPEC=15: (NSHFTS) The number of simultaneous shifts in
	//                       a multi-shift QR iteration.

	//             ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
	//                       following meanings.
	//                       0:  During the multi-shift QR sweep,
	//                           xLAQR5 does not accumulate reflections and
	//                           does not use matrix-matrix multiply to
	//                           update the far-from-diagonal matrix
	//                           entries.
	//                       1:  During the multi-shift QR sweep,
	//                           xLAQR5 and/or xLAQRaccumulates reflections and uses
	//                           matrix-matrix multiply to update the
	//                           far-from-diagonal matrix entries.
	//                       2:  During the multi-shift QR sweep.
	//                           xLAQR5 accumulates reflections and takes
	//                           advantage of 2-by-2 block structure during
	//                           matrix-matrix multiplies.
	//                       (If xTRMM is slower than xGEMM, then
	//                       IPARMQ(ISPEC=16)=1 may be more efficient than
	//                       IPARMQ(ISPEC=16)=2 despite the greater level of
	//                       arithmetic work implied by the latter choice.)

	//      N       (input) integer scalar
	//              N is the order of the Hessenberg matrix H.

	//      ILO     (input) INTEGER
	//      IHI     (input) INTEGER
	//              It is assumed that H is already upper triangular
	//              in rows and columns 1:ILO-1 and IHI+1:N.

	// Further Details
	// ===============

	//      Little is known about how best to choose these parameters.
	//      It is possible to use different values of the parameters
	//      for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.

	//      It is probably best to choose different parameters for
	//      different matrices and different parameters at different
	//      times during the iteration, but this has not been
	//      implemented --- yet.


	//      The best choices of most of the parameters depend
	//      in an ill-understood way on the relative execution
	//      rate of xLAQR3 and xLAQR5 and on the nature of each
	//      particular eigenvalue problem.  Experiment may be the
	//      only practical way to determine which choices are most
	//      effective.

	//      Following is a list of default values supplied by IPARMQ.
	//      These defaults may be adjusted in order to attain better
	//      performance in any particular computational environment.

	//      IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
	//                       Default: 75. (Must be at least 11.)

	//      IPARMQ(ISPEC=13) Recommended deflation window size.
	//                       This depends on ILO, IHI and NS, the
	//                       number of simultaneous shifts returned
	//                       by IPARMQ(ISPEC=15).  The default for
	//                       (IHI-ILO+1).LE.500 is NS.  The default
	//                       for (IHI-ILO+1).GT.500 is 3*NS/2.

	//      IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.

	//      IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
	//                       a multi-shift QR iteration.

	//                       If IHI-ILO+1 is ...

	//                       greater than      ...but less    ... the
	//                       or equal to ...      than        default is

	//                               0               30       NS =   2+
	//                              30               60       NS =   4+
	//                              60              150       NS =  10
	//                             150              590       NS =  **
	//                             590             3000       NS =  64
	//                            3000             6000       NS = 128
	//                            6000             infinity   NS = 256

	//                   (+)  By default matrices of this order are
	//                        passed to the implicit double shift routine
	//                        xLAHQR.  See IPARMQ(ISPEC=12) above.   These
	//                        values of NS are used only in case of a rare
	//                        xLAHQR failure.

	//                   (**) The asterisks (**) indicate an ad-hoc
	//                        function increasing from 10 to 64.

	//      IPARMQ(ISPEC=16) Select structured matrix multiply.
	//                       (See ISPEC=16 above for details.)
	//                       Default: 3.

	using namespace std;
	
	size_t nh, ns;
    if (ispec == 15 || ispec == 13 || ispec == 16) {
		// Set the number simultaneous shifts
		nh = ihi - ilo + 1;
		ns = 2;
		if (nh >= 30) {
			ns = 4;
		}
		if (nh >= 60) {
			ns = 10;
		}
		if (nh >= 150) {
			ns = max((size_t)10, nh / int(log((double)nh)/log(2.) + 0.5));
		}
		if (nh >= 590) {
			ns = 64;
		}
		if (nh >= 3000) {
			ns = 128;
		}
		if (nh >= 6000) {
			ns = 256;
		}
		ns = max((size_t)2, ns - ns % 2);
    }

	int ret_val;
    if (ispec == 12) {
		// Matrices of order smaller than NMIN get sent
		// to xLAHQR, the classic double shift algorithm.
		// This must be at least 11. ====
		ret_val = 75;
    } else if (ispec == 14) {
		// INIBL: skip a multi-shift qr iteration and
		// whenever aggressive early deflation finds
		// at least (NIBBLE*(window size)/100) deflations. ====
		ret_val = 14;
    } else if (ispec == 15) { // NSHFTS: The number of simultaneous shifts
		ret_val = ns;
    } else if (ispec == 13) { // NW: deflation window size.
		if (nh <= 500) {
			ret_val = ns;
		} else {
			ret_val = ns * 3 / 2;
		}
    } else if (ispec == 16) {
		// IACC22: Whether to accumulate reflections
		// before updating the far-from-diagonal elements
		// and whether to use 2-by-2 block structure while
		// doing it.  A small amount of work could be saved
		// by making this choice dependent also upon the
		// NH=IHI-ILO+1.
		ret_val = 0;
		/*
		if (ns >= 14) {
			ret_val = 1;
		}
		*/
		if (ns >= 14) {
			ret_val = 2;
		}
    } else {
		ret_val = -1;
    }
    return ret_val;
}

static size_t zlaqr0_(bool wantt, bool wantz, size_t n, size_t ilo, size_t ihi, std::complex<double> *_h__, size_t ldh, std::complex<double> *_w, size_t iloz, size_t ihiz, std::complex<double> *_z__, size_t ldz, std::complex<double> *_work, bool is_zero);

size_t zlahqr_(bool wantt, bool wantz, size_t n, size_t ilo, size_t ihi, std::complex<double> *h, int ldh, std::complex<double> *w, size_t iloz, size_t ihiz, std::complex<double> *z, size_t ldz){
	using namespace std;
    // 
    // -- LAPACK auxiliary routine (version 3.2) --
    // Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
    // November 2006
    // 
    // .. Scalar Arguments ..
    // ..
    // .. Array Arguments ..
    // ..
    // 
    // Purpose
    // =======
    // 
    // ZLAHQR is an auxiliary routine called by CHSEQR to update the
    // eigenvalues and Schur decomposition already computed by CHSEQR, by
    // dealing with the Hessenberg submatrix in rows and columns ILO to
    // IHI.
    // 
    // Arguments
    // =========
    // 
    // WANTT   (input) bool
    // = true : the full Schur form T is required;
    // = false: only eigenvalues are required.
    // 
    // WANTZ   (input) bool
    // = true : the matrix of Schur vectors Z is required;
    // = false: Schur vectors are not required.
    // 
    // N       (input) int
    // The order of the matrix H.  N >= 0.
    // 
    // ILO     (input) int
    // IHI     (input) int
    // It is assumed that H is already upper triangular in rows and
    // columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1).
    // ZLAHQR works primarily with the Hessenberg submatrix in rows
    // and columns ILO to IHI, but applies transformations to all of
    // H if WANTT is true.
    // 1 <= ILO <= max(1,IHI); IHI <= N.
    // 
    // H       (input/output) std::complex<double> array, dimension (LDH,N)
    // On entry, the upper Hessenberg matrix H.
    // On exit, if INFO is zero and if WANTT is true, then H
    // is upper triangular in rows and columns ILO:IHI.  If INFO
    // is zero and if WANTT is false, then the contents of H
    // are unspecified on exit.  The output state of H in case
    // INF is positive is below under the description of INFO.
    // 
    // LDH     (input) int
    // The leading dimension of the array H. LDH >= max(1,N).
    // 
    // W       (output) std::complex<double> array, dimension (N)
    // The computed eigenvalues ILO to IHI are stored in the
    // corresponding elements of W. If WANTT is true, the
    // eigenvalues are stored in the same order as on the diagonal
    // of the Schur form returned in H, with W(i) = H(i,i).
    // 
    // ILOZ    (input) int
    // IHIZ    (input) int
    // Specify the rows of Z to which transformations must be
    // applied if WANTZ is true.
    // 1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
    // 
    // Z       (input/output) std::complex<double> array, dimension (LDZ,N)
    // If WANTZ is true, on entry Z must contain the current
    // matrix Z of transformations accumulated by CHSEQR, and on
    // exit Z has been updated; transformations are applied only to
    // the submatrix Z(ILOZ:IHIZ,ILO:IHI).
    // If WANTZ is false, Z is not referenced.
    // 
    // LDZ     (input) int
    // The leading dimension of the array Z. LDZ >= max(1,N).
    // 
    // INFO    (output) int
    // =   0: successful exit
    // >  0: if INFO = i, ZLAHQR failed to compute all the
    // eigenvalues ILO to IHI in a total of 30 iterations
    // per eigenvalue; elements i+1:ihi of W contain
    // those eigenvalues which have been successfully
    // computed.
    // 
    // If INFO  >  0 and WANTT is false, then on exit,
    // the remaining unconverged eigenvalues are the
    // eigenvalues of the upper Hessenberg matrix
    // rows and columns ILO thorugh INFO of the final,
    // output value of H.
    // 
    // If INFO  >  0 and WANTT is true, then on exit
    // (*)       (initial value of H)*U  = U*(final value of H)
    // where U is an orthognal matrix.    The final
    // value of H is upper Hessenberg and triangular in
    // rows and columns INFO+1 through IHI.
    // 
    // If INFO  >  0 and WANTZ is true, then on exit
    // (final value of Z)  = (initial value of Z)*U
    // where U is the orthogonal matrix in (*)
    // (regardless of the value of WANTT.)
    // 
    // Further Details
    // ===============
    // 
    // 02-96 Based on modifications by
    // David Day, Sandia National Laboratory, USA
    // 
    // 12-04 Further modifications by
    // Ralph Byers, University of Kansas, USA
    // This is a modified version of ZLAHQR from LAPACK version 3.0.
    // It is (1) more robust against overflow and underflow and
    // (2) adopts the more conservative Ahues & Tisseur stopping
    // criterion (LAWN 122, 1997).
    // 
    // =========================================================
    // 
    // .. Parameters ..
	const size_t itmax = 30 ;
	// parameter          ( itmax = 30 )
	const std::complex<double> zero = std::complex<double>(0.0e0, 0.0e0);
	const std::complex<double> one = std::complex<double>(1.0e0, 0.0e0);
	// parameter          ( zero = ( 0.0e0, 0.0e0 ), one = ( 1.0e0, 0.0e0 ) )
	const double rzero = 0.0e0;
	const double half = 0.5e0 ;
	// parameter          ( rzero = 0.0e0, rone = 1.0e0, half = 0.5e0 )
	const double dat1 = 3.0e0 / 4.0e0 ;
	// parameter          ( dat1 = 3.0e0 / 4.0e0 )
	// ..

	size_t info = 0;
	// 
	// Quick return if possible
	// 
	if( n == 0 ) return 0;
	if( ilo == ihi ){
		w[ilo-1] = h[(ilo-1)+(ilo-1)*ldh];
		return 0;
	}
	// 
	// ==== clear out the trash ====
	for(size_t j = ilo; j+3 <= ihi; ++j){
		h[(( j+2)-1)+(j-1)*ldh] = zero;
		h[(( j+3)-1)+(j-1)*ldh] = zero;
	}
	if( ilo <= ihi-2 ) h[(ihi-1)+(( ihi-2 )-1)*ldh] = zero;
	// ==== ensure that subdiagonal entries are real ====
	size_t jhi;
	size_t jlo;
	if( wantt ){
		jlo = 1;
		jhi = n;
	}else{
		jlo = ilo;
		jhi = ihi;
	}
	for(size_t i = ilo + 1; i <= ihi; ++i){
		if( std::imag( h[(i-1)+(( i-1 )-1)*ldh] ) != rzero ){
			// ==== The following redundant normalization
			// .    avoids problems with both gradual and
			// .    sudden underflow in std::abs(H(I,I-1)) ====
			std::complex<double> sc = h[(i-1)+(( i-1 )-1)*ldh] / abs1( h[(i-1)+(( i-1 )-1)*ldh] );
			sc = std::conj( sc ) / std::abs( sc );
			h[(i-1)+(( i-1 )-1)*ldh] = std::abs( h[(i-1)+(( i-1 )-1)*ldh] );
			RNP::TBLAS::Scale( jhi-i+1, sc, &h[(i-1)+(i-1)*ldh], ldh );
			RNP::TBLAS::Scale( min( jhi, i+1 )-jlo+1, std::conj( sc ), &h[(jlo-1)+(i-1)*ldh], 1 );
			if( wantz ) RNP::TBLAS::Scale( ihiz-iloz+1, std::conj( sc ), &z[(iloz-1)+(i-1)*ldz], 1 );
		}
	}
	const size_t nh = ihi - ilo + 1;
	const size_t nz = ihiz - iloz + 1;
	
	// Set machine-dependent constants for the stopping criterion.
	const double safmin = std::numeric_limits<double>::min();
	//const double safmax = rone / safmin;
	const double ulp = 2*std::numeric_limits<double>::epsilon();
	const double smlnum = safmin*( double( nh ) / ulp );
	// 
	// I1 and I2 are the indices of the first row and last column of H
	// to which transformations must be applied. If eigenvalues only are
	// being computed, I1 and I2 are set inside the main loop.
	// 
	size_t i1, i2;
	if( wantt ){
		i1 = 1;
		i2 = n;
	}
	// 
	// The main loop begins here. I is the loop index and decreases from
	// IHI to ILO in steps of 1. Each iteration of the loop works
	// with the active submatrix in rows and columns L to I.
	// Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
	// H(L,L-1) is negligible so that the matrix splits.
	// 
	size_t i = ihi;
	while(1){
		if( i < ilo ){ break; }
		// 
		// Perform QR iterations on rows and columns ILO to I until a
		// submatrix of order 1 splits off at the bottom because a
		// subdiagonal element has become negligible.
		// 
		size_t l = ilo;
		bool converged = false;
		for(size_t its = 0; its < itmax; ++its){
			// Look for a single small subdiagonal element.
			int k;
			for(k = i; k >= (int)l + 1; --k){
				if( abs1( h[(k-1)+(( k-1 )-1)*ldh] ) <= smlnum ){ break; }
				double tst = abs1( h[(k-2)+(( k-1 )-1)*ldh] ) + abs1( h[(k-1)+(k-1)*ldh] );
				if( tst == zero ){
					if( (size_t)k-2 >= ilo ) tst = tst + std::abs( std::real( h[(k-2)+(( k-2 )-1)*ldh] ) );
					if( (size_t)k+1 <= ihi ) tst = tst + std::abs( std::real( h[k+(k-1)*ldh] ) );
				}
				// ==== The following is a conservative small subdiagonal
				// .    deflation criterion due to Ahues & Tisseur (LAWN 122,
				// .    1997). It has better mathematical foundation and
				// .    improves accuracy in some examples.  ====
				if( std::abs( std::real( h[(k-1)+(( k-1 )-1)*ldh] ) ) <= ulp*tst ){
					double ab = max( abs1( h[(k-1)+(( k-1 )-1)*ldh] ), abs1( h[(k-2)+(k-1)*ldh] ) );
					double ba = min( abs1( h[(k-1)+(( k-1 )-1)*ldh] ), abs1( h[(k-2)+(k-1)*ldh] ) );
					double aa = max( abs1( h[(k-1)+(k-1)*ldh] ), abs1( h[(k-2)+(( k-1 )-1)*ldh]-h[(k-1)+(k-1)*ldh] ) );
					double bb = min( abs1( h[(k-1)+(k-1)*ldh] ), abs1( h[(k-2)+(( k-1 )-1)*ldh]-h[(k-1)+(k-1)*ldh] ) );
					double s = aa + ab;
					if( ba*( ab / s ) <= max( smlnum, ulp*( bb*( aa / s ) ) ) ){ break; }
				}
			}
			l = k;
			if( l > ilo ){ // H(L,L-1) is negligible
				h[(l-1)+(( l-1 )-1)*ldh] = zero;
			}
			// 
			// Exit from loop if a submatrix of order 1 has split off.
			// 
			if( l >= i ){
				// H(I,I-1) is negligible: one eigenvalue has converged.
				w[i-1] = h[(i-1)+(i-1)*ldh];
				// return to start of the main loop with new value of I.
				i = l - 1;
				converged = true;
				break;
			}
			// 
			// Now the active submatrix is in rows and columns L to I. If
			// eigenvalues only are being computed, only the active submatrix
			// need be transformed.
			// 
			if(  !wantt ){
				i1 = l;
				i2 = i;
			}
			std::complex<double> t;
			if( its == 10 ){ // Exceptional shift.
				double s = dat1*std::abs( std::real( h[l+(l-1)*ldh] ) );
				t = s + h[(l-1)+(l-1)*ldh];
			}else if( its == 20 ){ // Exceptional shift.
				double s = dat1*std::abs( std::real( h[(i-1)+(( i-1 )-1)*ldh] ) );
				t = s + h[(i-1)+(i-1)*ldh];
			}else{ // Wilkinson's shift.
				t = h[(i-1)+(i-1)*ldh];
				std::complex<double> u = sqrt( h[(i-2)+(i-1)*ldh] )*sqrt( h[(i-1)+(( i-1 )-1)*ldh] );
				double s = abs1( u );
				if( s != rzero ){
					std::complex<double> x = half*( h[(i-2)+(( i-1 )-1)*ldh]-t );
					double sx = abs1( x );
					s = max( s, abs1( x ) );
					std::complex<double> y = s*sqrt( pow( x / s , 2)+pow( u / s , 2) );
					if( sx > rzero ){
						if( std::real( x / sx )*std::real( y )+std::imag( x / sx )* std::imag( y ) < rzero )y = -y;
					}
					t = t - u*RNP::TLASupport::ComplexDivide( u, ( x+y ) );
				}
			}
			// 
			// Look for two consecutive small subdiagonal elements.
			// 
			bool breaked_out = false;
			int m;
			std::complex<double> v[2];
			for(m = i - 1; m >= (int)l + 1; --m){
				// 
				// Determine the effect of starting the single-shift QR
				// iteration at row M, and see if this would make H(M,M-1)
				// negligible.
				// 
				std::complex<double> h11 = h[(m-1)+(m-1)*ldh];
				std::complex<double> h22 = h[m+(( m+1 )-1)*ldh];
				std::complex<double> h11s = h11 - t;
				double h21 = std::real( h[m+(m-1)*ldh] );
				double s = abs1( h11s ) + std::abs( h21 );
				h11s = h11s / s;
				h21 = h21 / s;
				v[0] = h11s;
				v[1] = h21;
				double h10 = std::real( h[(m-1)+(( m-1 )-1)*ldh] );
				if( std::abs( h10 )*std::abs( h21 ) <= ulp* ( abs1( h11s )*( abs1( h11 )+abs1( h22 ) ) ) ){
					breaked_out = true;
					break;
				}
			}
			if(!breaked_out){
				 std::complex<double> h11 = h[(l-1)+(l-1)*ldh];
				 std::complex<double> h22 = h[l+(( l+1 )-1)*ldh];
				 std::complex<double> h11s = h11 - t;
				 double h21 = std::real( h[l+(l-1)*ldh] );
				 double s = abs1( h11s ) + std::abs( h21 );
				 h11s = h11s / s;
				 h21 = h21 / s;
				 v[0] = h11s;
				 v[1] = h21;
			}
			// 
			// Single-shift QR step
			// 
			for(k = m; (size_t)k < i; ++k){
				// 
				// The first iteration of this loop determines a reflection G
				// from the vector V and applies it from left and right to H,
				// thus creating a nonzero bulge below the subdiagonal.
				// 
				// Each subsequent iteration determines a reflection G to
				// restore the Hessenberg form in the (K-1)th column, and thus
				// chases the bulge one step toward the bottom of the active
				// submatrix.
				// 
				// V(2) is always real before the call to ZLARFG, and hence
				// after the call T2 ( = T1*V(2) ) is also real.
				// 
				if( k > m ) RNP::TBLAS::Copy( 2, &h[(k-1)+(( k-1 )-1)*ldh], 1, v, 1 );
				std::complex<double> t1;
				RNP::TLASupport::GenerateElementaryReflector( 2, &v[0], &v[1], 1, &t1 );
				if( k > m ){
					h[(k-1)+(( k-1 )-1)*ldh] = v[0];
					h[k+(( k-1 )-1)*ldh] = zero;
				}
				std::complex<double> v2 = v[1];
				double t2 = std::real( t1*v2 );
				// 
				// Apply G from the left to transform the rows of the matrix
				// in columns K to I2.
				// 
				for(size_t j = k; j <= i2; ++j){
					std::complex<double> sum = std::conj( t1 )*h[(k-1)+(j-1)*ldh] + t2*h[k+(j-1)*ldh];
					h[(k-1)+(j-1)*ldh] = h[(k-1)+(j-1)*ldh] - sum;
					h[k+(j-1)*ldh] = h[k+(j-1)*ldh] - sum*v2;
				}
				// 
				// Apply G from the right to transform the columns of the
				// matrix in rows I1 to min(K+2,I).
				// 
				for(size_t j = i1; j <= min( (size_t)k+2, i ); ++j){
					std::complex<double> sum = t1*h[(j-1)+(k-1)*ldh] + t2*h[(j-1)+(( k+1 )-1)*ldh];
					h[(j-1)+(k-1)*ldh] = h[(j-1)+(k-1)*ldh] - sum;
					h[(j-1)+(( k+1 )-1)*ldh] = h[(j-1)+(( k+1 )-1)*ldh] - sum*std::conj( v2 );
				}
				// 
				if( wantz ){ // Accumulate transformations in the matrix Z
					for(size_t j = iloz; j <= ihiz; ++j){
						std::complex<double> sum = t1*z[(j-1)+(k-1)*ldz] + t2*z[(j-1)+(( k+1 )-1)*ldz];
						z[(j-1)+(k-1)*ldz] = z[(j-1)+(k-1)*ldz] - sum;
						z[(j-1)+(( k+1 )-1)*ldz] = z[(j-1)+(( k+1 )-1)*ldz] - sum*std::conj( v2 );
					}
				}
				// 
				if( k == m  &&  m > (int)l ){
					// 
					// If the QR step was started at row M > L because two
					// consecutive small subdiagonals were found, then extra
					// scaling must be performed to ensure that H(M,M-1) remains
					// real.
					// 
					std::complex<double> temp = one - t1;
					temp = temp / std::abs( temp );
					h[m+(m-1)*ldh] = h[m+(m-1)*ldh]*std::conj( temp );
					if( (size_t)m+2 <= i ) h[(( m+2)-1)+(( m+1 )-1)*ldh] = h[(( m+2)-1)+(( m+1 )-1)*ldh]*temp;
					for(size_t j = m; j <= i; ++j){
						if( j != (size_t)m+1 ){
							if( i2 > j ) RNP::TBLAS::Scale( i2-j, temp, &h[(j-1)+(( j+1 )-1)*ldh], ldh );
							RNP::TBLAS::Scale( j-i1, std::conj( temp ), &h[(i1-1)+(j-1)*ldh], 1 );
							if( wantz ){
								RNP::TBLAS::Scale( nz, std::conj( temp ), &z[(iloz-1)+(j-1)*ldz], 1 );
							}
						}
					}
				}
			}

			// Ensure that H(I,I-1) is real.
			std::complex<double> temp = h[(i-1)+(( i-1 )-1)*ldh];
			if( std::imag( temp ) != rzero ){
				double rtemp = std::abs( temp );
				h[(i-1)+(( i-1 )-1)*ldh] = rtemp;
				temp = temp / rtemp;
				if( i2 > i ) RNP::TBLAS::Scale( i2-i, std::conj( temp ), &h[(i-1)+(( i+1 )-1)*ldh], ldh );
				RNP::TBLAS::Scale( i-i1, temp, &h[(i1-1)+(i-1)*ldh], 1 );
				if( wantz ){
					RNP::TBLAS::Scale( nz, temp, &z[(iloz-1)+(i-1)*ldz], 1 );
				}
			}
		}
		// Failure to converge in remaining number of iterations
		if(!converged){
			return i;
		}
	}
	return info;
}

void ztrexc_(int n, std::complex<double> *t, int ldt, std::complex<double> *q, int ldq, int ifst, int ilst){
	using namespace std;
    // 
    // -- LAPACK routine (version 3.2) --
    // -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // November 2006
    // 
    // .. Scalar Arguments ..
    // ..
    // .. Array Arguments ..
    // ..
    // 
    // Purpose
    // =======
    // 
    // ZTREXC reorders the Schur factorization of a complex matrix
    // A = Q*T*Q**H, so that the diagonal element of T with row index IFST
    // is moved to row ILST.
    // 
    // The Schur form T is reordered by a unitary similarity transformation
    // Z**H*T*Z, and optionally the matrix Q of Schur vectors is updated by
    // postmultplying it with Z.
    // 
    // Arguments
    // =========
    // 
    // COMPQ   (input) char*1
    // = 'V':  update the matrix Q of Schur vectors;
    // = 'N':  do not update Q.
    // 
    // N       (input) int
    // The order of the matrix T. N >= 0.
    // 
    // T       (input/output) std::complex<double> array, dimension (LDT,N)
    // On entry, the upper triangular matrix T.
    // On exit, the reordered upper triangular matrix.
    // 
    // LDT     (input) int
    // The leading dimension of the array T. LDT >= max(1,N).
    // 
    // Q       (input/output) std::complex<double> array, dimension (LDQ,N)
    // On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
    // On exit, if COMPQ = 'V', Q has been postmultiplied by the
    // unitary transformation matrix Z which reorders T.
    // If COMPQ = 'N', Q is not referenced.
    // 
    // LDQ     (input) int
    // The leading dimension of the array Q.  LDQ >= max(1,N).
    // 
    // IFST    (input) int
    // ILST    (input) int
    // Specify the reordering of the diagonal elements of T:
    // The element with row index IFST is moved to row ILST by a
    // sequence of transpositions between adjacent elements.
    // 1 <= IFST <= N; 1 <= ILST <= N.
    // 
    // INFO    (output) int
    // = 0:  successful exit
    // < 0:  if INFO = -i, the i-th argument had an illegal value
    // 
    // =====================================================================
    // 
    // .. Local Scalars ..
      int m1;
      int m2;
      int m3;

      const bool wantq = ( NULL != q );
     
      // 
      // Quick return if possible
      // 
      if( n == 1  ||  ifst == ilst ) return;
      // 
      if( ifst < ilst ){
      // 
      // Move the IFST-th diagonal element forward down the diagonal.
      // 
         m1 = 0;
         m2 = -1;
         m3 = 1;
      }else{
      // 
      // Move the IFST-th diagonal element backward up the diagonal.
      // 
         m1 = -1;
         m2 = 0;
         m3 = -1;
      }
      // 
      for(int k = ifst + m1; ((m3 < 0) ? (k >= ilst + m2) : (k <= ilst + m2)); k += m3){
      // 
      // Interchange the k-th and (k+1)-th diagonal elements.
      // 
         std::complex<double> t11 = t[(k-1)+(k-1)*ldt];
         std::complex<double> t22 = t[k+(( k+1 )-1)*ldt];
         // 
         // Determine the transformation to perform the interchange.
         // 
      double cs;
      std::complex<double> sn, temp;
         RNP::TLASupport::GeneratePlaneRotation( t[(k-1)+(( k+1 )-1)*ldt], t22-t11, &cs, &sn, &temp );
         
         // Apply transformation to the matrix T.
         if( k+2 <= n ) RNP::TLASupport::ApplyPlaneRotation( n-k-1, &t[(k-1)+(( k+2 )-1)*ldt], ldt, &t[k+(( k+2 )-1)*ldt], ldt, cs, sn );
         RNP::TLASupport::ApplyPlaneRotation( k-1, &t[0+(k-1)*ldt], 1, &t[0+(( k+1 )-1)*ldt], 1, cs, std::conj( sn ) );
         // 
         t[(k-1)+(k-1)*ldt] = t22;
         t[k+(( k+1 )-1)*ldt] = t11;
         // 
         if( wantq ){ // Accumulate transformation in the matrix Q.
            RNP::TLASupport::ApplyPlaneRotation( n, &q[0+(k-1)*ldq], 1, &q[0+(( k+1 )-1)*ldq], 1, cs, std::conj( sn ) );
         }
      }
}

static void zlaqr3_(bool wantt, bool wantz, size_t n, size_t ktop, size_t kbot, size_t nw, std::complex<double> *_h__, size_t ldh, size_t iloz, size_t ihiz, std::complex<double> *_z__, size_t ldz, size_t *ns, size_t *nd, std::complex<double> *_sh, std::complex<double> *_v, size_t ldv, size_t nh, std::complex<double> *_t, size_t ldt, size_t nv, std::complex<double> *_wv, size_t ldwv, std::complex<double> *_work, bool is_three)
{
	std::complex<double> *h__ = (std::complex<double>*)_h__;
	std::complex<double> *v = (std::complex<double>*)_v;
	std::complex<double> *t = (std::complex<double>*)_t;
	std::complex<double> *sh = (std::complex<double>*)_sh;
	std::complex<double> *z__ = (std::complex<double>*)_z__;
	std::complex<double> *work = (std::complex<double>*)_work;
	
	size_t h_offset, t_offset, v_offset, z_offset;

	// ******************************************************************
	// Aggressive early deflation:

	// This subroutine accepts as input an upper Hessenberg matrix
	// H and performs an unitary similarity transformation
	// designed to detect and deflate fully converged eigenvalues from
	// a trailing principal submatrix.  On output H has been over-
	// written by a new Hessenberg matrix that is a perturbation of
	// an unitary similarity transformation of H.  It is to be
	// hoped that the final version of H has many zero subdiagonal
	// entries.

	// ******************************************************************
	// WANTT   (input) LOGICAL
	//      If .TRUE., then the Hessenberg matrix H is fully updated
	//      so that the triangular Schur factor may be
	//      computed (in cooperation with the calling subroutine).
	//      If .FALSE., then only enough of H is updated to preserve
	//      the eigenvalues.

	// WANTZ   (input) LOGICAL
	//      If .TRUE., then the unitary matrix Z is updated so
	//      so that the unitary Schur factor may be computed
	//      (in cooperation with the calling subroutine).
	//      If .FALSE., then Z is not referenced.

	// N       (input) INTEGER
	//      The order of the matrix H and (if WANTZ is .TRUE.) the
	//      order of the unitary matrix Z.

	// KTOP    (input) INTEGER
	//      It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
	//      KBOT and KTOP together determine an isolated block
	//      along the diagonal of the Hessenberg matrix.

	// KBOT    (input) INTEGER
	//      It is assumed without a check that either
	//      KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
	//      determine an isolated block along the diagonal of the
	//      Hessenberg matrix.

	// NW      (input) INTEGER
	//      Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).

	// H       (input/output) COMPLEX*16 array, dimension (LDH,N)
	//      On input the initial N-by-N section of H stores the
	//      Hessenberg matrix undergoing aggressive early deflation.
	//      On output H has been transformed by a unitary
	//      similarity transformation, perturbed, and the returned
	//      to Hessenberg form that (it is to be hoped) has some
	//      zero subdiagonal entries.

	// LDH     (input) integer
	//      Leading dimension of H just as declared in the calling
	//      subroutine.  N .LE. LDH

	// ILOZ    (input) INTEGER
	// IHIZ    (input) INTEGER
	//      Specify the rows of Z to which transformations must be
	//      applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.

	// Z       (input/output) COMPLEX*16 array, dimension (LDZ,N)
	//      IF WANTZ is .TRUE., then on output, the unitary
	//      similarity transformation mentioned above has been
	//      accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.
	//      If WANTZ is .FALSE., then Z is unreferenced.

	// LDZ     (input) integer
	//      The leading dimension of Z just as declared in the
	//      calling subroutine.  1 .LE. LDZ.

	// NS      (output) integer
	//      The number of unconverged (ie approximate) eigenvalues
	//      returned in SR and SI that may be used as shifts by the
	//      calling subroutine.

	// ND      (output) integer
	//      The number of converged eigenvalues uncovered by this
	//      subroutine.

	// SH      (output) COMPLEX*16 array, dimension KBOT
	//      On output, approximate eigenvalues that may
	//      be used for shifts are stored in SH(KBOT-ND-NS+1)
	//      through SR(KBOT-ND).  Converged eigenvalues are
	//      stored in SH(KBOT-ND+1) through SH(KBOT).

	// V       (workspace) COMPLEX*16 array, dimension (LDV,NW)
	//      An NW-by-NW work array.

	// LDV     (input) integer scalar
	//      The leading dimension of V just as declared in the
	//      calling subroutine.  NW .LE. LDV

	// NH      (input) integer scalar
	//      The number of columns of T.  NH.GE.NW.

	// T       (workspace) COMPLEX*16 array, dimension (LDT,NW)

	// LDT     (input) integer
	//      The leading dimension of T just as declared in the
	//      calling subroutine.  NW .LE. LDT

	// NV      (input) integer
	//      The number of rows of work array WV available for
	//      workspace.  NV.GE.NW.

	// WV      (workspace) COMPLEX*16 array, dimension (LDWV,NW)

	// LDWV    (input) integer
	//      The leading dimension of W just as declared in the
	//      calling subroutine.  NW .LE. LDV

	// WORK    (workspace) COMPLEX*16 array, dimension LWORK.
	//      On exit, WORK(1) is set to an estimate of the optimal value
	//      of LWORK for the given values of N, NW, KTOP and KBOT.

	// LWORK   (input) integer
	//      The dimension of the work array WORK.  LWORK = 2nw
	//      suffices, but greater efficiency may result from larger
	//      values of LWORK.

	//      If LWORK = -1, then a workspace query is assumed; we
	//      only estimate the optimal workspace size for the given
	//      values of N, NW, KTOP and KBOT.  The estimate is returned
	//      in WORK(1).  No error message related to LWORK is issued
	//      by XERBLA.  Neither H nor Z are accessed.

	// ================================================================
	// Based on contributions by
	//    Karen Braman and Ralph Byers, Department of Mathematics,
	//    University of Kansas, USA

	// ================================================================
	
	using namespace std;

	/* Parameter adjustments */
	h_offset = 1 + ldh;
	h__ -= h_offset;
	z_offset = 1 + ldz;
	z__ -= z_offset;
	--sh;
	v_offset = 1 + ldv;
	v -= v_offset;
	t_offset = 1 + ldt;
	t -= t_offset;
	--work;

	/* Function Body */
/* Computing MIN */
	size_t jw = min(nw, kbot - ktop + 1);

/*     ==== Nothing to do ... */
/*     ... for an empty active block ... ==== */
	*ns = 0;
	*nd = 0;
	if (ktop > kbot) {
		return;
	}
/*     ... nor for an empty deflation window. ==== */
	if (nw < 1) {
		return;
	}

/*     ==== Machine constants ==== */

	const double safmin = std::numeric_limits<double>::min();
	//const double safmax = 1. / safmin;
	const double ulp = 2*std::numeric_limits<double>::epsilon();
	const double smlnum = safmin * ((double) (n) / ulp);

/*     ==== Setup deflation window ==== */

/* Computing MIN */
	jw = min(nw ,kbot - ktop + 1);
	size_t kwtop = kbot - jw + 1;
	std::complex<double> s;
	if (kwtop == ktop) {
		s = 0.;
	} else {
		s = h__[kwtop+(kwtop-1)*ldh];
	}

	if (kbot == kwtop) { // 1-by-1 deflation window: not much to do
		sh[kwtop] = h__[kwtop+kwtop*ldh];
		*ns = 1;
		*nd = 0;
/* Computing MAX */
		if (abs1(s) <= max(smlnum,ulp * abs1(h__[kwtop+kwtop*ldh]))) {
			*ns = 0;
			*nd = 1;
			if (kwtop > ktop) {
				h__[kwtop+(kwtop-1)*ldh] = 0.;
			}
		}
		return;
	}

/*     ==== Convert to spike-triangular form.  (In case of a */
/*     .    rare QR failure, this routine continues to do */
/*     .    aggressive early deflation using that part of */
/*     .    the deflation window that converged using INFQR */
/*     .    here and there to keep track.) ==== */

	RNP::TBLAS::CopyMatrix<'U'>(jw, jw, &h__[kwtop + kwtop * ldh], ldh, _t, ldt);
	RNP::TBLAS::Copy(jw-1, &h__[kwtop + 1 + kwtop * ldh], ldh+1, &_t[1+0*ldt], ldt+1);

	RNP::TBLAS::SetMatrix<'A'>(jw, jw, 0., 1., _v, ldv);
	size_t infqr;
	if (jw > (size_t)iparmq_(12, jw, 1, jw) && is_three) {
		infqr = zlaqr0_(true, true, jw, 1, jw, _t, ldt, &sh[kwtop], 1, jw, _v, ldv, _work, false);
	} else {
		infqr = zlahqr_(true, true, jw, 1, jw, _t, ldt, &sh[kwtop], 1, jw, _v, ldv);
	}

/*     ==== Deflation detection loop ==== */

	size_t ifst, ilst;
	*ns = jw;
	ilst = infqr + 1;
	for (size_t knt = infqr + 1; knt <= jw; ++knt) {
		// Small spike tip deflation test
		double foo = abs1(t[(*ns)+(*ns)*ldt]);
		if (foo == 0.) {
			foo = abs1(s);
		}
		if(abs1(s)*abs1(v[1+(*ns)*ldv]) <= max(smlnum,ulp * foo)){
			// One more converged eigenvalue
			--(*ns);
		} else {
			// One undeflatable eigenvalue.  Move it up out of the
			//    way.   (ZTREXC can not fail in this case.)
			ifst = *ns;
			ztrexc_(jw, _t, ldt, _v, ldv, ifst, ilst);
			++ilst;
		}
	}

/*        ==== Return to Hessenberg form ==== */

	if (*ns == 0) {
		s = 0.;
	}

	if (*ns < jw) {
		// sorting the diagonal of T improves accuracy for graded matrices.
		for (size_t i = infqr; i < *ns; ++i) {
			ifst = i;
			for (size_t j = i; j < *ns; ++j) {
				if (abs1(_t[j+j*ldt]) > abs1(_t[ifst+ifst*ldt])) {
					ifst = j;
				}
			}
			ilst = i;
			if (ifst != ilst) {
				ztrexc_(jw, _t, ldt, _v, ldv, ifst+1, ilst+1);
			}
		}
	}

/*     ==== Restore shift/eigenvalue array from T ==== */

	for (size_t i = infqr; i < jw; ++i) {
		_sh[kwtop-1 + i] = _t[i+i*ldt];
	}


	if (*ns < jw || s == 0.) {
		if (*ns > 1 && (s != 0.)) {
			// Reflect spike back into lower triangle
			RNP::TBLAS::Copy(*ns, _v, ldv, _work, 1);
			for (size_t i = 0; i < *ns; ++i) {
				_work[i] = std::conj(_work[i]);
			}
			std::complex<double> beta = work[1];
			std::complex<double> tau;
			RNP::TLASupport::GenerateElementaryReflector(*ns, &beta, &_work[1], 1, &tau);
			_work[0] = 1.;

			RNP::TBLAS::SetMatrix<'L'>(jw-2, jw-2, 0., 0., &_t[2+0*ldt], ldt);

			RNP::TLASupport::ApplyElementaryReflector<'L'>(*ns, jw, _work, 1, std::conj(tau), _t, ldt, &_work[jw]);
			RNP::TLASupport::ApplyElementaryReflector<'R'>(*ns, *ns, _work, 1, tau, _t, ldt, &_work[jw]);
			RNP::TLASupport::ApplyElementaryReflector<'R'>(jw, *ns, _work, 1, tau, _v, ldv, &_work[jw]);

			RNP::TLASupport::HessenbergReduction(jw, 0, (*ns)-1, _t, ldt, _work, &_work[jw], jw);
		}

/*        ==== Copy updated reduced window into place ==== */

		if (kwtop > 1) {
			h__[kwtop+(kwtop-1)*ldh] = s*std::conj(_v[0+0*ldv]);
		}
		RNP::TBLAS::CopyMatrix<'U'>(jw, jw, _t, ldt, &h__[kwtop + kwtop * ldh], ldh);
		RNP::TBLAS::Copy(jw-1, &t[ldt + 2], ldt+1, &h__[kwtop + 1 + kwtop * ldh], ldh+1);

/*        ==== Accumulate orthogonal matrix in order update */
/*        .    H and Z, if requested.  ==== */

		if (*ns > 1 && (s != 0.)) {
			RNP::TLASupport::ApplyOrthognalMatrixFromElementaryReflectors<'R','N'>(jw, *ns-1, *ns-1, &_t[1+0*ldt], ldt, &_work[0], &_v[0+1*ldv], ldv, &_work[jw]);
		}

/*        ==== Update vertical slab in H ==== */

		size_t ltop;
		if (wantt) {
			ltop = 1;
		} else {
			ltop = ktop;
		}
		for (size_t krow = ltop; krow < kwtop; krow += nv) {
			size_t kln = min(nv,kwtop - krow);
			RNP::TBLAS::MultMM<'N','N'>(kln, jw, jw, 1., &h__[krow + kwtop * ldh], ldh, _v, ldv, 0., _wv, ldwv);
			RNP::TBLAS::CopyMatrix<'A'>(kln, jw, _wv, ldwv, &h__[krow + kwtop * ldh], ldh);
		}

/*        ==== Update horizontal slab in H ==== */

		if (wantt) {
			for (size_t kcol = kbot + 1; kcol <= n; kcol += nh) {
				size_t kln = min(nh, n - kcol + 1);
				RNP::TBLAS::MultMM<'C','N'>(jw, kln, jw, 1., _v, ldv, &h__[kwtop + kcol * ldh], ldh, 0., _t, ldt);
				RNP::TBLAS::CopyMatrix<'A'>(jw, kln, _t, ldt, &h__[kwtop + kcol * ldh], ldh);
			}
		}

/*        ==== Update vertical slab in Z ==== */

		if (wantz) {
			for (size_t krow = iloz; krow <= ihiz; krow += nv) {
				size_t kln = min(nv, ihiz - krow + 1);
				RNP::TBLAS::MultMM<'N','N'>(kln, jw, jw, 1., &z__[krow + kwtop * ldz], ldz, _v, ldv, 0., _wv, ldwv);
				RNP::TBLAS::CopyMatrix<'A'>(kln, jw, _wv, ldwv, &z__[krow + kwtop * ldz], ldz);
			}
		}
	}

	// Return the number of deflations ...
	*nd = jw - *ns;

	// ... and the number of shifts. (Subtracting
	// INFQR from the spike length takes care
	// of the case of a rare QR failure while
	// calculating eigenvalues of the deflation
	// window.)
	*ns -= infqr;
}

void zlaqr1_(size_t n, std::complex<double> *h, size_t ldh, const std::complex<double> &s1, const std::complex<double> &s2, std::complex<double> *v){
	using namespace std;

    // Given a 2-by-2 or 3-by-3 matrix H, ZLAQR1 sets v to a
    // scalar multiple of the first column of the product
    // 
    // (*)  K = (H - s1*I)*(H - s2*I)
    // 
    // scaling to avoid overflows and most underflows.
    // 
    // This is useful for starting double implicit shift bulges
    // in the QR algorithm.
    // 
    // 
    // N      (input) int
    // Order of the matrix H. N must be either 2 or 3.
    // 
    // H      (input) std::complex<double> array of dimension (LDH,N)
    // The 2-by-2 or 3-by-3 matrix H in (*).
    // 
    // LDH    (input) int
    // The leading dimension of H as declared in
    // the calling procedure.  LDH >= N
    // 
    // S1     (input) std::complex<double>
    // S2     S1 and S2 are the shifts defining K in (*) above.
    // 
    // V      (output) std::complex<double> array of dimension N
    // A scalar multiple of the first column of the
    // matrix K in (*).
    // 
    // ================================================================
    // Based on contributions by
    // Karen Braman and Ralph Byers, Department of Mathematics,
    // University of Kansas, USA
    // 
    // ================================================================
    // 
    // .. Parameters ..
      const std::complex<double> zero = std::complex<double>(0.0e0, 0.0e0);
      // parameter          ( zero = ( 0.0e0, 0.0e0 ) )
      const double rzero = 0.0e0 ;
      // parameter          ( rzero = 0.0e0 )
      // ..
      // .. Local Scalars ..
      std::complex<double> cdum;
      std::complex<double> h21s;
      std::complex<double> h31s;
      double s;

      if( n == 2 ){
         s = abs1( h[0+0*ldh]-s2 ) + abs1( h[1+0*ldh] );
         if( s == rzero ){
            v[1-1] = zero;
            v[2-1] = zero;
         }else{
            h21s = h[1+0*ldh] / s;
            v[1-1] = h21s*h[0+1*ldh] + ( h[0+0*ldh]-s1 )* ( ( h[0+0*ldh]-s2 ) / s );
            v[2-1] = h21s*( h[0+0*ldh]+h[1+1*ldh]-s1-s2 );
         }
      }else{
         s = abs1( h[0+0*ldh]-s2 ) + abs1( h[1+0*ldh] ) + abs1( h[(3-1)+0*ldh] );
         if( s == zero ){
            v[1-1] = zero;
            v[2-1] = zero;
            v[3-1] = zero;
         }else{
            h21s = h[1+0*ldh] / s;
            h31s = h[(3-1)+0*ldh] / s;
            v[1-1] = ( h[0+0*ldh]-s1 )*( ( h[0+0*ldh]-s2 ) / s ) + h[0+1*ldh]*h21s + h[0+(3-1)*ldh]*h31s;
            v[2-1] = h21s*( h[0+0*ldh]+h[1+1*ldh]-s1-s2 ) + h[1+(3-1)*ldh]*h31s;
            v[3-1] = h31s*( h[0+0*ldh]+h[(3-1)+(3-1)*ldh]-s1-s2 ) + h21s*h[(3-1)+1*ldh];
         }
      }
}

static inline double absv(std::complex<double> zdum){ return std::abs( std::real( zdum ) / 2.0 ) + std::abs( std::imag( zdum ) / 2.0 ); }
// Declarations need repairing; reference/pointers need to be replaced.

template <char uplo, char trans, char diag, char normin>
struct SolveTrV_Scaled{
	SolveTrV_Scaled(size_t n, std::complex<double> *a, size_t lda, std::complex<double> *x, double *scale, double *cnorm){
		using namespace std;

		// Purpose
		// =======
		// 
		// ZLATRS solves one of the triangular systems
		// 
		// A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
		// 
		// with scaling to prevent overflow.  Here A is an upper or lower
		// triangular matrix, A**T denotes the transpose of A, A**H denotes the
		// conjugate transpose of A, x and b are n-element vectors, and s is a
		// scaling factor, usually less than or equal to 1, chosen so that the
		// components of x will be less than the overflow threshold.  If the
		// unscaled problem will not cause overflow, the Level 2 BLAS routine
		// ZTRSV is called. If the matrix A is singular (A(j,j) = 0 for some j),
		// then s is set to 0 and a non-trivial solution to A*x = 0 is returned.
		// 
		// Arguments
		// =========
		// 
		// UPLO    (input) char*1
		// Specifies whether the matrix A is upper or lower triangular.
		// = 'U':  Upper triangular
		// = 'L':  Lower triangular
		// 
		// TRANS   (input) char*1
		// Specifies the operation applied to A.
		// = 'N':  Solve A * x = s*b     (No transpose)
		// = 'T':  Solve A**T * x = s*b  (Transpose)
		// = 'C':  Solve A**H * x = s*b  (Conjugate transpose)
		// 
		// DIAG    (input) char*1
		// Specifies whether or not the matrix A is unit triangular.
		// = 'N':  Non-unit triangular
		// = 'U':  Unit triangular
		// 
		// NORMIN  (input) char*1
		// Specifies whether CNORM has been set or not.
		// = 'Y':  CNORM contains the column norms on entry
		// = 'N':  CNORM is not set on entry.  On exit, the norms will
		// be computed and stored in CNORM.
		// 
		// N       (input) int
		// The order of the matrix A.  N >= 0.
		// 
		// A       (input) std::complex<double> array, dimension (LDA,N)
		// The triangular matrix A.  If UPLO = 'U', the leading n by n
		// upper triangular part of the array A contains the upper
		// triangular matrix, and the strictly lower triangular part of
		// A is not referenced.  If UPLO = 'L', the leading n by n lower
		// triangular part of the array A contains the lower triangular
		// matrix, and the strictly upper triangular part of A is not
		// referenced.  If DIAG = 'U', the diagonal elements of A are
		// also not referenced and are assumed to be 1.
		// 
		// LDA     (input) int
		// The leading dimension of the array A.  LDA >= max (1,N).
		// 
		// X       (input/output) std::complex<double> array, dimension (N)
		// On entry, the right hand side b of the triangular system.
		// On exit, X is overwritten by the solution vector x.
		// 
		// SCALE   (output) double
		// The scaling factor s for the triangular system
		// A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b.
		// If SCALE = 0, the matrix A is singular or badly scaled, and
		// the vector x is an exact or approximate solution to A*x = 0.
		// 
		// CNORM   (input or output) double array, dimension (N)
		// 
		// If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
		// contains the norm of the off-diagonal part of the j-th column
		// of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
		// to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
		// must be greater than or equal to the 1-norm.
		// 
		// If NORMIN = 'N', CNORM is an output argument and CNORM(j)
		// returns the 1-norm of the offdiagonal part of the j-th column
		// of A.
		// 
		// INFO    (output) int
		// = 0:  successful exit
		// < 0:  if INFO = -k, the k-th argument had an illegal value
		// 
		// Further Details
		// ======= =======
		// 
		// A rough bound on x is computed; if that is less than overflow, ZTRSV
		// is called, otherwise, specific code is used which checks for possible
		// overflow or divide-by-zero at every operation.
		// 
		// A columnwise scheme is used for solving A*x = b.  The basic algorithm
		// if A is lower triangular is
		// 
		// x[1:n] := b[1:n]
		// for j = 1, ..., n
		// x(j) := x(j) / A(j,j)
		// x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
		// end
		// 
		// Define bounds on the components of x after j iterations of the loop:
		// M(j) = bound on x[1:j]
		// G(j) = bound on x[j+1:n]
		// Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.
		// 
		// Then for iteration j+1 we have
		// M(j+1) <= G(j) / | A(j+1,j+1) |
		// G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
		// <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )
		// 
		// where CNORM(j+1) is greater than or equal to the infinity-norm of
		// column j+1 of A, not counting the diagonal.  Hence
		// 
		// G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
		// 1<=i<=j
		// and
		// 
		// |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
		// 1<=i< j
		// 
		// Since |x(j)| <= M(j), we use the Level 2 BLAS routine ZTRSV if the
		// reciprocal of the largest M(j), j=1,..,n, is larger than
		// max(underflow, 1/overflow).
		// 
		// The bound on x(j) is also used to determine when a step in the
		// columnwise method can be performed without fear of overflow.  If
		// the computed bound is greater than a large constant, x is scaled to
		// prevent overflow, but if the bound overflows, x is set to 0, x(j) to
		// 1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
		// 
		// Similarly, a row-wise scheme is used to solve A**T *x = b  or
		// A**H *x = b.  The basic algorithm for A upper triangular is
		// 
		// for j = 1, ..., n
		// x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)
		// end
		// 
		// We simultaneously compute two bounds
		// G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
		// M(j) = bound on x(i), 1<=i<=j
		// 
		// The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we
		// add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.
		// Then the bound on x(j) is
		// 
		// M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
		// 
		// <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
		// 1<=i<=j
		// 
		// and we can safely call ZTRSV if 1/M(n) and 1/G(n) are both greater
		// than max(underflow, 1/overflow).
		// 
		// =====================================================================
		// 
		// .. Parameters ..
		  const double zero = 0.0e+0;
		  const double half = 0.5e+0;
		  const double one = 1.0e+0;
		  const double two = 2.0e+0 ;
		  // parameter          ( zero = 0.0e+0, half = 0.5e+0, one = 1.0e+0, two = 2.0e+0 )
		  // ..
		  // .. Local Scalars ..
		  bool notran;
		  bool nounit;
		  bool upper;
		  int jfirst;
		  int jinc;
		  int jlast;
		  double bignum;
		  double grow;
		  double rec;
		  double smlnum;
		  double tjj;
		  double tmax;
		  double tscal;
		  double xbnd;
		  double xj;
		  double xmax;
		  std::complex<double> csumj;
		  std::complex<double> tjjs;
		  std::complex<double> uscal;
		  std::complex<double> zdum;

		  upper = ( uplo == 'U' );
		  notran = ( trans == 'N' );
		  nounit = ( diag == 'N' );

		  // 
		  // Quick return if possible
		  // 
		  if( n == 0 ) return;
		  // 
		  // Determine machine dependent parameters to control overflow.
		  // 
		  smlnum = std::numeric_limits<double>::min();
		  smlnum = smlnum / (2*std::numeric_limits<double>::epsilon());
		  bignum = one / smlnum;
		  *scale = one;
		  // 
		  if( normin == 'N'){
		  // 
		  // Compute the 1-norm of each column, not including the diagonal.
		  // 
			 if( upper ){
			 // 
			 // A is upper triangular.
			 // 
				for(size_t j = 0; j < n; ++j){
				   cnorm[j] = RNP::TBLAS::Asum( j, &a[0+j*lda], 1 );
				}
			 }else{
			 // 
			 // A is lower triangular.
			 // 
				for(size_t j = 0; j < n - 1; ++j){
				   cnorm[j] = RNP::TBLAS::Asum( n-j-1, &a[(j+1)+j*lda], 1 );
				}
				cnorm[n-1] = zero;
			 }
		  }
		  // 
		  // Scale the column norms by TSCAL if the maximum element in CNORM is
		  // greater than BIGNUM/2.
		  // 
		  size_t imax = RNP::TBLAS::MaximumIndex( n, cnorm, 1 );
		  tmax = cnorm[imax];
		  if( tmax <= bignum*half ){
			 tscal = one;
		  }else{
			 tscal = half / ( smlnum*tmax );
			 RNP::TBLAS::Scale( n, tscal, cnorm, 1 );
		  }
		  // 
		  // Compute a bound on the computed solution vector to see if the
		  // Level 2 BLAS routine ZTRSV can be used.
		  // 
		  xmax = zero;
		  for(size_t j = 0; j < n; ++j){
			 xmax = max( xmax, absv( x[j] ) );
		  }
		  xbnd = xmax;
		  // 
		  if( notran ){
		  // 
		  // Compute the growth in A * x = b.
		  // 
			 if( upper ){
				jfirst = n;
				jlast = 1;
				jinc = -1;
			 }else{
				jfirst = 1;
				jlast = n;
				jinc = 1;
			 }
			 // 
			 if( tscal != one ){
				grow = zero;
			 }else{
				 // 
				 if( nounit ){
				 // 
				 // A is non-unit triangular.
				 // 
				 // Compute GROW = 1/G(j) and XBND = 1/M(j).
				 // Initially, G(0) = max{x(i), i=1,...,n}.
				 // 
					grow = half / max( xbnd, smlnum );
					xbnd = grow;
					bool breaked_out = false;
					for(int j = jfirst; ((jinc < 0) ? (j >= jlast) : (j <= jlast)); j += jinc){
					// 
					// Exit the loop if the growth factor is too small.
					// 
					   if( grow <= smlnum ){
							breaked_out = true;
							break;
						}
					   // 
					   tjjs = a[(j-1)+(j-1)*lda];
					   tjj = abs1( tjjs );
					   // 
					   if( tjj >= smlnum ){
					   // 
					   // M(j) = G(j-1) / std::abs(A(j,j))
					   // 
						  xbnd = min( xbnd, min( one, tjj )*grow );
					   }else{
					   // 
					   // M(j) could overflow, set XBND to 0.
					   // 
						  xbnd = zero;
					   }
					   // 
					   if( tjj+cnorm[j-1] >= smlnum ){
					   // 
					   // G(j) = G(j-1)*( 1 + CNORM(j) / std::abs(A(j,j)) )
					   // 
						  grow = grow*( tjj / ( tjj+cnorm[j-1] ) );
					   }else{
					   // 
					   // G(j) could overflow, set GROW to 0.
					   // 
						  grow = zero;
					   }
					}
					if(!breaked_out){
						grow = xbnd;
					}
				 }else{
				 // 
				 // A is unit triangular.
				 // 
				 // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
				 // 
					grow = min( one, half / max( xbnd, smlnum ) );
					for(int j = jfirst; ((jinc < 0) ? (j >= jlast) : (j <= jlast)); j += jinc){
					// 
					// Exit the loop if the growth factor is too small.
					// 
					   if( grow <= smlnum ) break;
					   // 
					   // G(j) = G(j-1)*( 1 + CNORM(j) )
					   // 
					   grow = grow*( one / ( one+cnorm[j-1] ) );
					}
				 }
			}
		  }else{ // Compute the growth in A**T * x = b  or  A**H * x = b.
		  // 
			 if( upper ){
				jfirst = 1;
				jlast = n;
				jinc = 1;
			 }else{
				jfirst = n;
				jlast = 1;
				jinc = -1;
			 }
			 // 
			 if( tscal != one ){
				grow = zero;
			 }else{
				 // 
				 if( nounit ){
				 // 
				 // A is non-unit triangular.
				 // 
				 // Compute GROW = 1/G(j) and XBND = 1/M(j).
				 // Initially, M(0) = max{x(i), i=1,...,n}.
				 // 
					grow = half / max( xbnd, smlnum );
					xbnd = grow;
					bool breaked_out = false;
					for(int j = jfirst; ((jinc < 0) ? (j >= jlast) : (j <= jlast)); j += jinc){
					// 
					// Exit the loop if the growth factor is too small.
					// 
					   if( grow <= smlnum ){
						breaked_out = true;
							break;
						}
					   // 
					   // G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
					   // 
					   xj = one + cnorm[j-1];
					   grow = min( grow, xbnd / xj );
					   // 
					   tjjs = a[(j-1)+(j-1)*lda];
					   tjj = abs1( tjjs );
					   // 
					   if( tjj >= smlnum ){
					   // 
					   // M(j) = M(j-1)*( 1 + CNORM(j) ) / std::abs(A(j,j))
					   // 
						  if( xj > tjj ) xbnd = xbnd*( tjj / xj );
					   }else{
					   // 
					   // M(j) could overflow, set XBND to 0.
					   // 
						  xbnd = zero;
					   }
					}
					if(!breaked_out){
						grow = min( grow, xbnd );
					}
				 }else{
				 // 
				 // A is unit triangular.
				 // 
				 // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
				 // 
					grow = min( one, half / max( xbnd, smlnum ) );
					for(int j = jfirst; ((jinc < 0) ? (j >= jlast) : (j <= jlast)); j += jinc){
					// 
					// Exit the loop if the growth factor is too small.
					// 
					   if( grow <= smlnum ) break;
					   // 
					   // G(j) = ( 1 + CNORM(j) )*G(j-1)
					   // 
					   xj = one + cnorm[j-1];
					   grow = grow / xj;
					}
				 }
				}
		  }
		  // 
		  if( ( grow*tscal ) > smlnum ){
		  // 
		  // Use the Level 2 BLAS solve if the reciprocal of the bound on
		  // elements of X is not too small.
		  // 
			 RNP::TBLAS::SolveTrV<uplo,trans,diag>(n, a, lda, x, 1 );
		  }else{
		  // 
		  // Use a Level 1 BLAS solve, scaling intermediate results.
		  // 
			 if( xmax > bignum*half ){
			 // 
			 // Scale X so that its components are less than or equal to
			 // BIGNUM in absolute value.
			 // 
				*scale = ( bignum*half ) / xmax;
				RNP::TBLAS::Scale( n, *scale, x, 1 );
				xmax = bignum;
			 }else{
				xmax = xmax*two;
			 }
			 // 
			 if( notran ){
			 // 
			 // Solve A * x = b
			 // 
				for(int j = jfirst; ((jinc < 0) ? (j >= jlast) : (j <= jlast)); j += jinc){
				// 
				// Compute x(j) = b(j) / A(j,j), scaling x if necessary.
				// 
				   xj = abs1( x[j-1] );
				   if( nounit ){
					  tjjs = a[(j-1)+(j-1)*lda]*tscal;
				   }else{
					  tjjs = tscal;
				   }
				   if(nounit || (!nounit && tscal != one)){
					   tjj = abs1( tjjs );
					   if( tjj > smlnum ){
					   // 
					   // std::abs(A(j,j)) > SMLNUM:
					   // 
						  if( tjj < one ){
							 if( xj > tjj*bignum ){
							 // 
							 // Scale x by 1/b(j).
							 // 
								rec = one / xj;
								RNP::TBLAS::Scale( n, rec, x, 1 );
								*scale *= rec;
								xmax = xmax*rec;
							 }
						  }
						  x[j-1] = RNP::TLASupport::ComplexDivide( x[j-1], tjjs );
						  xj = abs1( x[j-1] );
					   }else if( tjj > zero ){
					   // 
					   // 0 < std::abs(A(j,j)) <= SMLNUM:
					   // 
						  if( xj > tjj*bignum ){
						  // 
						  // Scale x by (1/std::abs(x(j)))*std::abs(A(j,j))*BIGNUM
						  // to avoid overflow when dividing by A(j,j).
						  // 
							 rec = ( tjj*bignum ) / xj;
							 if( cnorm[j-1] > one ){
							 // 
							 // Scale by 1/CNORM(j) to avoid overflow when
							 // multiplying x(j) times column j.
							 // 
								rec = rec / cnorm[j-1];
							 }
							 RNP::TBLAS::Scale( n, rec, x, 1 );
							 *scale *= rec;
							 xmax = xmax*rec;
						  }
						  x[j-1] = RNP::TLASupport::ComplexDivide( x[j-1], tjjs );
						  xj = abs1( x[j-1] );
					   }else{
					   // 
					   // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
					   // scale = 0, and compute a solution to A*x = 0.
					   // 
						  for(size_t i = 1; i <= n; ++i){
							 x[i-1] = zero;
						  }
						  x[j-1] = one;
						  xj = one;
						  *scale = zero;
						  xmax = zero;
					   }
					  }
	// 
	// Scale x if necessary to avoid overflow when adding a
	// multiple of column j of A.
	// 
				   if( xj > one ){
					  rec = one / xj;
					  if( cnorm[j-1] > ( bignum-xmax )*rec ){
					  // 
					  // Scale x by 1/(2*std::abs(x(j))).
					  // 
						 rec = rec*half;
						 RNP::TBLAS::Scale( n, rec, x, 1 );
						 *scale *= rec;
					  }
				   }else if( xj*cnorm[j-1] > ( bignum-xmax ) ){
				   // 
				   // Scale x by 1/2.
				   // 
					  RNP::TBLAS::Scale( n, half, x, 1 );
					  *scale *= half;
				   }
				   // 
				   if( upper ){
					  if( j > 1 ){
					  // 
					  // Compute the update
					  // x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
					  // 
						 RNP::TBLAS::Axpy( j-1, -x[j-1]*tscal, &a[0+(j-1)*lda], 1, x, 1 );
						 size_t i = 1+RNP::TBLAS::MaximumIndex( j-1, x, 1 );
						 xmax = abs1( x[i-1] );
					  }
				   }else{
					  if((size_t)j < n ){
					  // 
					  // Compute the update
					  // x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
					  // 
						 RNP::TBLAS::Axpy( n-j, -x[j-1]*tscal, &a[j+(j-1)*lda], 1, &x[j], 1 );
						 size_t i = j + 1+RNP::TBLAS::MaximumIndex( n-j, &x[j], 1 );
						 xmax = abs1( x[i-1] );
					  }
				   }
				}
				// 
			 }else if( trans == 'T'){
			 // 
			 // Solve A**T * x = b
			 // 
				for(int j = jfirst; ((jinc < 0) ? (j >= jlast) : (j <= jlast)); j += jinc){
				// 
				// Compute x(j) = b(j) - sum A(k,j)*x(k).
				// k<>j
				// 
				   xj = abs1( x[j-1] );
				   uscal = tscal;
				   rec = one / max( xmax, one );
				   if( cnorm[j-1] > ( bignum-xj )*rec ){
				   // 
				   // If x(j) could overflow, scale x by 1/(2*XMAX).
				   // 
					  rec = rec*half;
					  if( nounit ){
						 tjjs = a[(j-1)+(j-1)*lda]*tscal;
					  }else{
						 tjjs = tscal;
					  }
					  tjj = abs1( tjjs );
					  if( tjj > one ){
					  // 
					  // Divide by A(j,j) when scaling x if A(j,j) > 1.
					  // 
						 rec = min( one, rec*tjj );
						 uscal = RNP::TLASupport::ComplexDivide( uscal, tjjs );
					  }
					  if( rec < one ){
						 RNP::TBLAS::Scale( n, rec, x, 1 );
						 *scale *= rec;
						 xmax = xmax*rec;
					  }
				   }
				   // 
				   csumj = zero;
				   if( uscal == std::complex<double>( one ) ){
				   // 
				   // If the scaling needed for A in the dot product is 1,
				   // call ZDOTU to perform the dot product.
				   // 
					  if( upper ){
						 csumj = RNP::TBLAS::Dot( j-1, &a[0+(j-1)*lda], 1, x, 1 );
					  }else if((size_t)j < n ){
						 csumj = RNP::TBLAS::Dot( n-j, &a[j+(j-1)*lda], 1, &x[j], 1 );
					  }
				   }else{
				   // 
				   // Otherwise, use in-line code for the dot product.
				   // 
					  if( upper ){
						 for(int i = 1; i < j; ++i){
							csumj += ( a[(i-1)+(j-1)*lda]*uscal )*x[i-1];
						 }
					  }else if((size_t)j < n ){
						 for(size_t i = j+1; i <= n; ++i){
							csumj += ( a[(i-1)+(j-1)*lda]*uscal )*x[i-1];
						 }
					  }
				   }
				   // 
				   if( uscal == std::complex<double>( tscal ) ){
				   // 
				   // Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
				   // was not used to scale the dotproduct.
				   // 
					  x[j-1] = x[j-1] - csumj;
					  xj = abs1( x[j-1] );
					  if( nounit ){
						 tjjs = a[(j-1)+(j-1)*lda]*tscal;
					  }else{
						 tjjs = tscal;
					  }
					  if(nounit || (!nounit && tscal != one)){
						  // 
						  // Compute x(j) = x(j) / A(j,j), scaling if necessary.
						  // 
						  tjj = abs1( tjjs );
						  if( tjj > smlnum ){
						  // 
						  // std::abs(A(j,j)) > SMLNUM:
						  // 
							 if( tjj < one ){
								if( xj > tjj*bignum ){
								// 
								// Scale X by 1/std::abs(x(j)).
								// 
								   rec = one / xj;
								   RNP::TBLAS::Scale( n, rec, x, 1 );
								   *scale *= rec;
								   xmax = xmax*rec;
								}
							 }
							 x[j-1] = RNP::TLASupport::ComplexDivide( x[j-1], tjjs );
						  }else if( tjj > zero ){
						  // 
						  // 0 < std::abs(A(j,j)) <= SMLNUM:
						  // 
							 if( xj > tjj*bignum ){
							 // 
							 // Scale x by (1/std::abs(x(j)))*std::abs(A(j,j))*BIGNUM.
							 // 
								rec = ( tjj*bignum ) / xj;
								RNP::TBLAS::Scale( n, rec, x, 1 );
								*scale *= rec;
								xmax = xmax*rec;
							 }
							 x[j-1] = RNP::TLASupport::ComplexDivide( x[j-1], tjjs );
						  }else{
						  // 
						  // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
						  // scale = 0 and compute a solution to A**T *x = 0.
						  // 
							 for(size_t i = 0; i < n; ++i){
								x[i] = zero;
							 }
							 x[j-1] = one;
							 *scale = zero;
							 xmax = zero;
						  }
						 }
				   }else{
				   // 
				   // Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
				   // product has already been divided by 1/A(j,j).
				   // 
					  x[j-1] = RNP::TLASupport::ComplexDivide( x[j-1], tjjs ) - csumj;
				   }
				   xmax = max( xmax, abs1( x[j-1] ) );
				}
				// 
			 }else{
			 // 
			 // Solve A**H * x = b
			 // 
				for(int j = jfirst; ((jinc < 0) ? (j >= jlast) : (j <= jlast)); j += jinc){
				// 
				// Compute x(j) = b(j) - sum A(k,j)*x(k).
				// k<>j
				// 
				   xj = abs1( x[j-1] );
				   uscal = tscal;
				   rec = one / max( xmax, one );
				   if( cnorm[j-1] > ( bignum-xj )*rec ){
				   // 
				   // If x(j) could overflow, scale x by 1/(2*XMAX).
				   // 
					  rec *= half;
					  if( nounit ){
						 tjjs = std::conj( a[(j-1)+(j-1)*lda] )*tscal;
					  }else{
						 tjjs = tscal;
					  }
					  tjj = abs1( tjjs );
					  if( tjj > one ){
					  // 
					  // Divide by A(j,j) when scaling x if A(j,j) > 1.
					  // 
						 rec = min( one, rec*tjj );
						 uscal = RNP::TLASupport::ComplexDivide( uscal, tjjs );
					  }
					  if( rec < one ){
						 RNP::TBLAS::Scale( n, rec, x, 1 );
						 *scale *= rec;
						 xmax = xmax*rec;
					  }
				   }
				   // 
				   csumj = zero;
				   if( uscal == std::complex<double>( one ) ){
				   // 
				   // If the scaling needed for A in the dot product is 1,
				   // call ZDOTC to perform the dot product.
				   // 
					  if( upper ){
						 csumj = RNP::TBLAS::ConjugateDot( j-1, &a[0+(j-1)*lda], 1, x, 1 );
					  }else if((size_t)j < n ){
						 csumj = RNP::TBLAS::ConjugateDot( n-j, &a[j+(j-1)*lda], 1, &x[j], 1 );
					  }
				   }else{
				   // 
				   // Otherwise, use in-line code for the dot product.
				   // 
					  if( upper ){
						 for(int i = 1; i < j; ++i){
							csumj += ( std::conj( a[(i-1)+(j-1)*lda] )*uscal )* x[i-1];
						 }
					  }else if((size_t)j < n ){
						 for(size_t i = j + 1; i <= n; ++i){
							csumj += ( std::conj( a[(i-1)+(j-1)*lda] )*uscal )* x[i-1];
						 }
					  }
				   }
				   // 
				   if( uscal == std::complex<double>( tscal ) ){
				   // 
				   // Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
				   // was not used to scale the dotproduct.
				   // 
					  x[j-1] = x[j-1] - csumj;
					  xj = abs1( x[j-1] );
					  if( nounit ){
						 tjjs = std::conj( a[(j-1)+(j-1)*lda] )*tscal;
					  }else{
						 tjjs = tscal;
					  }
					  if(nounit || (!nounit && tscal != one)){
						  // 
						  // Compute x(j) = x(j) / A(j,j), scaling if necessary.
						  // 
						  tjj = abs1( tjjs );
						  if( tjj > smlnum ){
						  // 
						  // std::abs(A(j,j)) > SMLNUM:
						  // 
							 if( tjj < one ){
								if( xj > tjj*bignum ){
								// 
								// Scale X by 1/std::abs(x(j)).
								// 
								   rec = one / xj;
								   RNP::TBLAS::Scale( n, rec, x, 1 );
								   *scale *= rec;
								   xmax = xmax*rec;
								}
							 }
							 x[j-1] = RNP::TLASupport::ComplexDivide( x[j-1], tjjs );
						  }else if( tjj > zero ){
						  // 
						  // 0 < std::abs(A(j,j)) <= SMLNUM:
						  // 
							 if( xj > tjj*bignum ){
							 // 
							 // Scale x by (1/std::abs(x(j)))*std::abs(A(j,j))*BIGNUM.
							 // 
								rec = ( tjj*bignum ) / xj;
								RNP::TBLAS::Scale( n, rec, x, 1 );
								*scale *= rec;
								xmax *= rec;
							 }
							 x[j-1] = RNP::TLASupport::ComplexDivide( x[j-1], tjjs );
						  }else{
						  // 
						  // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
						  // scale = 0 and compute a solution to A**H *x = 0.
						  // 
							 for(size_t i = 0; i < n; ++i){
								x[i] = zero;
							 }
							 x[j-1] = one;
							 *scale = zero;
							 xmax = zero;
						  }
						}
				   }else{
				   // 
				   // Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
				   // product has already been divided by 1/A(j,j).
				   // 
					  x[j-1] = RNP::TLASupport::ComplexDivide( x[j-1], tjjs ) - csumj;
				   }
				   xmax = max( xmax, abs1( x[j-1] ) );
				}
			 }
			 *scale /= tscal;
		  }
		  // 
		  // Scale the column norms by 1/TSCAL for return.
		  // 
		  if( tscal != one ){
			 RNP::TBLAS::Scale( n, one / tscal, cnorm, 1 );
		  }
	}
};




int ztrevc_(char howmny, bool *_select, size_t n, std::complex<double> *_t, size_t ldt, std::complex<double> *_vl, size_t ldvl, std::complex<double> *_vr, size_t ldvr, size_t mm, size_t *m, std::complex<double> *_work, double *rwork)
{
	bool *select = _select;
	std::complex<double> *t = (std::complex<double>*)_t;
	std::complex<double> *vl = (std::complex<double>*)_vl;
	std::complex<double> *vr = (std::complex<double>*)_vr;
	std::complex<double> *work = (std::complex<double>*)_work;

	/* System generated locals */
	size_t t_offset, vl_offset, vr_offset;

	/* Local variables */
	int ii, is;
	double scale;
	double remax;

/*  Purpose */
/*  ======= */

/*  ZTREVC computes some or all of the right and/or left eigenvectors of */
/*  a complex upper triangular matrix T. */
/*  Matrices of this type are produced by the Schur factorization of */
/*  a complex general matrix:  A = Q*T*Q**H, as computed by ZHSEQR. */

/*  The right eigenvector x and the left eigenvector y of T corresponding */
/*  to an eigenvalue w are defined by: */

/*               T*x = w*x,     (y**H)*T = w*(y**H) */

/*  where y**H denotes the conjugate transpose of the vector y. */
/*  The eigenvalues are not input to this routine, but are read directly */
/*  from the diagonal of T. */

/*  This routine returns the matrices X and/or Y of right and left */
/*  eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an */
/*  input matrix.  If Q is the unitary factor that reduces a matrix A to */
/*  Schur form T, then Q*X and Q*Y are the matrices of right and left */
/*  eigenvectors of A. */

/*  Arguments */
/*  ========= */

/*  HOWMNY  (input) CHARACTER*1 */
/*          = 'A':  compute all right and/or left eigenvectors; */
/*          = 'B':  compute all right and/or left eigenvectors, */
/*                  backtransformed using the matrices supplied in */
/*                  VR and/or VL; */
/*          = 'S':  compute selected right and/or left eigenvectors, */
/*                  as indicated by the logical array SELECT. */

/*  SELECT  (input) LOGICAL array, dimension (N) */
/*          If HOWMNY = 'S', SELECT specifies the eigenvectors to be */
/*          computed. */
/*          The eigenvector corresponding to the j-th eigenvalue is */
/*          computed if SELECT(j) = .TRUE.. */
/*          Not referenced if HOWMNY = 'A' or 'B'. */

/*  N       (input) INTEGER */
/*          The order of the matrix T. N >= 0. */

/*  T       (input/output) COMPLEX*16 array, dimension (LDT,N) */
/*          The upper triangular matrix T.  T is modified, but restored */
/*          on exit. */

/*  LDT     (input) INTEGER */
/*          The leading dimension of the array T. LDT >= max(1,N). */

/*  VL      (input/output) COMPLEX*16 array, dimension (LDVL,MM) */
/*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must */
/*          contain an N-by-N matrix Q (usually the unitary matrix Q of */
/*          Schur vectors returned by ZHSEQR). */
/*          On exit, if SIDE = 'L' or 'B', VL contains: */
/*          if HOWMNY = 'A', the matrix Y of left eigenvectors of T; */
/*          if HOWMNY = 'B', the matrix Q*Y; */
/*          if HOWMNY = 'S', the left eigenvectors of T specified by */
/*                           SELECT, stored consecutively in the columns */
/*                           of VL, in the same order as their */
/*                           eigenvalues. */
/*          Not referenced if SIDE = 'R'. */

/*  LDVL    (input) INTEGER */
/*          The leading dimension of the array VL.  LDVL >= 1, and if */
/*          SIDE = 'L' or 'B', LDVL >= N. */

/*  VR      (input/output) COMPLEX*16 array, dimension (LDVR,MM) */
/*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must */
/*          contain an N-by-N matrix Q (usually the unitary matrix Q of */
/*          Schur vectors returned by ZHSEQR). */
/*          On exit, if SIDE = 'R' or 'B', VR contains: */
/*          if HOWMNY = 'A', the matrix X of right eigenvectors of T; */
/*          if HOWMNY = 'B', the matrix Q*X; */
/*          if HOWMNY = 'S', the right eigenvectors of T specified by */
/*                           SELECT, stored consecutively in the columns */
/*                           of VR, in the same order as their */
/*                           eigenvalues. */
/*          Not referenced if SIDE = 'L'. */

/*  LDVR    (input) INTEGER */
/*          The leading dimension of the array VR.  LDVR >= 1, and if */
/*          SIDE = 'R' or 'B'; LDVR >= N. */

/*  MM      (input) INTEGER */
/*          The number of columns in the arrays VL and/or VR. MM >= M. */

/*  M       (output) INTEGER */
/*          The number of columns in the arrays VL and/or VR actually */
/*          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M */
/*          is set to N.  Each selected eigenvector occupies one */
/*          column. */

/*  WORK    (workspace) COMPLEX*16 array, dimension (2*N) */

/*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N) */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

/*  Further Details */
/*  =============== */

/*  The algorithm used in this program is basically backward (forward) */
/*  substitution, with scaling to make the the code robust against */
/*  possible overflow. */

/*  Each eigenvector is normalized so that the element of largest */
/*  magnitude has magnitude 1; here the magnitude of a complex number */
/*  (x,y) is taken to be |x| + |y|. */

	using namespace std;

	const bool rightv = (NULL != vr);
	const bool leftv = (NULL != vl);
	//const bool bothv = rightv && leftv;
	
	/* Parameter adjustments */
	--select;
	t_offset = 1 + ldt;
	t -= t_offset;
	vl_offset = 1 + ldvl;
	vl -= vl_offset;
	vr_offset = 1 + ldvr;
	vr -= vr_offset;
	--work;
	--rwork;

	/* Function Body */

	//const bool allv = (howmny == 'A');
	const bool over = (howmny == 'B');
	const bool somev = (howmny == 'S');

/*     Set M to the number of columns required to store the selected */
/*     eigenvectors. */

	if (somev) {
		*m = 0;
		for (size_t j = 0; j < n; ++j) {
			if (_select[j]) {
				++(*m);
			}
		}
	} else {
		*m = n;
	}

	if (n == 0) {
		return 0;
	}

/*     Set the constants to control overflow. */

	const double unfl = std::numeric_limits<double>::min();
	//const double ovfl = 1. / unfl;
	const double ulp = 2*std::numeric_limits<double>::epsilon();
	const double smlnum = unfl * (n / ulp);

/*     Store the diagonal elements of T in working array WORK. */

	for (size_t i = 0; i < n; ++i) {
		_work[i+n] = _t[i+i*ldt];
	}

/*     Compute 1-norm of each column of strictly upper triangular */
/*     part of T to control overflow in triangular solver. */

	rwork[1] = 0.;
	for (size_t j = 1; j < n; ++j) {
		rwork[j] = RNP::TBLAS::Asum(j, &_t[0+j*ldt], 1);
	}

	if (rightv) { // Compute right eigenvectors.
		is = *m;
		size_t ki = n+1; while(ki --> 1){
			if (somev && ! select[ki]) {
				continue;
			}
			double smin = max(ulp * abs1(t[ki + ki * ldt]),smlnum);

			_work[0] = 1.;

			// Form right-hand side.
			for (size_t k = 1; k < ki; ++k) {
				work[k] = -t[k+ki*ldt];
			}

			// Solve the triangular system:
			//    (T(1:KI-1,1:KI-1) - T(KI,KI))*X = SCALE*WORK.
			for (size_t k = 1; k < ki; ++k) {
				t[k+k*ldt] -= t[ki+ki*ldt];
				if (abs1(t[k + k * ldt]) < smin) {
					t[k+k*ldt] = smin;
				}
			}

			if (ki > 1) {
				SolveTrV_Scaled<'U','N','N','Y'>(ki-1, _t, ldt, _work, &scale, &rwork[1]);
				work[ki] = scale;
			}

			// Copy the vector x or Q*x to VR and normalize.
			if (! over) {
				RNP::TBLAS::Copy(ki, _work, 1, &vr[is * ldvr + 1], 1);

				ii = 1+RNP::TBLAS::MaximumIndex(ki, &vr[is * ldvr + 1], 1);
				remax = 1. / abs1(vr[ii + is * ldvr]);
				RNP::TBLAS::Scale(ki, remax, &vr[is * ldvr + 1], 1);

				for (size_t k = ki + 1; k <= n; ++k) {
					vr[k+is*ldvr] = 0.;
				}
			} else {
				if (ki > 1) {
					RNP::TBLAS::MultMV<'N'>(n, ki - 1, 1., _vr, ldvr, _work, 1, scale, &vr[ki * ldvr + 1], 1);
				}

				ii = 1+RNP::TBLAS::MaximumIndex(n, &vr[ki * ldvr + 1], 1);
				remax = 1. / abs1(vr[ii + ki * ldvr]);
				RNP::TBLAS::Scale(n, remax, &vr[ki * ldvr + 1], 1);
			}

			// Set back the original diagonal elements of T.
			for (size_t k = 1; k < ki; ++k) {
				t[k+k*ldt] = work[k+n];
			}

			--is;
		}
	}

	if (leftv) { // Compute left eigenvectors.

		is = 1;
		for (size_t ki = 1; ki <= n; ++ki) {
			if (somev && ! select[ki]) {
				continue;
			}
			double smin = max(ulp * abs1(t[ki + ki * ldt]), smlnum);

			_work[n-1] = 1.;

			// Form right-hand side.
			for (size_t k = ki + 1; k <= n; ++k) {
				work[k] = -std::conj(t[ki + k * ldt]);
			}

			// Solve the triangular system:
			//    (T(KI+1:N,KI+1:N) - T(KI,KI))'*X = SCALE*WORK.
			for (size_t k = ki + 1; k <= n; ++k) {
				t[k+k*ldt] -= t[ki+ki*ldt];
				if (abs1(t[k + k * ldt]) < smin) {
					t[k+k*ldt] = smin;
				}
			}

			if (ki < n) {
				SolveTrV_Scaled<'U','C','N','Y'>(n-ki, &t[ki + 1 + (ki + 1) * ldt], ldt, &work[ki + 1], &scale, &rwork[1]);
				work[ki] = scale;
			}

			// Copy the vector x or Q*x to VL and normalize.
			if (! over) {
				RNP::TBLAS::Copy(n - ki + 1, &work[ki], 1, &vl[ki + is * ldvl], 1);

				ii = RNP::TBLAS::MaximumIndex(n - ki + 1, &vl[ki + is * ldvl], 1) + ki;
				remax = 1. / abs1(vl[ii + is * ldvl]);
				RNP::TBLAS::Scale(n - ki + 1, remax, &vl[ki + is * ldvl], 1);

				for (size_t k = 1; k < ki; ++k) {
					vl[k+is*ldvl] = 0.;
				}
			} else {
				if (ki < n) {
					RNP::TBLAS::MultMV<'N'>(n, n - ki, 1., &vl[(ki + 1) * ldvl + 1], ldvl, &_work[ki], 1, scale, &vl[ki * ldvl + 1], 1);
				}

				ii = 1+RNP::TBLAS::MaximumIndex(n, &vl[ki * ldvl + 1], 1);
				remax = 1. / abs1(vl[ii + ki * ldvl]);
				RNP::TBLAS::Scale(n, remax, &vl[ki * ldvl + 1], 1);
			}

			// Set back the original diagonal elements of T.
			for (size_t k = ki + 1; k <= n; ++k) {
				t[k+k*ldt] = work[k+n];
			}

			++is;
		}
	}
	return 0;
} /* ztrevc_ */



static void zlaqr5(bool wantt, bool wantz, int kacc22, int n, int ktop, int kbot, int nshfts, std::complex<double> *s, std::complex<double> *h, int ldh, int iloz, int ihiz, std::complex<double> *z, int ldz, std::complex<double> *v, int ldv, std::complex<double> *u, int ldu, int nv, std::complex<double> *wv, int ldwv, int nh, std::complex<double> *wh, int ldwh){
	using namespace std;
	
    // 
    // This auxiliary subroutine called by ZLAQR0 performs a
    // single small-bulge multi-shift QR sweep.
    // 
    // WANTT  (input) bool scalar
    // WANTT = true if the triangular Schur factor
    // is being computed.  WANTT is set to false otherwise.
    // 
    // WANTZ  (input) bool scalar
    // WANTZ = true if the unitary Schur factor is being
    // computed.  WANTZ is set to false otherwise.
    // 
    // KACC22 (input) int with value 0, 1, or 2.
    // Specifies the computation mode of far-from-diagonal
    // orthogonal updates.
    // = 0: ZLAQR5 does not accumulate reflections and does not
    // use matrix-matrix multiply to update far-from-diagonal
    // matrix entries.
    // = 1: ZLAQR5 accumulates reflections and uses matrix-matrix
    // multiply to update the far-from-diagonal matrix entries.
    // = 2: ZLAQR5 accumulates reflections, uses matrix-matrix
    // multiply to update the far-from-diagonal matrix entries,
    // and takes advantage of 2-by-2 block structure during
    // matrix multiplies.
    // 
    // N      (input) int scalar
    // N is the order of the Hessenberg matrix H upon which this
    // subroutine operates.
    // 
    // KTOP   (input) int scalar
    // KBOT   (input) int scalar
    // These are the first and last rows and columns of an
    // isolated diagonal block upon which the QR sweep is to be
    // applied. It is assumed without a check that
    // either KTOP = 1  or   H(KTOP,KTOP-1) = 0
    // and
    // either KBOT = N  or   H(KBOT+1,KBOT) = 0.
    // 
    // NSHFTS (input) int scalar
    // NSHFTS gives the number of simultaneous shifts.  NSHFTS
    // must be positive and even.
    // 
    // S      (input/output) std::complex<double> array of size (NSHFTS)
    // S contains the shifts of origin that define the multi-
    // shift QR sweep.  On output S may be reordered.
    // 
    // H      (input/output) std::complex<double> array of size (LDH,N)
    // On input H contains a Hessenberg matrix.  On output a
    // multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied
    // to the isolated diagonal block in rows and columns KTOP
    // through KBOT.
    // 
    // LDH    (input) int scalar
    // LDH is the leading dimension of H just as declared in the
    // calling procedure.  LDH >= MAX(1,N).
    // 
    // ILOZ   (input) int
    // IHIZ   (input) int
    // Specify the rows of Z to which transformations must be
    // applied if WANTZ is true. 1  <=  ILOZ  <=  IHIZ  <=  N
    // 
    // Z      (input/output) std::complex<double> array of size (LDZ,IHI)
    // If WANTZ = true, then the QR Sweep unitary
    // similarity transformation is accumulated into
    // Z(ILOZ:IHIZ,ILO:IHI) from the right.
    // If WANTZ = false, then Z is unreferenced.
    // 
    // LDZ    (input) int scalar
    // LDA is the leading dimension of Z just as declared in
    // the calling procedure. LDZ >= N.
    // 
    // V      (workspace) std::complex<double> array of size (LDV,NSHFTS/2)
    // 
    // LDV    (input) int scalar
    // LDV is the leading dimension of V as declared in the
    // calling procedure.  LDV >= 3.
    // 
    // U      (workspace) std::complex<double> array of size
    // (LDU,3*NSHFTS-3)
    // 
    // LDU    (input) int scalar
    // LDU is the leading dimension of U just as declared in the
    // in the calling subroutine.  LDU >= 3*NSHFTS-3.
    // 
    // NH     (input) int scalar
    // NH is the number of columns in array WH available for
    // workspace. NH >= 1.
    // 
    // WH     (workspace) std::complex<double> array of size (LDWH,NH)
    // 
    // LDWH   (input) int scalar
    // Leading dimension of WH just as declared in the
    // calling procedure.  LDWH >= 3*NSHFTS-3.
    // 
    // NV     (input) int scalar
    // NV is the number of rows in WV agailable for workspace.
    // NV >= 1.
    // 
    // WV     (workspace) std::complex<double> array of size
    // (LDWV,3*NSHFTS-3)
    // 
    // LDWV   (input) int scalar
    // LDWV is the leading dimension of WV as declared in the
    // in the calling subroutine.  LDWV >= NV.
    // 
    //============================================================
    // Based on contributions by
    // Karen Braman and Ralph Byers, Department of Mathematics,
    // University of Kansas, USA
    // 
    //============================================================
    // Reference:
    // 
    // K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
    // Algorithm Part I: Maintaining Well Focused Shifts, and
    // Level 3 Performance, SIAM Journal of Matrix Analysis,
    // volume 23, pages 929--947, 2002.
    // 
    //============================================================
    // .. Parameters ..
	const std::complex<double> zero = std::complex<double>(0.0e0, 0.0e0);
	const std::complex<double> one = std::complex<double>(1.0e0, 0.0e0);
	// parameter          ( zero = ( 0.0e0, 0.0e0 ), one = ( 1.0e0, 0.0e0 ) )
	const double rzero = 0.0e0;
	const double rone = 1.0e0 ;
	// parameter          ( rzero = 0.0e0, rone = 1.0e0 )
	// ..
	// .. Local Scalars ..
	std::complex<double> alpha;
	std::complex<double> beta;
	std::complex<double> cdum;
	std::complex<double> refsum;
	double h11;
	double h12;
	double h21;
	double h22;
	double safmax;
	double safmin;
	double scl;
	double smlnum;
	double tst1;
	double tst2;
	double ulp;
	int i2;
	int i4;
	int incol;
	int j;
	int j2;
	int j4;
	int jbot;
	int jcol;
	int jlen;
	int jrow;
	int jtop;
	int k;
	int k1;
	int kdu;
	int kms;
	int knz;
	int krcol;
	int kzs;
	int m;
	int m22;
	int mbot;
	int mend;
	int mstart;
	int mtop;
	int nbmps;
	int ndcol;
	int ns;
	int nu;
	bool accum;
	bool blk22;
	bool bmp22;

	std::complex<double> vt[3];

	// If there are no shifts, then there is nothing to do.
	if( nshfts < 2 ) return;

	// If the active block is empty or 1-by-1, then there
	// .    is nothing to do.
	if( ktop >= kbot ) return;
	
	// NSHFTS is supposed to be even, but if it is odd,
	// .    then simply reduce it by one. 
	ns = nshfts - nshfts%2;

	// Machine constants for deflation
	safmin = std::numeric_limits<double>::min();
	safmax = rone / safmin;
	ulp = 2*std::numeric_limits<double>::epsilon();
	smlnum = safmin*( double(n) / ulp );

	// Use accumulated reflections to update far-from-diagonal
	// .    entries ?
	accum = ( kacc22 == 1 )  ||  ( kacc22 == 2 );

	// If so, exploit the 2-by-2 block structure?
	blk22 = ( ns > 2 )  &&  ( kacc22 == 2 );

	// clear trash
	if( ktop+2 <= kbot ) h[(( ktop+2)-1)+(ktop-1)*ldh] = zero;

	// NBMPS = number of 2-shift bulges in the chain
	nbmps = ns / 2;

	// KDU = width of slab
	kdu = 6*nbmps - 3;

	// Create and chase chains of NBMPS bulges
	for(incol = 3*( 1-nbmps ) + ktop - 1; ((3*nbmps - 2 < 0) ? (incol >= kbot - 2) : (incol <= kbot - 2)); incol += 3*nbmps - 2){
		ndcol = incol + kdu;
		if( accum ) RNP::TBLAS::SetMatrix<'A'>(kdu, kdu, 0., 1., u, ldu );

		// Near-the-diagonal bulge chase.  The following loop
		// .    performs the near-the-diagonal part of a small bulge
		// .    multi-shift QR sweep.  Each 6*NBMPS-2 column diagonal
		// .    chunk extends from column INCOL to column NDCOL
		// .    (including both column INCOL and column NDCOL). The
		// .    following loop chases a 3*NBMPS column long chain of
		// .    NBMPS bulges 3*NBMPS-2 columns to the right.  (INCOL
		// .    may be less than KTOP and and NDCOL may be greater than
		// .    KBOT indicating phantom columns from which to chase
		// .    bulges before they are actually introduced or to which
		// .    to chase bulges beyond column KBOT.) 
		for(krcol = incol; krcol <= min( incol+3*nbmps-3, kbot-2 ); ++krcol){
			// Bulges number MTOP to MBOT are active double implicit
			// .    shift bulges.  There may or may not also be small
			// .    2-by-2 bulge, if there is room.  The inactive bulges
			// .    (if any) must wait until the active bulges have moved
			// .    down the diagonal to make room.  The phantom matrix
			// .    paradigm described above helps keep track. 
			mtop = max( (int)1, ( ( ktop-1 )-krcol+2 ) / 3+1 );
			mbot = min( nbmps, ( kbot-krcol ) / 3 );
			m22 = mbot + 1;
			bmp22 = ( mbot < nbmps )  &&  ( krcol+3*( m22-1 ) ) ==  ( kbot-2 );

			// Generate reflections to chase the chain right
			// .    one column.  (The minimum value of K is KTOP-1.)
			for(m = mtop; m <= mbot; ++m){
				k = krcol + 3*( m-1 );
				if( k == ktop-1 ){
					zlaqr1_(3, &h[(ktop-1)+(ktop-1)*ldh], ldh, s[( 2*m-1 )-1], s[( 2*m )-1], &v[0+(m-1)*ldv] );
					alpha = v[0+(m-1)*ldv];
					RNP::TLASupport::GenerateElementaryReflector(3, &alpha, &v[1+(m-1)*ldv], 1, &v[0+(m-1)*ldv] );
				}else{
					beta = h[k+(k-1)*ldh];
					v[1+(m-1)*ldv] = h[(( k+2)-1)+(k-1)*ldh];
					v[(3-1)+(m-1)*ldv] = h[(( k+3)-1)+(k-1)*ldh];
					RNP::TLASupport::GenerateElementaryReflector( 3, &beta, &v[1+(m-1)*ldv], 1, &v[0+(m-1)*ldv] );
					// 
					// A Bulge may collapse because of vigilant
					// .    deflation or destructive underflow.  In the
					// .    underflow case, try the two-small-subdiagonals
					// .    trick to try to reinflate the bulge. 
					// 
					if( h[(( k+3)-1)+(k-1)*ldh] != zero  ||  h[(( k+3)-1)+(( k+1 )-1)*ldh] !=  zero  ||  h[(( k+3)-1)+(( k+2 )-1)*ldh] == zero ){
						// 
						// Typical case: not collapsed (yet).
						// 
						h[k+(k-1)*ldh] = beta;
						h[(( k+2)-1)+(k-1)*ldh] = zero;
						h[(( k+3)-1)+(k-1)*ldh] = zero;
					}else{
						// 
						// Atypical case: collapsed.  Attempt to
						// .    reintroduce ignoring H(K+1,K) and H(K+2,K).
						// .    If the fill resulting from the new
						// .    reflector is too large, then abandon it.
						// .    Otherwise, use the new one.
						// 
						zlaqr1_(3, &h[k+(( k+1 )-1)*ldh], ldh, s[( 2*m-1 )-1], s[( 2*m )-1], vt );
						alpha = vt[1-1];
						RNP::TLASupport::GenerateElementaryReflector( 3, &alpha, &vt[2-1], 1, &vt[1-1] );
						refsum = std::conj( vt[1-1] )* ( h[k+(k-1)*ldh]+std::conj( vt[2-1] )* h[(( k+2)-1)+(k-1)*ldh] );
						// 
						if( abs1( h[(( k+2)-1)+(k-1)*ldh]-refsum*vt[2-1] )+ abs1( refsum*vt[3-1] ) > ulp* ( abs1( h[(k-1)+(k-1)*ldh] )+abs1( h[k+(( k+1 )-1)*ldh] )+abs1( h[(( k+2)-1)+(( k+2 )-1)*ldh] ) ) ){
							// 
							// Starting a new bulge here would
							// .    create non-negligible fill.  Use
							// .    the old one with trepidation.
							// 
							h[k+(k-1)*ldh] = beta;
							h[(( k+2)-1)+(k-1)*ldh] = zero;
							h[(( k+3)-1)+(k-1)*ldh] = zero;
						}else{
							// 
							// Stating a new bulge here would
							// .    create only negligible fill.
							// .    Replace the old reflector with
							// .    the new one.
							// 
							h[k+(k-1)*ldh] = h[k+(k-1)*ldh] - refsum;
							h[(( k+2)-1)+(k-1)*ldh] = zero;
							h[(( k+3)-1)+(k-1)*ldh] = zero;
							v[0+(m-1)*ldv] = vt[1-1];
							v[1+(m-1)*ldv] = vt[2-1];
							v[(3-1)+(m-1)*ldv] = vt[3-1];
						}
					}
				}
			}

			// Generate a 2-by-2 reflection, if needed.
			k = krcol + 3*( m22-1 );
			if( bmp22 ){
				if( k == ktop-1 ){
					zlaqr1_(2, &h[k+(( k+1 )-1)*ldh], ldh, s[( 2*m22-1 )-1], s[( 2*m22 )-1], &v[0+(m22-1)*ldv] );
					beta = v[0+(m22-1)*ldv];
					RNP::TLASupport::GenerateElementaryReflector( 2, &beta, &v[1+(m22-1)*ldv], 1, &v[0+(m22-1)*ldv] );
				}else{
					beta = h[k+(k-1)*ldh];
					v[1+(m22-1)*ldv] = h[(( k+2)-1)+(k-1)*ldh];
					RNP::TLASupport::GenerateElementaryReflector( 2, &beta, &v[1+(m22-1)*ldv], 1, &v[0+(m22-1)*ldv] );
					h[k+(k-1)*ldh] = beta;
					h[(( k+2)-1)+(k-1)*ldh] = zero;
				}
			}

			// Multiply H by reflections from the left
			if( accum ){
				jbot = min( ndcol, kbot );
			}else if( wantt ){
				jbot = n;
			}else{
				jbot = kbot;
			}
			for(j = max( ktop, krcol ); j <= jbot; ++j){
				mend = min( mbot, ( j-krcol+2 ) / 3 );
				for(m = mtop; m <= mend; ++m){
					k = krcol + 3*( m-1 );
					refsum = std::conj( v[0+(m-1)*ldv] )* ( h[k+(j-1)*ldh]+std::conj( v[1+(m-1)*ldv] )* h[(( k+2)-1)+(j-1)*ldh]+std::conj( v[(3-1)+(m-1)*ldv] )*h[(( k+3)-1)+(j-1)*ldh] );
					h[k+(j-1)*ldh] = h[k+(j-1)*ldh] - refsum;
					h[(( k+2)-1)+(j-1)*ldh] = h[(( k+2)-1)+(j-1)*ldh] - refsum*v[1+(m-1)*ldv];
					h[(( k+3)-1)+(j-1)*ldh] = h[(( k+3)-1)+(j-1)*ldh] - refsum*v[(3-1)+(m-1)*ldv];
				}
			}
			if( bmp22 ){
				k = krcol + 3*( m22-1 );
				for(j = max( k+1, ktop ); j <= jbot; ++j){
					refsum = std::conj( v[0+(m22-1)*ldv] )* ( h[k+(j-1)*ldh]+std::conj( v[1+(m22-1)*ldv] )* h[(( k+2)-1)+(j-1)*ldh] );
					h[k+(j-1)*ldh] = h[k+(j-1)*ldh] - refsum;
					h[(( k+2)-1)+(j-1)*ldh] = h[(( k+2)-1)+(j-1)*ldh] - refsum*v[1+(m22-1)*ldv];
				}
			}

			// Multiply H by reflections from the right.
			// .    Delay filling in the last row until the
			// .    vigilant deflation check is complete.
			if( accum ){
				jtop = max( ktop, incol );
			}else if( wantt ){
				jtop = 1;
			}else{
				jtop = ktop;
			}
			for(m = mtop; m <= mbot; ++m){
				if( v[0+(m-1)*ldv] != zero ){
					k = krcol + 3*( m-1 );
					for(j = jtop; j <= min( kbot, k+3 ); ++j){
						refsum = v[0+(m-1)*ldv]*( h[(j-1)+(( k+1 )-1)*ldh]+v[1+(m-1)*ldv]* h[(j-1)+(( k+2 )-1)*ldh]+v[(3-1)+(m-1)*ldv]*h[(j-1)+(( k+3 )-1)*ldh] );
						h[(j-1)+(( k+1 )-1)*ldh] = h[(j-1)+(( k+1 )-1)*ldh] - refsum;
						h[(j-1)+(( k+2 )-1)*ldh] = h[(j-1)+(( k+2 )-1)*ldh] - refsum*std::conj( v[1+(m-1)*ldv] );
						h[(j-1)+(( k+3 )-1)*ldh] = h[(j-1)+(( k+3 )-1)*ldh] - refsum*std::conj( v[(3-1)+(m-1)*ldv] );
					}

					if( accum ){
						// Accumulate U. (If necessary, update Z later
						// .    with with an efficient matrix-matrix
						// .    multiply.)
						kms = k - incol;
						for(j = max( (int)1, ktop-incol ); j <= kdu; ++j){
							refsum = v[0+(m-1)*ldv]*( u[(j-1)+(( kms+1 )-1)*ldu]+v[1+(m-1)*ldv]* u[(j-1)+(( kms+2 )-1)*ldu]+v[(3-1)+(m-1)*ldv]*u[(j-1)+(( kms+3 )-1)*ldu] );
							u[(j-1)+(( kms+1 )-1)*ldu] = u[(j-1)+(( kms+1 )-1)*ldu] - refsum;
							u[(j-1)+(( kms+2 )-1)*ldu] = u[(j-1)+(( kms+2 )-1)*ldu] - refsum*std::conj( v[1+(m-1)*ldv] );
							u[(j-1)+(( kms+3 )-1)*ldu] = u[(j-1)+(( kms+3 )-1)*ldu] - refsum*std::conj( v[(3-1)+(m-1)*ldv] );
						}
					}else if( wantz ){
						// U is not accumulated, so update Z
						// .    now by multiplying by reflections
						// .    from the right.
						for(j = iloz; j <= ihiz; ++j){
							refsum = v[0+(m-1)*ldv]*( z[(j-1)+(( k+1 )-1)*ldz]+v[1+(m-1)*ldv]* z[(j-1)+(( k+2 )-1)*ldz]+v[(3-1)+(m-1)*ldv]*z[(j-1)+(( k+3 )-1)*ldz] );
							z[(j-1)+(( k+1 )-1)*ldz] = z[(j-1)+(( k+1 )-1)*ldz] - refsum;
							z[(j-1)+(( k+2 )-1)*ldz] = z[(j-1)+(( k+2 )-1)*ldz] - refsum*std::conj( v[1+(m-1)*ldv] );
							z[(j-1)+(( k+3 )-1)*ldz] = z[(j-1)+(( k+3 )-1)*ldz] - refsum*std::conj( v[(3-1)+(m-1)*ldv] );
						}
					}
				}
			}

			// Special case: 2-by-2 reflection (if needed)
			k = krcol + 3*( m22-1 );
			if( bmp22  &&  ( v[0+(m22-1)*ldv] != zero ) ){
				for(j = jtop; j <= min( kbot, k+3 ); ++j){
					refsum = v[0+(m22-1)*ldv]*( h[(j-1)+(( k+1 )-1)*ldh]+v[1+(m22-1)*ldv]* h[(j-1)+(( k+2 )-1)*ldh] );
					h[(j-1)+(( k+1 )-1)*ldh] = h[(j-1)+(( k+1 )-1)*ldh] - refsum;
					h[(j-1)+(( k+2 )-1)*ldh] = h[(j-1)+(( k+2 )-1)*ldh] - refsum*std::conj( v[1+(m22-1)*ldv] );
				}

				if( accum ){
					kms = k - incol;
					for(j = max( (int)1, ktop-incol ); j <= kdu; ++j){
						refsum = v[0+(m22-1)*ldv]*( u[(j-1)+(( kms+1 )-1)*ldu]+v[1+(m22-1)*ldv]* u[(j-1)+(( kms+2 )-1)*ldu] );
						u[(j-1)+(( kms+1 )-1)*ldu] = u[(j-1)+(( kms+1 )-1)*ldu] - refsum;
						u[(j-1)+(( kms+2 )-1)*ldu] = u[(j-1)+(( kms+2 )-1)*ldu] - refsum*std::conj( v[1+(m22-1)*ldv] );
					}
				}else if( wantz ){
					for(j = iloz; j <= ihiz; ++j){
						refsum = v[0+(m22-1)*ldv]*( z[(j-1)+(( k+1 )-1)*ldz]+v[1+(m22-1)*ldv]* z[(j-1)+(( k+2 )-1)*ldz] );
						z[(j-1)+(( k+1 )-1)*ldz] = z[(j-1)+(( k+1 )-1)*ldz] - refsum;
						z[(j-1)+(( k+2 )-1)*ldz] = z[(j-1)+(( k+2 )-1)*ldz] - refsum*std::conj( v[1+(m22-1)*ldv] );
					}
				}
			}

			// Vigilant deflation check
			mstart = mtop;
			if( krcol+3*( mstart-1 ) < ktop ) mstart = mstart + 1;
			mend = mbot;
			if( bmp22 ) mend = mend + 1;
			if( krcol == kbot-2 ) mend = mend + 1;
			for(m = mstart; m <= mend; ++m){
				k = min( kbot-1, krcol+3*( m-1 ) );

				// The following convergence test requires that
				// .    the tradition small-compared-to-nearby-diagonals
				// .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
				// .    criteria both be satisfied.  The latter improves
				// .    accuracy in some examples. Falling back on an
				// .    alternate convergence criterion when TST1 or TST2
				// .    is zero (as done here) is traditional but probably
				// .    unnecessary.
				if( h[k+(k-1)*ldh] != zero ){
					tst1 = abs1( h[(k-1)+(k-1)*ldh] ) + abs1( h[k+(( k+1 )-1)*ldh] );
					if( tst1 == rzero ){
						if( k >= ktop+1 ) tst1 = tst1 + abs1( h[(k-1)+(( k-1 )-1)*ldh] );
						if( k >= ktop+2 ) tst1 = tst1 + abs1( h[(k-1)+(( k-2 )-1)*ldh] );
						if( k >= ktop+3 ) tst1 = tst1 + abs1( h[(k-1)+(( k-3 )-1)*ldh] );
						if( k <= kbot-2 ) tst1 = tst1 + abs1( h[(( k+2)-1)+(( k+1 )-1)*ldh] );
						if( k <= kbot-3 ) tst1 = tst1 + abs1( h[(( k+3)-1)+(( k+1 )-1)*ldh] );
						if( k <= kbot-4 ) tst1 = tst1 + abs1( h[(( k+4)-1)+(( k+1 )-1)*ldh] );
					}
					if( abs1( h[k+(k-1)*ldh] ) <= max( smlnum, ulp*tst1 ) ){
						h12 = max( abs1( h[k+(k-1)*ldh] ), abs1( h[(k-1)+(( k+1 )-1)*ldh] ) );
						h21 = min( abs1( h[k+(k-1)*ldh] ), abs1( h[(k-1)+(( k+1 )-1)*ldh] ) );
						h11 = max( abs1( h[k+(( k+1 )-1)*ldh] ), abs1( h[(k-1)+(k-1)*ldh]-h[k+(( k+1 )-1)*ldh] ) );
						h22 = min( abs1( h[k+(( k+1 )-1)*ldh] ), abs1( h[(k-1)+(k-1)*ldh]-h[k+(( k+1 )-1)*ldh] ) );
						scl = h11 + h12;
						tst2 = h22*( h11 / scl );
						// 
						if( tst2 == rzero  ||  h21*( h12 / scl ) <=  max( smlnum, ulp*tst2 ) )h[k+(k-1)*ldh] = zero;
					}
				}
			}

			// Fill in the last row of each bulge.
			mend = min( nbmps, ( kbot-krcol-1 ) / 3 );
			for(m = mtop; m <= mend; ++m){
				k = krcol + 3*( m-1 );
				refsum = v[0+(m-1)*ldv]*v[(3-1)+(m-1)*ldv]*h[(( k+4)-1)+(( k+3 )-1)*ldh];
				h[(( k+4)-1)+(( k+1 )-1)*ldh] = -refsum;
				h[(( k+4)-1)+(( k+2 )-1)*ldh] = -refsum*std::conj( v[1+(m-1)*ldv] );
				h[(( k+4)-1)+(( k+3 )-1)*ldh] = h[(( k+4)-1)+(( k+3 )-1)*ldh] - refsum*std::conj( v[(3-1)+(m-1)*ldv] );
			}

			// End of near-the-diagonal bulge chase.
		}

		// Use U (if accumulated) to update far-from-diagonal
		// .    entries in H.  If required, use U to update Z as
		// .    well.
		if( accum ){
			if( wantt ){
				jtop = 1;
				jbot = n;
			}else{
				jtop = ktop;
				jbot = kbot;
			}
			if( (  !blk22 )  ||  ( incol < ktop )  ||  ( ndcol > kbot )  ||  ( ns <= 2 ) ){
				// Updates not exploiting the 2-by-2 block
				// .    structure of U.  K1 and NU keep track of
				// .    the location and size of U in the special
				// .    cases of introducing bulges and chasing
				// .    bulges off the bottom.  In these special
				// .    cases and in case the number of shifts
				// .    is NS = 2, there is no 2-by-2 block
				// .    structure to exploit. 
				k1 = max( (int)1, ktop-incol );
				nu = ( kdu-max( 0, ndcol-kbot ) ) - k1 + 1;
				// Horizontal Multiply
				for(jcol = min( ndcol, kbot ) + 1; ((nh < 0) ? (jcol >= jbot) : (jcol <= jbot)); jcol += nh){
					jlen = min( nh, jbot-jcol+1 );
					RNP::TBLAS::MultMM<'C','N'>(nu, jlen, nu, 1., &u[(k1-1)+(k1-1)*ldu], ldu, &h[(( incol+k1)-1)+(jcol-1)*ldh], ldh, 0., wh, ldwh );
					RNP::TBLAS::CopyMatrix<'A'>(nu, jlen, wh, ldwh, &h[(( incol+k1)-1)+(jcol-1)*ldh], ldh );
				}
				// Vertical multiply
				for(jrow = jtop; ((nv < 0) ? (jrow >= max( ktop, incol ) - 1) : (jrow <= max( ktop, incol ) - 1)); jrow += nv){
					jlen = min( nv, max( ktop, incol )-jrow );
					RNP::TBLAS::MultMM<'N','N'>(jlen, nu, nu, 1., &h[(jrow-1)+(( incol+k1 )-1)*ldh], ldh, &u[(k1-1)+(k1-1)*ldu], ldu, 0., wv, ldwv );
					RNP::TBLAS::CopyMatrix<'A'>(jlen, nu, wv, ldwv, &h[(jrow-1)+(( incol+k1 )-1)*ldh], ldh );
				}
				// Z multiply (also vertical)
				if( wantz ){
					for(jrow = iloz; ((nv < 0) ? (jrow >= ihiz) : (jrow <= ihiz)); jrow += nv){
						jlen = min( nv, ihiz-jrow+1 );
						RNP::TBLAS::MultMM<'N','N'>(jlen, nu, nu, 1., &z[(jrow-1)+(( incol+k1 )-1)*ldz], ldz, &u[(k1-1)+(k1-1)*ldu], ldu, 0., wv, ldwv );
						RNP::TBLAS::CopyMatrix<'A'>(jlen, nu, wv, ldwv, &z[(jrow-1)+(( incol+k1 )-1)*ldz], ldz );
					}
				}
			}else{
				// Updates exploiting U's 2-by-2 block structure.
				// .    (I2, I4, J2, J4 are the last rows and columns
				// .    of the blocks.)
				i2 = ( kdu+1 ) / 2;
				i4 = kdu;
				j2 = i4 - i2;
				j4 = kdu;
				// KZS and KNZ deal with the band of zeros
				// .    along the diagonal of one of the triangular
				// .    blocks.
				kzs = ( j4-j2 ) - ( ns+1 );
				knz = ns + 1;
				// Horizontal multiply
				for(jcol = min( ndcol, kbot ) + 1; ((nh < 0) ? (jcol >= jbot) : (jcol <= jbot)); jcol += nh){
					jlen = min( nh, jbot-jcol+1 );
					// Copy bottom of H to top+KZS of scratch
					// (The first KZS rows get multiplied by zero.)
					RNP::TBLAS::CopyMatrix<'A'>(knz, jlen, &h[(( incol+1+j2)-1)+(jcol-1)*ldh], ldh, &wh[kzs+0*ldwh], ldwh );
					// Multiply by U21'
					RNP::TBLAS::SetMatrix<'A'>(kzs, jlen, 0., 0., wh, ldwh );
					RNP::TBLAS::MultTrM<'L','U','C','N'>(knz, jlen, one, &u[j2+(( 1+kzs )-1)*ldu], ldu, &wh[kzs+0*ldwh], ldwh );
					// Multiply top of H by U11'
					RNP::TBLAS::MultMM<'C','N'>(i2, jlen, j2, 1., u, ldu, &h[incol+(jcol-1)*ldh], ldh, 1., wh, ldwh );
					// Copy top of H to bottom of WH
					RNP::TBLAS::CopyMatrix<'A'>(j2, jlen, &h[incol+(jcol-1)*ldh], ldh, &wh[i2+0*ldwh], ldwh );
					// Multiply by U21'
					RNP::TBLAS::MultTrM<'L','L','C','N'>(j2, jlen, one, &u[0+(( i2+1 )-1)*ldu], ldu, &wh[i2+0*ldwh], ldwh );
					// Multiply by U22
					RNP::TBLAS::MultMM<'C','N'>(i4-i2, jlen, j4-j2, 1., &u[j2+(( i2+1 )-1)*ldu], ldu, &h[(( incol+1+j2)-1)+(jcol-1)*ldh], ldh, 1., &wh[i2+0*ldwh], ldwh );
					// Copy it back
					RNP::TBLAS::CopyMatrix<'A'>(kdu, jlen, wh, ldwh, &h[incol+(jcol-1)*ldh], ldh );
				}
				// Vertical multiply
				for(jrow = jtop; ((nv < 0) ? (jrow >= max( incol, ktop ) - 1) : (jrow <= max( incol, ktop ) - 1)); jrow += nv){
					jlen = min( nv, max( incol, ktop )-jrow );
					// Copy right of H to scratch (the first KZS columns get multiplied by zero)
					RNP::TBLAS::CopyMatrix<'A'>(jlen, knz, &h[(jrow-1)+(( incol+1+j2 )-1)*ldh], ldh, &wv[0+(( 1+kzs )-1)*ldwv], ldwv );
					// Multiply by U21
					RNP::TBLAS::SetMatrix<'A'>(jlen, kzs, 0., 0., wv, ldwv );
					RNP::TBLAS::MultTrM<'R','U','N','N'>(jlen, knz, one, &u[j2+(( 1+kzs )-1)*ldu], ldu, &wv[0+(( 1+kzs )-1)*ldwv], ldwv );
					// Multiply by U11
					RNP::TBLAS::MultMM<'N','N'>(jlen, i2, j2, 1., &h[(jrow-1)+(( incol+1 )-1)*ldh], ldh, u, ldu, 1., wv, ldwv );
					// 
					// Copy left of H to right of scratch
					// 
					RNP::TBLAS::CopyMatrix<'A'>(jlen, j2, &h[(jrow-1)+(( incol+1 )-1)*ldh], ldh, &wv[0+(( 1+i2 )-1)*ldwv], ldwv );
					// Multiply by U21
					RNP::TBLAS::MultTrM<'R','L','N','N'>(jlen, i4-i2, one, &u[0+(( i2+1 )-1)*ldu], ldu, &wv[0+(( 1+i2 )-1)*ldwv], ldwv );
					// Multiply by U22
					RNP::TBLAS::MultMM<'N','N'>(jlen, i4-i2, j4-j2, 1., &h[(jrow-1)+(( incol+1+j2 )-1)*ldh], ldh, &u[j2+(( i2+1 )-1)*ldu], ldu, 1., &wv[0+(( 1+i2 )-1)*ldwv], ldwv );
					// Copy it back
					RNP::TBLAS::CopyMatrix<'A'>(jlen, kdu, wv, ldwv, &h[(jrow-1)+(( incol+1 )-1)*ldh], ldh );
				}

				// Multiply Z (also vertical)
				if( wantz ){
					for(jrow = iloz; ((nv < 0) ? (jrow >= ihiz) : (jrow <= ihiz)); jrow += nv){
						jlen = min( nv, ihiz-jrow+1 );

						// Copy right of Z to left of scratch (first
						// .     KZS columns get multiplied by zero)
						RNP::TBLAS::CopyMatrix<'A'>(jlen, knz, &z[(jrow-1)+(( incol+1+j2 )-1)*ldz], ldz, &wv[0+(( 1+kzs )-1)*ldwv], ldwv );
						// Multiply by U12
						RNP::TBLAS::SetMatrix<'A'>(jlen, kzs, 0., 0., wv, ldwv );
						RNP::TBLAS::MultTrM<'R','U','N','N'>(jlen, knz, one, &u[j2+(( 1+kzs )-1)*ldu], ldu, &wv[0+(( 1+kzs )-1)*ldwv], ldwv );
						// Multiply by U11
						RNP::TBLAS::MultMM<'N','N'>(jlen, i2, j2, 1., &z[(jrow-1)+(( incol+1 )-1)*ldz], ldz, u, ldu, 1., wv, ldwv );
						// Copy left of Z to right of scratch
						RNP::TBLAS::CopyMatrix<'A'>(jlen, j2, &z[(jrow-1)+(( incol+1 )-1)*ldz], ldz, &wv[0+(( 1+i2 )-1)*ldwv], ldwv );
						// Multiply by U21
						RNP::TBLAS::MultTrM<'R','L','N','N'>(jlen, i4-i2, one, &u[0+(( i2+1 )-1)*ldu], ldu, &wv[0+(( 1+i2 )-1)*ldwv], ldwv );
						// Multiply by U22
						RNP::TBLAS::MultMM<'N','N'>(jlen, i4-i2, j4-j2, 1., &z[(jrow-1)+(( incol+1+j2 )-1)*ldz], ldz, &u[j2+(( i2+1 )-1)*ldu], ldu, 1., &wv[0+(( 1+i2 )-1)*ldwv], ldwv );
						// Copy the result back to Z
						RNP::TBLAS::CopyMatrix<'A'>(jlen, kdu, wv, ldwv, &z[(jrow-1)+(( incol+1 )-1)*ldz], ldz );
					}
				}
			}
		}
	}
}

static size_t zlaqr0_(bool wantt, bool wantz, size_t n, size_t ilo, size_t ihi, std::complex<double> *_h__, size_t ldh, std::complex<double> *_w, size_t iloz, size_t ihiz, std::complex<double> *_z__, size_t ldz, std::complex<double> *_work, bool is_zero)
// is_zero is true if this should be zlaqr0, otherwise this is zlaqr4
{
	std::complex<double> *h__ = (std::complex<double>*)_h__;
	std::complex<double> *w = (std::complex<double>*)_w;
	std::complex<double> *z__ = (std::complex<double>*)_z__;
	std::complex<double> *work = (std::complex<double>*)_work;
	
	/* System generated locals */
	int h_offset, z_offset, i__2, i__3;

	// Purpose
	// =======

	// ZLAQR0 computes the eigenvalues of a Hessenberg matrix H
	// and, optionally, the matrices T and Z from the Schur decomposition
	// H = Z T Z**H, where T is an upper triangular matrix (the
	// Schur form), and Z is the unitary matrix of Schur vectors.

	// Optionally Z may be postmultiplied into an input unitary
	// matrix Q so that this routine can give the Schur factorization
	// of a matrix A which has been reduced to the Hessenberg form H
	// by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H.

	// Arguments
	// =========

	// WANTT   (input) LOGICAL
	//      = .TRUE. : the full Schur form T is required;
	//      = .FALSE.: only eigenvalues are required.

	// WANTZ   (input) LOGICAL
	//      = .TRUE. : the matrix of Schur vectors Z is required;
	//      = .FALSE.: Schur vectors are not required.

	// N     (input) int
	//       The order of the matrix H.  N .GE. 0.

	// ILO   (input) int
	// IHI   (input) int
	//       It is assumed that H is already upper triangular in rows
	//       and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1,
	//       H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
	//       previous call to ZGEBAL, and then passed to ZGEHRD when the
	//       matrix output by ZGEBAL is reduced to Hessenberg form.
	//       Otherwise, ILO and IHI should be set to 1 and N,
	//       respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
	//       If N = 0, then ILO = 1 and IHI = 0.

	// H     (input/output) COMPLEX*16 array, dimension (LDH,N)
	//       On entry, the upper Hessenberg matrix H.
	//       On exit, if INFO = 0 and WANTT is .TRUE., then H
	//       contains the upper triangular matrix T from the Schur
	//       decomposition (the Schur form). If INFO = 0 and WANT is
	//       .FALSE., then the contents of H are unspecified on exit.
	//       (The output value of H when INFO.GT.0 is given under the
	//       description of INFO below.)

	//       This subroutine may explicitly set H(i,j) = 0 for i.GT.j and
	//       j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.

	// LDH   (input) int
	//       The leading dimension of the array H. LDH .GE. max(1,N).

	// W        (output) COMPLEX*16 array, dimension (N)
	//       The computed eigenvalues of H(ILO:IHI,ILO:IHI) are stored
	//       in W(ILO:IHI). If WANTT is .TRUE., then the eigenvalues are
	//       stored in the same order as on the diagonal of the Schur
	//       form returned in H, with W(i) = H(i,i).

	// Z     (input/output) COMPLEX*16 array, dimension (LDZ,IHI)
	//       If WANTZ is .FALSE., then Z is not referenced.
	//       If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
	//       replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
	//       orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
	//       (The output value of Z when INFO.GT.0 is given under
	//       the description of INFO below.)

	// LDZ   (input) int
	//       The leading dimension of the array Z.  if WANTZ is .TRUE.
	//       then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1.

	// WORK  (workspace/output) COMPLEX*16 array, dimension LWORK
	//       On exit, if LWORK = -1, WORK(1) returns an estimate of
	//       the optimal value for LWORK.

	// LWORK (input) int
	//       The dimension of the array WORK.  LWORK .GE. max(1,N)
	//       is sufficient, but LWORK typically as large as 6n may
	//       be required for optimal performance.  A workspace query
	//       to determine the optimal workspace size is recommended.

	//       If LWORK = -1, then ZLAQR0 does a workspace query.
	//       In this case, ZLAQR0 checks the input parameters and
	//       estimates the optimal workspace size for the given
	//       values of N, ILO and IHI.  The estimate is returned
	//       in WORK(1).  No error message related to LWORK is
	//       issued by XERBLA.  Neither H nor Z are accessed.


	// INFO  (output) int
	//         =  0:  successful exit
	//       .GT. 0:  if INFO = i, ZLAQR0 failed to compute all of
	//            the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
	//            and WI contain those eigenvalues which have been
	//            successfully computed.  (Failures are rare.)

	//            If INFO .GT. 0 and WANT is .FALSE., then on exit,
	//            the remaining unconverged eigenvalues are the eigen-
	//            values of the upper Hessenberg matrix rows and
	//            columns ILO through INFO of the final, output
	//            value of H.

	//            If INFO .GT. 0 and WANTT is .TRUE., then on exit

	//       (*)  (initial value of H)*U  = U*(final value of H)

	//            where U is a unitary matrix.  The final
	//            value of  H is upper Hessenberg and triangular in
	//            rows and columns INFO+1 through IHI.

	//            If INFO .GT. 0 and WANTZ is .TRUE., then on exit

	//              (final value of Z(ILO:IHI,ILOZ:IHIZ)
	//               =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U

	//            where U is the unitary matrix in (*) (regard-
	//            less of the value of WANTT.)

	//            If INFO .GT. 0 and WANTZ is .FALSE., then Z is not
	//            accessed.

	// ================================================================
	// Based on contributions by
	//    Karen Braman and Ralph Byers, Department of Mathematics,
	//    University of Kansas, USA

	// ================================================================
	// References:
	//   K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
	//   Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
	//   Performance, SIAM Journal of Matrix Analysis, volume 23, pages
	//   929--947, 2002.

	//   K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
	//   Algorithm Part II: Aggressive Early Deflation, SIAM Journal
	//   of Matrix Analysis, volume 23, pages 948--973, 2002.

	// ================================================================
	// .. Parameters ..

	// ==== Matrices of order NTINY or smaller must be processed by
	// .    ZLAHQR because of insufficient subdiagonal scratch space.
	// .    (This is a hard limit.) ====

	// ==== Exceptional deflation windows:  try to cure rare
	// .    slow convergence by varying the size of the
	// .    deflation window after KEXNW iterations. ====

	// ==== Exceptional shifts: try to cure rare slow convergence
	// .    with ad-hoc exceptional shifts every KEXSH iterations.
	// .    ====

	// ==== The constant WILK1 is used to form the exceptional
	// .    shifts. ====

	using namespace std;

	/* Parameter adjustments */
	h_offset = 1 + ldh;
	h__ -= h_offset;
	--w;
	z_offset = 1 + ldz;
	z__ -= z_offset;
	--work;


	if(n == 0){
		return 0;
	}

	if(n <= 11){ // Tiny matrices must use ZLAHQR.
		return zlahqr_(wantt, wantz, n, ilo, ihi, _h__, ldh, _w, iloz, ihiz, _z__, ldz);
	}
	
	// Use small bulge multi-shift QR with aggressive early
	// deflation on larger-than-tiny matrices.
	size_t info = 0;

	// Set up job flags for ILAENV.

	// NWR = recommended deflation window size.  At this
	// point,  N .GT. NTINY = 11, so there is enough
	// subdiagonal workspace for NWR.GE.2 as required.
	// (In fact, there is enough subdiagonal space for
	// NWR.GE.3.)
	size_t nwr = iparmq_(13, n, ilo, ihi);
	nwr = max((size_t)2, nwr);
	nwr = min(min(ihi-ilo+1, (n-1)/3), nwr);

	// NSR = recommended number of simultaneous shifts.
	// At this point N .GT. NTINY = 11, so there is at
	// enough subdiagonal workspace for NSR to be even
	// and greater than or equal to two as required.
	size_t nsr = iparmq_(15, n, ilo, ihi);
	nsr = min(min(nsr, (n+6)/9), ihi-ilo);
	nsr = max((size_t)2, nsr - nsr % 2);

	// ZLAHQR/ZLAQR0 crossover point
	size_t nmin = iparmq_(12, n, ilo, ihi);
	nmin = max((size_t)11, nmin);

	// Nibble crossover point
	size_t nibble = iparmq_(14, n, ilo, ihi);

	// Accumulate reflections during ttswp?  Use block
	// 2-by-2 structure during matrix-matrix multiply?
	size_t kacc22 = iparmq_(16, n, ilo, ihi);
	kacc22 = max((size_t)0, kacc22);
	kacc22 = min((size_t)2, kacc22);

	// NWMAX = the largest possible deflation window for
	// which there is sufficient workspace.
	size_t nwmax = min((n-1)/3, n/2);
	size_t nw = nwmax;

	// NSMAX = the Largest number of simultaneous shifts
	// for which there is sufficient workspace.
	size_t nsmax = min((n+6)/9, (2*n)/3);
	nsmax -= nsmax % 2;

	// NDFL: an iteration count restarted at deflation.
	size_t ndfl = 1;
	int ndec;

	// ITMAX = iteration limit
	size_t itmax = max((size_t)10, ihi-ilo+1) * 30;

	// Last row and column in the active block
	size_t kbot = ihi;

	// Main Loop
	for(size_t it = 0; it < itmax; ++it){
		// Done when KBOT falls below ILO
		if(kbot < ilo){
			return info;
		}

		// Locate active block

		bool breaked_out = false;
		size_t k = kbot+1; while(k --> ilo+1){
			if(h__[k+(k-1)*ldh] == 0.){
				breaked_out = true;
				break;
			}
		}
		if(!breaked_out){
			k = ilo;
		}
		size_t ktop = k;

		// Select deflation window size:
		// Typical Case:
		//   If possible and advisable, nibble the entire
		//   active block.  If not, use size MIN(NWR,NWMAX)
		//   or MIN(NWR+1,NWMAX) depending upon which has
		//   the smaller corresponding subdiagonal entry
		//   (a heuristic).
		//
		// Exceptional Case:
		//   If there have been no deflations in KEXNW or
		//   more iterations, then vary the deflation window
		//   size.   At first, because, larger windows are,
		//   in general, more powerful than smaller ones,
		//   rapidly increase the window to the maximum possible.
		//   Then, gradually reduce the window size. ====

		size_t nh = kbot - ktop + 1;
		size_t nwupbd = min(nh,nwmax);
		if(ndfl < 5){
			nw = min(nwupbd,nwr);
		} else {
			i__2 = nwupbd, i__3 = 2*nw;
			nw = min(i__2,i__3);
		}
		if(nw < nwmax){
			if(nw >= nh - 1){
				nw = nh;
			} else {
				size_t kwtop = kbot - nw + 1;
				i__2 = kwtop + (kwtop - 1) * ldh;
				i__3 = kwtop - 1 + (kwtop - 2) * ldh;
				if(abs1(h__[kwtop + (kwtop - 1) * ldh]) > abs1(h__[kwtop - 1 + (kwtop - 2) * ldh])){
					++nw;
				}
			}
		}
		if(ndfl < 5){
			ndec = -1;
		} else if(ndec >= 0 || nw >= nwupbd){
			++ndec;
			if(nw - ndec < 2){
				ndec = 0;
			}
			nw -= ndec;
		}

		// Aggressive early deflation:
		// split workspace under the subdiagonal into
		//   - an nw-by-nw work array V in the lower
		//     left-hand-corner,
		//   - an NW-by-at-least-NW-but-more-is-better
		//     (NW-by-NHO) horizontal work array along
		//     the bottom edge,
		//   - an at-least-NW-but-more-is-better (NHV-by-NW)
		//     vertical work array along the left-hand-edge.
		size_t kv = n - nw + 1;
		size_t kt = nw + 1;
		size_t nho = n - nw - 1 - kt + 1;
		size_t kwv = nw + 2;
		size_t nve = n - nw - kwv + 1;

		// Aggressive early deflation
		size_t ls, ld;
		zlaqr3_(wantt, wantz, n, ktop, kbot, nw, &h__[h_offset], ldh, iloz, ihiz, &z__[z_offset], ldz, &ls, &ld, &w[1], &h__[kv + ldh], ldh, nho, &h__[kv + kt * ldh], ldh, nve, &h__[kwv + ldh], ldh, &work[1], is_zero);

		// Adjust KBOT accounting for new deflations.
		kbot -= ld;

		// KS points to the shifts.
		size_t ks = kbot - ls + 1;

		// Skip an expensive QR sweep if there is a (partly
		// heuristic) reason to expect that many eigenvalues
		// will deflate without it.  Here, the QR sweep is
		// skipped if many eigenvalues have just been deflated
		// or if the remaining active block is small.
		if(ld == 0 || (ld * 100 <= nw * nibble && kbot - ktop + 1 > min(nmin,nwmax))){
			// NS = nominal number of simultaneous shifts.
			// This may be lowered (slightly) if ZLAQR3
			// did not provide that many shifts.
			size_t ns = min(min(nsmax,nsr), max((size_t)2, kbot - ktop));
			ns -= ns % 2;

			// If there have been no deflations
			// in a multiple of KEXSH iterations,
			// then try exceptional shifts.
			// Otherwise use shifts provided by
			// ZLAQR3 above or from the eigenvalues
			// of a trailing principal submatrix.
			if(ndfl % 6 == 0){
				ks = kbot - ns + 1;
				for(int i = kbot; i >= (int)ks+1; i -= 2){
					w[i] = h__[i+i*ldh] + abs1(h__[i+(i-1)*ldh]) * .75;
					w[i-1] = w[i];
				}
			} else {
				// Got NS/2 or fewer shifts? Use ZLAQR4 or
				// ZLAHQR on a trailing principal submatrix to
				// get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
				// there is enough space below the subdiagonal
				// to fit an NS-by-NS scratch array.)
				if(kbot - ks + 1 <= ns / 2){
					ks = kbot - ns + 1;
					kt = n - ns + 1;
					RNP::TBLAS::CopyMatrix<'A'>(ns, ns, &h__[ks + ks * ldh], ldh, &h__[kt + ldh], ldh);
					size_t inf;
					if(ns > nmin && is_zero){
						inf = zlaqr0_(false, false, ns, 1, ns, &h__[kt + ldh], ldh, &w[ks], 1, 1, NULL, 1, _work, false);
					} else {
						inf = zlahqr_(false, false, ns, 1, ns, &h__[kt + ldh], ldh, &w[ks], 1, 1, NULL, 1);
					}
					ks += inf;

					// In case of a rare QR failure use
					// eigenvalues of the trailing 2-by-2
					// principal submatrix.  Scale to avoid
					// overflows, underflows and subnormals.
					// (The scale factor S can not be zero,
					// because H(KBOT,KBOT-1) is nonzero.)
					if(ks >= kbot){
						double s =
							abs1(h__[kbot-1+(kbot-1)*ldh]) +
							abs1(h__[kbot  +(kbot-1)*ldh]) +
							abs1(h__[kbot-1+ kbot   *ldh]) +
							abs1(h__[kbot  + kbot   *ldh]);
						std::complex<double> aa = h__[kbot-1+(kbot-1)*ldh] / s;
						std::complex<double> cc = h__[kbot  +(kbot-1)*ldh] / s;
						std::complex<double> bb = h__[kbot-1+ kbot   *ldh] / s;
						std::complex<double> dd = h__[kbot  + kbot   *ldh] / s;
						std::complex<double> tr2 = (aa+dd)/2.;
						std::complex<double> det = (aa-tr2)*(dd-tr2) - bb*cc;
						std::complex<double> rtdisc = sqrt(-det);
						w[kbot-1] = (tr2+rtdisc)*s;
						w[kbot] = (tr2-rtdisc)*s;
						ks = kbot - 1;
					}
				}

				if(kbot - ks + 1 > ns){
					// Sort the shifts (Helps a little)
					bool sorted = false;
					k = kbot+1; while(k --> ks+1){
						if(sorted){
							break;
						}
						sorted = true;
						for(size_t i = ks; i < k; ++i){
							if(abs1(w[i]) < abs1(w[i+1])){
								sorted = false;
								std::swap(w[i], w[i+1]);
							}
						}
					}
				}
			}

			// If there are only two shifts, then use only one.
			if(kbot - ks + 1 == 2){
				if(abs1(w[kbot]-h__[kbot+kbot*ldh]) < abs1(w[kbot-1]-h__[kbot+kbot*ldh])){
					w[kbot-1] = w[kbot];
				} else {
					w[kbot] = w[kbot-1];
				}
			}

			// Use up to NS of the the smallest magnatiude
			// shifts.  If there aren't NS shifts available,
			// then use them all, possibly dropping one to
			// make the number of shifts even.
			ns = min(ns, kbot - ks + 1);
			ns -= ns % 2;
			ks = kbot - ns + 1;

			// Small-bulge multi-shift QR sweep:
			// split workspace under the subdiagonal into
			// - a KDU-by-KDU work array U in the lower
			//   left-hand-corner,
			// - a KDU-by-at-least-KDU-but-more-is-better
			//   (KDU-by-NHo) horizontal work array WH along
			//   the bottom edge,
			// - and an at-least-KDU-but-more-is-better-by-KDU
			//   (NVE-by-KDU) vertical work WV arrow along
			//   the left-hand-edge.
			size_t kdu = ns * 3 - 3;
			size_t ku = n - kdu + 1;
			size_t kwh = kdu + 1;
			nho = n - kdu - 3 - (kdu + 1) + 1;
			kwv = kdu + 4;
			nve = n - kdu - kwv + 1;

			// Small-bulge multi-shift QR sweep
			zlaqr5(wantt, wantz, kacc22, n, ktop, kbot, ns, &w[ks], _h__, ldh, iloz, ihiz, _z__, ldz, _work, 3, &h__[ku + ldh], ldh, nve, &h__[kwv + ldh], ldh, nho, &h__[ku + kwh * ldh], ldh);
		}

		// Note progress (or the lack of it).
		if(ld > 0){
			ndfl = 1;
		} else {
			++ndfl;
		}

		// End of main loop
	}

	// Iteration limit exceeded.  Set INFO to show where
	// the problem occurred and exit.
	return kbot;
}

static int zhseqr_(char job, char compz, size_t n, size_t ilo, size_t ihi, std::complex<double> *h, size_t ldh, std::complex<double> *w, std::complex<double> *z, size_t ldz, std::complex<double> *_work)
{
	//   Purpose
	//   =======

	//   ZHSEQR computes the eigenvalues of a Hessenberg matrix H
	//   and, optionally, the matrices T and Z from the Schur decomposition
	//   H = Z T Z**H, where T is an upper triangular matrix (the
	//   Schur form), and Z is the unitary matrix of Schur vectors.

	//   Optionally Z may be postmultiplied into an input unitary
	//   matrix Q so that this routine can give the Schur factorization
	//   of a matrix A which has been reduced to the Hessenberg form H
	//   by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H.

	//   Arguments
	//   =========

	//   JOB   (input) CHARACTER*1
	//         = 'E':  compute eigenvalues only;
	//         = 'S':  compute eigenvalues and the Schur form T.

	//   COMPZ (input) CHARACTER*1
	//         = 'N':  no Schur vectors are computed;
	//         = 'I':  Z is initialized to the unit matrix and the matrix Z
	//                 of Schur vectors of H is returned;
	//         = 'V':  Z must contain an unitary matrix Q on entry, and
	//                 the product Q*Z is returned.

	//   N     The order of the matrix H.  N .GE. 0.

	//   ILO   (input) int
	//   IHI   (input) int
	//         It is assumed that H is already upper triangular in rows
	//         and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
	//         set by a previous call to ZGEBAL, and then passed to ZGEHRD
	//         when the matrix output by ZGEBAL is reduced to Hessenberg
	//         form. Otherwise ILO and IHI should be set to 1 and N
	//         respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
	//         If N = 0, then ILO = 1 and IHI = 0.

	//   H     (input/output) COMPLEX*16 array, dimension (LDH,N)
	//         On entry, the upper Hessenberg matrix H.
	//         On exit, if INFO = 0 and JOB = 'S', H contains the upper
	//         triangular matrix T from the Schur decomposition (the
	//         Schur form). If INFO = 0 and JOB = 'E', the contents of
	//         H are unspecified on exit.  (The output value of H when
	//         INFO.GT.0 is given under the description of INFO below.)

	//         Unlike earlier versions of ZHSEQR, this subroutine may
	//         explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1
	//         or j = IHI+1, IHI+2, ... N.

	//   LDH   The leading dimension of the array H. LDH .GE. max(1,N).

	//   W        (output) COMPLEX*16 array, dimension (N)
	//         The computed eigenvalues. If JOB = 'S', the eigenvalues are
	//         stored in the same order as on the diagonal of the Schur
	//         form returned in H, with W(i) = H(i,i).

	//   Z     (input/output) COMPLEX*16 array, dimension (LDZ,N)
	//         If COMPZ = 'N', Z is not referenced.
	//         If COMPZ = 'I', on entry Z need not be set and on exit,
	//         if INFO = 0, Z contains the unitary matrix Z of the Schur
	//         vectors of H.  If COMPZ = 'V', on entry Z must contain an
	//         N-by-N matrix Q, which is assumed to be equal to the unit
	//         matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,
	//         if INFO = 0, Z contains Q*Z.
	//         Normally Q is the unitary matrix generated by ZUNGHR
	//         after the call to ZGEHRD which formed the Hessenberg matrix
	//         H. (The output value of Z when INFO.GT.0 is given under
	//         the description of INFO below.)

	//   LDZ   The leading dimension of the array Z.  if COMPZ = 'I' or
	//         COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1.

	//   WORK  (workspace/output) COMPLEX*16 array, dimension N
	//         On exit, if INFO = 0, WORK(1) returns an estimate of
	//         the optimal value for LWORK.

	//   INFO  (output) int
	//           =  0:  successful exit
	//         .LT. 0:  if INFO = -i, the i-th argument had an illegal
	//                  value
	//         .GT. 0:  if INFO = i, ZHSEQR failed to compute all of
	//              the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
	//              and WI contain those eigenvalues which have been
	//              successfully computed.  (Failures are rare.)

	//              If INFO .GT. 0 and JOB = 'E', then on exit, the
	//              remaining unconverged eigenvalues are the eigen-
	//              values of the upper Hessenberg matrix rows and
	//              columns ILO through INFO of the final, output
	//              value of H.

	//              If INFO .GT. 0 and JOB   = 'S', then on exit

	//         (*)  (initial value of H)*U  = U*(final value of H)

	//              where U is a unitary matrix.  The final
	//              value of  H is upper Hessenberg and triangular in
	//              rows and columns INFO+1 through IHI.

	//              If INFO .GT. 0 and COMPZ = 'V', then on exit

	//                (final value of Z)  =  (initial value of Z)*U

	//              where U is the unitary matrix in (*) (regard-
	//              less of the value of JOB.)

	//              If INFO .GT. 0 and COMPZ = 'I', then on exit
	//                    (final value of Z)  = U
	//              where U is the unitary matrix in (*) (regard-
	//              less of the value of JOB.)

	//              If INFO .GT. 0 and COMPZ = 'N', then Z is not
	//              accessed.

	//   ================================================================
	//           Default values supplied by
	//           ILAENV(ISPEC,'ZHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).
	//           It is suggested that these defaults be adjusted in order
	//           to attain best performance in each particular
	//           computational environment.

	//          ISPEC=12: The ZLAHQR vs ZLAQR0 crossover point.
	//                    Default: 75. (Must be at least 11.)

	//          ISPEC=13: Recommended deflation window size.
	//                    This depends on ILO, IHI and NS.  NS is the
	//                    number of simultaneous shifts returned
	//                    by ILAENV(ISPEC=15).  (See ISPEC=15 below.)
	//                    The default for (IHI-ILO+1).LE.500 is NS.
	//                    The default for (IHI-ILO+1).GT.500 is 3*NS/2.

	//          ISPEC=14: Nibble crossover point. (See IPARMQ for
	//                    details.)  Default: 14% of deflation window
	//                    size.

	//          ISPEC=15: Number of simultaneous shifts in a multishift
	//                    QR iteration.

	//                    If IHI-ILO+1 is ...

	//                    greater than      ...but less    ... the
	//                    or equal to ...      than        default is

	//                         1               30          NS =   2(+)
	//                        30               60          NS =   4(+)
	//                        60              150          NS =  10(+)
	//                       150              590          NS =  **
	//                       590             3000          NS =  64
	//                      3000             6000          NS = 128
	//                      6000             infinity      NS = 256

	//                (+)  By default some or all matrices of this order
	//                     are passed to the implicit double shift routine
	//                     ZLAHQR and this parameter is ignored.  See
	//                     ISPEC=12 above and comments in IPARMQ for
	//                     details.

	//               (**)  The asterisks (**) indicate an ad-hoc
	//                     function of N increasing from 10 to 64.

	//          ISPEC=16: Select structured matrix multiply.
	//                    If the number of simultaneous shifts (specified
	//                    by ISPEC=15) is less than 14, then the default
	//                    for ISPEC=16 is 0.  Otherwise the default for
	//                    ISPEC=16 is 2.

	//   ================================================================
	//   Based on contributions by
	//      Karen Braman and Ralph Byers, Department of Mathematics,
	//      University of Kansas, USA

	//   ================================================================
	//   References:
	//     K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
	//     Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
	//     Performance, SIAM Journal of Matrix Analysis, volume 23, pages
	//     929--947, 2002.

	//     K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
	//     Algorithm Part II: Aggressive Early Deflation, SIAM Journal
	//     of Matrix Analysis, volume 23, pages 948--973, 2002.

	//   ================================================================
	//   .. Parameters ..

	//   ==== Matrices of order NTINY or smaller must be processed by
	//   .    ZLAHQR because of insufficient subdiagonal scratch space.
	//   .    (This is a hard limit.) ====

	//   ==== NL allocates some local workspace to help small matrices
	//   .    through a rare ZLAHQR failure.  NL .GT. NTINY = 11 is
	//   .    required and NL .LE. NMIN = ILAENV(ISPEC=12,...) is recom-
	//   .    mended.  (The default value of NMIN is 75.)  Using NL = 49
	//   .    allows up to six simultaneous shifts and a 16-by-16
	//   .    deflation window.  ====
	
	using namespace std;
	
	std::complex<double> hl[49*49];
	int kbot;
	std::complex<double> workl[49];

	/* Function Body */
	const bool wantt = (job == 'S');
	const bool initz = (compz == 'I');
	const bool wantz = initz || (compz == 'V');

	int info = 0;

	if (n == 0) {
		return 0;
	}

	// copy eigenvalues isolated by ZGEBAL
	if (ilo > 1) {
		RNP::TBLAS::Copy(ilo - 1, h, ldh+1, w, 1);
	}
	if (ihi < n) {
		RNP::TBLAS::Copy(n - ihi, &h[ihi+ihi*ldh], ldh + 1, &w[ihi], 1);
	}

	// Initialize Z, if requested
	if (initz) {
		RNP::TBLAS::SetMatrix<'A'>(n, n, 0., 1., z, ldz);
	}

	// Quick return if possible
	if(ilo == ihi){
		w[ilo-1] = h[(ilo-1)+(ilo-1)*ldh];
		return 0;
	}

	// ZLAHQR/ZLAQR0 crossover point
	size_t nmin = iparmq_(12, n, ilo, ihi);
	nmin = max((size_t)11,nmin);

	// ZLAQR0 for big matrices; ZLAHQR for small ones
	if (n > nmin) {
		info = zlaqr0_(wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi, z, ldz, _work, true);
		//info = zlaqr0_(wantt, wantz, n, ilo, ihi, (doublecomplex*)h, ldh, (doublecomplex*)w, ilo, ihi, (doublecomplex*)z, ldz, (doublecomplex*)_work, true);
	} else { // Small matrix
		info = zlahqr_(wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi, z, ldz);
		if (info > 0) {
			// A rare ZLAHQR failure! ZLAQR0 sometimes succeeds when ZLAHQR fails.
			kbot = info;

			if (n >= 49) {
				// Larger matrices have enough subdiagonal scratch
				// space to call ZLAQR0 directly.
				info = zlaqr0_(wantt, wantz, n, ilo, kbot, h, ldh, w, ilo, ihi, z, ldz, _work, true);
				//info = zlaqr0_(wantt, wantz, n, ilo, kbot, (doublecomplex*)h, ldh, (doublecomplex*)w, ilo, ihi, (doublecomplex*)z, ldz, (doublecomplex*)_work, true);
			} else {
				// Tiny matrices don't have enough subdiagonal
				// scratch space to benefit from ZLAQR0.  Hence,
				// tiny matrices must be copied into a larger
				// array before calling ZLAQR0.
				RNP::TBLAS::CopyMatrix<'A'>(n, n, h, ldh, hl, 49);
				hl[n+(n-1)*49] = 0;
				RNP::TBLAS::SetMatrix<'A'>(49, 49-n, 0., 0., &hl[0+n*49], 49);
				info = zlaqr0_(wantt, wantz, 49, ilo, kbot, hl, 49, w, ilo, ihi, z, ldz, workl, true);
				//info = zlaqr0_(wantt, wantz, 49, ilo, kbot, (doublecomplex*)hl, 49, (doublecomplex*)w, ilo, ihi, (doublecomplex*)z, ldz, (doublecomplex*)workl, true);
				if (wantt || info != 0) {
					RNP::TBLAS::CopyMatrix<'A'>(n, n, hl, 49, h, ldh);
				}
			}
		}
	}

	// Clear out the trash, if necessary.
	if ((wantt || info != 0) && n > 2) {
		RNP::TBLAS::SetMatrix<'L'>(n-2, n-2, 0., 0., &h[2+0*ldh], ldh);
	}

	return info;
}

static void zunghr_(size_t n, size_t ilo, size_t ihi, std::complex<double> *a, size_t lda, std::complex<double> *tau, std::complex<double> *work){
	using namespace std;

	// Purpose
	// =======
	// 
	// ZUNGHR generates a complex unitary matrix Q which is defined as the
	// product of IHI-ILO elementary reflectors of order N, as returned by
	// ZGEHRD:
	// 
	// Q = H(ilo) H(ilo+1) . . . H(ihi-1).
	// 
	// Arguments
	// =========
	// 
	// N       (input) int
	// The order of the matrix Q. N >= 0.
	// 
	// ILO     (input) int
	// IHI     (input) int
	// ILO and IHI must have the same values as in the previous call
	// of ZGEHRD. Q is equal to the unit matrix except in the
	// submatrix Q(ilo+1:ihi,ilo+1:ihi).
	// 1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
	// 
	// A       (input/output) std::complex<double> array, dimension (LDA,N)
	// On entry, the vectors which define the elementary reflectors,
	// as returned by ZGEHRD.
	// On exit, the N-by-N unitary matrix Q.
	// 
	// LDA     (input) int
	// The leading dimension of the array A. LDA >= max(1,N).
	// 
	// TAU     (input) std::complex<double> array, dimension (N-1)
	// TAU(i) must contain the scalar factor of the elementary
	// reflector H(i), as returned by ZGEHRD.
	// 
	// WORK    (workspace/output) std::complex<double> array, dimension (MAX(1,LWORK))
	// On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
	// 
	// LWORK   (input) int
	// The dimension of the array WORK. LWORK >= IHI-ILO.
	// For optimum performance LWORK >= (IHI-ILO)*NB, where NB is
	// the optimal blocksize.

	size_t nh = ihi - ilo;

	if(n == 0){ return; }
	// 
	// Shift the vectors which define the elementary reflectors one
	// column to the right, and set the first ilo and the last n-ihi
	// rows and columns to those of the unit matrix
	// 
	size_t j = ihi+1;
	while(j --> ilo+1){
		for(size_t i = 0; i < j; ++i){
			a[i+j*lda] = 0;
		}
		for(size_t i = j+1; i <= ihi; ++i){
			a[i+j*lda] = a[i+(j-1)*lda];
		}
		for(size_t i = ihi+1; i < n; ++i){
			a[i+j*lda] = 0;
		}
	}
	for(size_t j = 0; j <= ilo; ++j){
		for(size_t i = 0; i < n; ++i){
			a[i+j*lda] = 0;
		}
		a[j+j*lda] = 1;
	}
	for(size_t j = ihi+1; j < n; ++j){
		for(size_t i = 0; i < n; ++i){
			a[i+j*lda] = 0;
		}
		a[j+j*lda] = 1;
	}
	if(nh > 0){ // Generate Q(ilo+1:ihi,ilo+1:ihi)	
		RNP::TLASupport::GenerateOrthognalMatrixFromElementaryReflectors(nh, nh, nh, &a[ilo+1+(ilo+1)*lda], lda, &tau[ilo], work);
	}
}

int RNP::Eigensystem(size_t n, 
	std::complex<double> *a, size_t lda,
	std::complex<double> *w,
	std::complex<double> *vl, size_t ldvl, std::complex<double> *vr, size_t ldvr,
	std::complex<double> *work_, double *rwork_, size_t lwork_)
{
	using namespace std;

	// Purpose
	// =======

	// ZGEEV computes for an N-by-N complex nonsymmetric matrix A, the
	// eigenvalues and, optionally, the left and/or right eigenvectors.

	// The right eigenvector v(j) of A satisfies
	//                  A * v(j) = lambda(j) * v(j)
	// where lambda(j) is its eigenvalue.
	// The left eigenvector u(j) of A satisfies
	//               u(j)**H * A = lambda(j) * u(j)**H
	// where u(j)**H denotes the conjugate transpose of u(j).

	// The computed eigenvectors are normalized to have Euclidean norm
	// equal to 1 and largest component real.

	// Arguments
	// =========

	// N       (input) INTEGER
	//         The order of the matrix A. N >= 0.

	// A       (input/output) COMPLEX*16 array, dimension (LDA,N)
	//         On entry, the N-by-N matrix A.
	//         On exit, A has been overwritten.

	// LDA     The leading dimension of the array A.  LDA >= max(1,N).

	// W       (output) COMPLEX*16 array, dimension (N)
	//         W contains the computed eigenvalues.

	// VL      (output) COMPLEX*16 array, dimension (LDVL,N)
	//         If not NULL, the left eigenvectors u(j) are stored one
	//         after another in the columns of VL, in the same order
	//         as their eigenvalues.
	//         u(j) = VL(:,j), the j-th column of VL.

	// LDVL    The leading dimension of the array VL.  LDVL >= 1; if
	//         JOBVL = 'V', LDVL >= N.

	// VR      (output) COMPLEX*16 array, dimension (LDVR,N)
	//         If not NULL, the right eigenvectors v(j) are stored one
	//         after another in the columns of VR, in the same order
	//         as their eigenvalues.
	//         v(j) = VR(:,j), the j-th column of VR.

	// LDVR    The leading dimension of the array VR.  LDVR >= 1; if
	//         JOBVR = 'V', LDVR >= N.

	// WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,2*N))

	// RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)

	// return  = 0:  successful exit
	//         < 0:  if INFO = -i, the i-th argument had an illegal value.
	//         > 0:  if INFO = i, the QR algorithm failed to compute all the
	//               eigenvalues, and no eigenvectors have been computed;
	//               elements and i+1:N of W contain eigenvalues which have
	//               converged.

	const bool wantvl = (vl != NULL);
	const bool wantvr = (vr != NULL);

	if(n == 0) {
		return 0;
	}
	if((size_t)-1 == lwork_){
		work_[0] = (double)(2*n);
		return 0;
	}
	
	std::complex<double> *work = work_;
	double *rwork = rwork_;
	if(NULL == work_ || lwork_ < 2*n){
		work = new std::complex<double>[2*n];
	}
	if(NULL == rwork_){
		rwork = new double[2*n];
	}

	const double eps = 2*std::numeric_limits<double>::epsilon();
	const double smlnum = sqrt(std::numeric_limits<double>::min()) / eps;
	const double bignum = 1. / smlnum;

	// Scale A if max element outside range [SMLNUM,BIGNUM]
	double anrm;
	RNP::TLASupport::CheapMatrixNorm<'M'>(n, n, a, lda, &anrm);
	bool scalea = false;
	double cscale = 1;
	if(anrm > 0. && anrm < smlnum){
		scalea = true;
		cscale = smlnum;
	}else if(anrm > bignum){
		scalea = true;
		cscale = bignum;
	}
	if(scalea){
		RNP::TLASupport::RescaleMatrix<'G'>(0, 0, anrm, cscale, n, n, a, lda);
	}

	// Balance the matrix
	// (CWorkspace: none)
	// (RWorkspace: need N)

	const size_t ibal = 0;
	size_t ilo, ihi;
	zgebal_('B', n, a, lda, &ilo, &ihi, &rwork[ibal]);

	// Reduce to upper Hessenberg form
	// (CWorkspace: need 2*N, prefer N+N*NB)
	// (RWorkspace: none)

	const size_t itau = 0;
	size_t iwrk = itau + n;
	RNP::TLASupport::HessenbergReduction(n, ilo-1, ihi-1, a, lda, &work[itau], &work[iwrk], n);

	int info;
	if(wantvl){
		RNP::TBLAS::CopyMatrix<'L'>(n, n, a, lda, vl, ldvl);

		// Generate unitary matrix in VL
		// (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
		// (RWorkspace: none)

		zunghr_(n, ilo-1, ihi-1, vl, ldvl, &work[itau], &work[iwrk]);

		// Perform QR iteration, accumulating Schur vectors in VL
		// (CWorkspace: need 1, prefer HSWORK (see comments) )
		// (RWorkspace: none)

		iwrk = itau;
		info = zhseqr_('S', 'V', n, ilo, ihi, a, lda, w, vl, ldvl, &work[iwrk]);

		if (wantvr) {
			// Want left and right eigenvectors
			// Copy Schur vectors to VR
			RNP::TBLAS::CopyMatrix<'F'>(n, n, vl, ldvl, vr, ldvr);
		}
	}else if(wantvr){
		// Want right eigenvectors
		// Copy Householder vectors to VR
		RNP::TBLAS::CopyMatrix<'L'>(n, n, a, lda, vr, ldvr);
		// Generate unitary matrix in VR
		// (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
		// (RWorkspace: none)

		zunghr_(n, ilo-1, ihi-1, vr, ldvr, &work[itau], &work[iwrk]);

		// Perform QR iteration, accumulating Schur vectors in VR
		// (CWorkspace: need 1, prefer HSWORK (see comments) )
		// (RWorkspace: none)

		iwrk = itau;
		info = zhseqr_('S', 'V', n, ilo, ihi, a, lda, w, vr, ldvr, &work[iwrk]);
	}else{
		// Compute eigenvalues only
		// (CWorkspace: need 1, prefer HSWORK (see comments) )
		// (RWorkspace: none)
		iwrk = itau;
		zhseqr_('E', 'N', n, ilo, ihi, a, lda, w, vr, ldvr, &work[iwrk]);
	}

	if(info == 0){
		size_t irwork = ibal + n;
		if(wantvl || wantvr){
			// Compute left and/or right eigenvectors
			// (CWorkspace: need 2*N)
			// (RWorkspace: need 2*N)
			size_t nout;
			ztrevc_('B', NULL, n, a, lda, vl, ldvl, vr, ldvr, n, &nout, &work[iwrk], &rwork[irwork]);
		}

		if(wantvl){
			// Undo balancing of left eigenvectors
			// (CWorkspace: none)
			// (RWorkspace: need N)
			zgebak_('B', 'L', n, ilo-1, ihi-1, &rwork[ibal], n, vl, ldvl);

			// Normalize left eigenvectors and make largest component real
			for(size_t i = 0; i < n; ++i){
				double scl = 1. / RNP::TBLAS::Norm2(n, &vl[0+i*ldvl], 1);
				RNP::TBLAS::Scale(n, scl, &vl[0+i*ldvl], 1);
				for(size_t k = 0; k < n; ++k){
					rwork[irwork+k] = std::norm(vl[k+i*ldvl]);
				}
				size_t k = RNP::TBLAS::MaximumIndex(n, &rwork[irwork], 1);
				RNP::TBLAS::Scale(n, std::conj(vl[k+i*ldvl]) / sqrt(rwork[irwork+k]), &vl[0+i*ldvl], 1);
				vl[k+i*ldvl] = vl[k+i*ldvl].real();
			}
		}

		if(wantvr){
			// Undo balancing of right eigenvectors
			// (CWorkspace: none)
			// (RWorkspace: need N)
			zgebak_('B', 'R', n, ilo-1, ihi-1, &rwork[ibal], n, vr, ldvr);

			// Normalize right eigenvectors and make largest component real
			for(size_t i = 0; i < n; ++i){
				double scl = 1. / RNP::TBLAS::Norm2(n, &vr[0+i*ldvr], 1);
				RNP::TBLAS::Scale(n, scl, &vr[0+i*ldvr], 1);
				for(size_t k = 0; k < n; ++k){
					rwork[irwork + k] = std::norm(vr[k+i*ldvr]);
				}
				size_t k = RNP::TBLAS::MaximumIndex(n, &rwork[irwork], 1);
				RNP::TBLAS::Scale(n, std::conj(vr[k+i*ldvr]) / sqrt(rwork[irwork+k]), &vr[0+i*ldvr], 1);
				vr[k+i*ldvr] = vr[k+i*ldvr].real();
			}
		}
	}
	
	// Undo scaling if necessary
	if(scalea){
		RNP::TLASupport::RescaleMatrix<'G'>(0, 0, cscale, anrm, n - info, 1, &w[info], max(n-info,size_t(1)));
		if(info > 0){
			RNP::TLASupport::RescaleMatrix<'G'>(0, 0, cscale, anrm, ilo-1, 1, w, n);
		}
	}
	
	if(NULL == work_ || lwork_ < 2*n){
		delete [] work;
	}
	if(NULL == rwork_){
		delete [] rwork;
	}
	
	return info;
}


