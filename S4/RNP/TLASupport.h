#ifndef _RNP_TLASUPPORT_H_
#define _RNP_TLASUPPORT_H_

/* LAPACK license:

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

Many routines here are not copied from the official LAPACK distribution
and are in the public domain.

*/

#include <iostream>

#include <complex>
#include <cmath>
#include <limits>
#include "TBLAS.h"

/*
 * Preprocessor flags:
 *   RNP_NUMERIC_ROBUST
 */

// It's almost ALWAYS a bad idea to not define it
#define RNP_NUMERIC_ROBUST

namespace RNP{
namespace TLASupport{

template <class T>
inline T Pythag2(T x, T y){ // dlapy2
	using namespace std;
	// returns sqrt(x**2+y**2), taking care not to cause unnecessary overflow.
    T a = abs(x);
    T b = abs(y);
    if(a < b){ swap(a,b); }
    if(0 == b){
		return a;
    }else{
		b /= a;
		return a * sqrt(b*b + T(1));
    }
}

template <class T>
inline T Pythag3(T x, T y, T z){ // dlapy3
	using namespace std;
	// returns sqrt(x**2+y**2+z**2), taking care not to cause unnecessary overflow.
    T xabs = abs(x);
    T yabs = abs(y);
    T zabs = abs(z);
	T w = xabs;
	if(yabs > w){ w = yabs; }
	if(zabs > w){ w = zabs; }
    if(0 == w){
		// W can be zero for max(0,nan,0)
		// adding all three entries together will make sure
		// NaN will not disappear.
		return xabs + yabs + zabs;
    }else{
		xabs /= w;
		yabs /= w;
		zabs /= w;
		return w * sqrt(xabs*xabs + yabs*yabs + zabs*zabs);
    }
}

template <class T>
void ComplexDivide(const std::complex<T> &num, const std::complex<T> &den, std::complex<T> *result){ // zladiv, dladiv, cladiv, sladiv
	// The algorithm is due to Robert L. Smith and can be found
	// in D. Knuth, The art of Computer Programming, Vol.2, p.195
	using namespace std;
	if(abs(den.imag()) < abs(den.real())){
		T e = den.imag()/den.real();
		T f = den.real() + den.imag()*e;
		*result = std::complex<double>(
			(num.real() + num.imag()*e) / f,
			(num.imag() - num.real()*e) / f
			);
	}else{
		T e = den.real()/den.imag();
		T f = den.imag() + den.real()*e;
		*result = std::complex<double>(
			(num.imag() + num.real()*e) / f,
			(num.imag()*e - num.real()) / f
			);
	}
}
template <class T>
std::complex<T> ComplexDivide(const std::complex<T> &num, const std::complex<T> &den){ // zladiv, dladiv, cladiv, sladiv
	// The algorithm is due to Robert L. Smith and can be found
	// in D. Knuth, The art of Computer Programming, Vol.2, p.195
	using namespace std;
	if(abs(den.imag()) < abs(den.real())){
		T e = den.imag()/den.real();
		T f = den.real() + den.imag()*e;
		return std::complex<double>(
			(num.real() + num.imag()*e) / f,
			(num.imag() - num.real()*e) / f
			);
	}else{
		T e = den.real()/den.imag();
		T f = den.imag() + den.real()*e;
		return std::complex<double>(
			(num.imag() + num.real()*e) / f,
			(num.imag()*e - num.real()) / f
			);
	}
}

template <class T>
size_t LastNonzeroColumn(size_t m, size_t n, const T *a, size_t lda){ // ilazlc, iladlc, ilaclc, ilaslc
	// if n = 0, returns -1
	if(n < 1 || T(0) != a[0+(n-1)*lda] || T(0) != a[m-1+(n-1)*lda]){ return n-1; }
	size_t j = n;
	while(j --> 0){
		for(size_t i = 0; i < m; ++i){
			if(T(0) != a[i+j*lda]){ return j; }
		}
	}
	return j; // should be -1, but unsigned
}

template <class T>
size_t LastNonzeroRow(size_t m, size_t n, const T *a, size_t lda){ // ilazlr, iladlr, ilaclr, ilaslr
	// if m = 0, returns -1
	if(m < 1 || T(0) != a[m-1+0*lda] || T(0) != a[m-1+(n-1)*lda]){ return m-1; }
	size_t i = m;
	while(i --> 0){
		for(size_t j = 0; j < n; ++j){
			if(T(0) != a[i+j*lda]){ return i; }
		}
	}
	return i; // should be -1, but unsigned
}

template <char norm='M'>
struct CheapHessenbergNorm{ // zlanhs, dlanhs, clanhs, slanhs
	template <class T>
	CheapHessenbergNorm(size_t n, const T *a, size_t lda, typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type *result, typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type *work = NULL){
		typedef typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type real_type;
		*result = 0;
		if(n < 1){ return; }
		if('M' == norm){ // max(abs(A(i,j)))
			for(size_t j = 0; j < n; ++j){
				size_t ilimit = j+2; if(n < ilimit){ ilimit = n; }
				for(size_t i = 0; i < ilimit; ++i){
					real_type ca = RNP::TBLAS::_RealOrComplexChooser<T>::_abs(a[i+j*lda]);
					if(ca > *result){ *result = ca; }
				}
			}
		}else if('O' == norm || '1' == norm){ // max col sum
			for(size_t j = 0; j < n; ++j){
				size_t ilimit = j+2; if(n < ilimit){ ilimit = n; }
				real_type sum = 0;
				for(size_t i = 0; i < ilimit; ++i){
					sum += RNP::TBLAS::_RealOrComplexChooser<T>::_abs(a[i+j*lda]);
				}
				if(sum > *result){ *result = sum; }
			}
		}else if('I' == norm){ // max row sum
			if(NULL == work){ // can't accumulate row sums
				for(size_t i = 0; i < n; ++i){
					size_t jstart = 0; if(i > 0){ jstart = i-1; }
					real_type sum = 0;
					for(size_t j = 0; j < n; ++j){
						sum += RNP::TBLAS::_RealOrComplexChooser<T>::_abs(a[i+j*lda]);
					}
					if(sum > *result){ *result = sum; }
				}
			}else{ // accumulate row sums in a cache-friendlier traversal order
				for(size_t i = 0; i < n; ++i){ work[i] = 0; }
				for(size_t j = 0; j < n; ++j){
					size_t ilimit = j+2; if(n < ilimit){ ilimit = n; }
					for(size_t i = 0; i < ilimit; ++i){
						work[i] += RNP::TBLAS::_RealOrComplexChooser<T>::_abs(a[i+j*lda]);
					}
				}
				for(size_t i = 0; i < n; ++i){
					if(work[i] > *result){ *result = work[i]; }
				}
			}
		}else if('F' == norm || 'E' == norm){ // Frobenius norm
			using namespace std;
#ifdef RNP_NUMERIC_ROBUST
			real_type scale = 0;
			real_type sum = 1;
			for(size_t j = 0; j < n; ++j){
				size_t ilimit = j+2; if(n < ilimit){ ilimit = n; }
				real_type sum = 0;
				for(size_t i = 0; i < ilimit; ++i){
					real_type ca = RNP::TBLAS::_RealOrComplexChooser<T>::_abs(a[i+j*lda]);
					if(scale < ca){
						real_type r = scale/ca;
						sum = real_type(1) + sum*r*r;
						scale = ca;
					}else{
						real_type r = ca/scale;
						sum += r*r;
					}
				}
			}
			*result = scale*sqrt(sum);
#else
			real_type sum = 0;
			for(size_t j = 0; j < n; ++j){
				size_t ilimit = j+2; if(n < ilimit){ ilimit = n; }
				real_type sum = 0;
				for(size_t i = 0; i < ilimit; ++i){
					sum += RNP::TBLAS::_RealOrComplexChooser<T>::_abs2(a[i+j*lda]);
				}
			}
			*result = sqrt(sum);
#endif
		}
	}
};

template <char norm='M'>
struct CheapMatrixNorm{ // zlange, dlange, clange, slange
	template <class T>
	CheapMatrixNorm(size_t m, size_t n, const T *a, size_t lda, typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type *result, typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type *work = NULL){
		typedef typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type real_type;
		*result = 0;
		if(n < 1){ return; }
		if('M' == norm){ // max(abs(A(i,j)))
			for(size_t j = 0; j < n; ++j){
				for(size_t i = 0; i < m; ++i){
					real_type ca = RNP::TBLAS::_RealOrComplexChooser<T>::_abs(a[i+j*lda]);
					if(ca > *result){ *result = ca; }
				}
			}
		}else if('O' == norm || '1' == norm){ // max col sum
			for(size_t j = 0; j < n; ++j){
				real_type sum = 0;
				for(size_t i = 0; i < m; ++i){
					sum += RNP::TBLAS::_RealOrComplexChooser<T>::_abs(a[i+j*lda]);
				}
				if(sum > *result){ *result = sum; }
			}
		}else if('I' == norm){ // max row sum
			if(NULL == work){ // can't accumulate row sums
				for(size_t i = 0; i < m; ++i){
					size_t jstart = 0; if(i > 0){ jstart = i-1; }
					real_type sum = 0;
					for(size_t j = 0; j < n; ++j){
						sum += RNP::TBLAS::_RealOrComplexChooser<T>::_abs(a[i+j*lda]);
					}
					if(sum > *result){ *result = sum; }
				}
			}else{ // accumulate row sums in a cache-friendlier traversal order
				for(size_t i = 0; i < m; ++i){ work[i] = 0; }
				for(size_t j = 0; j < n; ++j){
					for(size_t i = 0; i < m; ++i){
						work[i] += RNP::TBLAS::_RealOrComplexChooser<T>::_abs(a[i+j*lda]);
					}
				}
				for(size_t i = 0; i < m; ++i){
					if(work[i] > *result){ *result = work[i]; }
				}
			}
		}else if('F' == norm || 'E' == norm){ // Frobenius norm
			using namespace std;
#ifdef RNP_NUMERIC_ROBUST
			real_type scale = 0;
			real_type sum = 1;
			for(size_t j = 0; j < n; ++j){
				real_type sum = 0;
				for(size_t i = 0; i < m; ++i){
					real_type ca = RNP::TBLAS::_RealOrComplexChooser<T>::_abs(a[i+j*lda]);
					if(scale < ca){
						real_type r = scale/ca;
						sum = real_type(1) + sum*r*r;
						scale = ca;
					}else{
						real_type r = ca/scale;
						sum += r*r;
					}
				}
			}
			*result = scale*sqrt(sum);
#else
			real_type sum = 0;
			for(size_t j = 0; j < n; ++j){
				real_type sum = 0;
				for(size_t i = 0; i < m; ++i){
					sum += RNP::TBLAS::_RealOrComplexChooser<T>::_abs2(a[i+j*lda]);
				}
			}
			*result = sqrt(sum);
#endif
		}
	}
};

template <char type='G'> // zlascl, dlascl, clascl, slascl
struct RescaleMatrix{
	template <class TS, class T>
	RescaleMatrix(size_t kl, size_t ku, const TS &cfrom, const TS &cto, size_t m, size_t n, T *a, size_t lda){
		//  Purpose
		//  =======

		//  ZLASCL multiplies the M by N complex matrix A by the real scalar
		//  CTO/CFROM.  This is done without over/underflow as long as the final
		//  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
		//  A may be full, upper triangular, lower triangular, upper Hessenberg,
		//  or banded.

		//  Arguments
		//  =========

		//  TYPE    (input) CHARACTER*1
		//      TYPE indices the storage type of the input matrix.
		//      = 'G':  A is a full matrix.
		//      = 'L':  A is a lower triangular matrix.
		//      = 'U':  A is an upper triangular matrix.
		//      = 'H':  A is an upper Hessenberg matrix.
		//      = 'B':  A is a symmetric band matrix with lower bandwidth KL
		//              and upper bandwidth KU and with the only the lower
		//              half stored.
		//      = 'Q':  A is a symmetric band matrix with lower bandwidth KL
		//              and upper bandwidth KU and with the only the upper
		//              half stored.
		//      = 'Z':  A is a band matrix with lower bandwidth KL and upper
		//              bandwidth KU.

		//  KL      (input) int
		//      The lower bandwidth of A.  Referenced only if TYPE = 'B',
		//      'Q' or 'Z'.

		//  KU      (input) int
		//      The upper bandwidth of A.  Referenced only if TYPE = 'B',
		//      'Q' or 'Z'.

		//  CFROM   (input) DOUBLE PRECISION
		//  CTO     (input) DOUBLE PRECISION
		//      The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
		//      without over/underflow if the final result CTO*A(I,J)/CFROM
		//      can be represented without over/underflow.  CFROM must be
		//      nonzero.

		//  M       (input) int
		//      The number of rows of the matrix A.  M >= 0.

		//  N       (input) int
		//      The number of columns of the matrix A.  N >= 0.

		//  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
		//      The matrix to be multiplied by CTO/CFROM.  See TYPE for the
		//      storage type.

		//  LDA     (input) int
		//      The leading dimension of the array A.  LDA >= max(1,M).

		if(n == 0 || m == 0){ return; }

		const TS smlnum = std::numeric_limits<TS>::min();
		const TS bignum = 1. / smlnum;

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
				}else if(std::abs(cfrom1) > std::abs(ctoc) && ctoc != 0.){
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

			if(type == 'G'){ // Full matrix
				for(size_t j = 0; j < n; ++j){
					for(size_t i = 0; i < m; ++i){
						a[i+j*lda] *= mul;
					}
				}
			}else if(type == 'L'){ // Lower triangular matrix
				for(size_t j = 0; j < n; ++j){
					for(size_t i = j; i < m; ++i){
						a[i+j*lda] *= mul;
					}
				}
			}else if(type == 'U'){
				// Upper triangular matrix
				for(size_t j = 0; j < n; ++j){
					size_t ilimit = j+1; if(m < ilimit){ ilimit = m; }
					for(size_t i = 0; i < ilimit; ++i){
						a[i+j*lda] *= mul;
					}
				}
			}else if(type == 'H'){
				// Upper Hessenberg matrix
				for(size_t j = 0; j < n; ++j){
					size_t ilimit = j+2; if(m < ilimit){ ilimit = m; };
					for(size_t i = 0; i < ilimit; ++i){
						a[i+j*lda] *= mul;
					}
				}
			}else if(type == 'B'){
				// Lower half of a symmetric band matrix
				for(size_t j = 0; j < n; ++j){
					size_t ilimit = n-j; if(kl+1 < ilimit){ ilimit = kl+1; }
					for(size_t i = 0; i < ilimit; ++i){
						a[i+j*lda] *= mul;
					}
				}
			}else if(type == 'Q'){
				// Upper half of a symmetric band matrix
				for(size_t j = 0; j < n; ++j){
					size_t istart = (ku > j) ? ku-j : 0;
					for(size_t i = istart; i <= ku; ++i){
						a[i+j*lda] *= mul;
					}
				}
			}else if(type == 'Z'){
				// Band matrix
				size_t k3 = 2*kl + ku + 1;
				for(size_t j = 0; j < n; ++j){
					size_t istart = kl+ku-j;
					if(kl > istart){ istart = kl; }
					size_t ilimit = kl + ku + m-j;
					if(k3 < ilimit){ ilimit = k3; }
					for(size_t i = istart; i < ilimit; ++i){
						a[i+j*lda] *= mul;
					}
				}
			}
		}while(!done);
	}
};

template <char side='L'>
struct ApplyElementaryReflector{ // zlarf, dlarf, clarf, slarf
	template <class T>
	ApplyElementaryReflector(size_t m, size_t n, const T *v, size_t incv, const T &tau, T *c, size_t ldc, T *work){
		// Applies an elementary reflector H to an M-by-N matrix C,
		// from either the left or the right. H is represented in the form
		//        H = I - tau * v * v'
		// where tau is a complex scalar and v is a complex vector.
		//
		// To apply H' (the conjugate transpose of H), supply conj(tau) instead

		// Arguments
		// =========
		//
		// SIDE    = 'L': form  H * C
		//         = 'R': form  C * H
		//
		// M       The number of rows of the matrix C.
		//
		// N       The number of columns of the matrix C.
		//
		// V                  (1 + (M-1)*abs(INCV)) if SIDE = 'L'
		//                 or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
		//         The vector v in the representation of H. V is not used if
		//         TAU = 0.
		//
		// INCV    The increment between elements of v. INCV <> 0.
		//
		// TAU     The value tau in the representation of H.
		//
		// C       (input/output) COMPLEX*16 array, dimension (LDC,N)
		//         On entry, the M-by-N matrix C.
		//         On exit, C is overwritten by the matrix H * C if SIDE = 'L',
		//         or C * H if SIDE = 'R'.
		//
		// LDC     The leading dimension of the array C. LDC >= max(1,M).
		//
		// WORK    (workspace) COMPLEX*16 array, dimension
		//                        (N) if SIDE = 'L'
		//                     or (M) if SIDE = 'R'

		size_t lenv = 0;
		size_t lenc = 0;
		if(T(0) != tau){
			// Set up variables for scanning V.  LASTV begins pointing to the end of V.
			if('L' == side){
				lenv = m;
			}else{
				lenv = n;
			}
			size_t i = 0;
			if(incv > 0){
				i = (lenv - 1) * incv;
			}
			// Look for the last non-zero row in V.
			while(lenv > 0 && (T(0) == v[i])){
				--lenv;
				i -= incv;
			}
			if('L' == side){ // Scan for the last non-zero column in C(1:lastv,:).
				lenc = 1+RNP::TLASupport::LastNonzeroColumn(lenv, n, c, ldc);
			}else{ // Scan for the last non-zero row in C(:,1:lastv).
				lenc = 1+RNP::TLASupport::LastNonzeroRow(m, lenv, c, ldc);
			}
		}
		
		if('L' == side){ // Form  H * C
			if(lenv > 0 && lenc > 0){
				// w(1:lastc,1) := C(1:lastv,1:lastc)' * v(1:lastv,1)
				RNP::TBLAS::MultMV<'C'>(lenv, lenc, T(1), c, ldc, v, incv, T(0), work, 1);
				// C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)'
				RNP::TBLAS::ConjugateRank1Update(lenv, lenc, -tau, v, incv, work, 1, c, ldc);
			}
		}else{ // Form  C * H
			if(lenv > 0 && lenc > 0){
				// w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
				RNP::TBLAS::MultMV<'N'>(lenc, lenv, T(1), c, ldc, v, incv, T(0), work, 1);
				// C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)'
				RNP::TBLAS::ConjugateRank1Update(lenc, lenv, -tau, work, 1, v, incv, c, ldc);
			}
		}
	}
};

template <char side, char trans, char direct, char storev>
struct ApplyElementaryReflectorBlocked{ // zlarfb, dlarfb, clarfb, slarfb
	template <class T>
	ApplyElementaryReflectorBlocked(size_t m, size_t n, size_t k, const T *v, size_t ldv, const T *t, size_t ldt, T *c, size_t ldc, T *work, size_t ldwork){
		//  ZLARFB applies a complex block reflector H or its transpose H' to a
		//  complex M-by-N matrix C, from either the left or the right.
		//
		//  Arguments
		//  =========
		//
		//  SIDE    (input) CHARACTER*1
		//          = 'L': apply H or H' from the Left
		//          = 'R': apply H or H' from the Right
		//
		//  TRANS   (input) CHARACTER*1
		//          = 'N': apply H (No transpose)
		//          = 'C': apply H' (Conjugate transpose)
		//
		//  DIRECT  (input) CHARACTER*1
		//          Indicates how H is formed from a product of elementary
		//          reflectors
		//          = 'F': H = H(1) H(2) . . . H(k) (Forward)
		//          = 'B': H = H(k) . . . H(2) H(1) (Backward)
		//
		//  STOREV  (input) CHARACTER*1
		//          Indicates how the vectors which define the elementary
		//          reflectors are stored:
		//          = 'C': Columnwise
		//          = 'R': Rowwise
		//
		//  M       (input) INTEGER
		//          The number of rows of the matrix C.
		//
		//  N       (input) INTEGER
		//          The number of columns of the matrix C.
		//
		//  K       (input) INTEGER
		//          The order of the matrix T (= the number of elementary
		//          reflectors whose product defines the block reflector).
		//
		//  V       (input) COMPLEX*16 array, dimension
		//                                (LDV,K) if STOREV = 'C'
		//                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
		//                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
		//          The matrix V. See further details.
		//
		//  LDV     (input) INTEGER
		//          The leading dimension of the array V.
		//          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
		//          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
		//          if STOREV = 'R', LDV >= K.
		//
		//  T       (input) COMPLEX*16 array, dimension (LDT,K)
		//          The triangular K-by-K matrix T in the representation of the
		//          block reflector.
		//
		//  LDT     (input) INTEGER
		//          The leading dimension of the array T. LDT >= K.
		//
		//  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
		//          On entry, the M-by-N matrix C.
		//          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
		//
		//  LDC     (input) INTEGER
		//          The leading dimension of the array C. LDC >= max(1,M).
		//
		//  WORK    (workspace) COMPLEX*16 array, dimension (LDWORK,K)
		//
		//  LDWORK  (input) INTEGER
		//          The leading dimension of the array WORK.
		//          If SIDE = 'L', LDWORK >= max(1,N);
		//          if SIDE = 'R', LDWORK >= max(1,M).
		if(m < 1 || n < 1){ return; }
		
		// System generated locals
		int c_offset, t_offset, v_offset, work_offset;

		// Local variables
		int i, j;
		int lastc;
		int lastv;

		// Parameter adjustments
		v_offset = 1 + ldv;
		v -= v_offset;
		t_offset = 1 + ldt;
		t -= t_offset;
		c_offset = 1 + ldc;
		c -= c_offset;
		work_offset = 1 + ldwork;
		work -= work_offset;


		if(storev== 'C'){
			if(direct== 'F'){
				// Let  V =  ( V1 )    (first K rows)
				//           ( V2 )
				// where  V1  is unit lower triangular.
				if(side== 'L'){
					// Form  H * C  or  H' * C  where  C = ( C1 )
					//                                     ( C2 )

					lastv = 1+LastNonzeroColumn(m, k, &v[v_offset], ldv); if(k > lastv){ lastv = k; }
					lastc = 1+LastNonzeroColumn(lastv, n, &c[c_offset], ldc);

					// W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
					// W := C1'
					for(j = 1; j <= k; ++j){
						RNP::TBLAS::Copy(lastc, &c[j + ldc], ldc, &work[j * ldwork + 1], 1);
						RNP::TBLAS::Conjugate(lastc, &work[j * ldwork + 1], 1);
					}

					// W := W * V1
					RNP::TBLAS::MultTrM<'R','L','N','U'>(lastc, k, 1., &v[v_offset], ldv, &work[work_offset], ldwork);
					if(lastv > k){
						// W := W + C2'*V2
						RNP::TBLAS::MultMM<'C','N'>(lastc, k, lastv - k, 1., &c[k + 1 + ldc], ldc, &v[k + 1 + ldv], ldv, 1., &work[work_offset], ldwork);
					}

					// W := W * T'  or  W * T
					if('N' == trans){
						RNP::TBLAS::MultTrM<'R','U','C','N'>(lastc, k, 1., &t[t_offset], ldt, &work[work_offset], ldwork);
					}else{
						RNP::TBLAS::MultTrM<'R','U','N','N'>(lastc, k, 1., &t[t_offset], ldt, &work[work_offset], ldwork);
					}

					// C := C - V * W'
					if(m > k){
						// C2 := C2 - V2 * W'
						RNP::TBLAS::MultMM<'N','C'>(lastv - k, lastc, k, -1., &v[k + 1 + ldv], ldv, &work[work_offset], ldwork, 1., &c[k + 1 + ldc], ldc);
					}

					// W := W * V1'
					RNP::TBLAS::MultTrM<'R','L','C','U'>(lastc, k, 1., &v[v_offset], ldv, &work[work_offset], ldwork);

					// C1 := C1 - W'
					for(j = 1; j <= k; ++j){
						for(i = 1; i <= lastc; ++i){
							c[j + i * ldc] = c[j + i * ldc] - std::conj(work[i + j * ldwork]);
						}
					}
				}else if(side== 'R'){
					// Form  C * H  or  C * H'  where  C = ( C1  C2 )

					lastv = 1+LastNonzeroColumn(n, k, &v[v_offset], ldv); if(k > lastv){ lastv = k; }
					lastc = 1+LastNonzeroColumn(m, lastv, &c[c_offset], ldc);

					// W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
					// W := C1
					for(j = 1; j <= k; ++j){
						RNP::TBLAS::Copy(lastc, &c[j * ldc + 1], 1, &work[j * ldwork + 1], 1);
					}

					// W := W * V1
					RNP::TBLAS::MultTrM<'R','L','N','U'>(lastc, k, 1., &v[v_offset], ldv, &work[work_offset], ldwork);
					if(lastv > k){
						// W := W + C2 * V2
						RNP::TBLAS::MultMM<'N','N'>(lastc, k, lastv - k, 1., &c[(k + 1) * ldc + 1], ldc, &v[k + 1 + ldv], ldv, 1., &work[work_offset], ldwork);
					}

					// W := W * T  or  W * T'
					RNP::TBLAS::MultTrM<'R','U',trans,'N'>(lastc, k, 1., &t[t_offset], ldt, &work[work_offset], ldwork);

					// C := C - W * V'
					if(lastv > k){
						// C2 := C2 - W * V2'
						RNP::TBLAS::MultMM<'N','C'>(lastc, lastv - k, k, -1., &work[work_offset], ldwork, &v[k + 1 + ldv], ldv, 1., &c[(k + 1) * ldc + 1], ldc);
					}

					// W := W * V1'
					RNP::TBLAS::MultTrM<'R','L','C','U'>(lastc, k, 1., &v[v_offset], ldv, &work[work_offset], ldwork);

					// C1 := C1 - W
					for(j = 1; j <= k; ++j){
						for(i = 1; i <= lastc; ++i){
							c[i + j * ldc] = c[i + j * ldc] - work[i + j * ldwork];
						}
					}
				}
			}else{
				// Let  V =  ( V1 )
				//           ( V2 )    (last K rows)
				// where  V2  is unit upper triangular.
				if(side=='L'){
					// Form  H * C  or  H' * C  where  C = ( C1 )
					//                                     ( C2 )

					lastv = 1+LastNonzeroColumn(m, k, &v[v_offset], ldv); if(k > lastv){ lastv = k; }
					lastc = 1+LastNonzeroColumn(lastv, n, &c[c_offset], ldc);

					// W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
					// W := C2'
					for(j = 1; j <= k; ++j){
						RNP::TBLAS::Copy(lastc, &c[lastv - k + j + ldc], ldc, &work[j * ldwork + 1], 1);
						RNP::TBLAS::Conjugate(lastc, &work[j * ldwork + 1], 1);
					}

					// W := W * V2
					RNP::TBLAS::MultTrM<'R','U','N','U'>(lastc, k, 1., &v[lastv - k + 1 + ldv], ldv, &work[
						work_offset], ldwork);
					if(lastv > k){
						// W := W + C1'*V1
						RNP::TBLAS::MultMM<'C','N'>(lastc, k, lastv - k, 1., &c[c_offset], ldc, &v[v_offset], ldv, 1., &work[work_offset], ldwork);
					}

					// W := W * T'  or  W * T
					if('N' == trans){
						RNP::TBLAS::MultTrM<'R','L','C','N'>(lastc, k, 1., &t[t_offset], ldt, &work[work_offset], ldwork);
					}else{
						RNP::TBLAS::MultTrM<'R','L','N','N'>(lastc, k, 1., &t[t_offset], ldt, &work[work_offset], ldwork);
					}

					// C := C - V * W'
					if(lastv > k){
						// C1 := C1 - V1 * W'
						RNP::TBLAS::MultMM<'N','C'>(lastv - k, lastc, k, -1., &v[v_offset], ldv, &work[work_offset], ldwork, 1., &c[c_offset], ldc);
					}

					// W := W * V2'
					RNP::TBLAS::MultTrM<'R','U','C','U'>(lastc, k, 1., &v[lastv - k + 1 + ldv], ldv, &work[work_offset], ldwork);

					// C2 := C2 - W'
					for(j = 1; j <= k; ++j){
						for(i = 1; i <= lastc; ++i){
							c[lastv - k + j + i * ldc] = c[lastv - k + j + i * ldc] - std::conj(work[i + j * ldwork]);
						}
					}

				}else if(side=='R'){
					// Form  C * H  or  C * H'  where  C = ( C1  C2 )
					lastv = 1+LastNonzeroColumn(n, k, &v[v_offset], ldv); if(k > lastv){ lastv = k; }
					lastc = 1+LastNonzeroColumn(m, lastv, &c[c_offset], ldc);

					// W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
					// W := C2
					for(j = 1; j <= k; ++j){
						RNP::TBLAS::Copy(lastc, &c[(lastv - k + j) * ldc + 1], 1, &work[j * ldwork + 1], 1);
					}

					// W := W * V2
					RNP::TBLAS::MultTrM<'R','U','N','U'>(lastc, k, 1., &v[lastv - k + 1 + ldv], ldv, &work[
						work_offset], ldwork);
					if(lastv > k){

					// W := W + C1 * V1
						RNP::TBLAS::MultMM<'N','N'>(lastc, k, lastv - k, 1., &c[c_offset], ldc, &v[v_offset], ldv, 1., &work[work_offset], ldwork);
					}

					// W := W * T  or  W * T'
					RNP::TBLAS::MultTrM<'R','L',trans,'N'>(lastc, k, 1., &t[t_offset], ldt, &work[work_offset], ldwork);

					// C := C - W * V'
					if(lastv > k){
						// C1 := C1 - W * V1'
						RNP::TBLAS::MultMM<'N','C'>(lastc, lastv - k, k, -1., &work[work_offset], ldwork, &v[v_offset], ldv, 1., &c[c_offset], ldc);
					}

					// W := W * V2'
					RNP::TBLAS::MultTrM<'R','U','C','U'>(lastc, k, 1., &v[lastv - k + 1 + ldv], ldv, &work[work_offset], ldwork);

					// C2 := C2 - W
					for(j = 1; j <= k; ++j){
						for(i = 1; i <= lastc; ++i){
							c[i + (lastv - k + j) * ldc] = c[i + (lastv - k + j) * ldc] - work[i + j * ldwork];
						}
					}
				}
			}

		}else if(storev=='R'){
			if(direct=='F'){
				// Let  V =  ( V1  V2 )    (V1: first K columns)
				// where  V1  is unit upper triangular.

				if(side=='L'){
					// Form  H * C  or  H' * C  where  C = ( C1 )
					//                                     ( C2 )

					lastv = 1+LastNonzeroColumn(k, m, &v[v_offset], ldv); if(k > lastv){ lastv = k; }
					lastc = 1+LastNonzeroColumn(lastv, n, &c[c_offset], ldc);

					// W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
					// W := C1'
					for(j = 1; j <= k; ++j){
						RNP::TBLAS::Copy(lastc, &c[j + ldc], ldc, &work[j * ldwork + 1], 1);
						RNP::TBLAS::Conjugate(lastc, &work[j * ldwork + 1], 1);
					}

					// W := W * V1'
					RNP::TBLAS::MultTrM<'R','U','C','U'>(lastc, k, 1., &v[v_offset], ldv, &work[work_offset], ldwork);
					if(lastv > k){

						// W := W + C2'*V2'
						RNP::TBLAS::MultMM<'C','C'>(lastc, k, lastv - k, 1., &c[k + 1 + ldc], ldc, &v[(k + 1) * ldv + 1], ldv, 1., &work[work_offset], ldwork);
					}

					// W := W * T'  or  W * T
					if('N' == trans){
						RNP::TBLAS::MultTrM<'R','U','C','N'>(lastc, k, 1., &t[t_offset], ldt, &work[work_offset], ldwork);
					}else{
						RNP::TBLAS::MultTrM<'R','U','N','N'>(lastc, k, 1., &t[t_offset], ldt, &work[work_offset], ldwork);
					}

					// C := C - V' * W'
					if(lastv > k){
						// C2 := C2 - V2' * W'
						RNP::TBLAS::MultMM<'C','C'>(lastv - k, lastc, k, -1., &v[(k + 1) * ldv + 1], ldv, &work[work_offset], ldwork, 1., &c[k + 1 + ldc], ldc);
					}

					// W := W * V1
					RNP::TBLAS::MultTrM<'R','U','N','U'>(lastc, k, 1., &v[v_offset], ldv, &work[work_offset], ldwork);

					// C1 := C1 - W'
					for(j = 1; j <= k; ++j){
						for(i = 1; i <= lastc; ++i){
							c[j + i * ldc] = c[j + i * ldc] - std::conj(work[i + j * ldwork]);
						}
					}

				}else if(side=='R'){
					// Form  C * H  or  C * H'  where  C = ( C1  C2 )
					lastv = 1+LastNonzeroColumn(k, n, &v[v_offset], ldv); if(k > lastv){ lastv = k; }
					lastc = 1+LastNonzeroColumn(m, lastv, &c[c_offset], ldc);

					// W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
					// W := C1
					for(j = 1; j <= k; ++j){
						RNP::TBLAS::Copy(lastc, &c[j * ldc + 1], 1, &work[j * ldwork + 1], 1);
					}

					// W := W * V1'
					RNP::TBLAS::MultTrM<'R','U','C','U'>(lastc, k, 1., &v[v_offset], ldv, &work[work_offset], ldwork);
					if(lastv > k){
						// W := W + C2 * V2'
						RNP::TBLAS::MultMM<'N','C'>(lastc, k, lastv - k, 1., &c[(k + 1) * ldc + 1], ldc, &v[(k + 1) * ldv + 1], ldv, 1., &work[work_offset], ldwork);
					}

					// W := W * T  or  W * T'
					RNP::TBLAS::MultTrM<'R','U',trans,'N'>(lastc, k, 1., &t[t_offset], ldt, &work[work_offset], ldwork);

					// C := C - W * V
					if(lastv > k){
						// C2 := C2 - W * V2
						RNP::TBLAS::MultMM<'N','N'>(lastc, lastv - k, k, -1., &work[work_offset], ldwork, &v[(k + 1) * ldv + 1], ldv, 1., &c[(k + 1) * ldc + 1], ldc);
					}

					// W := W * V1
					RNP::TBLAS::MultTrM<'R','U','N','U'>(lastc, k, 1., &v[v_offset], ldv, &work[work_offset], ldwork);

					// C1 := C1 - W
					for(j = 1; j <= k; ++j){
						for(i = 1; i <= lastc; ++i){
							c[i + j * ldc] = c[i + j * ldc] - work[i + j * ldwork];
						}
					}
				}
			}else{
				// Let  V =  ( V1  V2 )    (V2: last K columns)
				// where  V2  is unit lower triangular.
				if(side=='L'){
					// Form  H * C  or  H' * C  where  C = ( C1 )
					// ( C2 )

					lastv = 1+LastNonzeroColumn(k, m, &v[v_offset], ldv); if(k > lastv){ lastv = k; }
					lastc = 1+LastNonzeroColumn(lastv, n, &c[c_offset], ldc);

					// W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
					// W := C2'
					for(j = 1; j <= k; ++j){
						RNP::TBLAS::Copy(lastc, &c[lastv - k + j + ldc], ldc, &work[j * ldwork + 1], 1);
						RNP::TBLAS::Conjugate(lastc, &work[j * ldwork + 1], 1);
					}

					// W := W * V2'
					RNP::TBLAS::MultTrM<'R','L','C','U'>(lastc, k, 1., &v[(lastv - k + 1) * ldv + 1], ldv, &work[work_offset], ldwork);
					if(lastv > k){
						// W := W + C1'*V1'
						RNP::TBLAS::MultMM<'C','C'>(lastc, k, lastv - k, 1., &c[c_offset], ldc, &v[v_offset], ldv, 1., &work[work_offset], ldwork);
					}

					// W := W * T'  or  W * T
					if('N' == trans){
						RNP::TBLAS::MultTrM<'R','L','C','N'>(lastc, k, 1., &t[t_offset], ldt, &work[work_offset], ldwork);
					}else{
						RNP::TBLAS::MultTrM<'R','L','N','N'>(lastc, k, 1., &t[t_offset], ldt, &work[work_offset], ldwork);
					}

					// C := C - V' * W'
					if(lastv > k){
						// C1 := C1 - V1' * W'
						RNP::TBLAS::MultMM<'C','C'>(lastv - k, lastc, k, -1., &v[v_offset], ldv, &work[work_offset], ldwork, 1., &c[c_offset], ldc);
					}

					// W := W * V2
					RNP::TBLAS::MultTrM<'R','L','N','U'>(lastc, k, 1., &v[(lastv - k + 1) * ldv + 1], ldv, &work[
						work_offset], ldwork);

					// C2 := C2 - W'
					for(j = 1; j <= k; ++j){
						for(i = 1; i <= lastc; ++i){
							c[lastv - k + j + i * ldc] = c[lastv - k + j + i * ldc] - std::conj(work[i + j * ldwork]);
						}
					}

				}else if(side=='R'){
					// Form  C * H  or  C * H'  where  C = ( C1  C2 )
					lastv = 1+LastNonzeroColumn(k, n, &v[v_offset], ldv); if(k > lastv){ lastv = k; }
					lastc = 1+LastNonzeroColumn(m, lastv, &c[c_offset], ldc);

					// W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
					// W := C2
					for(j = 1; j <= k; ++j){
						RNP::TBLAS::Copy(lastc, &c[(lastv - k + j) * ldc + 1], 1, &work[j * ldwork + 1], 1);
					}

					// W := W * V2'
					RNP::TBLAS::MultTrM<'R','L','C','U'>(lastc, k, 1., &v[(lastv - k + 1) * ldv + 1], ldv, &work[work_offset], ldwork);
					if(lastv > k){
						// W := W + C1 * V1'
						RNP::TBLAS::MultMM<'N','C'>(lastc, k, lastv - k, 1., &c[c_offset], ldc, &v[v_offset], ldv, 1., &work[work_offset], ldwork);
					}

					// W := W * T  or  W * T'
					RNP::TBLAS::MultTrM<'R','L',trans,'N'>(lastc, k, 1., &t[t_offset], ldt, &work[work_offset], ldwork);

					// C := C - W * V
					if(lastv > k){
						// C1 := C1 - W * V1
						RNP::TBLAS::MultMM<'N','N'>(lastc, lastv - k, k, -1., &work[work_offset], ldwork, &v[v_offset], ldv, 1., &c[c_offset], ldc);
					}

					// W := W * V2
					RNP::TBLAS::MultTrM<'R','L','N','U'>(lastc, k, 1., &v[(lastv - k + 1) * ldv + 1], ldv, &work[work_offset], ldwork);

					// C1 := C1 - W
					for(j = 1; j <= k; ++j){
						for(i = 1; i <= lastc; ++i){
							c[i + (lastv - k + j) * ldc] = c[i + (lastv - k + j) * ldc] - work[i + j * ldwork];
						}
					}
				}
			}
		}
	}
};


template <class T> // zlarfg, dlarfg, clarfg, slarfg
void GenerateElementaryReflector(size_t n, T *alpha, T *x, size_t incx, T *tau){
	typedef typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type real_type;

	// Purpose
	// =======
	// 
	// ZLARFG generates a complex elementary reflector H of order n, such
	// that
	// 
	// H' * ( alpha ) = ( beta ),   H' * H = I.
	// (   x   )   (   0  )
	// 
	// where alpha and beta are scalars, with beta real, and x is an
	// (n-1)-element complex vector. H is represented in the form
	// 
	// H = I - tau * ( 1 ) * ( 1 v' ) ,
	// ( v )
	// 
	// where tau is a complex scalar and v is a complex (n-1)-element
	// vector. Note that H is not hermitian.
	// 
	// If the elements of x are all zero and alpha is real, then tau = 0
	// and H is taken to be the unit matrix.
	// 
	// Otherwise  1 <= real(tau) <= 2  and  std::abs(tau-1) <= 1 .
	// 
	// Arguments
	// =========
	// 
	// N       The order of the elementary reflector.
	// 
	// ALPHA   (input/output) std::complex<double>
	// On entry, the value alpha.
	// On exit, it is overwritten with the value beta.
	// 
	// X       (input/output) std::complex<double> array, dimension
	// (1+(N-2)*std::abs(INCX))
	// On entry, the vector x.
	// On exit, it is overwritten with the vector v.
	// 
	// INCX    The increment between elements of X. INCX > 0.
	// 
	// TAU     (output) std::complex<double>
	// The value tau.
	// 
	using namespace std;

	if(n == 0){
		*tau = 0;
		return;
	}
	
	real_type xnorm = RNP::TBLAS::Norm2( n-1, x, incx );
	if(xnorm == 0 && RNP::TBLAS::_RealOrComplexChooser<T>::_imag(*alpha) == 0){ // H  =  I
		*tau = 0;
	}else{ // general case
		real_type beta = RNP::TLASupport::Pythag3( RNP::TBLAS::_RealOrComplexChooser<T>::_real(*alpha), RNP::TBLAS::_RealOrComplexChooser<T>::_imag(*alpha), xnorm );
		if(RNP::TBLAS::_RealOrComplexChooser<T>::_real(*alpha) > 0){ beta = -beta; };
		const real_type safmin = std::numeric_limits<real_type>::min() / (2*std::numeric_limits<real_type>::epsilon());
		const real_type rsafmn = real_type(1) / safmin;
		// 
		size_t knt = 0;
		if(abs(beta) < safmin){ // XNORM, BETA may be inaccurate; scale X and recompute them
			do{
				++knt;
				RNP::TBLAS::Scale( n-1, rsafmn, x, incx );
				beta *= rsafmn;
				*alpha *= rsafmn;
			}while(abs(beta) < safmin );
			// New BETA is at most 1, at least SAFMIN
			xnorm = RNP::TBLAS::Norm2( n-1, x, incx );
			beta = RNP::TLASupport::Pythag3( RNP::TBLAS::_RealOrComplexChooser<T>::_real(*alpha), RNP::TBLAS::_RealOrComplexChooser<T>::_imag(*alpha), xnorm );
			if(RNP::TBLAS::_RealOrComplexChooser<T>::_real(*alpha) > 0){ beta = -beta; }
		}
		*tau = (beta - *alpha) / beta;
		*alpha = real_type(1)/(*alpha-beta);
		RNP::TBLAS::Scale( n-1, *alpha, x, incx );
		// If ALPHA is subnormal, it may lose relative accuracy
		while(knt --> 0){ 
			beta *= safmin;
		}
		*alpha = beta;
	}
}

template <class T> // zlarfp, dlarfp, clarfp, slarfp
void GenerateElementaryReflectorPositive(size_t n, T *alpha, T *x, size_t incx, T *tau){
	typedef typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type real_type;
	
	using namespace std;

	// Purpose
	// =======

	// ZLARFP generates a complex elementary reflector H of order n, such
	// that
	//       H' * ( alpha ) = ( beta ),   H' * H = I.
	//            (   x   )   (   0  )
	// where alpha and beta are scalars, beta is real and non-negative, and
	// x is an (n-1)-element complex vector.  H is represented in the form
	//       H = I - tau * ( 1 ) * ( 1 v' )
	//                     ( v )
	// where tau is a complex scalar and v is a complex (n-1)-element
	// vector. Note that H is not hermitian.
	// If the elements of x are all zero and alpha is real, then tau = 0
	// and H is taken to be the unit matrix.
	// Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .

	// Arguments
	// =========

	// N       The order of the elementary reflector.

	// ALPHA   (input/output) COMPLEX*16
	//         On entry, the value alpha.
	//         On exit, it is overwritten with the value beta.

	// X       (input/output) COMPLEX*16 array, dimension
	//                        (1+(N-2)*abs(INCX))
	//         On entry, the vector x.
	//         On exit, it is overwritten with the vector v.

	// INCX    The increment between elements of X. INCX > 0.

	// TAU     The value tau.

	if(n < 1){
		*tau = 0;
		return;
	}

	real_type xnorm = RNP::TBLAS::Norm2(n-1, x, incx);

	if(xnorm == 0){
		// H  =  [1-alpha/abs(alpha) 0; 0 I], sign chosen so alpha >= 0.
		if(RNP::TBLAS::_RealOrComplexChooser<T>::_imag(*alpha) == 0){
			if(RNP::TBLAS::_RealOrComplexChooser<T>::_real(*alpha) >= 0){
				// When tau == 0, the vector is special-cased to be
				// all zeros in the application routines.  We do not need
				// to clear it.
				*tau = 0;
			}else{
				// However, the application routines rely on explicit
				// zero checks when tau != 0, and we must clear X.
				*tau = 2;
				--n; // n not needed anymore
				while(n --> 0){
					*x = 0;
					x += incx;
				}
				*alpha = -*alpha;
			}
		}else{
			// Only "reflecting" the diagonal entry to be real and non-negative.
			xnorm = RNP::TBLAS::_RealOrComplexChooser<T>::_abs(*alpha);
			*tau = real_type(1) - *alpha/xnorm;
			--n; // n not needed anymore
			while(n --> 0){
				*x = 0;
				x += incx;
			}
			*alpha = xnorm;
		}
	}else{ // general case
		real_type amag = RNP::TBLAS::_RealOrComplexChooser<T>::_abs(*alpha);
		real_type beta = RNP::TLASupport::Pythag2(amag, xnorm); if(RNP::TBLAS::_RealOrComplexChooser<T>::_real(*alpha) < 0){ beta = -beta; }
		const real_type safmin = numeric_limits<real_type>::min() / numeric_limits<real_type>::epsilon();
		const real_type rsafmn = 1. / safmin;

		size_t knt = 0;
		if(abs(beta) < safmin){
			// xnorm, beta may be inaccurate; scale x and recompute them
			do{
				++knt;
				RNP::TBLAS::Scale(n-1, rsafmn, x, incx);
				beta *= rsafmn;
				*alpha *= rsafmn;
			}while(abs(beta) < safmin);
			// New beta is at most 1, at least safmin
			xnorm = RNP::TBLAS::Norm2(n-1, x, incx);
			amag = RNP::TBLAS::_RealOrComplexChooser<T>::_abs(*alpha);
			beta = RNP::TLASupport::Pythag2(amag, xnorm); if(RNP::TBLAS::_RealOrComplexChooser<T>::_real(*alpha) < 0){ beta = -beta; }
		}
		*alpha += beta;
		if(beta < 0){
			beta = -beta;
		}else{
			// The following original 3 lines is specific to complex numbers
			//alphr = alphi * (alphi / alpha->real()) + xnorm * (xnorm / alpha->real());
			//*tau = std::complex<double>(alphr / beta, -alphi / beta);
			//*alpha = std::complex<double>(-alphr, alphi);
			
			// Equivalent to:
			// alphr = -( alphi * (alphi / alpha->real()) + xnorm * (xnorm / alpha->real()) );
			// alpha = std::complex<double>(alphr, alphi);
			// tau = -alpha/beta
			//
			// alphr = alpha->real();
			// temp = alphi * (alphi / alphr) + xnorm * (xnorm / alphr);
			// alpha = std::complex<double>(alphr-(temp+alphr), alphi);
			// tau = -alpha/beta
			//
			// Note that: temp+alphr = alphr*[(alpha/alphr)^2 + (xnorm/alphr)^2]
			// 
			real_type alphr = RNP::TBLAS::_RealOrComplexChooser<T>::_real(*alpha);
			real_type ana = RNP::TBLAS::_RealOrComplexChooser<T>::_abs((*alpha / alphr));
			xnorm *= (xnorm/alphr);
			*alpha -= (alphr*ana)*ana;
			*alpha -= xnorm;
		}
		*tau = -*alpha/beta;
		*alpha = RNP::TLASupport::ComplexDivide(T(1), *alpha);
		RNP::TBLAS::Scale(n-1, (*alpha), x, incx);

		// If beta is subnormal, it may lose relative accuracy
		while(knt --> 0){
			beta *= safmin;
		}
		*alpha = beta;
	}
}


template <char side='L', char trans='N'>
struct ApplyOrthognalMatrixFromElementaryReflectors{ // zunmqr, zunm2r, cunmqr, cunm2r
	template <class T>
	ApplyOrthognalMatrixFromElementaryReflectors(size_t m, size_t n, size_t k,
		T *a, int lda, T *tau, T *c, int ldc, T *work)
	{
		using namespace std;

		// Purpose
		// =======

		// ZUNM2R overwrites the general complex m-by-n matrix C with
		//       Q * C  if SIDE = 'L' and TRANS = 'N', or
		//       Q'* C  if SIDE = 'L' and TRANS = 'C', or
		//       C * Q  if SIDE = 'R' and TRANS = 'N', or
		//       C * Q' if SIDE = 'R' and TRANS = 'C',
		// where Q is a complex unitary matrix defined as the product of k
		// elementary reflectors
		//       Q = H(1) H(2) . . . H(k)
		// as returned by ZGEQRF. Q is of order m if SIDE = 'L' and of order n
		// if SIDE = 'R'.

		// Arguments
		// =========

		// SIDE    = 'L': apply Q or Q' from the Left
		//         = 'R': apply Q or Q' from the Right

		// TRANS   = 'N': apply Q  (No transpose)
		//         = 'C': apply Q' (Conjugate transpose)

		// M       The number of rows of the matrix C. M >= 0.

		// N       The number of columns of the matrix C. N >= 0.

		// K       The number of elementary reflectors whose product defines
		//         the matrix Q.
		//         If SIDE = 'L', M >= K >= 0;
		//         if SIDE = 'R', N >= K >= 0.

		// A       (input) COMPLEX*16 array, dimension (LDA,K)
		//         The i-th column must contain the vector which defines the
		//         elementary reflector H(i), for i = 1,2,...,k, as returned by
		//         ZGEQRF in the first k columns of its array argument A.
		//         A is modified by the routine but restored on exit.

		// LDA     The leading dimension of the array A.
		//         If SIDE = 'L', LDA >= max(1,M);
		//         if SIDE = 'R', LDA >= max(1,N).

		// TAU     TAU(i) must contain the scalar factor of the elementary
		//         reflector H(i), as returned by ZGEQRF.

		// C       (input/output) COMPLEX*16 array, dimension (LDC,N)
		//         On entry, the m-by-n matrix C.
		//         On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.

		// LDC     The leading dimension of the array C. LDC >= max(1,M).

		// WORK    (workspace) COMPLEX*16 array, dimension
		//                                  (N) if SIDE = 'L',
		//                                  (M) if SIDE = 'R'

		const bool left = ('L' == side);
		const bool notran = ('N' == trans);

		if(m < 1 || n < 1 || k < 1){ return; }

		int istart, iend, iinc;
		if((left && ! notran) || (! left && notran)){
			istart = 0;
			iend = k;
			iinc = 1;
		}else{
			istart = k-1;
			iend = 0;
			iinc = -1;
		}

		size_t mi, ni;
		size_t ic, jc;
		if(left){
			ni = n;
			jc = 0;
		}else{
			mi = m;
			ic = 0;
		}

		for(int i = istart; iinc < 0 ? i >= iend : i < iend; i += iinc){
			if(left){ // H(i) or H(i)' is applied to C(i:m,1:n)
				mi = m - i;
				ic = i;
			}else{ // H(i) or H(i)' is applied to C(1:m,i:n)
				ni = n - i;
				jc = i;
			}

			// Apply H(i) or H(i)'
			T taui;
			if(notran){
				taui = tau[i];
			}else{
				taui = RNP::TBLAS::_RealOrComplexChooser<T>::_conj(tau[i]);
			}
			T aii = a[i+i*lda];
			a[i+i*lda] = 1;
			RNP::TLASupport::ApplyElementaryReflector<side>(mi, ni, &a[i+i*lda], 1, taui, &c[ic+jc*ldc], ldc, work);
			a[i+i*lda] = aii;
		}
	}
};

template <class T> // zungqr, zung2r, cungqr, cung2r
void GenerateOrthognalMatrixFromElementaryReflectors(size_t m, size_t n, size_t k,
	T *a, int lda, const T *tau, T *work)
{
	using namespace std;

	// ZUNG2R generates an m by n complex matrix Q with orthonormal columns,
	// which is defined as the first n columns of a product of k elementary
	// reflectors of order m
	//       Q  =  H(1) H(2) . . . H(k)
	// as returned by ZGEQRF.

	// Arguments
	// =========

	// M       The number of rows of the matrix Q. M >= 0.

	// N       The number of columns of the matrix Q. M >= N >= 0.

	// K       The number of elementary reflectors whose product defines the
	//         matrix Q. N >= K >= 0.

	// A       (input/output) COMPLEX*16 array, dimension (LDA,N)
	//         On entry, the i-th column must contain the vector which
	//         defines the elementary reflector H(i), for i = 1,2,...,k, as
	//         returned by ZGEQRF in the first k columns of its array
	//         argument A.
	//         On exit, the m by n matrix Q.

	// LDA     (input) int
	//         The first dimension of the array A. LDA >= max(1,M).

	// TAU     (input) COMPLEX*16 array, dimension (K)
	//         TAU(i) must contain the scalar factor of the elementary
	//         reflector H(i), as returned by ZGEQRF.

	// WORK    (workspace) COMPLEX*16 array, dimension (N)

	if(n < 1){ return; }

	// Initialise columns k+1:n to columns of the unit matrix
	for(size_t j = k; j < n; ++j){
		for(size_t l = 0; l < m; ++l){
			a[l+j*lda] = 0;
		}
		a[j+j*lda] = 1;
	}

	for(int i = k-1; i >= 0; --i){
		// Apply H(i) to A(i:m,i:n) from the left
		if((size_t)i < n-1){
			a[i+i*lda] = 1;
			RNP::TLASupport::ApplyElementaryReflector<'L'>(m-i, n-i-1, &a[i+i*lda], 1, tau[i], &a[i+(i+1)*lda], lda, work);
		}
		if((size_t)i < m-1){
			RNP::TBLAS::Scale(m-i-1, -tau[i], &a[i+1+i*lda], 1);
		}
		a[i+i*lda] = T(1) - tau[i];

		// Set A(1:i-1,i) to zero
		for(int l = 0; l < i; ++l){
			a[l+i*lda] = 0;
		}
	}
}

template <class T> // zgeqrf, zgeqr2, dgeqrf, dgeqr2, cgeqrf, cgeqr2, sgeqrf, sgeqr2
void QRFactorization(size_t m, size_t n, T *a, size_t lda, T *tau, T *work){
	using namespace std;

	// ZGEQR2 computes a QR factorization of a complex m by n matrix A = Q * R.

	// Arguments
	// =========

	// M       The number of rows of the matrix A.  M >= 0.

	// N       The number of columns of the matrix A.  N >= 0.

	// A       (input/output) COMPLEX*16 array, dimension (LDA,N)
	//         On entry, the m by n matrix A.
	//         On exit, the elements on and above the diagonal of the array
	//         contain the min(m,n) by n upper trapezoidal matrix R (R is
	//         upper triangular if m >= n); the elements below the diagonal,
	//         with the array TAU, represent the unitary matrix Q as a
	//         product of elementary reflectors (see Further Details).

	// LDA     The leading dimension of the array A.  LDA >= max(1,M).

	// TAU     (output) COMPLEX*16 array, dimension (min(M,N))
	//         The scalar factors of the elementary reflectors (see Further
	//         Details).

	// WORK    (workspace) COMPLEX*16 array, dimension (N)

	// Further Details
	// ===============

	// The matrix Q is represented as a product of elementary reflectors
	//    Q = H(1) H(2) . . . H(k), where k = min(m,n).
	// Each H(i) has the form
	//    H(i) = I - tau * v * v'
	// where tau is a complex scalar, and v is a complex vector with
	// v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
	// and tau in TAU(i).

	size_t k = m; if(n < k){ k = n; }
	for(size_t i = 0; i < k; ++i){
		// Generate elementary reflector H(i) to annihilate A(i+1:m,i)
		size_t row = m-1; if(i+1 < row){ row = i+1; }
		RNP::TLASupport::GenerateElementaryReflector(m-i, &a[i+i*lda], &a[row+i*lda], 1, &tau[i]);
		if(i < n-1){
			// Apply H(i)' to A(i:m,i+1:n) from the left
			std::complex<double> alpha = a[i+i*lda];
			a[i+i*lda] = 1;
			RNP::TLASupport::ApplyElementaryReflector<'L'>(m-i, n-i-1, &a[i+i*lda], 1, std::conj(tau[i]), &a[i+(i+1)*lda], lda, work);
			a[i+i*lda] = alpha;
		}
	}
}

template <class T> // zlartg, dlartg, clartg, slartg
void GeneratePlaneRotation(const T &f, const T &g, typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type *cs, T *sn, T *r)
{
	typedef typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type real_type;
	using namespace std;

	// Purpose
	// =======

	// ZLARTG generates a plane rotation so that

	//    [  CS  SN  ]     [ F ]     [ R ]
	//    [  __      ]  .  [   ]  =  [   ]   where CS**2 + |SN|**2 = 1.
	//    [ -SN  CS  ]     [ G ]     [ 0 ]

	// This is a faster version of the BLAS1 routine ZROTG, except for
	// the following differences:
	//    F and G are unchanged on return.
	//    If G=0, then CS=1 and SN=0.
	//    If F=0, then CS=0 and SN is chosen so that R is real.

	// Arguments
	// =========

	// F       (input) COMPLEX*16
	//         The first component of vector to be rotated.

	// G       (input) COMPLEX*16
	//         The second component of vector to be rotated.

	// CS      (output) DOUBLE PRECISION
	//         The cosine of the rotation.

	// SN      (output) COMPLEX*16
	//         The sine of the rotation.

	// R       (output) COMPLEX*16
	//         The nonzero component of the rotated vector.

	// Further Details
	// ======= =======

	// 3-5-96 - Modified with a new algorithm by W. Kahan and J. Demmel

	const real_type safmin = std::numeric_limits<real_type>::min();
	const real_type eps = std::numeric_limits<real_type>::epsilon();
	const real_type safmn2 = sqrt(safmin/eps);
	const real_type safmx2 = real_type(1) / safmn2;
	
	real_type aif = RNP::TBLAS::_RealOrComplexChooser<T>::_absinf(f);
	real_type aig = RNP::TBLAS::_RealOrComplexChooser<T>::_absinf(g);
	real_type scale = ((aif > aig) ? aif : aig);
	T fs = f;
	T gs = g;
	int count = 0;
	// Edit by vkl: check for inf; otherwise never terminates
	if(scale >= safmx2){
		do{
			++count;
			fs *= safmn2;
			gs *= safmn2;
			scale *= safmn2;
		}while(scale >= safmx2 && scale <= std::numeric_limits<real_type>::max());
	}else if(scale <= safmn2){
		if(T(0) == g){
			*cs = real_type(1);
			*sn = 0;
			*r = f;
			return;
		}
		do{
			--count;
			fs *= safmx2;
			gs *= safmx2;
			scale *= safmx2;
		}while(scale <= safmn2 && scale <= std::numeric_limits<real_type>::max());
	}
	real_type f2 = RNP::TBLAS::_RealOrComplexChooser<T>::_abs2(fs);
	real_type g2 = RNP::TBLAS::_RealOrComplexChooser<T>::_abs2(gs);
	real_type f2limit = g2; if(1 > f2limit){ f2limit = 1; }
	if(f2 <= f2limit * safmin){
		// This is a rare case: F is very small. */
		if(T(0) == f){
			*cs = 0;
			*r = RNP::TBLAS::_RealOrComplexChooser<T>::_abs(g);
			// Do complex/real division explicitly with two real divisions */
			*sn = gs/abs(gs);
			return;
		}
		real_type f2s = abs(fs);
		// G2 and G2S are accurate */
		// G2 is at least SAFMIN, and G2S is at least SAFMN2 */
		real_type g2s = sqrt(g2);
		// Error in CS from underflow in F2S is at most */
		// UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS */
		// If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN, */
		// and so CS .lt. sqrt(SAFMIN) */
		// If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN */
		// and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS) */
		// Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S */
		*cs = f2s / g2s;
		// Make sure abs(FF) = 1 */
		// Do complex/real division explicitly with 2 real divisions */
		T ff(f);
		if(RNP::TBLAS::_RealOrComplexChooser<T>::_absinf(f) <= real_type(1)){
			ff *= safmx2;
		}
		ff /= RNP::TBLAS::_RealOrComplexChooser<T>::_abs(ff);
		*sn = ff*RNP::TBLAS::_RealOrComplexChooser<T>::_conj(gs/g2s);
		*r = *cs*f + *sn*g;
	}else{
		// This is the most common case. */
		// Neither F2 nor F2/G2 are less than SAFMIN */
		// F2S cannot overflow, and it is accurate */
		real_type f2s = sqrt(g2 / f2 + real_type(1));
		// Do the F2S(real)*FS(complex) multiply with two real multiplies */
		*r = f2s*fs;
		*cs = real_type(1) / f2s;
		// Do complex/real division explicitly with two real divisions */
		*sn = (*r/(f2 + g2)) * RNP::TBLAS::_RealOrComplexChooser<T>::_conj(gs);
		while(count > 0){
			*r *= safmx2;
			--count;
		}
		while(count < 0){
			*r *= safmn2;
			++count;
		}
	}
}

template <class T, class TC, class TS> // zrot, drot, crot, srot
void ApplyPlaneRotation(size_t n, T *cx, size_t incx, T *cy, size_t incy, const TC &c, const TS &s){

	// Purpose
	// =======

	// ZROT   applies a plane rotation, where the cos (C) is real and the
	// sin (S) is complex, and the vectors CX and CY are complex.

	// Arguments
	// =========

	// N       The number of elements in the vectors CX and CY.

	// CX      (input/output) COMPLEX*16 array, dimension (N)
	//         On input, the vector X.
	//         On output, CX is overwritten with C*X + S*Y.

	// INCX    The increment between successive values of CY.  INCX <> 0.

	// CY      (input/output) COMPLEX*16 array, dimension (N)
	//         On input, the vector Y.
	//         On output, CY is overwritten with -CONJG(S)*X + C*Y.

	// INCY    The increment between successive values of CY.  INCX <> 0.

	// C       (input) DOUBLE PRECISION
	// S       (input) COMPLEX*16
	//         C and S define a rotation
	//            [  C          S  ]
	//            [ -conjg(S)   C  ]
	//         where C*C + S*CONJG(S) = 1.0.

	for(size_t i = 0; i < n; ++i){
		T stemp = c*(*cx) + s*(*cy);
		*cy = c*(*cy) - std::conj(s)*(*cx);
		*cx = stemp;
		cx += incx;
		cy += incy;
	}
}

template <class T> // zgehd2, dgehd2, cgehd2, sgehd2
void HessenbergReductionUnblocked(size_t n, size_t ilo, size_t ihi, T *a, size_t lda, T *tau, T *work){
	// ZGEHD2 reduces a complex general matrix A to upper Hessenberg form H
	// by a unitary similarity transformation:  Q' * A * Q = H .

	// Arguments
	// =========

	// N       The order of the matrix A.  N >= 0.

	// ILO     (input) INTEGER
	// IHI     (input) INTEGER
	//         It is assumed that A is already upper triangular in rows
	//         and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
	//         set by a previous call to ZGEBAL; otherwise they should be
	//         set to 1 and N respectively. See Further Details.
	//         1 <= ILO <= IHI <= max(1,N).

	// A       (input/output) COMPLEX*16 array, dimension (LDA,N)
	//         On entry, the n by n general matrix to be reduced.
	//         On exit, the upper triangle and the first subdiagonal of A
	//         are overwritten with the upper Hessenberg matrix H, and the
	//         elements below the first subdiagonal, with the array TAU,
	//         represent the unitary matrix Q as a product of elementary
	//         reflectors. See Further Details.

	// LDA     The leading dimension of the array A.  LDA >= max(1,N).

	// TAU     (output) COMPLEX*16 array, dimension (N-1)
	//         The scalar factors of the elementary reflectors (see Further
	//         Details).

	// WORK    (workspace) COMPLEX*16 array, dimension (N)

	// Further Details
	// ===============

	// The matrix Q is represented as a product of (ihi-ilo) elementary
	// reflectors

	//    Q = H(ilo) H(ilo+1) . . . H(ihi-1).

	// Each H(i) has the form

	//    H(i) = I - tau * v * v'

	// where tau is a complex scalar, and v is a complex vector with
	// v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
	// exit in A(i+2:ihi,i), and tau in TAU(i).

	// The contents of A are illustrated by the following example, with
	// n = 7, ilo = 2 and ihi = 6:

	// on entry,                        on exit,

	// ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
	// (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
	// (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
	// (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
	// (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
	// (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
	// (                         a )    (                          a )

	// where a denotes an element of the original matrix A, h denotes a
	// modified element of the upper Hessenberg matrix H, and vi denotes an
	// element of the vector defining H(i).

	for(size_t i = ilo; i < ihi; ++i){
		// Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
		T alpha = a[i+1+i*lda];
		size_t row = i+2; if((int)n-1 < (int)row){ row = n-1; }
		RNP::TLASupport::GenerateElementaryReflector(ihi-i, &alpha, &a[row+i*lda], 1, &tau[i]);
		a[i+1+i*lda] = 1;
		// Apply H(i) to A(1:ihi,i+1:ihi) from the right */
		RNP::TLASupport::ApplyElementaryReflector<'R'>(ihi+1, ihi-i, &a[i+1+i*lda], 1, tau[i], &a[0+(i+1)*lda], lda, work);
		// Apply H(i)' to A(i+1:ihi,i+1:n) from the left */
		RNP::TLASupport::ApplyElementaryReflector<'L'>(ihi-i, n-i-1, &a[i+1+i*lda], 1, std::conj(tau[i]), &a[i+1+(i+1)*lda], lda, work);
		a[i+1+i*lda] = alpha;
	}
}

template <class T> // zlahr2, dlahr2, clahr2, slahr2
void HessenbergReductionBlocked(size_t n, size_t k, size_t nb, T *a, size_t lda, T *tau, T *t, size_t ldt, T *y, size_t ldy){
	if(n <= 1){ return; }
	T ei = 0;
	for(size_t i = 0; i < nb; ++i){
		if(i > 0){
			// Update A(K+1:N,I)
			// Update I-th column of A - Y * V'
			RNP::TBLAS::Conjugate(i, &a[k+i+0*lda], lda);
            RNP::TBLAS::MultMV<'N'>(n-k, i, -1., &y[k+0*ldy], ldy, &a[k+i+0*lda], lda, 1., &a[k+(i-1)*lda], 1);
            RNP::TBLAS::Conjugate(i, &a[k+i+0*lda], lda);
			// Apply I - V * T' * V' to this column (call it b) from the
			// left, using the last column of T as workspace
			//
			// Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
			//          ( V2 )             ( b2 )
			//
			// where V1 is unit lower triangular
			//
			// w := V1' * b1
			RNP::TBLAS::Copy(i, &a[k+i*lda], 1, &t[0+(nb-1)*ldt], 1);
			RNP::TBLAS::MultTrV<'L','C','U'>(i, &a[k+0*lda], lda, &t[0+(nb-1)], 1);
			// w := w + V2'*b2
			RNP::TBLAS::MultMV<'C'>(n-k-i, i, 1., &a[k+i+0*lda], lda, &a[k+i+i*lda], 1, 1., &t[0+(nb-1)*ldt], 1);
			// w := T'*w
			RNP::TBLAS::MultTrV<'U','C','N'>(i, t,ldt, &t[0+(nb-1)*ldt], 1);
			// b2 := b2 - V2*w
			RNP::TBLAS::MultMV<'N'>(n-k-i, i, -1., &a[k+i+0*lda], lda, &t[0+(nb-1)*ldt], 1, 1., &a[k+i+i*lda], 1);
			// b1 := b1 - V1*w
			RNP::TBLAS::MultTrV<'L','N','U'>(i, &a[k+0*lda], lda, &t[0+(nb-1)*ldt], 1);
			RNP::TBLAS::Axpy(i, -1., &t[0+(nb-1)*ldt], 1, &a[k+i*lda], 1);
			a[k+i-1+(i-1)*lda] = ei;
		}
		// Generate the elementary reflector H(I) to annihilate
		// A(K+I+1:N,I)
		size_t arow = (k+i+2 < n ? k+i+1 : n-1);
		GenerateElementaryReflector(n-k-i, &a[k+i+i*lda], &a[arow+i*lda], 1, &tau[i]);
		ei = a[k+i+i*lda];
		a[k+i+i*lda] = 1.;
		// Compute  Y(K+1:N,I)
		RNP::TBLAS::MultMV<'N'>(n-k, n-k-i, 1., &a[k+(i+1)*lda], lda, &a[k+i+i*lda], 1, 0., &y[k+i*ldy], 1);
		RNP::TBLAS::MultMV<'C'>(n-k-i, i, 1., &a[k+i+0*lda], lda, &a[k+i+i*lda], 1, 0., &t[0+i*ldt], 1);
		RNP::TBLAS::MultMV<'N'>(n-k, i, -1., &y[k+0*ldy], ldy, &t[0+i*ldt], 1, 1., &y[k+i*ldy], 1);
		RNP::TBLAS::Scale(n-k, tau[i], &y[k+i*ldy], 1);
		// Compute T(1:I,I)
		RNP::TBLAS::Scale(i, -tau[i], &t[0+i*ldt], 1);
		RNP::TBLAS::MultTrV<'U','N','N'>(i, t, ldt, &t[0+i*ldt], 1);
		t[i+i*ldt] = tau[i];
	}
	a[k+nb-1+(nb-1)*lda] = ei;
	// Compute Y(1:K,1:NB)
	RNP::TBLAS::CopyMatrix<'A'>(k, nb, &a[0+1*lda],lda, y,ldy);
	RNP::TBLAS::MultTrM<'R','L','N','U'>(k, nb, 1., &a[k+0*lda], lda, y, ldy);
	if(n > k+nb){
		RNP::TBLAS::MultMM<'N','N'>(k, nb, n-k-nb, 1., &a[0+(2+nb)*lda], lda, &a[k+nb+0*lda], lda, 1., y, ldy);
	}
	RNP::TBLAS::MultTrM<'R','U','N','N'>(k, nb, 1., t, ldt, y, ldy);
}

template <class T> // zgehrd, dgehrd, cgehrd, sgehrd
void HessenbergReduction(size_t n, size_t ilo, size_t ihi, T *a, size_t lda, T *tau, T *work, size_t lwork){
	//
	//  -- LAPACK routine (version 3.2.1)                                  --
	//  -- LAPACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//  -- April 2009                                                      --
	//
	//
	//  Purpose
	//  =======
	//
	//  ZGEHRD reduces a complex general matrix A to upper Hessenberg form H by
	//  an unitary similarity transformation:  Q' * A * Q = H .
	//
	//  Arguments
	//  =========
	//
	//  N       (input) INTEGER
	//          The order of the matrix A.  N >= 0.
	//
	//  ILO     (input) INTEGER
	//  IHI     (input) INTEGER
	//          It is assumed that A is already upper triangular in rows
	//          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
	//          set by a previous call to ZGEBAL; otherwise they should be
	//          set to 1 and N respectively. See Further Details.
	//          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
	//
	//  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
	//          On entry, the N-by-N general matrix to be reduced.
	//          On exit, the upper triangle and the first subdiagonal of A
	//          are overwritten with the upper Hessenberg matrix H, and the
	//          elements below the first subdiagonal, with the array TAU,
	//          represent the unitary matrix Q as a product of elementary
	//          reflectors. See Further Details.
	//
	//  LDA     (input) INTEGER
	//          The leading dimension of the array A.  LDA >= max(1,N).
	//
	//  TAU     (output) COMPLEX*16 array, dimension (N-1)
	//          The scalar factors of the elementary reflectors (see Further
	//          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to
	//          zero.
	//
	//  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
	//          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
	//
	//  LWORK   (input) INTEGER
	//          The length of the array WORK.  LWORK >= max(1,N).
	//          For optimum performance LWORK >= N*NB, where NB is the
	//          optimal blocksize.
	//
	//          If LWORK = -1, then a workspace query is assumed; the routine
	//          only calculates the optimal size of the WORK array, returns
	//          this value as the first entry of the WORK array, and no error
	//          message related to LWORK is issued by XERBLA.
	//
	//  INFO    (output) INTEGER
	//          = 0:  successful exit
	//          < 0:  if INFO = -i, the i-th argument had an illegal value.
	//
	//  Further Details
	//  ===============
	//
	//  The matrix Q is represented as a product of (ihi-ilo) elementary
	//  reflectors
	//
	//     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
	//
	//  Each H(i) has the form
	//
	//     H(i) = I - tau * v * v'
	//
	//  where tau is a complex scalar, and v is a complex vector with
	//  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
	//  exit in A(i+2:ihi,i), and tau in TAU(i).
	//
	//  The contents of A are illustrated by the following example, with
	//  n = 7, ilo = 2 and ihi = 6:
	//
	//  on entry,                        on exit,
	//
	//  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
	//  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
	//  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
	//  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
	//  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
	//  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
	//  (                         a )    (                          a )
	//
	//  where a denotes an element of the original matrix A, h denotes a
	//  modified element of the upper Hessenberg matrix H, and vi denotes an
	//  element of the vector defining H(i).
	//
	//  This file is a slight modification of LAPACK-3.0's DGEHRD
	//  subroutine incorporating improvements proposed by Quintana-Orti and
	//  Van de Geijn (2006). (See DLAHR2.)
	//
	
	static const size_t nbopt = 32; // block size
	static const size_t nbmin = 2; // minimum block size
	static const size_t nx_parm = 128; // crossover point from blocked to unblocked
	
	T t[nbopt*nbopt];
	const size_t ldt = nbopt;
	
	size_t nb = nbopt;
	
	work[0] = n*nb;
	if((size_t)-1 == lwork){ return; }
	
	for(size_t i = 0; i < ilo; ++i){
		tau[i] = 0.;
	}
	for(size_t i = ihi; i+1 < n; ++i){
		tau[i] = 0.;
	}
	const size_t nh = ihi - ilo + 1;
	if(nh <= 1){ work[0] = 1; return; }
	
	size_t iws = 0;
	
	const size_t nx = nb > nx_parm ? nb : nx_parm;
	if(1 < nb && nb < nh){
		// determine cross over point
		if(nx < nh){
			iws = n*nb;
			if(lwork < iws){ // not enough workspace
				if(lwork >= n*nbmin){
					nb = lwork/n;
				}else{
					nb = 1;
				}
			}
		}
	}

	size_t i = ilo;

	if(nb < nbmin || nb >= nh){
	}else{ // use blocked code
		for(i = ilo; i+nx+1 <= ihi; i += nb){
			size_t ib = (nb < (ihi-i) ? nb : (ihi-i));
			// Reduce columns i:i+ib-1 to Hessenberg form, returning the
			// matrices V and T of the block reflector H = I - V*T*V'
			// which performs the reduction, and also the matrix Y = A*V*T
			HessenbergReductionBlocked(ihi, i, ib, &a[0+i*lda],lda, &tau[i], t, ldt, work, n);
			// Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
			// right, computing  A := A - Y * V'. V(i+ib,ib-1) must be set to 1
			T ei = a[i+ib+(i+ib-1)*lda];
            a[i+ib+(i+ib-1)*lda] = 1.;
            RNP::TBLAS::MultMM<'N','C'>(ihi+1, ihi-i-ib+1+1, ib+1, -1., work,n, &a[i+ib+i*lda],lda, 1., &a[0+(i+ib)*lda],lda);
            a[i+ib+(i+ib-1)*lda] = ei;
			// Apply the block reflector H to A(1:i,i+1:i+ib-1) from the right
			RNP::TBLAS::MultTrM<'R','L','C','U'>(i+1, ib-1+1, 1., &a[i+1+i*lda],lda, work,n);
			for(size_t j = 0; j < ib-1; ++j){
				RNP::TBLAS::Axpy(i+1, -1., &work[n*j],1, &a[0+(i+j)*lda],1);
			}
			// Apply the block reflector H to A(i+1:ihi,i+ib:n) from the left
			ApplyElementaryReflectorBlocked<'L','C','F','C'>(ihi-i+1, n-i-ib+1, ib, &a[i+1+i*lda], lda, t, ldt, &a[i+1+(i+ib)*lda], lda, work, n);
		}
	}
	
	HessenbergReductionUnblocked(n, i, ihi, a, lda, tau, work);
	work[0] = iws;
}

template <char compq='N', char compz='N'>
struct GeneralizedHessenbergReduction{ // zgghrd, dgghrd, cgghrd, sgghrd
	template <class T>
	GeneralizedHessenbergReduction(size_t n, size_t ilo, size_t ihi,
		T *a, size_t lda, T *b, size_t ldb, T *q, size_t ldq, T *z, size_t ldz)
	{
		using namespace std;

		// Purpose
		// =======

		// ZGGHRD reduces a pair of complex matrices (A,B) to generalized upper
		// Hessenberg form using unitary transformations, where A is a
		// general matrix and B is upper triangular.  The form of the
		// generalized eigenvalue problem is
		//    A*x = lambda*B*x,
		// and B is typically made upper triangular by computing its QR
		// factorization and moving the unitary matrix Q to the left side
		// of the equation.

		// This subroutine simultaneously reduces A to a Hessenberg matrix H:
		//    Q**H*A*Z = H
		// and transforms B to another upper triangular matrix T:
		//    Q**H*B*Z = T
		// in order to reduce the problem to its standard form
		//    H*y = lambda*T*y
		// where y = Z**H*x.

		// The unitary matrices Q and Z are determined as products of Givens
		// rotations.  They may either be formed explicitly, or they may be
		// postmultiplied into input matrices Q1 and Z1, so that
		//      Q1 * A * Z1**H = (Q1*Q) * H * (Z1*Z)**H
		//      Q1 * B * Z1**H = (Q1*Q) * T * (Z1*Z)**H
		// If Q1 is the unitary matrix from the QR factorization of B in the
		// original equation A*x = lambda*B*x, then ZGGHRD reduces the original
		// problem to generalized Hessenberg form.

		// Arguments
		// =========

		// COMPQ   (input) CHARACTER*1
		//         = 'N': do not compute Q;
		//         = 'I': Q is initialized to the unit matrix, and the
		//                unitary matrix Q is returned;
		//         = 'V': Q must contain a unitary matrix Q1 on entry,
		//                and the product Q1*Q is returned.

		// COMPZ   (input) CHARACTER*1
		//         = 'N': do not compute Q;
		//         = 'I': Q is initialized to the unit matrix, and the
		//                unitary matrix Q is returned;
		//         = 'V': Q must contain a unitary matrix Q1 on entry,
		//                and the product Q1*Q is returned.

		// N       (input) INTEGER
		//         The order of the matrices A and B.  N >= 0.

		// ILO     (input) INTEGER
		// IHI     (input) INTEGER
		//         ILO and IHI mark the rows and columns of A which are to be
		//         reduced.  It is assumed that A is already upper triangular
		//         in rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI are
		//         normally set by a previous call to ZGGBAL; otherwise they
		//         should be set to 1 and N respectively.
		//         1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.

		// A       (input/output) COMPLEX*16 array, dimension (LDA, N)
		//         On entry, the N-by-N general matrix to be reduced.
		//         On exit, the upper triangle and the first subdiagonal of A
		//         are overwritten with the upper Hessenberg matrix H, and the
		//         rest is set to zero.

		// LDA     (input) INTEGER
		//         The leading dimension of the array A.  LDA >= max(1,N).

		// B       (input/output) COMPLEX*16 array, dimension (LDB, N)
		//         On entry, the N-by-N upper triangular matrix B.
		//         On exit, the upper triangular matrix T = Q**H B Z.  The
		//         elements below the diagonal are set to zero.

		// LDB     (input) INTEGER
		//         The leading dimension of the array B.  LDB >= max(1,N).

		// Q       (input/output) COMPLEX*16 array, dimension (LDQ, N)
		//         On entry, if COMPQ = 'V', the unitary matrix Q1, typically
		//         from the QR factorization of B.
		//         On exit, if COMPQ='I', the unitary matrix Q, and if
		//         COMPQ = 'V', the product Q1*Q.
		//         Not referenced if COMPQ='N'.

		// LDQ     (input) INTEGER
		//         The leading dimension of the array Q.
		//         LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.

		// Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)
		//         On entry, if COMPZ = 'V', the unitary matrix Z1.
		//         On exit, if COMPZ='I', the unitary matrix Z, and if
		//         COMPZ = 'V', the product Z1*Z.
		//         Not referenced if COMPZ='N'.

		// LDZ     (input) INTEGER
		//         The leading dimension of the array Z.
		//         LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.

		// Further Details
		// ===============

		// This routine reduces A to Hessenberg and B to triangular form by
		// an unblocked reduction, as described in _Matrix_Computations_,
		// by Golub and van Loan (Johns Hopkins Press).

		const bool ilq = ((compq == 'V') || (compq == 'I'));
		const bool ilz = ((compz == 'V') || (compz == 'I'));

		// Initialize Q and Z if desired.
		if(compq == 'I'){
			RNP::TBLAS::SetMatrix<'F'>(n, n, 0, 1, q, ldq);
		}
		if(compz == 'I'){
			RNP::TBLAS::SetMatrix<'F'>(n, n, 0, 1, z, ldz);
		}

		if(n <= 1){ return; }

		// Zero out lower triangle of B
		for(size_t jcol = 0; jcol < n-1; ++jcol){
			for(size_t jrow = jcol+1; jrow < n; ++jrow){
				b[jrow+jcol*ldb] = 0;
			}
		}

		// Reduce A and B
		for(size_t jcol = ilo; (int)jcol <= (int)ihi - 2; ++jcol){
			for(size_t jrow = ihi; jrow >= jcol + 2; --jrow){
				typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type c;
				T s, ctemp;
				
				// Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL)
				ctemp = a[jrow-1+jcol*lda];
				RNP::TLASupport::GeneratePlaneRotation(ctemp, a[jrow+jcol*lda], &c, &s, &a[jrow-1+jcol*lda]);
				a[jrow+jcol*lda] = 0;
				RNP::TLASupport::ApplyPlaneRotation(n-jcol-1, &a[jrow-1+(jcol+1) * lda], lda, &a[jrow+(jcol+1)*lda], lda, c, s);
				RNP::TLASupport::ApplyPlaneRotation(n+1-jrow, &b[jrow-1+(jrow-1)*ldb], ldb, &b[jrow+(jrow-1)*ldb], ldb, c, s);
				if(ilq){
					RNP::TLASupport::ApplyPlaneRotation(n, &q[0+(jrow-1)*ldq], 1, &q[0+jrow*ldq], 1, c, std::conj(s));
				}

				// Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1)
				ctemp = b[jrow+jrow*ldb];
				RNP::TLASupport::GeneratePlaneRotation(ctemp, b[jrow+(jrow-1)*ldb], &c, &s, &b[jrow+jrow*ldb]);
				b[jrow+(jrow-1)*ldb] = 0;
				RNP::TLASupport::ApplyPlaneRotation(ihi+1, &a[0+jrow*lda], 1, &a[0+(jrow-1)*lda], 1, c, s);
				RNP::TLASupport::ApplyPlaneRotation(jrow, &b[0+jrow*ldb], 1, &b[0+(jrow-1)*ldb], 1, c, s);
				if(ilz){
					RNP::TLASupport::ApplyPlaneRotation(n, &z[0+jrow*ldz], 1, &z[0+(jrow-1)*ldz], 1, c, s);
				}
			}
		}
	}
};

template <class T> // zgetrf, dgetrf, cgetrf, sgetrf
int LUDecomposition(size_t m, size_t n, T *a, size_t lda, size_t *pivots){
	int info = 0;
	size_t min_dim = (m < n ? m : n);
	for(size_t j = 0; j < min_dim; ++j){
		size_t jp = j + RNP::TBLAS::MaximumIndex(m-j, &a[j+j*lda], 1);
		pivots[j] = jp;
		if(T(0) != a[jp+j*lda]){
			if(jp != j){
				RNP::TBLAS::Swap(n, &a[j+0*lda], lda, &a[jp+0*lda], lda);
			}
			if(j < m){
				RNP::TBLAS::Scale(m-j-1, T(1)/a[j+j*lda], &a[j+1+j*lda], 1); // possible overflow when inverting A(j,j)
			}
		}else{
			info = j;
		}
		if(j < min_dim){
			RNP::TBLAS::Rank1Update(m-j-1, n-j-1, T(-1), &a[j+1+j*lda], 1, &a[j+(j+1)*lda], lda, &a[j+1+(j+1)*lda], lda);
		}
	}
	if(0 != info){ return info+1; }
	return 0;
}


template <char trans='N'>
struct LUSolve{
	template <class T>
	LUSolve(size_t n, size_t nRHS, T *a, size_t lda, size_t *ipiv, T *b, size_t ldb){
		if(0 == n || nRHS == 0){ return; }
		
		if(trans == 'N'){
			for(size_t i = 0; i < n; ++i){
				if(ipiv[i] != i){
					RNP::TBLAS::Swap(nRHS, &b[i+0*ldb], ldb, &b[ipiv[i]+0*ldb], ldb);
				}
			}
			RNP::TBLAS::SolveTrM<'L','L','N','U'>(n, nRHS, T(1), a, lda, b, ldb);
			RNP::TBLAS::SolveTrM<'L','U','N','N'>(n, nRHS, T(1), a, lda, b, ldb);
		}else{
			if(trans == 'T'){
				RNP::TBLAS::SolveTrM<'L','U','T','N'>(n, nRHS, T(1), a, lda, b, ldb);
				RNP::TBLAS::SolveTrM<'L','L','T','U'>(n, nRHS, T(1), a, lda, b, ldb);
			}else if(trans == 'C'){
				RNP::TBLAS::SolveTrM<'L','U','C','N'>(n, nRHS, T(1), a, lda, b, ldb);
				RNP::TBLAS::SolveTrM<'L','L','C','U'>(n, nRHS, T(1), a, lda, b, ldb);
			}
			for(int i = (int)n-1; i >= 0; --i){
				if((int)ipiv[i] != i){
					RNP::TBLAS::Swap(nRHS, &b[i+0*ldb], ldb, &b[ipiv[i]+0*ldb], ldb);
				}
			}
		}
	}
};

template <class T>
void Determinant(size_t n, T *a, size_t lda, T *mant, typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type *base, int *expo, size_t *pivots = NULL){
	typedef typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type real_type;
	static const real_type b = 16.;
	size_t *ipiv = pivots;
	if(NULL == pivots){
		ipiv = new size_t[n];
	}
	LUDecomposition(n,n, a, lda, ipiv);
	*base = b;
	*expo = 0;
	*mant = T(1);
	for(size_t i = 0; i < n; ++i){
		*mant *= a[i+i*lda];
		if(T(0) == *mant){ break; }
		while(RNP::TBLAS::_RealOrComplexChooser<T>::_abs(*mant) < 1){
			*mant *= b;
			(*expo)--;
		}
		while(RNP::TBLAS::_RealOrComplexChooser<T>::_abs(*mant) >= b){
			*mant /= b;
			(*expo)++;
		}
		if(ipiv[i] != i){ *mant = -(*mant); }
	}
	if(NULL == pivots){
		delete [] ipiv;
	}
}

// if inv is true, the permuations are applied in reverse order (applies the inverse permutation
template <char side, char inv>
struct ApplyPermutations{
	template <class T>
	ApplyPermutations(size_t m, size_t n, T *a, size_t lda, const size_t *pivots){
		if('L' == side){
			if('N' == inv){
				for(size_t i = 0; i < m; ++i){
					size_t ip = pivots[i];
					if(i == ip){ continue; }
					RNP::TBLAS::Swap(n, &a[i+0*lda],lda, &a[ip+0*lda],lda);
				}
			}else{
				size_t i = m;
				while(i --> 0){
					size_t ip = pivots[i];
					if(i == ip){ continue; }
					RNP::TBLAS::Swap(n, &a[i+0*lda],lda, &a[ip+0*lda],lda);
				}
			}
		}else{
			if('N' == inv){
				for(size_t j = 0; j < n; ++j){
					size_t jp = pivots[j];
					if(j == jp){ continue; }
					RNP::TBLAS::Swap(m, &a[0+j*lda],1, &a[0+jp*lda],1);
				}
			}else{
				size_t j = n;
				while(j --> 0){
					size_t jp = pivots[j];
					if(j == jp){ continue; }
					RNP::TBLAS::Swap(m, &a[0+j*lda],1, &a[0+jp*lda],1);
				}
			}
		}
	}
};


/*
template <class T>
void Transpose(size_t m, size_t n, T *a, size_t lda, size_t lwork = 0, size_t *work = NULL){
	size_t *move = work;
	size_t stack_move = 0;
	size_t nmove = lwork;
	if(0 == lwork || NULL == work){
		move = &stack_move;
		nmove = 1;
	}
	nmove *= 8*sizeof(size_t);
	
	if(m < 2 || n < 2){ return; }
	if(m == n){
		for(size_t i = 0; i < n-1; ++i){
			for(size_t j = i+1; j < n; ++j){
				std::swap(a[i+j*lda], a[j+i*lda]);
			}
		}
		return;
	}
	size_t ncount = 2;
	size_t k = m*n-1;
	for(size_t i = 0; i < nmove; ++i){
		move[i] = 0;
	}
}*/

}; // namespace TLASupport
}; // namespace RNP



/*
      IF (M.LT.3 .OR. N.LT.3) GO TO 30
C CALCULATE THE NUMBER OF FIXED POINTS, EUCLIDS ALGORITHM
C FOR GCD(M-1,N-1).
      IR2 = M - 1
      IR1 = N - 1
   20 IR0 = MOD(IR2,IR1)
      IR2 = IR1
      IR1 = IR0
      IF (IR0.NE.0) GO TO 20
      NCOUNT = NCOUNT + IR2 - 1
C SET INITIAL VALUES FOR SEARCH
   30 I = 1
      IM = M
C AT LEAST ONE LOOP MUST BE RE-ARRANGED
      GO TO 80
C SEARCH FOR LOOPS TO REARRANGE
   40 MAX = K - I
      I = I + 1
      IF (I.GT.MAX) GO TO 160
      IM = IM + M
      IF (IM.GT.K) IM = IM - K
      I2 = IM
      IF (I.EQ.I2) GO TO 40
      IF (I.GT.IWRK) GO TO 60
      IF (MOVE(I).EQ.0) GO TO 80
      GO TO 40
   50 I2 = M*I1 - K*(I1/N)
   60 IF (I2.LE.I .OR. I2.GE.MAX) GO TO 70
      I1 = I2
      GO TO 50
   70 IF (I2.NE.I) GO TO 40
C REARRANGE THE ELEMENTS OF A LOOP AND ITS COMPANION LOOP
   80 I1 = I
      KMI = K - I
      B = A(I1+1)
      I1C = KMI
      C = A(I1C+1)
   90 I2 = M*I1 - K*(I1/N)
      I2C = K - I2
      IF (I1.LE.IWRK) MOVE(I1) = 2
      IF (I1C.LE.IWRK) MOVE(I1C) = 2
      NCOUNT = NCOUNT + 2
      IF (I2.EQ.I) GO TO 110
      IF (I2.EQ.KMI) GO TO 100
      A(I1+1) = A(I2+1)
      A(I1C+1) = A(I2C+1)
      I1 = I2
      I1C = I2C
      GO TO 90
C FINAL STORE AND TEST FOR FINISHED
  100 D = B
      B = C
      C = D
  110 A(I1+1) = B
      A(I1C+1) = C
      IF (NCOUNT.LT.MN) GO TO 40
C NORMAL RETURN
  120 IOK = 0
      RETURN

      GO TO 120
C ERROR RETURNS.
  160 IOK = I
  170 RETURN
  180 IOK = -1
      GO TO 170
  190 IOK = -2
      GO TO 170
      END
*/
#endif // _RNP_TLASUPPORT_H_
