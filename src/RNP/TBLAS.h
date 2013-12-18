#ifndef _RNP_TBLAS_H_
#define _RNP_TBLAS_H_

/* BLAS license:

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

Note that much of the low level BLAS was not copied from the official
BLAS distribution, and is in the public domain.

*/

#define _USE_MATH_DEFINES
#include <complex>
#include <cmath>

/*
 * Preprocessor flags:
 *   RNP_HAVE_BLAS
 *   RNP_TBLAS_USE_RANDOM
 */

namespace RNP{
namespace TBLAS{

// Classes to differentiate complex types
template <class T>
struct _RealOrComplexChooser{
	typedef T value_type;
	typedef T real_type;

	inline static value_type _conj(const value_type &v){ return v; }
	inline static real_type _real(const value_type &v){ return v; }
	inline static real_type _imag(const value_type &v){ return 0; }
	// abs is euclidean norm, abs2 is euclidean norm squared
	// abs1 is 1-norm (sum of absolute values)
	// absinf is infinity-norm (max absolute value)
	inline static real_type _abs(const value_type &v){ using namespace std; return abs(v); }
	inline static real_type _abs2(const value_type &v){ return v*v; }
	inline static real_type _abs1(const value_type &v){ using namespace std; return abs(v); }
	inline static real_type _absinf(const value_type &v){ using namespace std; return abs(v); }
};
template <class T>
struct _RealOrComplexChooser<std::complex<T> >{
	typedef std::complex<T> value_type;
	typedef typename std::complex<T>::value_type real_type;

	inline static value_type _conj(const value_type &v){ using namespace std; return conj(v); }
	inline static real_type _real(const value_type &v){ return v.real(); }
	inline static real_type _imag(const value_type &v){ return v.imag(); }
	inline static real_type _abs(const value_type &v){ using namespace std; return abs(v); }
	inline static real_type _abs2(const value_type &v){ using namespace std; return norm(v); }
	inline static real_type _abs1(const value_type &v){ using namespace std; return abs(v.real()) + abs(v.imag()); }
	inline static real_type _absinf(const value_type &v){
		using namespace std;
		real_type ar = abs(v.real());
		real_type ai = abs(v.imag());
		return ((ar > ai) ? ar : ai);
	}
};

template <class TV, class T>
void Fill(size_t n, const TV &value, T *x, size_t incx){
	while(n --> 0){
		*x = value;
		x += incx;
	}
}

//// Level 1

template <class T>
void Swap(size_t n, T *x, size_t incx, T *y, size_t incy){
	while(n --> 0){
		std::swap(*x, *y);
		x += incx; y += incy;
	}
}

template <class S, class T>
void Scale(size_t n, const S &scale, T *x, size_t incx){
	while(n --> 0){
		(*x) *= scale;
		x += incx;
	}
}
template <class S, class T>
void Scale(size_t m, size_t n, const S &scale, T *a, size_t lda){
	for(size_t j = 0; j < n; ++j){
		Scale(m, scale, &a[0+j*lda], 1);
	}
}

template <class T>
void Copy(size_t n, const T *src, size_t incsrc, T *dst, size_t incdst){
	while(n --> 0){
		(*dst) = (*src);
		src += incsrc; dst += incdst;
	}
}

template <class S, class T>
void Axpy(size_t n, const S &alpha, const T *x, size_t incx, T *y, size_t incy){
	while(n --> 0){
		(*y) += alpha*(*x);
		x += incx; y += incy;
	}
}

template <class S, class T>
void Axpy(size_t m, size_t n, const S &alpha, const T *x, size_t ldx, T *y, size_t ldy){
	for(size_t j = 0; j < n; ++j){
		Axpy(m, alpha, &x[0+j*ldx], 1, &y[0+j*ldy], 1);
	}
}

template <class T>
T Dot(size_t n, const T *x, size_t incx, const T *y, size_t incy){
	T sum(0);
	while(n --> 0){
		sum += (*x)*(*y);
		x += incx; y += incy;
	}
	return sum;
}

template <class T>
T ConjugateDot(size_t n, const T *x, size_t incx, const T *y, size_t incy){
	T sum(0);
	while(n --> 0){
		sum += _RealOrComplexChooser<T>::_conj(*x)*(*y);
		x += incx; y += incy;
	}
	return sum;
}

template <class T>
typename _RealOrComplexChooser<T>::real_type Norm2(size_t n, const T *x, size_t incx){
	typedef typename _RealOrComplexChooser<T>::real_type real_type;
	using namespace std;
	real_type ssq(1), scale(0);
	while(n --> 0){
		if(0 != _RealOrComplexChooser<T>::_real(*x)){
			real_type temp = abs(_RealOrComplexChooser<T>::_real(*x));
			if(scale < temp){
				real_type r = scale/temp;
				ssq = ssq*r*r + real_type(1);
				scale = temp;
			}else{
				real_type r = temp/scale;
				ssq += r*r;
			}
		}
		if(0 != _RealOrComplexChooser<T>::_imag(*x)){
			real_type temp = abs(_RealOrComplexChooser<T>::_imag(*x));
			if(scale < temp){
				real_type r = scale/temp;
				ssq = ssq*r*r + real_type(1);
				scale = temp;
			}else{
				real_type r = temp/scale;
				ssq += r*r;
			}
		}
		x += incx;
	}
	using namespace std;
	return scale*sqrt(ssq);
}

template <class T>
typename _RealOrComplexChooser<T>::real_type Asum(size_t n, const T *x, size_t incx){
	typedef typename _RealOrComplexChooser<T>::real_type real_type;
	real_type sum(0);
	while(n --> 0){
		sum += _RealOrComplexChooser<T>::_abs1(*x);
		x += incx;
	}
	return sum;
}

template <class T>
size_t MaximumIndex(size_t n, const T *x, size_t incx){
	if(n < 1){ return 0; }
	typedef typename _RealOrComplexChooser<T>::real_type real_type;
	real_type mv = _RealOrComplexChooser<T>::_abs1(*x);
	size_t mi = 0;
	for(size_t i = 1; i < n; ++i){
		x += incx;
		real_type cv = _RealOrComplexChooser<T>::_abs1(*x);
		if(cv > mv){ mi = i; mv = cv; }
	}
	return mi;
}

//// Level 2

template <char trans>
struct MultMV{ // zgemv, dgemv, cgemv, sgemv
	template <class A, class B, class T>
	MultMV(size_t m, size_t n, const A &alpha, const T *a, size_t lda,
           const T *x, size_t incx,
           const B &beta, T *y, size_t incy)
	{
		if(m < 1 || n < 1 || (A(0) == alpha && B(1) == beta)){ return; }
		const bool noconj = (trans == 'T');
		const size_t leny = ((trans == 'N') ? m : n);
		
		if(B(1) != beta){
			T *yy = y;
			size_t l = leny;
			if(B(0) == beta){
				while(l --> 0){
					*yy = 0;
					yy += incy;
				}
			}else{
				while(l --> 0){
					*yy *= beta;
					yy += incy;
				}
			}
		}
		if(A(0) == alpha){ return; }
		if('N' == trans){ // y <- alpha*A*x + y
			const T *xx = x;
			for(size_t j = 0; j < n; ++j){
				if(T(0) != *xx){
					T temp(alpha*(*xx));
					T *yy = y;
					for(size_t i = 0; i < m; ++i){
						*yy += temp*a[i+j*lda];
						yy += incy;
					}
				}
				xx += incx;
			}
		}else{
			for(size_t j = 0; j < n; ++j){
				T temp(0);
				if(noconj){
					const T *xx = x;
					for(size_t i = 0; i < m; ++i){
						temp += a[i+j*lda]*(*xx);
						xx += incx;
					}
				}else{
					const T *xx = x;
					for(size_t i = 0; i < m; ++i){
						temp += _RealOrComplexChooser<T>::_conj(a[i+j*lda])*(*xx);
						xx += incx;
					}
				}
				*y += alpha*temp;
				y += incy;
			}
		}
	}
};

template <char uplo, char trans, char diag>
struct MultTrV{ // ztrmv, dtrmv, ctrmv, strmv
	template <class T>
	MultTrV(size_t n, const T *a, size_t lda, T *x, size_t incx){
		if(n < 1){ return; }
		const bool noconj = ('T' == trans);
		const bool nounit = ('N' == trans);
		
		if('N' == trans){ // x <- A*x
			T *xx = x;
			if('U' == uplo){
				for(size_t j = 0; j < n; ++j){
					if(T(0) != *xx){
						T temp(*xx);
						T *x2 = x;
						for(size_t i = 0; i < j; ++i){
							*x2 += temp*a[i+j*lda];
							x2 += incx;
						}
						if(nounit){ *xx *= a[j+j*lda]; }
					}
					xx += incx;
				}
			}else{
				size_t j = n;
				xx += (n-1)*incx;
				T *x2i = xx;
				while(j --> 0){
					if(T(0) != *xx){
						T temp(*xx);
						T *x2 = x2i;
						size_t i = n;
						while(i --> j+1){
							*x2 += temp*a[i+j*lda];
							x2 -= incx;
						}
						if(nounit){ *xx *= a[j+j*lda]; }
					}
					xx -= incx;
				}
			}
		}else{ // x <- A'*x or x <- conj(A')*x;
			if('U' == uplo){
				T *xx = x + (n-1)*incx;
				size_t j = n;
				while(j --> 0){
					T temp(*xx);
					T *x2 = xx;
					if(noconj){
						if(nounit){ temp *= a[j+j*lda]; }
						for(size_t i = (j > 0 ? j-1 : 0); i-- > 0; ){
							x2 -= incx;
							temp += a[i+j*lda]*(*x2);
						}
					}else{
						if(nounit){ temp *= _RealOrComplexChooser<T>::_conj(a[j+j*lda]); }
						for(size_t i = (j > 0 ? j-1 : 0); i-- > 0; ){
							x2 -= incx;
							temp += _RealOrComplexChooser<T>::_conj(a[i+j*lda])*(*x2);
						}
					}
					*xx = temp;
					xx -= incx;
				}
			}else{
				T *xx = x;
				for(size_t j = 0; j < n; ++j){
					T temp(*xx);
					T *x2 = xx;
					if(noconj){
						if(nounit){ temp *= a[j+j*lda]; }
						for(size_t i = j+1; i < n; ++i){
							x2 += incx;
							temp += a[i+j*lda]*(*x2);
						}
					}else{
						if(nounit){ temp *= _RealOrComplexChooser<T>::_conj(a[j+j*lda]); }
						for(size_t i = j+1; i < n; ++i){
							x2 += incx;
							temp += _RealOrComplexChooser<T>::_conj(a[i+j*lda])*(*x2);
						}
					}
					*xx = temp;
					xx += incx;
				}
			}
		}
	}
};


template <char uplo>
struct MultHermV{ // zhemv, chemv, dsymv, ssymv
	template <class TA, class TB, class T>
	MultHermV(size_t n, const TA &alpha, const T *a, size_t lda, const T *x, size_t incx, const TB &beta, T *y, size_t incy){
		// 
		// ZHEMV  performs the matrix-vector  operation
		//   y := alpha*A*x + beta*y,
		// where alpha and beta are scalars, x and y are n element vectors and
		// A is an n by n hermitian matrix.
		// 
		// Arguments
		// ==========
		// 
		// UPLO
		// On entry, UPLO specifies whether the upper or lower
		// triangular part of the array A is to be referenced as
		// follows:
		// 
		// UPLO = 'U' or 'u'   Only the upper triangular part of A
		// is to be referenced.
		// 
		// UPLO = 'L' or 'l'   Only the lower triangular part of A
		// is to be referenced.
		// 
		// N      - int.
		// On entry, N specifies the order of the matrix A.
		// N must be at least zero.
		// 
		// ALPHA  - std::complex<double>      .
		// On entry, ALPHA specifies the scalar alpha.
		// 
		// A      - std::complex<double>       array of DIMENSION ( LDA, n ).
		// Before entry with  UPLO = 'U' or 'u', the leading n by n
		// upper triangular part of the array A must contain the upper
		// triangular part of the hermitian matrix and the strictly
		// lower triangular part of A is not referenced.
		// Before entry with UPLO = 'L' or 'l', the leading n by n
		// lower triangular part of the array A must contain the lower
		// triangular part of the hermitian matrix and the strictly
		// upper triangular part of A is not referenced.
		// Note that the imaginary parts of the diagonal elements need
		// not be set and are assumed to be zero.
		// 
		// LDA    - int.
		// On entry, LDA specifies the first dimension of A as declared
		// in the calling (sub) program. LDA must be at least
		// max( 1, n ).
		// Unchanged on exit.
		// 
		// X      - std::complex<double>       array of dimension at least
		// ( 1 + ( n - 1 )*abs( INCX ) ).
		// Before entry, the incremented array X must contain the n
		// element vector x.
		// 
		// INCX   - int.
		// On entry, INCX specifies the increment for the elements of
		// X. INCX must not be zero.
		// 
		// BETA   - std::complex<double>      .
		// On entry, BETA specifies the scalar beta. When BETA is
		// supplied as zero then Y need not be set on input.
		// 
		// Y      - std::complex<double>       array of dimension at least
		// ( 1 + ( n - 1 )*abs( INCY ) ).
		// Before entry, the incremented array Y must contain the n
		// element vector y. On exit, Y is overwritten by the updated
		// vector y.
		// 
		// INCY   - int.
		// On entry, INCY specifies the increment for the elements of
		// Y. INCY must not be zero.
		// 
		if((n == 0) || ((alpha == TA(0)) && (beta == TB(1)))){ return; }
		if(beta != TB(1)){ // First form  y := beta*y.
			T *yy = y;
			if(beta == TB(0)){
				for(size_t i = 0; i < n; ++i){
					*yy = 0;
					yy += incy;
				}
			}else{
				for(size_t i = 0; i < n; ++i){
					*yy *= beta;
					yy += incy;
				}
			}
		}
		if(alpha == TA(0)) return;
		if(uplo == 'U'){ // Form y when A is stored in upper triangle.
			const T *xj = x;
			T *yj = y;
			for(size_t j =  0; j < n; ++j){
				T temp1 = alpha*(*xj);
				T temp2 = 0;
				const T *xi = x;
				T *yi = y;
				for(size_t i = 0; i < j; ++i){
					*yi += temp1*a[i+j*lda];
					temp2 += _RealOrComplexChooser<T>::_conj(a[i+j*lda])*(*xi);
					xi += incx;
					yi += incy;
				}
				*yj += temp1*_RealOrComplexChooser<T>::_real(a[j+j*lda]) + alpha*temp2;
				xj += incx;
				yj += incy;
			}
		}else{ // Form y when A is stored in lower triangle.
			const T *xj = x;
			T *yj = y;
			for(size_t j = 0; j < n; ++j){
				T temp1 = alpha*(*xj);
				T temp2 = 0;
				*yj += temp1*_RealOrComplexChooser<T>::_real(a[j+j*lda]);
				const T *xi = xj+1;
				T *yi = yj+1;
				for(size_t i = j+1; i < n; ++i){
					*yi += temp1*a[i+j*lda];
					temp2 += _RealOrComplexChooser<T>::_conj(a[i+j*lda])*(*xi);
					xi += incx;
					yi += incy;
				}
				*yj += alpha*temp2;
				xj += incx;
				yj += incy;
			}
		}
	}
};

template <char uplo, char trans, char diag>
struct SolveTrV{ // ztrsv, dtrsv, ctrsv, strsv
	template <class T>
	SolveTrV(size_t n, const T *a, size_t lda, T *x, size_t incx){
		// ZTRSV  solves one of the systems of equations
		//    A*x = b,   or   A'*x = b,   or   conjg( A' )*x = b,
		// where b and x are n element vectors and A is an n by n unit, or
		// non-unit, upper or lower triangular matrix.

		// No test for singularity or near-singularity is included in this
		// routine. Such tests must be performed before calling this routine.

		// Arguments
		// ==========

		// UPLO   - whether the matrix is an upper or
		//          lower triangular matrix as follows:
		//             UPLO = 'U' or 'u'   A is an upper triangular matrix.
		//             UPLO = 'L' or 'l'   A is a lower triangular matrix.

		// TRANS  - specifies the equations to be solved as follows:
		//             TRANS = 'N' or 'n'   A*x = b.
		//             TRANS = 'T' or 't'   A'*x = b.
		//             TRANS = 'C' or 'c'   conjg( A' )*x = b.

		// DIAG   - whether or not A is unit triangular as follows:
		//             DIAG = 'U' or 'u'   A is assumed to be unit triangular.
		//             DIAG = 'N' or 'n'   A is not assumed to be unit
		//                                 triangular.

		// N      - the order of the matrix A. N must be at least zero.

		// A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
		//          Before entry with  UPLO = 'U' or 'u', the leading n by n
		//          upper triangular part of the array A must contain the upper
		//          triangular matrix and the strictly lower triangular part of
		//          A is not referenced.
		//          Before entry with UPLO = 'L' or 'l', the leading n by n
		//          lower triangular part of the array A must contain the lower
		//          triangular matrix and the strictly upper triangular part of
		//          A is not referenced.
		//          Note that when  DIAG = 'U' or 'u', the diagonal elements of
		//          A are not referenced either, but are assumed to be unity.

		// LDA    - the first dimension of A as declared
		//          in the calling (sub) program. LDA must be at least
		//          max( 1, n ).

		// X      - COMPLEX*16       array of dimension at least
		//          ( 1 + ( n - 1 )*abs( INCX ) ).
		//          Before entry, the incremented array X must contain the n
		//          element right-hand side vector b. On exit, X is overwritten
		//          with the solution vector x.

		// INCX   - the increment for the elements of X. INCX must not be zero.
/*
		// Parameter adjustments
		int a_offset = 1 + lda;
		a -= a_offset;
		
		--x;
		int i, j, ix, jx, kx;

		if(n == 0){ return; }

		const bool noconj = (trans == 'T');
		const bool nounit = (diag == 'N');

		// Set up the start point in X if the increment is not unity. This
		// will be  ( N - 1 )*INCX  too small for descending loops.

		if(incx <= 0){
			kx = 1 - (n - 1) * incx;
		}else{
			kx = 1;
		}

		// Start the operations. In this version the elements of A are
		// accessed sequentially with one pass through A.

		if(trans == 'N'){ // Form  x := inv( A )*x.
			if(uplo == 'U'){
				jx = kx + (n - 1) * incx;
				for(j = n; j >= 1; --j){
					if(x[jx] != 0.){
						if(nounit){
							x[jx] /= a[j+j*lda];
						}
						T temp = x[jx];
						ix = jx;
						for(i = j - 1; i >= 1; --i){
							ix -= incx;
							x[ix] -= temp * a[i+j*lda];
						}
					}
					jx -= incx;
				}
			}else{
				jx = kx;
				for(j = 1; j <= n; ++j){
					if(x[jx] != 0.){
						if(nounit){
							x[jx] /= a[j+j*lda];
						}
						T temp = x[jx];
						ix = jx;
						for(i = j + 1; i <= n; ++i){
							ix += incx;
							x[ix] -= temp * a[i+j*lda];
						}
					}
					jx += incx;
				}
			}
		}else{ // Form  x := inv( A' )*x  or  x := inv( conjg( A' ) )*x.
			if(uplo == 'U'){
				jx = kx;
				for(j = 1; j <= n; ++j){
					ix = kx;
					T temp = x[jx];
					if(noconj){
						for(i = 1; i <= j - 1; ++i){
							temp -= a[i+j*lda] * x[ix];
							ix += incx;
						}
						if(nounit){
							temp /= a[j+j*lda];
						}
					}else{
						for(i = 1; i <= j - 1; ++i){
							temp -= _RealOrComplexChooser<T>::_conj(a[i+j*lda]) * x[ix];
							ix += incx;
						}
						if(nounit){
							temp /= _RealOrComplexChooser<T>::_conj(a[j+j*lda]);
						}
					}
					x[jx] = temp;
					jx += incx;
				}
			}else{
				kx += (n - 1) * incx;
				jx = kx;
				for(j = n; j >= 1; --j){
					ix = kx;
					T temp = x[jx];
					if(noconj){
						for(i = n; i >= j + 1; --i){
							temp -= a[i+j*lda] * x[ix];
							ix -= incx;
						}
						if(nounit){
							temp /= a[j+j*lda];
						}
					}else{
						for(i = n; i >= j + 1; --i){
							temp -= _RealOrComplexChooser<T>::_conj(a[i+j*lda]) * x[ix];
							ix -= incx;
						}
						if(nounit){
							temp /= _RealOrComplexChooser<T>::_conj(a[j+j*lda]);
						}
					}
					x[jx] = temp;
					jx -= incx;
				}
			}
		}
		*/
		
		if(n == 0){ return; }

		const bool noconj = (trans == 'T');
		const bool nounit = (diag == 'N');

		// Start the operations. In this version the elements of A are
		// accessed sequentially with one pass through A.

		if(trans == 'N'){ // Form  x := inv( A )*x.
			if(uplo == 'U'){
				T *xjx = x + (n-1)*incx;
				size_t j = n; while(j --> 0){
					if(*xjx != 0.){
						if(nounit){
							*xjx /= a[j+j*lda];
						}
						T temp = *xjx;
						T *xix = xjx;
						size_t i = j; while(i --> 0){
							xix -= incx;
							*xix -= temp * a[i+j*lda];
						}
					}
					xjx -= incx;
				}
			}else{
				T *xjx = x;
				for(size_t j = 0; j < n; ++j){
					if(*xjx != 0.){
						if(nounit){
							*xjx /= a[j+j*lda];
						}
						T temp = *xjx;
						T *xix = xjx;
						for(size_t i = j; i < n; ++i){
							xix += incx;
							*xix -= temp * a[i+j*lda];
						}
					}
					xjx += incx;
				}
			}
		}else{ // Form  x := inv( A' )*x  or  x := inv( conjg( A' ) )*x.
			if(uplo == 'U'){
				T *xjx = x;
				for(size_t j = 0; j < n; ++j){
					T* xix = x;
					T temp = *xjx;
					if(noconj){
						for(size_t i = 0; i < j; ++i){
							temp -= a[i+j*lda] * *xix;
							xix += incx;
						}
						if(nounit){
							temp /= a[j+j*lda];
						}
					}else{
						for(size_t i = 0; i < j; ++i){
							temp -= _RealOrComplexChooser<T>::_conj(a[i+j*lda]) * *xix;
							xix += incx;
						}
						if(nounit){
							temp /= _RealOrComplexChooser<T>::_conj(a[j+j*lda]);
						}
					}
					*xjx = temp;
					xjx += incx;
				}
			}else{
				T *xjx = x + (n-1)*incx;
				size_t j = n; while(j --> 0){
					T *xix = x + (n-1)*incx;
					T temp = *xjx;
					if(noconj){
						size_t i = n; while(i --> j+2){
							temp -= a[i+j*lda] * *xix;
							xix -= incx;
						}
						if(nounit){
							temp /= a[j+j*lda];
						}
					}else{
						size_t i = n; while(i --> j+2){
							temp -= _RealOrComplexChooser<T>::_conj(a[i+j*lda]) * *xix;
							xix -= incx;
						}
						if(nounit){
							temp /= _RealOrComplexChooser<T>::_conj(a[j+j*lda]);
						}
					}
					*xjx = temp;
					xjx -= incx;
				}
			}
		}
	}
};


template <class A, class T> // zgeru, dgeru, cgeru, sgeru
void Rank1Update(size_t m, size_t n, const A &alpha, const T *x, size_t incx, const T *y, size_t incy, T *a, size_t lda){
	// performs A += alpha*x*conj(y')
	if(m < 1 || n < 1 || A(0) == alpha){ return; }
	for(size_t j = 0; j < n; ++j){
		if(T(0) != *y){
			T temp = alpha*(*y);
			const T *xx = x;
			for(size_t i = 0; i < m; ++i){
				a[i+j*lda] += (*xx)*temp;
				xx += incx;
			}
		}
		y += incy;
	}
}

template <class A, class T> // zgerc, cgerc
void ConjugateRank1Update(size_t m, size_t n, const A &alpha, const T *x, size_t incx, const T *y, size_t incy, T *a, size_t lda){
	// performs A += alpha*x*conj(y')
	if(m < 1 || n < 1 || A(0) == alpha){ return; }
	for(size_t j = 0; j < n; ++j){
		if(T(0) != *y){
			T temp = alpha*_RealOrComplexChooser<T>::_conj(*y);
			const T *xx = x;
			for(size_t i = 0; i < m; ++i){
				a[i+j*lda] += (*xx)*temp;
				xx += incx;
			}
		}
		y += incy;
	}
}

//// Level 3

// MultMM (zgemm)
template <char transa='N', char transb='N'>
struct MultMM{
	template <class A, class B, class T>
	MultMM(size_t m, size_t n, size_t k, const A &alpha, const T *a, size_t lda,
		const T *b, size_t ldb, const B &beta, T *c, size_t ldc){

		// ZGEMM  performs one of the matrix-matrix operations
		//    C := alpha*op( A )*op( B ) + beta*C,
		// where  op( X ) is one of
		//    op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
		// alpha and beta are scalars, and A, B and C are matrices, with op( A )
		// an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

		// Arguments
		// ==========

		// TRANSA - On entry, TRANSA specifies the form of op( A ) to be used in
		//          the matrix multiplication as follows:
		//             TRANSA = 'N' or 'n',  op( A ) = A.
		//             TRANSA = 'T' or 't',  op( A ) = A'.
		//             TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).

		// TRANSB - On entry, TRANSB specifies the form of op( B ) to be used in
		//          the matrix multiplication as follows:
		//             TRANSB = 'N' or 'n',  op( B ) = B.
		//             TRANSB = 'T' or 't',  op( B ) = B'.
		//             TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).

		// M      - On entry,  M  specifies  the number  of rows  of the  matrix
		//          op( A )  and of the  matrix  C.

		// N      - On entry,  N  specifies the number  of columns of the matrix
		//          op( B ) and the number of columns of the matrix C.

		// K      - On entry,  K  specifies  the number of columns of the matrix
		//          op( A ) and the number of rows of the matrix op( B ).

		// ALPHA  - On entry, ALPHA specifies the scalar alpha.

		// A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
		//          k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
		//          Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
		//          part of the array  A  must contain the matrix  A,  otherwise
		//          the leading  k by m  part of the array  A  must contain  the
		//          matrix A.

		// LDA    - On entry, LDA specifies the first dimension of A as declared
		//          in the calling (sub) program. When  TRANSA = 'N' or 'n' then
		//          LDA must be at least  max( 1, m ), otherwise  LDA must be at
		//          least  max( 1, k ).

		// B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
		//          n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
		//          Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
		//          part of the array  B  must contain the matrix  B,  otherwise
		//          the leading  n by k  part of the array  B  must contain  the
		//          matrix B.
		//          Unchanged on exit.

		// LDB    - On entry, LDB specifies the first dimension of B as declared
		//          in the calling (sub) program. When  TRANSB = 'N' or 'n' then
		//          LDB must be at least  max( 1, k ), otherwise  LDB must be at
		//          least  max( 1, n ).

		// BETA   - On entry,  BETA  specifies the scalar  beta.  When  BETA  is
		//          supplied as zero then C need not be set on input.

		// C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
		//          Before entry, the leading  m by n  part of the array  C must
		//          contain the matrix  C,  except when  beta  is zero, in which
		//          case C need not be set on entry.
		//          On exit, the array  C  is overwritten by the  m by n  matrix
		//          ( alpha*op( A )*op( B ) + beta*C ).

		// LDC    - On entry, LDC specifies the first dimension of C as declared
		//          in  the  calling  (sub)  program.   LDC  must  be  at  least
		//          max( 1, m ).

		//    Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
		//    conjugated or transposed, set  CONJA and CONJB  as true if  A  and
		//    B  respectively are to be  transposed but  not conjugated  and set
		//    NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
		//    and the number of rows of  B  respectively.

		const bool nota = (transa == 'N');
		const bool notb = (transb == 'N');
		const bool conja = (transa == 'C');
		const bool conjb = (transb == 'C');

		if(m == 0 || n == 0 || ((A(0) == alpha || k == 0) && (B(1) == beta))){ return; }

		if(A(0) == alpha){
			if(B(0) == beta){
				for(size_t j = 0; j < n; ++j){
					for(size_t i = 0; i < m; ++i){
						c[i+j*ldc] = 0;
					}
				}
			}else{
				for(size_t j = 0; j < n; ++j){
					for(size_t i = 0; i < m; ++i){
						c[i+j*ldc] *= beta;
					}
				}
			}
			return;
		}

		if(notb){
			if(nota){ // Form  C := alpha*A*B + beta*C.
				for(size_t j = 0; j < n; ++j){
					if(B(0) == beta){
						for(size_t i = 0; i < m; ++i){
							c[i+j*ldc] = 0;
						}
					}else if(B(1) != beta){
						for(size_t i = 0; i < m; ++i){
							c[i+j*ldc] *= beta;
						}
					}
					for(size_t l = 0; l < k; ++l){
						if(T(0) != b[l+j*ldb]){
							T temp = alpha * b[l+j*ldb];
							for(size_t i = 0; i < m; ++i){
								c[i+j*ldc] += temp * a[i+l*lda];
							}
						}
					}
				}
			}else if(conja){ // Form  C := alpha*conjg( A' )*B + beta*C.
				for(size_t j = 0; j < n; ++j){
					for(size_t i = 0; i < m; ++i){
						T temp = 0;
						for(size_t l = 0; l < k; ++l){
							temp += _RealOrComplexChooser<T>::_conj(a[l+i*lda]) * b[l+j*ldb];
						}
						if(B(0) == beta){
							c[i+j*ldc] = alpha * temp;
						}else{
							c[i+j*ldc] = alpha * temp + beta * c[i+j*ldc];
						}
					}
				}
			}else{ // Form  C := alpha*A'*B + beta*C
				for(size_t j = 0; j < n; ++j){
					for(size_t i = 0; i < m; ++i){
						T temp = 0;
						for(size_t l = 0; l < k; ++l){
							temp += a[l+i*lda] * b[l+j*ldb];
						}
						if(B(0) == beta){
							c[i+j*ldc] = alpha * temp;
						}else{
							c[i+j*ldc] = alpha * temp + beta * c[i+j*ldc];
						}
					}
				}
			}
		}else if(nota){
			if(conjb){ // Form  C := alpha*A*conjg( B' ) + beta*C.
				for(size_t j = 0; j < n; ++j){
					if(B(0) == beta){
						for(size_t i = 0; i < m; ++i){
							c[i+j*ldc] = 0;
						}
					}else if(B(1) != beta){
						for(size_t i = 0; i < m; ++i){
							c[i+j*ldc] *= beta;
						}
					}
					for(size_t l = 0; l < k; ++l){
						if(T(0) != b[j+l*ldb]){
							T temp = alpha * _RealOrComplexChooser<T>::_conj(b[j+l*ldb]);
							for(size_t i = 0; i < m; ++i){
								c[i+j*ldc] += temp * a[i+l*lda];
							}
						}
					}
				}
			}else{ // Form  C := alpha*A*B' + beta*C
				for(size_t j = 0; j < n; ++j){
					if(B(0) == beta){
						for(size_t i = 0; i < m; ++i){
							c[i+j*ldc] = 0;
						}
					}else if(B(1) != beta){
						for(size_t i = 0; i < m; ++i){
							c[i+j*ldc] *= beta;
						}
					}
					for(size_t l = 0; l < k; ++l){
						if(T(0) != b[j+l*ldb]){
							T temp = alpha * b[j+l*ldb];
							for(size_t i = 0; i < m; ++i){
								c[i+j*ldc] += temp * a[i+l*lda];
							}
						}
					}
				}
			}
		}else if(conja){
			if(conjb){ // Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
				for(size_t j = 0; j < n; ++j){
					for(size_t i = 0; i < m; ++i){
						T temp = 0;
						for(size_t l = 0; l < k; ++l){
							temp += _RealOrComplexChooser<T>::_conj(a[l+i*lda]) * _RealOrComplexChooser<T>::_conj(b[j+l*ldb]);
						}
						if(B(0) == beta){
							c[i+j*ldc] = alpha * temp;
						}else{
							c[i+j*ldc] = alpha * temp + beta * c[i+j*ldc];
						}
					}
				}
			}else{ // Form  C := alpha*conjg( A' )*B' + beta*C
				for(size_t j = 0; j < n; ++j){
					for(size_t i = 0; i < m; ++i){
						T temp = 0;
						for(size_t l = 0; l < k; ++l){
							temp += _RealOrComplexChooser<T>::_conj(a[l+i*lda]) * b[j+l*ldb];
						}
						if(B(0) == beta){
							c[i+j*ldc] = alpha * temp;
						}else{
							c[i+j*ldc] = alpha * temp + beta * c[i+j*ldc];
						}
					}
				}
			}
		}else{
			if(conjb){ // Form  C := alpha*A'*conjg( B' ) + beta*C
				for(size_t j = 0; j < n; ++j){
					for(size_t i = 0; i < m; ++i){
						T temp = 0;
						for(size_t l = 0; l < k; ++l){
							temp += a[l + i * lda] * _RealOrComplexChooser<T>::_conj(b[j+l*ldb]);
						}
						if(B(0) == beta){
							c[i+j*ldc] = alpha * temp;
						}else{
							c[i+j*ldc] = alpha * temp + beta * c[i+j*ldc];
						}
					}
				}
			}else{ // Form  C := alpha*A'*B' + beta*C
				for(size_t j = 0; j < n; ++j){
					for(size_t i = 0; i < m; ++i){
						T temp = 0;
						for(size_t l = 0; l < k; ++l){
							temp += a[l+i*lda] * b[j+l*ldb];
						}
						if(B(0) == beta){
							c[i+j*ldc] = alpha * temp;
						}else{
							c[i+j*ldc] = alpha * temp + beta * c[i+j*ldc];
						}
					}
				}
			}
		}
	}
};

template <char side, char uplo, char transa, char diag>
struct MultTrM{ // ztrmm, dtrmm, ctrmm, strmm
	template <class T, class TA>
	MultTrM(size_t m, size_t n, const TA &alpha, const T *a, size_t lda, T *b, size_t ldb){
		// Purpose
		// =======
		// 
		// ZTRMM  performs one of the matrix-matrix operations
		//   B := alpha*op( A )*B,   or   B := alpha*B*op( A )
		// where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
		// non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
		//   op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
		// 
		// Arguments
		// ==========
		// 
		// SIDE   - char*1.
		// On entry,  SIDE specifies whether  op( A ) multiplies B from
		// the left or right as follows:
		//   SIDE = 'L' or 'l'   B := alpha*op( A )*B.
		//   SIDE = 'R' or 'r'   B := alpha*B*op( A ).
		// 
		// UPLO   - char*1.
		// On entry, UPLO specifies whether the matrix A is an upper or
		// lower triangular matrix as follows:
		//   UPLO = 'U' or 'u'   A is an upper triangular matrix.
		//   UPLO = 'L' or 'l'   A is a lower triangular matrix.
		// 
		// TRANSA - char*1.
		// On entry, TRANSA specifies the form of op( A ) to be used in
		// the matrix multiplication as follows:
		//   TRANSA = 'N' or 'n'   op( A ) = A.
		//   TRANSA = 'T' or 't'   op( A ) = A'.
		//   TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
		// 
		// DIAG   - char*1.
		// On entry, DIAG specifies whether or not A is unit triangular
		// as follows:
		//   DIAG = 'U' or 'u'   A is assumed to be unit triangular.
		//   DIAG = 'N' or 'n'   A is not assumed to be unit triangular.
		// 
		// M      - int.
		// On entry, M specifies the number of rows of B. M must be at
		// least zero.
		// 
		// N      - int.
		// On entry, N specifies the number of columns of B.  N must be
		// at least zero.
		// 
		// ALPHA  - std::complex<double>      .
		// On entry,  ALPHA specifies the scalar  alpha. When  alpha is
		// zero then  A is not referenced and  B need not be set before
		// entry.
		// 
		// A      - std::complex<double>       array of DIMENSION ( LDA, k ), where k is m
		// when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
		// Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
		// upper triangular part of the array  A must contain the upper
		// triangular matrix  and the strictly lower triangular part of
		// A is not referenced.
		// Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
		// lower triangular part of the array  A must contain the lower
		// triangular matrix  and the strictly upper triangular part of
		// A is not referenced.
		// Note that when  DIAG = 'U' or 'u',  the diagonal elements of
		// A  are not referenced either,  but are assumed to be  unity.
		// Unchanged on exit.
		// 
		// LDA    - int.
		// On entry, LDA specifies the first dimension of A as declared
		// in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
		// LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
		// then LDA must be at least max( 1, n ).
		// 
		// B      - std::complex<double>       array of DIMENSION ( LDB, n ).
		// Before entry,  the leading  m by n part of the array  B must
		// contain the matrix  B,  and  on exit  is overwritten  by the
		// transformed matrix.
		// 
		// LDB    - int.
		// On entry, LDB specifies the first dimension of B as declared
		// in  the  calling  (sub)  program.   LDB  must  be  at  least
		// max( 1, m ).
		// 

		const bool lside = (side == 'L');
		const bool noconj = (transa == 'T');
		const bool nounit = (diag == 'N');
		const bool upper = (uplo == 'U');
	 
		if(m == 0 || n == 0){ return; }

		if(alpha == TA(0)){
			for(size_t j = 0; j < n; ++j){
				for(size_t i = 0; i < m; ++i){
					b[i+j*ldb] = T(0);
				}
				}
			return;
		}

		if(lside){
			if(transa == 'N'){ // Form  B := alpha*A*B.
				if(upper){
					for(size_t j = 0; j < n; ++j){
						for(size_t k = 0; k < m; ++k){
							if(b[k+j*ldb] != T(0)){
								T temp = alpha*b[k+j*ldb];
								for(size_t i = 0; i < k; ++i){
									b[i+j*ldb] += temp*a[i+k*lda];
								}
								if(nounit) temp *= a[k+k*lda];
								b[k+j*ldb] = temp;
							}
						}
					}
				}else{
					for(size_t j = 0; j < n; ++j){
						for(size_t k = m-1; (int)k >= 0; --k){
							if(b[k+j*ldb] != T(0)){
								T temp = alpha*b[k+j*ldb];
								b[k+j*ldb] = temp;
								if(nounit) b[k+j*ldb] *= a[k+k*lda];
								for(size_t i = k+1; i < m; ++i){
									b[i+j*ldb] += temp*a[i+k*lda];
								}
							}
						}
					}
				}
			}else{ // Form  B := alpha*A'*B   or   B := alpha*conjg( A' )*B.
				if(upper){
					for(size_t j = 0; j < n; ++j){
						for(size_t i = m-1; (int)i >= 0; --i){
							T temp = b[i+j*ldb];
							if(noconj){
								if(nounit) temp *= a[i+i*lda];
								for(size_t k = 0; k < i; ++k){
									temp += a[k+i*lda]*b[k+j*ldb];
								}
							}else{
								if(nounit) temp *= std::conj(a[i+i*lda]);
								for(size_t k = 0; k < i; ++k){
									temp += std::conj(a[k+i*lda])*b[k+j*ldb];
								}
							}
							b[i+j*ldb] = alpha*temp;
						}
					}
				}else{
					for(size_t j = 0; j < n; ++j){
						for(size_t i = 0; i < m; ++i){
							T temp = b[i+j*ldb];
							if(noconj){
								if(nounit) temp *= a[i+i*lda];
								for(size_t k = i+1; k < m; ++k){
									temp += a[k+i*lda]*b[k+j*ldb];
								}
							}else{
								if(nounit) temp *= std::conj(a[i+i*lda]);
								for(size_t k = i+1; k < m; ++k){
									temp += std::conj(a[k+i*lda])*b[k+j*ldb];
								}
							}
							b[i+j*ldb] = alpha*temp;
						}
					}
				}
			}
		}else{
			if(transa == 'N'){ // Form  B := alpha*B*A.
				if(upper){
					for(size_t j = n-1; (int)j >= 0; --j){
						T temp = alpha;
						if(nounit) temp *= a[j+j*lda];
						for(size_t i = 0; i < m; ++i){
							b[i+j*ldb] *= temp;
						}
						for(size_t k = 0; k < j; ++k){
							if(a[k+j*lda] != T(0)){
								temp = alpha*a[k+j*lda];
								for(size_t i = 0; i < m; ++i){
									b[i+j*ldb] += temp*b[i+k*ldb];
								}
							}
						}
					}
				}else{
					for(size_t j = 0; j < n; ++j){
						T temp = alpha;
						if(nounit) temp *= a[j+j*lda];
						for(size_t i = 0; i < m; ++i){
							b[i+j*ldb] *= temp;
						}
						for(size_t k = j+1; k < n; ++k){
							if(a[k+j*lda] != T(0)){
								temp = alpha*a[k+j*lda];
								for(size_t i = 0; i < m; ++i){
									b[i+j*ldb] += temp*b[i+k*ldb];
								}
							}
						}
					}
				}
			}else{ // Form  B := alpha*B*A'   or   B := alpha*B*conjg( A' ).
				if(upper){
					for(size_t k = 0; k < n; ++k){
						T temp;
						for(size_t j = 0; j < k; ++j){
							if(a[j+k*lda] != T(0)){
								if(noconj){
									temp = alpha*a[j+k*lda];
								}else{
									temp = alpha*std::conj(a[j+k*lda]);
								}
								for(size_t i = 0; i < m; ++i){
									b[i+j*ldb] += temp*b[i+k*ldb];
								}
							}
						}
						temp = alpha;
						if(nounit){
							if(noconj){
								temp *= a[k+k*lda];
							}else{
								temp *= std::conj(a[k+k*lda]);
							}
						}
						if(temp != T(1)){
							for(size_t i = 0; i < m; ++i){
								b[i+k*ldb] *= temp;
							}
						}
					}
				}else{
					for(size_t k = n-1; (int)k >= 0; --k){
						T temp;
						for(size_t j = k+1; j < n; ++j){
							if(a[j+k*lda] != T(0)){
								if(noconj){
									temp = alpha*a[j+k*lda];
								}else{
									temp = alpha*std::conj(a[j+k*lda]);
								}
								for(size_t i = 0; i < m; ++i){
									b[i+j*ldb] += temp*b[i+k*ldb];
								}
							}
						}
						temp = alpha;
						if(nounit){
							if(noconj){
								temp *= a[k+k*lda];
							}else{
								temp = temp*std::conj(a[k+k*lda]);
							}
						}
						if(temp != T(1)){
							for(size_t i = 0; i < m; ++i){
								b[i+k*ldb] *= temp;
							}
						}
					}
				}
			}
		}
	}
};


template <char side='L', char uplo='L', char transa='N', char diag='N'>
struct SolveTrM{ // ztrsm, dtrsm, ctrsm, strsm
	template <class TA, class T>
	SolveTrM(size_t m, size_t n, const TA &alpha, const T *a, size_t lda, T *b, size_t ldb){
		// ZTRSM  solves one of the matrix equations
		//    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
		// where alpha is a scalar, X and B are m by n matrices, A is a unit, or
		// non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
		//    op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
		// The matrix X is overwritten on B.

		// Arguments
		// ==========

		// SIDE   - whether op( A ) appears on the left or right of X as follows:
		//             SIDE = 'L' or 'l'   op( A )*X = alpha*B.
		//             SIDE = 'R' or 'r'   X*op( A ) = alpha*B.

		// UPLO   - whether the matrix A is an upper or lower triangular matrix as follows:
		//             UPLO = 'U' or 'u'   A is an upper triangular matrix.
		//             UPLO = 'L' or 'l'   A is a lower triangular matrix.

		// TRANSA - the form of op( A ) to be used in the matrix multiplication as follows:
		//             TRANSA = 'N' or 'n'   op( A ) = A.
		//             TRANSA = 'T' or 't'   op( A ) = A'.
		//             TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).

		// DIAG   - whether or not A is unit triangular as follows:
		//             DIAG = 'U' or 'u'   A is assumed to be unit triangular.
		//             DIAG = 'N' or 'n'   A is not assumed to be unit triangular.

		// M      - number of rows of B.
		// N      - number of columns of B.

		// ALPHA  - the scalar  alpha. When  alpha is
		//          zero then  A is not referenced and  B need not be set before
		//          entry.

		// A      - array of DIMENSION ( LDA, k ), where k is m
		//          when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
		//          Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
		//          upper triangular part of the array  A must contain the upper
		//          triangular matrix  and the strictly lower triangular part of
		//          A is not referenced.
		//          Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
		//          lower triangular part of the array  A must contain the lower
		//          triangular matrix  and the strictly upper triangular part of
		//          A is not referenced.
		//          Note that when  DIAG = 'U' or 'u',  the diagonal elements of
		//          A  are not referenced either,  but are assumed to be  unity.

		// LDA    - the first dimension of A as declared
		//          in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
		//          LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
		//          then LDA must be at least max( 1, n ).

		// B      - array of DIMENSION ( LDB, n ).
		//          Before entry,  the leading  m by n part of the array  B must
		//          contain  the  right-hand  side  matrix  B,  and  on exit  is
		//          overwritten by the solution matrix  X.

		// LDB    - the first dimension of B as declared
		//          in  the  calling  (sub)  program.   LDB  must  be  at  least
		//          max( 1, m ).
				
		const bool lside = (side == 'L');
		const bool noconj = (transa == 'T');
		const bool nounit = (diag == 'N');
		const bool upper = (uplo == 'U');

		if(m == 0 || n == 0){ return; }

		if(TA(0) == alpha){
			for(size_t j = 0; j < n; ++j){
				for(size_t i = 0; i < m; ++i){
					b[i+j*ldb] = 0;
				}
			}
			return;
		}

		if(lside){
			if(transa == 'N'){ // Form  B := alpha*inv( A )*B.
				if(upper){
					for(size_t j = 0; j < n; ++j){
						if(TA(1) != alpha){
							for(size_t i = 0; i < m; ++i){
								b[i+j*ldb] = alpha * b[i+j*ldb];
							}
						}
						for(int k = (int)m-1; k >= 0; --k){
							if(T(0) != b[k+j*ldb]){
								if(nounit){
									b[k+j*ldb] /= a[k+k*lda];
								}
								for(size_t i = 0; i < (size_t)k; ++i){
									b[i+j*ldb] -= b[k+j*ldb] * a[i+k*lda];
								}
							}
						}
					}
				}else{
					for(size_t j = 0; j < n; ++j){
						if(TA(1) != alpha){
							for(size_t i = 0; i < m; ++i){
								b[i+j*ldb] = alpha * b[i+j*ldb];
							}
						}
						for(size_t k = 0; k < m; ++k){
							if(T(0) != b[k+j*ldb]){
								if(nounit){
									b[k+j*ldb] /= a[k+k*lda];
								}
								for(size_t i = k+1; i < m; ++i){
									b[i+j*ldb] -= b[k+j*ldb] * a[i+k*lda];
								}
							}
						}
					}
				}
			}else{ // Form  B := alpha*inv( A' )*B or B := alpha*inv( conjg( A' ) )*B
				if(upper){
					for(size_t j = 0; j < n; ++j){
						for(size_t i = 0; i < m; ++i){
							T temp = alpha * b[i+j*ldb];
							if(noconj){
								for(size_t k = 0; k < i; ++k){
									temp -= a[k+i*lda] * b[k+j*ldb];
								}
								if(nounit){
									temp /= a[i+i*lda];
								}
							}else{
								for(size_t k = 0; k < i; ++k){
									temp -= _RealOrComplexChooser<T>::_conj(a[k+i*lda]) * b[k+j*ldb];
								}
								if(nounit){
									temp /= _RealOrComplexChooser<T>::_conj(a[i+i*lda]);
								}
							}
							b[i+j*ldb] = temp;
						}
					}
				}else{
					for(size_t j = 0; j < n; ++j){
						for(int i = (int)m-1; i >= 0; --i){
							T temp = alpha * b[i+j*ldb];
							if(noconj){
								for(size_t k = i+1; k < m; ++k){
									temp -= a[k+i*lda] * b[k+j*ldb];
								}
								if(nounit){
									temp /= a[i+i*lda];
								}
							}else{
								for(size_t k = i+1; k < m; ++k){
									temp -= _RealOrComplexChooser<T>::_conj(a[k+i*lda]) * b[k+j*ldb];
								}
								if(nounit){
									temp /= _RealOrComplexChooser<T>::_conj(a[i+i*lda]);
								}
							}
							b[i+j*ldb] = temp;
						}
					}
				}
			}
		}else{
			if(transa == 'N'){ // Form  B := alpha*B*inv( A )
				if(upper){
					for(size_t j = 0; j < n; ++j){
						if(TA(1) != alpha){
							for(size_t i = 0; i < m; ++i){
								b[i+j*ldb] = alpha * b[i+j*ldb];
							}
						}
						for(size_t k = 0; k < j; ++k){
							if(T(0) != a[k+j*lda]){
								for(size_t i = 0; i < m; ++i){
									b[i+j*ldb] -= a[k+j*lda] * b[i+k*ldb];
								}
							}
						}
						if(nounit){
							T temp = T(1) / a[j+j*lda];
							for(size_t i = 0; i < m; ++i){
								b[i+j*ldb] = temp * b[i+j*ldb];
							}
						}
					}
				}else{
					for(int j = (int)n-1; j >= 0; --j){
						if(TA(1) != alpha){
							for(size_t i = 0; i < m; ++i){
								b[i+j*ldb] = alpha * b[i+j*ldb];
							}
						}
						for(size_t k = j+1; k < n; ++k){
							if(T(0) != a[k+j*lda]){
								for(size_t i = 0; i < m; ++i){
									b[i+j*ldb] -= a[k+j*lda] * b[i+k*ldb];
								}
							}
						}
						if(nounit){
							T temp = T(1) / a[j+j*lda];
							for(size_t i = 0; i < m; ++i){
								b[i+j*ldb] = temp * b[i+j*ldb];
							}
						}
					}
				}
			}else{ // Form  B := alpha*B*inv( A' ) or B := alpha*B*inv( conjg( A' ) )
				if(upper){
					for(int k = (int)n-1; k >= 0; --k){
						T temp;
						if(nounit){
							if(noconj){
								temp = T(1) / a[k+k*lda];
							}else{
								temp = T(1) / _RealOrComplexChooser<T>::_conj(a[k+k*lda]);
							}
							for(size_t i = 0; i < m; ++i){
								b[i+k*ldb] = temp * b[i+k*ldb];
							}
						}
						for(size_t j = 0; j < (size_t)k; ++j){
							if(T(0) != a[j+k*lda]){
								if(noconj){
									temp = a[j+k*lda];
								}else{
									temp = _RealOrComplexChooser<T>::_conj(a[j+k*lda]);
								}
								for(size_t i = 0; i < m; ++i){
									b[i+j*ldb] -= temp * b[i+k*ldb];
								}
							}
						}
						if(T(1) != alpha){
							for(size_t i = 0; i < m; ++i){
								b[i+k*ldb] = alpha * b[i+k*ldb];
							}
						}
					}
				}else{
					for(size_t k = 0; k < n; ++k){
						T temp;
						if(nounit){
							if(noconj){
								temp = T(1) / a[k+k*lda];
							}else{
								temp = T(1) / _RealOrComplexChooser<T>::_conj(a[k+k*lda]);
							}
							for(size_t i = 0; i < m; ++i){
								b[i+k*ldb] = temp * b[i+k*ldb];
							}
						}
						for(size_t j = k+1; j < n; ++j){
							if(T(0) != a[j+k*lda]){
								if(noconj){
									temp = a[j+k*lda];
								}else{
									temp = _RealOrComplexChooser<T>::_conj(a[j+k*lda]);
								}
								for(size_t i = 0; i < m; ++i){
									b[i+j*ldb] -= temp * b[i+k*ldb];
								}
							}
						}
						if(T(1) != alpha){
							for(size_t i = 0; i < m; ++i){
								b[i+k*ldb] = alpha * b[i+k*ldb];
							}
						}
					}
				}
			}
		}
	}
};


//// Common/simple LAPACK routines

template <class T>
void Conjugate(size_t n, T *x, size_t incx){ // zlacgv, clacgv
	while(n --> 0){
		_RealOrComplexChooser<T>::_conj(*x);
		x += incx;
	}
}

template <char uplo='A'>
struct CopyMatrix{ // zlacpy, dlacpy, clacpy, slacpy
	template <class T>
	CopyMatrix(size_t m, size_t n, const T *a, size_t lda, T *b, size_t ldb){
		if('U' == uplo){
			for(size_t j = 0; j < n; ++j){
				size_t ilimit = j+1; if(m < ilimit){ ilimit = m; }
				for(size_t i = 0; i < ilimit; ++i){
					b[i+j*ldb] = a[i+j*lda];
				}
			}
		}else if('L' == uplo){
			for(size_t j = 0; j < n; ++j){
				for(size_t i = j; i < m; ++i){
					b[i+j*ldb] = a[i+j*lda];
				}
			}
		}else{
			for(size_t j = 0; j < n; ++j){
				for(size_t i = 0; i < m; ++i){
					b[i+j*ldb] = a[i+j*lda];
				}
			}
		}
	}
};

template <char uplo='F'>
struct SetMatrix{ // zlaset, dlaset, claset, slaset
	template <class TA, class TB, class T>
	SetMatrix(size_t m, size_t n, const TA &offdiag, const TB &diag, T *a, size_t lda){
		const size_t minmn = ((m < n) ? m : n);
		if('U' == uplo){
			for(size_t j = 1; j < n; ++j){
				size_t ilimit = m; if(j < ilimit){ ilimit = j; }
				for(size_t i = 0; i < ilimit; ++i){
					a[i+j*lda] = offdiag;
				}
			}
		}else if('L' == uplo){
			for(size_t j = 0; j < minmn; ++j){
				for(size_t i = j+1; i < m; ++i){
					a[i+j*lda] = offdiag;
				}
			}
		}else{
			for(size_t j = 0; j < n; ++j){
				for(size_t i = 0; i < m; ++i){
					a[i+j*lda] = offdiag;
				}
			}
		}
		for(size_t i = 0; i < minmn; ++i){
			a[i+i*lda] = diag;
		}
	}
};

#if defined(RNP_TBLAS_USE_RANDOM)
#include "Random.h"

template <class T>
struct _RealOrComplexRand{
	typedef T value_type;
	typedef T real_type;

	inline static void randvec(size_t n, value_type *x, int iseed[4] = NULL){
		RNP::Random::RandomRealsUniform01(n, x, iseed);
	}
	inline static size_t numreals(){ return 1; }
};
template <class T>
struct _RealOrComplexRand<std::complex<T> >{
	typedef std::complex<T> value_type;
	typedef typename std::complex<T>::value_type real_type;

	inline static void randvec(size_t n, value_type *x, int iseed[4] = NULL){
		RNP::Random::RandomRealsUniform01((sizeof(value_type)/sizeof(real_type))*n, reinterpret_cast<real_type*>(x), iseed);
	}
	inline static size_t numreals(){ return 2; }
};

// dist specifies the distribution of the random numbers:
//      = 1:  real and imaginary parts each uniform (0,1)
//      = 2:  real and imaginary parts each uniform (-1,1)
//      = 3:  real and imaginary parts each normal (0,1)
//      = 4:  uniformly distributed on the disc abs(z) < 1
//      = 5:  uniformly distributed on the circle abs(z) = 1
template <int dist=2>
struct RandomVector{
	template <class T>
	RandomVector(size_t n, T *x, int iseed[4] = NULL){
		using namespace std;
		static const size_t chunksize = 128;
		for(size_t i = 0; i < n; i += chunksize){
			size_t chunk = n-i; if(chunksize < chunk){ chunk = chunksize; }
			_RealOrComplexRand<T>::randvec(chunk, &x[i], iseed);
			if(1 == dist){
			}else if(2 == dist){
				for(size_t j = 0; j < chunk; ++j){
					x[i+j] = 2*x[i+j]-1;
				}
			}
		}
	}
	template <class T>
	RandomVector(size_t n, std::complex<T> *x, int iseed[4] = NULL){
		using namespace std;
		typedef std::complex<T> value_type;
		static const size_t chunksize = 64;
		typedef typename _RealOrComplexRand<value_type>::real_type real_type;
		real_type buf[2*chunksize];
		for(size_t i = 0; i < n; i += chunksize){
			size_t chunk = n-i; if(chunksize < chunk){ chunk = chunksize; }
			_RealOrComplexRand<real_type>::randvec(2*chunk, buf, iseed);
			if(1 == dist){
				for(size_t j = 0; j < chunk; ++j){
					x[i+j] = value_type(buf[2*j+0], buf[2*j+1]);
				}
			}else if(2 == dist){
				for(size_t j = 0; j < chunk; ++j){
					x[i+j] = value_type(2*buf[2*j+0]-1, 2*buf[2*j+1]-1);
				}
			}else if(3 == dist){
				for(size_t j = 0; j < chunk; ++j){
					x[i+j] = sqrt(-2*log(buf[2*j+0])) * exp(value_type(0,2*M_PI*buf[2*j+1]));
				}
			}else if(4 == dist){
				for(size_t j = 0; j < chunk; ++j){
					x[i+j] = sqrt(buf[2*j+0]) * exp(value_type(0,2*M_PI*buf[2*j+1]));
				}
			}else if(5 == dist){
				for(size_t j = 0; j < chunk; ++j){
					x[i+j] = exp(value_type(0,2*M_PI*buf[2*j+0]));
				}
			}
		}
	}
};
#endif // RNP_TBLAS_USE_RANDOM

}; // namespace TBLAS
}; // namespace RNP

#ifdef RNP_HAVE_BLAS
# include "TBLAS_ext.h"
#endif

#endif // _RNP_TBLAS_H_
