#ifndef _RNP_TBLAS_EXT_H_
#define _RNP_TBLAS_EXT_H_

#define _USE_MATH_DEFINES
#include <complex>

namespace RNP{
typedef int integer;
namespace TBLAS{

//#define RNP_FORTRAN_NAME(LCASE,UCASE) LCASE ## _
#ifndef RNP_FORTRAN_NAME
# define RNP_FORTRAN_NAME(LCASE,UCASE) F77_FUNC(LCASE,UCASE)
#endif

//// Level 1
extern "C" void RNP_FORTRAN_NAME(sswap,SSWAP)(const integer &n, float *sx, const integer &incx, float *sy, const integer &incy);
extern "C" void RNP_FORTRAN_NAME(dswap,DSWAP)(const integer &n, double *sx, const integer &incx, double *sy, const integer &incy);
extern "C" void RNP_FORTRAN_NAME(cswap,CSWAP)(const integer &n, std::complex<float> *sx, const integer &incx, std::complex<float> *sy, const integer &incy);
extern "C" void RNP_FORTRAN_NAME(zswap,ZSWAP)(const integer &n, std::complex<double> *sx, const integer &incx, std::complex<double> *sy, const integer &incy);
template <>
inline void Swap<float>(size_t n, float *x, size_t incx, float *y, size_t incy){
	RNP_FORTRAN_NAME(sswap,SSWAP)(n, x, incx, y, incy);
}
template <>
inline void Swap<double>(size_t n, double *x, size_t incx, double *y, size_t incy){
	RNP_FORTRAN_NAME(dswap,DSWAP)(n, x, incx, y, incy);
}
template <>
inline void Swap<std::complex<float> >(size_t n, std::complex<float> *x, size_t incx, std::complex<float> *y, size_t incy){
	RNP_FORTRAN_NAME(cswap,CSWAP)(n, x, incx, y, incy);
}
template <>
inline void Swap<std::complex<double> >(size_t n, std::complex<double> *x, size_t incx, std::complex<double> *y, size_t incy){
	RNP_FORTRAN_NAME(zswap,ZSWAP)(n, x, incx, y, incy);
}

extern "C" void RNP_FORTRAN_NAME(sscal,SSCAL)(const integer &n, const float &sa, float *sx, const integer &incx);
extern "C" void RNP_FORTRAN_NAME(dscal,DSCAL)(const integer &n, const double &sa, double *sx, const integer &incx);
extern "C" void RNP_FORTRAN_NAME(cscal,CSCAL)(const integer &n, const std::complex<float> &sa, std::complex<float> *sx, const integer &incx);
extern "C" void RNP_FORTRAN_NAME(zscal,ZSCAL)(const integer &n, const std::complex<double> &sa, std::complex<double> *sx, const integer &incx);
extern "C" void RNP_FORTRAN_NAME(csscal,CSSCAL)(const integer &n, const float &sa, std::complex<float> *cx, const integer &incx);
extern "C" void RNP_FORTRAN_NAME(zdscal,ZDSCAL)(const integer &n, const double &sa, std::complex<double> *cx, const integer &incx);
template <>
inline void Scale<float,float>(size_t n, const float &scale, float *x, size_t incx){
	RNP_FORTRAN_NAME(sscal,SSCAL)(n, scale, x, incx);
}
template <>
inline void Scale<double,double>(size_t n, const double &scale, double *x, size_t incx){
	RNP_FORTRAN_NAME(dscal,DSCAL)(n, scale, x, incx);
}
template <>
inline void Scale<std::complex<float>,std::complex<float> >(size_t n, const std::complex<float> &scale, std::complex<float> *x, size_t incx){
	RNP_FORTRAN_NAME(cscal,CSCAL)(n, scale, x, incx);
}
template <>
inline void Scale<std::complex<double>,std::complex<double> >(size_t n, const std::complex<double> &scale, std::complex<double> *x, size_t incx){
	RNP_FORTRAN_NAME(zscal,ZSCAL)(n, scale, x, incx);
}
template <>
inline void Scale<float,std::complex<float> >(size_t n, const float &scale, std::complex<float> *x, size_t incx){
	RNP_FORTRAN_NAME(csscal,CSSCAL)(n, scale, x, incx);
}
template <>
inline void Scale<double,std::complex<double> >(size_t n, const double &scale, std::complex<double> *x, size_t incx){
	RNP_FORTRAN_NAME(zdscal,ZDSCAL)(n, scale, x, incx);
}

extern "C" void RNP_FORTRAN_NAME(scopy,SCOPY)(const integer &n, const float *sx, const integer &incx, float *sy, const integer &incy);
extern "C" void RNP_FORTRAN_NAME(dcopy,DCOPY)(const integer &n, const double *sx, const integer &incx, double *sy, const integer &incy);
extern "C" void RNP_FORTRAN_NAME(ccopy,CCOPY)(const integer &n, const std::complex<float> *sx, const integer &incx, std::complex<float> *sy, const integer &incy);
extern "C" void RNP_FORTRAN_NAME(zcopy,ZCOPY)(const integer &n, const std::complex<double> *sx, const integer &incx, std::complex<double> *sy, const integer &incy);
template <>
inline void Copy<float>(size_t n, const float *src, size_t incsrc, float *dst, size_t incdst){
	RNP_FORTRAN_NAME(scopy,SCOPY)(n, src, incsrc, dst, incdst);
}
template <>
inline void Copy<double>(size_t n, const double *src, size_t incsrc, double *dst, size_t incdst){
	RNP_FORTRAN_NAME(dcopy,DCOPY)(n, src, incsrc, dst, incdst);
}
template <>
inline void Copy<std::complex<float> >(size_t n, const std::complex<float> *src, size_t incsrc, std::complex<float> *dst, size_t incdst){
	RNP_FORTRAN_NAME(ccopy,CCOPY)(n, src, incsrc, dst, incdst);
}
template <>
inline void Copy<std::complex<double> >(size_t n, const std::complex<double> *src, size_t incsrc, std::complex<double> *dst, size_t incdst){
	RNP_FORTRAN_NAME(zcopy,ZCOPY)(n, src, incsrc, dst, incdst);
}

extern "C" void RNP_FORTRAN_NAME(saxpy,SAXPY)(const integer &n, const float &sa, const float *sx, const integer &incx, float *sy, const integer &incy);
extern "C" void RNP_FORTRAN_NAME(daxpy,dAXPY)(const integer &n, const double &sa, const double *sx, const integer &incx, double *sy, const integer &incy);
extern "C" void RNP_FORTRAN_NAME(caxpy,CAXPY)(const integer &n, const std::complex<float> &sa, const std::complex<float> *sx, const integer &incx, std::complex<float> *sy, const integer &incy);
extern "C" void RNP_FORTRAN_NAME(zaxpy,ZAXPY)(const integer &n, const std::complex<double> &sa, const std::complex<double> *sx, const integer &incx, std::complex<double> *sy, const integer &incy);
template <>
inline void Axpy<float,float>(size_t n, const float &alpha, const float *x, size_t incx, float *y, size_t incy){
	RNP_FORTRAN_NAME(saxpy,SAXPY)(n, alpha, x, incx, y, incy);
}
template <>
inline void Axpy<double,double>(size_t n, const double &alpha, const double *x, size_t incx, double *y, size_t incy){
	RNP_FORTRAN_NAME(daxpy,dAXPY)(n, alpha, x, incx, y, incy);
}
template <>
inline void Axpy<std::complex<float>,std::complex<float> >(size_t n, const std::complex<float> &alpha, const std::complex<float> *x, size_t incx, std::complex<float> *y, size_t incy){
	RNP_FORTRAN_NAME(caxpy,CAXPY)(n, alpha, x, incx, y, incy);
}
template <>
inline void Axpy<std::complex<double>,std::complex<double> >(size_t n, const std::complex<double> &alpha, const std::complex<double> *x, size_t incx, std::complex<double> *y, size_t incy){
	RNP_FORTRAN_NAME(zaxpy,ZAXPY)(n, alpha, x, incx, y, incy);
}
/*
template <class T>
T Dot(size_t n, const T *x, size_t incx, T *y, size_t incy){
	T sum(0);
	while(n --> 0){
		sum += (*x)*(*y);
		x += incx; y += incy;
	}
	return sum;
}

template <class T>
T ConjugateDot(size_t n, const T *x, size_t incx, T *y, size_t incy){
	T sum(0);
	while(n --> 0){
		sum += _RealOrComplexChooser<T>::_conj(*x)*(*y);
		x += incx; y += incy;
	}
	return sum;
}
*/
extern "C" float RNP_FORTRAN_NAME(snrm2,SNRM2)(const integer &n, const float *x, const integer &incx);
extern "C" double RNP_FORTRAN_NAME(dnrm2,DNRM2)(const integer &n, const double *x, const integer &incx);
extern "C" float RNP_FORTRAN_NAME(scnrm2,SCNRM2)(const integer &n, const std::complex<float> *x, const integer &incx);
extern "C" double RNP_FORTRAN_NAME(dznrm2,DZNRM2)(const integer &n, const std::complex<double> *x, const integer &incx);
template <>
inline float Norm2<float>(size_t n, const float *x, size_t incx){
	return RNP_FORTRAN_NAME(snrm2,SNRM2)(n, x, incx);
}
template <>
inline double Norm2<double>(size_t n, const double *x, size_t incx){
	return RNP_FORTRAN_NAME(dnrm2,DNRM2)(n, x, incx);
}
template <>
inline float Norm2<std::complex<float> >(size_t n, const std::complex<float> *x, size_t incx){
	return RNP_FORTRAN_NAME(scnrm2,SCNRM2)(n, x, incx);
}
template <>
inline double Norm2<std::complex<double> >(size_t n, const std::complex<double> *x, size_t incx){
	return RNP_FORTRAN_NAME(dznrm2,DZNRM2)(n, x, incx);
}

extern "C" float RNP_FORTRAN_NAME(sasum,SASUM)(const integer &n, const float *x, const integer &incx);
extern "C" double RNP_FORTRAN_NAME(dasum,DASUM)(const integer &n, const double *x, const integer &incx);
extern "C" float RNP_FORTRAN_NAME(scasum,SCASUM)(const integer &n, const std::complex<float> *x, const integer &incx);
extern "C" double RNP_FORTRAN_NAME(dzasum,DZASUM)(const integer &n, const std::complex<double> *x, const integer &incx);
template <>
inline float Asum<float>(size_t n, const float *x, size_t incx){
	return RNP_FORTRAN_NAME(sasum,SASUM)(n, x, incx);
}
template <>
inline double Asum<double>(size_t n, const double *x, size_t incx){
	return RNP_FORTRAN_NAME(dasum,DASUM)(n, x, incx);
}
template <>
inline float Asum<std::complex<float> >(size_t n, const std::complex<float> *x, size_t incx){
	return RNP_FORTRAN_NAME(scasum,SCASUM)(n, x, incx);
}
template <>
inline double Asum<std::complex<double> >(size_t n, const std::complex<double> *x, size_t incx){
	return RNP_FORTRAN_NAME(dzasum,DZASUM)(n, x, incx);
}

extern "C" integer RNP_FORTRAN_NAME(isamax,ISAMAX)(const integer &n, const float *sx, const integer &incx);
extern "C" integer RNP_FORTRAN_NAME(idamax,IDAMAX)(const integer &n, const double *sx, const integer &incx);
extern "C" integer RNP_FORTRAN_NAME(icamax,ICAMAX)(const integer &n, const std::complex<float> *sx, const integer &incx);
extern "C" integer RNP_FORTRAN_NAME(izamax,IZAMAX)(const integer &n, const std::complex<double> *sx, const integer &incx);
template <>
inline size_t MaximumIndex<float>(size_t n, const float *x, size_t incx){
	return RNP_FORTRAN_NAME(isamax,ISAMAX)(n, x, incx) - 1;
}
template <>
inline size_t MaximumIndex<double>(size_t n, const double *x, size_t incx){
	return RNP_FORTRAN_NAME(idamax,IDAMAX)(n, x, incx) - 1;
}
template <>
inline size_t MaximumIndex<std::complex<float> >(size_t n, const std::complex<float> *x, size_t incx){
	return RNP_FORTRAN_NAME(icamax,ICAMAX)(n, x, incx) - 1;
}
template <>
inline size_t MaximumIndex<std::complex<double> >(size_t n, const std::complex<double> *x, size_t incx){
	return RNP_FORTRAN_NAME(izamax,IZAMAX)(n, x, incx) - 1;
}

//// Level 2

extern "C" void RNP_FORTRAN_NAME(zgemv,ZGEMV)(const char *trans, const integer &m, const integer &n, 
	const std::complex<double> &alpha, const std::complex<double> *a, const integer &lda, const std::complex<double> *x,
	const integer &incx, const std::complex<double> &beta, std::complex<double> *y, const integer &incy);
template <>
template <>
inline MultMV<'N'>::MultMV(size_t m, size_t n, const std::complex<double> &alpha, const std::complex<double> *a, size_t lda,
       const std::complex<double> *x, size_t incx,
       const std::complex<double> &beta, std::complex<double> *y, size_t incy)
{
	RNP_FORTRAN_NAME(zgemv,ZGEMV)("N", m, n, alpha, a, lda, x, incx, beta, y, incy);
}

/*
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
				const T *x2i = xx;
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
						size_t i = j-1;
						while(i --> 0){
							x2 -= incx;
							temp += a[i+j*lda]*(*x2);
						}
					}else{
						if(nounit){ temp *= _RealOrComplexChooser<T>::_conj(a[j+j*lda]); }
						size_t i = j-1;
						while(i --> 0){
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
*/
//// Level 3

extern "C" void RNP_FORTRAN_NAME(zgemm,ZGEMM)(const char *transa, const char *transb, const integer &m, const integer &
	n, const integer &k, const std::complex<double> &alpha, const std::complex<double> *a, const integer &lda, 
	const std::complex<double> *b, const integer &ldb, const std::complex<double> &beta, std::complex<double> *c, const integer &ldc);
	
template <>
template <>
inline MultMM<'N','N'>::MultMM(size_t m, size_t n, size_t k, const double &alpha, const std::complex<double> *a, size_t lda,
	const std::complex<double> *b, size_t ldb, const double &beta, std::complex<double> *c, size_t ldc){
	RNP_FORTRAN_NAME(zgemm,ZGEMM)("N", "N", m,n,k, std::complex<double>(alpha),a,lda, b,ldb, std::complex<double>(beta),c,ldc);
}


/*

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
*/


}; // namespace TBLAS
}; // namespace RNP

#endif // _RNP_TBLAS_EXT_H_
