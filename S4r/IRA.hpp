#ifndef IRA_HPP_INCLUDED
#define IRA_HPP_INCLUDED

#include <complex>
#include <cstring>

// Implicitly Restarted Arnoldi method (ARPACK)

// Solves the generalized complex non-symmetric eigenvalue problem
//    A*x = lambda*B*x
// for a few (lambda,x) eigenpairs. The choice of which part of the
// spectrum to obtain is determined by the EigenvalueCompare function.
// The regular non-shift-invert mode only works for extremal parts
// of the spectrum. For interior eigenvalues, it will not converge,
// and shift-invert mode must be used. When using shift-invert mode,
// the comparison function must prioritize large values (e.g., to
// compute eigenvalues nearest a shift, one should use LargestMagnitude.
//
// If shift_invert is set to false, then ApplyOp should do:
//   y <- inv(B)*A*x
// otherwise, when shift_invert is true, then ApplyOp should do:
//   y <- inv(A-shift*B)*x
// In either case, ApplyB should do y <- B*x

namespace IRA{

typedef double real_type;
typedef std::complex<double> complex_type;

class ComplexEigensystemImpl;
class ComplexEigensystem{
	ComplexEigensystemImpl *impl;
public:
	struct Params{
		size_t max_iterations;   // actual number of iterations is returned here
		double tol;              // actual tolerance used is returned here
		bool invariant_subspace; // if true, then v contains invariant subspace instead of eigenvectors
		Params(){
			max_iterations = 1000;
			tol = 0;                 // default to machine precision
			invariant_subspace = false;
		}
	};
	ComplexEigensystem(
		size_t n, size_t n_wanted,
		size_t n_arnoldi, // n_wanted+1 <= n_arnoldi <= n
		bool shift_invert, const complex_type &shift,
		const Params &params,
		complex_type *resid0 = NULL // if not NULL, we commandeer the memory for storing the residual
	);
	virtual ~ComplexEigensystem();
	
	// Solution routines (trigger a solve)
	size_t GetConvergedCount();
	const complex_type *GetEigenvalues();
	const complex_type *GetEigenvectors();
	
	// Common interface routines
	virtual bool IsOpInPlace() const = 0;
	// If IsAInPlace returns true, then y contains the input
	// vector rather than x, and x is NULL.
	virtual void ApplyOp(size_t n, const complex_type *x, complex_type *y) = 0;
	// Returns whether Bop is the identity.
	virtual bool IsBIdentity() const{ return true; }
	// True if the B operator can be applied in place: y <- Bop*y
	// instead of the usual y <- Bop*x
	virtual bool IsBInPlace() const = 0;
	// If IsBIdentity returns true, then ApplyB will never be called.
	// If IsBInPlace returns true, then y contains the input
	// vector rather than x, and x is NULL.
	virtual void ApplyB(size_t n, const complex_type *x, complex_type *y) = 0;
	
	// Return true to use user specified shifts (typically this is not used).
	virtual bool CanSupplyShifts() const{ return false; }
	// Requesting n shifts returned in s. The array s0 contains the current
	// Ritz eigenvalue estimates.
	virtual void GetShifts(size_t n, const complex_type *s0, complex_type *s) const = 0;
	
	// If a is more desireable than b, then Compare(a,b) should return false.
	virtual bool EigenvalueCompare(const complex_type &a, const complex_type &b) const = 0;
	
	static bool LargestMagnitude(const complex_type &a, const complex_type &b){
		return std::abs(a) < std::abs(b);
	}
	static bool SmallestMagnitude(const complex_type &a, const complex_type &b){
		return std::abs(a) < std::abs(b);
	}
	static bool LargestReal(const complex_type &a, const complex_type &b){
		return a.real() < b.real();
	}
};

namespace BLAS{

template <typename TV, typename T>
void Set(size_t m, size_t n, const TV &offdiag, const TV &diag, T* a, size_t lda){
	for(size_t j = 0; j < n; ++j){
		for(size_t i = 0; i < m; ++i){
			a[i+j*lda] = (i == j ? diag : offdiag);
		}
	}
}

template <typename T>
void Copy(size_t m, size_t n, const T* src, size_t ldsrc, T* dst, size_t lddst){
	if(m == ldsrc && m == lddst){
		memcpy(dst, src, sizeof(T) * m*n);
	}else{
		for(size_t j = 0; j < n; ++j){
			memcpy(&dst[j*lddst], &src[j*lddst], sizeof(T) * m);
		}
	}
}

template <typename TS, typename T>
void Copy(size_t n, const TS* src, size_t incsrc, T* dst, size_t incdst){
	if(sizeof(T) == sizeof(TS) && 1 == incsrc && 1 == incdst){
		memcpy(dst, src, sizeof(T) * n);
	}else{
		while(n --> 0){
			*dst = *src;
			src += incsrc; dst += incdst;
		}
	}
}

template <typename T>
T ConjugateDot(size_t n, const T* x, size_t incx, const T* y, size_t incy){
	// Don't call BLAS due to return arg problem
	T sum = 0.f;
	while(n --> 0){
		sum += (std::conj(*x)*(*y));
		x += incx; y += incy;
	}
	return sum;
}

template <typename T>
typename T::value_type Norm2(size_t n, const T* x, size_t incx){
	typedef typename T::value_type real_type;
	static const real_type rzero(0);
	static const real_type rone(1);
	real_type ssq(1), scale(0);
	while(n --> 0){
		if(rzero != std::real(*x)){
			real_type temp = std::abs(std::real(*x));
			if(scale < temp){
				real_type r = scale/temp;
				ssq = ssq*r*r + rone;
				scale = temp;
			}else{
				real_type r = temp/scale;
				ssq += r*r;
			}
		}
		if(rzero != std::imag(*x)){
			real_type temp = std::abs(std::imag(*x));
			if(scale < temp){
				real_type r = scale/temp;
				ssq = ssq*r*r + rone;
				scale = temp;
			}else{
				real_type r = temp/scale;
				ssq += r*r;
			}
		}
		x += incx;
	}
	return scale*std::sqrt(ssq);
}

template <typename TS, typename T>
void Scale(size_t n, const TS &alpha, T* x, size_t incx){
	while(n --> 0){
		*x *= alpha;
		x += incx;
	}
}

template <class S, class T>
void Axpy(size_t n, const S &alpha, const T *x, size_t incx, T *y, size_t incy){
	while(n --> 0){
		(*y) += alpha*(*x);
		x += incx; y += incy;
	}
}

} // namespace BLAS

} // namespace IRA

#endif // IRA_HPP_INCLUDED
