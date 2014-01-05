#ifndef _RNP_LINEAR_SOLVE_H_
#define _RNP_LINEAR_SOLVE_H_

#include "TBLAS.h"
#include "TLASupport.h"

namespace RNP{

template <char trans='N'>
struct LinearSolve{
	template <class T>
	LinearSolve(size_t n, size_t nRHS, T *a, size_t lda, T *b, size_t ldb, int *info = NULL, size_t *pivots = NULL){
		if(NULL != info){ *info = 0; }
		if(0 == n || nRHS == 0){ return; }
		
		if(pivots){} // prevent unused parameter warning
		
		int iinfo = 0;
		{ // LU decomposition
			for(size_t j = 0; j < n; ++j){
				size_t jp = j + RNP::TBLAS::MaximumIndex(n-j, &a[j+j*lda], 1);
				if(T(0) != a[jp+j*lda]){
					if(jp != j){
						RNP::TBLAS::Swap(n, &a[j+0*lda], lda, &a[jp+0*lda], lda);
						RNP::TBLAS::Swap(nRHS, &b[j+0*lda], ldb, &b[jp+0*lda], ldb);
					}
					if(j < n){
						RNP::TBLAS::Scale(n-j-1, T(1)/a[j+j*lda], &a[j+1+j*lda], 1); // possible overflow when inverting A(j,j)
					}
				}else{
					iinfo = j+1;
				}
				if(j < n){
					RNP::TBLAS::Rank1Update(n-j-1, n-j-1, T(-1), &a[j+1+j*lda], 1, &a[j+(j+1)*lda], lda, &a[j+1+(j+1)*lda], lda);
				}
			}
		}
		if(0 == iinfo){
			if(trans == 'T'){
				RNP::TBLAS::SolveTrM<'L','U','T','N'>(n, nRHS, T(1), a, lda, b, ldb);
				RNP::TBLAS::SolveTrM<'L','L','T','U'>(n, nRHS, T(1), a, lda, b, ldb);
			}else if(trans == 'C'){
				RNP::TBLAS::SolveTrM<'L','U','C','N'>(n, nRHS, T(1), a, lda, b, ldb);
				RNP::TBLAS::SolveTrM<'L','L','C','U'>(n, nRHS, T(1), a, lda, b, ldb);
			}else{
				RNP::TBLAS::SolveTrM<'L','L','N','U'>(n, nRHS, T(1), a, lda, b, ldb);
				RNP::TBLAS::SolveTrM<'L','U','N','N'>(n, nRHS, T(1), a, lda, b, ldb);
			}
		}
		if(NULL != info){ *info = iinfo; }
	}
};

}; // namespace RNP

#ifdef RNP_HAVE_LAPACK
# include "LinearSolve_lapack.h"
#endif

#endif // _RNP_LINEAR_SOLVE_H_
