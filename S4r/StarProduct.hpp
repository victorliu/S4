#ifndef S4R_STAR_PRODUCT_HPP_INCLUDED
#define S4R_STAR_PRODUCT_HPP_INCLUDED

#include <S4r/Types.hpp>
#include <S4r/Pseudoinverse.hpp>

namespace S4r{

// We will refer to the following block diagram
//   n0                 n1
//   a0 --> [     ] --> a1
//          [     ]
//   b0 <-- [     ] <-- b1

// Given a transfer matrix
//                  n1   n1
//   [ a0 ]   n0 [ T00  T01 ] [ a1 ]
//   [    ] =    [          ] [    ]
//   [ b0 ]   n0 [ T10  T11 ] [ b1 ]
// Converts it into a scattering matrix
//                  n0   n1
//   [ a1 ]   n1 [ S00  S01 ] [ a0 ]
//   [    ] =    [          ] [    ]
//   [ b0 ]   n0 [ S10  S11 ] [ b1 ]
//
// We have
//   S00 =           inv(T00)
//   S01 =          -inv(T00) T01
//   S10 =       T10 inv(T00)
//   S11 = T11 - T10 inv(T00) T01
// 
template <typename MatrixType>
MatrixType T2Smat(typename MatrixType::Index n0, typename MatrixType::Index n1, const MatrixType &T){
	assert(T.rows() == 2*n0 && T.cols() == 2*n1);
	MatrixType S(n0+n1, n0+n1);
	Pseudoinverse(
		n0, n1, T.block(0,0, n0,n1).data(), T.block(0,0, n0,n1).outerstride(),
		S.block(0,0, n1,n0).data(), S.block(0,0,n1,n0).outerStride()
	);
	S.block(n1,0, n0,n0) = T.block(n0,0, n0,n1) * S.block(0,0,n1,n0);
	S.block(0,n0, n1,n1) = -S.block(0,0,n1,n0) * T.block(0,n1, n0,n1);
	S.block(n1,n0, n0,n1) = T.block(n0,n1, n0,n1) + T.block(n0,0, n0,n1) * S.block(0,n0, n1,n1);
	return S;
}

// Same as T2Smat except that T is of the form
//                 n1  n1
//   [ a0 ]   n0 [ Ta  Tb ] [ a1 ]
//   [    ] =    [        ] [    ]
//   [ b0 ]   n0 [ Tb  Ta ] [ b1 ]
template <typename MatrixType>
MatrixType T2Sblocks(typename MatrixType::Index n0, typename MatrixType::Index n1, const MatrixType &Ta, const MatrixType &Tb){
	assert(Ta.rows() == n0 && Tb.rows() == n0 && Ta.cols() == n1 && Tb.cols() == n1);
	MatrixType S(n0+n1, n0+n1);
	
	Pseudoinverse(
		n0, n1, Ta.data(), Ta.outerStride(),
		S.block(0,0, n1,n0).data(), S.block(0,0,n1,n0).outerStride()
	);
	
	S.block(n1,0, n0,n0) = Tb * S.block(0,0,n1,n0);
	S.block(0,n0, n1,n1) = -S.block(0,0,n1,n0) * Tb;
	S.block(n1,n0, n0,n1) = Ta + Tb * S.block(0,n0, n1,n1);
	return S;
}


//                  n0   n1
//   [ a1 ]   n1 [ A00  A01 ] [ a0 ]
//   [    ] =    [          ] [    ]
//   [ b0 ]   n0 [ A10  A11 ] [ b1 ]

//                  n1   n2
//   [ a2 ]   n2 [ B00  B01 ] [ a1 ]
//   [    ] =    [          ] [    ]
//   [ b1 ]   n1 [ B10  B11 ] [ b2 ]

// When we cascade two such matrices A and B, the net result C is given by
//   C00 = B00 inv(I - A01 B10) A00
//   C01 = B00 inv(I - A01 B10) A01 B11 + B01
//   C10 = A10 + A11 inv(I - B10 A01) B10 A00
//   C11 =       A11 inv(I - B10 A01) B11
template <typename MatrixType>
MatrixType StarProduct(
	typename MatrixType::Index n0, typename MatrixType::Index n1, typename MatrixType::Index n2,
	const MatrixType &A, const MatrixType &B
){
	assert(A.rows() == n0+n1 && A.cols() == n0+n1);
	assert(B.rows() == n1+n2 && B.cols() == n1+n2);
	
	MatrixType C(n0+n2, n0+n2);
	
	{
		MatrixType T(n1,n1);
		T.setIdentity();
		T -= A.block(0,n0, n1,n1) * B.block(n2,0, n1,n1);
		Eigen::PartialPivLU<MatrixType> lu(T);
		C.block(0, 0, n2,n0) = B.block(0,0, n2,n1) * lu.solve(A.block(0,0, n1,n0));
		C.block(0,n0, n2,n2) = B.block(0,0, n2,n1) * lu.solve(A.block(0,n0, n1,n1) * B.block(n2,n1,n1,n2)) + B.block(0,n1, n2,n2);
	}
	{
		MatrixType T(n1,n1);
		T.setIdentity();
		T -= B.block(n2,0, n1,n1) * A.block(0,n0, n1,n1);
		Eigen::PartialPivLU<MatrixType> lu(T);
		C.block(n2, 0, n0,n0) = A.block(n1,n0, n0,n1) * lu.solve(B.block(n2,0, n1,n1) * A.block(0,0, n1,n0)) + A.block(n1,0, n0,n0);
		C.block(n2,n0, n0,n2) = A.block(n1,n0, n0,n1) * lu.solve(B.block(n2,n1, n1,n2));
	}
	return C;
}


} // namespace S4r

#endif // S4R_STAR_PRODUCT_HPP_INCLUDED
