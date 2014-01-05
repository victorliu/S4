#ifndef S4R_EIGENSYSTEMS_HPP_INCLUDED
#define S4R_EIGENSYSTEMS_HPP_INCLUDED

#include <S4r/Types.hpp>

namespace S4r{

// These routines compute `cols` eigenpairs of the matrix A,
// returning the eigenvalues in vals, and the eigenvectors
// in vecs.

int Eigensystem(
	SpMat::Index cols, const SpMat &A,
	CVec &vals, CMat &vecs
);

int HermitianEigensystem(
	SpMat::Index cols, const SpMat &A,
	CVec &vals, CMat &vecs
);

} // namespace S4r

#endif // S4R_EIGENSYSTEMS_HPP_INCLUDED
