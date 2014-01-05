#ifndef S4R_TYPES_HPP_INCLUDED
#define S4R_TYPES_HPP_INCLUDED

#include <complex>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace S4r{

//// Fundamental types
typedef std::complex<double> doublecomplex;
typedef unsigned int bitfield;

//// Linear algebra objects
// Real and complex matrices
typedef Eigen::Matrix<double       ,Eigen::Dynamic,Eigen::Dynamic> RMat;
typedef Eigen::Matrix<doublecomplex,Eigen::Dynamic,Eigen::Dynamic> CMat;
// Real and complex vectors
typedef Eigen::Matrix<double,       Eigen::Dynamic,1> RVec;
typedef Eigen::Matrix<doublecomplex,Eigen::Dynamic,1> CVec;
// 2D/3D Geometric objects
typedef Eigen::Matrix<double,        2, 1> Vec2;
typedef Eigen::Matrix<double,        2, 2> Mat2;
typedef Eigen::Matrix<double,        3, 1> Vec3;
typedef Eigen::Matrix<doublecomplex, 2, 1> CVec2;
typedef Eigen::Matrix<doublecomplex, 3, 1> CVec3;
// Constitutive tensor
typedef Eigen::Matrix<doublecomplex, 2, 2> CTensor2;
typedef Eigen::Matrix<doublecomplex, 3, 3> CTensor3;
// Sparse matrices
typedef Eigen::SparseMatrix<doublecomplex>          SpMat;
typedef Eigen::Triplet<doublecomplex, SpMat::Index> Triplet;

} // namespace S4r

#endif // S4R_TYPES_HPP_INCLUDED
