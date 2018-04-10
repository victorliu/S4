/* Copyright (C) 2009-2011, Stanford University
 * This file is part of S4
 * Written by Victor Liu (vkl@stanford.edu)
 *
 * S4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * S4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef _RCWA_H_
#define _RCWA_H_

#include <cstddef> // for size_t
#include <complex>

// Possible values for epstype:
#define EPSILON2_TYPE_FULL            0
#define EPSILON2_TYPE_BLKDIAG1        1
#define EPSILON2_TYPE_BLKDIAG2        2
#define EPSILON2_TYPE_BLKDIAG1_SCALAR 5
#define EPSILON2_TYPE_BLKDIAG2_SCALAR 6


// The rcwa.h/rcwa.cpp files contain the core RCWA and S-Matrix logic.
// The details of computing the Fourier expansions of the dielectric
// constant is done elsewhere.

// A layer is described by a thickness and a list of Fourier
// coefficients of its dielectric constant.

// To compute S-matrices, we require the eigenvector and
// eigenvalue matrices.


// Purpose
// =======
// Computes the layer band structure (waveguide modes) of a layer,
// which is assumed to be uniform in the z-direction. The layer's
// dielectric constant is expressed in Fourier space. The eigen-
// modes have z-dependence exp(iqz), where q is the propagation
// constant in the z-direction. Assuming the electric and magnetic
// fields have Fourier expansion coefficients ex, ey, hx, hy, the
// eigenmode amplitudes are related to the field components by
//   [ -ex ]   [   A    -A  ] [ a ]
//   [  ey ] = [            ] [   ] , A = kp phi inv(diag(q))/omega
//   [  hx ]   [   B     B  ] [ b ]   B = phi
//   [  hy ]   [            ] [   ]
// a and b are the mode amplitudes at a particular z-location.
//
// Arguments
// =========
// omega       - (INPUT) The operating angular frequency.
// n           - (INPUT) The number of Fourier orders used in the
//               dielectric expansion. We will denote 2*n as 2n,
//               and 4*n as 4n.
// kx, ky      - (INPUT) Arrays of length n. The x- and y-components
//               of the Fourier k-vectors of the dielectric expansion.
// Epsilon_inv - (INPUT) Matrix of size n x n. The inverse of the
//               Fourier dielectric constant coupling matrix. The non-
//               inverse version is denoted Epsilon:
//                 Epsilon_{ij} = \int_{unit cell}
//                   \epsilon(x,y) \exp[ i(kx[i]-kx[j])x
//                                     + i(ky[i]-ky[j])y] dx dy
// Epsilon2    - (INPUT) Matrix of size 2n x 2n. The in-plane
//               Fourier dielectric constant coupling matrix.
//               This matrix should approximately be
//                 [ Epsilon          ]
//                 [          Epsilon ]
//               but deviates from it when taking into account proper
//               Fourier factorization rules for the normal and
//               tangential electric field continuity relations. The
//               exact definition is that for Fourier expansion
//               coefficients of the E and D fields,
//                 [ -dy ] = Epsilon2 [ -ey ]
//                 [  dx ]            [  ex ]
// q           - (OUTPUT) Array of length 2n. The eigenvalues
//               corresponding to the eigenmodes. Physically, these
//               are the z-direction propagation constants.
// kp          - (OUTPUT) Matrix of size 2n x 2n. The matrix
//                 omega^2 - [  ky ] Epsilon_inv [ ky  -kx ]
//                           [ -kx ]
//               where kx and ky are to be interpreted as diagonal
//               matrices. This quantity is used often enough to
//               merit returning it and storing it.
//               If NULL, then this parameter is not accessed.
// phi         - (OUTPUT) Matrix of size 2n x 2n. The eigenvector
//               matrix corresponding to the eigensystem. These are
//               the eigenvectors to the master equation
//                 Epsilon2 kp - [ kx ] [ kx  ky ] phi = phi q^2
//                               [ ky ]
//               where kx, ky, and q are to be interpreted as
//               diagonal matrices.
// work        - (IN/OUT) Workspace. If NULL or lwork == 0, then the
//               space is internally allocated. If lwork == -1, then
//               the desired workspace is returned in work[0].real().
// rwork       - (WORK) Real valued workspace. Should be length 4n.
//               If NULL, then the space is internally allocated.
// iwork       - (WORK) Integer workspace. Should be length n. If NULL,
//               then the space is internally allocated.
// lwork       - (INPUT) The length of the work array. If -1, then a
//               workspace query is performed and the optimal lwork
//               is returned in work[0].real().
void SolveLayerEigensystem(
	std::complex<double> omega,
	size_t n,
	const double *kx,
	const double *ky,
	const std::complex<double> *Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	const std::complex<double> *Epsilon2, // size (2*glist.n)^2 (dielectric/normal-field matrix)
	int epstype,
	std::complex<double> *q, // length 2*glist.n
	std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix) (optional)
	std::complex<double> *phi, // size (2*glist.n)^2
	std::complex<double> *work = NULL, // length lwork
	double *rwork = NULL, // length 4*n
	size_t lwork = 0 // set to -1 for query into work[0], at least 4*n*n+2*n
);

// Purpose
// =======
// Same as SolveLayerEigensystem, but assumes Epsilon is eps*I.
// It assumed that Epsilon2 = diag(Epsilon, Epsilon).
// See the documentation for SolveLayerEigensystem for argument
// descriptions.
void SolveLayerEigensystem_uniform(
	std::complex<double> omega,
	size_t n,
	const double *kx,
	const double *ky,
	std::complex<double> eps,
	std::complex<double> *q, // length 2*glist.n
	std::complex<double> *kp = NULL, // size (2*glist.n)^2 (k-parallel matrix)
	std::complex<double> *phi = NULL // size (2*glist.n)^2 (optional)
);


// Purpose
// =======
// Multiplies a matrix by the k-parallel matrix (kp), as passed
// to SolveLayerEigensystem.
//
// Arguments
// =========
// trans       - (INPUT) If "C", multiplies by the conjugate transpose.
//                       If "N", multiplies by the original matrix.
// omega       - (INPUT) The operating angular frequency.
// n           - (INPUT) The number of Fourier orders used in the
//               dielectric expansion. We will denote 2*n as 2n,
//               and 4*n as 4n.
// kx, ky      - (INPUT) Arrays of length n. The x- and y-components
//               of the Fourier k-vectors of the dielectric expansion.
// Epsilon_inv - (INPUT) Matrix of size n x n. The inverse of the
//               Fourier dielectric constant coupling matrix. The non-
//               inverse version is denoted Epsilon:
//                 Epsilon_{ij} = \int_{unit cell}
//                   \epsilon(x,y) \exp[ i(kx[i]-kx[j])x
//                                     + i(ky[i]-ky[j])y] dx dy
// Epsilon2_type - (INPUT) Ignored if kp is not NULL. When kp is NULL,
//               and Epsilon2_type is EPSILON2_TYPE_BLKDIAG1_SCALAR or
//               EPSILON2_TYPE_BLKDIAG2_SCALAR, Epsilon_inv is assumed
//               to be a scalar value and only its first element is
//               accessed.
// kp          - (INPUT) Matrix of size 2n x 2n. The matrix
//                 omega^2 - [  ky ] Epsilon_inv [ ky  -kx ]
//                           [ -kx ]
//               where kx and ky are to be interpreted as diagonal
//               matrices.
//               If NULL, then this parameter is not accessed (See
//               Epsilon2_type).
// ncols       - (INPUT) Number of columns of the target matrix.
// a           - (INPUT) Pointer to the first element of the target
//               matrix. Number of rows should be 2*n.
// lda         - (INPUT) Leading dimension of matrix a, >= 2*n.
// y           - (OUTPUT) Pointer to first element of the resulting
//               matrix. (y = kp.a) Dimension should be 2*n by ncols.
// ldy         - (INPUT) Leading dimension of matrix y.

void MultKPMatrix(
	const char *trans,
	const std::complex<double> omega,
	const size_t n,
	const double *kx,
	const double *ky,
	const std::complex<double> *Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	const int Epsilon2_type,
	const std::complex<double> *kp, // kp if exists. Can be NULL
	const size_t ncols,
	const std::complex<double> *a,
	const size_t lda,
	std::complex<double> *y, // same size as a, y <- kp*a
	const size_t ldy
);

// Purpose
// =======
// Initialize an S-Matrix to the identity matrix.
//
// Arguments
// =========
// n - (INPUT) Number of Fourier orders.
// S - (OUTPUT) Matrix of size 4n x 4n. Set to the identity matrix.
void InitSMatrix(
	size_t n,
	std::complex<double> *S // size (4*n)^2
);

// Purpose
// =======
// Computes the S-Matrix for a layer stack.
//
// Arguments
// =========
// nlayers   - (INPUT) Number of layers in the stack.
// n         - (INPUT) Number of Fourier orders.
// kx, ky    - (INPUT) Arrays of length n. The x- and y-components
//             of the Fourier k-vectors of the dielectric expansion.
// thickness - (INPUT) Array of length nlayers. Thickness of each
//             layer, in order.
// q         - (INPUT) Array of length nlayers. Pointers to the q
//             arrays of each layer (computed by SolveLayerEigensystem*
//             functions), in the same order as thickness.
// Epsilon_inv - (INPUT) size n^2 matrix, inverse of dielectric Fourier
//             coupling matrix.
// epstype   - (INPUT) Type code of the Epsilon2 matrix (see above).
// kp        - (INPUT) Array of length nlayers. Pointers to the kp
//             matrices of each layer, in the same order as thickness.
// phi       - (INPUT) Array of length nlayers. Pointers to the phi
//             matrices of each layer, in the same order as thickness.
// S         - (OUTPUT) Matrix of size 4n x 4n. The output S-matrix.
// work      - (IN/OUT) Workspace. If NULL or lwork == 0, then the
//             space is internally allocated. If lwork == -1, then the
//             desired workspace is returned in work[0].real().
// iwork     - (WORK) Integer workspace, of length 2n. If NULL, the the
//             space is internally allocated.
// lwork     - (INPUT) The length of the work array. If -1, then a
//             workspace query is performed and the optimal lwork is
//             returned in work[0].real().
void GetSMatrix( // appends the layers to an existing S matrix
	size_t nlayers,
	size_t n, // glist.n
	const double *kx, const double *ky,
	std::complex<double> omega,
	const double *thickness, // list of thicknesses
	const std::complex<double> **q, // list of q vectors
	const std::complex<double> **Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	int *epstype,
	const std::complex<double> **kp,
	const std::complex<double> **phi,
	std::complex<double> *S, // size (4*n)^2
	std::complex<double> *work = NULL, // length lwork
	size_t *iwork = NULL, // length n2
	size_t lwork = 0 // set to -1 for query into work[0], at least 4*n*(4*n+1)
);


int SolveAll(
	size_t nlayers,
	size_t n, // glist.n
	const double *kx, const double *ky,
	std::complex<double> omega,
	const double *thickness, // list of thicknesses
	const std::complex<double> **q, // list of q vectors
	const std::complex<double> **Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	int *epstype,
	const std::complex<double> **kp,
	const std::complex<double> **phi,
	std::complex<double> *ab, // length 2*n*nlayers
	std::complex<double> *work_ = NULL, // length lwork
	size_t *iwork = NULL, // length n2*nlayers
	size_t lwork = 0 // set to -1 for query into iwork[0], at least 6*n2^2*nlayers
);

// Purpose
// =======
// Given input amplitudes to a layer stack, computes the mode
// amplitudes within any layer. The forward mode amplitudes a are
// assumed to have the phase origin at the start (z=0) of each layer,
// while the backward mode amplitudes b are assumed to have the phase
// origin at the end of the layer (z=thickness).
//
// Arguments
// =========
// nlayers     - (INPUT) Number of layers in the stack.
// which_layer - (INPUT) Offset of which layer in the stack the output
//               ab should be computed. 0 <= which_layer < nlayers.
// n           - (INPUT) Number of Fourier orders.
// kx, ky      - (INPUT) Arrays of length n. The x- and y-components
//               of the Fourier k-vectors of the dielectric expansion.
// thickness   - (INPUT) List of thicknesses. (See GetSMatrix).
// q, kp, phi  - (INPUT) List of layer properties. (See GetSMatrix).
// Epsilon_inv - (INPUT) size n^2 matrix, inverse of dielectric Fourier
//             coupling matrix.
// epstype     - (INPUT) Type code of the Epsilon2 matrix (see above).
// a0, bN      - (INPUT) Input mode amplitudes, each of length 2n. The
//               forward amplitudes a0 are incident from the front of
//               the stack (layer offset 0), while the backward
//               amplitudes are incident from the back (layer offset
//               nlayers-1). If either is NULL, it is assumed to be
//               a zero vector.
// ab          - (OUTPUT) Output mode amplitudes, length 4n. The
//               forward in the first 2n locations, and the backward
//               mode amplitudes are stored in the remaining 2n
//               locations. The forward amplitudes are relative to the
//               front of the layer, while the backward amplitudes are
//               relative to the back of the layer.
// work        - (IN/OUT) Workspace. If NULL or lwork == 0, then the
//               space is internally allocated. If lwork == -1, then
//               the desired workspace is returned in work[0].real().
// iwork       - (WORK) Integer workspace, of length 2n. If NULL, the the
//               space is internally allocated.
// lwork       - (INPUT) The length of the work array. If -1, then a
//               workspace query is performed and the optimal lwork is
//               returned in work[0].real().
int SolveInterior(
	size_t nlayers,
	size_t which_layer,
	size_t n, // glist.n
	const double *kx, const double *ky,
	std::complex<double> omega,
	const double *thickness, // list of thicknesses
	const std::complex<double> **q, // list of q vectors
	const std::complex<double> **Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	int *epstype,
	const std::complex<double> **kp,
	const std::complex<double> **phi,
	std::complex<double> *a0, // length 2*n
	std::complex<double> *bN, // length 2*n
	std::complex<double> *ab, // length 4*n
	std::complex<double> *work_ = NULL, // length lwork
	size_t *iwork = NULL, // length n2
	size_t lwork = 0 // set to -1 for query into work[0], at least 2*(4*n)^2 + 2*(2*n) + 4*n*(4*n+1)
);

//////////////////////// Solution manipulators ////////////////////////

// Purpose
// =======
// Transforms a mode amplitude vector from SolveInterior to a
// particular z-coordinate within the layer.
//
// Arguments
// =========
// n         - (INPUT) Number of Fourier orders.
// q         - (INPUT) Length 2n. The layer eigenvalues.
// thickness - (INPUT) The thickness of the layer.
// dz        - (INPUT) The z-offset from the start of the layer where
//             the vector ab should be translated.
// ab        - (OUTPUT) Length 4n. On entry, the solution vector with
//             mode amplitudes relative to the ends of the layer. On
//             exit, the mode amplitudes at the specified offset.
void TranslateAmplitudes(
	size_t n, // glist.n
	const std::complex<double> *q, // length 2*glist.n
	double thickness,
	double dz,
	std::complex<double> *ab
);

// Purpose
// =======
// Returns the plane integral of the Poynting flux over the unit cell
// over a surface normal to the z-direction. The result is not time
// averaged (lacks a factor of 0.5).
//
// Arguments
// =========
// n          - (INPUT) The number of Fourier orders.
// kx, ky     - (INPUT) Arrays of length n. The x- and y-components
//              of the Fourier k-vectors of the dielectric expansion.
// omega      - (INPUT) The operating angular frequency.
// q, kp, phi - (INPUT) Layer band structure. (See above)
// Epsilon_inv - (INPUT) size n^2 matrix, inverse of dielectric Fourier
//              coupling matrix.
// epstype    - (INPUT) Type code of the Epsilon2 matrix (see above).
// ab         - (INPUT) Length 4n. The mode amplitudes within the layer
//              at a particular z-offset, e.g. returned by
//              TranslateAmplitudes.
// forward,   - (OUTPUT) The complex integrated Poynting flux over the
// backward     the unit cell.
// work       - (WORK) Length 4*2n. If NULL, then the space is
//              internally allocated.
void GetZPoyntingFlux(
	size_t n, // glist.n
	const double *kx, const double *ky,
	std::complex<double> omega,
	const std::complex<double> *q, // length 2*glist.n
	const std::complex<double> *Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	int epstype,
	const std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	const std::complex<double> *phi, // size (2*glist.n)^2
	const std::complex<double> *ab, // length 4*glist.n
	std::complex<double> *forward,
	std::complex<double> *backward,
	std::complex<double> *work = NULL // length 4*n2 or NULL
);

// Purpose
// =======
// Returns the plane integral of the Poynting flux over the unit cell
// over a surface normal to the z-direction for each diffracted order.
// The result is not time averaged (lacks a factor of 0.5).
//
// Arguments
// =========
// n          - (INPUT) The number of Fourier orders.
// kx, ky     - (INPUT) Arrays of length n. The x- and y-components
//              of the Fourier k-vectors of the dielectric expansion.
// omega      - (INPUT) The operating angular frequency.
// q, kp, phi - (INPUT) Layer band structure. (See above)
// Epsilon_inv - (INPUT) size n^2 matrix, inverse of dielectric Fourier
//              coupling matrix.
// epstype    - (INPUT) Type code of the Epsilon2 matrix (see above).
// ab         - (INPUT) Length 4n. The mode amplitudes within the layer
//              at a particular z-offset, e.g. returned by
//              TranslateAmplitudes.
// forward,   - (OUTPUT) Length 2n. The complex integrated Poynting
// backward     flux over the the unit cell for each G-vector.
// work       - (WORK) Length 4*2n. If NULL, then the space is
//              internally allocated.
void GetZPoyntingFluxComponents(
	size_t n, // glist.n
	const double *kx, const double *ky,
	std::complex<double> omega,
	const std::complex<double> *q, // length 2*glist.n
	const std::complex<double> *Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	int epstype,
	const std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	const std::complex<double> *phi, // size (2*glist.n)^2
	const std::complex<double> *ab, // length 4*glist.n
	std::complex<double> *forward,
	std::complex<double> *backward,
	std::complex<double> *work = NULL // length 4*n2 or NULL
);

// Purpose
// =======
// Returns the electric and/or magnetic field at a particular point
// by reconstructing the Fourier series representation of the field.
//
// Arguments
// =========
// n           - (INPUT) The number of Fourier orders.
// kx, ky      - (INPUT) Arrays of length n. The x- and y-components
//               of the Fourier k-vectors of the dielectric expansion.
// omega       - (INPUT) The operating angular frequency.
// q, kp, phi  - (INPUT) Layer band structure. (See above)
// epsilon_inv - (INPUT) Inverse of Epsilon matrix. (See above)
// epstype     - (INPUT) Type code of the Epsilon2 matrix (see above).
// ab          - (INPUT) Length 4n. The mode amplitudes within the layer
//               at the desired z-offset within the layer,
//               e.g. returned by TranslateAmplitudes.
// r           - (INPUT) The in-plane x- and y-coordinates of the point
//               at which the field should be obtained.
// efield,     - (OUTPUT) The complex field at the desired point. Note
// hfield        no time averaging is performed (no factor of 0.5).
// work        - (WORK) Length 8*2n. If NULL, then the space is
//               internally allocated.
void GetFieldAtPoint(
	size_t n, // glist.n
	const double *kx, const double *ky,
	std::complex<double> omega,
	const std::complex<double> *q, // length 2*glist.n
	const std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	const std::complex<double> *phi, // size (2*glist.n)^2
	const std::complex<double> *epsilon_inv, // size (glist.n)^2, non NULL for efield != NULL || kp == NULL
	int epstype,
	const std::complex<double> *ab, // length 4*glist.n
	const double r[2], // coordinates within layer
	std::complex<double> efield[3],
	std::complex<double> hfield[3],
	std::complex<double> *work = NULL
);
void GetFieldOnGrid(
	size_t n, // glist.n
	int *G, // length 2*glist.n, pairs of uv coordinates of Lk
	const double *kx, const double *ky,
	std::complex<double> omega,
	const std::complex<double> *q, // length 2*glist.n
	const std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	const std::complex<double> *phi, // size (2*glist.n)^2
	const std::complex<double> *epsilon_inv, // size (glist.n)^2, non NULL for efield != NULL || kp == NULL
	int epstype,
	const std::complex<double> *ab, // length 4*glist.n
	const size_t nxy[2], // number of points per lattice direction
	const double *xy0, // origin of grid
	std::complex<double> *efield,
	std::complex<double> *hfield
);

// Purpose
// =======
// Returns the plane integral of the stress tensor over the unit cell
// over a surface normal to the z-direction. The result is not time
// averaged (lacks a factor of 0.5).
//
// Arguments
// =========
// n           - (INPUT) The number of Fourier orders.
// kx, ky      - (INPUT) Arrays of length n. The x- and y-components
//                of the Fourier k-vectors of the dielectric expansion.
// omega       - (INPUT) The operating angular frequency.
// q, kp, phi  - (INPUT) Layer band structure. (See above)
// epsilon_inv - (INPUT) Inverse of Epsilon matrix. (See above)
// epsilon2    - (INPUT) Epsilon2 matrix. (See above)
// epstype     - (INPUT) Type code of the Epsilon2 matrix (see above).
// ab          - (INPUT) Length 4n. The mode amplitudes within the layer
//               at a particular z-offset, e.g. returned by
//               TranslateAmplitudes.
// integral    - (OUTPUT) The complex integrated stress tensor over the
//               the unit cell.
// work        - (WORK) Length 8*2n. If NULL, then the space is
//               internally allocated.
void GetZStressTensorIntegral(
	size_t n, // glist.n
	const double *kx, const double *ky,
	std::complex<double> omega,
	const std::complex<double> *q, // length 2*glist.n
	const std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	const std::complex<double> *phi, // size (2*glist.n)^2
	const std::complex<double> *epsilon_inv, // size (glist.n)^2
	const std::complex<double> *epsilon2, // size (2*glist.n)^2
	int epstype,
	const std::complex<double> *ab, // length 4*glist.n
	std::complex<double> integral[3],
	std::complex<double> *work = NULL // length 8*2*glist.n
);

// Purpose
// =======
// Returns the volume integral of a density over the unit cell
// throughout the entire thickness of the cell. The result is not time
// averaged (lacks a factor of 0.5).
//
// Arguments
// =========
// which       - (INPUT) Computes one of:
//                 U : Total energy density:
//                       epsilon |E|^2 + |H|^2
//                 E : Electric energy density:
//                       epsilon |E|^2
//                 E : Magnetic energy density:
//                       |H|^2
//                 e : Electric field density:
//                       |E|^2
//               Any other value is invalid.
// n           - (INPUT) The number of Fourier orders.
// kx, ky      - (INPUT) Arrays of length n. The x- and y-components
//                of the Fourier k-vectors of the dielectric expansion.
// omega       - (INPUT) The operating angular frequency.
// thickness   - (INPUT) Thickness of the layer.
// q, kp, phi  - (INPUT) Layer band structure. (See above)
// epsilon_inv - (INPUT) Inverse of Epsilon matrix. (See above)
// epsilon2    - (INPUT) Epsilon2 matrix. (See above)
// epstype     - (INPUT) Type code of the Epsilon2 matrix (see above).
// ab          - (INPUT) Length 4n. The mode amplitudes within the layer
//               at a particular z-offset, e.g. returned by
//               TranslateAmplitudes.
// integral    - (OUTPUT) The volume integral of the specified quantity.
// work        - (WORK) Length 4n x 4n. If NULL, then the space is
//               internally allocated.
void GetLayerVolumeIntegral(
	char which,
	size_t n, // glist.n
	const double *kx, const double *ky,
	std::complex<double> omega,
	double thickness,
	const std::complex<double> *q, // length 2*glist.n
	const std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	const std::complex<double> *phi, // size (2*glist.n)^2
	const std::complex<double> *epsilon_inv, // size (glist.n)^2
	const std::complex<double> *epsilon2, // size (2*glist.n)^2
	int epstype,
	const std::complex<double> *ab, // length 4*glist.n
	std::complex<double> *integral,
	std::complex<double> *work
);


// Purpose
// =======
// Returns the line integral of the square magnitudes of electric and
// magnetic field components in the depth (z) direction. The result
// is not time averaged (lacks a factor of 0.5).
//
// Arguments
// =========
// n           - (INPUT) The number of Fourier orders.
// kx, ky      - (INPUT) Arrays of length n. The x- and y-components
//                of the Fourier k-vectors of the dielectric expansion.
// omega       - (INPUT) The operating angular frequency.
// thickness   - (INPUT) Thickness of the layer.
// r           - (INPUT) x and y coordinates of location to integrate.
// q, kp, phi  - (INPUT) Layer band structure. (See above)
// epsilon_inv - (INPUT) Inverse of Epsilon matrix. (See above)
// epsilon2    - (INPUT) Epsilon2 matrix. (See above)
// epstype     - (INPUT) Type code of the Epsilon2 matrix (see above).
// ab          - (INPUT) Length 4n. The mode amplitudes within the layer
//               at a particular z-offset, e.g. returned by
//               TranslateAmplitudes.
// integral    - (OUTPUT) The integrals, in order, of
//               |Ex|^2, |Ey|^2, |Ez|^2, |Hx|^2, |Hy|^2, |Hz|^2
// work        - (WORK) Length 12 x 4n. If NULL, then the space is
//               internally allocated.
void GetLayerZIntegral(
	size_t n, // glist.n
	const double *kx,
	const double *ky,
	std::complex<double> omega,
	const double &thickness, const double r[2],
	const std::complex<double> *q, // length 2*glist.n
	const std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	const std::complex<double> *phi, // size (2*glist.n)^2
	const std::complex<double> *epsilon_inv, // size (glist.n)^2
	const std::complex<double> *epsilon2, // size (2*glist.n)^2
	int epstype,
	const std::complex<double> *ab, // length 4*glist.n
	double integral[6],
	std::complex<double> *work
);

#endif // _RCWA_H_
