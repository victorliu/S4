/* Copyright (C) 2009-2016, Victor Liu
 * This file is part of S4
 * Written by Victor Liu (victorliu@alumni.stanford.edu)
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

#ifndef S4_H_INCLUDED
#define S4_H_INCLUDED

/************************* Fundamental types *************************/
typedef double S4_real;

#ifdef __cplusplus
#include <complex>
typedef std::complex<S4_real> S4_complex;
#endif

/***************************** S4 types ******************************/
typedef int S4_MaterialID;
typedef int S4_LayerID;
typedef struct S4_Simulation_ S4_Simulation;

typedef struct S4_Options_{
	int use_discretized_epsilon;

	// Set use_subpixel_smoothing to nonzero if subpixel smoothing should
	// be done on the Fourier transform of the dielectric functions.
	// This will cause the dielectric function to be discretized instead
	// of using the analytic Fourier transforms of specified shapes.
	int use_subpixel_smoothing;

	// Set use_Lanczos_smoothing to nonzero if Lanczos smoothing should
	// be done on the Fourier transform of the dielectric functions.
	int use_Lanczos_smoothing;

	// Set use_polarization_basis to nonzero if a decomposition into the
	// structure conforming polarization basis should be performed. By
	// default, the basis is not everywhere-unit-normalized and purely
	// real unless one of the following two options are used.
	int use_polarization_basis;

	// Set use_jones_vector_basis to nonzero if a complex Jones
	// vector polarization basis should be used. This setting has no
	// effect if use_polarization_basis is zero.
	int use_jones_vector_basis;

	// Set use_normal_vector_basis to nonzero if a normalized normal
	// vector polarization basis should be used. This setting has no
	// effect if use_polarization_basis is zero.
	int use_normal_vector_basis;

	// Set use_normal_vector_field to nonzero if a normal vector field
	// should be generated. The field is everywhere normal to material
	// interfaces and forced to have unit magnitude everywhere.
	int use_normal_vector_field;

	// Set resolution to be the multiple of the reciprocal lattice order
	// extent by which the polarization basis vector field should be
	// discretized. For example, if the reciprocal lattice orders go up
	// to G = {2,1}, then the vector field will be generated on a grid
	// of dimension 2*f x 1*f, where f is this value.
	// It is preferable to set this value to have small integer factors
	// in order to make the FFT more efficient. Powers of 2 are best.
	// The default value is 64, and must be at least 2 to satisfy the
	// Nyquist sampling criterion.
	int resolution;

	// Set vector_field_dump_filename_prefix to non-NULL if the
	// automatically generated vector field in the polarization
	// decomposition should be output to files whose prefix is the
	// specified string (suffix is layer name).
	char *vector_field_dump_filename_prefix;

	// Specifies how G vectors should be chosen (currently not used)
	// 0 = Circular truncation
	// 1 = Parallelogramic truncation
	int lattice_truncation;

	// Specifies how much simulation output progress to output to stdout
	int verbosity;

	// Set use_experimental_fmm to nonzero if the experimental FMM
	// formulation should be used. (For experimental/debug purposes)
	int use_experimental_fmm;

	// Set use_less_memory to nonzero if intermediate matrices should not
	// be stored in memory when possible.
	int use_less_memory;

	S4_real lanczos_smoothing_width;
	int lanczos_smoothing_power;
} S4_Options;

#define S4_MSG_ERROR    1
#define S4_MSG_WARNING  3
#define S4_MSG_INFO     5
#define S4_MSG_STATUS   9
typedef int (*S4_message_handler)(
	void *data, const char *fname, int level, const char *msg
);

/**************************** Public API *****************************/

#ifdef __cplusplus
extern "C" {
#endif

/********************************************************************/
/* Simulation object constructor, destructor, and copy constructor. */
/********************************************************************/
S4_Simulation* S4_Simulation_New(const S4_real *Lr, unsigned int nG, int *G);
void S4_Simulation_Destroy(S4_Simulation *S);
S4_Simulation* S4_Simulation_Clone(const S4_Simulation *S);

S4_message_handler S4_Simulation_SetMessageHandler(
	S4_Simulation *S, S4_message_handler handler, void *data
);

unsigned int S4_Lattice_Count(const S4_real *Lr, unsigned int nG);
int S4_Lattice_Reciprocate(const S4_real *Lr, S4_real *Lk);

/**********************************/
/* Simulation getters and setters */
/**********************************/
int S4_Simulation_SetLattice(S4_Simulation *S, const S4_real *Lr);
int S4_Simulation_GetLattice(const S4_Simulation *S, S4_real *Lr);
int S4_Simulation_SetBases(S4_Simulation *S, unsigned int nG, int *G);
int S4_Simulation_GetBases(const S4_Simulation *S, int *G);
int S4_Simulation_SetFrequency(S4_Simulation *S, const S4_real *freq_complex);
int S4_Simulation_GetFrequency(const S4_Simulation *S, S4_real *freq_complex);
int S4_Simulation_LayerCount(const S4_Simulation *S);
int S4_Simulation_TotalThickness(const S4_Simulation *S, S4_real *thickness);

/******************************/
/* Material related functions */
/******************************/
#define S4_MATERIAL_TYPE_SCALAR_REAL      2
#define S4_MATERIAL_TYPE_SCALAR_COMPLEX   3
#define S4_MATERIAL_TYPE_XYTENSOR_REAL    4
#define S4_MATERIAL_TYPE_XYTENSOR_COMPLEX 5
S4_MaterialID S4_Simulation_SetMaterial(
	S4_Simulation *S, S4_MaterialID M, const char *name, int type, const S4_real *eps
);
S4_MaterialID S4_Simulation_GetMaterialByName(
	const S4_Simulation *S, const char *name
);
int S4_Material_GetName(
	const S4_Simulation *S, S4_MaterialID M, const char **name
);
int S4_Material_GetEpsilon(
	const S4_Simulation *S, S4_MaterialID M, S4_real *eps
);

/***************************/
/* Layer related functions */
/***************************/
S4_LayerID S4_Simulation_SetLayer(
	S4_Simulation *S, S4_LayerID L, const char *name, const S4_real *thickness,
	S4_LayerID copy, S4_MaterialID material
); /*
if NULL == L:
	Adds a new layer with optional name and thickness
if NULL != L:
	if NULL == name, then the name is not changed.
*/
S4_LayerID S4_Simulation_GetLayerByName(
	const S4_Simulation *S, const char *name
);
int S4_Layer_GetName(
	const S4_Simulation *S, S4_LayerID L, const char **name
);
int S4_Layer_GetThickness(
	const S4_Simulation *S, S4_LayerID L, S4_real *thickness
);


/**********************************/
/* Layer region related functions */
/**********************************/
int S4_Layer_ClearRegions(
	S4_Simulation *S, S4_LayerID L
);
#define S4_REGION_TYPE_INTERVAL   11
#define S4_REGION_TYPE_RECTANGLE  12
#define S4_REGION_TYPE_ELLIPSE    13
#define S4_REGION_TYPE_CIRCLE     14
int S4_Layer_SetRegionHalfwidths(
	S4_Simulation *S, S4_LayerID L, S4_MaterialID M,
	int type, const S4_real *halfwidths,
	const S4_real *center, const S4_real *angle_frac
);

#define S4_REGION_TYPE_POLYGON    21
int S4_Layer_SetRegionVertices(
	S4_Simulation *S, S4_LayerID L, S4_MaterialID M,
	int type, int nv, const S4_real *v,
	const S4_real *center, const S4_real *angle_frac
);
int S4_Layer_IsCopy(S4_Simulation *S, S4_LayerID L);

/********************************/
/* Excitation related functions */
/********************************/

int S4_Simulation_ExcitationPlanewave(
	S4_Simulation *S, const S4_real *kdir, const S4_real *udir,
	const S4_real *amp_u, const S4_real *amp_v
);
int S4_Simulation_ExcitationExterior(S4_Simulation *S, int n, const int *exg, const double *ex);
int S4_Simulation_ExcitationDipole(S4_Simulation *S, const double k[2], const char *layer, const double pos[2], const double moment[6]);

/***********************************/
/* Solution hint related functions */
/***********************************/
int S4_Simulation_SolveLayer(S4_Simulation *S, S4_LayerID L);

/****************************/
/* Output related functions */
/****************************/

int S4_Simulation_GetPowerFlux(
	S4_Simulation *S, S4_LayerID layer, const S4_real *offset,
	S4_real *power
);
int S4_Simulation_GetPowerFluxes(
	S4_Simulation *S, S4_LayerID layer, const S4_real *offset,
	S4_real *power
);
// waves should be size 2*11*S->n_G
// Each wave is length 11:
//   { kx, ky, kzr, kzi, ux, uy, uz, cur, cui, cvr, cvi }
// The first n_G waves are forward propagating, the second are backward.
// Each set of n_G waves are ordered in the basis ordering.
int S4_Simulation_GetWaves(S4_Simulation *S, S4_LayerID layer, S4_real *wave);

int S4_Simulation_GetFieldPlane(
	S4_Simulation *S, const int nxy[2], const S4_real *xyz0,
	S4_real *E, S4_real *H
);

int S4_Simulation_GetEpsilon(
	S4_Simulation *S, int nxy[2], const S4_real *xyz0, S4_real *eps
); // eps is {real,imag}

// Returns a solution error code
// Tint is a vector of time averaged stress tensor integral
int S4_Simulation_GetStressTensorIntegral(
	S4_Simulation *S, S4_LayerID layer, const S4_real *offset,
	S4_real *Tint
);

// Returns a solution error code
// which can be 'U', 'E', 'H', 'e'
// 'E' is epsilon*|E|^2, 'H' is |H|^2, 'e' is |E|^2, 'U' is 'E'+'H'
#define S4_VOLUME_INTEGRAL_ENERGY_E  'E'
#define S4_VOLUME_INTEGRAL_ENERGY_H  'H'
#define S4_VOLUME_INTEGRAL_ENERGY    'U'
#define S4_VOLUME_INTEGRAL_E_SQUARED 'e'
int S4_Simulation_GetLayerVolumeIntegral(
	S4_Simulation *S, S4_LayerID layer, int which, S4_real *integral
);
int S4_Simulation_GetLayerZIntegral(
	S4_Simulation *S, S4_LayerID layer, const S4_real *r,
	S4_real *integral
);

/***************************************/
/* Mode/band-solving related functions */
/***************************************/
// Determinant is (rmant[0]+i*rmant[1])*base^expo
int S4_Simulation_GetSMatrixDeterminant(
	S4_Simulation *S, const S4_real *k,
	S4_real *rmant, S4_real *base, int *expo
);


#ifdef __cplusplus
} /* extern "C" */
#endif

#include "S4_internal.h"

#endif /* S4_H_INCLUDED */
