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

#ifndef _S4_H_
#define _S4_H_

#include <stdio.h>

// Note: the interface described in this file must be kept C-compatible.

#ifdef __cplusplus
#include <complex>
extern "C" {
#endif
#include "pattern/pattern.h"
#ifdef __cplusplus
}
#endif


#ifdef ENABLE_S4_TRACE
# ifdef __cplusplus
#  include <cstdio>
# else
#  include <stdio.h>
# endif
# define S4_TRACE(...) fprintf(stderr, __VA_ARGS__)
# define S4_CHECK if(1)
#else
# define S4_TRACE(...)
# define S4_CHECK if(0)
#endif




#define S4_VERB(verb,...) \
	do{ \
		if(S->options.verbosity >= verb){ \
			fprintf(stdout, "[%d] ", verb); \
			fprintf(stdout, __VA_ARGS__); \
		} \
	}while(0)

#ifdef __cplusplus
# include <cstdlib>
#else
# include <stdlib.h>
#endif

void* S4_malloc(size_t size);
void S4_free(void *ptr);

typedef struct Material_{
	char *name;    // name of material
	struct Material_ *next; // linked-list next pointer
	int type; // 0 = scalar epsilon, 1 = tensor
	union{
		double s[2]; // real and imaginary parts of epsilon
		double abcde[10];
		// [ a b 0 ]
		// [ c d 0 ]
		// [ 0 0 e ]
	} eps;
} Material;

typedef struct Layer_{
	char *name;       // name of layer
	double thickness; // thickness of layer
	char *material;   // name of background material
	Pattern pattern;  // See pattern.h
	char *copy;       // See below.
	struct Layer_ *next; // linked-list next pointer
} Layer;
// If a layer is a copy, then `copy' is the name of the layer that should
// be copied, and `material' and `pattern' are inherited, and so they can
// be arbitrary. For non-copy layers, copy should be NULL.

typedef struct Options_{
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

	double lanczos_smoothing_width;
	int lanczos_smoothing_power;
} Options;

typedef struct Solution_{
	int *G; // length 2*n_G; uv coordinates of Lk
	void **layer_bands;    // opaque pointer to the layer band structure
	void **layer_solution; // opaque pointer to the layer solution vector
	double *kx, *ky; // each length n_G
} Solution;

struct FieldCache;

typedef struct Excitation_Planewave_{
	double hx[2],hy[2]; // re,im components of H_x,H_y field
	size_t order;
} Excitation_Planewave;

typedef struct Excitation_Dipole_{
	double pos[2];
	double moment[6]; // dipole moment (jxr, jxi, jyr, jyi, jzr, jzi)
	// The source is at the layer interface after the layer specified by ex_layer.
} Excitation_Dipole;

typedef struct Excitation_Exterior_{
	size_t n;
	int *Gindex1; /* 1-based G index, negative means backwards (length n) */
	double *coeff; /* complex coefficients of each G index (length 2n) */
} Excitation_Exterior;

typedef struct Excitation_{
	union{
		Excitation_Planewave planewave; // type 0
		Excitation_Dipole dipole; // type 1
		Excitation_Exterior exterior; // type 2
	} sub;
	int type;
	char *layer; // name of layer after which excitation is applied
} Excitation;

typedef struct Simulation_{
	double Lr[4]; // real space lattice:
	              //  {Lr[0],Lr[1]} is the first basis vector's x and y coords.
	              //  {Lr[2],Lr[3]} is the second basis vector's x and y coords.
	double Lk[4]; // reciprocal lattice:
	              //  2pi*{Lk[0],Lk[1]} is the first basis vector's x and y coords.
	              //  2pi*{Lk[2],Lk[3]} is the second basis vector's x and y coords.
	              // Computed by taking the inverse of Lr as a 2x2 column-major matrix.
	int n_G;            // Number of G-vectors (Fourier planewave orders).
	Material *material; // linked list of materials
	Layer *layer;       // linked list of layers

	// Excitation
	double omega[2]; // real and imaginary parts of omega
	double k[2]; // xy components of k vector, as fraction of k0 = omega
	Excitation exc;

	Solution *solution; // The solution object is not allocated until needed.
	Options options;

	struct FieldCache *field_cache; // Internal cache of vector field FT when using polarization bases
} Simulation;

#ifdef __cplusplus
extern "C" {
#endif



//// Layer methods
void Layer_Init(Layer *L, const char *name, double thickness, const char *material, const char *copy);
void Layer_Destroy(Layer *L);

//// Material methods
void Material_Init(Material *M, const char *name, const double eps[2]);
void Material_InitTensor(Material *M, const char *name, const double abcde[10]);
void Material_Destroy(Material *M);

//// Simulation methods
void Simulation_Init(Simulation *S);
void Simulation_Destroy(Simulation *S);
void Simulation_Clone(const Simulation *S, Simulation *T);

void Simulation_DestroySolution(Simulation *S);
void Simulation_DestroyLayerSolutions(Simulation *S);

// Destroys the solution belonging to a given simulation and sets
// S->solution to NULL.

///////////////// Simulation parameter specifications /////////////////
// These will invalidate any existing partially computed solutions.

// Assumes that S->Lr is filled in, and computes Lk.
// For 1D lattices, S->Lr[2] = S->Lr[3] = 0, and the corresponding
// Lk[2] and Lk[3] are set to 0.
// Returns:
//  -1 if S is NULL
//   0 on success
//   1 if basis vectors are degenerate (rank Lr < 2)
//   2 if basis vectors are zero (Lr == 0)
int Simulation_MakeReciprocalLattice(Simulation *S);

Material* Simulation_AddMaterial(Simulation *S);
// Adds a blank material to the end of the linked list in S and returns
// a pointer to it.

Layer* Simulation_AddLayer(Simulation *S);
// Adds a blank layer to the end of the linked list in S and returns
// a pointer to it.

int Simulation_SetNumG(Simulation *S, int n);
int Simulation_GetNumG(const Simulation *S, int **G);
double Simulation_GetUnitCellSize(const Simulation *S);

// Returns NULL if the layer name could not be found, or on error.
// If index is non-NULL, the offset of the layer in the list is returned.
Layer* Simulation_GetLayerByName(const Simulation *S, const char *name, int *index);

// Returns NULL if the material name could not be found, or on error.
// If index is non-NULL, the offset of the material in the list is returned.
Material* Simulation_GetMaterialByName(const Simulation *S, const char *name, int *index);
Material* Simulation_GetMaterialByIndex(const Simulation *S, int i);

// Returns -n if the n-th argument is invalid.
// S and L should be valid pointers, material is the offset into the material list (not bounds checked).
// angle should be in radians
// vert should be of length 2*nvert
int Simulation_AddLayerPatternCircle   (Simulation *S, Layer *layer, int material, const double center[2], double radius);
int Simulation_AddLayerPatternEllipse  (Simulation *S, Layer *layer, int material, const double center[2], double angle, const double halfwidths[2]);
int Simulation_AddLayerPatternRectangle(Simulation *S, Layer *layer, int material, const double center[2], double angle, const double halfwidths[2]);
int Simulation_AddLayerPatternPolygon  (Simulation *S, Layer *layer, int material, const double center[2], double angle, int nvert, const double *vert);
int Simulation_RemoveLayerPatterns(Simulation *S, Layer *layer);
int Simulation_ChangeLayerThickness(Simulation *S, Layer *layer, const double *thickness);

// Returns 14 if no layers present
// exg is length 2*n, pairs of G index (1-based index, negative for backwards modes), and 0,1 polarization (E field x,y)
// ex is length 2*n, pairs of re,im complex coefficients
int Simulation_MakeExcitationExterior(Simulation *S, int n, const int *exg, const double *ex);
int Simulation_MakeExcitationPlanewave(Simulation *S, const double angle[2], const double pol_s[2], const double pol_p[2], size_t order);
int Simulation_MakeExcitationDipole(Simulation *S, const double k[2], const char *layer, const double pos[2], const double moment[6]);


// Internal functions
#ifdef __cplusplus
// Field cache manipulation
void Simulation_InvalidateFieldCache(Simulation *S);
std::complex<double>* Simulation_GetCachedField(const Simulation *S, const Layer *layer);
void Simulation_AddFieldToCache(Simulation *S, const Layer *layer, size_t n, const std::complex<double> *P, size_t Plen);
#endif

//////////////////////// Simulation solutions ////////////////////////
// These functions try to compute the minimum amount necessary to
// obtain the desired result. All partial solutions are stored so
// that later requerying of the same results returns the stored
// result.

// Solution error codes:
// -n - n-th argument is invalid
// 0  - successful
// 1  - allocation error
// 9  - lattice n_G not set
// 10 - copy referenced an unknown layer name
// 11 - a copy of a copy was made
// 12 - duplicate layer name found
// 13 - excitation layer name not found
// 14 - no layers
// 15 - material not found


// Sets up S->solution; if one already exists, destroys it and allocates a new one
// Possibly reduces S->n_G
// Returns a solution error code
int Simulation_InitSolution(Simulation *S);

// Returns a solution error code
int Simulation_SolveLayer(Simulation *S, Layer *layer);

// Returns a solution error code
// powers[0] - 0.5 real forw
// powers[1] - 0.5 real back
// powers[2] - imag forw
// powers[3] - imag back
int Simulation_GetPoyntingFlux(Simulation *S, Layer *layer, double offset, double powers[4]);
int Simulation_GetPoyntingFluxByG(Simulation *S, Layer *layer, double offset, double *powers);

// Returns a list of S->n_G complex numbers of mode propagation constants
// q should be length 2*S->n_G
int Simulation_GetPropagationConstants(Simulation *S, Layer *layer, double *q);

// Returns lists of 2*S->n_G complex numbers of forward and backward amplitudes
// forw and back should each be length 4*S->n_G
int Simulation_GetAmplitudes(Simulation *S, Layer *layer, double offset, double *forw, double *back);
// waves should be size 2*11*S->n_G
// Each wave is:
//   { kx, ky, kzr, kzi, ux, uy, uz, cur, cui, cvr, cvi }
int Simulation_GetWaves(Simulation *S, Layer *layer, double *wave);

// Returns a solution error code
// Tint is a vector of time averaged stress tensor integral
int Simulation_GetStressTensorIntegral(Simulation *S, Layer *layer, double offset, double Tint[6]);

// Returns a solution error code
// which can be 'U', 'E', 'H', 'e'
// 'E' is epsilon*|E|^2, 'H' is |H|^2, 'e' is |E|^2, 'U' is 'E'+'H'
int Simulation_GetLayerVolumeIntegral(Simulation *S, Layer *layer, char which, double integral[2]);
int Simulation_GetLayerZIntegral(Simulation *S, Layer *layer, const double r[2], double integral[6]);

// Outputs a POV-Ray render script of a unit cell of the structure.
// Return value can be ignored for valid inputs.
int Simulation_OutputStructurePOVRay(Simulation *S, FILE *fp);

// Outputs a PostScript rendering of the layer pattern to stdout.
// Return value can be ignored for valid inputs.
int Simulation_OutputLayerPatternDescription(Simulation *S, Layer *layer, FILE *fp);

// Returns a solution error code
// Outputs the Fourier reconstruction of the layer pattern to stdout in Gnuplot splot format.
// The unit cell is discretized into nu and nv cells in the lattice directions.
int Simulation_OutputLayerPatternRealization(Simulation *S, Layer *layer, int nu, int nv, FILE *fp);

// Returns a solution error code
// E field is stored as {Exr,Eyr,Ezr,Exi,Eyi,Ezi}
int Simulation_GetField(Simulation *S, const double r[3], double fE[6], double fH[6]);
int Simulation_GetFieldPlane(Simulation *S, int nxy[2], double z, double *E, double *H);

// Returns a solution error code
int Simulation_GetEpsilon(Simulation *S, const double r[3], double eps[2]); // eps is {real,imag}

// Returns a solution error code
// Determinant is (rmant[0]+i*rmant[1])*base^expo
int Simulation_GetSMatrixDeterminant(Simulation *S, double rmant[2], double *base, int *expo);


#ifdef __cplusplus
}
#endif

#endif // _S4_H_
