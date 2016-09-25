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
typedef struct S4_Material_ S4_Material;
typedef struct S4_Layer_ S4_Layer;
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

/**************************** Public API *****************************/



#include "S4_internal.h"

#endif /* S4_H_INCLUDED */
