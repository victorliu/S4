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

#include "config.h"

#define _USE_MATH_DEFINES
#include "S4.h"
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include <float.h>
#include <TBLAS.h>
#ifdef HAVE_BLAS
# include <TBLAS_ext.h>
#endif
#include <TLASupport.h>
#include <LinearSolve.h>
#ifdef HAVE_LAPACK
# include <LinearSolve_lapack.h>
#endif
#include "rcwa.h"
#include "fmm/fmm.h"
extern "C" {
#include "gsel.h"
}

//#include <IO.h>
#include <cstdio>

#include <numalloc.h>


void* S4_malloc(size_t size){ // for debugging
	void* ret = malloc_aligned(size, 16);
	return ret;
}
void S4_free(void *ptr){
	free_aligned(ptr);
}
void* S4_realloc(void *ptr, size_t size){ // for debugging
	void* ret = realloc_aligned(ptr, size, 16);
	return ret;
}


static double geom_norm3d(const double v[3]){
	double a[3] = {std::abs(v[0]),std::abs(v[1]),std::abs(v[2])};
	double w = a[0];
	if(a[1] > w){ w = a[1]; }
	if(a[2] > w){ w = a[2]; }
	if(0 == w){
		return a[0] + a[1] + a[2];
	}else{
		a[0] /= w; a[1] /= w; a[2] /= w;
		w *= std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
		return w;
	}
}
static double geom_dot3d(const double *u, const double *v){
	return (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]);
}
static double geom_normalize3d(double *v){
	const double n = geom_norm3d(v);
	v[0] /= n; v[1] /= n; v[2] /= n;
	return n;
}
static void geom_cross3d(const double *a, const double *b, double *result){
	result[0] = a[1]*b[2] - a[2]*b[1];
	result[1] = a[2]*b[0] - a[0]*b[2];
	result[2] = a[0]*b[1] - a[1]*b[0];
}


struct LayerModes{
	std::complex<double> *q; // length 2*glist.n
	std::complex<double> *kp; // size (2*glist.n)^2 (k-parallel matrix)
	std::complex<double> *phi; // size (2*glist.n)^2
	std::complex<double> *Epsilon2; // size (2*glist.n)^2 (dielectric/normal-field matrix)
	std::complex<double> *Epsilon_inv; // size (glist.n)^2 inverse of usual dielectric Fourier coupling matrix
	// max total size needed: 2n+13nn
	int epstype;
};
struct Solution_{
	std::complex<double> *ab;
	int *solved;
};

// This structure caches the Fourier transform of the polarization basis
// field, allowing a substantial increase in speed when doing frequency
// scans. Invalidates on SetNumG, and layer patterning
struct FieldCache{
	int n;
	const S4_Layer *layer;
	std::complex<double> *P; // 2n x 2n matrix, allocated along with this structure itself
	FieldCache *next;
};

// Private functions

// If layer_solution is null, only computes layer_modes if it is non-NULL.
// If layer_solution is not null, all layer's modes are computed.
int Simulation_GetLayerSolution(S4_Simulation *S, S4_Layer *layer, LayerModes **layer_modes, std::complex<double> **layer_solution);

// These two assume S->solution is set up already
int Simulation_ComputeLayerSolution(S4_Simulation *S, S4_Layer *L, LayerModes **layer_modes, std::complex<double> **layer_solution);
int Simulation_ComputeLayerModes(S4_Simulation *S, S4_Layer *L, LayerModes **layer_modes);
void Simulation_SetExcitationType(S4_Simulation *S, int type);
void Simulation_CopyExcitation(const S4_Simulation *from, S4_Simulation *to);
int Simulation_GetSMatrix(S4_Simulation *S, int layer_from, int layer_to, std::complex<double> *M);

// Field cache manipulation
void Simulation_InvalidateFieldCache(S4_Simulation *S);
std::complex<double>* Simulation_GetCachedField(S4_Simulation *S, const S4_Layer *layer);
void Simulation_AddFieldToCache(S4_Simulation *S, const S4_Layer *layer, size_t n, const std::complex<double> *P, size_t Plen);

void Layer_Destroy(S4_Layer *L){
	S4_TRACE("> Layer_Destroy(L=%p)\n", L);
	if(NULL == L){
		S4_TRACE("< Layer_Destroy (failed; L == NULL)\n");
		return;
	}
	free(L->name); L->name = NULL;
	if(NULL != L->pattern.shapes){
		for(int i = 0; i < L->pattern.nshapes; ++i){
			if(POLYGON == L->pattern.shapes[i].type && NULL != L->pattern.shapes[i].vtab.polygon.vertex){
				S4_free(L->pattern.shapes[i].vtab.polygon.vertex);
			}
		}
		free(L->pattern.shapes);
		L->pattern.shapes = NULL;
	}
	if(NULL != L->pattern.parent){ free(L->pattern.parent); L->pattern.parent = NULL; }
	Simulation_DestroyLayerModes(L);
	S4_TRACE("< Layer_Destroy\n");
}
void Material_Destroy(S4_Material *M){
	S4_TRACE("> Material_Destroy(M=%p)\n", M);
	if(NULL == M){
		S4_TRACE("< Material_Destroy (failed; M == NULL)\n");
		return;
	}
	S4_TRACE("< Material_Destroy\n");
}


S4_Simulation* S4_Simulation_New(const S4_real *Lr, unsigned int nG, int *G){
	S4_TRACE("> S4_Simulation_New(Lr=%p, nG=%u, G=%p)\n", Lr, nG, G);
	S4_Simulation *S = (S4_Simulation*)malloc(sizeof(S4_Simulation));
	S->solution = NULL; // needed by other initializations

	if(NULL == Lr){
		S->Lr[0] = 1;
		S->Lr[1] = 0;
		S->Lr[2] = 0;
		S->Lr[3] = 1;
	}else{
		memcpy(S->Lr, Lr, sizeof(S4_real) * 4);
	}
	int ret = S4_Lattice_Reciprocate(S->Lr, S->Lk);
	/*
	if(1 == ret){
		if(NULL != S->msg){
			S->msg(S->msgdata, "S4_Lattice_Reciprocate", S4_MSG_ERROR, "Degenerate lattice basis");
		}
	}else if(2 == ret){
		if(NULL != S->msg){
			S->msg(S->msgdata, "S4_Lattice_Reciprocate", S4_MSG_ERROR, "Both lattice vectors are zero");
		}
	}
	*/
	if(nG < 1){ nG = 1; }
	S->n_G = nG;
	S->G = (int*)S4_malloc(sizeof(int) * 2*nG);
	S->n_materials = 0;
	S->n_materials_alloc = 4;
	S->material = (S4_Material*)malloc(sizeof(S4_Material) * S->n_materials_alloc);
	S->n_layers = 0;
	S->n_layers_alloc = 4;
	S->layer = (S4_Layer*)malloc(sizeof(S4_Layer) * S->n_layers_alloc);
	S->omega[0] = 1;
	S->omega[1] = 0;
	S->k[0] = 0;
	S->k[1] = 0;
	S->exc.type = 0;
	S->exc.layer = NULL;
	S->exc.sub.planewave.hx[0] = 0;
	S->exc.sub.planewave.hx[1] = 0;
	S->exc.sub.planewave.hy[0] = 0;
	S->exc.sub.planewave.hy[1] = 0;
	S->exc.sub.planewave.order = 0;
	S->exc.sub.planewave.backwards = 0;

	S->options.use_discretized_epsilon = 0;
	S->options.use_subpixel_smoothing = 0;
	S->options.use_Lanczos_smoothing = 0;
	S->options.use_polarization_basis = 0;
	S->options.use_jones_vector_basis = 0;
	S->options.use_normal_vector_basis = 0;
	S->options.resolution = 8;
	S->options.use_normal_vector_field = 0;
	S->options.vector_field_dump_filename_prefix = NULL;
	S->options.lattice_truncation = 0;
	S->options.verbosity = 0;
	S->options.use_experimental_fmm = 0;
	S->options.use_less_memory = 0;

	S->options.lanczos_smoothing_width = 1.0;
	S->options.lanczos_smoothing_power = 1;

	S->field_cache = NULL;
	
	S->msg = NULL;
	S->msgdata = NULL;

	if(NULL != G){
		memcpy(S->G, G, sizeof(int) * 2*S->n_G);
	}else{
		// Get G vectors
		if(0 != S->Lr[2] || 0 != S->Lr[3]){
			unsigned int NG = S->n_G;
			G_select(S->options.lattice_truncation, &NG, S->Lk, S->G);
			S->n_G = NG;
		}else{
			// 1D lattice
			S->G[0] = 0; S->G[1] = 0;
			int remaining = (S->n_G-1)/2;
			S->n_G = 1+2*remaining;
			for(int i = 0; i < remaining; ++i){
				S->G[2+4*i+0] = i+1;
				S->G[2+4*i+1] = 0;
				S->G[2+4*i+2] = -(i+1);
				S->G[2+4*i+3] = 0;
			}
		}
		S4_VERB(1, "Using %d G-vectors\n", S->n_G);
	}
	S->kx = (double*)S4_malloc(sizeof(double)*2*S->n_G);
	S->ky = S->kx + S->n_G;

	S4_TRACE("< S4_Simulation_New\n");
	return S;
}

unsigned int S4_Lattice_Count(const S4_real *Lr, unsigned int nG){
	S4_real Lk[4];
	S4_Lattice_Reciprocate(Lr, Lk);
	// Get G vectors
	if(nG <= 1){ return nG; }
	if(0 != Lr[2] || 0 != Lr[3]){
		unsigned int NG = nG;
		G_select(0, &NG, Lk, NULL);
		return NG;
	}else{
		if(0 == (nG % 2)){
			nG--;
		}
		return nG;
	}
}

void S4_Simulation_Destroy(S4_Simulation *S){
	S4_TRACE("> S4_Simulation_Destroy(S=%p) [omega=%f]\n", S, S->omega[0]);
	if(NULL == S){
		S4_TRACE("< S4_Simulation_Destroy (failed; S == NULL) [omega=%f]\n", S->omega[0]);
		return;
	}
	if(NULL != S->solution){
		Simulation_DestroySolution(S);
	}
	for(int i = 0; i < S->n_layers; ++i){
		Layer_Destroy(&S->layer[i]);
	}
	free(S->layer);
	for(int i = 0; i < S->n_materials; ++i){
		Material_Destroy(&S->material[i]);
	}
	free(S->material);
	Simulation_SetExcitationType(S, -1);
	Simulation_InvalidateFieldCache(S);
	if(NULL != S->options.vector_field_dump_filename_prefix){
		free(S->options.vector_field_dump_filename_prefix);
		S->options.vector_field_dump_filename_prefix = NULL;
	}
	S4_free(S->kx);
	S4_free(S->G);
	free(S);
	S4_TRACE("< S4_Simulation_Destroy [omega=%f]\n", S->omega[0]);
}
S4_Simulation* S4_Simulation_Clone(const S4_Simulation *S){
	S4_TRACE("> S4_Simulation_Clone(S=%p) [omega=%f]\n", S, S->omega[0]);
	S4_Simulation *T = (S4_Simulation*)malloc(sizeof(S4_Simulation));
	if(NULL == S || NULL == T){
		S4_TRACE("< S4_Simulation_Clone (failed; S == NULL or T == NULL) [omega=%f]\n", S->omega[0]);
		if(NULL != S->msg){
			S->msg(S->msgdata, "S4_Simulation_Clone", S4_MSG_ERROR, "S == NULL or T == NULL");
		}
		return NULL;
	}

	memcpy(T, S, sizeof(S4_Simulation));

	T->n_materials_alloc = S->n_materials_alloc;
	T->material = (S4_Material*)malloc(sizeof(S4_Material) * T->n_materials_alloc);
	for(int i = 0; i < S->n_materials; ++i){
		const S4_Material *M = &(S->material[i]);
		S4_Simulation_SetMaterial(T, -1, M->name, M->type, &M->eps.abcde[0]);
	}

	T->n_layers_alloc = S->n_layers_alloc;
	T->layer = (S4_Layer*)malloc(sizeof(S4_Layer) * T->n_layers_alloc);
	for(int i = 0; i < S->n_layers; ++i){
		const S4_Layer *L = &(S->layer[i]);
		S4_LayerID id = S4_Simulation_SetLayer(T, -1, L->name, &L->thickness, L->copy, L->material);
		S4_Layer *L2 = &T->layer[id];
		// Copy pattern
		L2->pattern.nshapes = L->pattern.nshapes;
		L2->pattern.shapes = (shape*)malloc(sizeof(shape)*L->pattern.nshapes);
		memcpy(L2->pattern.shapes, L->pattern.shapes, sizeof(shape)*L->pattern.nshapes);
		L2->pattern.parent = NULL;
		L2->modes = NULL;
	}

	Simulation_CopyExcitation(S, T);

	T->solution = NULL;
	T->field_cache = NULL;

	S4_TRACE("< S4_Simulation_Clone [omega=%f]\n", S->omega[0]);
	return T;
}

S4_message_handler S4_Simulation_SetMessageHandler(
	S4_Simulation *S, S4_message_handler handler, void *data
){
	if(NULL == S){ return NULL; }
	S->msg = handler;
	S->msgdata = data;
}

int S4_Simulation_SetLattice(S4_Simulation *S, const S4_real *Lr){
	if(NULL == S){ return -1; }
	if(NULL == Lr){ return -2; }
	Simulation_DestroySolution(S);
	memcpy(S->Lr, Lr, sizeof(S4_real) * 4);
	int ret = S4_Lattice_Reciprocate(S->Lr, S->Lk);
	if(1 == ret){
		if(NULL != S->msg){
			S->msg(S->msgdata, "S4_Lattice_Reciprocate", S4_MSG_ERROR, "Degenerate lattice basis");
		}
	}else if(2 == ret){
		if(NULL != S->msg){
			S->msg(S->msgdata, "S4_Lattice_Reciprocate", S4_MSG_ERROR, "Both lattice vectors are zero");
		}
	}
	return ret;
}
int S4_Simulation_GetLattice(const S4_Simulation *S, S4_real *Lr){
	if(NULL == S){ return -1; }
	if(NULL == Lr){ return -2; }
	memcpy(Lr, S->Lr, sizeof(S4_real) * 4);
	return 0;
}
int S4_Simulation_SetBases(S4_Simulation *S, unsigned int nG, int *G){
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(nG < 1){ ret = -2; }
	if(0 != ret){
		return ret;
	}

	Simulation_DestroySolution(S);
	Simulation_InvalidateFieldCache(S);

	S->n_G = nG;
	S->G = (int*)S4_realloc(S->G, sizeof(int)*2*S->n_G);
	S->kx = (double*)S4_realloc(S->kx, sizeof(double)*2*S->n_G);
	S->ky = S->kx + S->n_G;

	if(NULL != G){
		memcpy(S->G, G, sizeof(int) * 2 * S->n_G);
	}else{
		if(0 != S->Lr[2] || 0 != S->Lr[3]){
			unsigned int NG = S->n_G;
			G_select(S->options.lattice_truncation, &NG, S->Lk, S->G);
			S->n_G = NG;
		}else{
			// 1D lattice
			S->G[0] = 0; S->G[1] = 0;
			int remaining = (S->n_G-1)/2;
			S->n_G = 1+2*remaining;
			for(int i = 0; i < remaining; ++i){
				S->G[2+4*i+0] = i+1;
				S->G[2+4*i+1] = 0;
				S->G[2+4*i+2] = -(i+1);
				S->G[2+4*i+3] = 0;
			}
		}
	}
	if(NULL != S->msg && nG != (unsigned int)S->n_G){
		char buffer[64];
		snprintf(buffer, 64, "Using %u bases", S->n_G);
		S->msg(S->msgdata, "S4_Simulation_SetBases", S4_MSG_INFO, buffer);
	}
	return 0;
}
int S4_Simulation_GetBases(const S4_Simulation *S, int *G){
	if(NULL == S){ return -1; }
	if(NULL != G){
		memcpy(G, S->G, sizeof(int) * 2 * S->n_G);
	}
	return S->n_G;
}
int S4_Simulation_SetFrequency(S4_Simulation *S, const S4_real *freq_complex){
	if(NULL == S){ return -1; }
	if(NULL == freq_complex){ return -2; }
	S4_TRACE("> S4_Simulation_SetFrequency(S=%p, freq=(%f,%f))\n", S, freq_complex[0], freq_complex[1]);
	S4_Simulation_DestroyLayerModes(S, -1);
	Simulation_DestroySolution(S);
	Simulation_InvalidateFieldCache(S);
	S->omega[0] = 2*M_PI*freq_complex[0];
	S->omega[1] = 2*M_PI*freq_complex[1];
	S4_TRACE("< S4_Simulation_SetFrequency\n");
	return 0;
}
int S4_Simulation_GetFrequency(const S4_Simulation *S, S4_real *freq_complex){
	if(NULL == S){ return -1; }
	if(NULL == freq_complex){ return -2; }
	freq_complex[0] = S->omega[0] / (2*M_PI);
	freq_complex[1] = S->omega[1] / (2*M_PI);
	return 0;
}
int S4_Simulation_LayerCount(const S4_Simulation *S){
	if(NULL == S){ return -1; }
	return S->n_layers;
}
int S4_Simulation_TotalThickness(const S4_Simulation *S, S4_real *thickness){
	if(NULL == S){ return -1; }
	if(NULL == thickness){ return -2; }
	*thickness = 0;
	for(int i = 0; i < S->n_layers; ++i){
		*thickness += S->layer[i].thickness;
	}
	return 0;
}

int S4_Lattice_Reciprocate(const S4_real *Lr, S4_real *Lk){
	S4_TRACE("> S4_Lattice_Reciprocate(Lr=%f,%f, %f,%f)\n", Lr[0], Lr[1], Lr[2], Lr[3]);
	double d;

	d = Lr[0]*Lr[3] - Lr[1]*Lr[2];
	if(0 == d){
		if(0 != Lr[2] || 0 != Lr[3]){ // 1D lattice has zero vector for second basis vector
			S4_TRACE("< S4_Lattice_Reciprocate (failed; degenerate lattice basis)\n");
			return 1; // degenerate lattice basis vectors
		}
		d = hypot(Lr[0], Lr[1]);
		if(0 == d){
			S4_TRACE("< S4_Lattice_Reciprocate (failed; both lattice vectors are zero)\n");
			return 2; // both basis vectors are zero
		}
		d = 1./(d*d);
		Lk[0] = Lr[0] * d;
		Lk[1] = Lr[1] * d;
		Lk[2] = 0; Lk[3] = 0;
	}else{
		/*if(d < 0){
			double t;
			t = Lr[0]; Lr[0] = Lr[2]; Lr[2] = t;
			t = Lr[1]; Lr[1] = Lr[3]; Lr[3] = t;
			d = -d;
		}*/
		d = 1./d;
		Lk[0] =  d*Lr[3];
		Lk[1] = -d*Lr[2];
		Lk[2] = -d*Lr[1];
		Lk[3] =  d*Lr[0];
	}
	S4_TRACE("< S4_Lattice_Reciprocate(Lk=%f,%f, %f,%f)\n", Lk[0], Lk[1], Lk[2], Lk[3]);
	return 0;
}

S4_MaterialID S4_Simulation_SetMaterial(
	S4_Simulation *S, S4_MaterialID id, const char *name, int type, const S4_real *eps
){
	int newmat = 0;
	S4_TRACE(
		"> S4_Simulation_SetMaterial(S=%p, id=%d, name=%p (%s), type=%d, eps=%p (%f))\n",
		S, id, name, (NULL != name ? name : ""), type, eps, (NULL != eps ? eps[0] : 0)
	);
	if(NULL == S){
		S4_TRACE("< S4_Simulation_SetMaterial (failed; S == NULL)\n");
		return -1;
	}
	S4_Material *M = NULL;
	if(id < 0){
		if(NULL != S->msg){
			S->msg(S->msgdata, "S4_Simulation_SetMaterial", S4_MSG_INFO, "Adding new material");
		}
		if(S->n_materials >= S->n_materials_alloc){
			 S->n_materials_alloc *= 2;
			 S->material = (S4_Material*)realloc(S->material, sizeof(S4_Material) * S->n_materials_alloc);
		}
		id = S->n_materials;
		M = &S->material[id];
		M->name = NULL;
		M->type = 0;
		M->eps.s[0] = 1;
		M->eps.s[1] = 0;

		S->n_materials++;
		newmat = 1;
	}else{
		if(NULL != S->msg){
			S->msg(S->msgdata, "S4_Simulation_SetMaterial", S4_MSG_INFO, "Updating existing material");
		}
		if(id >= S->n_materials){ return -1; }
		M = &S->material[id];
	}
	if(NULL != name){
		if(NULL != M->name){
			free(M->name);
		}
		M->name = strdup(name);
	}
	if(NULL != eps){
		switch(type){
		case S4_MATERIAL_TYPE_SCALAR_REAL:
			M->eps.s[0] = eps[0];
			M->eps.s[1] = 0;
			M->type = 0;
			break;
		case S4_MATERIAL_TYPE_SCALAR_COMPLEX:
			M->eps.s[0] = eps[0];
			M->eps.s[1] = eps[1];
			M->type = 0;
			break;
		case S4_MATERIAL_TYPE_XYTENSOR_REAL:
			for(int i = 0; i < 5; ++i){
				M->eps.abcde[2*i] = eps[i];
			}
			M->type = 1;
			break;
		case S4_MATERIAL_TYPE_XYTENSOR_COMPLEX:
			for(int i = 0; i < 10; ++i){
				M->eps.abcde[i] = eps[i];
			}
			M->type = 1;
			break;
		default:
			if(newmat){
				S->n_materials--;
			}
			M = NULL;
			break;
		}
	}
	S4_TRACE("< S4_Simulation_SetMaterial (returning id=%d)\n", id);
	return id;
}
S4_MaterialID S4_Simulation_GetMaterialByName(
	const S4_Simulation *S, const char *name
){
	S4_TRACE("> S4_Simulation_GetMaterialByName(S=%p, name=%p (%s)) [omega=%f]\n",
		S, name, (NULL == name ? "" : name), S->omega[0]);
	if(NULL == name || NULL == S){ return -1; }
	for(int i = 0; i < S->n_materials; ++i){
		if(0 == strcmp(S->material[i].name, name)){
			S4_TRACE("< S4_Simulation_GetMaterialByName returning %s [omega=%f]\n", S->material[i].name, S->omega[0]);
			return i;
		}
	}
	if(NULL != S->msg){
		char buffer[64];
		snprintf(buffer, 64, "Material not found: %s", name);
		S->msg(S->msgdata, "S4_Simulation_GetMaterialByName", S4_MSG_ERROR, buffer);
	}
	S4_TRACE("< S4_Simulation_GetMaterialByName (failed; material name not found) [omega=%f]\n", S->omega[0]);
	return -1;
}
int S4_Material_GetName(
	const S4_Simulation *S, S4_MaterialID id, const char **name
){
	if(NULL == S){ return -1; }
	if(id < 0 || id >= S->n_materials){ return -2; }
	*name = S->material[id].name;
	return 0;
}
int S4_Material_GetEpsilon(
	const S4_Simulation *S, S4_MaterialID id, S4_real *eps
){
	if(NULL == S){ return -1; }
	if(id < 0 || id >= S->n_materials){ return -2; }
	const S4_Material *M = &S->material[id];
	for(int i = 0; i < 18; ++i){
		eps[i] = 0;
	}
	if(0 == M->type){
		for(int i = 0; i < 3; ++i){
			eps[8*i+0] = M->eps.s[0];
			eps[8*i+1] = M->eps.s[1];
		}
	}else{
		eps[0] = M->eps.abcde[0];
		eps[1] = M->eps.abcde[1];
		eps[2] = M->eps.abcde[4];
		eps[3] = M->eps.abcde[5];
		eps[6] = M->eps.abcde[2];
		eps[7] = M->eps.abcde[3];
		eps[8] = M->eps.abcde[6];
		eps[9] = M->eps.abcde[7];
		eps[16] = M->eps.abcde[8];
		eps[17] = M->eps.abcde[9];
	}
	return 0;
}

S4_LayerID S4_Simulation_SetLayer(
	S4_Simulation *S, S4_LayerID id, const char *name, const S4_real *thickness,
	S4_LayerID copy, S4_MaterialID material
){
	S4_TRACE(
		"> S4_Simulation_SetLayer(S=%p, id=%d, name=%p (%s), thickness=%f, copy=%d, material=%d)\n",
		S, id, name, (NULL != name ? name : ""), (NULL != thickness ? *thickness : 0), copy, material
	);
	if(NULL == S){
		S4_TRACE("< S4_Simulation_SetLayer (failed; S == NULL)\n");
		return -1;
	}
	if(copy < 0 && material < 0){ return -1; }
	if(copy >= S->n_layers){
		if(NULL != S->msg){
			S->msg(S->msgdata, "S4_Simulation_SetLayer", S4_MSG_ERROR, "Attempt to copy non-existent layer");
		}
		return -1;
	}
	if(copy >= 0 && S->layer[copy].copy >= 0){
		if(NULL != S->msg){
			S->msg(S->msgdata, "S4_Simulation_SetLayer", S4_MSG_ERROR, "Attempt to copy a copied layer; de-referencing");
		}
		copy = S->layer[copy].copy;
	}
	Simulation_DestroySolution(S);
	S4_Layer *L = NULL;
	if(id < 0){
		if(NULL != S->msg){
			S->msg(S->msgdata, "S4_Simulation_SetLayer", S4_MSG_INFO, "Adding new layer");
		}
		if(S->n_layers >= S->n_layers_alloc){
			S->n_layers_alloc *= 2;
			S->layer = (S4_Layer*)realloc(S->layer, sizeof(S4_Layer) * S->n_layers_alloc);
		}
		id = S->n_layers;
		L = &S->layer[id];
		S->n_layers++;

		L->name = NULL;
		L->thickness = 0;
		L->material = -1;
		L->copy = -1;
		L->pattern.nshapes = 0;
		L->pattern.shapes = NULL;
		L->pattern.parent = NULL;
		L->modes = NULL;
	}else{
		if(NULL != S->msg){
			S->msg(S->msgdata, "S4_Simulation_SetLayer", S4_MSG_INFO, "Updating existing layer");
		}
		if(id >= S->n_layers){ return -1; }
		L = &S->layer[id];
	}

	if(NULL != name){
		if(NULL != L->name){
			free(L->name);
		}
		L->name = strdup(name);
	}
	if(copy >= 0){
		if(L->copy >= 0){
			S4_Layer_ClearRegions(S, id);
		}
		L->copy = copy;
		L->material = -1;
	}else if(material >= 0){
		if(L->material >= 0){
			S4_Layer_ClearRegions(S, id);
		}
		L->material = material;
		L->copy = -1;
	}
	if(NULL != thickness){
		L->thickness = *thickness;
	}
	S4_TRACE("< S4_Simulation_SetLayer (returning id=%d)\n", id);
	return id;
}
S4_LayerID S4_Simulation_GetLayerByName(
	const S4_Simulation *S, const char *name
){
	S4_TRACE("> S4_Simulation_GetLayerByName(S=%p, name=%p (%s)) [omega=%f]\n",
		S, name, (NULL == name ? "" : name), S->omega[0]);
	if(NULL == S || NULL == name){ return -1; }
	for(int i = 0; i < S->n_layers; ++i){
		if(0 == strcmp(S->layer[i].name, name)){
			S4_TRACE("< S4_Simulation_GetLayerByName returning %s [omega=%f]\n", S->layer[i].name, S->omega[0]);
			return i;
		}
	}
	if(NULL != S->msg){
		char buffer[64];
		snprintf(buffer, 64, "Layer not found: %s", name);
		S->msg(S->msgdata, "S4_Simulation_GetLayerByName", S4_MSG_ERROR, buffer);
	}
	S4_TRACE("< S4_Simulation_GetLayerByName (failed; material name not found) [omega=%f]\n", S->omega[0]);
	return -1;
}
int S4_Layer_GetName(
	const S4_Simulation *S, S4_LayerID id, const char **name
){
	if(NULL == S){ return -1; }
	if(id < 0 || id >= S->n_layers){ return -2; }
	*name = S->layer[id].name;
	return 0;
}
int S4_Layer_GetThickness(
	const S4_Simulation *S, S4_LayerID id, S4_real *thickness
){
	if(NULL == S){ return -1; }
	if(id < 0 || id >= S->n_layers){ return -2; }
	*thickness = S->layer[id].thickness;
	return 0;
}
int S4_Layer_ClearRegions(
	S4_Simulation *S, S4_LayerID id
){
	if(NULL == S){ return -1; }
	if(id < 0 || id >= S->n_layers){ return -2; }
	S4_Layer *L = &S->layer[id];
	L->pattern.nshapes = 0;
	if(NULL != L->pattern.shapes){
		free(L->pattern.shapes);
		L->pattern.shapes = NULL;
	}
	if(NULL != L->pattern.parent){
		free(L->pattern.parent);
		L->pattern.parent = NULL;
	}
	Simulation_DestroyLayerModes(L);
	return 0;
}

int S4_Layer_SetRegionHalfwidths(
	S4_Simulation *S, S4_LayerID Lid, S4_MaterialID Mid,
	int type, const S4_real *halfwidths,
	const S4_real *center, const S4_real *angle_frac
){
	S4_TRACE("> S4_Layer_SetRegionHalfwidths(S=%p, layer=%d, material=%d, halfwidths=%p (%f,%f), center=%p (%f,%f), angle_frac=%f)\n",
		S, Lid, Mid,
		halfwidths, (NULL == halfwidths ? 0 : halfwidths[0]), (NULL == halfwidths ? 0 : halfwidths[1]),
		center, (NULL == center ? 0 : center[0]), (NULL == center ? 0 : center[1]),
		(NULL == angle_frac ? 0 : *angle_frac));
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(Lid < 0 || Lid >= S->n_layers){ ret = -2; }
	if(Mid < 0 || Mid >= S->n_materials){ ret = -3; }
	if(NULL == halfwidths){ ret = -5; }
	if(0 != ret){
		S4_TRACE("< S4_Layer_SetRegionHalfwidths (failed; ret = %d)\n", ret);
		return ret;
	}
	S4_real hw[2] = { halfwidths[0], halfwidths[1] };
	{
		S4_real Lr[4];
		S4_Simulation_GetLattice(S, Lr);
		int lattice1d = 0;
		if(0 == Lr[1] && 0 == Lr[2] && 0 == Lr[3]){
			lattice1d = 1;
		}
		if(lattice1d && (0 != hw[1])){
			if(NULL != S->msg){
				S->msg(S->msgdata, "S4_Layer_SetRegionHalfwidths", S4_MSG_WARNING, "1D lattice has nonzero second halfwidth; ignoring");
			}
			hw[1] = 0;
		}
		if(lattice1d && (S4_REGION_TYPE_INTERVAL != type)){
			if(NULL != S->msg){
				S->msg(S->msgdata, "S4_Layer_SetRegionHalfwidths", S4_MSG_ERROR, "1D lattice may only have interval regions");
			}
			ret = -4;
			S4_TRACE("< S4_Layer_SetRegionHalfwidths (failed; ret = %d)\n", ret);
			return ret;
		}
		if((hw[0] < 0) || (!lattice1d && hw[1] < 0)){
			if(NULL != S->msg){
				S->msg(S->msgdata, "S4_Layer_SetRegionHalfwidths", S4_MSG_WARNING, "Halfwidth is negative; taking absolute value");
			}
			hw[0] = -hw[0];
		}
	}
	S4_Layer *L = &S->layer[Lid];

	Simulation_DestroyLayerModes(L);
	Simulation_DestroySolution(S);
	Simulation_InvalidateFieldCache(S);

	int n = L->pattern.nshapes++;
	L->pattern.shapes = (shape*)realloc(L->pattern.shapes, sizeof(shape)*L->pattern.nshapes);
	if(NULL == L->pattern.shapes){ return 1; }
	shape *sh = &L->pattern.shapes[n];

	if(NULL != center){
		sh->center[0] = center[0];
		sh->center[1] = center[1];
	}
	if(NULL != angle_frac){
		sh->angle = 2*M_PI*(*angle_frac);
	}
	sh->tag = Mid;
	switch(type){
	case S4_REGION_TYPE_INTERVAL:
		sh->type = RECTANGLE;
		sh->vtab.rectangle.halfwidth[0] = hw[0];
		sh->vtab.rectangle.halfwidth[1] = 0;
		break;
	case S4_REGION_TYPE_RECTANGLE:
		sh->type = RECTANGLE;
		sh->vtab.rectangle.halfwidth[0] = hw[0];
		sh->vtab.rectangle.halfwidth[1] = hw[1];
		break;
	case S4_REGION_TYPE_ELLIPSE:
		if(halfwidths[0] == halfwidths[1]){
			sh->type = CIRCLE;
			sh->vtab.circle.radius = hw[0];
		}else{
			sh->type = ELLIPSE;
			sh->vtab.ellipse.halfwidth[0] = hw[0];
			sh->vtab.ellipse.halfwidth[1] = hw[1];
		}
		break;
	case S4_REGION_TYPE_CIRCLE:
		sh->type = CIRCLE;
		sh->vtab.circle.radius = hw[0];
		break;
	default:
		L->pattern.nshapes--;
		S4_TRACE("  S4_Layer_SetRegionHalfwidths unknown type: %d\n", type);
		break;
	}

	S4_TRACE("< S4_Layer_SetRegionHalfwidths\n");
	return 0;
}

int S4_Layer_SetRegionVertices(
	S4_Simulation *S, S4_LayerID Lid, S4_MaterialID Mid,
	int type, int nv, const S4_real *v,
	const S4_real *center, const S4_real *angle_frac
){
	S4_TRACE("> S4_Layer_SetRegionVertices(S=%p, layer=%d, material=%d, nv=%d, v=%p, center=%p (%f,%f), angle_frac=%f)\n",
		S, Lid, Mid,
		nv, v,
		center, (NULL == center ? 0 : center[0]), (NULL == center ? 0 : center[1]),
		(NULL == angle_frac ? 0 : *angle_frac));
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(Lid < 0 || Lid >= S->n_layers){ ret = -2; }
	if(Mid < 0 || Mid >= S->n_materials){ ret = -3; }
	if(nv < 2){ ret = -5; }
	if(NULL == v){ ret = -7; }
	if(0 != ret){
		S4_TRACE("< S4_Layer_SetRegionVertices (failed; ret = %d)\n", ret);
		return ret;
	}
	S4_Layer *L = &S->layer[Lid];

	Simulation_DestroyLayerModes(L);
	Simulation_DestroySolution(S);
	Simulation_InvalidateFieldCache(S);

	int n = L->pattern.nshapes++;
	L->pattern.shapes = (shape*)realloc(L->pattern.shapes, sizeof(shape)*L->pattern.nshapes);
	if(NULL == L->pattern.shapes){ return 1; }
	shape *sh = &L->pattern.shapes[n];

	if(NULL != center){
		sh->center[0] = center[0];
		sh->center[1] = center[1];
	}
	if(NULL != angle_frac){
		sh->angle = 2*M_PI* (*angle_frac);
	}
	sh->tag = Mid;
	switch(type){
	case S4_REGION_TYPE_POLYGON:
		sh->type = POLYGON;
		sh->vtab.polygon.n_vertices = nv;
		sh->vtab.polygon.vertex = (double*)S4_malloc(sizeof(double)*nv*2);
		for(int i = 0; i < nv; ++i){
			sh->vtab.polygon.vertex[2*i+0] = v[2*i+0];
			sh->vtab.polygon.vertex[2*i+1] = v[2*i+1];
		}
		{
			double area2 = 0;
			int p, q;
			for(p=n-1, q=0; q < n; p = q++){
				area2 += v[2*p+0]*v[2*q+1] - v[2*q+0]*v[2*p+1];
			}
			if(area2 < 0){
				if(NULL != S->msg){
					S->msg(S->msgdata, "S4_Layer_SetRegionVertices", S4_MSG_WARNING, "Polygon has negative orientation; reversing orientation");
				}
				for(int i = 0, j = nv-1; i < nv; ++i, --j){
					sh->vtab.polygon.vertex[2*i+0] = v[2*j+0];
					sh->vtab.polygon.vertex[2*i+1] = v[2*j+1];
				}
			}
		}
		break;
	default:
		L->pattern.nshapes--;
		break;
	}

	S4_TRACE("< S4_Layer_SetRegionVertices\n");
	return 0;
}
int S4_Layer_IsCopy(S4_Simulation *S, S4_LayerID L){
	if(NULL == S){ return -1; }
	if(L < 0 || L >= S->n_layers){ return -2; }
	return (S->layer[L].copy >= 0) ? 1 : 0;
}

void Simulation_DestroyLayerModes(S4_Layer *layer){
	if(NULL != layer->modes){
		if(NULL != layer->modes->q){ S4_free(layer->modes->q); }
		layer->modes->q = NULL;
		free(layer->modes);
		layer->modes = NULL;
	}
}
void S4_Simulation_DestroyLayerModes(S4_Simulation *S, S4_LayerID id){
	if(id < 0){
		for(int i = 0; i < S->n_layers; ++i){
			Simulation_DestroyLayerModes(&S->layer[i]);
		}
	}else if(id < S->n_layers){
		Simulation_DestroyLayerModes(&S->layer[id]);
	}
}

int Simulation_RemoveLayerPatterns(S4_Simulation *S, S4_Layer *layer){
	S4_TRACE("> Simulation_RemoveLayerPatterns(S=%p, layer=%p) [omega=%f]\n",
		S, layer, S->omega[0]);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(0 != ret){
		S4_TRACE("< Simulation_RemoveLayerPatterns (failed; ret = %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	Simulation_DestroySolution(S);

	layer->pattern.nshapes = 0;
	if(NULL != layer->pattern.shapes){
		free(layer->pattern.shapes);
		layer->pattern.shapes = NULL;
	}
	if(NULL != layer->pattern.parent){
		free(layer->pattern.parent);
		layer->pattern.parent = NULL;
	}
	Simulation_DestroyLayerModes(layer);

	S4_TRACE("< Simulation_RemoveLayerPatterns [omega=%f]\n", S->omega[0]);
	return 0;
}
int Simulation_ChangeLayerThickness(S4_Simulation *S, S4_Layer *layer, const double *thick){
	S4_TRACE("> Simulation_ChangeLayerThickness(S=%p, layer=%p, thick=%g) [omega=%f]\n",
		S, layer, *thick, S->omega[0]);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(NULL == thick || *thick < 0){ ret = -3; }
	if(0 != ret){
		S4_TRACE("< Simulation_RemoveLayerPatterns (failed; ret = %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	layer->thickness = *thick;
	Simulation_DestroyLayerSolutions(S);

	S4_TRACE("< Simulation_ChangeLayerThickness [omega=%f]\n", S->omega[0]);
	return 0;
}
int Simulation_SetNumG(S4_Simulation *S, int n){
	int ret = 0;
	S4_TRACE("> Simulation_SetNumG(S=%p, n=%d) [omega=%f]\n", S, n, S->omega[0]);
	if(NULL == S){ ret = -1; }
	if(n < 1){ ret = -2; }
	if(0 != ret){
		S4_TRACE("< Simulation_SetNumG (failed; ret = %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	Simulation_DestroySolution(S);
	Simulation_InvalidateFieldCache(S);

	S->n_G = n;
	S->G = (int*)S4_realloc(S->G, sizeof(int)*2*S->n_G);
	S->kx = (double*)S4_realloc(S->kx, sizeof(double)*2*S->n_G);
	S->ky = S->kx + S->n_G;

	if(0 != S->Lr[2] || 0 != S->Lr[3]){
		unsigned int NG = S->n_G;
		G_select(S->options.lattice_truncation, &NG, S->Lk, S->G);
		S->n_G = NG;
	}else{
		// 1D lattice
		S->G[0] = 0; S->G[1] = 0;
		int remaining = (S->n_G-1)/2;
		S->n_G = 1+2*remaining;
		for(int i = 0; i < remaining; ++i){
			S->G[2+4*i+0] = i+1;
			S->G[2+4*i+1] = 0;
			S->G[2+4*i+2] = -(i+1);
			S->G[2+4*i+3] = 0;
		}
	}
	S4_TRACE("< Simulation_SetNumG [omega=%f]\n", S->omega[0]);
	return 0;
}
int Simulation_GetNumG(const S4_Simulation *S, int **G){
	int ret = 0;
	S4_TRACE("> Simulation_GetNumG(S=%p, G=%p) [omega=%f]\n", S, G, S->omega[0]);
	if(NULL == S){ ret = -1; }
	if(0 != ret){
		S4_TRACE("< Simulation_GetNumG (failed; ret = %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	ret = S->n_G;

	if(NULL != G){
		*G = S->G;
	}

	S4_TRACE("< Simulation_GetNumG [omega=%f]\n", S->omega[0]);
	return ret;
}
int Simulation_AddLayerPatternCircle(
	S4_Simulation *S,
	S4_Layer *layer, int material,
	const double center[2],
	double radius
){
	S4_TRACE("> Simulation_AddLayerPatternCircle(S=%p, layer=%p, material=%d, center=%p (%f,%f), radius=%f)\n",
		S, layer,
		material,
		center, (NULL == center ? 0.0 : center[0]), (NULL == center ? 0.0 : center[1]),
		radius);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(material < 0){ ret = -3; }
	if(NULL == center){ ret = -4; }
	if(radius < 0){ ret = -5; }
	if(0 != ret){
		S4_TRACE("< Simulation_AddLayerPatternCircle (failed; ret = %d)\n", ret);
		return ret;
	}

	Simulation_DestroyLayerModes(layer);
	Simulation_DestroySolution(S);
	Simulation_InvalidateFieldCache(S);

	int n = layer->pattern.nshapes++;
	layer->pattern.shapes = (shape*)realloc(layer->pattern.shapes, sizeof(shape)*layer->pattern.nshapes);
	if(NULL == layer->pattern.shapes){ return 1; }
	shape *sh = &layer->pattern.shapes[n];
	sh->type = CIRCLE;
	sh->center[0] = center[0];
	sh->center[1] = center[1];
	sh->angle = 0;
	sh->vtab.circle.radius = radius;
	sh->tag = material;

	S4_TRACE("< Simulation_AddLayerPatternCircle\n");
	return 0;
}
int Simulation_AddLayerPatternEllipse(
	S4_Simulation *S,
	S4_Layer *layer, int material,
	const double center[2],
	double angle,
	const double halfwidths[2]
){
	S4_TRACE("> Simulation_AddLayerPatternEllipse(S=%p, layer=%p, material=%d, center=%p (%f,%f), angle=%f, halfwidths=%p (%f,%f))\n",
		S, layer,
		material,
		center, (NULL == center ? 0 : center[0]), (NULL == center ? 0 : center[1]),
		angle,
		halfwidths, (NULL == halfwidths ? 0 : halfwidths[0]), (NULL == halfwidths ? 0 : halfwidths[1]));
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(material < 0){ ret = -3; }
	if(NULL == center){ ret = -4; }
	if(NULL == halfwidths){ ret = -6; }
	if(0 != ret){
		S4_TRACE("< Simulation_AddLayerPatternEllipse (failed; ret = %d)\n", ret);
		return ret;
	}

	Simulation_DestroyLayerModes(layer);
	Simulation_DestroySolution(S);
	Simulation_InvalidateFieldCache(S);

	int n = layer->pattern.nshapes++;
	layer->pattern.shapes = (shape*)realloc(layer->pattern.shapes, sizeof(shape)*layer->pattern.nshapes);
	if(NULL == layer->pattern.shapes){ return 1; }
	shape *sh = &layer->pattern.shapes[n];
	sh->type = ELLIPSE;

	sh->center[0] = center[0];
	sh->center[1] = center[1];
	sh->angle = angle;
	sh->vtab.ellipse.halfwidth[0] = halfwidths[0];
	sh->vtab.ellipse.halfwidth[1] = halfwidths[1];
	sh->tag = material;

	S4_TRACE("< Simulation_AddLayerPatternEllipse\n");
	return 0;
}

int Simulation_AddLayerPatternRectangle(
	S4_Simulation *S,
	S4_Layer *layer, int material,
	const double center[2],
	double angle,
	const double halfwidths[2]
){
	S4_TRACE("> Simulation_AddLayerPatternRectangle(S=%p, layer=%p, material=%d, center=%p (%f,%f), angle=%f, halfwidths=%p (%f,%f))\n",
		S, layer,
		material,
		center, (NULL == center ? 0 : center[0]), (NULL == center ? 0 : center[1]),
		angle,
		halfwidths, (NULL == halfwidths ? 0 : halfwidths[0]), (NULL == halfwidths ? 0 : halfwidths[1]));
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(material < 0){ ret = -3; }
	if(NULL == center){ ret = -4; }
	if(NULL == halfwidths){ ret = -6; }
	if(0 != ret){
		S4_TRACE("< Simulation_AddLayerPatternRectangle (failed; ret = %d)\n", ret);
		return ret;
	}

	Simulation_DestroyLayerModes(layer);
	Simulation_DestroySolution(S);
	Simulation_InvalidateFieldCache(S);

	int n = layer->pattern.nshapes++;
	layer->pattern.shapes = (shape*)realloc(layer->pattern.shapes, sizeof(shape)*layer->pattern.nshapes);
	if(NULL == layer->pattern.shapes){ return 1; }
	shape *sh = &layer->pattern.shapes[n];
	sh->type = RECTANGLE;
	sh->center[0] = center[0];
	sh->center[1] = center[1];
	sh->angle = angle;
	sh->vtab.rectangle.halfwidth[0] = halfwidths[0];
	sh->vtab.rectangle.halfwidth[1] = halfwidths[1];
	sh->tag = material;

	S4_TRACE("< Simulation_AddLayerPatternRectangle\n");
	return 0;
}
int Simulation_AddLayerPatternPolygon(
	S4_Simulation *S,
	S4_Layer *layer, int material,
	const double center[2],
	double angle,
	int nvert,
	const double *vert
){
	S4_TRACE("> Simulation_AddLayerPatternPolygon(S=%p, layer=%p, material=%d, center=%p (%f,%f), angle=%f, nvert=%d, vert=%p)\n",
		S, layer,
		material,
		center, (NULL == center ? 0 : center[0]), (NULL == center ? 0 : center[1]),
		angle,
		nvert,
		vert);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(material < 0){ ret = -3; }
	if(NULL == center){ ret = -4; }
	if(nvert < 3){ ret = -6; }
	if(NULL == vert){ ret = -7; }
	if(0 != ret){
		S4_TRACE("< Simulation_AddLayerPatternPolygon (failed; ret = %d)\n", ret);
		return ret;
	}

	Simulation_DestroyLayerModes(layer);
	Simulation_DestroySolution(S);
	Simulation_InvalidateFieldCache(S);

	int n = layer->pattern.nshapes++;
	layer->pattern.shapes = (shape*)realloc(layer->pattern.shapes, sizeof(shape)*layer->pattern.nshapes);
	if(NULL == layer->pattern.shapes){ return 3; }
	shape *sh = &layer->pattern.shapes[n];
	sh->type = POLYGON;
	sh->center[0] = center[0];
	sh->center[1] = center[1];
	sh->angle = angle;
	sh->vtab.polygon.n_vertices = nvert;
	sh->vtab.polygon.vertex = (double*)S4_malloc(sizeof(double)*nvert*2);
	for(int i = 0; i < nvert; ++i){
		sh->vtab.polygon.vertex[2*i+0] = vert[2*i+0];
		sh->vtab.polygon.vertex[2*i+1] = vert[2*i+1];
	}
	sh->tag = material;

	S4_TRACE("< Simulation_AddLayerPatternPolygon\n");
	return 0;
}

// Returns:
// -n if n-th argument is invalid
// 1  - allocation error
// 9  - Lattice n_G not set
// 10 - copy referenced an unknown layer name
// 11 - a copy of a copy was made
// 12 - duplicate layer name found
// 13 - excitation layer name not found
// 14 - no layers
// 15 - material not found
// 16 - invalid 1D layer patterning
int Simulation_InitSolution(S4_Simulation *S){
	S4_TRACE("> Simulation_InitSolution(S=%p) [omega=%f]\n", S, S->omega[0]);

	if(NULL == S){
		S4_TRACE("< Simulation_InitSolution (failed; S == NULL) [omega=%f]\n", S->omega[0]);
		return -1;
	}
	if(S->n_G < 1){
		if(NULL != S->msg){
			S->msg(S->msgdata, "Simulation_InitSolution", S4_MSG_ERROR, "Basis set contains fewer than 1 element");
		}
		S4_TRACE("< Simulation_InitSolution (failed; S->n_G < 1) [omega=%f]\n", S->omega[0]);
		return 9;
	}

	// Check that every material referenced exists // no, we are sure due to way we added patterns

	// Check that layer copies and the excitation layer names exist
	bool found_ex_layer = (NULL == S->exc.layer);
	for(int i = 0; i < S->n_layers; ++i){
		S4_Layer *L = &(S->layer[i]);
		if(L->copy >= 0){
			bool found = false;
			for(int j = 0; j < S->n_layers; ++j){
				const S4_Layer *L2 = &(S->layer[j]);
				if(j == L->copy){
					if(L2->copy >= 0){
						S4_TRACE("< Simulation_InitSolution (failed; layer %s is referenced by a copy but is also a copy) [omega=%f]\n", L2->name, S->omega[0]);
						return 11;
					}
					found = true;
					break;
				}
			}
			if(!found){
				S4_TRACE("< Simulation_InitSolution (failed; could not find layer %d is referenced by a copy) [omega=%f]\n", L->copy, S->omega[0]);
				return 10;
			}
		}else{
			/*
			// check that no duplicate names exist
			for(int j = 0; j < i; ++j){
				const S4_Layer *L2 = &(S->layer[j]);
				if(0 == strcmp(L2->name, L->name)){
					S4_TRACE("< Simulation_InitSolution (failed; layer name %s appears more than once) [omega=%f]\n", L->name, S->omega[0]);
					return 12;
				}
			}
			*/
		}
		if(!found_ex_layer && L == S->exc.layer){
			found_ex_layer = true;
		}
		// check that if we have a 1D pattern, the only shapes are rectangles
		if(0 == S->Lr[2] && 0 == S->Lr[3]){
			for(int k = 0; k < L->pattern.nshapes; ++k){
				if(RECTANGLE != L->pattern.shapes[k].type){
					return 16;
				}
			}
		}
		// Initialize the layer pattern
		if(NULL != L->pattern.parent){
			free(L->pattern.parent);
		}
		L->pattern.parent = (int*)malloc(sizeof(int)*L->pattern.nshapes);
		int error = Pattern_GetContainmentTree(&L->pattern);
		if(0 != error){
			S4_TRACE("< Simulation_InitSolution (failed; Pattern_GetContainmentTree returned %d for layer %s) [omega=%f]\n", error, L->name, S->omega[0]);
			return error;
		}
	}
	if(S->n_layers < 1){
		S4_TRACE("< Simulation_InitSolution (failed; fewer than one layer found) [omega=%f]\n", S->omega[0]);
		if(NULL != S->msg){
			S->msg(S->msgdata, "Simulation_InitSolution", S4_MSG_ERROR, "Fewer than one layer found");
		}
		return 14;
	}
	if(!found_ex_layer){
		S4_TRACE("< Simulation_InitSolution (failed; excitation layer not found) [omega=%f]\n", S->omega[0]);
		if(NULL != S->msg){
			S->msg(S->msgdata, "Simulation_InitSolution", S4_MSG_ERROR, "Excitation layer not found");
		}
		return 13;
	}

	if(NULL != S->solution){
		Simulation_DestroySolution(S);
	}
	S->solution = (Solution_*)malloc(sizeof(Solution_));
	S->solution->ab = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>) * S->n_layers * 4 * S->n_G);
	S->solution->solved = (int*)malloc(sizeof(int) * S->n_layers);
	if(NULL == S->solution || NULL == S->solution->ab || NULL == S->solution->solved){
		S4_TRACE("< Simulation_InitSolution (failed; could not allocate S->solution) [omega=%f]\n", S->omega[0]);
		if(NULL != S->msg){
			S->msg(S->msgdata, "Simulation_InitSolution", S4_MSG_ERROR, "Could not allocate solution memory");
		}
		return 1;
	}
	memset(S->solution->solved, 0, sizeof(int) * S->n_layers);

	S4_TRACE("I  Simulation_InitSolution G: (%d) [omega=%f]\n", S->n_G, S->omega[0]);

	for(int i = 0; i < S->n_G; ++i){
		S->kx[i] = S->k[0]*S->omega[0] + 2*M_PI*(S->Lk[0]*S->G[2*i+0] + S->Lk[2]*S->G[2*i+1]);
		S->ky[i] = S->k[1]*S->omega[0] + 2*M_PI*(S->Lk[1]*S->G[2*i+0] + S->Lk[3]*S->G[2*i+1]);
	}

	S4_TRACE("< Simulation_InitSolution [omega=%f]\n", S->omega[0]);
	return 0;
}

// realistically, can return -2 if layer is not found, or an InitSolution code
int Simulation_GetLayerSolution(S4_Simulation *S, S4_Layer *layer, LayerModes **layer_modes, std::complex<double> **layer_solution){
	S4_TRACE("> Simulation_GetLayerSolution(S=%p, layer=%p, layer_modes=%p (%p), layer_solution=%p (%p)) [omega=%f]\n",
		S, layer,
		layer_modes, (NULL != layer_modes ? *layer_modes : NULL), layer_solution, (NULL != layer_solution ? *layer_solution : NULL), S->omega[0]);
	int error;
	Solution_ *sol;
	if(NULL == S){
		S4_TRACE("< Simulation_GetLayerSolution (failed; S == NULL) [omega=%f]\n", S->omega[0]);
		return -1;
	}
	if(NULL == layer){
		S4_TRACE("< Simulation_GetLayerSolution (failed; layer == NULL) [omega=%f]\n", S->omega[0]);
		return -2;
	}

	if(NULL == S->solution){
		error = Simulation_InitSolution(S);
		if(0 != error){
			S4_TRACE("< Simulation_GetLayerSolution (failed; Simulation_InitSolution returned %d) [omega=%f]\n", error, S->omega[0]);
			return error;
		}
	}
	sol = S->solution;

	for(int i = 0; i < S->n_layers; ++i){
		S4_Layer *L = &(S->layer[i]);
		if(L == layer){
			if((NULL == L->modes && L->copy < 0) || (L->copy >= 0 && NULL == S->layer[L->copy].modes) || !sol->solved[i]){
				error = Simulation_ComputeLayerSolution(S, L, layer_modes, layer_solution);
				if(0 != error){ // should never happen
					S4_TRACE("< Simulation_GetLayerSolution (failed; Simulation_ComputeLayerSolution returned %d) [omega=%f]\n", error, S->omega[0]);
					return 2;
				}
				if(L->copy < 0){
					L->modes = *layer_modes;
				}
			}else{
				if(L->copy >= 0){
					*layer_modes = S->layer[L->copy].modes;
				}else{
					*layer_modes = L->modes;
				}
				*layer_solution = &sol->ab[i*4*S->n_G];
			}

			S4_TRACE("< Simulation_GetLayerSolution [omega=%f]\n", S->omega[0]);
			return 0;
		}
	}

	S4_TRACE("< Simulation_GetLayerSolution (failed; layer not found) [omega=%f]\n", S->omega[0]);
	return -2;
}

int S4_Simulation_SolveLayer(S4_Simulation *S, S4_LayerID layer){
	S4_TRACE("> Simulation_SolveLayer(S=%p, layer=%d) [omega=%f]\n", S, layer, S->omega[0]);
	if(NULL == S){
		S4_TRACE("< Simulation_GetLayerSolution (failed; S == NULL) [omega=%f]\n", S->omega[0]);
		return -1;
	}
	if(layer < 0 || layer >= S->n_layers){
		S4_TRACE("< Simulation_GetLayerSolution (failed; layer == NULL) [omega=%f]\n", S->omega[0]);
		return -2;
	}

	LayerModes *Lmodes = NULL;
	std::complex<double> *Lsoln = NULL;
	int ret = Simulation_GetLayerSolution(S, &S->layer[layer], &Lmodes, &Lsoln);

	S4_TRACE("< Simulation_SolveLayer [omega=%f]\n", S->omega[0]);
	return ret;
}

int Simulation_ComputeLayerSolution(S4_Simulation *S, S4_Layer *L, LayerModes **layer_modes, std::complex<double> **layer_solution){
	S4_TRACE("> Simulation_ComputeLayerSolution(S=%p, L=%p (%s), layer_modes=%p (%p), LayerSolution=%p (%p)) [omega=%f]\n",
		S, L, (NULL != L && NULL != L->name ? L->name : ""), layer_modes, (NULL != layer_modes ? *layer_modes : NULL), layer_solution, (NULL != layer_solution ? *layer_solution : NULL), S->omega[0]);
	if(NULL != layer_modes && NULL == layer_solution){
		// only need to compute modes for L
		Simulation_ComputeLayerModes(S, L, layer_modes);
		S4_TRACE("< Simulation_ComputeLayerSolution [omega=%f]\n", S->omega[0]);
		return 0;
	}
	if(NULL == layer_solution){
		S4_TRACE("< Simulation_ComputeLayerSolution [omega=%f]\n", S->omega[0]);
		return 0;
	}
	// At this point, layer_solution != NULL, we need to get all modes

	S4_VERB(1, "Computing solution in layer: %s\n", NULL != L->name ? L->name : "");

	const size_t n = S->n_G;
	const size_t n2 = 2*n;
	const size_t n4 = 2*n2;

	// compute all modes and then get solution
	Solution_ *sol = S->solution;
	int which_layer = 0;
	bool found_layer = false;
	for(int i = 0; i < S->n_layers; ++i){
		S4_Layer *SL = &(S->layer[i]);
		if(NULL == SL->modes && SL->copy < 0){
			Simulation_ComputeLayerModes(S, SL, &SL->modes);
		}
		if(L == SL){
			found_layer = true;
			which_layer = i;
			if(SL->copy >= 0){
				*layer_modes = S->layer[SL->copy].modes;
			}else{
				*layer_modes = SL->modes;
			}
			*layer_solution = &sol->ab[i*n4];
		}
	}
	if(!found_layer){
		S4_TRACE("< Simulation_ComputeLayerSolution (failed; could not find layer) [omega=%f]\n", S->omega[0]);
		return -1;
	}

	// Make arrays of q, kp, and phi
	double *lthick = (double*)S4_malloc(sizeof(double)*S->n_layers);
	int *lepstype = (int*)S4_malloc(sizeof(int)*S->n_layers);
	const std::complex<double> **lq   = (const std::complex<double> **)S4_malloc(sizeof(const std::complex<double> *)*S->n_layers*4);
	if(NULL == lthick || NULL == lq){
		S4_TRACE("< Simulation_ComputeLayerSolution (failed; could not allocate work arrays) [omega=%f]\n", S->omega[0]);
		return 1;
	}
	const std::complex<double> **lepsinv  = lq  + S->n_layers;
	const std::complex<double> **lkp  = lepsinv  + S->n_layers;
	const std::complex<double> **lphi = lkp + S->n_layers;

	for(int i = 0; i < S->n_layers; ++i){
		const S4_Layer *SL = &(S->layer[i]);
		const S4_Layer *SLmodes = SL;
		if(SL->copy >= 0){ SLmodes = &S->layer[SL->copy]; }
		lthick[i] = SL->thickness;
		lq  [i] = SLmodes->modes->q;
		lepsinv [i] = SLmodes->modes->Epsilon_inv;
		lepstype[i] = SLmodes->modes->epstype;
		lkp [i] = SLmodes->modes->kp;
		lphi[i] = SLmodes->modes->phi;
	}

	// Compose the RCWA solution
	int error = 0;
	if(0 == S->exc.type){
		// Front incidence by planewave
		const size_t order = S->exc.sub.planewave.order;
		const bool inc_back = (0 != S->exc.sub.planewave.backwards);
		const size_t ind_fb = (inc_back ? S->n_layers-1 : 0);
		const size_t phicopy_size = (NULL == lphi[ind_fb] ? 0 : n2*n2);
		std::complex<double> *ab0 = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>)*(n2+phicopy_size));
		std::complex<double> *phicopy = (NULL == lphi[ind_fb] ? NULL : ab0 + n2);
		RNP::TBLAS::Fill(n2, 0., ab0,1);
		if(order < n){
			ab0[order+0] = std::complex<double>(S->exc.sub.planewave.hx[0], S->exc.sub.planewave.hx[1]);
			ab0[order+n] = std::complex<double>(S->exc.sub.planewave.hy[0], S->exc.sub.planewave.hy[1]);
		}
		// [ kp.phi.inv(q) -kp.phi.inv(q) ] [ a ] = [-ey;ex ]
		// [     phi            phi       ] [ b ]   [ hx;hy ]
		// We assume b = 0 or a = 0.
		// ab = inv(phi)*[ hx;hy ]
		if(NULL != lphi[ind_fb]){
			RNP::TBLAS::CopyMatrix<'A'>(n2,n2, lphi[ind_fb],n2, phicopy,n2);
			RNP::LinearSolve<'N'>(n2,1, phicopy,n2, ab0,n2, NULL, NULL);
		}

		if(S->options.use_less_memory){
			S4_TRACE("I  Calling SolveInterior(layer_count=%d, which_layer=%d, n=%d, lthick,lq,lkp,lphi={\n", S->n_layers, which_layer, S->n_G);
			for(int i = 0; i < S->n_layers; ++i){
				S4_TRACE("I    %f, %p (0,0=%f,%f), %p (0,0=%f,%f), %p (0,0=%f,%f)\n", lthick[i],
					lq[i], lq[i][0].real(), lq[i][0].imag(),
					lkp[i], NULL != lkp[i] ? lkp[i][0].real() : 0., NULL != lkp[i] ? lkp[i][0].imag() : 0.,
					lphi[i], NULL != lphi[0] ? lphi[i][0].real() : 1., NULL != lphi[0] ? lphi[i][0].imag() : 0.);
			}
			S4_TRACE("I   }, a0[0]=%f,%f, a0[n]=%f,%f, ...) [omega=%f]\n", ab0[0].real(), ab0[0].imag(), ab0[S->n_G].real(), ab0[S->n_G].imag(), S->omega[0]);

			error = SolveInterior(
				S->n_layers, which_layer,
				S->n_G,
				S->kx, S->ky,
				std::complex<double>(S->omega[0], S->omega[1]),
				lthick, lq, lepsinv, lepstype, lkp, lphi,
				inc_back ? NULL : ab0, // length 2*n
				inc_back ? ab0 : NULL, // bN
				(*layer_solution));
		}else{
			// Solve all at once
			std::complex<double> *pab = sol->ab;
			memset(pab, 0, sizeof(std::complex<double>) * S->n_layers * n4);
			if(!inc_back){
				memcpy(pab, ab0, sizeof(std::complex<double>) * n2);
			}else{
				memcpy(&pab[S->n_layers*n4 - n2], ab0, sizeof(std::complex<double>) * n2);
			}
			const size_t lwork = 6*S->n_layers*n2*n2;
			std::complex<double> *work = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>) * lwork);
			size_t *iwork = (size_t*)S4_malloc(sizeof(size_t) * S->n_layers*n2);
			SolveAll(
				S->n_layers, S->n_G, S->kx, S->ky,
				std::complex<double>(S->omega[0], S->omega[1]),
				lthick, lq, lepsinv, lepstype, lkp, lphi,
				pab,
				work, iwork, lwork
			);
			S4_free(iwork);
			S4_free(work);
			for(size_t i = 0; i < S->n_layers; ++i){
				sol->solved[i] = 1;
			}
		}
		S4_free(ab0);
	}else if(2 == S->exc.type){
		Excitation_Exterior *ext = &(S->exc.sub.exterior);
		// Front incidence by planewave
		size_t phicopy_size = (NULL == lphi[0] ? 0 : n2*n2);
		std::complex<double> *a0 = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>)*(n4+phicopy_size));
		std::complex<double> *bN = a0 + n2;
		RNP::TBLAS::Fill(n4, 0., a0,1);
		for(size_t i = 0; i < ext->n; ++i){
			int gindex = ext->Gindex1[2*i+0];
			const int pol = ext->Gindex1[2*i+1];
			if(gindex > 0){
				gindex--;
				if(0 == pol){ // pol is for E field
					a0[gindex+n] =  std::complex<double>(ext->coeff[2*i+0], ext->coeff[2*i+1]);
				}else{
					a0[gindex+0] = -std::complex<double>(ext->coeff[2*i+0], ext->coeff[2*i+1]);
				}
			}else if(gindex < 0){
				gindex++;
				gindex = -gindex;
				if(0 == pol){ // pol is for E field
					bN[gindex+n] = -std::complex<double>(ext->coeff[2*i+0], ext->coeff[2*i+1]);
				}else{
					bN[gindex+0] =  std::complex<double>(ext->coeff[2*i+0], ext->coeff[2*i+1]);
				}
			}
		}

		S4_TRACE("I  Calling SolveInterior(layer_count=%d, which_layer=%d, n=%d, lthick,lq,lkp,lphi={\n", S->n_layers, which_layer, S->n_G);
		for(int i = 0; i < S->n_layers; ++i){
			S4_TRACE("I    %f, %p (0,0=%f,%f), %p (0,0=%f,%f), %p (0,0=%f,%f)\n", lthick[i],
				lq[i], lq[i][0].real(), lq[i][0].imag(),
				lkp[i], NULL != lkp[i] ? lkp[i][0].real() : 0., NULL != lkp[i] ? lkp[i][0].imag() : 0.,
				lphi[i], NULL != lphi[0] ? lphi[i][0].real() : 1., NULL != lphi[0] ? lphi[i][0].imag() : 0.);
		}
		S4_TRACE("I   }, a0[0]=%f,%f, a0[n]=%f,%f, ...) [omega=%f]\n", a0[0].real(), a0[0].imag(), a0[S->n_G].real(), a0[S->n_G].imag(), S->omega[0]);

		error = SolveInterior(
			S->n_layers, which_layer,
			S->n_G,
			S->kx, S->ky,
			std::complex<double>(S->omega[0], S->omega[1]),
			lthick, lq, lepsinv, lepstype, lkp, lphi,
			a0, // length 2*n
			bN, // bN
			(*layer_solution));
		S4_free(a0);
	}else if(1 == S->exc.type){
		S4_Layer *l[2];
		l[0] = S->exc.layer;
		const int li = (l[0] - &S->layer[0]);
		if(NULL == l[0]){ return 13; }
		l[1] = l[0]+1;
		std::complex<double> *ab = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>)*(n4+n4*n4+n2*n2));
		std::complex<double> *work4 = ab + n4;
		std::complex<double> *work2 = work4 + n4*n4;
		std::complex<double> J0[3] = {
			std::complex<double>(S->exc.sub.dipole.moment[0],S->exc.sub.dipole.moment[1]),
			std::complex<double>(S->exc.sub.dipole.moment[2],S->exc.sub.dipole.moment[3]),
			std::complex<double>(S->exc.sub.dipole.moment[4],S->exc.sub.dipole.moment[5])
		};
		for(size_t i = 0; i < n; ++i){
			const double phaseangle = -(S->kx[i] * S->exc.sub.dipole.pos[0] + S->ky[i] * S->exc.sub.dipole.pos[1]);
			const std::complex<double> phase(cos(phaseangle), sin(phaseangle));
			ab[0*n+i] = J0[2]*phase; // -ky eta jz
			//ab[1*n+i] = ; //  kx eta jz
			ab[2*n+i] =  J0[1]*phase; //   jy
			ab[3*n+i] = -J0[0]*phase; //  -jx
		}
		// Make eta*jz
		RNP::TBLAS::MultMV<'N'>(n,n, 1.,l[0]->modes->Epsilon_inv,n, &ab[0*n],1, 0.,&ab[1*n],1);
		RNP::TBLAS::Copy(n, &ab[1*n],1, &ab[0*n],1);
		// finish ab[0*n] and ab[1*n]
		for(size_t i = 0; i < n; ++i){
			ab[0*n+i] *= -S->ky[i];
			ab[1*n+i] *=  S->kx[i];
		}
		// At this point ab is (p_z, p_par)
		// Solve:
		// [ (omega^2 - Kappa_{l+1}) Phi_{l+1} q_{l+1}^{-1} [ 1 - f_{l+1} S_{21}(l+1,N) ,
		//              (omega^2 - Kappa_l) Phi_l q_l^{-1} [ 1 - f_l S_{12}(0,l) ;
		//   Phi_{l+1} [ 1 + f_{l+1} S_{21}(l+1,N) , -Phi_l [ 1 + f_l S_{12}(0,l)         ] [ a_{l+1} ; b_l ]
		// == [ p_z ; p_par ]
		//
		// In code:
		// [factorization] [RHS] == [pzp]
		//

		// First solve for the outgoing waves immediately adjacent to the dipole
		// Get S matrix portions first
		Simulation_GetSMatrix(S, 0, li, work4);
		RNP::TBLAS::CopyMatrix<'A'>(n2,n2, &work4[0+n2*n4],n4, work2,n2);
		Simulation_GetSMatrix(S, li+1, -1, work4);
		// Make upper right
		for(size_t i = 0; i < n2; ++i){ // first scale work2 to make -f_l*S12(0,l)
			std::complex<double> f = -std::exp(lq[li][i] * std::complex<double>(0,lthick[li]));
			RNP::TBLAS::Scale(n2, f, &work2[i+0*n2],n2);
		}
		RNP::TBLAS::CopyMatrix<'A'>(n2,n2, work2,n2, &work4[0+n2*n4],n4);
		for(size_t i = 0; i < n2; ++i){
			work4[i+(n2+i)*n4] += 1.;
			RNP::TBLAS::Scale(n2, 1./lq[li][i], &work4[i+n2*n4],n4);
		}
		if(NULL != lphi[li]){ // use lower right as buffer
			RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,lphi[li],n2, &work4[0+n2*n4],n4, 0.,&work4[n2+n2*n4],n4);
		}else{
			RNP::TBLAS::CopyMatrix<'A'>(n2,n2, &work4[0+n2*n4],n4, &work4[n2+n2*n4],n4);
		}
		MultKPMatrix("N",
			std::complex<double>(S->omega[0], S->omega[1]),
			S->n_G,
			S->kx, S->ky,
			lepsinv[li], lepstype[li], lkp[li],
			n2, &work4[n2+n2*n4],n4,
			&work4[0+n2*n4],n4
		); // upper right complete

		// For lower right, -f_l S12(0,li) is still in work2
		for(size_t i = 0; i < n2; ++i){
			work2[i+i*n4] -= 1.;
		}
		if(NULL != lphi[li]){ // use lower right as buffer
			RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,lphi[li],n2, work2,n2, 0.,&work4[n2+n2*n4],n4);
		}else{
			RNP::TBLAS::CopyMatrix<'A'>(n2,n2, work2,n2, &work4[n2+n2*n4],n4);
		} // lower right complete

		// For upper left, make f_l+1 S21(l+1,N) in lower left first
		for(size_t i = 0; i < n2; ++i){ // make f_l+1*S12(l+1,N)
			std::complex<double> f = std::exp(lq[li+1][i] * std::complex<double>(0,lthick[li+1]));
			RNP::TBLAS::Scale(n2, f, &work4[n2+i+0*n2],n4);
		}
		RNP::TBLAS::CopyMatrix<'A'>(n2,n2, &work4[n2+0*n2],n4, &work4[0+0*n4],n4); // copy to upper left
		for(size_t i = 0; i < n2; ++i){
			work4[0+i+(0+i)*n4] -= 1.;
		}
		for(size_t i = 0; i < n2; ++i){
			RNP::TBLAS::Scale(n2, -1./lq[li+1][i], &work4[(0+i)+0*n4],n4);
		}
		if(NULL != lphi[li]){ // use work2 as buffer
			RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,lphi[li+1],n2, &work4[0+0*n4],n4, 0.,work2,n2);
		}else{
			RNP::TBLAS::CopyMatrix<'A'>(n2,n2, &work4[0+0*n4],n4, work2,n2);
		}
		MultKPMatrix("N",
			std::complex<double>(S->omega[0], S->omega[1]),
			S->n_G,
			S->kx, S->ky,
			lepsinv[li], lepstype[li], lkp[li],
			n2, work2,n2,
			&work4[0+0*n4],n4
		); // upper left complete

		// Lower left
		for(size_t i = 0; i < n2; ++i){
			work4[n2+i+(0+i)*n4] += 1.;
		}
		if(NULL != lphi[li+1]){
			RNP::TBLAS::CopyMatrix<'A'>(n2,n2, &work4[n2+0*n4],n4, work2,n2);
			RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,lphi[li+1],n2, work2,n2, 0.,&work4[n2+0*n4],n4);
		}

		// Solve for a_l+1, bl
		RNP::LinearSolve<'N'>(n4,1, work4,n4, ab,n4);

		// Now solve for the waves in the layers that have been requested
		if(which_layer <= li){
			error = SolveInterior(
				li, which_layer,
				S->n_G,
				S->kx, S->ky,
				std::complex<double>(S->omega[0], S->omega[1]),
				lthick, lq, lepsinv, lepstype, lkp, lphi,
				NULL, // length 2*n
				&ab[n2], // bN
				(*layer_solution));
		}else{
			error = SolveInterior(
				S->n_layers-li, which_layer-li,
				S->n_G,
				S->kx, S->ky,
				std::complex<double>(S->omega[0], S->omega[1]),
				lthick+li, lq+li, lepsinv+li, lepstype+li, lkp+li, lphi+li,
				&ab[0], // length 2*n
				NULL, // bN
				(*layer_solution));
		}
		S4_free(ab);
	}
	sol->solved[which_layer] = 1;

	S4_TRACE("I  ab[0] = %f,%f [omega=%f]\n", (*layer_solution)[0].real(), (*layer_solution)[0].imag(), S->omega[0]);

	S4_free(lq);
	S4_free(lepstype);
	S4_free(lthick);

	if(0 != error){
		S4_TRACE("< Simulation_ComputeLayerSolution (failed; SolveInterior returned %d) [omega=%f]\n", error, S->omega[0]);
		return error;
	}

	S4_TRACE("< Simulation_ComputeLayerSolution [omega=%f]\n", S->omega[0]);
	return 0;
}

double Simulation_GetUnitCellSize(const S4_Simulation *S){
	if(0 == S->Lr[2] && 0 == S->Lr[3]){
		return hypot(S->Lr[0], S->Lr[1]);
	}else{
		return fabs(S->Lr[0]*S->Lr[3] - S->Lr[1]*S->Lr[2]);
	}
}

int Simulation_ComputeLayerModes(S4_Simulation *S, S4_Layer *L, LayerModes **layer_modes){
	S4_TRACE("> Simulation_ComputeLayerModes(S=%p, L=%p (%s), modes=%p (%p)) [omega=%f]\n", S, L, (NULL != L && NULL != L->name ? L->name : ""), layer_modes, (NULL != layer_modes ? *layer_modes : NULL), S->omega[0]);
	if(NULL == S){
		S4_TRACE("< Simulation_ComputeLayerModes (failed; S == NULL) [omega=%f]\n", S->omega[0]);
		return -1;
	}
	if(NULL == L){
		S4_TRACE("< Simulation_ComputeLayerModes (failed; L == NULL) [omega=%f]\n", S->omega[0]);
		return -2;
	}
	if(NULL == layer_modes){
		S4_TRACE("< Simulation_ComputeLayerModes (early exit; modes == NULL) [omega=%f]\n", S->omega[0]);
		return 0;
	}

	S4_VERB(1, "Computing modes of layer: %s\n", NULL != L->name ? L->name : "");

	*layer_modes = (LayerModes*)malloc(sizeof(LayerModes));
	LayerModes *pB = *layer_modes;
	const int n = S->n_G;
	const int n2 = 2*n;
	const int nn = n*n;
	const int n2n2 = n2*n2;

	size_t kp_size = n2n2;
	size_t phi_size = n2n2;
	size_t Epsilon_inv_size = nn;
	size_t Epsilon2_size = n2n2;
	pB->epstype = EPSILON2_TYPE_FULL;
	if(0 == L->pattern.nshapes){
		const S4_Material *M;
		if(L->copy < 0){
			M = &S->material[L->material];
		}else{
			M = &S->material[S->layer[L->copy].material];
		}
		if(0 == M->type){
			phi_size = 0;
			pB->epstype = EPSILON2_TYPE_BLKDIAG1_SCALAR;
		}
	}
	if(S->options.use_less_memory){
		kp_size = 0;
	}

	pB->q = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>)*(
		2*n + // for q
		kp_size + phi_size + Epsilon_inv_size + Epsilon2_size
		));
	pB->kp = pB->q + n2;
	pB->phi = pB->kp + kp_size;
	pB->Epsilon_inv = pB->phi + phi_size;
	pB->Epsilon2 = pB->Epsilon_inv + Epsilon_inv_size;
	//pB->Epsilon = pB->Epsilon2 + n2n2;

	if(0 == phi_size){
		pB->phi = NULL;
	}
	if(0 == kp_size){
		pB->kp = NULL;
	}

	// Outline of the epsilon matrix generation code below:
	//
	// If no shapes
	//   If background material is scalar epsilon
	//     Solve uniform layer eigensystem
	//   Else
	//     Generate the 4 quadrants of Epsilon2 and solve the generic layer eigensystem
	// Else -- not uniform layer
	//   If discretize epsilon
	//     If using subpixel smoothing
	//       Apply Kottke
	//     ElseIf polarization basis
	//       If use complex basis
	//         Apply PolBasisJones
	//       Else
	//         Apply PolBasisNV
	//     Else
	//       Apply FFT
	//   Else
	//     Apply ClosedForm
	if(0 == L->pattern.nshapes){
		const S4_Material *M;
		if(L->copy < 0){
			//eps_scalar = Simulation_GetEpsilonByName(S, L->material);
			M = &S->material[L->material];
		}else{
			//eps_scalar = Simulation_GetEpsilonByName(S, Simulation_GetLayerByName(S, L->copy, NULL)->material);
			M = &S->material[S->layer[L->copy].material];
		}
		if(0 == M->type){
			std::complex<double> eps_scalar(M->eps.s[0], M->eps.s[1]);
			RNP::TBLAS::SetMatrix<'A'>(n,n,0.,1./eps_scalar,pB->Epsilon_inv, n);
			RNP::TBLAS::SetMatrix<'A'>(n2,n2,0.,eps_scalar,pB->Epsilon2, n2);
			SolveLayerEigensystem_uniform(
				std::complex<double>(S->omega[0],S->omega[1]), n, S->kx, S->ky,
				eps_scalar, pB->q, pB->kp, pB->phi);
		}else{
			RNP::TBLAS::SetMatrix<'A'>(n,n,0.,1./std::complex<double>(M->eps.abcde[8],M->eps.abcde[9]),pB->Epsilon_inv, n);
			RNP::TBLAS::SetMatrix<'A'>(n,n,0.,std::complex<double>(M->eps.abcde[0],M->eps.abcde[1]),&pB->Epsilon2[0+0*n2], n2);
			RNP::TBLAS::SetMatrix<'A'>(n,n,0.,std::complex<double>(M->eps.abcde[4],M->eps.abcde[5]),&pB->Epsilon2[n+0*n2], n2);
			RNP::TBLAS::SetMatrix<'A'>(n,n,0.,std::complex<double>(M->eps.abcde[2],M->eps.abcde[3]),&pB->Epsilon2[0+n*n2], n2);
			RNP::TBLAS::SetMatrix<'A'>(n,n,0.,std::complex<double>(M->eps.abcde[6],M->eps.abcde[7]),&pB->Epsilon2[n+n*n2], n2);

			S4_VERB(1, "Solving eigensystem of layer: %s\n", NULL != L->name ? L->name : "");
			SolveLayerEigensystem(
				std::complex<double>(S->omega[0],S->omega[1]), n, S->kx, S->ky,
				pB->Epsilon_inv, pB->Epsilon2, pB->epstype, pB->q, pB->kp, pB->phi);
		}
	}else{ // not a uniform layer
		S4_VERB(1, "Generating epsilon matrix of layer: %s\n", NULL != L->name ? L->name : "");
		if(S->options.use_experimental_fmm){
			FMMGetEpsilon_Experimental(S, L, n, pB->Epsilon2, pB->Epsilon_inv);
		}else{
			if(S->options.use_discretized_epsilon){
				if(S->options.use_subpixel_smoothing){
					FMMGetEpsilon_Kottke(S, L, n, pB->Epsilon2, pB->Epsilon_inv);
				}else{ // not using subpixel smoothing
					FMMGetEpsilon_FFT(S, L, n, pB->Epsilon2, pB->Epsilon_inv);
					if(S->options.use_polarization_basis){
						if(S->options.use_jones_vector_basis){
							FMMGetEpsilon_PolBasisJones(S, L, n, pB->Epsilon2, pB->Epsilon_inv);
						}else if(S->options.use_normal_vector_basis){
							FMMGetEpsilon_PolBasisNV(S, L, n, pB->Epsilon2, pB->Epsilon_inv);
						}else{
							FMMGetEpsilon_PolBasisVL(S, L, n, pB->Epsilon2, pB->Epsilon_inv);
						}
					}
				}
			}else{
				FMMGetEpsilon_ClosedForm(S, L, n, pB->Epsilon2, pB->Epsilon_inv);
				if(S->options.use_polarization_basis){
					if(S->options.use_jones_vector_basis){
						FMMGetEpsilon_PolBasisJones(S, L, n, pB->Epsilon2, pB->Epsilon_inv);
					}else if(S->options.use_normal_vector_basis){
						FMMGetEpsilon_PolBasisNV(S, L, n, pB->Epsilon2, pB->Epsilon_inv);
					}else{
						FMMGetEpsilon_PolBasisVL(S, L, n, pB->Epsilon2, pB->Epsilon_inv);
					}
				}
			}
		}
//std::cerr << pB->Epsilon2[0] << "\t" << pB->Epsilon2[1] << "\t" << pB->Epsilon_inv[0] << "\t" << pB->Epsilon_inv[1] << std::endl;
		S4_VERB(1, "Solving eigensystem of layer: %s\n", NULL != L->name ? L->name : "");
		{
			size_t lwork = (size_t)-1;
			double *rwork = (double*)S4_malloc(sizeof(double) * 4*n);
			std::complex<double> *work = NULL;
			std::complex<double> dum;
			SolveLayerEigensystem(
				std::complex<double>(S->omega[0],S->omega[1]), n,
				S->kx, S->ky,
				pB->Epsilon_inv, pB->Epsilon2, pB->epstype,
				pB->q, pB->kp, pB->phi,
				&dum, rwork, lwork
			);
			lwork = (int)(dum.real() + 0.5);
			work = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>) * lwork);
			SolveLayerEigensystem(
				std::complex<double>(S->omega[0],S->omega[1]), n,
				S->kx, S->ky,
				pB->Epsilon_inv, pB->Epsilon2, pB->epstype,
				pB->q, pB->kp, pB->phi,
				work, rwork, lwork
			);
			S4_free(work);
			S4_free(rwork);
		}
	}
	S4_TRACE("I  q[0] = %f,%f [omega=%f]\n", pB->q[0].real(), pB->q[0].imag(), S->omega[0]);

	S4_TRACE("< Simulation_ComputeLayerModes [omega=%f]\n", S->omega[0]);
	return 0;
}

// realistically, can only return an InitSolution code
int S4_Simulation_GetPowerFlux(S4_Simulation *S, S4_LayerID id, const double *offset, double *powers){
	S4_TRACE("> S4_Simulation_GetPowerFlux(S=%p, layer=%d, offset=%f, powers=%p) [omega=%f]\n",
		S, id, (NULL == offset ? 0 : *offset), powers, S->omega[0]);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(id < 0 || id >= S->n_layers){ return -2; }
	if(NULL == powers){ ret = -4; }
	if(0 != ret){
		S4_TRACE("< S4_Simulation_GetPowerFlux (failed; ret = %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}
	S4_Layer *layer = &S->layer[id];

	LayerModes *Lmodes;
	std::complex<double> *Lsoln;
	const double off = (NULL != offset ? *offset : 0);

	ret = Simulation_GetLayerSolution(S, layer, &Lmodes, &Lsoln);
	if(0 != ret){
		S4_TRACE("< S4_Simulation_GetPowerFlux (failed; Simulation_GetLayerSolution returned %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	const int n = S->n_G;
	const int n2 = 2*n;
	const int n4 = 2*n2;

	std::complex<double> *ab = (std::complex<double> *)S4_malloc(sizeof(std::complex<double>) * (n4+4*n2));
	if(NULL == ab){
		S4_TRACE("< S4_Simulation_GetPowerFlux (failed; allocation failed) [omega=%f]\n", S->omega[0]);
		return 1;
	}
	std::complex<double> *work = ab + n4;

	memcpy(ab, Lsoln, sizeof(std::complex<double>) * n4);
	TranslateAmplitudes(n, Lmodes->q, layer->thickness, off, ab);

	std::complex<double> forw, back;
	GetZPoyntingFlux(n, S->kx, S->ky, std::complex<double>(S->omega[0],S->omega[1]), Lmodes->q, Lmodes->Epsilon_inv, Lmodes->epstype, Lmodes->kp, Lmodes->phi, ab, &forw, &back, work);
	powers[0] = forw.real();
	powers[1] = back.real();
	powers[2] = forw.imag();
	powers[3] = back.imag();

	S4_free(ab);
	S4_TRACE("< S4_Simulation_GetPowerFlux returning %f, %f, %f, %f [omega=%f]\n", powers[0], powers[1], powers[2], powers[3], S->omega[0]);
	return 0;
}

int Simulation_GetPoyntingFluxByG(S4_Simulation *S, S4_Layer *layer, double offset, double *powers){
	S4_TRACE("> Simulation_GetPoyntingFluxByG(S=%p, layer=%p, offset=%f, powers=%p) [omega=%f]\n",
		S, layer, offset, powers, S->omega[0]);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(NULL == powers){ ret = -4; }
	if(0 != ret){
		S4_TRACE("< Simulation_GetPoyntingFluxByG (failed; ret = %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	LayerModes *Lmodes;
	std::complex<double> *Lsoln;

	ret = Simulation_GetLayerSolution(S, layer, &Lmodes, &Lsoln);
	if(0 != ret){
		S4_TRACE("< Simulation_GetPoyntingFluxByG (failed; Simulation_GetLayerSolution returned %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	const int n = S->n_G;
	const int n2 = 2*n;
	const int n4 = 2*n2;

	std::complex<double> *ab = (std::complex<double> *)S4_malloc(sizeof(std::complex<double>) * (n2+n4+4*n2));
	if(NULL == ab){
		S4_TRACE("< Simulation_GetPoyntingFluxByG (failed; allocation failed) [omega=%f]\n", S->omega[0]);
		return 1;
	}
	std::complex<double> *forw = ab + n4;
	std::complex<double> *back = forw + n;
	std::complex<double> *work = back + n;

	memcpy(ab, Lsoln, sizeof(std::complex<double>) * n4);
	TranslateAmplitudes(n, Lmodes->q, layer->thickness, offset, ab);

	GetZPoyntingFluxComponents(n, S->kx, S->ky, std::complex<double>(S->omega[0],S->omega[1]), Lmodes->q, Lmodes->Epsilon_inv, Lmodes->epstype, Lmodes->kp, Lmodes->phi, ab, forw, back, work);
	for(int i = 0; i < n; ++i){
		powers[4*i+0] = forw[i].real();
		powers[4*i+1] = back[i].real();
		powers[4*i+2] = forw[i].imag();
		powers[4*i+3] = back[i].imag();
	}

	S4_free(ab);
	S4_TRACE("< Simulation_GetPoyntingFluxByG [omega=%f]\n", S->omega[0]);
	return 0;
}

int Simulation_GetPropagationConstants(S4_Simulation *S, S4_Layer *L, double *q){
	S4_TRACE("> Simulation_GetPropagationConstants(S=%p, layer=%p, q=%p) [omega=%f]\n",
		S, L, q, S->omega[0]);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == L){ ret = -2; }
	if(0 != ret){
		S4_TRACE("< Simulation_GetPropagationConstants (failed; ret = %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	const int n = S->n_G;

	// compute all modes and then get solution
	LayerModes *layer_modes;
	bool found_layer = false;
	for(int i = 0; i < S->n_layers; ++i){
		S4_Layer *SL = &(S->layer[i]);
		if(NULL == SL->modes && SL->copy < 0){
			Simulation_ComputeLayerModes(S, SL, &SL->modes);
		}
		if(L == SL){
			found_layer = true;
			if(SL->copy >= 0){
				layer_modes = S->layer[SL->copy].modes;
			}else{
				layer_modes = SL->modes;
			}
		}
	}
	if(!found_layer){
		S4_TRACE("< Simulation_GetPropagationConstants (failed; could not find layer) [omega=%f]\n", S->omega[0]);
		return -1;
	}

	for(int i = 0; i < n; ++i){
		q[2*i+0] = layer_modes->q[i].real();
		q[2*i+1] = layer_modes->q[i].imag();
	}

	S4_TRACE("< Simulation_GetPropagationConstants [omega=%f]\n", S->omega[0]);
	return ret;
}

int Simulation_GetAmplitudes(S4_Simulation *S, S4_Layer *layer, double offset, double *forw, double *back){
	S4_TRACE("> Simulation_GetAmplitudes(S=%p, layer=%p, offset=%f, forw=%p, back=%p) [omega=%f]\n",
		S, layer, offset, forw, back, S->omega[0]);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(0 != ret){
		S4_TRACE("< Simulation_GetAmplitudes (failed; ret = %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	LayerModes *Lmodes;
	std::complex<double> *Lsoln;

	ret = Simulation_GetLayerSolution(S, layer, &Lmodes, &Lsoln);
	if(0 != ret){
		S4_TRACE("< Simulation_GetAmplitudes (failed; Simulation_GetLayerSolution returned %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	const int n = S->n_G;
	const int n2 = 2*n;
	const int n4 = 2*n2;

	std::complex<double> *ab = (std::complex<double> *)S4_malloc(sizeof(std::complex<double>) * n4);
	if(NULL == ab){
		S4_TRACE("< Simulation_GetAmplitudes (failed; allocation failed) [omega=%f]\n", S->omega[0]);
		return 1;
	}

	memcpy(ab, Lsoln, sizeof(std::complex<double>) * n4);
	TranslateAmplitudes(n, Lmodes->q, layer->thickness, offset, ab);

	if(NULL != forw){
		for(int i = 0; i < n2; ++i){
			forw[2*i+0] = ab[i].real();
			forw[2*i+1] = ab[i].imag();
		}
	}
	if(NULL != back){
		for(int i = 0; i < n2; ++i){
			back[2*i+0] = ab[n2+i].real();
			back[2*i+1] = ab[n2+i].imag();
		}
	}

	S4_free(ab);
	S4_TRACE("< Simulation_GetAmplitudes [omega=%f]\n", S->omega[0]);
	return 0;
}


int S4_Simulation_GetWaves(S4_Simulation *S, S4_LayerID id, S4_real *wave){
	S4_TRACE("> Simulation_GetWaves(S=%p, layer=%d, wave=%p) [omega=%f]\n",
		S, id, wave, S->omega[0]);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(id < 0 || id >= S->n_layers){ ret = -2; }
	if(NULL == wave){ ret = -3; }
	if(0 != ret){
		S4_TRACE("< Simulation_GetWaves (failed; ret = %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	S4_Layer *layer = &S->layer[id];

	const double w = S->omega[0];

	LayerModes *Lmodes;
	std::complex<double> *Lsoln;
	ret = Simulation_GetLayerSolution(S, layer, &Lmodes, &Lsoln);
	if(0 != ret){
		S4_TRACE("< Simulation_GetWaves (failed; Simulation_GetLayerSolution returned %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	const int n = S->n_G;
	const int n2 = 2*n;
	const int n4 = 2*n2;

	const std::complex<double> *ab = Lsoln;

	for(int i = 0; i < n; ++i){
		const double kx = S->kx[i];
		const double ky = S->ky[i];
		const std::complex<double> qi = Lmodes->q[i];
		for(int j = 0; j < 2; ++j){
			double *curwave = &wave[(j*n+i)*11];
			curwave[0] = kx;
			curwave[1] = ky;
			curwave[2] = Lmodes->q[i].real();
			curwave[3] = Lmodes->q[i].imag();
			if(0 != j){
				curwave[2] = -curwave[2];
				curwave[3] = -curwave[3];
			}
			double rd[3] = {
				curwave[0],
				curwave[1],
				curwave[2]
			};
			geom_normalize3d(rd);
/*
if(qi.imag() == 0){
	std::cout << "i = " << i << ", j = " << j << ", ab = " << ab[i+j*n2] << ", " << ab[i+n+j*n2] << std::endl;
	std::cout << "   kx = " << kx << ", ky = " << ky << std::endl;
	std::cout << "   qi = " << qi << std::endl;
}*/

{
			const double qsign(0 == j ? 1 : -1);
			const std::complex<double> Hx = ab[i  +j*n2];
			const std::complex<double> Hy = ab[i+n+j*n2];
			const std::complex<double> H[2][3] = {
				{ Hx, 0, -(kx*Hx)/(qsign*qi.real()) },
				{ 0, Hy, -(ky*Hy)/(qsign*qi.real()) }
			};
			// The E field is Cross[H,{kx,ky,(qsign*qi.real())}] / leps
			std::complex<double> E[2][3];
			const double sgn(qsign);
			const S4_Material *M = NULL;
			if(layer->material >= 0){
				M = &S->material[layer->material];
			}else{
				M = &S->material[S->layer[layer->copy].material];
			}
			const double leps = (NULL != M ? M->eps.s[0] : 1);
			E[0][2] = ky*Hx / leps;
			E[0][0] =  kx*E[0][2]       / (sgn*qi.real());
			E[0][1] = (ky*E[0][2] - Hx*w) / (sgn*qi.real());
			E[1][2] = (-kx*Hy) / leps;
			E[1][0] = (Hy*w + kx*E[1][2]) / (sgn*qi.real());
			E[1][1] = (ky*E[1][2]) / (sgn*qi.real());
/*if(qi.imag() == 0){
	std::cout << "   leps = " << leps << std::endl;
	std::cout << "   E[0] = " << E[0][0] << ", " << E[0][1] << ", " << E[0][2] << std::endl;
	std::cout << "   E[1] = " << E[1][0] << ", " << E[1][1] << ", " << E[1][2] << std::endl;
}*/
			double ru[3] = { 1, 0, -kx/(qsign*Lmodes->q[i].real()) };
			geom_normalize3d(ru);
			double v[3];
			geom_cross3d(rd, ru, v);
			std::complex<double> c0(ru[0]*(E[0][0]+E[1][0]) + ru[1]*(E[0][1]+E[1][1]) + ru[2]*(E[0][2]+E[1][2]));
			std::complex<double> c1( v[0]*(E[0][0]+E[1][0]) +  v[1]*(E[0][1]+E[1][1]) +  v[2]*(E[0][2]+E[1][2]));
			// Re-polarize so u is the s-polarization vector
			{
				// Don't touch d, only perform a rotation of u around d.
				double p[3] = {1,0,0};
				const double dxy2 = rd[0]*rd[0]+rd[1]*rd[1];
				const double tol = 1e-15;
				if(dxy2 > tol*tol){
					if(fabs(rd[2]) > tol){
						p[0] = rd[0];
						p[1] = rd[1];
						p[2] = -dxy2/rd[2];
					}else{
						p[0] = 0;
						p[1] = 0;
						p[2] = 1;
					}
				}
				// Orthogonalize p to be sure
				double pdd(geom_dot3d(p, rd));
				p[0] -= pdd*rd[0];
				p[1] -= pdd*rd[1];
				p[2] -= pdd*rd[2];
				geom_normalize3d(p);
				// Now perform rotations
				const double cs = geom_dot3d(p, ru);
				const double sn = geom_dot3d(p, v);
				const std::complex<double> c0p = cs*c0 + sn*c1;
				c1 = cs*c1 - sn*c0;
				c0 = c0p;

				curwave[ 4] = p[0];
				curwave[ 5] = p[1];
				curwave[ 6] = p[2];
				curwave[ 7] = c0.real();
				curwave[ 8] = c0.imag();
				curwave[ 9] = c1.real();
				curwave[10] = c1.imag();
			}
}
		}
	}

	S4_TRACE("< Simulation_GetWaves [omega=%f]\n", S->omega[0]);
	return 0;
}

void Simulation_DestroySolution(S4_Simulation *S){
	S4_TRACE("> Simulation_DestroySolution(S=%p) [omega=%f]\n", S, S->omega[0]);
	Solution_ *sol;
	if(NULL == S){
		S4_TRACE("< Simulation_DestroySolution (early exit; S == NULL) [omega=%f]\n", S->omega[0]);
		return;
	}
	sol = S->solution;
	if(NULL == sol){
		S4_TRACE("< Simulation_DestroySolution (early exit; sol == NULL) [omega=%f]\n", S->omega[0]);
		return;
	}
	if(NULL != sol->ab){
		S4_free(sol->ab);
		sol->ab = NULL;
	}
	if(NULL != sol->solved){
		free(sol->solved);
		sol->solved = NULL;
	}
	free(S->solution); S->solution = NULL;

	S4_TRACE("< Simulation_DestroySolution [omega=%f]\n", S->omega[0]);
}

void Simulation_DestroyLayerSolutions(S4_Simulation *S){
	S4_TRACE("> Simulation_DestroyLayerSolutions(S=%p) [omega=%f]\n", S, S->omega[0]);
	Solution_ *sol;
	if(NULL == S){
		S4_TRACE("< Simulation_DestroyLayerSolutions (early exit; S == NULL) [omega=%f]\n", S->omega[0]);
		return;
	}
	sol = S->solution;
	if(NULL == sol){
		S4_TRACE("< Simulation_DestroyLayerSolutions (early exit; sol == NULL) [omega=%f]\n", S->omega[0]);
		return;
	}

	for(int i = 0; i < S->n_layers; ++i){
		sol->solved[i] = 0;
	}

	S4_TRACE("< Simulation_DestroyLayerSolutions [omega=%f]\n", S->omega[0]);
}

S4_Material* Simulation_GetMaterialByName(const S4_Simulation *S, const char *name, int *index){
	S4_TRACE("> Simulation_GetMaterialByName(S=%p, name=%p (%s)) [omega=%f]\n",
		S, name, (NULL == name ? "" : name), S->omega[0]);
	for(int i = 0; i < S->n_materials; ++i){
		if(0 == strcmp(S->material[i].name, name)){
			if(NULL != index){ *index = i; }
			S4_TRACE("< Simulation_GetMaterialByName returning %s [omega=%f]\n", S->material[i].name, S->omega[0]);
			return &(S->material[i]);
		}
	}
	S4_TRACE("< Simulation_GetMaterialByName (failed; material name not found) [omega=%f]\n", S->omega[0]);
	return NULL;
}

S4_Material* Simulation_GetMaterialByIndex(const S4_Simulation *S, int i){
	S4_TRACE("> Simulation_GetMaterialByIndex(S=%p, i=%d) [omega=%f]\n", S, i, S->omega[0]);
	if(i < 0 || i >= S->n_materials){
		S4_TRACE("< Simulation_GetMaterialByIndex (failed; index out of bounds) [omega=%f]\n", S->omega[0]);
		return NULL;
	}
	S4_TRACE("< Simulation_GetMaterialByIndex returning %s [omega=%f]\n", S->material[i].name, S->omega[0]);
	return &(S->material[i]);
}
S4_Layer* Simulation_GetLayerByName(const S4_Simulation *S, const char *name, int *index){
	S4_TRACE("> Simulation_GetLayerByName(S=%p, name=%p (%s), index=%p) [omega=%f]\n",
		S, name, (NULL == name ? "" : name), index, S->omega[0]);
	if(NULL == S || NULL == name){ return NULL; }
	for(int i = 0; i < S->n_layers; ++i){
		if(0 == strcmp(S->layer[i].name, name)){
			if(NULL != index){ *index = i; }
			S4_TRACE("< Simulation_GetLayerByName returning %p [omega=%f]\n", &(S->layer[i]), S->omega[0]);
			return &(S->layer[i]);
		}
	}
	S4_TRACE("< Simulation_GetLayerByName (failed; name not found) [omega=%f]\n", S->omega[0]);
	return NULL;
}

static void output_povray_shape(FILE *f, shape *s, double t){
	switch(s->type){
	case CIRCLE:
		fprintf(f, "cylinder{\n");
		fprintf(f, "\t<0,0,0>, <0,0,%f>, %f\n", t, s->vtab.circle.radius);
		break;
	case ELLIPSE:
		fprintf(f, "cylinder{\n");
		fprintf(f, "\t<0,0,0>, <0,0,%f>, 1\n", t);
		fprintf(f, "\tscale +x*%f\n", s->vtab.ellipse.halfwidth[0]);
		fprintf(f, "\tscale +y*%f\n", s->vtab.ellipse.halfwidth[1]);
		break;
	case RECTANGLE:
		fprintf(f, "box{\n");
		fprintf(f, "\t<-1,-1,0>, <1,1,%f>\n", t);
		fprintf(f, "\tscale +x*%f\n", s->vtab.rectangle.halfwidth[0]);
		fprintf(f, "\tscale +y*%f\n", s->vtab.rectangle.halfwidth[1]);
		break;
	case POLYGON:
		fprintf(f, "prism{\n");
		fprintf(f, "\tlinear_spline\n");
		fprintf(f, "\t0, %f, %d,\n", t, s->vtab.polygon.n_vertices);
		for(int i = 0; i < s->vtab.polygon.n_vertices; ++i){
			fprintf(f, "\t<%f,%f>%c\n",
				s->vtab.polygon.vertex[2*i+0],
				s->vtab.polygon.vertex[2*i+1],
				i+1 == s->vtab.polygon.n_vertices ? ' ' : ','
			);
		}
		break;
	default:
		break;
	}
	fprintf(f, "\trotate +z*%f\n", s->angle * 180./M_PI);
	fprintf(f, "\ttranslate +x*%f\n", s->center[0]);
	fprintf(f, "\ttranslate +y*%f\n", s->center[1]);
	fprintf(f, "}\n");
}

int Simulation_OutputStructurePOVRay(S4_Simulation *S, FILE *fp){
	S4_TRACE("> Simulation_OutputStructurePOVRay(S=%p)\n", S);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(0 != ret){
		S4_TRACE("< Simulation_OutputStructurePOVRay (failed; ret = %d)\n", ret);
		return ret;
	}

	if(NULL == S->solution){
		ret = Simulation_InitSolution(S);
		if(0 != ret){
			S4_TRACE("< Simulation_OutputStructurePOVRay (failed; Simulation_InitSolution returned %d)\n", ret);
			return ret;
		}
	}

	// characteristic unit cell size
	double charsize = 1;
	{
		double u = hypot(S->Lr[0], S->Lr[1]);
		double v = hypot(S->Lr[2], S->Lr[3]);
		charsize = (u > v ? u : v);
	}

	// get Voronoi defining points
	double vorpts[8];
	vorpts[0] = S->Lr[0];
	vorpts[1] = S->Lr[1];
	vorpts[2] = S->Lr[2];
	vorpts[3] = S->Lr[3];
	vorpts[4] = S->Lr[0] + S->Lr[2];
	vorpts[5] = S->Lr[1] + S->Lr[3];
	vorpts[6] = S->Lr[0] - S->Lr[2];
	vorpts[7] = S->Lr[1] - S->Lr[3];
	if(hypot(vorpts[6],vorpts[7]) < hypot(vorpts[4],vorpts[5])){
		vorpts[4] = vorpts[6];
		vorpts[5] = vorpts[7];
	}

	FILE *f = fp;
	if(NULL == fp){ f = stdout; }

	// output preample
	fprintf(f,
		"// -w320 -h240\n"
		"\n"
		"#version 3.6;\n"
		"\n"
		"#include \"colors.inc\"\n"
		"#include \"textures.inc\"\n"
		"#include \"shapes.inc\"\n"
		"\n"
		"global_settings {max_trace_level 5 assumed_gamma 1.0}\n"
		"\n"
		"camera {\n"
		"\tlocation <%f, %f, %f>\n"
		"\tdirection <0, 0,  2.25>\n"
		"\tright x*1.33\n"
		"\tlook_at <0,0,0>\n"
		"}\n"
		"\n",
		-3.0*charsize, 6.0*charsize, -9.0*charsize
	);
	fprintf(f,
		"#declare Dist=80.0;\n"
		"light_source {< -25, 50, -50> color White\n"
		"\tfade_distance Dist fade_power 2\n"
		"}\n"
		"light_source {< 50, 10,  -4> color Gray30\n"
		"\tfade_distance Dist fade_power 2\n"
		"}\n"
		"light_source {< 0, 100,  0> color Gray30\n"
		"\tfade_distance Dist fade_power 2\n"
		"}\n"
		"\n"
		"sky_sphere {\n"
		"\tpigment {\n"
		"\t\tgradient y\n"
		"\t\tcolor_map {\n"
		"\t\t\t[0, 1  color White color White]\n"
		"\t\t}\n"
		"\t}\n"
		"}\n"
		"\n"
		"#declare Xaxis = union{\n"
		"\tcylinder{\n"
		"\t\t<0,0,0>,<0.8,0,0>,0.05\n"
		"\t}\n"
		"\tcone{\n"
		"\t\t<0.8,0,0>, 0.1, <1,0,0>, 0\n"
		"\t}\n"
		"\ttexture { pigment { color Red } }\n"
		"}\n"
		"#declare Yaxis = union{\n"
		"\tcylinder{\n"
		"\t\t<0,0,0>,<0,0.8,0>,0.05\n"
		"\t}\n"
		"\tcone{\n"
		"\t\t<0,0.8,0>, 0.1, <0,1,0>, 0\n"
		"\t}\n"
		"\ttexture { pigment { color Green } }\n"
		"}\n"
		"#declare Zaxis = union{\n"
		"\tcylinder{\n"
		"\t<0,0,0>,<0,0,0.8>,0.05\n"
		"\t}\n"
		"\tcone{\n"
		"\t\t<0,0,0.8>, 0.1, <0,0,1>, 0\n"
		"\t}\n"
		"\ttexture { pigment { color Blue } }\n"
		"}\n"
		"#declare Axes = union{\n"
		"\tobject { Xaxis }\n"
		"\tobject { Yaxis }\n"
		"\tobject { Zaxis }\n"
		"}\n"
	);

	// Output material texture definitions
	for(int i = 0; i < S->n_materials; ++i){
		const S4_Material *M = &(S->material[i]);
		if(NULL != M->name){
			if(0 == strcasecmp("air", M->name) || 0 == strcasecmp("vacuum", M->name)){
				fprintf(f, "#declare Material_%s = texture{ pigment{ color transmit 1.0 } }\n", M->name);
			}else{
				fprintf(f, "#declare Material_%s = texture{ pigment{ rgb <%f,%f,%f> } }\n",
					M->name,
					(double)rand()/(double)RAND_MAX,
					(double)rand()/(double)RAND_MAX,
					(double)rand()/(double)RAND_MAX
				);
			}
		}
	}

	unsigned int layer_name_counter = 0;
	double layer_offset = 0;
	for(int i = 0; i < S->n_layers; ++i){
		const S4_Layer *L = &(S->layer[i]);
		const char *name = L->name;
		if(NULL == name){
			fprintf(f, "#declare Layer_%u = union{\n", layer_name_counter++);
		}else{
			fprintf(f, "#declare Layer_%s = union{\n", name);
		}

		bool comment_out = false;
		if(L->material >= 0){
			if(0 == strcasecmp("air", S->material[L->material].name) || 0 == strcasecmp("vacuum", S->material[L->material].name)){
				comment_out = true;
			}
		}

		if(comment_out){
			fprintf(f, "/*\n");
		}
		fprintf(f, "\tdifference{\n");

		// output base uniform layer
		fprintf(f, "\t\tintersection{\n");
		for(unsigned int i = 0; i < 3; ++i){
			double dist = 0.5*hypot(vorpts[2*i+0], vorpts[2*i+1]);
			fprintf(f, "\t\t\tplane{ <%f,%f,0>, %f }\n",  vorpts[2*i+0],  vorpts[2*i+1], dist);
			fprintf(f, "\t\t\tplane{ <%f,%f,0>, %f }\n", -vorpts[2*i+0], -vorpts[2*i+1], dist);
		}
		fprintf(f, "\t\t\tplane{ <0,0,-1>, 0 }\n");
		fprintf(f, "\t\t\tplane{ <0,0,1>, %f }\n", L->thickness);
		fprintf(f, "\t\t}\n");


		fprintf(f, "// nshapes = %d\n", L->pattern.nshapes);


		for(int i = 0; i < L->pattern.nshapes; ++i){
			if(-1 == L->pattern.parent[i]){
				output_povray_shape(f, &L->pattern.shapes[i], L->thickness);
			}
		}

		if(L->copy >= 0){
			const S4_Layer *Lcopy = &S->layer[L->copy];
			if(NULL != S->material[Lcopy->material].name){
				fprintf(f, "\t\ttexture { Material_%s }\n", S->material[Lcopy->material].name);
			}else{
				fprintf(f, "\t\ttexture { Material_ }\n");
			}
		}else{
			if(NULL != S->material[L->material].name){
				fprintf(f, "\t\ttexture { Material_%s }\n", S->material[L->material].name);
			}
		}
		fprintf(f, "\t}\n");
		if(comment_out){
			fprintf(f, "*/\n");
		}

		for(int i = 0; i < L->pattern.nshapes; ++i){
			comment_out = false;

			const S4_Material *m = &S->material[L->pattern.shapes[i].tag];
			if(NULL != m && NULL != m->name){
				if(0 == strcasecmp("air", m->name) || 0 == strcasecmp("vacuum", m->name)){
					comment_out = true;
				}
			}

			if(comment_out){
				fprintf(f, "/*\n");
			}

			fprintf(f, "\tdifference{\n");
			fprintf(f, "\t\tintersection{\n");
			output_povray_shape(f, &L->pattern.shapes[i], L->thickness);
			for(int j = i+1; j < L->pattern.nshapes; ++j){
				if(i == L->pattern.parent[j]){
					output_povray_shape(f, &L->pattern.shapes[j], L->thickness);
				}
			}
			fprintf(f, "\t\t\tplane{ <0,0,-1>, 0 }\n");
			fprintf(f, "\t\t\tplane{ <0,0,1>, %f }\n", L->thickness);
			fprintf(f, "\t\t}\n");


			fprintf(f, "\t\ttexture { Material_%s }\n", m != NULL ? m->name : "");
			fprintf(f, "\t}\n");
			if(comment_out){
				fprintf(f, "*/\n");
			}
		}

		fprintf(f, "\ttranslate +z*%f\n", layer_offset);
		fprintf(f, "}\n");

		layer_offset += L->thickness;
	}

	// Output postamble
	fprintf(f,
		"#declare Layers = union {\n"
	);
	layer_name_counter = 0;
	for(int i = 0; i < S->n_layers; ++i){
		const S4_Layer *L = &(S->layer[i]);
		const char *name = L->name;
		const char *comment = (0 < i && i+1 < S->n_layers) ? "" : "//";
		if(NULL == name){
			fprintf(f, "\t%sobject{ Layer_%u }\n", comment, layer_name_counter++);
		}else{
			fprintf(f, "\t%sobject{ Layer_%s }\n", comment, name);
		}
	}
	fprintf(f,
		"}\n"
		"\n"
		"Axes\n"
		"Layers\n"
	);
	S4_TRACE("< Simulation_OutputStructurePOVRay\n");
	return 0;
}

int Simulation_OutputLayerPatternDescription(S4_Simulation *S, S4_Layer *layer, FILE *fp){
	S4_TRACE("> Simulation_OutputLayerPatternDescription(S=%p, layer=%p)\n",
		S, layer);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(0 != ret){
		S4_TRACE("< Simulation_OutputLayerPatternDescription (failed; ret = %d)\n", ret);
		return ret;
	}

	const double *L = &S->Lr[0];

	FILE *f = fp;
	if(NULL == fp){ f = stdout; }

	double scale[2] = {4*72, 4*72};
	fprintf(f, "%f %f scale\n", scale[0], scale[1]);
	fprintf(f, "%f setlinewidth\n", 2./scale[0]);
	fprintf(f, "%f %f translate\n", 8.5*0.5*72/scale[0], 11*0.5*72/scale[1]);
	// Draw the origin cross
	fprintf(f, "newpath %f %f moveto %f %f lineto stroke\n", -0.05, 0.00, 0.05, 0.00);
	fprintf(f, "newpath %f %f moveto %f %f lineto stroke\n", 0.00, -0.05, 0.00, 0.05);

	// Draw the clip marks
	fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto\nstroke\n",
		0.55*L[0]+0.50*L[2], 0.55*L[1]+0.50*L[3],
		0.50*L[0]+0.50*L[2], 0.50*L[1]+0.50*L[3],
		0.50*L[0]+0.55*L[2], 0.50*L[1]+0.55*L[3]);
	fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto\nstroke\n",
		-0.50*L[0]+0.55*L[2], -0.50*L[1]+0.55*L[3],
		-0.50*L[0]+0.50*L[2], -0.50*L[1]+0.50*L[3],
		-0.55*L[0]+0.50*L[2], -0.55*L[1]+0.50*L[3]);
	fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto\nstroke\n",
		-0.55*L[0]-0.50*L[2], -0.55*L[1]-0.50*L[3],
		-0.50*L[0]-0.50*L[2], -0.50*L[1]-0.50*L[3],
		-0.50*L[0]-0.55*L[2], -0.50*L[1]-0.55*L[3]);
	fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto\nstroke\n",
		0.50*L[0]-0.55*L[2], 0.50*L[1]-0.55*L[3],
		0.50*L[0]-0.50*L[2], 0.50*L[1]-0.50*L[3],
		0.55*L[0]-0.50*L[2], 0.55*L[1]-0.50*L[3]);

	// Draw the patterns
	for(int i = 0; i < layer->pattern.nshapes; ++i){
		shape *s = &layer->pattern.shapes[i];
		fprintf(f, "gsave %f %f translate %f rotate\n", s->center[0], s->center[1], s->angle * 180./M_PI);
		switch(s->type){
		case CIRCLE:
			fprintf(f, "newpath 0 0 %f 0 360 arc closepath stroke\n", s->vtab.circle.radius);
			break;
		case ELLIPSE:
			fprintf(f, "gsave 1 %f scale ", s->vtab.ellipse.halfwidth[1]/s->vtab.ellipse.halfwidth[0]);
			fprintf(f, "newpath 0 0 %f 0 360 arc closepath stroke ", s->vtab.ellipse.halfwidth[0]);
			fprintf(f, "grestore\n");
			break;
		case RECTANGLE:
			fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath stroke\n",
				-s->vtab.rectangle.halfwidth[0], -s->vtab.rectangle.halfwidth[1],
				 s->vtab.rectangle.halfwidth[0], -s->vtab.rectangle.halfwidth[1],
				 s->vtab.rectangle.halfwidth[0],  s->vtab.rectangle.halfwidth[1],
				-s->vtab.rectangle.halfwidth[0],  s->vtab.rectangle.halfwidth[1]);
			break;
		case POLYGON:
			if(s->vtab.polygon.n_vertices >= 3){
				fprintf(f, "newpath %f %f moveto ", s->vtab.polygon.vertex[0], s->vtab.polygon.vertex[1]);
				for(int j = 1; j < s->vtab.polygon.n_vertices; ++j){
					fprintf(f, "%f %f lineto ", s->vtab.polygon.vertex[2*j+0], s->vtab.polygon.vertex[2*j+1]);
				}
				fprintf(f, "closepath stroke\n");
			}
			break;
		default:
			break;
		}
		fprintf(f, "grestore\n");
	}

	fprintf(f, "showpage\n");


	S4_TRACE("< Simulation_OutputLayerPatternDescription\n");
	return 0;
}

int Simulation_OutputLayerPatternRealization(S4_Simulation *S, S4_Layer *layer, int nx, int ny, FILE *fp){
	S4_TRACE("> Simulation_OutputLayerPatternRealization(S=%p, layer=%p, nx=%d, ny=%d)\n",
		S, layer, nx, ny);
	if(NULL == S){
		S4_TRACE("< Simulation_OutputLayerPatternRealization (failed; S == NULL)\n");
		return -1;
	}
	if(NULL == layer){
		S4_TRACE("< Simulation_OutputLayerPatternRealization (failed; layer == NULL)\n");
		return -2;
	}

	if(NULL == S->solution){
		int error = Simulation_InitSolution(S);
		if(0 != error){
			S4_TRACE("< Simulation_OutputLayerPatternRealization (failed; Simulation_InitSolution returned %d)\n", error);
			return error;
		}
	}

	FILE *f = fp;
	if(NULL == fp){ f = stdout; }

	double *values = (double*)S4_malloc(sizeof(double)*2*(layer->pattern.nshapes+1));
	for(int i = -1; i < layer->pattern.nshapes; ++i){
		const S4_Material *M;
		if(-1 == i){
			M = &S->material[layer->material];
		}else{
			M = &S->material[layer->pattern.shapes[i].tag];
		}
		if(0 == M->type){
			values[2*(i+1)+0] = M->eps.s[0];
			values[2*(i+1)+1] = M->eps.s[1];
		}else{
			values[2*(i+1)+0] = M->eps.abcde[8];
			values[2*(i+1)+1] = M->eps.abcde[9];
		}
	}

	const double unit_cell_size = Simulation_GetUnitCellSize(S);
	const int ndim = (0 == S->Lr[2] && 0 == S->Lr[3]) ? 1 : 2;
	for(int j = 0; j < ny; ++j){
		for(int i = 0; i < nx; ++i){
			double ruv[2] = {
				-0.5 + ((double)i+0.5)/(double)nx,
				-0.5 + ((double)j+0.5)/(double)ny };
			double r[2] = {
				ruv[0] * S->Lr[0] + ruv[1] * S->Lr[2],
				ruv[0] * S->Lr[1] + ruv[1] * S->Lr[3] };
			double z[2] = {0,0};
			for(int g = 0; g < S->n_G; ++g){
				double f[2] = {
					S->G[2*g+0] * S->Lk[0] + S->G[2*g+1] * S->Lk[2],
					S->G[2*g+0] * S->Lk[1] + S->G[2*g+1] * S->Lk[3]
				};
				double ft[2];
				Pattern_GetFourierTransform(&layer->pattern, values, f, ndim, unit_cell_size, ft);

				//double theta = (S->k[0] + f[0])*r[0] + (S->k[1] + f[1])*r[1];
				double theta = (f[0])*r[0] + (f[1])*r[1];
				double ca = cos(2*M_PI*theta);
				double sa = sin(2*M_PI*theta);
				//z += std::complex<double>(ft[0],ft[1]) * std::exp(std::complex<double>(0,theta));
				z[0] += ft[0]*ca - ft[1]*sa;
				z[1] += ft[0]*sa + ft[1]*ca;
			}
			fprintf(f, "%d\t%d\t%f\t%f\n", i, j, z[0], z[1]);
		}
		fprintf(f, "\n");
	}

	S4_free(values);

	S4_TRACE("< Simulation_OutputLayerPatternRealization\n");
	return 0;
}
int Simulation_GetField(S4_Simulation *S, const double r[3], double fE[6], double fH[6]){
	S4_TRACE("> Simulation_GetField(S=%p, r=%p (%f,%f,%f), fE=%p, fH=%p)\n",
		S, r, (NULL == r ? 0 : r[0]), (NULL == r ? 0 : r[1]), (NULL == r ? 0 : r[2]), fE, fH);
	if(NULL == S){
		S4_TRACE("< Simulation_GetField (failed; S == NULL)\n");
		return -1;
	}
	if(NULL == r){
		S4_TRACE("< Simulation_GetField (failed; r == NULL)\n");
		return -2;
	}
	if(NULL == fE && NULL == fH){
		S4_TRACE("< Simulation_GetField (early exit; fE and fH are NULL)\n");
		return 0;
	}

	const size_t n2 = 2*S->n_G;
	const size_t n4 = 2*n2;

	S4_Layer *L = NULL;
	double dz = r[2];
	{
		double z = 0;
		int i;
		for(i = 0; i < S->n_layers && r[2] > z+S->layer[i].thickness; ++i){
			z += S->layer[i].thickness;
			if(i+1 == S->n_layers){ break; }
			dz -= S->layer[i].thickness;
		}
		L = &(S->layer[i]);
	}
	if(NULL == L){
		S4_TRACE("< Simulation_GetField (failed; no layers found)\n");
		return 14;
	}
//fprintf(stderr, "(%f,%f,%f) in %s: dz = %f\n", r[0], r[1], r[2], L->name, dz);

	LayerModes *Lmodes;
	std::complex<double> *Lsoln;
	int ret = Simulation_GetLayerSolution(S, L, &Lmodes, &Lsoln);
	if(0 != ret){
		S4_TRACE("< Simulation_GetField (failed; Simulation_GetLayerSolution returned %d)\n", ret);
		return ret;
	}

	std::complex<double> *ab = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>) * (n4+8*n2));
	if(NULL == ab){
		S4_TRACE("< Simulation_GetField (failed; allocation failed)\n");
		return 1;
	}
	std::complex<double> *work = ab + n4;

	RNP::TBLAS::Copy(n4, Lsoln,1, ab,1);
	//RNP::IO::PrintVector(n4, ab, 1);
	TranslateAmplitudes(S->n_G, Lmodes->q, L->thickness, dz, ab);
	std::complex<double> efield[3], hfield[3];
	GetFieldAtPoint(
		S->n_G, S->kx, S->ky, std::complex<double>(S->omega[0],S->omega[1]),
		Lmodes->q, Lmodes->kp, Lmodes->phi, Lmodes->Epsilon_inv, Lmodes->epstype,
		ab, r, (NULL != fE ? efield : NULL) , (NULL != fH ? hfield : NULL), work);
	if(NULL != fE){
		fE[0] = efield[0].real();
		fE[1] = efield[1].real();
		fE[2] = efield[2].real();
		fE[3] = efield[0].imag();
		fE[4] = efield[1].imag();
		fE[5] = efield[2].imag();
	}
	if(NULL != fH){
		fH[0] = hfield[0].real();
		fH[1] = hfield[1].real();
		fH[2] = hfield[2].real();
		fH[3] = hfield[0].imag();
		fH[4] = hfield[1].imag();
		fH[5] = hfield[2].imag();
	}
	S4_free(ab);

	S4_TRACE("< Simulation_GetField\n");
	return 0;
}
int Simulation_GetFieldPlane(S4_Simulation *S, int nxy[2], double zz, double *E, double *H){
	S4_TRACE("> Simulation_GetFieldPlane(S=%p, nxy=%p (%d,%d), z=%f, E=%p, H=%p)\n",
		S, nxy, (NULL == nxy ? 0 : nxy[0]), (NULL == nxy ? 0 : nxy[1]), zz, E, H);
	if(NULL == S){
		S4_TRACE("< Simulation_GetFieldPlane (failed; S == NULL)\n");
		return -1;
	}
	if(NULL == nxy){
		S4_TRACE("< Simulation_GetFieldPlane (failed; r == NULL)\n");
		return -2;
	}
	if(NULL == E || NULL == H){
		S4_TRACE("< Simulation_GetFieldPlane (early exit; E or H are NULL)\n");
		return 0;
	}

	const size_t n2 = 2*S->n_G;
	const size_t n4 = 2*n2;

	S4_Layer *L = NULL;
	double dz = zz;
	{
		double z = 0;
		int i;
		for(i = 0; i < S->n_layers && zz > z+S->layer[i].thickness; ++i){
			z += S->layer[i].thickness;
			if(i+1 == S->n_layers){ break; }
			dz -= S->layer[i].thickness;
		}
		L = &(S->layer[i]);
	}
	if(NULL == L){
		S4_TRACE("< Simulation_GetFieldPlane (failed; no layers found)\n");
		return 14;
	}
//fprintf(stderr, "(%f,%f,%f) in %s: dz = %f\n", r[0], r[1], r[2], L->name, dz);

	LayerModes *Lmodes;
	std::complex<double> *Lsoln;
	int ret = Simulation_GetLayerSolution(S, L, &Lmodes, &Lsoln);
	if(0 != ret){
		S4_TRACE("< Simulation_GetFieldPlane (failed; Simulation_GetLayerSolution returned %d)\n", ret);
		return ret;
	}

	std::complex<double> *ab = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>) * (n4+8*n2));
	if(NULL == ab){
		S4_TRACE("< Simulation_GetFieldPlane (failed; allocation failed)\n");
		return 1;
	}
	std::complex<double> *work = ab + n4;

	RNP::TBLAS::Copy(n4, Lsoln,1, ab,1);
	//RNP::IO::PrintVector(n4, ab, 1);
	TranslateAmplitudes(S->n_G, Lmodes->q, L->thickness, dz, ab);
	size_t snxy[2] = { (size_t)nxy[0], (size_t)nxy[1] };
	GetFieldOnGrid(
		S->n_G, S->G, S->kx, S->ky, std::complex<double>(S->omega[0],S->omega[1]),
		Lmodes->q, Lmodes->kp, Lmodes->phi, Lmodes->Epsilon_inv, Lmodes->epstype,
		ab, snxy, NULL,
		reinterpret_cast<std::complex<double>*>(E),
		reinterpret_cast<std::complex<double>*>(H)
	);
	S4_free(ab);

	S4_TRACE("< Simulation_GetFieldPlane\n");
	return 0;
}

int Simulation_GetEpsilon(S4_Simulation *S, const double r[3], double eps[2]){
	S4_TRACE("> Simulation_GetEpsilon(S=%p, r=%p (%f,%f,%f), eps=%p)\n",
		S, r, (NULL == r ? 0 : r[0]), (NULL == r ? 0 : r[1]), (NULL == r ? 0 : r[2]), eps);
	if(NULL == S){
		S4_TRACE("< Simulation_GetEpsilon (failed; S == NULL)\n");
		return -1;
	}
	if(NULL == r){
		S4_TRACE("< Simulation_GetEpsilon (failed; r == NULL)\n");
		return -2;
	}
	if(NULL == S->solution){
		int error = Simulation_InitSolution(S);
		if(0 != error){
			S4_TRACE("< Simulation_OutputLayerPatternRealization (failed; Simulation_InitSolution returned %d)\n", error);
			return error;
		}
	}

	const S4_Layer *L = NULL;
	{
		double z = 0;
		int i;
		for(i = 0; i < S->n_layers && r[2] > z+S->layer[i].thickness; ++i){
			z += S->layer[i].thickness;
		}
		if(i >= S->n_layers){ i = S->n_layers-1; }
		L = &(S->layer[i]);
	}
	if(NULL == L){
		S4_TRACE("< Simulation_GetEpsilon (failed; no layers found)\n");
		return 14;
	}

	if(L->copy >= 0){ L = &S->layer[L->copy]; }

	double *values = (double*)S4_malloc(sizeof(double)*2*(L->pattern.nshapes+1));
	for(int i = -1; i < L->pattern.nshapes; ++i){
		const S4_Material *M;
		if(-1 == i){
			M = &S->material[L->material];
		}else{
			M = &(S->material[L->pattern.shapes[i].tag]);
		}
		if(0 == M->type){
			values[2*(i+1)+0] = M->eps.s[0];
			values[2*(i+1)+1] = M->eps.s[1];
		}else{
			values[2*(i+1)+0] = M->eps.abcde[8];
			values[2*(i+1)+1] = M->eps.abcde[9];
		}
	}

	double mp1 = 0;
	int pwr = S->options.lanczos_smoothing_power;
	if(S->options.use_Lanczos_smoothing){
		mp1 = GetLanczosSmoothingOrder(S);
		mp1 *= S->options.lanczos_smoothing_width;
	}
	const double unit_cell_size = Simulation_GetUnitCellSize(S);
	const int ndim = (0 == S->Lr[2] && 0 == S->Lr[3]) ? 1 : 2;
	eps[0] = 0;
	eps[1] = 0;
	for(int g = 0; g < S->n_G; ++g){
		double f[2] = {
			S->G[2*g+0] * S->Lk[0] + S->G[2*g+1] * S->Lk[2],
			S->G[2*g+0] * S->Lk[1] + S->G[2*g+1] * S->Lk[3]
		};

		double ft[2];
		Pattern_GetFourierTransform(&L->pattern, values, f, ndim, unit_cell_size, ft);
		if(S->options.use_Lanczos_smoothing){
			double sigma = GetLanczosSmoothingFactor(mp1, pwr, f);
			ft[0] *= sigma;
			ft[1] *= sigma;
		}
		//printf("g = %d, ft = %g, %g\n", g, ft[0], ft[1]);

		//double theta = (S->k[0] + f[0])*r[0] + (S->k[1] + f[1])*r[1];
		double theta = (f[0])*r[0] + (f[1])*r[1];
		double ca = cos(2*M_PI*theta);
		double sa = sin(2*M_PI*theta);
		//z += std::complex<double>(ft[0],ft[1]) * std::exp(std::complex<double>(0,theta));

		eps[0] += ft[0]*ca - ft[1]*sa;
		eps[1] += ft[0]*sa + ft[1]*ca;
	}
	S4_free(values);

	S4_TRACE("< Simulation_GetEpsilon\n");
	return 0;
}

int Simulation_GetSMatrixDeterminant(S4_Simulation *S, double rmant[2], double *rbase, int *rexpo){
	S4_TRACE("> Simulation_GetSMatrixDeterminant(S=%p, rmand=%p, rbase=%p, rexpo=%p)\n", S, rmant, rbase, rexpo);
	if(NULL == S){ return -1; }
	if(NULL == rmant){ return -2; }
	if(NULL == rbase){ return -3; }
	if(NULL == rexpo){ return -4; }

	const size_t n4 = 4*S->n_G;
	std::complex<double> *M = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>)*n4*n4);
	int ret = Simulation_GetSMatrix(S, 0, -1, M);
	if(0 != ret){
		S4_TRACE("< Simulation_GetSMatrixDeterminant (failed; Simulation_GetSMatrix returned %d)\n", ret);
		return ret;
	}
	std::complex<double> mant;
	double base;
	int expo;
	RNP::TLASupport::Determinant(n4, M, n4, &mant, &base, &expo, NULL);
	S4_free(M);

	rmant[0] = mant.real();
	rmant[1] = mant.imag();
	*rbase = base;
	*rexpo = expo;

	S4_TRACE("< Simulation_GetSMatrixDeterminant\n");
	return 0;
}


// Returns a solution error code
// Tint is a vector of time averaged stress tensor integral
int Simulation_GetStressTensorIntegral(S4_Simulation *S, S4_Layer *layer, double offset, double Tint[6]){
	S4_TRACE("> Simulation_GetStressTensorIntegral(S=%p, layer=%p, offset=%f, Tint=%p)\n", S, layer, offset, Tint);

	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(NULL == Tint){ ret = -4; }
	if(0 != ret){
		S4_TRACE("< Simulation_GetStressTensorIntegral (failed; ret = %d)\n", ret);
		return ret;
	}

	LayerModes *Lmodes;
	std::complex<double> *Lsoln;

	ret = Simulation_GetLayerSolution(S, layer, &Lmodes, &Lsoln);
	if(0 != ret){
		S4_TRACE("< Simulation_GetStressTensorIntegral (failed; Simulation_GetLayerSolution returned %d)\n", ret);
		return ret;
	}

	const int n = S->n_G;
	const int n2 = 2*n;
	const int n4 = 2*n2;

	std::complex<double> *ab = (std::complex<double> *)S4_malloc(sizeof(std::complex<double>) * (n4+8*n2));
	if(NULL == ab){
		S4_TRACE("< Simulation_GetStressTensorIntegral (failed; allocation failed)\n");
		return 1;
	}
	std::complex<double> *work = ab + n4;

	memcpy(ab, Lsoln, sizeof(std::complex<double>) * n4);
	TranslateAmplitudes(n, Lmodes->q, layer->thickness, offset, ab);

	std::complex<double> integral[3];
	GetZStressTensorIntegral(
		n, S->kx, S->ky,
		std::complex<double>(S->omega[0],S->omega[1]),
		Lmodes->q, Lmodes->kp, Lmodes->phi, Lmodes->Epsilon_inv, Lmodes->Epsilon2, Lmodes->epstype, ab, integral, work);
	Tint[0] = integral[0].real();
	Tint[1] = integral[1].real();
	Tint[2] = integral[2].real();
	Tint[3] = integral[0].imag();
	Tint[4] = integral[1].imag();
	Tint[5] = integral[2].imag();

	S4_free(ab);

	S4_TRACE("< Simulation_GetStressTensorIntegral\n");
	return 0;
}

// Returns a solution error code
// which can be 'U', 'E', 'H', 'e'
// 'E' is epsilon*|E|^2, 'H' is |H|^2, 'e' is |E|^2, 'U' is 'E'+'H'
int Simulation_GetLayerVolumeIntegral(S4_Simulation *S, S4_Layer *layer, char which, double integral[2]){
	S4_TRACE("> Simulation_GetLayerVolumeIntegral(S=%p, layer=%p, which=%c, integral=%p)\n", S, layer, which, integral);

	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(which != 'U' && which != 'E' && which != 'H' && which != 'e'){ ret = -3; }
	if(NULL == integral){ ret = -4; }
	if(0 != ret){
		S4_TRACE("< Simulation_GetLayerVolumeIntegral (failed; ret = %d)\n", ret);
		return ret;
	}

	LayerModes *Lmodes;
	std::complex<double> *Lsoln;

	ret = Simulation_GetLayerSolution(S, layer, &Lmodes, &Lsoln);
	if(0 != ret){
		S4_TRACE("< Simulation_GetLayerVolumeIntegral (failed; Simulation_GetLayerSolution returned %d)\n", ret);
		return ret;
	}

	const int n = S->n_G;
	const int n2 = 2*n;
	const int n4 = 2*n2;

	std::complex<double> *work = (std::complex<double> *)S4_malloc(sizeof(std::complex<double>) * (n4*n4));
	if(NULL == work){
		S4_TRACE("< Simulation_GetLayerVolumeIntegral (failed; allocation failed)\n");
		return 1;
	}

	std::complex<double> zintegral;
	GetLayerVolumeIntegral(which,
		n, S->kx, S->ky,
		std::complex<double>(S->omega[0],S->omega[1]),
		layer->thickness, Lmodes->q, Lmodes->kp, Lmodes->phi, Lmodes->Epsilon_inv, Lmodes->Epsilon2, Lmodes->epstype, Lsoln, &zintegral, work);

	S4_free(work);

	integral[0] = zintegral.real();
	integral[1] = zintegral.imag();

	S4_TRACE("< Simulation_GetLayerVolumeIntegral\n");
	return 0;
}

// Returns a solution error code
int Simulation_GetLayerZIntegral(S4_Simulation *S, S4_Layer *layer, const double r[2], double integral[6]){
	S4_TRACE("> Simulation_GetLayerZIntegral(S=%p, layer=%p, r=%g,%g, integral=%p)\n", S, layer, r[0], r[1], integral);

	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(NULL == r){ ret = -3; }
	if(NULL == integral){ ret = -4; }
	if(0 != ret){
		S4_TRACE("< Simulation_GetLayerZIntegral (failed; ret = %d)\n", ret);
		return ret;
	}

	LayerModes *Lmodes;
	std::complex<double> *Lsoln;

	ret = Simulation_GetLayerSolution(S, layer, &Lmodes, &Lsoln);
	if(0 != ret){
		S4_TRACE("< Simulation_GetLayerZIntegral (failed; Simulation_GetLayerSolution returned %d)\n", ret);
		return ret;
	}

	const int n = S->n_G;
	const int n2 = 2*n;
	const int n4 = 2*n2;

	std::complex<double> *work = (std::complex<double> *)S4_malloc(sizeof(std::complex<double>) * (12*n4));
	if(NULL == work){
		S4_TRACE("< Simulation_GetLayerZIntegral (failed; allocation failed)\n");
		return 1;
	}

	GetLayerZIntegral(
		n, S->kx, S->ky,
		std::complex<double>(S->omega[0],S->omega[1]),
		layer->thickness, r, Lmodes->q, Lmodes->kp, Lmodes->phi, Lmodes->Epsilon_inv, Lmodes->Epsilon2, Lmodes->epstype, Lsoln, integral, work);

	S4_free(work);


	S4_TRACE("< Simulation_GetLayerZIntegral\n");
	return 0;
}

int Simulation_MakeExcitationPlanewave(S4_Simulation *S, const double angle[2], const double pol_s[2], const double pol_p[2], size_t order){
	S4_TRACE("> Simulation_MakeExcitationPlanewave(S=%p, angle=%p (%f,%f), pol_s=%p (%f,%f), pol_p=%p (%f,%f))\n", S,
		angle, (NULL != angle) ? angle[0] : 0, (NULL != angle) ? angle[1] : 0,
		pol_s, (NULL != pol_s) ? pol_s[0] : 0, (NULL != pol_s) ? pol_s[1] : 0,
		pol_p, (NULL != pol_p) ? pol_p[0] : 0, (NULL != pol_p) ? pol_p[1] : 0);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == angle){ ret = -2; }
	if(NULL == pol_s){ ret = -3; }
	if(NULL == pol_p){ ret = -4; }
	if(0 != ret){
		S4_TRACE("< Simulation_MakeExcitationPlanewave (failed; ret = %d)\n", ret);
		return ret;
	}

	if(NULL == S->layer){
		S4_TRACE("< Simulation_MakeExcitationPlanewave (failed; no layers)\n");
		return 14;
	}

	Simulation_DestroySolution(S);
	S->exc.type = 0;

	const S4_Material *M = &S->material[S->layer[0].material];
	if(NULL == M){
		S4_TRACE("< Simulation_MakeExcitationPlanewave (failed; material %d not defined)\n", S->layer[0].material);
		return 15;
	}
	std::complex<double> layer_eps;
	if(0 == M->type){
		layer_eps = std::complex<double>(M->eps.s[0], M->eps.s[1]);
	}else{
		layer_eps = std::complex<double>(M->eps.abcde[8], M->eps.abcde[9]);
	}
	double root_eps = std::sqrt(layer_eps).real();

	// At normal incidence, we assume
	// {0,0,k} is the k-vector
	// {e_s, e_p, 0} is the electric field
	// e_s is the s-polarization complex ampltiude
	// so then the magnetic field is
	// {-e_p, e_s, 0}
	// Now we apply a rotation in the xz-plane by angle[0], then in the xy-plane by angle[1]
	// [ cos1 -sin1  0 ] [  cos0  0  sin0 ] [ x ]   [  x*cos0*cos1 - y*sin1 + z*cos1*sin0 ]
	// [ sin1  cos1  0 ] [   0    1   0   ] [ y ] = [  x*cos0*sin1 + y*cos1 + z*sin1*sin0 ]
	// [  0     0    1 ] [ -sin0  0  cos0 ] [ z ]   [ -x*sin0               + z*cos0      ]
	const double c0 = cos(angle[0]);
	const double s0 = sin(angle[0]);
	const double c1 = cos(angle[1]);
	const double s1 = sin(angle[1]);
	const double c_s = cos(pol_s[1]);
	const double s_s = sin(pol_s[1]);
	const double c_p = cos(pol_p[1]);
	const double s_p = sin(pol_p[1]);
	S->k[0] = c1*s0*root_eps;
	S->k[1] = s1*s0*root_eps;

	S->exc.sub.planewave.hx[0] = -c0*c1*pol_s[0]*c_s - s1*pol_p[0]*c_p;
	S->exc.sub.planewave.hx[1] = -c0*c1*pol_s[0]*s_s - s1*pol_p[0]*s_p;
	S->exc.sub.planewave.hy[0] = -c0*s1*pol_s[0]*c_s + c1*pol_p[0]*c_p;
	S->exc.sub.planewave.hy[1] = -c0*s1*pol_s[0]*s_s + c1*pol_p[0]*s_p;
	S->exc.sub.planewave.order = order;
	S->exc.sub.planewave.backwards = false;
	/*
	S->ex[0] = ;
	S->ex[1] = ;
	S->ey[0] = ;
	S->ey[1] = ;
	*/
	/*
	std::cout << "k: " << S->k[0] << "\t" << S->k[1] << std::endl;
	std::cout << "hx: " << S->exc.sub.planewave.hx[0] << "\t" << S->exc.sub.planewave.hx[1] << std::endl;
	std::cout << "hy: " << S->exc.sub.planewave.hy[0] << "\t" << S->exc.sub.planewave.hy[1] << std::endl;
	//*/

	S4_TRACE("< Simulation_MakeExcitationPlanewave\n");
	return 0;
}

int S4_Simulation_ExcitationDipole(S4_Simulation *S, const double k[2], const char *layer, const double pos[2], const double moment[6]){
	S4_TRACE("> Simulation_MakeExcitationDipole(S=%p, k=%p (%f,%f), moment=%p (%f,%f,%f,%f,%f,%f))\n", S,
		pos, (NULL != pos) ? pos[0] : 0, (NULL != pos) ? pos[1] : 0,
		moment, (NULL != moment) ? moment[0] : 0, (NULL != moment) ? moment[1] : 0, (NULL != moment) ? moment[2] : 0,
		(NULL != moment) ? moment[3] : 0, (NULL != moment) ? moment[4] : 0, (NULL != moment) ? moment[5] : 0);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == pos){ ret = -2; }
	if(NULL == moment){ ret = -3; }
	if(0 != ret){
		S4_TRACE("< Simulation_MakeExcitationDipole (failed; ret = %d)\n", ret);
		return ret;
	}

	if(NULL == S->layer){
		S4_TRACE("< Simulation_MakeExcitationDipole (failed; no layers)\n");
		return 14;
	}

	Simulation_DestroySolution(S);
	S->exc.type = 1;

	S->k[0] = k[0];
	S->k[1] = k[1];
	S->exc.sub.dipole.pos[0] = pos[0];
	S->exc.sub.dipole.pos[1] = pos[1];
	for(size_t i = 0; i < 6; ++i){
		S->exc.sub.dipole.moment[i] = moment[i];
	}
	S->exc.layer = Simulation_GetLayerByName(S, layer, NULL);

	S4_TRACE("< Simulation_MakeExcitationDipole\n");
	return 0;
}

int S4_Simulation_ExcitationExterior(S4_Simulation *S, int n, const int *exg, const double *ex){
	S4_TRACE("> Simulation_MakeExcitationExterior(S=%p, n=%d, exg=%p, ex=%p)\n", S,
		n, exg, ex);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(0 == n){ ret = -2; }
	if(NULL == exg){ ret = -3; }
	if(NULL == ex){ ret = -4; }
	if(0 != ret){
		S4_TRACE("< Simulation_MakeExcitationExterior (failed; ret = %d)\n", ret);
		return ret;
	}

	Simulation_DestroySolution(S);
	S->exc.type = 2;

	Excitation_Exterior *ext = &(S->exc.sub.exterior);

	ext->n = n;
	ext->Gindex1 = (int*)malloc(sizeof(int) * 2*n);
	memcpy(ext->Gindex1, exg, sizeof(int) * 2*n);
	ext->coeff = (double*)malloc(sizeof(double) * 2*n);
	memcpy(ext->coeff, ex, sizeof(double) * 2*n);

	S4_TRACE("< Simulation_MakeExcitationExterior\n");
	return 0;
}

void Simulation_InvalidateFieldCache(S4_Simulation *S){
	S4_TRACE("> Simulation_InvalidateFieldCache(S=%p) [omega=%f]\n", S, S->omega[0]);
	while(NULL != S->field_cache){
		FieldCache *t = S->field_cache;
		S->field_cache = S->field_cache->next;
		// The P pointer is just t+1, so don't S4_free it!
		S4_free(t);
	}
	S->field_cache = NULL;
	S4_TRACE("< Simulation_InvalidateFieldCache [omega=%f]\n", S->omega[0]);
}
std::complex<double>* Simulation_GetCachedField(const S4_Simulation *S, const S4_Layer *layer){
	S4_TRACE("> Simulation_GetCachedField(S=%p, layer=%p) [omega=%f]\n", S, layer, S->omega[0]);
	std::complex<double> *P = NULL;
	FieldCache *f = S->field_cache;
	while(NULL != f){
		if(layer == f->layer && S->n_G == f->n){
			P = f->P;
			break;
		}
		f = f->next;
	}
	S4_TRACE("< Simulation_GetCachedField returning P = %p [omega=%f]\n", P, S->omega[0]);
	return P;
}
void Simulation_AddFieldToCache(S4_Simulation *S, const S4_Layer *layer, size_t n, const std::complex<double> *P, size_t Plen){
	S4_TRACE("> Simulation_AddFieldToCache(S=%p, layer=%p, n=%d, P=%p) [omega=%f]\n", S, layer, (int)n, P, S->omega[0]);
	FieldCache *f = (FieldCache*)S4_malloc(sizeof(FieldCache)+sizeof(std::complex<double>)*Plen);
	f->P = (std::complex<double>*)(f+1);
	memcpy(f->P, P, sizeof(std::complex<double>)*Plen);
	f->layer = layer;
	f->n = n;
	f->next = S->field_cache;
	S->field_cache = f;
	S4_TRACE("< Simulation_AddFieldToCache [omega=%f]\n", S->omega[0]);
}

void Simulation_SetExcitationType(S4_Simulation *S, int type){
	S4_TRACE("> Simulation_SetExcitationType(S=%p, type=%d\n", S, type);
	if(1 == S->exc.type){
		if(NULL != S->exc.layer){
			free(S->exc.layer);
		}
	}else if(2 == S->exc.type){
		if(NULL != S->exc.sub.exterior.Gindex1){ free(S->exc.sub.exterior.Gindex1); }
		if(NULL != S->exc.sub.exterior.coeff){ free(S->exc.sub.exterior.coeff); }
	}
	S->exc.type = type;
	S4_TRACE("< Simulation_SetExcitationType\n");
}
void Simulation_CopyExcitation(const S4_Simulation *from, S4_Simulation *to){
	S4_TRACE("> Simulation_CopyExcitation(from=%p, to=%p\n", from, to);
	memcpy(&(to->exc), &(from->exc), sizeof(Excitation));
	if(1 == to->exc.type){
		to->exc.layer = to->layer + (from->exc.layer - from->layer);
	}else if(2 == to->exc.type){
		const int n = from->exc.sub.exterior.n;
		to->exc.sub.exterior.Gindex1 = (int*)malloc(sizeof(int) * 2*n);
		memcpy(to->exc.sub.exterior.Gindex1, from->exc.sub.exterior.Gindex1, sizeof(int) * 2*n);
		to->exc.sub.exterior.coeff = (double*)malloc(sizeof(double) * 2*n);
		memcpy(to->exc.sub.exterior.coeff, from->exc.sub.exterior.coeff, sizeof(double) * 2*n);
	}
	S4_TRACE("< Simulation_CopyExcitation\n");
}

int Simulation_GetSMatrix(S4_Simulation *S, int from, int to, std::complex<double> *M){
	S4_TRACE("> Simulation_GetSMatrix(S=%p, from=%d, to=%d)\n", S, from, to);

	if(-1 != to && to < from){ return -3; }

	if(NULL == S->solution){
		int error = Simulation_InitSolution(S);
		if(0 != error){
			S4_TRACE("< Simulation_GetSMatrix (failed; Simulation_InitSolution returned %d)\n", error);
			return error;
		}
	}

	// compute all modes
	for(int i = 0; i < S->n_layers; ++i){
		S4_Layer *SL = &(S->layer[i]);
		if(from <= i && (-1 == to || i <= to)){
			if(NULL == SL->modes){
				Simulation_ComputeLayerModes(S, SL, &SL->modes);
			}
		}
	}

	// Make arrays of q, kp, and phi
	double *lthick = (double*)S4_malloc(sizeof(double)*S->n_layers);
	int *lepstype = (int*)S4_malloc(sizeof(int)*S->n_layers);
	const std::complex<double> **lq   = (const std::complex<double> **)S4_malloc(sizeof(const std::complex<double> *)*S->n_layers*4);
	const std::complex<double> **lepsinv  = lq  + S->n_layers;
	const std::complex<double> **lkp  = lepsinv  + S->n_layers;
	const std::complex<double> **lphi = lkp + S->n_layers;

	for(int i = 0; i < S->n_layers; ++i){
		S4_Layer *SL = &(S->layer[i]);
		if(from <= i && (-1 == to || i <= to)){
			lthick[i] = SL->thickness;
			lq  [i] = SL->modes->q;
			lepsinv[i] = SL->modes->Epsilon_inv;
			lepstype[i] = SL->modes->epstype;
			lkp [i] = SL->modes->kp;
			lphi[i] = SL->modes->phi;
		}
	}

	GetSMatrix(S->n_layers, S->n_G, S->kx, S->ky, std::complex<double>(S->omega[0], S->omega[1]), lthick, lq, lepsinv, lepstype, lkp, lphi, M);

	S4_free(lq);
	S4_free(lepstype);
	S4_free(lthick);

	S4_TRACE("< Simulation_GetSMatrix\n");
	return 0;
}





int S4_Simulation_ExcitationPlanewave(
	S4_Simulation *S, const S4_real *kdir, const S4_real *udir,
	const S4_real *amp_u, const S4_real *amp_v
){
	S4_TRACE("> S4_Simulation_ExcitationPlanewave(S=%p, kdir=(%f,%f,%f), udir=(%f,%f,%f), amp_u=(%f,%f), amp_v=(%f,%f))\n", S,
		kdir[0], kdir[1], kdir[2],
		udir[0], udir[1], udir[2],
		amp_u[0], amp_u[1],
		amp_v[0], amp_v[1]
	);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == kdir){ ret = -2; }
	if(NULL == udir){ ret = -3; }
	if(NULL == amp_u){ ret = -4; }
	if(NULL == amp_v){ ret = -5; }
	if(0 != ret){
		S4_TRACE("< S4_Simulation_ExcitationPlanewave (failed; ret = %d)\n", ret);
		return ret;
	}

	if(S->n_layers <= 0){
		S4_TRACE("< S4_Simulation_ExcitationPlanewave (failed; no layers)\n");
		return 14;
	}
	Simulation_DestroySolution(S);
	S->exc.type = 0;

	const S4_Material *M;
	if(kdir[2] < 0){
		M = &S->material[S->layer[S->n_layers-1].material];
	}else{
		M = &S->material[S->layer[0].material];
	}

	if(NULL == M){
		S4_TRACE("< S4_Simulation_ExcitationPlanewave (failed; material %d not defined)\n", S->layer[0].material);
		return 15;
	}
	std::complex<double> layer_eps;
	if(0 == M->type){
		layer_eps = std::complex<double>(M->eps.s[0], M->eps.s[1]);
	}else{
		layer_eps = std::complex<double>(M->eps.abcde[8], M->eps.abcde[9]);
	}
	double root_eps = std::sqrt(layer_eps).real();

	S4_real il = 1./sqrt(kdir[0]*kdir[0] + kdir[1]*kdir[1] + kdir[2]*kdir[2]);
	const S4_real kn[3] = {
		il*kdir[0], il*kdir[1], il*kdir[2]
	};
	S4_real uk = kn[0]*udir[0] + kn[1]*udir[1] + kn[2]*udir[2];
	S4_real un[3] = {
		udir[0] - uk*kdir[0],
		udir[1] - uk*kdir[1],
		udir[2] - uk*kdir[2]
	};
	il = 1./sqrt(udir[0]*udir[0] + udir[1]*udir[1] + udir[2]*udir[2]);
	un[0] *= il;
	un[1] *= il;
	un[2] *= il;
	const S4_real vn[3] = {
		kn[1]*un[2]-kn[2]*un[1],
		kn[2]*un[0]-kn[0]*un[2],
		kn[0]*un[1]-kn[1]*un[0]
	};
	// E field is amp_u * u + amp_v * v
	// H field is amp_u * v - amp_v * u

	S4_real k_new[2] = { root_eps * kn[0], root_eps * kn[1] };

	if(k_new[0] != S->k[0] || k_new[1] != S->k[1]){
		S4_Simulation_DestroyLayerModes(S, -1);
		S->k[0] = k_new[0];
		S->k[1] = k_new[1];
	}

	S->exc.sub.planewave.hx[0] = root_eps*(amp_u[0]*vn[0] - amp_v[0]*un[0]);
	S->exc.sub.planewave.hx[1] = root_eps*(amp_u[1]*vn[0] - amp_v[1]*un[0]);
	S->exc.sub.planewave.hy[0] = root_eps*(amp_u[0]*vn[1] - amp_v[0]*un[1]);
	S->exc.sub.planewave.hy[1] = root_eps*(amp_u[1]*vn[1] - amp_v[1]*un[1]);
	S->exc.sub.planewave.order = 0;
	if(kn[2] < 0){
		S->exc.sub.planewave.backwards = 1;
	}

	S4_TRACE("< S4_Simulation_ExcitationPlanewave\n");
	return 0;
}

int S4_Simulation_GetFieldPlane(S4_Simulation *S, const int nxy[2], const S4_real *xyz0, S4_real *E, S4_real *H){
	S4_TRACE("> S4_Simulation_GetFieldPlane(S=%p, nxy=%p (%d,%d), r0=(%f,%f,%f), E=%p, H=%p)\n",
		S, nxy, (NULL == nxy ? 0 : nxy[0]), (NULL == nxy ? 0 : nxy[1]), xyz0[0], xyz0[1], xyz0[2], E, H);
	if(NULL == S){
		S4_TRACE("< S4_Simulation_GetFieldPlane (failed; S == NULL)\n");
		return -1;
	}
	if(NULL == nxy){
		S4_TRACE("< S4_Simulation_GetFieldPlane (failed; nxy == NULL)\n");
		return -2;
	}
	if(NULL == xyz0){
		S4_TRACE("< S4_Simulation_GetFieldPlane (failed; xyz0 == NULL)\n");
		return -3;
	}
	if(NULL == E && NULL == H){
		S4_TRACE("< S4_Simulation_GetFieldPlane (early exit; E and H both NULL)\n");
		return 0;
	}
	
	if(1 == nxy[0] && 1 == nxy[1]){
		int ret = Simulation_GetField(S, xyz0, E, H);
		if(0 == ret){
			std::swap(E[1], E[3]);
			std::swap(E[2], E[4]);
			std::swap(E[2], E[3]);
			std::swap(H[1], H[3]);
			std::swap(H[2], H[4]);
			std::swap(H[2], H[3]);
		}
		return ret;
	}

	const size_t n2 = 2*S->n_G;
	const size_t n4 = 2*n2;

	S4_Layer *L = NULL;
	double dz = xyz0[2];
	{
		double z = 0;
		int i;
		for(i = 0; i < S->n_layers && xyz0[2] > z+S->layer[i].thickness; ++i){
			z += S->layer[i].thickness;
			if(i+1 == S->n_layers){ break; }
			dz -= S->layer[i].thickness;
		}
		L = &(S->layer[i]);
	}
	if(NULL == L){
		S4_TRACE("< S4_Simulation_GetFieldPlane (failed; no layers found)\n");
		return 14;
	}
//fprintf(stderr, "(%f,%f,%f) in %s: dz = %f\n", r[0], r[1], r[2], L->name, dz);

	LayerModes *Lmodes;
	std::complex<double> *Lsoln;
	int ret = Simulation_GetLayerSolution(S, L, &Lmodes, &Lsoln);
	if(0 != ret){
		S4_TRACE("< S4_Simulation_GetFieldPlane (failed; Simulation_GetLayerSolution returned %d)\n", ret);
		return ret;
	}

	std::complex<double> *ab = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>) * (n4+8*n2));
	if(NULL == ab){
		S4_TRACE("< S4_Simulation_GetFieldPlane (failed; allocation failed)\n");
		return 1;
	}
	std::complex<double> *work = ab + n4;

	RNP::TBLAS::Copy(n4, Lsoln,1, ab,1);
	//RNP::IO::PrintVector(n4, ab, 1);
	TranslateAmplitudes(S->n_G, Lmodes->q, L->thickness, dz, ab);
	const size_t snxy[2] = { (size_t)nxy[0], (size_t)nxy[1] };
	const double xy0[2] = { xyz0[0], xyz0[1] };
	GetFieldOnGrid(
		S->n_G, S->G, S->kx, S->ky, std::complex<double>(S->omega[0],S->omega[1]),
		Lmodes->q, Lmodes->kp, Lmodes->phi, Lmodes->Epsilon_inv, Lmodes->epstype,
		ab, snxy, xy0,
		reinterpret_cast<std::complex<double>*>(E),
		reinterpret_cast<std::complex<double>*>(H)
	);
	S4_free(ab);

	S4_TRACE("< S4_Simulation_GetFieldPlane\n");
	return 0;
}