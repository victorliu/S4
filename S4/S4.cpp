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

#include <IO.h>
#include <cstdio>

#include <numalloc.h>


void* S4_malloc(size_t size){ // for debugging
	void* ret = malloc_aligned(size, 16);
	//memset(ret, 0x0, size);
	return ret;
}
void S4_free(void *ptr){
	free_aligned(ptr);
}

struct LayerBands{
	std::complex<double> *q; // length 2*glist.n
	std::complex<double> *kp; // size (2*glist.n)^2 (k-parallel matrix)
	std::complex<double> *phi; // size (2*glist.n)^2
	std::complex<double> *Epsilon2; // size (2*glist.n)^2 (dielectric/normal-field matrix)
	std::complex<double> *Epsilon_inv; // size (glist.n)^2 inverse of usual dielectric Fourier coupling matrix
	// max total size needed: 2n+13nn
	int epstype;
};
struct LayerSolution{
	std::complex<double> *ab; // length 2*glist.n
	// total size needed: 4n
};

// This structure caches the Fourier transform of the polarization basis
// field, allowing a substantial increase in speed when doing frequency
// scans. Invalidates on SetNumG, and layer patterning
struct FieldCache{
	int n;
	const Layer *layer;
	std::complex<double> *P; // 2n x 2n matrix, allocated along with this structure itself
	FieldCache *next;
};

// Private functions

// If layer_solution is null, only computes layer_bands if it is non-NULL.
// If layer_solution is not null, all layer's bands are computed.
int Simulation_GetLayerSolution(Simulation *S, Layer *layer, LayerBands **layer_bands, LayerSolution **layer_solution);

// These two assume S->solution is set up already
int Simulation_ComputeLayerSolution(Simulation *S, Layer *L, LayerBands **layer_bands, LayerSolution **layer_solution);
int Simulation_ComputeLayerBands(Simulation *S, Layer *L, LayerBands **bands);
void Simulation_SetExcitationType(Simulation *S, int type);
void Simulation_CopyExcitation(const Simulation *from, Simulation *to);
int Simulation_GetSMatrix(Simulation *S, int layer_from, int layer_to, std::complex<double> *M);

// Field cache manipulation
void Simulation_InvalidateFieldCache(Simulation *S);
std::complex<double>* Simulation_GetCachedField(Simulation *S, const Layer *layer);
void Simulation_AddFieldToCache(Simulation *S, const Layer *layer, size_t n, const std::complex<double> *P, size_t Plen);

void Layer_Init(Layer *L, const char *name, double thickness, const char *material, const char *copy){
	S4_TRACE("> Layer_Init(L=%p, name=%p (%s), thickness=%f, material=%p (%s), copy=%p (%s))\n",
		L, name, (NULL == name ? "" : name), thickness,
		material, (NULL == material ? "" : material),
		copy, (NULL == copy ? "" : copy));
	if(NULL == L){
		S4_TRACE("< Layer_Init (failed; L == NULL)\n");
		return;
	}
	if(NULL != copy){
		L->copy = strdup(copy);
		L->material = NULL;
	}else{
		if(NULL != material){ L->material = strdup(material); }else{ L->material = NULL; }
		L->copy = NULL;
	}
	if(NULL != name){ L->name = strdup(name); }else{ L->name = NULL; }
	L->thickness = thickness;
	L->pattern.nshapes = 0;
	L->pattern.shapes = NULL;
	L->pattern.parent = NULL;
	S4_TRACE("< Layer_Init\n");
}
void Layer_Destroy(Layer *L){
	S4_TRACE("> Layer_Destroy(L=%p)\n", L);
	if(NULL == L){
		S4_TRACE("< Layer_Destroy (failed; L == NULL)\n");
		return;
	}
	free(L->name); L->name = NULL;
	free(L->material); L->material = NULL;
	if(NULL != L->copy){ free(L->copy); L->copy = NULL; }
	if(NULL != L->pattern.shapes){ free(L->pattern.shapes); L->pattern.shapes = NULL; }
	if(NULL != L->pattern.parent){ free(L->pattern.parent); L->pattern.parent = NULL; }
	S4_TRACE("< Layer_Destroy\n");
}
void Material_Init(Material *M, const char *name, const double eps[2]){
	S4_TRACE("> Material_Init(M=%p, name=%p (%s), eps=%p (%f,%f))\n",
		M, name, (NULL == name ? "" : name),
		eps, (NULL == eps ? 0 : eps[0]), (NULL == eps ? 0 : eps[1]));
	if(NULL == M){
		S4_TRACE("< Material_Init (failed; M == NULL)\n");
		return;
	}
	if(NULL != name){ M->name = strdup(name); }else{ M->name = NULL; }
	M->type = 0;
	if(NULL != eps){
		M->eps.s[0] = eps[0];
		M->eps.s[1] = eps[1];
	}else{ M->eps.s[0] = 1; M->eps.s[1] = 0; }
	S4_TRACE("< Material_Init\n");
}
void Material_InitTensor(Material *M, const char *name, const double abcde[10]){
	S4_TRACE("> Material_InitTensor(M=%p, name=%p (%s), eps=%p (%f,%f))\n",
		M, name, (NULL == name ? "" : name),
		abcde, (NULL == abcde ? 0 : abcde[8]), (NULL == abcde ? 0 : abcde[9]));
	if(NULL == M){
		S4_TRACE("< Material_InitTensor (failed; M == NULL)\n");
		return;
	}
	if(NULL != name){ M->name = strdup(name); }else{ M->name = NULL; }
	M->type = 1;
	if(NULL != abcde){
		for(size_t i = 0; i < 10; ++i){	M->eps.abcde[i] = abcde[i]; }
	}else{
		for(size_t i = 0; i < 10; ++i){	M->eps.abcde[i] = 0; }
		M->eps.abcde[0] = 1;
		M->eps.abcde[6] = 1;
		M->eps.abcde[8] = 1;
	}
	S4_TRACE("< Material_InitTensor\n");
}
void Material_Destroy(Material *M){
	S4_TRACE("> Material_Destroy(M=%p)\n", M);
	if(NULL == M){
		S4_TRACE("< Material_Destroy (failed; M == NULL)\n");
		return;
	}
	free(M->name); M->name = NULL;
	S4_TRACE("< Material_Destroy\n");
}


void Simulation_Init(Simulation *S){
	S4_TRACE("> Simulation_Init(S=%p)\n", S);
	S->solution = NULL; // needed by other initializations

	S->Lr[0] = 1;
	S->Lr[1] = 0;
	S->Lr[2] = 0;
	S->Lr[3] = 1;
	Simulation_MakeReciprocalLattice(S);
	S->n_G = 0;
	S->material = NULL;
	S->layer = NULL;
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

	S4_TRACE("< Simulation_Init\n");
}

void Simulation_Destroy(Simulation *S){
	S4_TRACE("> Simulation_Destroy(S=%p) [omega=%f]\n", S, S->omega[0]);
	if(NULL == S){
		S4_TRACE("< Simulation_Destroy (failed; S == NULL) [omega=%f]\n", S->omega[0]);
		return;
	}
	if(NULL != S->solution){
		Simulation_DestroySolution(S);
	}
	while(NULL != S->layer){
		Layer *l = S->layer;
		S->layer = S->layer->next;
		Layer_Destroy(l);
		S4_free(l);
	}
	while(NULL != S->material){
		Material *m = S->material;
		S->material = S->material->next;
		Material_Destroy(m);
		S4_free(m);
	}
	Simulation_SetExcitationType(S, -1);
	Simulation_InvalidateFieldCache(S);
	if(NULL != S->options.vector_field_dump_filename_prefix){
		free(S->options.vector_field_dump_filename_prefix);
		S->options.vector_field_dump_filename_prefix = NULL;
	}
	S4_TRACE("< Simulation_Destroy [omega=%f]\n", S->omega[0]);
}
void Simulation_Clone(const Simulation *S, Simulation *T){
	S4_TRACE("> Simulation_Clone(S=%p, T=%p) [omega=%f]\n", S, T, S->omega[0]);
	if(NULL == S || NULL == T){
		S4_TRACE("< Simulation_Clone (failed; S == NULL or T == NULL) [omega=%f]\n", S->omega[0]);
		return;
	}

	memcpy(T, S, sizeof(Simulation));

	Layer *L = S->layer;
	Layer **L2 = &T->layer;
	while(NULL != L){
		*L2 = (Layer*)S4_malloc(sizeof(Layer));
		Layer_Init(*L2, L->name, L->thickness, L->material, L->copy);

		// Copy pattern
		(*L2)->pattern.nshapes = L->pattern.nshapes;
		(*L2)->pattern.shapes = (shape*)malloc(sizeof(shape)*L->pattern.nshapes);
		memcpy((*L2)->pattern.shapes, L->pattern.shapes, sizeof(shape)*L->pattern.nshapes);
		(*L2)->pattern.parent = NULL;

		(*L2)->next = NULL;
		L2 = &((*L2)->next);
		L = L->next;
	}

	Material *M = S->material;
	Material **M2 = &T->material;
	while(NULL != M){
		*M2 = (Material*)S4_malloc(sizeof(Material));
		if(0 == M->type){
			Material_Init(*M2, M->name, M->eps.s);
		}else{
			Material_InitTensor(*M2, M->name, M->eps.abcde);
		}

		(*M2)->next = NULL;
		M2 = &((*M2)->next);
		M = M->next;
	}

	Simulation_CopyExcitation(S, T);

	T->solution = NULL;
	T->field_cache = NULL;

	S4_TRACE("< Simulation_Clone [omega=%f]\n", S->omega[0]);
}

int Simulation_MakeReciprocalLattice(Simulation *S){
	S4_TRACE("> Simulation_MakeReciprocalLattice(S=%p) [omega=%f]\n", S, S->omega[0]);
	double d;
	if(NULL == S){
		S4_TRACE("< Simulation_MakeReciprocalLattice (failed; S == NULL) [omega=%f]\n", S->omega[0]);
		return -1;
	}

	Simulation_DestroySolution(S);

	d = S->Lr[0]*S->Lr[3] - S->Lr[1]*S->Lr[2];
	if(0 == d){
		if(0 != S->Lr[2] || 0 != S->Lr[3]){ // 1D lattice has zero vector for second basis vector
			S4_TRACE("< Simulation_MakeReciprocalLattice (failed; degenerate lattice basis)\n");
			return 1; // degenerate lattice basis vectors
		}
		d = hypot(S->Lr[0], S->Lr[1]);
		if(0 == d){
			S4_TRACE("< Simulation_MakeReciprocalLattice (failed; both lattice vectors are zero)\n");
			return 2; // both basis vectors are zero
		}
		d = 1./(d*d);
		S->Lk[0] = S->Lr[0] * d;
		S->Lk[1] = S->Lr[1] * d;
		S->Lk[2] = 0; S->Lk[3] = 0;
	}else{
		if(d < 0){
			double t;
			t = S->Lr[0]; S->Lr[0] = S->Lr[2]; S->Lr[2] = t;
			t = S->Lr[1]; S->Lr[1] = S->Lr[3]; S->Lr[3] = t;
			d = -d;
		}
		d = 1./d;
		S->Lk[0] =  d*S->Lr[3];
		S->Lk[1] = -d*S->Lr[2];
		S->Lk[2] = -d*S->Lr[1];
		S->Lk[3] =  d*S->Lr[0];
	}
	S4_TRACE("< Simulation_MakeReciprocalLattice [omega=%f]\n", S->omega[0]);
	return 0;
}

Material* Simulation_AddMaterial(Simulation *S){
	S4_TRACE("> Simulation_AddMaterial(S=%p) [omega=%f]\n", S, S->omega[0]);
	if(NULL == S){
		S4_TRACE("< Simulation_AddMaterial (failed; S == NULL) [omega=%f]\n", S->omega[0]);
		return NULL;
	}
	Material *M;
	if(NULL == S->material){
		S->material = (Material *)S4_malloc(sizeof(Material));
		S->material->next = NULL;
		M = S->material;
	}else{
		M = S->material;
		while(NULL != M->next){
			M = M->next;
		}
		M->next = (Material *)S4_malloc(sizeof(Material));
		M = M->next;
		M->next = NULL;
	}
	Material_Init(M, NULL, NULL);
	S4_TRACE("< Simulation_AddMaterial (returning M=%p) [omega=%f]\n", M, S->omega[0]);
	return M;
}

Layer* Simulation_AddLayer(Simulation *S){
	S4_TRACE("> Simulation_AddLayer(S=%p) [omega=%f]\n", S, S->omega[0]);
	if(NULL == S){
		S4_TRACE("< Simulation_AddLayer (failed; S == NULL) [omega=%f]\n", S->omega[0]);
		return NULL;
	}
	Simulation_DestroySolution(S);
	Layer *L;
	if(NULL == S->layer){
		S->layer = (Layer *)S4_malloc(sizeof(Layer));
		S->layer->next = NULL;
		L = S->layer;
	}else{
		L = S->layer;
		while(NULL != L->next){
			L = L->next;
		}
		L->next = (Layer *)S4_malloc(sizeof(Layer));
		L = L->next;
		L->next = NULL;
	}
	Layer_Init(L, NULL, 0, NULL, NULL);
	S4_TRACE("< Simulation_AddLayer (returning L=%p) [omega=%f]\n", L, S->omega[0]);
	return L;
}
int Simulation_RemoveLayerPatterns(Simulation *S, Layer *layer){
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

	S4_TRACE("< Simulation_RemoveLayerPatterns [omega=%f]\n", S->omega[0]);
	return 0;
}
int Simulation_ChangeLayerThickness(Simulation *S, Layer *layer, const double *thick){
	S4_TRACE("> Simulation_ChangeLayerThickness(S=%p, layer=%p, thick=%g) [omega=%f]\n",
		S, layer, thick, S->omega[0]);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(thick < 0){ ret = -3; }
	if(0 != ret){
		S4_TRACE("< Simulation_RemoveLayerPatterns (failed; ret = %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	layer->thickness = *thick;
	Simulation_DestroyLayerSolutions(S);

	S4_TRACE("< Simulation_ChangeLayerThickness [omega=%f]\n", S->omega[0]);
	return 0;
}
int Simulation_SetNumG(Simulation *S, int n){
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

	S4_TRACE("< Simulation_SetNumG [omega=%f]\n", S->omega[0]);
	return 0;
}
int Simulation_GetNumG(const Simulation *S, int **G){
	int ret = 0;
	S4_TRACE("> Simulation_GetNumG(S=%p, G=%p) [omega=%f]\n", S, G, S->omega[0]);
	if(NULL == S){ ret = -1; }
	if(0 != ret){
		S4_TRACE("< Simulation_GetNumG (failed; ret = %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	ret = S->n_G;

	if(NULL != G){
		if(NULL != S->solution){
			*G = S->solution->G;
		}
	}

	S4_TRACE("< Simulation_GetNumG [omega=%f]\n", S->omega[0]);
	return ret;
}
int Simulation_AddLayerPatternCircle(
	Simulation *S,
	Layer *layer, int material,
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
	Simulation *S,
	Layer *layer, int material,
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
	Simulation *S,
	Layer *layer, int material,
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
	Simulation *S,
	Layer *layer, int material,
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
int Simulation_InitSolution(Simulation *S){
	S4_TRACE("> Simulation_InitSolution(S=%p) [omega=%f]\n", S, S->omega[0]);
	int layer_count, i;
	Solution *sol;
	if(NULL == S){
		S4_TRACE("< Simulation_InitSolution (failed; S == NULL) [omega=%f]\n", S->omega[0]);
		return -1;
	}
	if(S->n_G < 1){
		S4_TRACE("< Simulation_InitSolution (failed; S->n_G < 1) [omega=%f]\n", S->omega[0]);
		return 9;
	}

	// Check that every material referenced exists // no, we are sure due to way we added patterns

	// Count layers
	// Check that layer copies and the excitation layer names exist
	Layer *L = S->layer;
	bool found_ex_layer = (NULL == S->exc.layer);
	layer_count = 0;
	while(NULL != L){
		++layer_count;
		if(NULL != L->copy){
			const Layer *L2 = S->layer;
			bool found = false;
			while(NULL != L2){
				if(0 == strcmp(L2->name, L->copy)){
					if(NULL != L2->copy){
						S4_TRACE("< Simulation_InitSolution (failed; layer %s is referenced by a copy but is also a copy) [omega=%f]\n", L2->name, S->omega[0]);
						return 11;
					}
					found = true;
					break;
				}
				L2 = L2->next;
			}
			if(!found){
				S4_TRACE("< Simulation_InitSolution (failed; could not find layer %s is referenced by a copy) [omega=%f]\n", L->copy, S->omega[0]);
				return 10;
			}
		}else{
			// check that no duplicate names exist
			const Layer *L2 = S->layer;
			while(L2 != L){
				if(0 == strcmp(L2->name, L->name)){
					S4_TRACE("< Simulation_InitSolution (failed; layer name %s appears more than once) [omega=%f]\n", L->name, S->omega[0]);
					return 12;
				}
				L2 = L2->next;
			}
		}
		if(!found_ex_layer && NULL != S->exc.layer && 0 == strcmp(S->exc.layer, L->name)){
			found_ex_layer = true;
		}
		// check that if we have a 1D pattern, the only shapes are rectangles
		if(0 == S->Lr[2] && 0 == S->Lr[3]){
			for(i = 0; i < L->pattern.nshapes; ++i){
				if(RECTANGLE != L->pattern.shapes[i].type){
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

		L = L->next;
	}
	if(layer_count < 1){
		S4_TRACE("< Simulation_InitSolution (failed; less than one layer found) [omega=%f]\n", S->omega[0]);
		return 14;
	}
	if(!found_ex_layer){
		S4_TRACE("< Simulation_InitSolution (failed; excitation layer not found) [omega=%f]\n", S->omega[0]);
		return 13;
	}

	if(NULL != S->solution){
		Simulation_DestroySolution(S);
	}
	S->solution = (Solution*)S4_malloc(sizeof(Solution));
	if(NULL == S->solution){
		S4_TRACE("< Simulation_InitSolution (failed; could not allocate S->solution) [omega=%f]\n", S->omega[0]);
		return 1;
	}
	sol = S->solution;
	sol->G = (int*)S4_malloc(sizeof(int)*2*S->n_G);
	if(NULL == sol->G){
		S4_TRACE("< Simulation_InitSolution (failed; could not allocate sol->G) [omega=%f]\n", S->omega[0]);
		return 1;
	}
	sol->layer_bands = (void**)S4_malloc(sizeof(void*)*2*layer_count);
	if(NULL == sol->layer_bands){
		S4_TRACE("< Simulation_InitSolution (failed; could not allocate sol->layer_bands) [omega=%f]\n", S->omega[0]);
		return 1;
	}
	sol->layer_solution = sol->layer_bands + layer_count;
	for(i = 0; i < 2*layer_count; ++i){
		sol->layer_bands[i] = NULL;
	}

	// Get G vectors
	if(0 != S->Lr[2] || 0 != S->Lr[3]){
		unsigned int NG = S->n_G;
		G_select(S->options.lattice_truncation, &NG, S->Lk, sol->G);
		S->n_G = NG;
	}else{
		// 1D lattice
		sol->G[0] = 0; sol->G[1] = 0;
		int remaining = (S->n_G-1)/2;
		S->n_G = 1+2*remaining;
		for(i = 0; i < remaining; ++i){
			sol->G[2+4*i+0] = i+1;
			sol->G[2+4*i+1] = 0;
			sol->G[2+4*i+2] = -(i+1);
			sol->G[2+4*i+3] = 0;
		}
	}
	S4_VERB(1, "Using %d G-vectors\n", S->n_G);
	S4_TRACE("I  Simulation_InitSolution G: (%d) [omega=%f]\n", S->n_G, S->omega[0]);

	sol->kx = (double*)S4_malloc(sizeof(double)*2*S->n_G);
	if(NULL == sol->kx){
		S4_TRACE("< Simulation_InitSolution (failed; could not allocate sol->kx) [omega=%f]\n", S->omega[0]);
		return 1;
	}
	sol->ky = sol->kx+S->n_G;

	for(i = 0; i < S->n_G; ++i){
		sol->kx[i] = S->k[0]*S->omega[0] + 2*M_PI*(S->Lk[0]*sol->G[2*i+0] + S->Lk[2]*sol->G[2*i+1]);
		sol->ky[i] = S->k[1]*S->omega[0] + 2*M_PI*(S->Lk[1]*sol->G[2*i+0] + S->Lk[3]*sol->G[2*i+1]);
	}

	S4_TRACE("< Simulation_InitSolution [omega=%f]\n", S->omega[0]);
	return 0;
}

// realistically, can return -2 if layer is not found, or an InitSolution code
int Simulation_GetLayerSolution(Simulation *S, Layer *layer, LayerBands **layer_bands, LayerSolution **layer_solution){
	S4_TRACE("> Simulation_GetLayerSolution(S=%p, layer=%p, layer_bands=%p (%p), layer_solution=%p (%p)) [omega=%f]\n",
		S, layer,
		layer_bands, (NULL != layer_bands ? *layer_bands : NULL), layer_solution, (NULL != layer_solution ? *layer_solution : NULL), S->omega[0]);
	int error;
	Solution *sol;
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

	LayerBands** Lbands = (LayerBands**)sol->layer_bands;
	LayerSolution** Lsoln = (LayerSolution**)sol->layer_solution;
	Layer *L = S->layer;
	while(NULL != L){
		if(L == layer){
			if(NULL == *Lsoln){
				error = Simulation_ComputeLayerSolution(S, L, layer_bands, layer_solution);
				if(0 != error){ // should never happen
					S4_TRACE("< Simulation_GetLayerSolution (failed; Simulation_ComputeLayerSolution returned %d) [omega=%f]\n", error, S->omega[0]);
					return 2;
				}
				if(NULL == L->copy){
					*Lbands = *layer_bands;
				}
				*Lsoln = *layer_solution;
			}else{
				if(NULL != L->copy){
					int refi = 0;
					Simulation_GetLayerByName(S, L->copy, &refi);
					*layer_bands = (LayerBands*)sol->layer_bands[refi];
				}else{
					*layer_bands = *Lbands;
				}
				*layer_solution = *Lsoln;
			}

			S4_TRACE("< Simulation_GetLayerSolution [omega=%f]\n", S->omega[0]);
			return 0;
		}
		++Lbands;
		++Lsoln;
		L = L->next;
	}

	S4_TRACE("< Simulation_GetLayerSolution (failed; layer not found) [omega=%f]\n", S->omega[0]);
	return -2;
}

int Simulation_SolveLayer(Simulation *S, Layer *layer){
	S4_TRACE("> Simulation_SolveLayer(S=%p, L=%p (%s)) [omega=%f]\n", S, layer, (NULL != layer && NULL != layer->name ? layer->name : ""), S->omega[0]);
	if(NULL == S){
		S4_TRACE("< Simulation_GetLayerSolution (failed; S == NULL) [omega=%f]\n", S->omega[0]);
		return -1;
	}
	if(NULL == layer){
		S4_TRACE("< Simulation_GetLayerSolution (failed; layer == NULL) [omega=%f]\n", S->omega[0]);
		return -2;
	}
	S4_CHECK{
		Layer *L = S->layer;
		size_t i = 0;
		while(NULL != L){
			if(L == layer){
				S4_TRACE("I  Simulation_SolveLayer: solve layer index %d [omega=%f]\n", i, S->omega[0]);
				break;
			}
			L = L->next;
			++i;
		}
	}

	LayerBands *Lbands = NULL;
	LayerSolution *Lsoln = NULL;
	int ret = Simulation_GetLayerSolution(S, layer, &Lbands, &Lsoln);

	S4_TRACE("< Simulation_SolveLayer [omega=%f]\n", S->omega[0]);
	return ret;
}

int Simulation_ComputeLayerSolution(Simulation *S, Layer *L, LayerBands **layer_bands, LayerSolution **layer_solution){
	S4_TRACE("> Simulation_ComputeLayerSolution(S=%p, L=%p (%s), layer_bands=%p (%p), LayerSolution=%p (%p)) [omega=%f]\n",
		S, L, (NULL != L && NULL != L->name ? L->name : ""), layer_bands, (NULL != layer_bands ? *layer_bands : NULL), layer_solution, (NULL != layer_solution ? *layer_solution : NULL), S->omega[0]);
	if(NULL != layer_bands && NULL == layer_solution){
		// only need to compute bands for L
		Simulation_ComputeLayerBands(S, L, layer_bands);

		S4_TRACE("< Simulation_ComputeLayerSolution [omega=%f]\n", S->omega[0]);
		return 0;
	}
	if(NULL == layer_solution){
		S4_TRACE("< Simulation_ComputeLayerSolution [omega=%f]\n", S->omega[0]);
		return 0;
	}
	// At this point, layer_solution != NULL, we need to get all bands

	S4_VERB(1, "Computing solution in layer: %s\n", NULL != L->name ? L->name : "");

	const size_t n = S->n_G;
	const size_t n2 = 2*n;
	const size_t n4 = 2*n2;

	// compute all bands and then get solution
	Layer *SL = S->layer;
	Solution *sol = S->solution;
	LayerBands **Lbands = (LayerBands**)sol->layer_bands;
	LayerSolution **Lsoln = (LayerSolution**)sol->layer_solution;
	int layer_count = 0;
	int which_layer = 0;
	bool found_layer = false;
	while(NULL != SL){
		if(NULL == *Lbands && NULL == SL->copy){
			Simulation_ComputeLayerBands(S, SL, Lbands);
		}
		if(L == SL){
			found_layer = true;
			which_layer = layer_count;
			*Lsoln = (LayerSolution*)S4_malloc(sizeof(LayerSolution));
			if(NULL != SL->copy){
				int refi = 0;
				Simulation_GetLayerByName(S, SL->copy, &refi);
				*layer_bands = (LayerBands*)sol->layer_bands[refi];
			}else{
				*layer_bands = *Lbands;
			}
			*layer_solution = *Lsoln;
			(*layer_solution)->ab = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>)*n4);
		}
		++Lbands;
		++Lsoln;
		SL = SL->next;
		++layer_count;
	}
	if(!found_layer){
		S4_TRACE("< Simulation_ComputeLayerSolution (failed; could not find layer) [omega=%f]\n", S->omega[0]);
		return -1;
	}

	// Make arrays of q, kp, and phi
	double *lthick = (double*)S4_malloc(sizeof(double)*layer_count);
	int *lepstype = (int*)S4_malloc(sizeof(int)*layer_count);
	const std::complex<double> **lq   = (const std::complex<double> **)S4_malloc(sizeof(const std::complex<double> *)*layer_count*4);
	if(NULL == lthick || NULL == lq){
		S4_TRACE("< Simulation_ComputeLayerSolution (failed; could not allocate work arrays) [omega=%f]\n", S->omega[0]);
		return 1;
	}
	const std::complex<double> **lepsinv  = lq  + layer_count;
	const std::complex<double> **lkp  = lepsinv  + layer_count;
	const std::complex<double> **lphi = lkp + layer_count;

	int i;
	Lbands = (LayerBands**)sol->layer_bands;
	for(i = 0, SL = S->layer; NULL != SL; SL = SL->next, i++){
		int refi = i;
		lthick[i] = SL->thickness;
		if(NULL != SL->copy){
			Simulation_GetLayerByName(S, SL->copy, &refi);
		}
		lq  [i] = Lbands[refi]->q;
		lepsinv [i] = Lbands[refi]->Epsilon_inv;
		lepstype[i] = Lbands[refi]->epstype;
		lkp [i] = Lbands[refi]->kp;
		lphi[i] = Lbands[refi]->phi;
	}

	// Compose the RCWA solution
	int error = 0;
	if(0 == S->exc.type){
		// Front incidence by planewave
		size_t order = S->exc.sub.planewave.order;
		size_t phicopy_size = (NULL == lphi[0] ? 0 : n2*n2);
		std::complex<double> *a0 = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>)*(n2+phicopy_size));
		std::complex<double> *phicopy = (NULL == lphi[0] ? NULL : a0 + n2);
		RNP::TBLAS::Fill(n2, 0., a0,1);
		if(order < n){
			a0[order+0] = std::complex<double>(S->exc.sub.planewave.hx[0], S->exc.sub.planewave.hx[1]);
			a0[order+n] = std::complex<double>(S->exc.sub.planewave.hy[0], S->exc.sub.planewave.hy[1]);
		}
		// [ kp.phi.inv(q) -kp.phi.inv(q) ] [ a ] = [-ey;ex ]
		// [     phi            phi       ] [ b ]   [ hx;hy ]
		// We assume b = 0.
		// Then a = inv(phi)*[ hx;hy ]
		if(NULL != lphi[0]){
			RNP::TBLAS::CopyMatrix<'A'>(n2,n2, lphi[0],n2, phicopy,n2);
			RNP::LinearSolve<'N'>(n2,1, phicopy,n2, a0,n2, NULL, NULL);
		}

		S4_TRACE("I  Calling SolveInterior(layer_count=%d, which_layer=%d, n=%d, lthick,lq,lkp,lphi={\n", layer_count, which_layer, S->n_G);
		for(i = 0; i < layer_count; ++i){
			S4_TRACE("I    %f, %p (0,0=%f,%f), %p (0,0=%f,%f), %p (0,0=%f,%f)\n", lthick[i],
				lq[i], lq[i][0].real(), lq[i][0].imag(),
				lkp[i], NULL != lkp[i] ? lkp[i][0].real() : 0., NULL != lkp[i] ? lkp[i][0].imag() : 0.,
				lphi[i], NULL != lphi[0] ? lphi[i][0].real() : 1., NULL != lphi[0] ? lphi[i][0].imag() : 0.);
		}
		S4_TRACE("I   }, a0[0]=%f,%f, a0[n]=%f,%f, ...) [omega=%f]\n", a0[0].real(), a0[0].imag(), a0[S->n_G].real(), a0[S->n_G].imag(), S->omega[0]);

		error = SolveInterior(
			layer_count, which_layer,
			S->n_G,
			S->solution->kx, S->solution->ky,
			std::complex<double>(S->omega[0], S->omega[1]),
			lthick, lq, lepsinv, lepstype, lkp, lphi,
			a0, // length 2*n
			NULL, // bN
			(*layer_solution)->ab);
		S4_free(a0);
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

		S4_TRACE("I  Calling SolveInterior(layer_count=%d, which_layer=%d, n=%d, lthick,lq,lkp,lphi={\n", layer_count, which_layer, S->n_G);
		for(size_t i = 0; i < layer_count; ++i){
			S4_TRACE("I    %f, %p (0,0=%f,%f), %p (0,0=%f,%f), %p (0,0=%f,%f)\n", lthick[i],
				lq[i], lq[i][0].real(), lq[i][0].imag(),
				lkp[i], NULL != lkp[i] ? lkp[i][0].real() : 0., NULL != lkp[i] ? lkp[i][0].imag() : 0.,
				lphi[i], NULL != lphi[0] ? lphi[i][0].real() : 1., NULL != lphi[0] ? lphi[i][0].imag() : 0.);
		}
		S4_TRACE("I   }, a0[0]=%f,%f, a0[n]=%f,%f, ...) [omega=%f]\n", a0[0].real(), a0[0].imag(), a0[S->n_G].real(), a0[S->n_G].imag(), S->omega[0]);

		error = SolveInterior(
			layer_count, which_layer,
			S->n_G,
			S->solution->kx, S->solution->ky,
			std::complex<double>(S->omega[0], S->omega[1]),
			lthick, lq, lepsinv, lepstype, lkp, lphi,
			a0, // length 2*n
			bN, // bN
			(*layer_solution)->ab);
		S4_free(a0);
	}else if(1 == S->exc.type){
		Layer *l[2];
		int li;
		l[0] = Simulation_GetLayerByName(S, S->exc.layer, &li);
		if(NULL == l[0]){ return 13; }
		l[1] = l[0]->next;
		if(NULL == l[1]){ return 13; }
		std::complex<double> *ab = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>)*(n4+n4*n4+n2*n2));
		std::complex<double> *work4 = ab + n4;
		std::complex<double> *work2 = work4 + n4*n4;
		std::complex<double> J0[3] = {
			std::complex<double>(S->exc.sub.dipole.moment[0],S->exc.sub.dipole.moment[1]),
			std::complex<double>(S->exc.sub.dipole.moment[2],S->exc.sub.dipole.moment[3]),
			std::complex<double>(S->exc.sub.dipole.moment[4],S->exc.sub.dipole.moment[5])
		};
		for(i = 0; i < n; ++i){
			const double phaseangle = -(S->solution->kx[i] * S->exc.sub.dipole.pos[0] + S->solution->ky[i] * S->exc.sub.dipole.pos[1]);
			const std::complex<double> phase(cos(phaseangle), sin(phaseangle));
			ab[0*n+i] = J0[2]*phase; // -ky eta jz
			//ab[1*n+i] = ; //  kx eta jz
			ab[2*n+i] =  J0[1]*phase; //   jy
			ab[3*n+i] = -J0[0]*phase; //  -jx
		}
		// Make eta*jz
		RNP::TBLAS::MultMV<'N'>(n,n, 1.,Lbands[li]->Epsilon_inv,n, &ab[0*n],1, 0.,&ab[1*n],1);
		RNP::TBLAS::Copy(n, &ab[1*n],1, &ab[0*n],1);
		// finish ab[0*n] and ab[1*n]
		for(i = 0; i < n; ++i){
			ab[0*n+i] *= -S->solution->ky[i];
			ab[1*n+i] *=  S->solution->kx[i];
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
		for(i = 0; i < n2; ++i){ // first scale work2 to make -f_l*S12(0,l)
			std::complex<double> f = -std::exp(lq[li][i] * std::complex<double>(0,lthick[li]));
			RNP::TBLAS::Scale(n2, f, &work2[i+0*n2],n2);
		}
		RNP::TBLAS::CopyMatrix<'A'>(n2,n2, work2,n2, &work4[0+n2*n4],n4);
		for(i = 0; i < n2; ++i){
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
			S->solution->kx, S->solution->ky,
			lepsinv[li], lepstype[li], lkp[li],
			n2, &work4[n2+n2*n4],n4,
			&work4[0+n2*n4],n4
		); // upper right complete

		// For lower right, -f_l S12(0,li) is still in work2
		for(i = 0; i < n2; ++i){
			work2[i+i*n4] -= 1.;
		}
		if(NULL != lphi[li]){ // use lower right as buffer
			RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,lphi[li],n2, work2,n2, 0.,&work4[n2+n2*n4],n4);
		}else{
			RNP::TBLAS::CopyMatrix<'A'>(n2,n2, work2,n2, &work4[n2+n2*n4],n4);
		} // lower right complete

		// For upper left, make f_l+1 S21(l+1,N) in lower left first
		for(i = 0; i < n2; ++i){ // make f_l+1*S12(l+1,N)
			std::complex<double> f = std::exp(lq[li+1][i] * std::complex<double>(0,lthick[li+1]));
			RNP::TBLAS::Scale(n2, f, &work4[n2+i+0*n2],n4);
		}
		RNP::TBLAS::CopyMatrix<'A'>(n2,n2, &work4[n2+0*n2],n4, &work4[0+0*n4],n4); // copy to upper left
		for(i = 0; i < n2; ++i){
			work4[0+i+(0+i)*n4] -= 1.;
		}
		for(i = 0; i < n2; ++i){
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
			S->solution->kx, S->solution->ky,
			lepsinv[li], lepstype[li], lkp[li],
			n2, work2,n2,
			&work4[0+0*n4],n4
		); // upper left complete

		// Lower left
		for(i = 0; i < n2; ++i){
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
				S->solution->kx, S->solution->ky,
				std::complex<double>(S->omega[0], S->omega[1]),
				lthick, lq, lepsinv, lepstype, lkp, lphi,
				NULL, // length 2*n
				&ab[n2], // bN
				(*layer_solution)->ab);
		}else{
			error = SolveInterior(
				layer_count-li, which_layer-li,
				S->n_G,
				S->solution->kx, S->solution->ky,
				std::complex<double>(S->omega[0], S->omega[1]),
				lthick+li, lq+li, lepsinv+li, lepstype+li, lkp+li, lphi+li,
				&ab[0], // length 2*n
				NULL, // bN
				(*layer_solution)->ab);
		}
		S4_free(ab);
	}

	S4_TRACE("I  ab[0] = %f,%f [omega=%f]\n", (*layer_solution)->ab[0].real(), (*layer_solution)->ab[0].imag(), S->omega[0]);

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

double Simulation_GetUnitCellSize(const Simulation *S){
	if(0 == S->Lr[2] && 0 == S->Lr[3]){
		return hypot(S->Lr[0], S->Lr[1]);
	}else{
		return fabs(S->Lr[0]*S->Lr[3] - S->Lr[1]*S->Lr[2]);
	}
}

int Simulation_ComputeLayerBands(Simulation *S, Layer *L, LayerBands **bands){
	S4_TRACE("> Simulation_ComputeLayerBands(S=%p, L=%p (%s), bands=%p (%p)) [omega=%f]\n", S, L, (NULL != L && NULL != L->name ? L->name : ""), bands, (NULL != bands ? *bands : NULL), S->omega[0]);
	if(NULL == S){
		S4_TRACE("< Simulation_ComputeLayerBands (failed; S == NULL) [omega=%f]\n", S->omega[0]);
		return -1;
	}
	if(NULL == L){
		S4_TRACE("< Simulation_ComputeLayerBands (failed; L == NULL) [omega=%f]\n", S->omega[0]);
		return -2;
	}
	if(NULL == bands){
		S4_TRACE("< Simulation_ComputeLayerBands (early exit; bands == NULL) [omega=%f]\n", S->omega[0]);
		return 0;
	}

	S4_VERB(1, "Computing bands of layer: %s\n", NULL != L->name ? L->name : "");

	*bands = (LayerBands*)S4_malloc(sizeof(LayerBands));
	LayerBands *pB = *bands;
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
		const Material *M;
		if(NULL == L->copy){
			M = Simulation_GetMaterialByName(S, L->material, NULL);
		}else{
			M = Simulation_GetMaterialByName(S, Simulation_GetLayerByName(S, L->copy, NULL)->material, NULL);
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
		const Material *M;
		if(NULL == L->copy){
			//eps_scalar = Simulation_GetEpsilonByName(S, L->material);
			M = Simulation_GetMaterialByName(S, L->material, NULL);
		}else{
			//eps_scalar = Simulation_GetEpsilonByName(S, Simulation_GetLayerByName(S, L->copy, NULL)->material);
			M = Simulation_GetMaterialByName(S, Simulation_GetLayerByName(S, L->copy, NULL)->material, NULL);
		}
		if(0 == M->type){
			std::complex<double> eps_scalar(M->eps.s[0], M->eps.s[1]);
			RNP::TBLAS::SetMatrix<'A'>(n,n,0.,1./eps_scalar,pB->Epsilon_inv, n);
			RNP::TBLAS::SetMatrix<'A'>(n2,n2,0.,eps_scalar,pB->Epsilon2, n2);
			SolveLayerEigensystem_uniform(
				std::complex<double>(S->omega[0],S->omega[1]), n, S->solution->kx, S->solution->ky,
				eps_scalar, pB->q, pB->kp, pB->phi);
		}else{
			RNP::TBLAS::SetMatrix<'A'>(n,n,0.,1./std::complex<double>(M->eps.abcde[8],M->eps.abcde[9]),pB->Epsilon_inv, n);
			RNP::TBLAS::SetMatrix<'A'>(n,n,0.,std::complex<double>(M->eps.abcde[0],M->eps.abcde[1]),&pB->Epsilon2[0+0*n2], n2);
			RNP::TBLAS::SetMatrix<'A'>(n,n,0.,std::complex<double>(M->eps.abcde[4],M->eps.abcde[5]),&pB->Epsilon2[n+0*n2], n2);
			RNP::TBLAS::SetMatrix<'A'>(n,n,0.,std::complex<double>(M->eps.abcde[2],M->eps.abcde[3]),&pB->Epsilon2[0+n*n2], n2);
			RNP::TBLAS::SetMatrix<'A'>(n,n,0.,std::complex<double>(M->eps.abcde[6],M->eps.abcde[7]),&pB->Epsilon2[n+n*n2], n2);

			S4_VERB(1, "Solving eigensystem of layer: %s\n", NULL != L->name ? L->name : "");
			SolveLayerEigensystem(
				std::complex<double>(S->omega[0],S->omega[1]), n, S->solution->kx, S->solution->ky,
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
				S->solution->kx, S->solution->ky,
				pB->Epsilon_inv, pB->Epsilon2, pB->epstype,
				pB->q, pB->kp, pB->phi,
				&dum, rwork, lwork
			);
			lwork = (int)(dum.real() + 0.5);
			work = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>) * lwork);
			SolveLayerEigensystem(
				std::complex<double>(S->omega[0],S->omega[1]), n,
				S->solution->kx, S->solution->ky,
				pB->Epsilon_inv, pB->Epsilon2, pB->epstype,
				pB->q, pB->kp, pB->phi,
				work, rwork, lwork
			);
			S4_free(work);
			S4_free(rwork);
		}
	}
	S4_TRACE("I  q[0] = %f,%f [omega=%f]\n", pB->q[0].real(), pB->q[0].imag(), S->omega[0]);

	S4_TRACE("< Simulation_ComputeLayerBands [omega=%f]\n", S->omega[0]);
	return 0;
}

// realistically, can only return an InitSolution code
int Simulation_GetPoyntingFlux(Simulation *S, Layer *layer, double offset, double powers[4]){
	S4_TRACE("> Simulation_GetPoyntingFlux(S=%p, layer=%p, offset=%f, powers=%p) [omega=%f]\n",
		S, layer, offset, powers, S->omega[0]);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(NULL == powers){ ret = -4; }
	if(0 != ret){
		S4_TRACE("< Simulation_GetPoyntingFlux (failed; ret = %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	LayerBands *Lbands;
	LayerSolution *Lsoln;

	ret = Simulation_GetLayerSolution(S, layer, &Lbands, &Lsoln);
	if(0 != ret){
		S4_TRACE("< Simulation_GetPoyntingFlux (failed; Simulation_GetLayerSolution returned %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	const int n = S->n_G;
	const int n2 = 2*n;
	const int n4 = 2*n2;

	std::complex<double> *ab = (std::complex<double> *)S4_malloc(sizeof(std::complex<double>) * (n4+4*n2));
	if(NULL == ab){
		S4_TRACE("< Simulation_GetPoyntingFlux (failed; allocation failed) [omega=%f]\n", S->omega[0]);
		return 1;
	}
	std::complex<double> *work = ab + n4;

	memcpy(ab, Lsoln->ab, sizeof(std::complex<double>) * n4);
	TranslateAmplitudes(n, Lbands->q, layer->thickness, offset, ab);

	std::complex<double> forw, back;
	GetZPoyntingFlux(n, S->solution->kx, S->solution->ky, std::complex<double>(S->omega[0],S->omega[1]), Lbands->q, Lbands->Epsilon_inv, Lbands->epstype, Lbands->kp, Lbands->phi, ab, &forw, &back, work);
	powers[0] = forw.real();
	powers[1] = back.real();
	powers[2] = forw.imag();
	powers[3] = back.imag();

	S4_free(ab);
	S4_TRACE("< Simulation_GetPoyntingFlux returning %f, %f, %f, %f [omega=%f]\n", powers[0], powers[1], powers[2], powers[3], S->omega[0]);
	return 0;
}

int Simulation_GetPoyntingFluxByG(Simulation *S, Layer *layer, double offset, double *powers){
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

	LayerBands *Lbands;
	LayerSolution *Lsoln;

	ret = Simulation_GetLayerSolution(S, layer, &Lbands, &Lsoln);
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

	memcpy(ab, Lsoln->ab, sizeof(std::complex<double>) * n4);
	TranslateAmplitudes(n, Lbands->q, layer->thickness, offset, ab);

	GetZPoyntingFluxComponents(n, S->solution->kx, S->solution->ky, std::complex<double>(S->omega[0],S->omega[1]), Lbands->q, Lbands->Epsilon_inv, Lbands->epstype, Lbands->kp, Lbands->phi, ab, forw, back, work);
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

int Simulation_GetPropagationConstants(Simulation *S, Layer *L, double *q){
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

	Solution *sol;
	if(NULL == S->solution){
		int error = Simulation_InitSolution(S);
		if(0 != error){
			S4_TRACE("< Simulation_GetPropagationConstants (failed; Simulation_InitSolution returned %d) [omega=%f]\n", error, S->omega[0]);
			return error;
		}
	}
	sol = S->solution;

	// compute all bands and then get solution
	Layer *SL = S->layer;
	LayerBands **Lbands = (LayerBands**)sol->layer_bands;
	LayerBands *layer_bands;
	bool found_layer = false;
	while(NULL != SL){
		if(NULL == *Lbands && NULL == SL->copy){
			Simulation_ComputeLayerBands(S, SL, Lbands);
		}
		if(L == SL){
			found_layer = true;
			if(NULL != SL->copy){
				int refi = 0;
				Simulation_GetLayerByName(S, SL->copy, &refi);
				layer_bands = (LayerBands*)sol->layer_bands[refi];
			}else{
				layer_bands = *Lbands;
			}
		}
		++Lbands;
		SL = SL->next;
	}
	if(!found_layer){
		S4_TRACE("< Simulation_GetPropagationConstants (failed; could not find layer) [omega=%f]\n", S->omega[0]);
		return -1;
	}

	for(int i = 0; i < n; ++i){
		q[2*i+0] = layer_bands->q[i].real();
		q[2*i+1] = layer_bands->q[i].imag();
	}

	S4_TRACE("< Simulation_GetPropagationConstants [omega=%f]\n", S->omega[0]);
	return ret;
}

int Simulation_GetAmplitudes(Simulation *S, Layer *layer, double offset, double *forw, double *back){
	S4_TRACE("> Simulation_GetAmplitudes(S=%p, layer=%p, offset=%f, forw=%p, back=%p) [omega=%f]\n",
		S, layer, offset, forw, back, S->omega[0]);
	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(0 != ret){
		S4_TRACE("< Simulation_GetAmplitudes (failed; ret = %d) [omega=%f]\n", ret, S->omega[0]);
		return ret;
	}

	LayerBands *Lbands;
	LayerSolution *Lsoln;

	ret = Simulation_GetLayerSolution(S, layer, &Lbands, &Lsoln);
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

	memcpy(ab, Lsoln->ab, sizeof(std::complex<double>) * n4);
	TranslateAmplitudes(n, Lbands->q, layer->thickness, offset, ab);

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

void Simulation_DestroySolution(Simulation *S){
	S4_TRACE("> Simulation_DestroySolution(S=%p) [omega=%f]\n", S, S->omega[0]);
	Solution *sol;
	if(NULL == S){
		S4_TRACE("< Simulation_DestroySolution (early exit; S == NULL) [omega=%f]\n", S->omega[0]);
		return;
	}
	sol = S->solution;
	if(NULL == sol){
		S4_TRACE("< Simulation_DestroySolution (early exit; sol == NULL) [omega=%f]\n", S->omega[0]);
		return;
	}
	if(NULL != sol->G){ S4_free(sol->G); sol->G = NULL; }
	if(NULL != sol->kx){ S4_free(sol->kx); sol->kx = NULL; }
	if(NULL != sol->layer_bands){
		void **Lbands = sol->layer_bands;
		void **Lsoln = sol->layer_solution;
		Layer *L = S->layer;
		while(NULL != L){
			// release Lbands
			if(NULL != Lbands){
				LayerBands *pB = (LayerBands*)(*Lbands);
				if(NULL != pB){
					if(NULL != pB->q){
						S4_free(pB->q);
					}
					S4_free(pB);
				}
			}
			// release Lsoln
			if(NULL != Lsoln){
				LayerSolution *pS = (LayerSolution*)(*Lsoln);
				if(NULL != pS){
					if(NULL != pS->ab){
						S4_free(pS->ab);
					}
					S4_free(pS);
				}
			}
			++Lbands;
			++Lsoln;
			L = L->next;
		}
		S4_free(sol->layer_bands);
		sol->layer_bands = NULL;
		sol->layer_solution = NULL;
	}
	S4_free(S->solution); S->solution = NULL;

	S4_TRACE("< Simulation_DestroySolution [omega=%f]\n", S->omega[0]);
}

void Simulation_DestroyLayerSolutions(Simulation *S){
	S4_TRACE("> Simulation_DestroyLayerSolutions(S=%p) [omega=%f]\n", S, S->omega[0]);
	Solution *sol;
	if(NULL == S){
		S4_TRACE("< Simulation_DestroyLayerSolutions (early exit; S == NULL) [omega=%f]\n", S->omega[0]);
		return;
	}
	sol = S->solution;
	if(NULL == sol){
		S4_TRACE("< Simulation_DestroyLayerSolutions (early exit; sol == NULL) [omega=%f]\n", S->omega[0]);
		return;
	}

	if(NULL != sol->layer_bands){
		void **Lsoln = sol->layer_solution;
		Layer *L = S->layer;
		while(NULL != L){
			// release Lsoln
			if(NULL != Lsoln){
				LayerSolution *pS = (LayerSolution*)(*Lsoln);
				if(NULL != pS){
					if(NULL != pS->ab){
						S4_free(pS->ab);
					}
					S4_free(pS);
				}
			}
			*Lsoln = NULL;
			++Lsoln;
			L = L->next;
		}
	}

	S4_TRACE("< Simulation_DestroyLayerSolutions [omega=%f]\n", S->omega[0]);
}

Material* Simulation_GetMaterialByName(const Simulation *S, const char *name, int *index){
	S4_TRACE("> Simulation_GetMaterialByName(S=%p, name=%p (%s)) [omega=%f]\n",
		S, name, (NULL == name ? "" : name), S->omega[0]);
	Material *M = S->material;
	int i = 0;
	while(NULL != M){
		if(0 == strcmp(M->name, name)){
			if(NULL != index){ *index = i; }
			S4_TRACE("< Simulation_GetMaterialByName returning %s [omega=%f]\n", M->name, S->omega[0]);
			return M;
		}
		M = M->next;
		++i;
	}
	S4_TRACE("< Simulation_GetMaterialByName (failed; material name not found) [omega=%f]\n", S->omega[0]);
	return NULL;
}

Material* Simulation_GetMaterialByIndex(const Simulation *S, int i){
	S4_TRACE("> Simulation_GetMaterialByIndex(S=%p, i=%d) [omega=%f]\n", S, i, S->omega[0]);
	Material *M = S->material;
	while(NULL != M){
		if(0 == i){
			S4_TRACE("< Simulation_GetMaterialByIndex returning %s [omega=%f]\n", M->name, S->omega[0]);
			return M;
		}
		--i;
		M = M->next;
	}
	S4_TRACE("< Simulation_GetMaterialByIndex (failed; index out of bounds) [omega=%f]\n", S->omega[0]);
	return NULL;
}
Layer* Simulation_GetLayerByName(const Simulation *S, const char *name, int *index){
	S4_TRACE("> Simulation_GetLayerByName(S=%p, name=%p (%s), index=%p) [omega=%f]\n",
		S, name, (NULL == name ? "" : name), index, S->omega[0]);
	if(NULL == S || NULL == name){ return NULL; }
	Layer *L = S->layer;
	int i = 0;
	while(NULL != L){
		if(0 == strcmp(L->name, name)){
			if(NULL != index){ *index = i; }
			S4_TRACE("< Simulation_GetLayerByName returning %p [omega=%f]\n", L, S->omega[0]);
			return L;
		}
		L = L->next;
		++i;
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

int Simulation_OutputStructurePOVRay(Simulation *S, FILE *fp){
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
	Material *M = S->material;
	while(NULL != M){
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
		M = M->next;
	}

	Layer *L = S->layer;
	unsigned int layer_name_counter = 0;
	double layer_offset = 0;
	int num_layers = 0;
	while(NULL != L){
		const char *name = L->name;
		if(NULL == name){
			fprintf(f, "#declare Layer_%u = union{\n", layer_name_counter++);
		}else{
			fprintf(f, "#declare Layer_%s = union{\n", name);
		}

		bool comment_out = false;
		if(NULL != L->material){
			if(0 == strcasecmp("air", L->material) || 0 == strcasecmp("vacuum", L->material)){
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

		if(NULL != L->copy){
			Layer *Lcopy = Simulation_GetLayerByName(S, L->copy, NULL);
			if(NULL != Lcopy){
				fprintf(f, "\t\ttexture { Material_%s }\n", Lcopy->material);
			}else{
				fprintf(f, "\t\ttexture { Material_ }\n");
			}
		}else{
			fprintf(f, "\t\ttexture { Material_%s }\n", L->material);
		}
		fprintf(f, "\t}\n");
		if(comment_out){
			fprintf(f, "*/\n");
		}

		for(int i = 0; i < L->pattern.nshapes; ++i){
			comment_out = false;

			Material *m = Simulation_GetMaterialByIndex(S, L->pattern.shapes[i].tag);
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
		num_layers++;
		L = L->next;
	}

	// Output postamble
	fprintf(f,
		"#declare Layers = union {\n"
	);
	L = S->layer;
	layer_name_counter = 0;
	int layer_counter = 0;
	while(NULL != L){
		const char *name = L->name;
		const char *comment = (0 < layer_counter && layer_counter+1 < num_layers) ? "" : "//";
		if(NULL == name){
			fprintf(f, "\t%sobject{ Layer_%u }\n", comment, layer_name_counter++);
		}else{
			fprintf(f, "\t%sobject{ Layer_%s }\n", comment, name);
		}
		layer_counter++;
		L = L->next;
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

int Simulation_OutputLayerPatternDescription(Simulation *S, Layer *layer, FILE *fp){
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

int Simulation_OutputLayerPatternRealization(Simulation *S, Layer *layer, int nx, int ny, FILE *fp){
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
		const Material *M;
		if(-1 == i){
			M = Simulation_GetMaterialByName(S, layer->material, NULL);
		}else{
			M = Simulation_GetMaterialByIndex(S, layer->pattern.shapes[i].tag);
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
					S->solution->G[2*g+0] * S->Lk[0] + S->solution->G[2*g+1] * S->Lk[2],
					S->solution->G[2*g+0] * S->Lk[1] + S->solution->G[2*g+1] * S->Lk[3]
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
int Simulation_GetField(Simulation *S, const double r[3], double fE[6], double fH[6]){
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

	Layer *L = S->layer;
	double dz = r[2];
	{
		double z = 0;
		while(NULL != L && r[2] > z+L->thickness){
			z += L->thickness;
			dz -= L->thickness;
			if(NULL == L->next){ break; }
			L = L->next;
		}
	}
	if(NULL == L){
		S4_TRACE("< Simulation_GetField (failed; no layers found)\n");
		return 14;
	}
//fprintf(stderr, "(%f,%f,%f) in %s: dz = %f\n", r[0], r[1], r[2], L->name, dz);

	LayerBands *Lbands;
	LayerSolution *Lsoln;
	int ret = Simulation_GetLayerSolution(S, L, &Lbands, &Lsoln);
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

	RNP::TBLAS::Copy(n4, Lsoln->ab,1, ab,1);
	//RNP::IO::PrintVector(n4, ab, 1);
	TranslateAmplitudes(S->n_G, Lbands->q, L->thickness, dz, ab);
	std::complex<double> efield[3], hfield[3];
	GetFieldAtPoint(
		S->n_G, S->solution->kx, S->solution->ky, std::complex<double>(S->omega[0],S->omega[1]),
		Lbands->q, Lbands->kp, Lbands->phi, Lbands->Epsilon_inv, Lbands->epstype,
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
int Simulation_GetFieldPlane(Simulation *S, int nxy[2], double zz, double *E, double *H){
	S4_TRACE("> Simulation_GetFieldPlane(S=%p, nxy=%p (%d,%d), z=%f, E=%p, H=%p)\n",
		S, (NULL == nxy ? 0 : nxy[0]), (NULL == nxy ? 0 : nxy[1]), zz, E, H);
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

	Layer *L = S->layer;
	double dz = zz;
	{
		double z = 0;
		while(NULL != L && zz > z+L->thickness){
			z += L->thickness;
			dz -= L->thickness;
			if(NULL == L->next){ break; }
			L = L->next;
		}
	}
	if(NULL == L){
		S4_TRACE("< Simulation_GetFieldPlane (failed; no layers found)\n");
		return 14;
	}
//fprintf(stderr, "(%f,%f,%f) in %s: dz = %f\n", r[0], r[1], r[2], L->name, dz);

	LayerBands *Lbands;
	LayerSolution *Lsoln;
	int ret = Simulation_GetLayerSolution(S, L, &Lbands, &Lsoln);
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

	RNP::TBLAS::Copy(n4, Lsoln->ab,1, ab,1);
	//RNP::IO::PrintVector(n4, ab, 1);
	TranslateAmplitudes(S->n_G, Lbands->q, L->thickness, dz, ab);
	size_t snxy[2] = { nxy[0], nxy[1] };
	GetFieldOnGrid(
		S->n_G, S->solution->G, S->solution->kx, S->solution->ky, std::complex<double>(S->omega[0],S->omega[1]),
		Lbands->q, Lbands->kp, Lbands->phi, Lbands->Epsilon_inv, Lbands->epstype,
		ab, snxy,
		reinterpret_cast<std::complex<double>*>(E),
		reinterpret_cast<std::complex<double>*>(H)
	);
	S4_free(ab);

	S4_TRACE("< Simulation_GetFieldPlane\n");
	return 0;
}

int Simulation_GetEpsilon(Simulation *S, const double r[3], double eps[2]){
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

	Layer *L = S->layer;
	{
		double z = 0;
		while(NULL != L && r[2] > z+L->thickness){
			z += L->thickness;
			if(NULL == L->next){ break; }
			L = L->next;
		}
	}
	if(NULL == L){
		S4_TRACE("< Simulation_GetField (failed; no layers found)\n");
		return 14;
	}

	if(NULL != L->copy){ L = Simulation_GetLayerByName(S, L->copy, NULL); }

	double *values = (double*)S4_malloc(sizeof(double)*2*(L->pattern.nshapes+1));
	for(int i = -1; i < L->pattern.nshapes; ++i){
		const Material *M;
		if(-1 == i){
			M = Simulation_GetMaterialByName(S, L->material, NULL);
		}else{
			M = Simulation_GetMaterialByIndex(S, L->pattern.shapes[i].tag);
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
			S->solution->G[2*g+0] * S->Lk[0] + S->solution->G[2*g+1] * S->Lk[2],
			S->solution->G[2*g+0] * S->Lk[1] + S->solution->G[2*g+1] * S->Lk[3]
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

int Simulation_GetSMatrixDeterminant(Simulation *S, double rmant[2], double *rbase, int *rexpo){
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
int Simulation_GetStressTensorIntegral(Simulation *S, Layer *layer, double offset, double Tint[6]){
	S4_TRACE("> Simulation_GetStressTensorIntegral(S=%p, layer=%p, offset=%f, Tint=%p)\n", S, layer, offset, Tint);

	int ret = 0;
	if(NULL == S){ ret = -1; }
	if(NULL == layer){ ret = -2; }
	if(NULL == Tint){ ret = -4; }
	if(0 != ret){
		S4_TRACE("< Simulation_GetStressTensorIntegral (failed; ret = %d)\n", ret);
		return ret;
	}

	LayerBands *Lbands;
	LayerSolution *Lsoln;

	ret = Simulation_GetLayerSolution(S, layer, &Lbands, &Lsoln);
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

	memcpy(ab, Lsoln->ab, sizeof(std::complex<double>) * n4);
	TranslateAmplitudes(n, Lbands->q, layer->thickness, offset, ab);

	std::complex<double> integral[3];
	GetZStressTensorIntegral(
		n, S->solution->kx, S->solution->ky,
		std::complex<double>(S->omega[0],S->omega[1]),
		Lbands->q, Lbands->kp, Lbands->phi, Lbands->Epsilon_inv, Lbands->Epsilon2, Lbands->epstype, ab, integral, work);
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
int Simulation_GetLayerVolumeIntegral(Simulation *S, Layer *layer, char which, double integral[2]){
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

	LayerBands *Lbands;
	LayerSolution *Lsoln;

	ret = Simulation_GetLayerSolution(S, layer, &Lbands, &Lsoln);
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
		n, S->solution->kx, S->solution->ky,
		std::complex<double>(S->omega[0],S->omega[1]),
		layer->thickness, Lbands->q, Lbands->kp, Lbands->phi, Lbands->Epsilon_inv, Lbands->Epsilon2, Lbands->epstype, Lsoln->ab, &zintegral, work);

	S4_free(work);

	integral[0] = zintegral.real();
	integral[1] = zintegral.imag();

	S4_TRACE("< Simulation_GetLayerVolumeIntegral\n");
	return 0;
}

// Returns a solution error code
int Simulation_GetLayerZIntegral(Simulation *S, Layer *layer, const double r[2], double integral[6]){
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

	LayerBands *Lbands;
	LayerSolution *Lsoln;

	ret = Simulation_GetLayerSolution(S, layer, &Lbands, &Lsoln);
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
		n, S->solution->kx, S->solution->ky,
		std::complex<double>(S->omega[0],S->omega[1]),
		layer->thickness, r, Lbands->q, Lbands->kp, Lbands->phi, Lbands->Epsilon_inv, Lbands->Epsilon2, Lbands->epstype, Lsoln->ab, integral, work);

	S4_free(work);


	S4_TRACE("< Simulation_GetLayerZIntegral\n");
	return 0;
}

int Simulation_MakeExcitationPlanewave(Simulation *S, const double angle[2], const double pol_s[2], const double pol_p[2], size_t order){
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

	const Material *M = Simulation_GetMaterialByName(S, S->layer->material, NULL);
	if(NULL == M){
		S4_TRACE("< Simulation_MakeExcitationPlanewave (failed; material %s not defined)\n", S->layer->material);
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
	/*
	S->ex[0] = ;
	S->ex[1] = ;
	S->ey[0] = ;
	S->ey[1] = ;
	*/
	/*
	std::cout << "k: " << S->k[0] << "\t" << S->k[1] << std::endl;
	std::cout << "hx: " << S->hx[0] << "\t" << S->hx[1] << std::endl;
	std::cout << "hy: " << S->hy[0] << "\t" << S->hy[1] << std::endl;
	//*/

	S4_TRACE("< Simulation_MakeExcitationPlanewave\n");
	return 0;
}

int Simulation_MakeExcitationDipole(Simulation *S, const double k[2], const char *layer, const double pos[2], const double moment[6]){
	S4_TRACE("> Simulation_MakeExcitationDipole(S=%p, k=%p (%f,%f), moment=%p (%f,%f,%f), ampphase=%p (%f,%f))\n", S,
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
	S->exc.layer = strdup(layer);

	S4_TRACE("< Simulation_MakeExcitationDipole\n");
	return 0;
}

int Simulation_MakeExcitationExterior(Simulation *S, int n, const int *exg, const double *ex){
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

void Simulation_InvalidateFieldCache(Simulation *S){
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
std::complex<double>* Simulation_GetCachedField(const Simulation *S, const Layer *layer){
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
void Simulation_AddFieldToCache(Simulation *S, const Layer *layer, size_t n, const std::complex<double> *P, size_t Plen){
	S4_TRACE("> Simulation_AddFieldToCache(S=%p, layer=%p, n=%d, P=%p) [omega=%f]\n", S, layer, n, P, S->omega[0]);
	FieldCache *f = (FieldCache*)S4_malloc(sizeof(FieldCache)+sizeof(std::complex<double>)*Plen);
	f->P = (std::complex<double>*)(f+1);
	memcpy(f->P, P, sizeof(std::complex<double>)*Plen);
	f->layer = layer;
	f->n = n;
	f->next = S->field_cache;
	S->field_cache = f;
	S4_TRACE("< Simulation_AddFieldToCache [omega=%f]\n", S->omega[0]);
}

void Simulation_SetExcitationType(Simulation *S, int type){
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
void Simulation_CopyExcitation(const Simulation *from, Simulation *to){
	S4_TRACE("> Simulation_CopyExcitation(from=%p, to=%p\n", from, to);
	memcpy(&(to->exc), &(from->exc), sizeof(Excitation));
	if(1 == to->exc.type){
		to->exc.layer = strdup(from->exc.layer);
	}else if(2 == to->exc.type){
		const int n = from->exc.sub.exterior.n;
		to->exc.sub.exterior.Gindex1 = (int*)malloc(sizeof(int) * 2*n);
		memcpy(to->exc.sub.exterior.Gindex1, from->exc.sub.exterior.Gindex1, sizeof(int) * 2*n);
		to->exc.sub.exterior.coeff = (double*)malloc(sizeof(double) * 2*n);
		memcpy(to->exc.sub.exterior.coeff, from->exc.sub.exterior.coeff, sizeof(double) * 2*n);
	}
	S4_TRACE("< Simulation_CopyExcitation\n");
}

int Simulation_GetSMatrix(Simulation *S, int from, int to, std::complex<double> *M){
	S4_TRACE("> Simulation_GetSMatrix(S=%p, from=%d, to=%d)\n", S, from, to);

	if(-1 != to && to < from){ return -3; }

	if(NULL == S->solution){
		int error = Simulation_InitSolution(S);
		if(0 != error){
			S4_TRACE("< Simulation_GetSMatrix (failed; Simulation_InitSolution returned %d)\n", error);
			return error;
		}
	}

	const size_t n2 = 2*S->n_G;
	const size_t n4 = 2*n2;

	// compute all bands
	Layer *SL = S->layer;
	Solution *sol = S->solution;
	LayerBands **Lbands = (LayerBands**)sol->layer_bands;
	int layer_count = 0;
	int layer_index = 0;
	while(NULL != SL){
		if(from <= layer_index && (-1 == to || layer_index <= to)){
			if(NULL == *Lbands){
				Simulation_ComputeLayerBands(S, SL, Lbands);
			}
			++Lbands;
			layer_count++;
		}
		SL = SL->next;
		++layer_index;
	}

	// Make arrays of q, kp, and phi
	double *lthick = (double*)S4_malloc(sizeof(double)*layer_count);
	int *lepstype = (int*)S4_malloc(sizeof(int)*layer_count);
	const std::complex<double> **lq   = (const std::complex<double> **)S4_malloc(sizeof(const std::complex<double> *)*layer_count*4);
	const std::complex<double> **lepsinv  = lq  + layer_count;
	const std::complex<double> **lkp  = lepsinv  + layer_count;
	const std::complex<double> **lphi = lkp + layer_count;

	Lbands = (LayerBands**)sol->layer_bands;
	for(layer_index = 0, SL = S->layer; NULL != SL; SL = SL->next, layer_index++){
		if(from <= layer_index && (-1 == to || layer_index <= to)){
			int i = layer_index;
			lthick[i] = SL->thickness;
			lq  [i] = Lbands[i]->q;
			lepsinv[i] = Lbands[i]->Epsilon_inv;
			lepstype[i] = Lbands[i]->epstype;
			lkp [i] = Lbands[i]->kp;
			lphi[i] = Lbands[i]->phi;
		}
	}

	GetSMatrix(layer_count, S->n_G, S->solution->kx, S->solution->ky, std::complex<double>(S->omega[0], S->omega[1]), lthick, lq, lepsinv, lepstype, lkp, lphi, M);

	S4_free(lq);
	S4_free(lepstype);
	S4_free(lthick);

	S4_TRACE("< Simulation_GetSMatrix\n");
	return 0;
}
