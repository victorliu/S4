#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <stdlib.h>
#include <S4.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#define S4V2_SIMULATION_TYPENAME "S4v2::Simulation"
#define S4V2_MATERIAL_TYPENAME   "S4v2::Material"
#define S4V2_LAYER_TYPENAME      "S4v2::Layer"

static S4_Simulation *R_S4v2_Simulation_check(SEXP sim){
	// TODO: check tag
	S4_Simulation *S = (S4_Simulation*)R_ExternalPtrAddr(sim);
	if(NULL == S){
		error("Encountered NULL S4 Simulation pointer");
	}
	return S;
}
static S4_MaterialID R_S4v2_Material_check(SEXP matid){
	// TODO: check tag
	S4_MaterialID id = (S4_MaterialID)R_ExternalPtrAddr(matid);
	return id-1;
}
static S4_LayerID R_S4v2_Layer_check(SEXP layid){
	// TODO: check tag
	S4_LayerID id = (S4_LayerID)R_ExternalPtrAddr(layid);
	return id-1;
}

static void R_S4v2_finalizer(SEXP ptr){
	S4_Simulation *S = R_S4v2_Simulation_check(ptr);
	S4_Simulation_Destroy(S);
	R_ClearExternalPtr(ptr); /* not really needed */
}

/*
lattice: double 2x2 matrix or single number
bases: integer scalar, vector, or nx2 matrix
*/
SEXP R_S4v2_Simulation_New(SEXP lattice, SEXP bases){
	S4_real Lr[4] = {0, 0, 0, 0};
	int nG = 1;
	int *G = NULL;
	SEXP ret = R_NilValue;

	if(isMatrix(lattice)){
		const double *p_lattice = REAL(lattice);
		Lr[0] = p_lattice[0];
		Lr[1] = p_lattice[1];
		Lr[2] = p_lattice[2];
		Lr[3] = p_lattice[3];
	}else{
		const double *p_lattice = REAL(lattice);
		Lr[0] = p_lattice[0];
	}

	if(isMatrix(bases)){
		nG = ncols(bases);
	}else{
		int nbases = length(bases);
		if(1 == nbases){
			nG = *INTEGER(bases);
		}else{
			nG = nbases;
		}
	}
	//printf("Lr = %f, %f, %f, %f, bases = %d\n", Lr[0], Lr[1], Lr[2], Lr[3], nG);

	S4_Simulation *S = S4_Simulation_New(Lr, nG, G);
	if(NULL == S){
		error_return("Failed to create a new S4 Simulation object");
	}
	ret = R_MakeExternalPtr(S, install(S4V2_SIMULATION_TYPENAME), R_NilValue);
	PROTECT(ret);
	R_RegisterCFinalizerEx(ret, R_S4v2_finalizer, TRUE);
	UNPROTECT(1);

	return ret;
}

SEXP R_S4v2_Simulation_Clone(SEXP sim){
	S4_Simulation *S = R_S4v2_Simulation_check(sim);
	S4_Simulation *T = S4_Simulation_Clone(S);
	SEXP ret = R_MakeExternalPtr(T, install("S4v2::Simulation"), R_NilValue);
	PROTECT(ret);
	R_RegisterCFinalizerEx(ret, R_S4v2_finalizer, TRUE);
	UNPROTECT(1);
	return ret;
}

SEXP R_S4v2_Simulation_SetFrequency(SEXP sim, SEXP freq){
	S4_Simulation *S = R_S4v2_Simulation_check(sim);
	const Rcomplex *f = COMPLEX(freq);
	S4_real fri[2] = { f->r, f->i };
	//printf("freq = %f, %f\n", fri[0], fri[1]);
	S4_Simulation_SetFrequency(S, fri);
	return R_NilValue;
}

SEXP R_S4v2_Simulation_AddMaterial(SEXP sim, SEXP name, SEXP epsilon){
	S4_Simulation *S = R_S4v2_Simulation_check(sim);
	const char *namestr = NULL;
	S4_real eps[18] = { 0 };
	int i, j;
	int type = 0;

	if(isString(name)){
		namestr = CHAR(STRING_ELT(name, 0));
	}

	if(isMatrix(epsilon)){
		type = 1;
		if(3 != nrows(epsilon) || 3 != ncols(epsilon)){
			error_return("Matrix epsilon must be 3x3");
		}
		if(isComplex(epsilon)){
			const Rcomplex *cplx = COMPLEX(epsilon);
			for(j = 0; j < 3; ++j){
				for(i = 0; i < 3; ++i){
					eps[2*(i+j*3)+0] = cplx[i+3*j].r;
					eps[2*(i+j*3)+1] = cplx[i+3*j].i;
				}
			}
		}else{
			const double *re = REAL(epsilon);
			for(j = 0; j < 3; ++j){
				for(i = 0; i < 3; ++i){
					eps[i+j*3] = re[i+3*j];
				}
			}
		}
		{
			S4_real abcde[10] = {
				eps[0], eps[1], eps[6], eps[7],
				eps[2], eps[3], eps[8], eps[9],
				eps[16], eps[17]
			};
			for(j = 0; j < 10; ++j){
				eps[j] = abcde[j];
			}
		}
		type = S4_MATERIAL_TYPE_XYTENSOR_COMPLEX;
	}else{
		if(isComplex(epsilon)){
			const Rcomplex *cplx = COMPLEX(epsilon);
			for(j = 0; j < 3; ++j){
				eps[8*j+0] = cplx[0].r;
				eps[8*j+1] = cplx[0].i;
			}
			type = S4_MATERIAL_TYPE_SCALAR_COMPLEX;
		}else{
			const double *re = REAL(epsilon);
			for(j = 0; j < 3; ++j){
				eps[8*j+0] = re[0];
			}
			type = S4_MATERIAL_TYPE_SCALAR_REAL;
		}
	}
	//printf("type = %d\n", type); for(j = 0; j < 18; ++j){ printf("\t%f\n", eps[j]); }
	S4_MaterialID M = S4_Simulation_SetMaterial(S, -1, namestr, type, eps);
	SEXP ret = R_MakeExternalPtr((void*)(M+1), install(S4V2_MATERIAL_TYPENAME), R_NilValue);
	return ret;
	/*
	PROTECT(ret);
	R_RegisterCFinalizerEx(ret, R_S4v2_finalizer, TRUE);

	SEXP result = PROTECT(allocVector(REALSXP, 1));
	REAL(result)[0] = asReal(a) + asReal(b);
	UNPROTECT(1);
	return
	*/
}
SEXP R_S4v2_Simulation_AddLayer(SEXP sim, SEXP name, SEXP thickness, SEXP copy, SEXP material){
	S4_Simulation *S = R_S4v2_Simulation_check(sim);
	const char *namestr = NULL;
	S4_MaterialID matid = R_S4v2_Material_check(material);
	S4_LayerID copyid = -1;
	S4_real t = *REAL(thickness);

	if(isString(name)){
		namestr = CHAR(STRING_ELT(name, 0));
	}
	if(!isNull(copy)){
		copyid = R_S4v2_Layer_check(copy);
	}
	//printf("namestr = %s, t = %f, copy = %d, matid = %d\n", namestr, t, copyid, matid);
	S4_LayerID L = S4_Simulation_SetLayer(S, -1, namestr, &t, copyid, matid);
	SEXP ret = R_MakeExternalPtr((void*)(L+1), install(S4V2_LAYER_TYPENAME), R_NilValue);
	return ret;
}
SEXP R_S4v2_Layer_SetRegionHalfwidths(SEXP sim, SEXP layer, SEXP material, SEXP shape, SEXP center, SEXP halfwidths, SEXP angle_frac){
	S4_Simulation *S = R_S4v2_Simulation_check(sim);
	S4_MaterialID matid = R_S4v2_Material_check(material);
	S4_LayerID layid = R_S4v2_Layer_check(layer);
	const double *vc = REAL(center);
	const double *vh = REAL(halfwidths);
	const S4_real c[2] = { vc[0], length(center) > 1 ? vc[1] : 0 };
	const S4_real h[2] = { vh[0], length(halfwidths) > 1 ? vh[1] : 0 };
	S4_real angle = *REAL(angle_frac);
	const char *shapestr = CHAR(STRING_ELT(shape, 0));
	int type;
	if(0 == strcmp("rectangle", shapestr)){
		type = S4_REGION_TYPE_RECTANGLE;
	}else if(0 == strcmp("circle", shapestr)){
		type = S4_REGION_TYPE_CIRCLE;
	}else if(0 == strcmp("interval", shapestr)){
		type = S4_REGION_TYPE_ELLIPSE;
	}else if(0 == strcmp("ellipse", shapestr)){
		type = S4_REGION_TYPE_INTERVAL;
	}else{
		error_return("Unrecognized shape");
	}
	//printf("L = %d, M = %d, type = %d, h = %f, %f, c = %f, %f, angle_frac = %f\n", layid, matid, type, h[0], h[1], c[0], c[1], angle);
	S4_Layer_SetRegionHalfwidths(S, layid, matid, type, h, c, &angle);
	return R_NilValue;
}
SEXP R_S4v2_Layer_SetRegionVertices(SEXP sim, SEXP layer, SEXP material, SEXP shape, SEXP center, SEXP vertices, SEXP angle_frac){
	S4_Simulation *S = R_S4v2_Simulation_check(sim);
	S4_MaterialID matid = R_S4v2_Material_check(material);
	S4_LayerID layid = R_S4v2_Layer_check(layer);
	const double *vc = REAL(center);
	const double *vh = REAL(vertices);
	const S4_real c[2] = { vc[0], length(center) > 1 ? vc[1] : 0 };
	int nv = ncols(vertices);
	S4_real angle = *REAL(angle_frac);
	const char *shapestr = CHAR(STRING_ELT(shape, 0));
	int type;
	if(0 == strcmp("polygon", shapestr)){
		type = S4_REGION_TYPE_POLYGON;
	}else{
		error_return("Unrecognized shape");
	}
	S4_Layer_SetRegionVertices(S, layid, matid, type, nv, vh, c, &angle);
	return R_NilValue;
}
SEXP R_S4v2_Simulation_ExcitationPlanewave(SEXP sim, SEXP k, SEXP u, SEXP cu, SEXP cv){
	S4_Simulation *S = R_S4v2_Simulation_check(sim);
	const double *vk = REAL(k);
	const double *vu = REAL(u);
	const Rcomplex *pcu = COMPLEX(cu);
	const Rcomplex *pcv = COMPLEX(cv);
	const double vcu[2] = { pcu->r, pcu->i };
	const double vcv[2] = { pcv->r, pcv->i };
	//printf("k = %f, %f, %f, u = %f, %f, %f, cu = %f, %f, cv = %f, %f\n", vk[0], vk[1], vk[2], vu[0], vu[1], vu[2], vcu[0], vcu[1], vcv[0], vcv[1]);
	S4_Simulation_ExcitationPlanewave(S, vk, vu, vcu, vcv);
	return R_NilValue;
}

SEXP R_S4v2_Simulation_GetEpsilon(SEXP sim, SEXP pt){
	S4_Simulation *S = R_S4v2_Simulation_check(sim);
	const double *r = REAL(pt);
	S4_real eps[2];
	Simulation_GetEpsilon(S, r, eps);
	Rcomplex reps;
	reps.r = eps[0];
	reps.i = eps[1];
	return ScalarComplex(reps);
}
SEXP R_S4v2_Layer_GetPowerFlux(SEXP sim, SEXP layer, SEXP offset){
	S4_Simulation *S = R_S4v2_Simulation_check(sim);
	S4_LayerID layid = R_S4v2_Layer_check(layer);
	S4_real off = *REAL(offset);
	S4_real power[4];
	S4_Simulation_GetPowerFlux(S, layid, &off, power);
	//printf("Power: %f, %f, %f, %f\n", power[0], power[1], power[2], power[3]);
	SEXP ret = PROTECT(allocVector(CPLXSXP, 2));
	Rcomplex *v = COMPLEX(ret);
	v[0].r = power[0]; v[0].i = power[2];
	v[1].r = power[1]; v[1].i = power[3];
	UNPROTECT(1);
	return ret;
}
SEXP R_S4v2_Layer_GetWaves(SEXP sim, SEXP layer){
	S4_Simulation *S = R_S4v2_Simulation_check(sim);
	S4_LayerID layid = R_S4v2_Layer_check(layer);
	S4_real *waves = NULL;
	int n, i, j;

	n = S4_Simulation_GetBases(S, NULL);
	const int n2 = 2*n;

	waves = (S4_real*)R_alloc(11*n2, sizeof(S4_real));
	S4_Simulation_GetWaves(S, layid, waves);
	SEXP ret = PROTECT(allocMatrix(REALSXP, 11, n2));
	double *pr = REAL(ret);
	for(j = 0; j < n2; ++j){
		for(i = 0; i < 11; ++i){
			pr[i+j*11] = waves[i+j*11];
		}
	}
	UNPROTECT(1);
	return ret;
}

