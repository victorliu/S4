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

#include "Python.h"
#include "config.h"

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <stdlib.h>
#include "S4.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void fft_init(void);
void fft_destroy(void);


static int CheckPyNumber(PyObject *obj){
	return PyFloat_Check(obj) || PyLong_Check(obj)
#if PY_MAJOR_VERSION < 3
		|| PyInt_Check(obj)
#endif
	;
}
static double AsNumberPyNumber(PyObject *obj){
	if(PyFloat_Check(obj)){
		return PyFloat_AsDouble(obj);
	}else if(PyLong_Check(obj)){
		return (double)PyLong_AsLong(obj);
	}
#if PY_MAJOR_VERSION < 3
	else if(PyInt_Check(obj)){
		return (double)PyInt_AsLong(obj);
	}
#endif
	return -1.0;
}
static int CheckPyComplex(PyObject *obj){
	return CheckPyNumber(obj) || PyComplex_Check(obj);
}
static void AsComplexPyComplex(PyObject *obj, double *re, double *im){
	if(CheckPyNumber(obj)){
		*re = AsNumberPyNumber(obj);
		*im = 0;
	}else if(PyComplex_Check(obj)){
		*re = PyComplex_RealAsDouble(obj);
		*im = PyComplex_ImagAsDouble(obj);
	}else{
		*re = -1.;
		*im = 0.;
	}
}

void HandleSolutionErrorCode(const char *fname, int code){
	static const char def[] = "An unknown error occurred";
	static const char* errstr[] = {
		def, /* 0 */
		"A memory allocation error occurred", /* 1 */
		def, /* 2 */
		def, /* 3 */
		def, /* 4 */
		def, /* 5 */
		def, /* 6 */
		def, /* 7 */
		def, /* 8 */
		"NumG was not set", /* 9 */
		"A layer copy referenced an unknown layer", /* 10 */
		"A layer copy referenced another layer copy", /* 11 */
		"A duplicate layer name was found", /* 12 */
		"Excitation layer name not found", /* 13 */
		"No layers exist in the structure", /* 14 */
		"A material name was not found", /* 15 */
		"Invalid patterning for 1D lattice", /* 16 */
		def
	};
	const char *str = def;
	if(0 < code && code <= 16){
		str = errstr[code];
		PyErr_Format(PyExc_RuntimeError, "%s: %s", fname, str);
	}else{
		PyErr_Format(PyExc_RuntimeError, "%s: %s, error code: %d", fname, str, code);
	}
}


#ifdef S4_DEBUG
# include "debug.h"
#endif

void threadsafe_init(void){
	extern void exactinit(void);
#ifdef HAVE_LAPACK
	extern float slamch_(const char *);
	extern double dlamch_(const char *);
#endif
	exactinit();
#ifdef HAVE_LAPACK
	slamch_("E");
	dlamch_("E");
#endif
	fft_init();
}
void threadsafe_destroy(void){
	fft_destroy();
}


struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

typedef struct{
	PyObject_HEAD
	Simulation S;
} S4Sim;
static PyTypeObject S4Sim_Type;

int lattice_converter(PyObject *obj, double *Lr){
	if(CheckPyNumber(obj)){
		Lr[0] = AsNumberPyNumber(obj);
		Lr[1] = 0;
		Lr[2] = 0;
		Lr[3] = 0;
		return 1;
	}else if(PyTuple_Check(obj) && (PyTuple_Size(obj) == 2)){
		unsigned i, j;
		for(j = 0; j < 2; ++j){
			PyObject *pj = PyTuple_GetItem(obj, j);
			if(!PyTuple_Check(pj) || (PyTuple_Size(pj) != 2)){
				PyErr_SetString(PyExc_TypeError, "2D lattice must be of the form ((ux,uy),(vx,vy))");
				return 0;
			}
			for(i = 0; i < 2; ++i){
				PyObject *pi = PyTuple_GetItem(pj, i);
				if(CheckPyNumber(pi)){
					Lr[i+j*2] = AsNumberPyNumber(pi);
				}else{
					PyErr_SetString(PyExc_TypeError, "2D lattice must be of the form ((ux,uy), (vx,vy))");
					return 0;
				}
			}
		}
		return 1;
	}
    PyErr_SetString(PyExc_TypeError, "Expected a number or 2-tuple");
	return 0;
}
struct epsilon_converter_data{
	int type; /* 0 = scalar, 1 = 3x3 tensor */
	double eps[18];
};
int epsilon_converter(PyObject *obj, struct epsilon_converter_data *data){
	double *eps = data->eps;
	if(CheckPyComplex(obj)){
		data->type = 0;
		AsComplexPyComplex(obj, &eps[0], &eps[1]);
		eps[2] = 0; eps[3] = 0;
		eps[4] = 0; eps[5] = 0;
		eps[6] = 0; eps[7] = 0;
		eps[8] = eps[0]; eps[9] = eps[1];
		eps[10] = 0; eps[11] = 0;
		eps[12] = 0; eps[13] = 0;
		eps[14] = 0; eps[15] = 0;
		eps[16] = eps[0]; eps[17] = eps[1];
		return 1;
	}else if(PyTuple_Check(obj) && (PyTuple_Size(obj) == 3)){
		unsigned i, j;
		data->type = 1;
		for(i = 0; i < 3; ++i){
			PyObject *pi = PyTuple_GetItem(obj, i);
			if(PyTuple_Check(pi) && (PyTuple_Size(pi) == 3)){
				for(j = 0; j < 3; ++j){
					PyObject *pj = PyTuple_GetItem(pi, j);
					if(CheckPyComplex(pj)){
						AsComplexPyComplex(pj, &eps[2*(i+j*3)+0], &eps[2*(i+j*3)+1]);
					}else{
						PyErr_SetString(PyExc_TypeError, "Material tensor must be a 3x3 matrix");
						return 0;
					}
				}
			}else{
				PyErr_SetString(PyExc_TypeError, "Material tensor must be a 3x3 matrix");
				return 0;
			}
		}
		return 1;
	}
    PyErr_SetString(PyExc_TypeError, "Expected a number or a 3x3-tuple");
	return 0;
}
struct polygon_converter_data{
	int nvert;
	double *vert;
};
int polygon_converter(PyObject *obj, struct polygon_converter_data *data){
	int i;
	if(!PyTuple_Check(obj)){
		return 0;
	}
	data->nvert = PyTuple_Size(obj);
	data->vert = (double*)malloc(sizeof(double) * 2 * data->nvert);
	for(i = 0; i < data->nvert; ++i){
		PyObject *pi = PyTuple_GetItem(obj, i);
		if(PyTuple_Check(pi) && (PyTuple_Size(pi) == 2)){
			unsigned j;
			for(j = 0; j < 2; ++j){
				PyObject *pj = PyTuple_GetItem(pi, j);
				if(CheckPyNumber(pj)){
					data->vert[2*i+j] = AsNumberPyNumber(pj);
				}else{
					free(data->vert);
					PyErr_SetString(PyExc_TypeError, "Polygon tensor must be a list of coordinate pairs");
					return 0;
				}
			}
		}else{
			free(data->vert);
			PyErr_SetString(PyExc_TypeError, "Polygon tensor must be a list of coordinate pairs");
			return 0;
		}
	}
	return 1;
}

static PyObject *S4Sim_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
	S4Sim *self;
	double Lr[4];
	Py_ssize_t nbasis;
	
	static char *kwlist[] = { "Lattice", "NumBasis", NULL };
	
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O&n:New", kwlist, &lattice_converter, &(Lr[0]), &nbasis)){ return NULL; }
	self = (S4Sim*)type->tp_alloc(type, 0);
	if(self != NULL){
		Simulation_Init(&(self->S));
		self->S.Lr[0] = Lr[0];
		self->S.Lr[1] = Lr[1];
		self->S.Lr[2] = Lr[2];
		self->S.Lr[3] = Lr[3];
		Simulation_MakeReciprocalLattice(&(self->S));
		Simulation_SetNumG(&(self->S), nbasis);
	}

	return (PyObject*)self;
}

static void S4Sim_dealloc(S4Sim* self){
	Simulation_Destroy(&(self->S));
	Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *S4Sim_Clone(S4Sim *self, PyObject *args){
	S4Sim *cpy;

	cpy = (S4Sim*)S4Sim_Type.tp_alloc(&S4Sim_Type, 0);
	if(cpy != NULL){
		Simulation_Clone(&(self->S), &(cpy->S));
	}
	return (PyObject*)cpy;
}
static PyObject *S4Sim_SetMaterial(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "Name", "Epsilon", NULL };
	const char *name;
	struct epsilon_converter_data epsdata;
	Material *M;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sO&:SetMaterial", kwlist, &name, &epsilon_converter, &epsdata)){ return NULL; }
	M = Simulation_GetMaterialByName(&(self->S), name, NULL);
	if(NULL == M){
		M = Simulation_AddMaterial(&(self->S));
		if(NULL == M){
			PyErr_Format(PyExc_MemoryError, "SetMaterial: There was a problem allocating the material named '%s'.", name);
			return NULL;
		}
		Material_Init(M, name, NULL);
	}
	
	if(0 == epsdata.type){
		M->eps.s[0] = epsdata.eps[0];
		M->eps.s[1] = epsdata.eps[1];
	}else{
		/* [ a b c ]    [ a b   ]
		 * [ d e f ] -> [ d e   ]
		 * [ g h i ]    [     i ]
		 */
		M->eps.abcde[0] = epsdata.eps[ 0]; M->eps.abcde[1] = epsdata.eps[ 1];
		M->eps.abcde[2] = epsdata.eps[ 2]; M->eps.abcde[3] = epsdata.eps[ 3];
		M->eps.abcde[4] = epsdata.eps[ 6]; M->eps.abcde[5] = epsdata.eps[ 7];
		M->eps.abcde[6] = epsdata.eps[ 8]; M->eps.abcde[7] = epsdata.eps[ 9];
		M->eps.abcde[8] = epsdata.eps[16]; M->eps.abcde[9] = epsdata.eps[17];
	}
	
	Py_RETURN_NONE;
}
static PyObject *S4Sim_AddLayer(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "Name", "Thickness", "Material", NULL };
	Layer *layer;
	const char *name;
	double thickness;
	const char *matname;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sds:AddLayer", kwlist, &name, &thickness, &matname)){ return NULL; }
	
	layer = Simulation_AddLayer(&(self->S));
	if(NULL == layer){
		PyErr_Format(PyExc_MemoryError, "AddLayer: There was a problem allocating the layer named '%s'.", name);
		return NULL;
	}
	Layer_Init(layer, name, thickness, matname, NULL);
	
	Py_RETURN_NONE;
}
static PyObject *S4Sim_AddLayerCopy(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "Name", "Thickness", "Layer", NULL };
	Layer *layer;
	const char *name;
	double thickness;
	const char *layername;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sds:AddLayerCopy", kwlist, &name, &thickness, &layername)){ return NULL; }
	
	layer = Simulation_AddLayer(&(self->S));
	if(NULL == layer){
		PyErr_Format(PyExc_MemoryError, "AddLayerCopy: There was a problem allocating the layer named '%s'.", name);
		return NULL;
	}
	Layer_Init(layer, name, thickness, NULL, layername);
	
	Py_RETURN_NONE;
}
static PyObject *S4Sim_SetLayerThickness(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "Layer", "Thickness", NULL };
	Layer *layer;
	const char *name;
	double thickness;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sd:SetLayerThickness", kwlist, &name, &thickness)){ return NULL; }
	
	layer = Simulation_GetLayerByName(&(self->S), name, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "SetLayerThickness: Layer named '%s' not found.", name);
		return NULL;
	}else{
		if(thickness < 0){
			PyErr_Format(PyExc_RuntimeError, "SetLayerThickness: Thickness must be non-negative.");
			return NULL;
		}
		Simulation_ChangeLayerThickness(&(self->S), layer, &thickness);
	}
	Py_RETURN_NONE;
}
static PyObject *S4Sim_RemoveLayerPatterns(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "Layer", "Thickness", "Material", NULL };
	Layer *layer;
	const char *name;
	
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s:ResetLayer", kwlist, &name)){ return NULL; }
	
	layer = Simulation_GetLayerByName(&(self->S), name, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "RemoveLayerPatterns: Layer named '%s' not found.", name);
		return NULL;
	}else{
		Simulation_RemoveLayerPatterns(&(self->S), layer);
	}
	Py_RETURN_NONE;
}
static PyObject *S4Sim_SetLayerPatternCircle(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "Layer", "Material", "Center", "Radius", NULL };
	Layer *layer;
	Material *M;
	const char *layername;
	const char *matname;
	int material_index;
	double center[2], radius;
	int ret;
	
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss(dd)d:SetLayerPatternCircle", kwlist, &layername, &matname, &center[0], &center[1], &radius)){ return NULL; }
	layer = Simulation_GetLayerByName(&(self->S), layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "SetLayerPatternCircle: Layer named '%s' not found.", layername);
		return NULL;
	}
	if(NULL != layer->copy){
		PyErr_Format(PyExc_RuntimeError, "SetLayerPatternCircle: Cannot pattern a layer copy.");
		return NULL;
	}
	M = Simulation_GetMaterialByName(&(self->S), matname, &material_index);
	if(NULL == M){
		PyErr_Format(PyExc_RuntimeError, "SetLayerPatternCircle: Material named '%s' not found.", matname);
		return NULL;
	}
	ret = Simulation_AddLayerPatternCircle(&(self->S), layer, material_index, center, radius);
	if(0 != ret){
		PyErr_Format(PyExc_MemoryError, "SetLayerPatternCircle: There was a problem allocating the pattern.");
		return NULL;
	}
	Py_RETURN_NONE;
}
static PyObject *S4Sim_SetLayerPatternEllipse(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "Layer", "Material", "Center", "Angle", "Halfwidths", NULL };
	Layer *layer;
	Material *M;
	const char *layername;
	const char *matname;
	int material_index;
	double center[2], tilt, halfwidths[2];
	int ret;
	
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss(dd)d(dd):SetLayerPatternEllipse", kwlist, &layername, &matname, &center[0], &center[1], &tilt, &halfwidths[0], &halfwidths[1])){ return NULL; }
	layer = Simulation_GetLayerByName(&(self->S), layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "SetLayerPatternEllipse: Layer named '%s' not found.", layername);
		return NULL;
	}
	if(NULL != layer->copy){
		PyErr_Format(PyExc_RuntimeError, "SetLayerPatternEllipse: Cannot pattern a layer copy.");
		return NULL;
	}
	M = Simulation_GetMaterialByName(&(self->S), matname, &material_index);
	if(NULL == M){
		PyErr_Format(PyExc_RuntimeError, "SetLayerPatternEllipse: Material named '%s' not found.", matname);
		return NULL;
	}
	ret = Simulation_AddLayerPatternEllipse(&(self->S), layer, material_index, center, (M_PI/180.)*tilt, halfwidths);
	if(0 != ret){
		PyErr_Format(PyExc_MemoryError, "SetLayerPatternEllipse: There was a problem allocating the pattern.");
		return NULL;
	}
	Py_RETURN_NONE;
}
static PyObject *S4Sim_SetLayerPatternRectangle(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "Layer", "Material", "Center", "Angle", "Halfwidths", NULL };
	Layer *layer;
	Material *M;
	const char *layername;
	const char *matname;
	int material_index;
	double center[2], tilt, halfwidths[2];
	int ret;
	
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss(dd)d(dd):SetLayerPatternRectangle", kwlist, &layername, &matname, &center[0], &center[1], &tilt, &halfwidths[0], &halfwidths[1])){ return NULL; }
	layer = Simulation_GetLayerByName(&(self->S), layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "SetLayerPatternRectangle: Layer named '%s' not found.", layername);
		return NULL;
	}
	if(NULL != layer->copy){
		PyErr_Format(PyExc_RuntimeError, "SetLayerPatternRectangle: Cannot pattern a layer copy.");
		return NULL;
	}
	M = Simulation_GetMaterialByName(&(self->S), matname, &material_index);
	if(NULL == M){
		PyErr_Format(PyExc_RuntimeError, "SetLayerPatternRectangle: Material named '%s' not found.", matname);
		return NULL;
	}
	ret = Simulation_AddLayerPatternRectangle(&(self->S), layer, material_index, center, (M_PI/180.)*tilt, halfwidths);
	if(0 != ret){
		PyErr_Format(PyExc_MemoryError, "SetLayerPatternRectangle: There was a problem allocating the pattern.");
		return NULL;
	}
	Py_RETURN_NONE;
}
static PyObject *S4Sim_SetLayerPatternPolygon(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "Layer", "Material", "Center", "Angle", "Vertices", NULL };
	Layer *layer;
	Material *M;
	const char *layername;
	const char *matname;
	int material_index;
	double center[2], tilt;
	struct polygon_converter_data polydata;
	int ret;
	
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss(dd)dO&:SetLayerPatternPolygon", kwlist, &layername, &matname, &center[0], &center[1], &tilt, &polygon_converter, &polydata)){ return NULL; }
	layer = Simulation_GetLayerByName(&(self->S), layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "SetLayerPatternPolygon: Layer named '%s' not found.", layername);
		return NULL;
	}
	if(NULL != layer->copy){
		PyErr_Format(PyExc_RuntimeError, "SetLayerPatternPolygon: Cannot pattern a layer copy.");
		return NULL;
	}
	M = Simulation_GetMaterialByName(&(self->S), matname, &material_index);
	if(NULL == M){
		PyErr_Format(PyExc_RuntimeError, "SetLayerPatternPolygon: Material named '%s' not found.", matname);
		return NULL;
	}
	ret = Simulation_AddLayerPatternPolygon(&(self->S), layer, material_index, center, (M_PI/180.)*tilt, polydata.nvert, polydata.vert);
	free(polydata.vert);
	if(0 != ret){
		PyErr_Format(PyExc_MemoryError, "SetLayerPatternPolygon: There was a problem allocating the pattern.");
		return NULL;
	}
	Py_RETURN_NONE;
}
static PyObject *S4Sim_SetExcitationPlanewave(S4Sim *self, PyObject *args, PyObject *kwds){
	int ret;
	static char *kwlist[] = { "IncidenceAngles", "sAmplitude", "pAmplitude", "Order", NULL };
	double angle[2];
	double pol_s[2], pol_p[2];
	Py_complex cs, cp;
	Py_ssize_t order = 0;
	cs.real = 0; cs.imag = 0;
	cp.real = 0; cp.imag = 0;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "(dd)|DDn:SetExcitationPlanewave", kwlist, &angle[0], &angle[1], &cs, &cp, &order)){ return NULL; }
	
	pol_s[0] = cs.real; pol_s[1] = cs.imag;
	pol_p[0] = cp.real; pol_p[1] = cp.imag;
	angle[0] *= (M_PI/180.);
	angle[1] *= (M_PI/180.);
	ret = Simulation_MakeExcitationPlanewave(&(self->S), angle, pol_s, pol_p, order);
	if(0 != ret){
		HandleSolutionErrorCode("SetExcitationPlanewave", ret);
		return NULL;
	}
	Py_RETURN_NONE;
}
static PyObject *S4Sim_SetFrequency(S4Sim *self, PyObject *args){
	Py_complex f;
	if(!PyArg_ParseTuple(args, "D:SetFrequency", &f)){ return NULL; }
	
	Simulation_DestroySolution(&(self->S));

	self->S.omega[0] = 2*M_PI*f.real;
	self->S.omega[1] = 2*M_PI*f.imag;
	if(self->S.omega[0] <= 0){
		PyErr_Warn(PyExc_RuntimeWarning, "A non-positive frequency was specified.");
	}
	if(self->S.omega[1] > 0){
		PyErr_Warn(PyExc_RuntimeWarning, "A frequency with positive imaginary part was specified.");
	}
	Py_RETURN_NONE;
}
static PyObject *S4Sim_GetReciprocalLattice(S4Sim *self, PyObject *args){
	if(!PyArg_ParseTuple(args, ":GetReciprocalLattice")){ return NULL; }
	return PyTuple_Pack(2,
		PyTuple_Pack(2,
			PyFloat_FromDouble(self->S.Lk[0]), PyFloat_FromDouble(self->S.Lk[1])
		),
		PyTuple_Pack(2,
			PyFloat_FromDouble(self->S.Lk[2]), PyFloat_FromDouble(self->S.Lk[3])
		)
	);
}
static PyObject *S4Sim_GetEpsilon(S4Sim *self, PyObject *args){
	int ret;
	double r[3], feps[2];
	if(!PyArg_ParseTuple(args, "ddd:GetEpsilon", &r[0], &r[1], &r[2])){ return NULL; }
	ret = Simulation_GetEpsilon(&(self->S), r, feps);
	if(0 != ret){
		HandleSolutionErrorCode("GetEpsilon", ret);
	}
	return PyComplex_FromDoubles(feps[0], feps[1]);
}
static PyObject *S4Sim_OutputLayerPatternPostscript(S4Sim *self, PyObject *args, PyObject *kwds){
	int ret;
	static char *kwlist[] = { "Layer", "Filename", NULL };
	const char *layername;
	const char *filename = NULL;
	Layer *layer;
	
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|s:OutputLayerPatternPostscript", kwlist, &layername, &filename)){ return NULL; }
	
	layer = Simulation_GetLayerByName(&(self->S), layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "OutputLayerPatternPostscript: Layer named '%s' not found.", layername);
		return NULL;
	}
	{
		FILE *fp = stdout;
		if(NULL != filename){
			fp = fopen(filename, "wb");
		}
		
		ret = Simulation_OutputLayerPatternDescription(&(self->S), layer, fp);
		if(0 != ret){
			HandleSolutionErrorCode("OutputLayerPatternDescription", ret);
			return NULL;
		}
		
		if(NULL != filename){
			fclose(fp);
		}
	}
	Py_RETURN_NONE;
}

static PyObject *S4Sim_OutputStructurePOVRay(S4Sim *self, PyObject *args, PyObject *kwds){
	PyErr_SetNone(PyExc_NotImplementedError); return NULL;
}

static PyObject *S4Sim_GetBasisSet(S4Sim *self, PyObject *args){
	int *G;
	int n, i, ret;
	PyObject *rv;
	if(!PyArg_ParseTuple(args, ":GetBasisSet")){ return NULL; }
		
	ret = Simulation_InitSolution(&(self->S));
	if(0 != ret){
		HandleSolutionErrorCode("GetBasisSet", ret);
		return NULL;
	}
	
	n = Simulation_GetNumG(&(self->S), &G);
	if(NULL == G){
		HandleSolutionErrorCode("GetBasisSet", 0);
		return NULL;
	}
	rv = PyTuple_New(n);
	for(i = 0; i < n; ++i){
		PyTuple_SetItem(rv, i, PyTuple_Pack(2, G[2*i+0], G[2*i+1]));
	}
	return rv;
}
static PyObject *S4Sim_GetAmplitudes(S4Sim *self, PyObject *args, PyObject *kwds){
	PyErr_SetNone(PyExc_NotImplementedError); return NULL;
}
static PyObject *S4Sim_GetPowerFlux(S4Sim *self, PyObject *args, PyObject *kwds){
	int ret;
	static char *kwlist[] = { "Layer", "zOffset", NULL };
	const char *layername;
	double offset = 0;
	double power[4];
	Layer *layer;
	
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|d:GetPowerFlux", kwlist, &layername, &offset)){ return NULL; }
	
	layer = Simulation_GetLayerByName(&(self->S), layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "GetPowerFlux: Layer named '%s' not found.", layername);
		return NULL;
	}
	ret = Simulation_GetPoyntingFlux(&(self->S), layer, offset, power);
	if(0 != ret){
		HandleSolutionErrorCode("GetPowerFlux", ret);
		return NULL;
	}
	
	return PyTuple_Pack(2,
		PyComplex_FromDoubles(power[0], power[2]),
		PyComplex_FromDoubles(power[1], power[3])
	);
}
static PyObject *S4Sim_GetPowerFluxByOrder(S4Sim *self, PyObject *args, PyObject *kwds){
	PyErr_SetNone(PyExc_NotImplementedError); return NULL;
}
static PyObject *S4Sim_GetStressTensorIntegral(S4Sim *self, PyObject *args, PyObject *kwds){
	PyErr_SetNone(PyExc_NotImplementedError); return NULL;
}
static PyObject *S4Sim_GetLayerVolumeIntegral(S4Sim *self, PyObject *args, PyObject *kwds){
	PyErr_SetNone(PyExc_NotImplementedError); return NULL;
}
static PyObject *S4Sim_GetLayerZIntegral(S4Sim *self, PyObject *args, PyObject *kwds){
	PyErr_SetNone(PyExc_NotImplementedError); return NULL;
}
static PyObject *S4Sim_GetFields(S4Sim *self, PyObject *args, PyObject *kwds){
	int ret;
	double r[3], fE[6],fH[6];
	if(!PyArg_ParseTuple(args, "ddd:GetFields", &r[0], &r[1], &r[2])){ return NULL; }
	
	ret = Simulation_GetField(&(self->S), r, fE, fH);
	if(0 != ret){
		HandleSolutionErrorCode("GetFields", ret);
		return NULL;
	}
	return PyTuple_Pack(2,
		PyTuple_Pack(3,
			PyComplex_FromDoubles(fE[0], fE[3]),
			PyComplex_FromDoubles(fE[1], fE[4]),
			PyComplex_FromDoubles(fE[2], fE[5])
		),
		PyTuple_Pack(3,
			PyComplex_FromDoubles(fH[0], fH[3]),
			PyComplex_FromDoubles(fH[1], fH[4]),
			PyComplex_FromDoubles(fH[2], fH[5])
		)
	);
}
static PyObject *S4Sim_GetFieldsOnGrid(S4Sim *self, PyObject *args, PyObject *kwds){
	int i, j, ret;
	static char *kwlist[] = { "zOffset", "NumSamples", "Format", "BaseFilename", NULL };
	Py_ssize_t nxy[2];
	double z;
	int len;
	double *Efields, *Hfields;
	const char *fmt;
	const char *fbasename = "field";
	char *filename;
	FILE *fp;
	int snxy[2];
	
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "d(nn)s|s:GetFieldsOnGrid", kwlist, &z, &nxy[0], &nxy[1], &fmt, &fbasename)){ return NULL; }
	len = strlen(fbasename);
	
	filename = (char*)malloc(sizeof(char) * (len+3));
	strcpy(filename, fbasename);
	filename[len+0] = '.';
	filename[len+2] = '\0';
	
	Efields = (double*)malloc(sizeof(double) * 2*3 * nxy[0] * nxy[1]);
	Hfields = (double*)malloc(sizeof(double) * 2*3 * nxy[0] * nxy[1]);
	
	snxy[0] = nxy[0];
	snxy[1] = nxy[1];
	ret = Simulation_GetFieldPlane(&(self->S), snxy, z, Efields, Hfields);
	if(0 != ret){
		HandleSolutionErrorCode("GetFields", ret);
		return NULL;
	}
	
	ret = 0;
	if(0 == strcmp("FileWrite", fmt)){
		filename[len+1] = 'E';
		fp = fopen(filename, "wb");
		for(i = 0; i < nxy[0]; ++i){
			for(j = 0; j < nxy[1]; ++j){
				fprintf(fp,
					"%d\t%d\t%.14g\t%.14g\t%.14g\t%.14g\t%.14g\t%.14g\n",
					i, j,
					Efields[2*(3*(i+j*nxy[0])+0)+0],
					Efields[2*(3*(i+j*nxy[0])+0)+1],
					Efields[2*(3*(i+j*nxy[0])+1)+0],
					Efields[2*(3*(i+j*nxy[0])+1)+1],
					Efields[2*(3*(i+j*nxy[0])+2)+0],
					Efields[2*(3*(i+j*nxy[0])+2)+1]
				);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		filename[len+1] = 'H';
		fp = fopen(filename, "wb");
		for(i = 0; i < nxy[0]; ++i){
			for(j = 0; j < nxy[1]; ++j){
				fprintf(fp,
					"%d\t%d\t%.14g\t%.14g\t%.14g\t%.14g\t%.14g\t%.14g\n",
					i, j,
					Hfields[2*(3*(i+j*nxy[0])+0)+0],
					Hfields[2*(3*(i+j*nxy[0])+0)+1],
					Hfields[2*(3*(i+j*nxy[0])+1)+0],
					Hfields[2*(3*(i+j*nxy[0])+1)+1],
					Hfields[2*(3*(i+j*nxy[0])+2)+0],
					Hfields[2*(3*(i+j*nxy[0])+2)+1]
				);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		free(Hfields);
		free(Efields);
		free(filename);
		Py_RETURN_NONE;
	}else if(0 == strcmp("FileAppend", fmt)){
		filename[len+1] = 'E';
		fp = fopen(filename, "ab");
		for(i = 0; i < nxy[0]; ++i){
			for(j = 0; j < nxy[1]; ++j){
				fprintf(fp,
					"%d\t%d\t%.14g\t%.14g\t%.14g\t%.14g\t%.14g\t%.14g\t%.14g\n",
					i, j, z,
					Efields[2*(3*(i+j*nxy[0])+0)+0],
					Efields[2*(3*(i+j*nxy[0])+0)+1],
					Efields[2*(3*(i+j*nxy[0])+1)+0],
					Efields[2*(3*(i+j*nxy[0])+1)+1],
					Efields[2*(3*(i+j*nxy[0])+2)+0],
					Efields[2*(3*(i+j*nxy[0])+2)+1]
				);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
		fclose(fp);
		filename[len+1] = 'H';
		fp = fopen(filename, "ab");
		for(i = 0; i < nxy[0]; ++i){
			for(j = 0; j < nxy[1]; ++j){
				fprintf(fp,
					"%d\t%d\t%.14g\t%.14g\t%.14g\t%.14g\t%.14g\t%.14g\t%.14g\n",
					i, j, z,
					Hfields[2*(3*(i+j*nxy[0])+0)+0],
					Hfields[2*(3*(i+j*nxy[0])+0)+1],
					Hfields[2*(3*(i+j*nxy[0])+1)+0],
					Hfields[2*(3*(i+j*nxy[0])+1)+1],
					Hfields[2*(3*(i+j*nxy[0])+2)+0],
					Hfields[2*(3*(i+j*nxy[0])+2)+1]
				);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
		fclose(fp);
		free(Hfields);
		free(Efields);
		free(filename);
		Py_RETURN_NONE;
	}else{ /* Array */
		unsigned k, i3;
		double *F[2] = { Efields, Hfields };
		PyObject *rv = PyTuple_New(2);
		
		for(k = 0; k < 2; ++k){
			PyObject *pk = PyTuple_New(nxy[0]);
			PyTuple_SetItem(rv, k, pk);
			for(i = 0; i < nxy[0]; ++i){
				PyObject *pi = PyTuple_New(nxy[1]);
				PyTuple_SetItem(pk, i, pi);
				for(j = 0; j < nxy[1]; ++j){
					PyObject *pj = PyTuple_New(3);
					PyTuple_SetItem(pi, j, pj);
					for(i3 = 0; i3 < 3; ++i3){
						PyTuple_SetItem(pj, i3, PyComplex_FromDoubles(
							F[k][2*(3*(i+j*nxy[0])+i3)+0],
							F[k][2*(3*(i+j*nxy[0])+i3)+1]
						));
					}
				}
			}
		}
		free(Hfields);
		free(Efields);
		free(filename);
		return rv;
	}
}
static PyObject *S4Sim_GetSMatrixDeterminant(S4Sim *self, PyObject *args, PyObject *kwds){
	PyErr_SetNone(PyExc_NotImplementedError); return NULL;
}

static PyMethodDef S4Sim_methods[] = {
	{"Clone"            , (PyCFunction)S4Sim_Clone, METH_NOARGS, PyDoc_STR("Clone() -> S4.Simulation")},
	/* Specification */
	{"SetMaterial"      , (PyCFunction)S4Sim_SetMaterial, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetMaterial(name,eps) -> None")},
	{"AddLayer"         , (PyCFunction)S4Sim_AddLayer, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("AddLayer(name,thickness,matname) -> None")},
	{"AddLayerCopy"     , (PyCFunction)S4Sim_AddLayerCopy, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("AddLayerCopy(name,thickness,layer) -> None")},
	{"SetLayerThickness", (PyCFunction)S4Sim_SetLayerThickness, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetLayerThickness(layer,thickness) -> None")},
	{"RemoveLayerPatterns", (PyCFunction)S4Sim_RemoveLayerPatterns, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("RemoveLayerPatterns(layer) -> None")},
	{"SetLayerPatternCircle", (PyCFunction)S4Sim_SetLayerPatternCircle, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetLayerPatternCircle(layer,matname,center,radius) -> None")},
	{"SetLayerPatternEllipse", (PyCFunction)S4Sim_SetLayerPatternEllipse, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetLayerPatternEllipse(layer,matname,center,angle,halfwidths) -> None")},
	{"SetLayerPatternRectangle", (PyCFunction)S4Sim_SetLayerPatternRectangle, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetLayerPatternRectangle(layer,matname,center,angle,halfwidths) -> None")},
	{"SetLayerPatternPolygon", (PyCFunction)S4Sim_SetLayerPatternPolygon, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetLayerPatternPolygon(layer,matname,center,angle,vertices) -> None")},
	{"SetExcitationPlanewave", (PyCFunction)S4Sim_SetExcitationPlanewave, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetExcitationPlanewave(angles,s_amp,p_amp) -> None")},
	{"SetFrequency", (PyCFunction)S4Sim_SetFrequency, METH_VARARGS, PyDoc_STR("SetFrequency(freq) -> None")},
	/* Outputs requiring no solutions */
	{"GetReciprocalLattice", (PyCFunction)S4Sim_GetReciprocalLattice, METH_NOARGS, PyDoc_STR("GetReciprocalLattice() -> ((px,py),(qx,qy))")},
	{"GetEpsilon", (PyCFunction)S4Sim_GetEpsilon, METH_VARARGS, PyDoc_STR("GetEpsilon(x,y,z) -> Complex")},
	{"OutputLayerPatternPostscript", (PyCFunction)S4Sim_OutputLayerPatternPostscript, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("OutputLayerPatternPostscript(layer,filename) -> None")},
	/* Outputs requiring solutions */
	{"OutputStructurePOVRay", (PyCFunction)S4Sim_OutputStructurePOVRay, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("OutputStructurePOVRay(filename) -> None")},
	{"GetBasisSet", (PyCFunction)S4Sim_GetBasisSet, METH_NOARGS, PyDoc_STR("GetBasisSet() -> Tuple")},
	{"GetAmplitudes", (PyCFunction)S4Sim_GetAmplitudes, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetAmplitudes(layer,zoffset) -> Tuple")},
	{"GetPowerFlux", (PyCFunction)S4Sim_GetPowerFlux, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetPowerFlux(layer,zoffset) -> (forw,back)")},
	{"GetPowerFluxByOrder", (PyCFunction)S4Sim_GetPowerFluxByOrder, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetPowerFluxByOrder(layer,zoffset) -> Tuple")},
	{"GetStressTensorIntegral", (PyCFunction)S4Sim_GetStressTensorIntegral, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetStressTensorIntegral(layer,zoffset) -> Complex")},
	{"GetLayerVolumeIntegral", (PyCFunction)S4Sim_GetLayerVolumeIntegral, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetLayerVolumeIntegral(layer,which) -> Complex")},
	{"GetLayerZIntegral", (PyCFunction)S4Sim_GetLayerZIntegral, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetLayerZIntegral(layer,which,pos) -> Complex")},
	{"GetFields", (PyCFunction)S4Sim_GetFields, METH_VARARGS, PyDoc_STR("GetFields(x,y,z) -> (Tuple,Tuple)")},
	{"GetFieldsOnGrid", (PyCFunction)S4Sim_GetFieldsOnGrid, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetFieldsOnGrid(z,nsamples,format,filename) -> Tuple")},
	{"GetSMatrixDeterminant", (PyCFunction)S4Sim_GetSMatrixDeterminant, METH_NOARGS, PyDoc_STR("GetSMatrixDeterminant() -> Tuple")},
	{NULL, NULL}
};
static PyTypeObject S4Sim_Type = {
	/* The ob_type field must be initialized in the module init function
	 * to be portable to Windows without using C++. */
	PyVarObject_HEAD_INIT(NULL, 0)
	"S4.Simulation",    /*tp_name*/
	sizeof(S4Sim),      /*tp_basicsize*/
	0,                  /*tp_itemsize*/
	/* methods */
	(destructor)S4Sim_dealloc, /*tp_dealloc*/
	0,                  /*tp_print*/
	0,                  /*tp_getattr*/
	0,                  /*tp_setattr*/
	0,                  /*tp_reserved*/
	0,                  /*tp_repr*/
	0,                  /*tp_as_number*/
	0,                  /*tp_as_sequence*/
	0,                  /*tp_as_mapping*/
	0,                  /*tp_hash*/
	0,                  /*tp_call*/
	0,                  /*tp_str*/
	0,                  /*tp_getattro*/
	0,                  /*tp_setattro*/
	0,                  /*tp_as_buffer*/
	Py_TPFLAGS_DEFAULT, /*tp_flags*/
	0,                  /*tp_doc*/
	0,                  /*tp_traverse*/
	0,                  /*tp_clear*/
	0,                  /*tp_richcompare*/
	0,                  /*tp_weaklistoffset*/
	0,                  /*tp_iter*/
	0,                  /*tp_iternext*/
	S4Sim_methods,      /*tp_methods*/
	0,                  /*tp_members*/
	0,                  /*tp_getset*/
	0,                  /*tp_base*/
	0,                  /*tp_dict*/
	0,                  /*tp_descr_get*/
	0,                  /*tp_descr_set*/
	0,                  /*tp_dictoffset*/
	0,                  /*tp_init*/
	0,                  /*tp_alloc*/
	S4Sim_new,          /*tp_new*/
	0,                  /*tp_free*/
	0,                  /*tp_is_gc*/
};

static PyObject *S4_new(PyObject *self, PyObject *args, PyObject *kwds){
	return (PyObject*)S4Sim_new(&S4Sim_Type, args, kwds);
}

static PyMethodDef S4_funcs[] = {
	{"New", (PyCFunction)S4_new, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("New() -> new S4 simulation object")},
	{NULL , NULL} /* sentinel */
};


/* Initialization function for the module (*must* be called PyInit_FunctionSampler1D) */
PyDoc_STRVAR(module_doc, "Stanford Stratified Structure Solver (S4): Fourier Modal Method.");

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef S4_module = {
	PyModuleDef_HEAD_INIT,
	"S4",
	module_doc,
	sizeof(struct module_state),
	S4_funcs,
	NULL,
	NULL,
	NULL,
	NULL
};
#define INITERROR return NULL
PyObject * PyInit_S4(void)
#else
#define INITERROR return
PyMODINIT_FUNC initS4(void)
#endif
{
	PyObject *m = NULL;
	
	//threadsafe_init();

	/* Finalize the type object including setting type of the new type
	 * object; doing it here is required for portability, too. */
	if(PyType_Ready(&S4Sim_Type) < 0){ INITERROR; }

	/* Create the module and add the functions */
#if PY_MAJOR_VERSION >= 3
	m = PyModule_Create(&S4_module);
#else
	m = Py_InitModule3("S4", S4_funcs, module_doc);
#endif
	if(m == NULL){ INITERROR; }

	struct module_state *st = GETSTATE(m);

	st->error = PyErr_NewException("S4.Error", NULL, NULL);
	if(st->error == NULL){
		Py_DECREF(m);
		INITERROR;
	}
	Py_INCREF(st->error);
	PyModule_AddObject(m, "Error", st->error);
	
#if PY_MAJOR_VERSION >= 3
	return m;
#endif
}

