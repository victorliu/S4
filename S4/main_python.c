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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#ifdef WIN32
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <stdlib.h>
#include "S4.h"
#include "convert.h"
#include "SpectrumSampler.h"
#include "cubature.h"
#include "Interpolator.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef Bool
#define Bool unsigned char
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
static PyObject *FromIntPyDefInt(int i){
#if PY_MAJOR_VERSION < 3
	return PyInt_FromLong(i);
#else
	return PyLong_FromLong(i);
#endif
}
static int CheckPyInt(PyObject *obj){
#if PY_MAJOR_VERSION < 3
	return PyInt_Check(obj);
#else
	return PyLong_Check(obj);
#endif
}
static long GetPyInt(PyObject *obj){
#if PY_MAJOR_VERSION < 3
	return PyInt_AsLong(obj);
#else
	return PyLong_AsLong(obj);
#endif
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
	S4_Simulation *S;
} S4Sim;

typedef struct
{
	PyObject_HEAD
	Interpolator I;
}S4Interpolator;

typedef struct
{
	PyObject_HEAD
	SpectrumSampler SpecS;
}S4SpectrumSampler;

/*
Description: this structure is defined to be used for argument convert function.
			 to make it understanding, the members' name keep th same with the C
			 API.
Note: the 'xy' pointer's resource is heap allocated, remember to free it.
*/
typedef struct
{
	int n;
	int ny;
	double *xy;
}S4Interpolator_Data;

/*
Description: the use of this structure is similar as structure 'S4Interpolator_Data'.
Note: both 'axg' and 'ex' is heap alloated, remember to free it.
*/
typedef struct
{
	int n;
	int *exg;
	double *ex;
}S4Excitation_Data;

/*
Description: the use is same as structure 'S4Interpolator_Data'.
			 But the members' name is same as the arguments of
			 python API(ignore case).
*/
typedef struct
{
	double freqStart;
	double freqEnd;
	int initialNumPoints;
	double rangeThreshold;
	double maxBend;
	double minimumSpacing;
	Bool parallelize;
}S4SpectrumSampler_Data;

static PyTypeObject S4Sim_Type;
static PyTypeObject S4Interpolator_Type;
static PyTypeObject S4SpectrumSampler_Type;

int bool_converter(PyObject *obj, int *b){
	if(PyBool_Check(obj)){
		if(Py_True == obj){
			*b = 1;
		}else{
			*b = 0;
		}
		return 1;
	}
    PyErr_SetString(PyExc_TypeError, "Expected a boolean");
	return 0;
}

struct lanczos_smoothing_settings{
	int set; /* whether the user set this option at all */
	int use;
	int set_power;
	int power;
	int set_width;
	double width;
};

int lanczos_converter(PyObject *obj, struct lanczos_smoothing_settings *s){
	s->set = 1;
	s->set_power = 0;
	s->set_width = 0;
	if(PyBool_Check(obj)){
		s->use = (Py_True == obj);
		return 1;
	}else if(PyDict_Check(obj)){
		PyObject *val;
		s->use = 1;
		if((val = PyDict_GetItemString(obj, "Width"))){
			s->set_width = 1;
			if(CheckPyNumber(val)){
				s->width = AsNumberPyNumber(val);
				if(s->width <= 0){
					PyErr_SetString(PyExc_ValueError, "Width must be positive");
					return 0;
				}
			}else{
				PyErr_SetString(PyExc_TypeError, "Width must be a positive number");
				return 0;
			}
		}
		if((val = PyDict_GetItemString(obj, "Power"))){
			s->set_power = 1;
			if(CheckPyInt(val)){
				s->power = GetPyInt(val);
				if(s->power <= 0){
					PyErr_SetString(PyExc_ValueError, "Power must be positive");
					return 0;
				}
			}else{
				PyErr_SetString(PyExc_TypeError, "Power must be a positive integer");
				return 0;
			}
		}
		return 1;
	}
    PyErr_SetString(PyExc_TypeError, "Expected a boolean or dictionary");
	return 0;
}

/*
Descritpion: used by function S4Sim_SetExcitationExterior() to parse arguments.
Parameters:
	data: the structure to store the arguments.
return :
	0: failed.
	1: success.
*/
int excitation_converter(PyObject *obj, S4Excitation_Data *data)
{
	if(!PyTuple_Check(obj))
	{
		PyErr_SetString(PyExc_TypeError, "parameter must be a tuple.");
		return 0;
	}

	//the first calling, exg and ex should to setted to NULL
	if(NULL == data->exg || NULL == data->ex)
	{
		data->n = PyTuple_Size(obj);	//return the size info needed to malloc.
		return 1;
	}

	for(int i = 0; i < data->n; i++)
	{
		PyObject *pi = PyTuple_GetItem(obj, i);
		char *pol;
		Py_ssize_t polLen;
		PyObject *pj;
		if(!PyTuple_Check(pi))
		{
			PyErr_SetString(PyExc_TypeError, "the tuple item must be a tuple.");
			return 0;
		}

		//get G index
		pj = PyTuple_GetItem(pi, 0);
		if(!CheckPyInt(pj))
		{
			PyErr_SetString(PyExc_TypeError, "the G index must be a integer.");
			return 0;
		}
		data->exg[2 * i + 0] = PyInt_AsLong(pj);

		//get polarization: 'x' or 'y'
		pj = PyTuple_GetItem(pi, 1);
		if(!PyString_Check(pj))
		{
			PyErr_SetString(PyExc_TypeError, "polalization should be specified by 'x' or 'y'.");
			return 0;
		}
		PyString_AsStringAndSize(pj, &pol, &polLen);
		if(1 != polLen || ('x' != pol[0] && 'y' != pol[0]))
		{
			PyErr_SetString(PyExc_TypeError, "polalization should be specified by 'x' or 'y'.");
			return 0;
		}
		if('x' == pol[0])
			data->exg[2 * i + 1] = 0;
		else
			data->exg[2 * i + 1] = 1;

		//get the complex coeffcient
		pj = PyTuple_GetItem(pi, 2);
		if(!CheckPyComplex(pj))
		{
			PyErr_Format(PyExc_TypeError, "cofficent should be complex");
			return 0;
		}
		AsComplexPyComplex(pj, &data->ex[2 * i + 0], &data->ex[2 * i + 1]);
	}
	return 1;
}

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
						AsComplexPyComplex(pj, &eps[2*(3*i+j)+0], &eps[2*(3*i+j)+1]);
					}else{
						PyErr_SetString(PyExc_TypeError, "S4_Material tensor must be a 3x3 matrix");
						return 0;
					}
				}
			}else{
				PyErr_SetString(PyExc_TypeError, "S4_Material tensor must be a 3x3 matrix");
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

/*
Description: the use is similar to excitation_converter().
*/
static int interpolator_table_converter(PyObject *args, S4Interpolator_Data *data)
{
	PyObject *pi, *pj;
	if(NULL == data)
		return 0;
	if(!PyTuple_Check(args))
	{
		PyErr_SetString(PyExc_TypeError, "the 'Table' argument isn't a tuple.");
		return 0;
	}

	data->n = PyTuple_Size(args);
	if(0 == data->n)
	{
		PyErr_SetString(PyExc_TypeError, "the 'Table' argument can't be empty");
		return 0;
	}
	for(int i = 0, ld = data->ny + 1; i < data->n; i++)
	{
		pi = PyTuple_GetItem(args, i);
		if(2 != PyTuple_Size(pi))
		{
			PyErr_SetString(PyExc_TypeError, "a the 'Table' should be like:((x, (y1, y2,...)), (x, (y1, y2,...)),...)");
			return 0;
		}
		pj = PyTuple_GetItem(pi, 1);
		if(!PyTuple_Check(pj))
		{
			PyErr_SetString(PyExc_TypeError, "the 'Table' should be like:((x, (y1, y2,...)), (x, (y1, y2,...)),...)");
			return 0;
		}

		if(NULL == data->xy)
		{
			data->ny = PyTuple_Size(pj);
			return 1;
		}
		data->xy[i*ld + 0] = PyFloat_AsDouble(PyTuple_GetItem(pi, 0));
		for(int j = 0; j < data->ny; j++)
			data->xy[i*ld + j + 1] = PyFloat_AsDouble(PyTuple_GetItem(pj, j));
	}
	return 1;
}

static PyObject *S4Interpolator_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	static char *kwlist[] = {"Type", "Table", NULL};
	const char *typeName;
	S4Interpolator *self;
	S4Interpolator_Data interData = {0, 0, NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sO&:interpolator_new", kwlist, &typeName, &interpolator_table_converter, &interData))
		return NULL;
	interData.xy = (double*)malloc(sizeof(double) * interData.n * (interData.ny + 1));
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sO&:interpolator_new", kwlist, &typeName, &interpolator_table_converter, &interData))
	{
		free(interData.xy); interData.xy = NULL;
		return NULL;
	}

	self = (S4Interpolator*)type->tp_alloc(type, 0);
	if(NULL != self)
	{
		Interpolator_type inter_type;
		if(0 == strcmp("linear", typeName))
			inter_type = Interpolator_LINEAR;
		else if(0 == strcmp("cublic spline", typeName))
			inter_type = Interpolator_CUBIC_SPLINE;
		else if(0 == strcmp("cubic hermite spline", typeName))
			inter_type = Interpolator_CUBIC_HERMITE_SPLINE;
		else
		{
			PyErr_SetString(PyExc_TypeError, "the 'type' should be 'linear'/'cubic spline'/'cubic hermite spline'.");
			free(interData.xy); interData.xy = NULL;
			return NULL;
		}
		self->I = Interpolator_New(interData.n, interData.ny, interData.xy, inter_type);
	}
	free(interData.xy); interData.xy = NULL;
	return (PyObject*)self;
}

static PyObject *S4Interpolator_Get(S4Interpolator *self, PyObject *args, PyObject *kwds)
{
	static char *kwlist[] = {"X", NULL};
	double x;
	double *ys;
	int ny;
	PyObject *ret;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "d:Get", kwlist, &x))
		return NULL;
	ys = Interpolator_Get(self->I, x, &ny);
	if(NULL == ys)
		Py_RETURN_NONE;
	ret = PyTuple_New(ny);
	for(int i = 0; i < ny; i++)
		PyTuple_SetItem(ret, i, Py_BuildValue("d", ys[i]));
	return ret;
}

static PyObject *S4SpectrumSampler_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	static char * kwlist[] = {"FreqStart", "FreqEnd", "InitialNumPoints", "RangeThreshold", \
		"MaxBend", "MinimumSpacing", "Parallelize", NULL};
	double x0, x1;
	SpectrumSampler_Options options = {33, 0.001, 10, 1e-6, 0};
	PyObject *py_expectBool = NULL;
	S4SpectrumSampler *self;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "dd|i|d|d|d|O!:SpectrumSampler_New", \
		kwlist, &x0, &x1, &options.initial_num_points, &options.range_threshold,\
		&options.max_bend, &options.min_dx, &PyBool_Type, &py_expectBool))
		return NULL;
	if(NULL != py_expectBool)
		options.parallelize = PyObject_IsTrue(py_expectBool);

	self = (S4SpectrumSampler*)type->tp_alloc(type, 0);
	if(NULL == self)
		return NULL;
	self->SpecS = SpectrumSampler_New(x0, x1, &options);
	return (PyObject*)self;
}

static PyObject *S4_NewSpectrumSampler(PyObject *self, PyObject *args, PyObject *kwds)
{
	return S4SpectrumSampler_new(&S4SpectrumSampler_Type, args, kwds);
}

static PyObject *S4Sim_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
	S4Sim *self;
	double Lr[4];
	Py_ssize_t nbasis;
	static char *kwlist[] = { "Lattice", "NumBasis", NULL };

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O&n:New", kwlist, &lattice_converter, &(Lr[0]), &nbasis)){ return NULL; }
	self = (S4Sim*)type->tp_alloc(type, 0);
	if(self != NULL){
		self->S = S4_Simulation_New(Lr, nbasis, NULL);
	}

	return (PyObject*)self;
}

static void S4Sim_dealloc(S4Sim* self){
	S4_Simulation_Destroy(self->S);
	Py_TYPE(self)->tp_free((PyObject*)self);
}

static void S4Interpolator_dealloc(S4Interpolator *self)
{
	Interpolator_Destroy(self->I);
	Py_TYPE(self)->tp_free((PyObject*) self);
}

static void S4SpectrumSampler_dealloc(S4SpectrumSampler *self)
{
	SpectrumSampler_Destroy(self->SpecS);
	Py_TYPE(self)->tp_free((PyObject*) self);
}

static PyObject *S4Sim_Clone(S4Sim *self, PyObject *args){
	S4Sim *cpy;

	cpy = (S4Sim*)S4Sim_Type.tp_alloc(&S4Sim_Type, 0);
	if(cpy != NULL){
		cpy->S = S4_Simulation_Clone(self->S);
	}
	return (PyObject*)cpy;
}

static PyObject *S4Sim_ConvertUnits(S4Sim *self, PyObject *args)
{
	double value;
	const char *from_units = NULL;
	const char *to_units = NULL;
	if(!PyArg_ParseTuple(args, "dss", &value, from_units, to_units))
		return NULL;
	if(0 == convert_units(&value, from_units, to_units))
		return Py_BuildValue("d", value);
	return NULL;
}

static PyObject *S4Sim_SetMaterial(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "Name", "Epsilon", NULL };
	const char *name;
	struct epsilon_converter_data epsdata;
	S4_Material *M;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sO&:SetMaterial", kwlist, &name, &epsilon_converter, &epsdata)){ return NULL; }
	M = S4_Simulation_GetMaterialByName(self->S, name);
	if(NULL == M){
		M = Simulation_AddMaterial(self->S);
		if(NULL == M){
			PyErr_Format(PyExc_MemoryError, "SetMaterial: There was a problem allocating the material named '%s'.", name);
			return NULL;
		}
		if(0 == epsdata.type){
			Material_Init(M, name, NULL);
		}else{
			Material_InitTensor(M, name, NULL);
		}
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

static PyObject *S4Sim_AddMaterial(S4Sim *self, PyObject *args, PyObject *kwds)
{
	return S4Sim_SetMaterial(self, args, kwds);
}

static PyObject *S4Sim_AddLayer(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "Name", "Thickness", "S4_Material", NULL };
	S4_Layer *layer;
	const char *name;
	double thickness;
	const char *matname;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sds:AddLayer", kwlist, &name, &thickness, &matname)){ return NULL; }

	layer = Simulation_AddLayer(self->S);
	if(NULL == layer){
		PyErr_Format(PyExc_MemoryError, "AddLayer: There was a problem allocating the layer named '%s'.", name);
		return NULL;
	}
	Layer_Init(layer, name, thickness, matname, NULL);

	Py_RETURN_NONE;
}

static PyObject *S4Sim_SetLayer(S4Sim *self, PyObject *args, PyObject *kwds)
{
	static char *kwlist[] = { "Name", "Thickness", "S4_Material", NULL };
	const char *name, *material = NULL;
	double thickness;
	S4_Layer *layer;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sd|s:SetLayer", kwlist, &name, &thickness, &material))
		return NULL;
	layer = Simulation_GetLayerByName(self->S, name, NULL);
	if(NULL == layer)
		S4Sim_AddLayer(self, args, kwds);
	else
	{
		layer->thickness = thickness;
		if(NULL != material)
			layer->material = strdup(material);
		Simulation_RemoveLayerPatterns(self->S, layer);
	}
	Py_RETURN_NONE;
}

static PyObject *S4Sim_AddLayerCopy(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "Name", "Thickness", "S4_Layer", NULL };
	S4_Layer *layer;
	const char *name;
	double thickness;
	const char *layername;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sds:AddLayerCopy", kwlist, &name, &thickness, &layername)){ return NULL; }

	layer = Simulation_AddLayer(self->S);
	if(NULL == layer){
		PyErr_Format(PyExc_MemoryError, "AddLayerCopy: There was a problem allocating the layer named '%s'.", name);
		return NULL;
	}
	Layer_Init(layer, name, thickness, NULL, layername);

	Py_RETURN_NONE;
}
static PyObject *S4Sim_SetLayerThickness(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "S4_Layer", "Thickness", NULL };
	S4_Layer *layer;
	const char *name;
	double thickness;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sd:SetLayerThickness", kwlist, &name, &thickness)){ return NULL; }

	layer = Simulation_GetLayerByName(self->S, name, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "SetLayerThickness: S4_Layer named '%s' not found.", name);
		return NULL;
	}else{
		if(thickness < 0){
			PyErr_Format(PyExc_RuntimeError, "SetLayerThickness: Thickness must be non-negative.");
			return NULL;
		}
		Simulation_ChangeLayerThickness(self->S, layer, &thickness);
	}
	Py_RETURN_NONE;
}
static PyObject *S4Sim_RemoveLayerRegions(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "S4_Layer", NULL };
	S4_Layer *layer;
	const char *name;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s:RemoveLayerRegions", kwlist, &name)){ return NULL; }

	layer = Simulation_GetLayerByName(self->S, name, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "RemoveLayerRegions: S4_Layer named '%s' not found.", name);
		return NULL;
	}else{
		Simulation_RemoveLayerPatterns(self->S, layer);
	}
	Py_RETURN_NONE;
}
static PyObject *S4Sim_SetRegionCircle(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "S4_Layer", "S4_Material", "Center", "Radius", NULL };
	S4_Layer *layer;
	S4_Material *M;
	const char *layername;
	const char *matname;
	int material_index;
	double center[2], radius;
	int ret;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss(dd)d:SetRegionCircle", kwlist, &layername, &matname, &center[0], &center[1], &radius)){ return NULL; }
	layer = Simulation_GetLayerByName(self->S, layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "SetRegionCircle: S4_Layer named '%s' not found.", layername);
		return NULL;
	}
	if(NULL != layer->copy){
		PyErr_Format(PyExc_RuntimeError, "SetRegionCircle: Cannot pattern a layer copy.");
		return NULL;
	}
	M = Simulation_GetMaterialByName(self->S, matname, &material_index);
	if(NULL == M){
		PyErr_Format(PyExc_RuntimeError, "SetRegionCircle: S4_Material named '%s' not found.", matname);
		return NULL;
	}
	ret = Simulation_AddLayerPatternCircle(self->S, layer, material_index, center, radius);
	if(0 != ret){
		PyErr_Format(PyExc_MemoryError, "SetRegionCircle: There was a problem allocating the pattern.");
		return NULL;
	}
	Py_RETURN_NONE;
}
static PyObject *S4Sim_SetRegionEllipse(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "S4_Layer", "S4_Material", "Center", "Angle", "Halfwidths", NULL };
	S4_Layer *layer;
	S4_Material *M;
	const char *layername;
	const char *matname;
	int material_index;
	double center[2], tilt, halfwidths[2];
	int ret;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss(dd)d(dd):SetRegionEllipse", kwlist, &layername, &matname, &center[0], &center[1], &tilt, &halfwidths[0], &halfwidths[1])){ return NULL; }
	layer = Simulation_GetLayerByName(self->S, layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "SetRegionEllipse: S4_Layer named '%s' not found.", layername);
		return NULL;
	}
	if(NULL != layer->copy){
		PyErr_Format(PyExc_RuntimeError, "SetRegionEllipse: Cannot pattern a layer copy.");
		return NULL;
	}
	M = Simulation_GetMaterialByName(self->S, matname, &material_index);
	if(NULL == M){
		PyErr_Format(PyExc_RuntimeError, "SetRegionEllipse: S4_Material named '%s' not found.", matname);
		return NULL;
	}
	ret = Simulation_AddLayerPatternEllipse(self->S, layer, material_index, center, (M_PI/180.)*tilt, halfwidths);
	if(0 != ret){
		PyErr_Format(PyExc_MemoryError, "SetRegionEllipse: There was a problem allocating the pattern.");
		return NULL;
	}
	Py_RETURN_NONE;
}
static PyObject *S4Sim_SetRegionRectangle(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "S4_Layer", "S4_Material", "Center", "Angle", "Halfwidths", NULL };
	S4_Layer *layer;
	S4_Material *M;
	const char *layername;
	const char *matname;
	int material_index;
	double center[2], tilt, halfwidths[2];
	int ret;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss(dd)d(dd):SetRegionRectangle", kwlist, &layername, &matname, &center[0], &center[1], &tilt, &halfwidths[0], &halfwidths[1])){ return NULL; }
	layer = Simulation_GetLayerByName(self->S, layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "SetRegionRectangle: S4_Layer named '%s' not found.", layername);
		return NULL;
	}
	if(NULL != layer->copy){
		PyErr_Format(PyExc_RuntimeError, "SetRegionRectangle: Cannot pattern a layer copy.");
		return NULL;
	}
	M = Simulation_GetMaterialByName(self->S, matname, &material_index);
	if(NULL == M){
		PyErr_Format(PyExc_RuntimeError, "SetRegionRectangle: S4_Material named '%s' not found.", matname);
		return NULL;
	}
	ret = Simulation_AddLayerPatternRectangle(self->S, layer, material_index, center, (M_PI/180.)*tilt, halfwidths);
	if(0 != ret){
		PyErr_Format(PyExc_MemoryError, "SetRegionRectangle: There was a problem allocating the pattern.");
		return NULL;
	}
	Py_RETURN_NONE;
}
static PyObject *S4Sim_SetRegionPolygon(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = { "S4_Layer", "S4_Material", "Center", "Angle", "Vertices", NULL };
	S4_Layer *layer;
	S4_Material *M;
	const char *layername;
	const char *matname;
	int material_index;
	double center[2], tilt;
	struct polygon_converter_data polydata;
	int ret;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss(dd)dO&:SetRegionPolygon", kwlist, &layername, &matname, &center[0], &center[1], &tilt, &polygon_converter, &polydata)){ return NULL; }
	layer = Simulation_GetLayerByName(self->S, layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "SetRegionPolygon: S4_Layer named '%s' not found.", layername);
		return NULL;
	}
	if(NULL != layer->copy){
		PyErr_Format(PyExc_RuntimeError, "SetRegionPolygon: Cannot pattern a layer copy.");
		return NULL;
	}
	M = Simulation_GetMaterialByName(self->S, matname, &material_index);
	if(NULL == M){
		PyErr_Format(PyExc_RuntimeError, "SetRegionPolygon: S4_Material named '%s' not found.", matname);
		return NULL;
	}
	ret = Simulation_AddLayerPatternPolygon(self->S, layer, material_index, center, (M_PI/180.)*tilt, polydata.nvert, polydata.vert);
	free(polydata.vert);
	if(0 != ret){
		PyErr_Format(PyExc_MemoryError, "SetRegionPolygon: There was a problem allocating the pattern.");
		return NULL;
	}
	Py_RETURN_NONE;
}

static PyObject *S4Sim_SetExcitationExterior(S4Sim *self, PyObject *args, PyObject *kwds)
{
	static char *kwlist[] = {"Excitations", NULL};
	S4Excitation_Data exciData = {0, NULL, NULL};	//set exg or ex NULL to get the size of tuple.
	int err;
	/*double calls to one function with some flag to get size info,
	so heap location variable can be initialized. learn from WIN32 API :) */
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O&:S4Sim_SetExcitationExterior", kwlist, &excitation_converter, &exciData))
		return NULL;
	exciData.exg = (int*)malloc(sizeof(int) * 2 * exciData.n);
	exciData.ex = (double*)malloc(sizeof(double)* 2 * exciData.n);
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O&:S4Sim_SetExcitationExterior", kwlist, &excitation_converter, &exciData))
	{
		free(exciData.exg); exciData.exg = NULL;
		free(exciData.ex); exciData.ex = NULL;
		return NULL;
	}

	err = Simulation_MakeExcitationExterior(self->S, exciData.n, exciData.exg, exciData.ex);
	free(exciData.exg); exciData.exg = NULL;
	free(exciData.ex); exciData.ex = NULL;
	if(0 != err)
	{
		HandleSolutionErrorCode("S4Sim_SetExcitationExterior", err);
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

	pol_s[0] = sqrt(cs.real*cs.real + cs.imag*cs.imag); pol_s[1] = atan2(cs.imag,cs.real);
	pol_p[0] = sqrt(cp.real*cp.real + cp.imag*cp.imag); pol_p[1] = atan2(cp.imag,cp.real);
	angle[0] *= (M_PI/180.);
	angle[1] *= (M_PI/180.);
	ret = Simulation_MakeExcitationPlanewave(self->S, angle, pol_s, pol_p, order);
	if(0 != ret){
		HandleSolutionErrorCode("SetExcitationPlanewave", ret);
		return NULL;
	}
	Py_RETURN_NONE;
}

static PyObject *S4Sim_SetFrequency(S4Sim *self, PyObject *args){
	Py_complex f;
	if(!PyArg_ParseTuple(args, "D:SetFrequency", &f)){ return NULL; }

	Simulation_DestroySolution(self->S);

	self->S->omega[0] = 2*M_PI*f.real;
	self->S->omega[1] = 2*M_PI*f.imag;
	if(self->S->omega[0] <= 0){
		PyErr_Warn(PyExc_RuntimeWarning, "A non-positive frequency was specified.");
	}
	if(self->S->omega[1] > 0){
		PyErr_Warn(PyExc_RuntimeWarning, "A frequency with positive imaginary part was specified.");
	}
	Py_RETURN_NONE;
}

static PyObject *S4Sim_GetReciprocalLattice(S4Sim *self, PyObject *args){
	return PyTuple_Pack(2,
		PyTuple_Pack(2,
			PyFloat_FromDouble(self->S->Lk[0]), PyFloat_FromDouble(self->S->Lk[1])
		),
		PyTuple_Pack(2,
			PyFloat_FromDouble(self->S->Lk[2]), PyFloat_FromDouble(self->S->Lk[3])
		)
	);
}

static PyObject *S4Sim_GetEpsilon(S4Sim *self, PyObject *args){
	int ret;
	double r[3], feps[2];
	if(!PyArg_ParseTuple(args, "ddd:GetEpsilon", &r[0], &r[1], &r[2])){ return NULL; }
	ret = Simulation_GetEpsilon(self->S, r, feps);
	if(0 != ret){
		HandleSolutionErrorCode("GetEpsilon", ret);
	}
	return PyComplex_FromDoubles(feps[0], feps[1]);
}

static PyObject *S4Sim_OutputLayerPatternRealization(S4Sim *self, PyObject *args, PyObject *kwds)
{
	static char *kwlist[] = { "S4_Layer", "Nu", "Nv", "Filename", NULL };
	const char *layerName;
	const char *fileName = NULL;
	int Nu, Nv;
	S4_Layer *layer;
	FILE *fp;
	int err;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sii|s:OutputLayerPatternRealization", kwlist, &layerName, &Nu, &Nv, &fileName))
		return NULL;

	layer = Simulation_GetLayerByName(self->S, layerName, NULL);
	if(NULL == layer)
	{
		PyErr_Format(PyExc_RuntimeError, "OutputLayerPatternRealization: S4_Layer named '%s' not Found.", layerName);
		return NULL;
	}

	fp = stdout;
	if(NULL != fileName)
		fp = fopen(fileName, "wb");
	if(NULL == fp)
	{
		PyErr_Format(PyExc_IOError, "OutputLayerPatternRealization: open file '%s' failed.", fileName);
		return NULL;
	}

	err = Simulation_OutputLayerPatternRealization(self->S, layer, Nu, Nv, fp);
	if(0 != err)
	{
		HandleSolutionErrorCode("OutputLayerPatternRealization", err);
		return NULL;
	}
	if(NULL != fp)
		fclose(fp);

	Py_RETURN_NONE;
 }

static PyObject *S4Sim_OutputLayerPatternPostscript(S4Sim *self, PyObject *args, PyObject *kwds){
	int ret;
	static char *kwlist[] = { "S4_Layer", "Filename", NULL };
	const char *layername;
	const char *filename = NULL;
	S4_Layer *layer;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|s:OutputLayerPatternPostscript", kwlist, &layername, &filename)){ return NULL; }

	layer = Simulation_GetLayerByName(self->S, layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "OutputLayerPatternPostscript: S4_Layer named '%s' not found.", layername);
		return NULL;
	}
	{
		FILE *fp = stdout;
		if(NULL != filename){
			fp = fopen(filename, "wb");
		}

		ret = Simulation_OutputLayerPatternDescription(self->S, layer, fp);
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
	int ret;
	static char *kwlist[] = { "Filename", NULL };
	const char *filename = NULL;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "|s:OutputStructurePOVRay", kwlist, &filename)){ return NULL; }

	{
		FILE *fp = stdout;
		if(NULL != filename){
			fp = fopen(filename, "wb");
		}

		ret = Simulation_OutputStructurePOVRay(self->S, fp);
		if(0 != ret){
			HandleSolutionErrorCode("OutputStructurePOVRay", ret);
			return NULL;
		}

		if(NULL != filename){
			fclose(fp);
		}
	}
	Py_RETURN_NONE;
}

static PyObject *S4Sim_GetBasisSet(S4Sim *self, PyObject *args){
	int *G;
	int n, i, ret;
	PyObject *rv;

	ret = Simulation_InitSolution(self->S);
	if(0 != ret){
		HandleSolutionErrorCode("GetBasisSet", ret);
		return NULL;
	}

	n = Simulation_GetNumG(self->S, &G);
	if(NULL == G){
		HandleSolutionErrorCode("GetBasisSet", 0);
		return NULL;
	}
	rv = PyTuple_New(n);
	for(i = 0; i < n; ++i){
		PyTuple_SetItem(rv, i,
			PyTuple_Pack(2,
				FromIntPyDefInt(G[2*i+0]),
				FromIntPyDefInt(G[2*i+1])
			)
		);
	}
	return rv;
}

static PyObject *S4Sim_GetAmplitudes(S4Sim *self, PyObject *args, PyObject *kwds){
	int ret, n, i, j;
	int *G;
	static char *kwlist[] = { "S4_Layer", "zOffset", NULL };
	const char *layername;
	double offset = 0;
	double *amp;
	S4_Layer *layer;
	PyObject *rv, *rventry;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|d:GetAmplitudes", kwlist, &layername, &offset)){ return NULL; }

	layer = Simulation_GetLayerByName(self->S, layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "GetAmplitudes: S4_Layer named '%s' not found.", layername);
		return NULL;
	}
	n = Simulation_GetNumG(self->S, &G);
	amp = (double*)malloc(sizeof(double)*8*n);
	ret = Simulation_GetAmplitudes(self->S, layer, offset, amp, &amp[4*n]);
	if(0 != ret){
		HandleSolutionErrorCode("GetAmplitudes", ret);
		return NULL;
	}

	rv = PyTuple_New(2);
	for(j = 0; j < 2; ++j){
		rventry = PyTuple_New(2*n);
		PyTuple_SetItem(rv, j, rventry);
		for(i = 0; i < 2*n; ++i){
			PyTuple_SetItem(rventry, i,
				PyComplex_FromDoubles(amp[4*n*j+2*i+0], amp[4*n*j+2*i+1])
			);
		}
	}
	free(amp);
	return rv;
}
static PyObject *S4Sim_GetPowerFlux(S4Sim *self, PyObject *args, PyObject *kwds){
	int ret;
	static char *kwlist[] = { "S4_Layer", "zOffset", NULL };
	const char *layername;
	double offset = 0;
	double power[4];
	S4_Layer *layer;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|d:GetPowerFlux", kwlist, &layername, &offset)){ return NULL; }

	layer = Simulation_GetLayerByName(self->S, layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "GetPowerFlux: S4_Layer named '%s' not found.", layername);
		return NULL;
	}
	ret = Simulation_GetPoyntingFlux(self->S, layer, offset, power);
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
	int ret, n, i;
	int *G;
	static char *kwlist[] = { "S4_Layer", "zOffset", NULL };
	const char *layername;
	double offset = 0;
	double *power;
	S4_Layer *layer;
	PyObject *rv;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|d:GetPowerFluxByOrder", kwlist, &layername, &offset)){ return NULL; }

	layer = Simulation_GetLayerByName(self->S, layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "GetPowerFluxByOrder: S4_Layer named '%s' not found.", layername);
		return NULL;
	}
	n = Simulation_GetNumG(self->S, &G);
	power = (double*)malloc(sizeof(double)*4*n);
	ret = Simulation_GetPoyntingFluxByG(self->S, layer, offset, power);
	if(0 != ret){
		HandleSolutionErrorCode("GetPowerFluxByOrder", ret);
		return NULL;
	}

	rv = PyTuple_New(n);
	for(i = 0; i < n; ++i){
		PyTuple_SetItem(rv, i,
			PyTuple_Pack(2,
				PyComplex_FromDoubles(power[4*i+0], power[4*i+2]),
				PyComplex_FromDoubles(power[4*i+1], power[4*i+3])
			)
		);
	}
	free(power);
	return rv;
}
static PyObject *S4Sim_GetStressTensorIntegral(S4Sim *self, PyObject *args, PyObject *kwds){
	int ret;
	static char *kwlist[] = { "S4_Layer", "zOffset", NULL };
	const char *layername;
	double offset = 0;
	double Tint[6];
	S4_Layer *layer;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|d:GetStressTensorIntegral", kwlist, &layername, &offset)){ return NULL; }

	layer = Simulation_GetLayerByName(self->S, layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "GetStressTensorIntegral: S4_Layer named '%s' not found.", layername);
		return NULL;
	}
	ret = Simulation_GetStressTensorIntegral(self->S, layer, offset, Tint);
	if(0 != ret){
		HandleSolutionErrorCode("GetStressTensorIntegral", ret);
		return NULL;
	}

	return PyTuple_Pack(3,
		PyComplex_FromDoubles(Tint[0], Tint[3]),
		PyComplex_FromDoubles(Tint[1], Tint[4]),
		PyComplex_FromDoubles(Tint[2], Tint[5])

	);
}
static PyObject *S4Sim_GetLayerVolumeIntegral(S4Sim *self, PyObject *args, PyObject *kwds){
	int ret;
	static char *kwlist[] = { "S4_Layer", "Quantity", NULL };
	const char *layername;
	const char *strwhat;
	double integral[2];
	char which = 'U';
	S4_Layer *layer;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss:GetLayerVolumeIntegral", kwlist, &layername, &strwhat)){ return NULL; }

	layer = Simulation_GetLayerByName(self->S, layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "GetLayerVolumeIntegral: S4_Layer named '%s' not found.", layername);
		return NULL;
	}

	if(0 == strcmp("U", strwhat)){
		which = 'U';
	}else if(0 == strcmp("E", strwhat)){
		which = 'E';
	}else if(0 == strcmp("H", strwhat)){
		which = 'H';
	}else if(0 == strcmp("e", strwhat)){
		which = 'e';
	}else{
		PyErr_Format(PyExc_RuntimeError, "GetLayerVolumeIntegral: Quantity must be one of: 'U', 'E', 'H', 'e'");
		return NULL;
	}

	ret = Simulation_GetLayerVolumeIntegral(self->S, layer, which, integral);
	if(0 != ret){
		HandleSolutionErrorCode("GetLayerVolumeIntegral", ret);
		return NULL;
	}

	return PyComplex_FromDoubles(integral[0], integral[1]);
}
static PyObject *S4Sim_GetLayerZIntegral(S4Sim *self, PyObject *args, PyObject *kwds){
	int ret;
	static char *kwlist[] = { "S4_Layer", "xy", NULL };
	const char *layername;
	double integral[6], r[2];
	S4_Layer *layer;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s(dd):GetLayerZIntegral", kwlist, &layername, &r[0], &r[1])){ return NULL; }

	layer = Simulation_GetLayerByName(self->S, layername, NULL);
	if(NULL == layer){
		PyErr_Format(PyExc_RuntimeError, "GetLayerZIntegral: S4_Layer named '%s' not found.", layername);
		return NULL;
	}

	ret = Simulation_GetLayerZIntegral(self->S, layer, r, integral);
	if(0 != ret){
		HandleSolutionErrorCode("GetLayerZIntegral", ret);
		return NULL;
	}
	return PyTuple_Pack(2,
		PyTuple_Pack(3,
			PyFloat_FromDouble(integral[0]),
			PyFloat_FromDouble(integral[1]),
			PyFloat_FromDouble(integral[2])
		),
		PyTuple_Pack(3,
			PyFloat_FromDouble(integral[3]),
			PyFloat_FromDouble(integral[4]),
			PyFloat_FromDouble(integral[5])
		)
	);
}
static PyObject *S4Sim_GetFields(S4Sim *self, PyObject *args, PyObject *kwds){
	int ret;
	double r[3], fE[6],fH[6];
	if(!PyArg_ParseTuple(args, "ddd:GetFields", &r[0], &r[1], &r[2])){ return NULL; }

	ret = Simulation_GetField(self->S, r, fE, fH);
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
	static char *kwlist[] = { "z", "NumSamples", "Format", "BaseFilename", NULL };
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
	ret = Simulation_GetFieldPlane(self->S, snxy, z, Efields, Hfields);
	if(0 != ret){
		HandleSolutionErrorCode("GetFieldsOnGrid", ret);
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

static PyObject *S4Sim_GetSMatrixDeterminant(S4Sim *self, PyObject *args){
	int ret;
	double mant[2], base;
	int expo;

	ret = Simulation_GetSMatrixDeterminant(self->S, mant, &base, &expo);
	if(0 != ret){
		HandleSolutionErrorCode("GetSMatrixDeterminant", ret);
		return NULL;
	}
	return PyTuple_Pack(3,
		PyComplex_FromDoubles(mant[0], mant[1]),
		PyFloat_FromDouble(base),
		FromIntPyDefInt(expo)
	);
}

static PyObject *S4Sim_SetVerbosity(S4Sim *self, PyObject *args, PyObject *kwds)
{
	static char *kwlist[] = {"Level", NULL};
	int level;
	if(PyArg_ParseTupleAndKeywords(args, kwds, "i:SetVerbosity", kwlist, &level))
		return NULL;

	if(level < 0 || level > 9)
	{
		PyErr_SetString(PyExc_TypeError, "the level should be specifiled in [0, 9].");
		return NULL;
	}
	self->S->options.verbosity = level;
	Py_RETURN_NONE;
}

static PyObject *S4Sim_SetOptions(S4Sim *self, PyObject *args, PyObject *kwds){
	static char *kwlist[] = {
		"Verbosity",                 /* int */
		"LatticeTruncation",         /* str */
		"DiscretizedEpsilon",        /* bool */
		"DiscretizationResolution",  /* int */
		"PolarizationDecomposition", /* bool */
		"PolarizationBasis",         /* str */
		"LanczosSmoothing",          /* obj */
		"SubpixelSmoothing",         /* bool */
		"ConserveMemory",            /* bool */
		NULL
	};
	int verbosity = -1;
	int discretized_epsilon = -1;
	int discretization_resolution = -1;
	int polarization_decomp = -1;
	struct lanczos_smoothing_settings lanczos_smoothing;
	int subpixel_smoothing = -1;
	int conserve_memory = -1;

	lanczos_smoothing.set = 0;

	const char *lattice_truncation = NULL;
	const char *polarization_basis = NULL;

	if(!PyArg_ParseTupleAndKeywords(
		args, kwds, "|isO&iO&sO&O&O&:SetOptions", kwlist,
		&verbosity,
		&lattice_truncation,
		&bool_converter, &discretized_epsilon,
		&discretization_resolution,
		&bool_converter, &polarization_decomp,
		&polarization_basis,
		&bool_converter, &lanczos_smoothing,
		&bool_converter, &subpixel_smoothing,
		&bool_converter, &conserve_memory
	)){ return NULL; }
	if(verbosity >= 0){
		if(verbosity > 9){ verbosity = 9; }
		self->S->options.verbosity = verbosity;
	}
	if(NULL != lattice_truncation){
		if(0 == strcmp("Circular", lattice_truncation)){
			self->S->options.lattice_truncation = 0;
		}else if(0 == strcmp("Parallelogramic", lattice_truncation)){
			self->S->options.lattice_truncation = 1;
		}else{
			PyErr_SetString(PyExc_ValueError, "LatticeTruncation must be one of: 'Circular', 'Parallelogramic'");
			return NULL;
		}
	}
	if(discretized_epsilon >= 0){
		self->S->options.use_discretized_epsilon = discretized_epsilon;
	}
	if(discretization_resolution >= 0){
		if(discretization_resolution < 2){
			PyErr_SetString(PyExc_ValueError, "DiscretizationResolution must be at least 2");
			return NULL;
		}
		self->S->options.resolution = discretization_resolution;
	}
	if(polarization_decomp >= 0){
		self->S->options.use_polarization_basis = polarization_decomp;
	}
	if(NULL != polarization_basis){
		if(0 == strcmp("Default", polarization_basis)){
			self->S->options.use_normal_vector_basis = 0;
			self->S->options.use_jones_vector_basis = 0;
		}else if(0 == strcmp("Normal", polarization_basis)){
			self->S->options.use_normal_vector_basis = 1;
			self->S->options.use_jones_vector_basis = 0;
		}else if(0 == strcmp("Jones", polarization_basis)){
			self->S->options.use_normal_vector_basis = 0;
			self->S->options.use_jones_vector_basis = 1;
		}else{
			PyErr_SetString(PyExc_ValueError, "PolarizationBasis must be one of: 'Default', 'Normal', 'Jones'");
			return NULL;
		}
	}
	if(lanczos_smoothing.set){
		self->S->options.use_Lanczos_smoothing = lanczos_smoothing.use;
		if(lanczos_smoothing.set_power){
			self->S->options.lanczos_smoothing_power = lanczos_smoothing.power;
		}
		if(lanczos_smoothing.set_width){
			self->S->options.lanczos_smoothing_width = lanczos_smoothing.width;
		}
	}
	if(subpixel_smoothing >= 0){
		self->S->options.use_subpixel_smoothing = subpixel_smoothing;
	}
	if(conserve_memory >= 0){
		self->S->options.use_less_memory = conserve_memory;
	}
	Py_RETURN_NONE;
}

static PyObject *S4SpectrumSampler_IsDone(S4SpectrumSampler *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args, ":IsDone"))
		return NULL;
	if(SpectrumSampler_IsDone(self->SpecS))
		Py_RETURN_TRUE;
	Py_RETURN_FALSE;
}

static PyObject *S4SpectrumSampler_IsParallelized(S4SpectrumSampler *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args, ":IsParallelized"))
		return NULL;
	if(SpectrumSampler_IsParallelized(self->SpecS))
		Py_RETURN_TRUE;
	Py_RETURN_FALSE;
}

static PyObject *S4SpectrumSampler_GetFrequency(S4SpectrumSampler *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args, ":GetFrequency"))
		return NULL;
	if(SpectrumSampler_IsParallelized(self->SpecS))
	{
		PyErr_SetString(PyExc_RuntimeError, "call 'GetFrequency' is illegal when is parallelized.");
		return NULL;
	}
	return Py_BuildValue("d", SpectrumSampler_GetFrequency(self->SpecS));
}

static PyObject *S4SpectrumSampler_GetFrequencies(S4SpectrumSampler *self, PyObject *args)
{
	Py_ssize_t nf;
	const double *freqs;
	PyObject *retObj;
	if(!PyArg_ParseTuple(args, ":GetFrequencies"))
		return NULL;
	if(!SpectrumSampler_IsParallelized(self->SpecS))
	{
		PyErr_SetString(PyExc_RuntimeError, "call 'GetFrequencies' is illegal when is not parallelized.");
		return NULL;
	}
	nf = SpectrumSampler_GetFrequencies(self->SpecS, &freqs);
	retObj = PyTuple_New(nf);
	for(int i = 0; i < nf; i++)
		PyTuple_SetItem(retObj, i, Py_BuildValue("d", freqs[i]));
	return retObj;
}

static PyObject *S4SpectrumSampler_SubmitResult(S4SpectrumSampler *self, PyObject *args, PyObject *kwds)
{
	static char *kwlist[] = {"Result", NULL};
	double y;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "d:SubmitResult", kwlist, &y))
		return NULL;
	if(SpectrumSampler_IsParallelized(self->SpecS))
	{
		PyErr_SetString(PyExc_RuntimeError, "call 'SubmitResult' is illegal when is parallelized.");
		return NULL;
	}
	SpectrumSampler_SubmitResult(self->SpecS, y);
	Py_RETURN_NONE;
}

static PyObject *S4SpectrumSampler_SubmitResults(S4SpectrumSampler *self, PyObject *args, PyObject *kwds)
{
	static char *kwlist[] = { "Results", NULL };
	int ny;
	double *y;
	PyObject *tupleObj = NULL;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!:SubmitResults", kwlist, &PyTuple_Type, &tupleObj))
		return NULL;
	if(!SpectrumSampler_IsParallelized(self->SpecS))
	{
		PyErr_SetString(PyExc_RuntimeError, "call 'SubmitResult' is illegal when is not parallelized.");
		return NULL;
	}
	ny = SpectrumSampler_GetSubmissionBuffer(self->SpecS, &y);
	if(ny != PyTuple_Size(tupleObj))
	{
		PyErr_SetString(PyExc_TypeError, "the length of the results is not equal to the buffer.");
		return NULL;
	}
	for(int i = 0; i < ny; i++)
		y[i] = PyFloat_AsDouble(PyTuple_GetItem(tupleObj, i));
	SpectrumSampler_SubmitResults(self->SpecS);
	Py_RETURN_NONE;
}

static PyObject *S4SpectrumSampler_GetSpectrum(S4SpectrumSampler *self, PyObject *args)
{
	int n;
	double pt[2];
	PyObject *retObj;
	SpectrumSampler_Enumerator e;
	if(!PyArg_ParseTuple(args, ":GetSpectrum"))
		return NULL;
	n = SpectrumSampler_GetNumPoints(self->SpecS);
	retObj = PyTuple_New(n);
	e = SpectrumSampler_GetPointEnumerator(self->SpecS);
	for(int i = 0; i < n; i++)
	{
		SpectrumSampler_Enumerator_Get(e, pt);
		PyTuple_SetItem(retObj, i, PyTuple_Pack(2, pt[0], pt[1]));
	}
	return retObj;
}

static PyMethodDef S4SpectrumSampler_methods[] =
{
	{"IsDone"			, (PyCFunction)S4SpectrumSampler_IsDone, METH_VARARGS, PyDoc_STR("IsDone() -> bool")},
	{"IsParallelized"	, (PyCFunction)S4SpectrumSampler_IsParallelized, METH_VARARGS, PyDoc_STR("IsParallelized() -> bool")},
	{"GetFrequency"		, (PyCFunction)S4SpectrumSampler_GetFrequency, METH_VARARGS, PyDoc_STR("GetFrequency() -> freq")},
	{"GetFrequencies"	, (PyCFunction)S4SpectrumSampler_GetFrequencies, METH_VARARGS, PyDoc_STR("GetFrequencies() -> tuple")},
	{"SubmitResult"		, (PyCFunction)S4SpectrumSampler_SubmitResult, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SubmitResult(result) -> None")},
	{"SubmitResults"	, (PyCFunction)S4SpectrumSampler_SubmitResults, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SubmitResults(results) -> None")},
	{"GetSpectrum"		, (PyCFunction)S4SpectrumSampler_GetSpectrum, METH_VARARGS, PyDoc_STR("GetSpectrum() -> tuple")},
	{NULL, NULL, 0, NULL}
};

static PyMethodDef	S4Interpolator_methods[] =
{
	{ "Get", (PyCFunction)S4Interpolator_Get, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("Get(x) -> Tuple") },

	{ NULL, NULL, 0, NULL }
};

static PyMethodDef S4Sim_methods[] = {
	{"Clone"            , (PyCFunction)S4Sim_Clone, METH_NOARGS, PyDoc_STR("Clone() -> S4.S4_Simulation")},
	/* Specification */
	{"AddMaterial"				, (PyCFunction)S4Sim_AddMaterial, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("AddMaterial(name,eps) -> None")},
	{"SetMaterial"				, (PyCFunction)S4Sim_SetMaterial, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetMaterial(name,eps) -> None")},
	{"AddLayer"					, (PyCFunction)S4Sim_AddLayer, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("AddLayer(name,thickness,matname) -> None")},
	{"AddLayerCopy"				, (PyCFunction)S4Sim_AddLayerCopy, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("AddLayerCopy(name,thickness,layer) -> None")},
	{"SetLayer"					, (PyCFunction)S4Sim_SetLayer, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetLayer(name,thickness,material) -> None")},
	{"SetLayerThickness"		, (PyCFunction)S4Sim_SetLayerThickness, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetLayerThickness(layer,thickness) -> None")},
	{"SetVerbosity"				, (PyCFunction)S4Sim_SetVerbosity, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetVerbosity(level) -> None")},
	{"RemoveLayerRegions"		, (PyCFunction)S4Sim_RemoveLayerRegions, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("RemoveLayerPatterns(layer) -> None")},
	{"SetRegionCircle"			, (PyCFunction)S4Sim_SetRegionCircle, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetLayerPatternCircle(layer,matname,center,radius) -> None")},
	{"SetRegionEllipse"			, (PyCFunction)S4Sim_SetRegionEllipse, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetLayerPatternEllipse(layer,matname,center,angle,halfwidths) -> None")},
	{"SetRegionRectangle"		, (PyCFunction)S4Sim_SetRegionRectangle, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetLayerPatternRectangle(layer,matname,center,angle,halfwidths) -> None")},
	{"SetRegionPolygon"			, (PyCFunction)S4Sim_SetRegionPolygon, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetLayerPatternPolygon(layer,matname,center,angle,vertices) -> None")},
	{"SetExcitationPlanewave"	, (PyCFunction)S4Sim_SetExcitationPlanewave, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetExcitationPlanewave(angles,s_amp,p_amp) -> None")},
	{"SetExcitationExterior"	, (PyCFunction)S4Sim_SetExcitationExterior, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetExcitationExterior(Excitations) -> None")},
	{"SetFrequency"				, (PyCFunction)S4Sim_SetFrequency, METH_VARARGS, PyDoc_STR("SetFrequency(freq) -> None")},
	/* Outputs requiring no solutions */
	{"GetReciprocalLattice"		, (PyCFunction)S4Sim_GetReciprocalLattice, METH_NOARGS, PyDoc_STR("GetReciprocalLattice() -> ((px,py),(qx,qy))")},
	{"GetEpsilon"				, (PyCFunction)S4Sim_GetEpsilon, METH_VARARGS, PyDoc_STR("GetEpsilon(x,y,z) -> Complex")},
	{"OutputLayerPatternPostscript", (PyCFunction)S4Sim_OutputLayerPatternPostscript, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("OutputLayerPatternPostscript(layer,filename) -> None")},
	{ "OutputLayerPatternRealization", (PyCFunction)S4Sim_OutputLayerPatternRealization, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("OutputLayerPatternRealization(layer, nu, nv, filename) -> None")},
	/* Outputs requiring solutions */
	{"OutputStructurePOVRay"	, (PyCFunction)S4Sim_OutputStructurePOVRay, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("OutputStructurePOVRay(filename) -> None")},
	{"GetBasisSet"				, (PyCFunction)S4Sim_GetBasisSet, METH_NOARGS, PyDoc_STR("GetBasisSet() -> Tuple")},
	{"GetAmplitudes"			, (PyCFunction)S4Sim_GetAmplitudes, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetAmplitudes(layer,zoffset) -> Tuple")},
	{"GetPowerFlux"				, (PyCFunction)S4Sim_GetPowerFlux, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetPowerFlux(layer,zoffset) -> (forw,back)")},
	{"GetPowerFluxByOrder"		, (PyCFunction)S4Sim_GetPowerFluxByOrder, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetPowerFluxByOrder(layer,zoffset) -> Tuple")},
	{"GetStressTensorIntegral"	, (PyCFunction)S4Sim_GetStressTensorIntegral, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetStressTensorIntegral(layer,zoffset) -> Complex")},
	{"GetLayerVolumeIntegral"	, (PyCFunction)S4Sim_GetLayerVolumeIntegral, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetLayerVolumeIntegral(layer,which) -> Complex")},
	{"GetLayerZIntegral"		, (PyCFunction)S4Sim_GetLayerZIntegral, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetLayerZIntegral(layer,which,pos) -> Complex")},
	/*
	{"GetEField"				, (PyCFunction)S4Sim_GetEField, METH_VARARGS, PyDoc_STR("GetEField(x,y,z) -> (Tuple)")},
	{"GetHField"				, (PyCFunction)S4Sim_GetHField, METH_VARARGS, PyDoc_STR("GetHField(x,y,z) -> (Tuple)")},
	*/
	{"GetFields"				, (PyCFunction)S4Sim_GetFields, METH_VARARGS, PyDoc_STR("GetFields(x,y,z) -> (Tuple,Tuple)")},
	{"GetFieldsOnGrid"			, (PyCFunction)S4Sim_GetFieldsOnGrid, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetFieldsOnGrid(z,nsamples,format,filename) -> Tuple")},
	{"GetSMatrixDeterminant"	, (PyCFunction)S4Sim_GetSMatrixDeterminant, METH_NOARGS, PyDoc_STR("GetSMatrixDeterminant() -> Tuple")},
	/*
	{"GetDiffractionOrder"		, (PyCFunction)S4Sim_GetDiffractionOrder, METH_VARARGS, PyDoc_STR("GetDiffractionOrder(m,n) -> order")},
	*/
	{"SetOptions"				, (PyCFunction)S4Sim_SetOptions, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetOptions() -> None")},
	{"GetPoyntingFlux"			, (PyCFunction)S4Sim_GetPowerFlux, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetPoyntingFlux(layer,zoffset) -> (forw, back)")},
	{"GetPoyntingFluxByOrder"	, (PyCFunction)S4Sim_GetPowerFluxByOrder, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("GetPoyntingFluxByOrder(layer,zoffset) -> Tuple")},
	/*
	{"GetGList"					, (PyCFunction)S4Sim_GetGList, METH_VARARGS, PyDoc_STR("GetGList() -> Tuple")},
	{"GetNumG"					, (PyCFunction)S4Sim_GetNumG, METH_VARARGS, PyDoc_STR("GetNumG() -> G num")},
	*/
	/*options*/
	/*
	{"SetBasisFieldDumpPrefix"	, (PyCFunction)S4Sim_SetBasisFieldDumpPrefix, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetBasisFieldDumpPrefix(prefix) -> None")},
	{"SetLatticeTruncation"		, (PyCFunction)S4Sim_SetLatticeTruncation, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SetLatticeTruncation(Trunc) -> NOne")},
	*/
	{NULL, NULL}
};

static PyTypeObject S4Sim_Type = {
	/* The ob_type field must be initialized in the module init function
	 * to be portable to Windows without using C++. */
	PyVarObject_HEAD_INIT(NULL, 0)
	"S4.S4_Simulation",    /*tp_name*/
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

static PyTypeObject S4Interpolator_Type = {
	/* The ob_type field must be initialized in the module init function
	* to be portable to Windows without using C++. */
	PyVarObject_HEAD_INIT(NULL, 0)
	"S4.Interpolator",    /*tp_name*/
	sizeof(S4Interpolator),      /*tp_basicsize*/
	0,                  /*tp_itemsize*/
	/* methods */
	(destructor)S4Interpolator_dealloc, /*tp_dealloc*/
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
	S4Interpolator_methods,      /*tp_methods*/
	0,                  /*tp_members*/
	0,                  /*tp_getset*/
	0,                  /*tp_base*/
	0,                  /*tp_dict*/
	0,                  /*tp_descr_get*/
	0,                  /*tp_descr_set*/
	0,                  /*tp_dictoffset*/
	0,                  /*tp_init*/
	0,                  /*tp_alloc*/
	S4Interpolator_new,          /*tp_new*/
	0,                  /*tp_free*/
	0,                  /*tp_is_gc*/
};

static PyTypeObject S4SpectrumSampler_Type = {
	/* The ob_type field must be initialized in the module init function
	* to be portable to Windows without using C++. */
	PyVarObject_HEAD_INIT(NULL, 0)
	"S4.SpectrumSampler",    /*tp_name*/
	sizeof(S4SpectrumSampler),      /*tp_basicsize*/
	0,                  /*tp_itemsize*/
	/* methods */
	(destructor)S4SpectrumSampler_dealloc, /*tp_dealloc*/
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
	S4SpectrumSampler_methods,      /*tp_methods*/
	0,                  /*tp_members*/
	0,                  /*tp_getset*/
	0,                  /*tp_base*/
	0,                  /*tp_dict*/
	0,                  /*tp_descr_get*/
	0,                  /*tp_descr_set*/
	0,                  /*tp_dictoffset*/
	0,                  /*tp_init*/
	0,                  /*tp_alloc*/
	S4SpectrumSampler_new,          /*tp_new*/
	0,                  /*tp_free*/
	0,                  /*tp_is_gc*/
};

static PyObject *S4_new(PyObject *self, PyObject *args, PyObject *kwds){
	return (PyObject*)S4Sim_new(&S4Sim_Type, args, kwds);
}

static PyObject *S4_NewInterpolator(PyObject *self, PyObject *args, PyObject *kwds)
{
	return (PyObject*)S4Interpolator_new(&S4Interpolator_Type, args, kwds);
}

//didn't finished yet
static PyObject *S4_SolveInParallel(PyObject *Self, PyObject *args, PyObject *kwds)
{
	static char *kwlist[] = {"S4_Layer", "Simulations", NULL};
	const char *layerName;
	//S4_solve_in
	Py_RETURN_NONE;
}

static PyMethodDef S4_funcs[] = {
	{"New"				, (PyCFunction)S4_new, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("New() -> new S4 simulation object")},
	{"SolveInParallel"	, (PyCFunction)S4_SolveInParallel, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("SolveInParallel(layer, sim_obj) -> None")},
	{ "NewInterpolator"	, (PyCFunction)S4_NewInterpolator, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("NewInterpolator(type, table) -> new S4 interpolator object") },
	//{"PrintTuple"		, (PyCFunction)S4_PrintTuple, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("PrintTuple(tuple) -> None")},
	{ "NewSpectrumSampler", (PyCFunction)S4_NewSpectrumSampler, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("NewSpectrumSampler() -> new S4 spectrum sampler onject")},
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
	if(PyType_Ready(&S4Interpolator_Type) < 0){ INITERROR; }
	if(PyType_Ready(&S4SpectrumSampler_Type) < 0){ INITERROR; }

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
