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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <stdlib.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include "S4.h"
#include "SpectrumSampler.h"
#include "cubature.h"
#include "Interpolator.h"
#include "convert.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void fft_init();
void fft_destroy();

#ifdef S4_DEBUG
# include "debug.h"
# ifdef HAVE_LIBPTHREAD
pthread_mutex_t g_mutex = PTHREAD_MUTEX_INITIALIZER;
int count = 0;
pthread_cond_t g_cond = PTHREAD_COND_INITIALIZER;
# endif
#endif

static char interactive_key = 'k';
static char thread_key = 't';
static char thread_count_key = 'c';

unsigned int get_max_threads();
void S4_threads_init(lua_State *L, unsigned int nthreads);
void S4_threads_destroy(lua_State *L);
void S4_threads_join(lua_State *L, int n);
void S4_threads_run(lua_State *L, int i, void* (*func)(void*), void *data);

#ifdef HAVE_LIBPTHREAD
#define PTW32_STATIC_LIB
#include <pthread.h>

#ifdef WIN32
# define WIN32_LEAN_AND_MEAN
# include <windows.h>
#elif defined(__APPLE__)
# include <sys/sysctl.h>
#else
# ifdef HAVE_UNISTD_H
#  include <unistd.h>
# endif
#endif

unsigned int get_max_threads(){
	unsigned int max_threads;
#ifdef WIN32
	DWORD ProcessorCount; 
	SYSTEM_INFO SystemInfo; 
	GetSystemInfo(&SystemInfo); 
	ProcessorCount = SystemInfo.dwNumberOfProcessors;
	max_threads = ProcessorCount;
#elif defined(__APPLE__)
	int cpuCount;
	size_t countSize = sizeof(cpuCount);
	sysctlbyname("hw.logicalcpu", &cpuCount, &countSize, NULL, 0);
	max_threads = cpuCount;
#else
	int nprocs = sysconf(_SC_NPROCESSORS_ONLN);
	max_threads = (nprocs>=1)?nprocs:1;
#endif
	return max_threads;
}

void S4_threads_init(lua_State *L, unsigned int nthreads){
#ifdef WIN32
	if(!pthread_win32_process_attach_np()){
		fprintf(stderr, "Pthread initialization failed.\n");
		return;
	}
#endif
	pthread_t *thread;
	lua_pushlightuserdata(L, (void *)&thread_count_key);
	lua_pushinteger(L, nthreads);
	lua_settable(L, LUA_REGISTRYINDEX);
	
	thread = (pthread_t*)malloc(sizeof(pthread_t)*nthreads);
	lua_pushlightuserdata(L, (void *)&thread_key);
	lua_pushlightuserdata(L, (void*)thread);
	lua_settable(L, LUA_REGISTRYINDEX);
}
void S4_threads_destroy(lua_State *L){
	int  nthreads;
	pthread_t *thread;
	lua_pushlightuserdata(L, (void *)&thread_count_key);
	lua_gettable(L, LUA_REGISTRYINDEX);
	nthreads = lua_tointeger(L, -1);
	lua_pop(L, 1);
	lua_pushlightuserdata(L, (void *)&thread_key);
	lua_gettable(L, LUA_REGISTRYINDEX);
	thread = (pthread_t*)lua_topointer(L, -1);
	lua_pop(L, 1);
	free(thread);
#ifdef WIN32
	pthread_win32_process_detach_np();
#endif
}
void S4_threads_join(lua_State *L, int n){
	int i, nthreads;
	pthread_t *thread;
	lua_pushlightuserdata(L, (void *)&thread_count_key);
	lua_gettable(L, LUA_REGISTRYINDEX);
	nthreads = lua_tointeger(L, -1);
	lua_pop(L, 1);
	lua_pushlightuserdata(L, (void *)&thread_key);
	lua_gettable(L, LUA_REGISTRYINDEX);
	thread = (pthread_t*)lua_topointer(L, -1);
	lua_pop(L, 1);
	if(n < nthreads){ nthreads = n; }
	for(i = 0; i < nthreads; ++i){
		pthread_join(thread[i], NULL);
	}
}
void S4_threads_run(lua_State *L, int i, void* (*func)(void*), void *data){
	int nthreads;
	pthread_t *thread;
	lua_pushlightuserdata(L, (void *)&thread_count_key);
	lua_gettable(L, LUA_REGISTRYINDEX);
	nthreads = lua_tointeger(L, -1);
	lua_pop(L, 1);
	lua_pushlightuserdata(L, (void *)&thread_key);
	lua_gettable(L, LUA_REGISTRYINDEX);
	thread = (pthread_t*)lua_topointer(L, -1);
	lua_pop(L, 1);
	
	if(i >= nthreads){
		i = i%nthreads;
		pthread_join(thread[i], NULL);
	}
	pthread_create(&thread[i], NULL, func, data);
}
#else /* HAVE_LIBPTHREAD */
unsigned int get_max_threads(){ return 1; }
void S4_threads_init(lua_State *L, unsigned int nthreads){}
void S4_threads_destroy(lua_State *L){}
void S4_threads_join(lua_State *L, int n){}
void S4_threads_run(lua_State *L, int i, void* (*func)(void*), void *data){
	func(data);
}
#endif /* HAVE_LIBPTHREAD */

void S4L_set_interactive(lua_State *L, int i){
	lua_pushlightuserdata(L, (void *)&interactive_key);
	lua_pushnumber(L, i);
	lua_settable(L, LUA_REGISTRYINDEX);
}
int S4L_get_interactive(lua_State *L){
	int n;
	lua_pushlightuserdata(L, (void *)&interactive_key);
	lua_gettable(L, LUA_REGISTRYINDEX);
	n = lua_tointeger(L, -1);
	lua_pop(L, 1);
	return n;
}

void S4L_error(lua_State *L, const char *fmt, ...){
	va_list argp;
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	fputc('\n', stderr);
	va_end(argp);
	
	if(!S4L_get_interactive(L)){
		lua_close(L);
		exit(EXIT_FAILURE);
	}
}

void HandleSolutionErrorCode(lua_State *L, const char *fname, int code){
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
		S4L_error(L, "%s: %s.", fname, str);
	}else{
		S4L_error(L, "%s: %s. Error code: %d", fname, str, code);
	}
}

typedef struct S4_solve_in_parallel_data_{
	Simulation *S;
	Layer *layer;
} S4_solve_in_parallel_data;

void* S4_solve_in_parallel(void *p){
	S4_solve_in_parallel_data *data = (S4_solve_in_parallel_data*)p;
	Simulation_SolveLayer(data->S, data->layer);
	return NULL;
}
/* Expected stack:
 *   1 layer name string
 *   2 Simulation object
 *   : Simulation objects
 */
static int S4L_SolveInParallel(lua_State *L){
	int i, n;
	S4_solve_in_parallel_data *data;
	const char *layer_name = luaL_checklstring(L, 1, NULL);
	
	n = lua_gettop(L);
	data = (S4_solve_in_parallel_data*)malloc(sizeof(S4_solve_in_parallel_data)*(n-1));
	for(i = 2; i <= n; ++i){
		Layer *layer;
		Simulation *S = (Simulation *)luaL_checkudata(L, i, "S4.Simulation");
		luaL_argcheck(L, S != NULL, i, "SolveInParallel: 'Simulation' object expected.");
		layer = Simulation_GetLayerByName(S, layer_name, NULL);
		if(NULL == layer){
			S4L_error(L, "SolveInParallel: Layer named '%s' not found.", layer_name);
		}
		data[i-2].S = S;
		data[i-2].layer = layer;

		S4_threads_run(L, i-2, &S4_solve_in_parallel, (void*)&data[i-2]);
	}
	S4_threads_join(L, n-1);
	free(data);
	
	return 0;
}

static int S4L_ConvertUnits(lua_State *L){
	double value;
	const char *from_units;
	const char *to_units;
	
	value = luaL_checknumber(L, 1);
	from_units = luaL_checklstring(L, 2, NULL);
	to_units = luaL_checklstring(L, 3, NULL);
	
	if(0 == convert_units(&value, from_units, to_units)){
		lua_pushnumber(L, value);
		return 1;
	}else{
		return 0;
	}
}

static void IntegrateFunction(
	unsigned ndim, const double *x, void *fdata,
	unsigned fdim, double *fval
){
	unsigned i;
	lua_State *L = (lua_State*)fdata;
	
	fval[0] = 0;
	
	lua_pushvalue(L, -1); /* push a copy of the function */
	
	for(i = 0; i < ndim; ++i){
		lua_pushnumber(L, x[i]);
	}
	if(0 != lua_pcall(L, ndim, fdim, 0)){
		S4L_error(L, "Error running Integrate function: %s", lua_tostring(L, -1));
        lua_pop(L, 1);
        return;
	}
	
	if(!lua_isnumber(L, -1)){
		lua_pop(L, 1);
		S4L_error(L, "Integrate function must return a number");
		return;
	}
	
	fval[0] = lua_tonumber(L, -1);
	lua_pop(L, 1);
}
static void IntegrateFunctionVectorized(
	unsigned ndim, unsigned npts, const double *x, void *fdata,
	unsigned fdim, double *fval
){
	lua_State *L = (lua_State*)fdata;
	unsigned i, j;
	unsigned table_size;
	
	if(fdim){} /* prevent unused parameter warning */
	
	fval[0] = 0;
	
	lua_pushvalue(L, -1); /* push a copy of the function */
	
	lua_createtable(L, npts, 0);
	
	for(i = 0; i < npts; ++i){
		lua_pushinteger(L, i+1);
		lua_createtable(L, ndim, 0);
		for(j = 0; j < ndim; ++j){
			lua_pushinteger(L, j+1);
			lua_pushnumber(L, x[i*ndim+j]);
			lua_settable(L, -3);
		}
		lua_settable(L, -3);
	}
	if(0 != lua_pcall(L, 1, 1, 0)){
		S4L_error(L, "Error running Integrate function: %s", lua_tostring(L, -1));
        lua_pop(L, 1);
        return;
	}
	
	if(!lua_istable(L, -1)){
		lua_pop(L, 1);
		S4L_error(L, "Integrate function must return a list of results");
		return;
	}

	table_size = lua_rawlen(L, -1);
	
	if(npts != table_size){
		lua_pop(L, 1);
		S4L_error(L, "Integrate function returned a list of incorrect length");
		return;
	}
	
	for(i = 0; i < npts; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, -2);
		fval[i] = lua_tonumber(L, -1);
		lua_pop(L, 1);
	}
	lua_pop(L, 1);
}

static int S4L_Integrate(lua_State *L){
	/* syntax:
	 *  S4.Integrate(f, [0, 1], [0, 2], opts) will integrate f(x,y) on x[0,1], y[0,2]
	 */
	int i, nargs, ndim = 0;
	
	double *xmin;
	double *xmax;
	
	unsigned MaxEval = 1000000;
	double AbsoluteError = 0;
	double RelativeError = 1e-6;
	int Parallelize = 0;
	
	double val[1], err[1];
	
	nargs = lua_gettop(L);
	luaL_checktype(L, 1, LUA_TFUNCTION);
	
	xmin = (double*)malloc(sizeof(double) * 2*nargs);
	xmax = xmin+nargs;
	
	for(i = 2; i <= nargs; ++i){
		int table_size;
		if(!lua_istable(L, i)){
			S4L_error(L, "Expected table at argument %d to Integrate\n", i);
			return 0;
		}
		
		table_size = lua_rawlen(L, i);
		
		if(i < nargs){
			if(2 != table_size){
				S4L_error(L, "Argument %d to Integrate must be table of length 2", i);
				return 0;
			}else{
				lua_pushinteger(L, 1);
				lua_gettable(L, i);
				xmin[ndim] = lua_tonumber(L, -1);
				lua_pop(L, 1);
				
				lua_pushinteger(L, 2);
				lua_gettable(L, i);
				xmax[ndim] = lua_tonumber(L, -1);
				lua_pop(L, 1);
				
				++ndim;
			}
		}else if(i == nargs){
			/* table is in the stack at index 't' */
			lua_pushnil(L);  /* first key */
			while(lua_next(L, nargs) != 0){
				if(LUA_TSTRING == lua_type(L, -2)){
					const char *key = lua_tostring(L, -2);
					if(0 == strcmp(key, "MaxEval")){
						MaxEval = lua_tointeger(L, -1);
					}else if(0 == strcmp(key, "AbsoluteError")){
						AbsoluteError = lua_tonumber(L, -1);
					}else if(0 == strcmp(key, "RelativeError")){
						RelativeError = lua_tonumber(L, -1);
					}else if(0 == strcmp(key, "Parallelize")){
						Parallelize = lua_toboolean(L, -1);
					}else{
						S4L_error(L, "Unrecognized option to Integrate: %s", key);
						return 0;
					}
				}else{
					if(2 != table_size){
						S4L_error(L, "Argument %d to Integrate must be table of length 2 or options table", i);
						return 0;
					}
					lua_pushinteger(L, 1);
					lua_gettable(L, i);
					xmin[ndim] = lua_tonumber(L, -1);
					lua_pop(L, 1);
					
					lua_pushinteger(L, 2);
					lua_gettable(L, i);
					xmax[ndim] = lua_tonumber(L, -1);
					lua_pop(L, 1);
					
					++ndim;
					lua_pop(L, 1);
					break;
				}
				lua_pop(L, 1);
			}
		}
	}
	
	/*for(i = 0; i < ndim; ++i){ printf("%f\t%f\n", xmin[i], xmax[i]); } */
	
	/* pop off everything but the function */
	lua_pop(L, nargs-1);
	
	if(Parallelize){
		adapt_integrate_v(
			1, &IntegrateFunctionVectorized, (void*)L,
			ndim, xmin, xmax,
			MaxEval, AbsoluteError, RelativeError,
			val, err);
	}else{
		adapt_integrate(
			1, &IntegrateFunction, (void*)L,
			ndim, xmin, xmax,
			MaxEval, AbsoluteError, RelativeError,
			val, err);
	}
	lua_pop(L, 1); /* pop off the function */
	
	lua_pushnumber(L, val[0]);
	lua_pushnumber(L, err[0]);
	
	return 2;
}

static int S4L_NewInterpolator(lua_State *L){
	Interpolator *I;
	Interpolator_type type;
	int n, ny = 0, i, j, ld = 0;
	double *xy = NULL;
	const char *type_name;
	
	I = (Interpolator *)lua_newuserdata(L, sizeof(Interpolator));
	luaL_getmetatable(L, "S4.Interpolator");
	lua_setmetatable(L, -2);
	
	type_name = luaL_checklstring(L, 1, NULL);
	if(0 == strcmp(type_name, "cubic hermite spline")){
		type = Interpolator_CUBIC_HERMITE_SPLINE;
	}else if(1 || 0 == strcmp(type_name, "linear")){
		type = Interpolator_LINEAR;
	}
	
	if(!lua_istable(L, 2)){
		S4L_error(L, "NewInterpolator: Table expected for argument 2.");
		goto S4L_NewInterpolator_error;
	}
	
	n = lua_rawlen(L, 2);
	if(n < 2){
		S4L_error(L, "NewInterpolator: Table must be of length 2 or more.");
		goto S4L_NewInterpolator_error;
	}
	for(i = 0; i < n; ++i){
		double x;
		lua_pushinteger(L, i+1);
		lua_gettable(L, 2); /* {x, {y1, y2, ... }} */
		if(!lua_istable(L, -1)){
			S4L_error(L, "NewInterpolator: Table must contain tables.");
			goto S4L_NewInterpolator_error;
		}
		lua_pushinteger(L, 1);
		lua_gettable(L, -2); /* get x */
		x = lua_tonumber(L, -1); lua_pop(L, 1);
		
		lua_pushinteger(L, 2);
		lua_gettable(L, -2); /* {y1, y2, ... } */
		if(!lua_istable(L, -1)){
			S4L_error(L, "NewInterpolator: Table must contain tables of form {x, {y1, y2, ...}.");
			goto S4L_NewInterpolator_error;
		}
		if(0 == i){
			ny = lua_rawlen(L, -1);
			ld = 1+ny;
			xy = (double*)malloc(sizeof(double)*n*ld);
		}
		xy[i*ld] = x;
		for(j = 0; j < ny; ++j){
			lua_pushinteger(L, j+1);
			lua_gettable(L, -2);
			xy[1+j+i*ld] = lua_tonumber(L, -1);
			lua_pop(L, 1);
		}
		lua_pop(L, 2);
	}
	
	*I = Interpolator_New(n, ny, xy, type);
	free(xy);
	return 1;
S4L_NewInterpolator_error:
	*I = NULL;
	return 1;
}
static int S4L_Interpolator__gc(lua_State *L){
	Interpolator *I = (Interpolator *)luaL_checkudata(L, 1, "S4.Interpolator");
	Interpolator_Destroy(*I);
	return 0;
}
static int S4L_Interpolator_Get(lua_State *L){
	double *y;
	int ny, i;
	Interpolator *I = (Interpolator *)luaL_checkudata(L, 1, "S4.Interpolator");
	luaL_argcheck(L, I != NULL, 1, "Get: 'Interpolator' object expected.");
	
	y = Interpolator_Get(*I, luaL_checknumber(L, 2), &ny);
	for(i = 0; i < ny; ++i){
		lua_pushnumber(L, y[i]);
	}
	return ny;
}
static int S4L_NewSpectrumSampler(lua_State *L){
	SpectrumSampler *sampler;
	SpectrumSampler_Options options;
	double x0, x1;
	
	options.initial_num_points = 33;
	options.range_threshold = 0.001;
	options.max_bend = cos(10.0*M_PI/180.0);
	options.min_dx = 1e-6;
	options.parallelize = 0;
	
	x0 = luaL_checknumber(L, 1);
	x1 = luaL_checknumber(L, 2);
	if(lua_gettop(L) > 2){
		lua_pushstring(L, "InitialNumPoints");
		lua_gettable(L, 3);
		if(!lua_isnil(L, -1)){
			options.initial_num_points = lua_tointeger(L, -1);
		}
		lua_pop(L, 1);
		
		lua_pushstring(L, "RangeThreshold");
		lua_gettable(L, 3);
		if(!lua_isnil(L, -1)){
			options.range_threshold = lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
		
		lua_pushstring(L, "MaxBend");
		lua_gettable(L, 3);
		if(!lua_isnil(L, -1)){
			options.max_bend = lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
		
		lua_pushstring(L, "MinimumSpacing");
		lua_gettable(L, 3);
		if(!lua_isnil(L, -1)){
			options.min_dx = lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
		
		lua_pushstring(L, "Parallelize");
		lua_gettable(L, 3);
		if(!lua_isnil(L, -1)){
			options.parallelize = lua_toboolean(L, -1);
		}
		lua_pop(L, 1);
	}
	
	sampler = (SpectrumSampler *)lua_newuserdata(L, sizeof(SpectrumSampler));
	luaL_getmetatable(L, "S4.SpectrumSampler");
	lua_setmetatable(L, -2);
	
	*sampler = SpectrumSampler_New(x0, x1, &options);
	return 1;
}
static int S4L_SpectrumSampler__gc(lua_State *L){
	SpectrumSampler *sampler = (SpectrumSampler *)luaL_checkudata(L, 1, "S4.SpectrumSampler");
	SpectrumSampler_Destroy(*sampler);
	return 0;
}
static int S4L_SpectrumSampler_IsDone(lua_State *L){
	SpectrumSampler *sampler = (SpectrumSampler *)luaL_checkudata(L, 1, "S4.SpectrumSampler");
	luaL_argcheck(L, sampler != NULL, 1, "IsDone: 'SpectrumSampler' object expected.");
	
	lua_pushboolean(L, SpectrumSampler_IsDone(*sampler));
	return 1;
}
static int S4L_SpectrumSampler_GetFrequency(lua_State *L){
	SpectrumSampler *sampler = (SpectrumSampler *)luaL_checkudata(L, 1, "S4.SpectrumSampler");
	luaL_argcheck(L, sampler != NULL, 1, "GetFrequency: 'SpectrumSampler' object expected.");
	luaL_argcheck(L, !SpectrumSampler_IsParallelized(*sampler), 1, "GetFrequency can only be used when Parallelize = false");
	
	lua_pushnumber(L, SpectrumSampler_GetFrequency(*sampler));
	return 1;
}
static int S4L_SpectrumSampler_SubmitResult(lua_State *L){
	double y;
	SpectrumSampler *sampler = (SpectrumSampler *)luaL_checkudata(L, 1, "S4.SpectrumSampler");
	luaL_argcheck(L, sampler != NULL, 1, "SubmitResult: 'SpectrumSampler' object expected.");
	luaL_argcheck(L, !SpectrumSampler_IsParallelized(*sampler), 1, "SubmitResult can only be used when Parallelize = false");
	
	y = luaL_checknumber(L, 2);
	lua_pushboolean(L, SpectrumSampler_SubmitResult(*sampler, y));
	return 1;
}
static int S4L_SpectrumSampler_GetFrequencies(lua_State *L){
	int i, nf;
	const double *f;
	SpectrumSampler *sampler = (SpectrumSampler *)luaL_checkudata(L, 1, "S4.SpectrumSampler");
	luaL_argcheck(L, sampler != NULL, 1, "GetFrequencies: 'SpectrumSampler' object expected.");
	luaL_argcheck(L, SpectrumSampler_IsParallelized(*sampler), 1, "GetFrequencies can only be used when Parallelize = true");
	
	nf = SpectrumSampler_GetFrequencies(*sampler, &f);
	if(nf <= 0){
		return 0;
	}
	
	lua_createtable(L, nf, 0);
	for(i = 0; i < nf; ++i){
		lua_pushinteger(L, i+1);
		lua_pushnumber(L, f[i]);
		lua_settable(L, -3);
	}
	
	return 1;
}
static int S4L_SpectrumSampler_SubmitResults(lua_State *L){
	int i, ny;
	double *y;
	SpectrumSampler *sampler = (SpectrumSampler *)luaL_checkudata(L, 1, "S4.SpectrumSampler");
	luaL_argcheck(L, sampler != NULL, 1, "SubmitResults: 'SpectrumSampler' object expected.");
	luaL_argcheck(L, SpectrumSampler_IsParallelized(*sampler), 1, "SubmitResults can only be used when Parallelize = true");
	
	ny = SpectrumSampler_GetSubmissionBuffer(*sampler, &y);
	luaL_argcheck(L, ny > 0, 1, "Internal buffer error");
	
	luaL_argcheck(L, lua_istable(L, 2), 2, "Expected table of values");
	luaL_argcheck(L, (int)lua_rawlen(L, 2) == ny, 2, "Length of table of values did not match expected length");
	
	for(i = 0; i < ny; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, -2);
		y[i] = luaL_checknumber(L, -1);
		lua_pop(L, 1);
	}
	
	lua_pushboolean(L, SpectrumSampler_SubmitResults(*sampler));
	return 1;
}
static int S4L_SpectrumSampler_GetSpectrum(lua_State *L){
	int i, n;
	SpectrumSampler_Enumerator e;
	double pt[2];
	SpectrumSampler *sampler = (SpectrumSampler *)luaL_checkudata(L, 1, "S4.SpectrumSampler");
	luaL_argcheck(L, sampler != NULL, 1, "SubmitResult: 'SpectrumSampler' object expected.");
	
	n = SpectrumSampler_GetNumPoints(*sampler);
	lua_createtable(L, n, 0);
	
	e = SpectrumSampler_GetPointEnumerator(*sampler);
	i = 1;
	while(SpectrumSampler_Enumerator_Get(e, pt)){
		lua_pushinteger(L, i++);
		lua_createtable(L, 2, 0);
		lua_pushinteger(L, 1);
		lua_pushnumber(L, pt[0]);
		lua_settable(L, -3);
		lua_pushinteger(L, 2);
		lua_pushnumber(L, pt[1]);
		lua_settable(L, -3);
		lua_settable(L, -3);
	}
	return 1;
}
static int S4L_NewSimulation(lua_State *L){
	Simulation *S;
	lua_gc(L, LUA_GCCOLLECT, 0);
	S = (Simulation *)lua_newuserdata(L, sizeof(Simulation));
	luaL_getmetatable(L, "S4.Simulation");
	lua_setmetatable(L, -2);
	Simulation_Init(S);
	return 1;
}
static int S4L_Simulation__gc(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	Simulation_Destroy(S);
	return 0;
}
static int S4L_Simulation_Clone(lua_State *L){
	Simulation *T;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "Clone: 'Simulation' object expected.");
	
	T = (Simulation *)lua_newuserdata(L, sizeof(Simulation));
	luaL_getmetatable(L, "S4.Simulation");
	lua_setmetatable(L, -2);
	
	Simulation_Clone(S, T);
	return 1;
}

static int S4L_Simulation_SetLattice(lua_State *L){
	int i,j;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "SetLattice: 'Simulation' object expected.");
	
	if(lua_isnumber(L, 2)){
		double period = lua_tonumber(L, 2);
		luaL_argcheck(L, period > 0, 2, "1D lattice constant must be positive");
		S->Lr[0] = period;
		S->Lr[1] = 0;
		S->Lr[2] = 0;
		S->Lr[3] = 0;
		Simulation_MakeReciprocalLattice(S);
	}else{
		luaL_checktype(L, 2, LUA_TTABLE);
		luaL_checktype(L, 3, LUA_TTABLE);
		
		for(i = 0; i < 2; ++i){
			for(j = 0; j < 2; ++j){
				lua_pushinteger(L, 1+j);
				lua_gettable(L, 2+i);
				if(!lua_isnumber(L, -1)){
					S4L_error(L, "SetLattice: Lattice coordinates must be numeric.");
				}else{
					S->Lr[2*i+j] = (double)lua_tonumber(L, -1);
				}
				lua_pop(L, 1);
			}
		}
		i = Simulation_MakeReciprocalLattice(S);
		if(0 != i){
			switch(i){
			case 1: /* degenerate */
				S4L_error(L, "SetLattice: Lattice vectors are degenerate (for 1D lattice, set second vector to zero).");
				break;
			case 2: /* both zero */
				S4L_error(L, "SetLattice: Lattice vectors are both zero.");
				break;
			default:
				break;
			}
		}
	}
	return 0;
}
static int S4L_Simulation_GetReciprocalLattice(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "GetReciprocalLattice: 'Simulation' object expected.");
	
	lua_createtable(L, 2, 0);    /* {}                   */
	lua_createtable(L, 2, 0);    /* {} {}                */
	lua_pushnumber(L, S->Lk[0]); /* {} {} L0             */
	lua_rawseti(L, -2, 1);       /* {} {L0}              */
	lua_pushnumber(L, S->Lk[1]); /* {} {L0} L1           */
	lua_rawseti(L, -2, 2);       /* {} {L0, L1}          */
	lua_rawseti(L, -2, 1);       /* {{L0, L1}}           */
	lua_createtable(L, 2, 0);    /* {{L0, L1}} {}        */
	lua_pushnumber(L, S->Lk[2]); /* {{L0, L1}} {} L2     */
	lua_rawseti(L, -2, 1);       /* {{L0, L1}} {L2}      */
	lua_pushnumber(L, S->Lk[3]); /* {{L0, L1}} {L2} L3   */
	lua_rawseti(L, -2, 2);       /* {{L0, L1}} {L2, L3}  */
	lua_rawseti(L, -2, 2);       /* {{L0, L1}, {L2, L3}} */
	
	return 1;
}

static int S4L_Simulation_SetNumG(lua_State *L){
	int n;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "SetNumG: 'Simulation' object expected.");
	n = luaL_checkint(L, 2);
	if(n < 1){
		S4L_error(L, "SetNumG: Must have at least 1 G-vector.");
	}
	Simulation_SetNumG(S, n);
	return 0;
}
static int S4L_Simulation_GetNumG(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "GetNumG: 'Simulation' object expected.");
	lua_pushinteger(L, Simulation_GetNumG(S, NULL));
	return 1;
}
static int S4L_Simulation_SetResolution(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "SetResolution: 'Simulation' object expected.");
	S->options.resolution = luaL_checkint(L, 2);
	if(S->options.resolution < 2){
		S4L_error(L, "SetResolution: Resolution must be at least 2.");
	}
	return 0;
}
static int S4L_Simulation_AddMaterial(lua_State *L){
	int i, j;
	double eps[18];
	Material *M;
	const char *name;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "AddMaterial: 'Simulation' object expected.");
	
	M = Simulation_AddMaterial(S);
	name = luaL_checklstring(L, 2, NULL);
	if(NULL == M){
		S4L_error(L, "AddMaterial: There was a problem allocating the material named '%s'.", name);
		return 0;
	}
	
	if(2 == lua_rawlen(L, 3)){
		for(j = 0; j < 2; ++j){
			lua_pushinteger(L, 1+j);
			lua_gettable(L, 3);
			if(!lua_isnumber(L, -1)){
				S4L_error(L, "AddMaterial: Material epsilon must be numeric.");
			}else{
				eps[j] = (double)lua_tonumber(L, -1);
			}
			lua_pop(L, 1);
		}
		Material_Init(M,
			name,
			eps);
	}else if(9 == lua_rawlen(L, 3)){
		for(i = 0; i < 9; ++i){
			lua_pushinteger(L, 1+i);
			lua_gettable(L, 3);
			for(j = 0; j < 2; ++j){
				lua_pushinteger(L, 1+j);
				lua_gettable(L, 4);
				if(!lua_isnumber(L, -1)){
					S4L_error(L, "AddMaterial: Material epsilon must be numeric.");
				}else{
					eps[2*i+j] = (double)lua_tonumber(L, -1);
				}
				lua_pop(L, 1);
			}
			lua_pop(L, 1);
		}
		/* [ a b c ]    [ a b   ]
		 * [ d e f ] -> [ d e   ]
		 * [ g h i ]    [     i ]
		 */
		eps[4] = eps[6]; eps[5] = eps[7];
		eps[6] = eps[8]; eps[7] = eps[9];
		eps[8] = eps[16]; eps[9] = eps[17];
		Material_InitTensor(M,
			name,
			eps);
	}else{
		S4L_error(L, "AddMaterial: Expected either a scalar or tensor value for material %s.", name);
	}
	
	return 0;
}
static int S4L_Simulation_SetMaterial(lua_State *L){
	int i, j;
	double eps[18];
	const char *name;
	Material *M;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "SetMaterial: 'Simulation' object expected.");
	
	name = luaL_checklstring(L, 2, NULL);
	M = Simulation_GetMaterialByName(S, name, NULL);
	if(NULL == M){
		M = Simulation_AddMaterial(S);
		if(NULL == M){
			S4L_error(L, "SetMaterial: There was a problem allocating the material named '%s'.", name);
			return 0;
		}
	}
	
	if(2 == lua_rawlen(L, 3)){
		for(j = 0; j < 2; ++j){
			lua_pushinteger(L, 1+j);
			lua_gettable(L, 3);
			if(!lua_isnumber(L, -1)){
				S4L_error(L, "SetMaterial: Material epsilon must be numeric.");
			}else{
				eps[j] = (double)lua_tonumber(L, -1);
			}
			lua_pop(L, 1);
		}
		M->eps.s[0] = eps[0];
		M->eps.s[1] = eps[1];
	}else if(9 == lua_rawlen(L, 3)){
		for(i = 0; i < 9; ++i){
			lua_pushinteger(L, 1+i);
			lua_gettable(L, 3);
			for(j = 0; j < 2; ++j){
				lua_pushinteger(L, 1+j);
				lua_gettable(L, 4);
				if(!lua_isnumber(L, -1)){
					S4L_error(L, "SetMaterial: Material epsilon must be numeric.");
				}else{
					eps[2*i+j] = (double)lua_tonumber(L, -1);
				}
				lua_pop(L, 1);
			}
			lua_pop(L, 1);
		}
		/* [ a b c ]    [ a b   ]
		 * [ d e f ] -> [ d e   ]
		 * [ g h i ]    [     i ]
		 */
		eps[4] = eps[6]; eps[5] = eps[7];
		eps[6] = eps[8]; eps[7] = eps[9];
		eps[8] = eps[16]; eps[9] = eps[17];
		for(i = 0; i < 10; ++i){
			M->eps.abcde[i] = eps[i];
		}
	}else{
		S4L_error(L, "SetMaterial: Expected either a scalar or tensor value for material %s.", name);
	}
	
	return 0;
}
static int S4L_Simulation_AddLayer(lua_State *L){
	Layer *layer;
	const char *name;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "AddLayer: 'Simulation' object expected.");
	
	layer = Simulation_AddLayer(S);
	name = luaL_checklstring(L, 2, NULL);
	if(NULL == layer){
		S4L_error(L, "AddLayer: There was a problem allocating the layer named '%s'.", name);
		return 0;
	}
	Layer_Init(layer,
		name,
		luaL_checknumber(L, 3),
		luaL_checklstring(L, 4, NULL),
		NULL);

	return 0;
}
static int S4L_Simulation_SetLayer(lua_State *L){
	Layer *layer;
	const char *name;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "SetLayer: 'Simulation' object expected.");
	
	name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, name, NULL);
	if(NULL == layer){
		layer = Simulation_AddLayer(S);
		if(NULL == layer){
			S4L_error(L, "SetLayer: There was a problem allocating the layer named '%s'.", name);
			return 0;
		}
		Layer_Init(layer,
			name,
			luaL_checknumber(L, 3),
			luaL_checklstring(L, 4, NULL),
			NULL);
	}else{
		Simulation_RemoveLayerPatterns(S, layer);
		layer->thickness = luaL_checknumber(L, 3);
	}
	return 0;
}
static int S4L_Simulation_SetLayerThickness(lua_State *L){
	Layer *layer;
	const char *name;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "SetLayerThickness: 'Simulation' object expected.");
	
	name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, name, NULL);
	if(NULL == layer){
		S4L_error(L, "SetLayerThickness: Layer named '%s' not found.", layer);
	}else{
		double thick = luaL_checknumber(L, 3);
		if(thick < 0){
			S4L_error(L, "SetLayerThickness: Thickness must be non-negative.");
		}
		Simulation_ChangeLayerThickness(S, layer, &thick);
	}
	return 0;
}

static int S4L_Simulation_AddLayerCopy(lua_State *L){
	Layer *layer;
	const char *name;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "AddLayerCopy: 'Simulation' object expected.");
	
	layer = Simulation_AddLayer(S);
	name = luaL_checklstring(L, 2, NULL);
	if(NULL == layer){
		S4L_error(L, "AddLayerCopy: There was a problem allocating the layer named '%s'.", name);
		return 0;
	}
	Layer_Init(layer,
		name,
		luaL_checknumber(L, 3),
		NULL,
		luaL_checklstring(L, 4, NULL));

	return 0;
}

/* Expected stack:
 *   1 Simulation
 *   2 layer name string
 *   3 circle material string
 *   4 circle center (numeric table of length 2)
 *   5 circle radius
 */
static int S4L_Simulation_SetLayerPatternCircle(lua_State *L){
	int j, ret;
	const char *layer_name;
	const char *material_name;
	int material_index;
	Layer *layer;
	Material *M;
	double center[2];
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "SetLayerPatternCircle: 'Simulation' object expected.");
	
	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "SetLayerPatternCircle: Layer named '%s' not found.", layer_name);
		return 0;
	}
	if(NULL != layer->copy){
		S4L_error(L, "SetLayerPatternCircle: Cannot pattern a layer copy.");
		return 0;
	}
	material_name = luaL_checklstring(L, 3, NULL);
	M = Simulation_GetMaterialByName(S, material_name, &material_index);
	if(NULL == M){
		S4L_error(L, "SetLayerPatternCircle: Material named '%s' not found.", material_name);
		return 0;
	}
	for(j = 0; j < 2; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 4);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "SetLayerPatternCircle: Circle center must be numeric.");
			return 0;
		}else{
			center[j] = (double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	ret = Simulation_AddLayerPatternCircle(S, layer, material_index, center, luaL_checknumber(L, 5));
	if(0 != ret){
		S4L_error(L, "SetLayerPatternCircle: There was a problem allocating the pattern.");
		return 0;
	}
	return 0;
}

/* Expected stack:
 *   1 Simulation
 *   2 layer name string
 *   3 ellipse material string
 *   4 ellipse center (numeric table of length 2)
 *   5 ellipse tilt angle (degrees)
 *   6 ellipse halfwidths (numeric table of length 2)
 */
static int S4L_Simulation_SetLayerPatternEllipse(lua_State *L){
	int j, ret;
	const char *layer_name;
	const char *material_name;
	int material_index;
	Layer *layer;
	Material *M;
	double center[2], halfwidths[2];
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "SetLayerPatternEllipse: 'Simulation' object expected");
	
	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "SetLayerPatternEllipse: Layer named '%s' not found.", layer_name);
		return 0;
	}
	if(NULL != layer->copy){
		S4L_error(L, "SetLayerPatternEllipse: Cannot pattern a layer copy.");
		return 0;
	}
	material_name = luaL_checklstring(L, 3, NULL);
	M = Simulation_GetMaterialByName(S, material_name, &material_index);
	if(NULL == M){
		S4L_error(L, "SetLayerPatternEllipse: Material named '%s' not found.", material_name);
		return 0;
	}
	for(j = 0; j < 2; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 4);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "SetLayerPatternEllipse: Ellipse center must be numeric.");
			return 0;
		}else{
			center[j] = (double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	for(j = 0; j < 2; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 6);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "SetLayerPatternEllipse: Halfwidths must be numeric.");
			return 0;
		}else{
			halfwidths[j] = (double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	ret = Simulation_AddLayerPatternEllipse(S, layer, material_index,
		center,
		(M_PI/180.)*luaL_checknumber(L, 5),
		halfwidths);
	if(0 != ret){
		S4L_error(L, "SetLayerPatternEllipse: There was a problem allocating the pattern.");
		return 0;
	}
	return 0;
}

/* Expected stack:
 *   1 Simulation
 *   2 layer name string
 *   3 rectangle material string
 *   4 rectangle center (numeric table of length 2)
 *   5 rectangle tilt angle (degrees)
 *   6 rectangle halfwidths (numeric table of length 2)
 */
static int S4L_Simulation_SetLayerPatternRectangle(lua_State *L){
	int j, ret;
	const char *layer_name;
	const char *material_name;
	int material_index;
	Layer *layer;
	Material *M;
	double center[2], halfwidths[2];
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "SetLayerPatternRectangle: 'Simulation' object expected.");
	
	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "SetLayerPatternRectangle: Layer named '%s' not found.", layer_name);
		return 0;
	}
	if(NULL != layer->copy){
		S4L_error(L, "SetLayerPatternRectangle: Cannot pattern a layer copy.");
		return 0;
	}
	material_name = luaL_checklstring(L, 3, NULL);
	
	M = Simulation_GetMaterialByName(S, material_name, &material_index);
	if(NULL == M){
		S4L_error(L, "SetLayerPatternRectangle: Material named '%s' not found.", material_name);
		return 0;
	}
	for(j = 0; j < 2; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 4);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "SetLayerPatternRectangle: Rectangle center must be numeric.");
			return 0;
		}else{
			center[j] = (double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	for(j = 0; j < 2; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 6);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "SetLayerPatternRectangle: Halfwidths must be numeric.");
			return 0;
		}else{
			halfwidths[j] = (double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	ret = Simulation_AddLayerPatternRectangle(S,
		layer, material_index,
		center,
		(M_PI/180.)*luaL_checknumber(L, 5),
		halfwidths);
	if(0 != ret){
		S4L_error(L, "SetLayerPatternRectangle: There was a problem allocating the pattern.");
		return 0;
	}
	return 0;
}

/* Expected stack:
 *   1 Simulation
 *   2 layer name string
 *   3 polygon material string
 *   4 polygon center (numeric table of length 2)
 *   5 polygon tilt angle (degrees)
 *   6 polygon vertices (numeric table of length 2*nverts)
 */
static int S4L_Simulation_SetLayerPatternPolygon(lua_State *L){
	int i,j, ret, nvert;
	double center[2];
	double *vert;
	const char *layer_name;
	const char *material_name;
	int material_index;
	Layer *layer;
	Material *M;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "SetLayerPatternPolygon: 'Simulation' object expected.");
	
	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "SetLayerPatternPolygon: Layer named '%s' not found.", layer_name);
		return 0;
	}
	if(NULL != layer->copy){
		S4L_error(L, "SetLayerPatternPolygon: Cannot pattern a layer copy.");
		return 0;
	}
	material_name = luaL_checklstring(L, 3, NULL);
	M = Simulation_GetMaterialByName(S, material_name, &material_index);
	if(NULL == M){
		S4L_error(L, "SetLayerPatternPolygon: Material named '%s' not found.", material_name);
		return 0;
	}
	for(j = 0; j < 2; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 4);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "SetLayerPatternPolygon: Polygon center must be numeric.");
			return 0;
		}else{
			center[j] = (double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	nvert = lua_rawlen(L, 6) / 2;
	vert = (double*)malloc(sizeof(double)*2*nvert);
	for(j = 0; j < nvert; ++j){
		for(i = 0; i < 2; ++i){
			lua_pushinteger(L, 1+2*j+i);
			lua_gettable(L, 6);
			if(!lua_isnumber(L, -1)){
				S4L_error(L, "SetLayerPatternPolygon: Polygon vertexes must be numeric.");
				return 0;
			}else{
				vert[2*j+i] = (double)lua_tonumber(L, -1);
			}
			lua_pop(L, 1);
		}
	}
	ret = Simulation_AddLayerPatternPolygon(S,
		layer, material_index,
		center,
		(M_PI/180.)*luaL_checknumber(L, 5),
		nvert,
		vert);
	free(vert);
	if(0 != ret){
		S4L_error(L, "SetLayerPatternPolygon: There was a problem allocating the pattern.");
		return 0;
	}
	return 0;
}


/* Expected stack:
 *   1 Simulation
 *   2 table
 * First table:
 *   Each entry a table of three entries, first is integer, second is 'x' or 'y', third is table of two entries re and im
 *   If integer is negative, indicates a backward propagating mode.
 */
static int S4L_Simulation_SetExcitationExterior(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	double *ex;
	int *exg;
	const char *pol;
	int i, j, n;
	luaL_argcheck(L, lua_istable(L, 2), 2, "SetExcitationExterior: table of order coefficients expected.");
	
	n = lua_rawlen(L, 2);
	ex = (double*)malloc(sizeof(double) * 2*n);
	exg = (int*)malloc(sizeof(int) * 2*n);
	
	for(i = 0; i < n; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, 2);
		if(!lua_istable(L, -1)){
			S4L_error(L, "SetExcitationExterior: List items must be {G-index, pol, {re, im}}.");
		}
		/* Get G index */
		lua_rawgeti(L, -1, 1);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "SetExcitationExterior: List items must be {G-index, pol, {re, im}}.");
		}
		exg[2*i+0] = lua_tointeger(L, -1);
		if(0 == exg[2*i+0]){
			S4L_error(L, "SetExcitationExterior: G-index cannot be zero.");
		}
		lua_pop(L, 1);
		
		/* Get polarization */
		lua_rawgeti(L, -1, 2);
		if(!lua_isstring(L, -1)){
			S4L_error(L, "SetExcitationExterior: List items must be {G-index, pol, {re, im}}.");
		}
		pol = lua_tostring(L, -1);
		if(0 == strcmp("x", pol)){
			exg[2*i+1] = 0;
		}else if(0 == strcmp("y", pol)){
			exg[2*i+1] = 1;
		}else{
			S4L_error(L, "SetExcitationExterior: polarization must be 'x' or 'y'.");
		}
		lua_pop(L, 1);
		
		/* get {re,im} */
		lua_rawgeti(L, -1, 3);
		if(!lua_istable(L, -1)){
			S4L_error(L, "SetExcitationExterior: List items must be {G-index, pol, {re, im}}.");
		}
		for(j = 0; j < 2; ++j){
			lua_rawgeti(L, -1, j+1);
			if(!lua_isnumber(L, -1)){
				S4L_error(L, "SetExcitationExterior: List items must be {G-index, pol, {re, im}}.");
			}
			ex[2*i+j] = lua_tonumber(L, -1);
			lua_pop(L, 1);
		}
		lua_pop(L, 1);
		
		lua_pop(L, 1);
	}
	Simulation_MakeExcitationExterior(S, n, exg, ex);
	return 0;
}

static int S4L_Simulation_SetExcitationInterior(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	S4L_error(L, "SetExcitationInterior: Not implemented.");
	return 0;
}

/* Expected stack:
 *   1 Simulation
 *   2 layer name string
 *   3 k (numeric table of length 2)
 *   4 position in plane (numeric table of length 2)
 *   5 dipole moment (numeric table of length 3)
 *   6 amplitude and phase (numeric table of length 2; phase in degrees)
 */
static int S4L_Simulation_SetExcitationDipole(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	const char *layer;
	int j, ret;
	double k[2], pos[2], moment[6], ampphase[2];
	luaL_argcheck(L, S != NULL, 1, "SetExcitationDipole: 'Simulation' object expected.");
	layer = luaL_checkstring(L, 2);
	luaL_argcheck(L, lua_istable(L, 3), 3, "SetExcitationDipole: k-vector expected.");
	luaL_argcheck(L, lua_istable(L, 4), 4, "SetExcitationDipole: position expected.");
	luaL_argcheck(L, lua_istable(L, 5), 5, "SetExcitationDipole: dipole moment expected.");
	luaL_argcheck(L, lua_isnumber(L, 6) || lua_istable(L, 6), 6, "SetExcitationDipole: amplitude/phase expected.");
	
	for(j = 0; j < 2; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 3);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "SetExcitationDipole: k-vector must be numeric.");
			return 0;
		}else{
			k[j] = (double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	for(j = 0; j < 2; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 4);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "SetExcitationDipole: position must be numeric.");
			return 0;
		}else{
			pos[j] = (double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	for(j = 0; j < 3; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 5);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "SetExcitationDipole: dipole moment must be numeric.");
			return 0;
		}else{
			moment[j] = (double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	if(lua_isnumber(L, 6)){
		ampphase[0] = (double)lua_tonumber(L, 6);
		ampphase[1] = 0;
	}else{
		for(j = 0; j < 2; ++j){
			lua_pushinteger(L, 1+j);
			lua_gettable(L, 6);
			if(!lua_isnumber(L, -1)){
				S4L_error(L, "SetExcitationDipole: amplitude/phase must be numeric.");
				return 0;
			}else{
				ampphase[j] = (double)lua_tonumber(L, -1);
			}
			lua_pop(L, 1);
		}
	}

	{
		double imomentlen = 1. / hypot(hypot(moment[0], moment[1]), moment[2]);
		moment[4] = (moment[2] * imomentlen) * ampphase[0];
		moment[5] = (moment[2] * imomentlen) * ampphase[1];
		moment[2] = (moment[1] * imomentlen) * ampphase[0];
		moment[3] = (moment[1] * imomentlen) * ampphase[1];
		moment[0] = (moment[0] * imomentlen) * ampphase[0];
		moment[1] = (moment[0] * imomentlen) * ampphase[1];
	}
	ret = Simulation_MakeExcitationDipole(S, k, layer, pos, moment);
	if(0 != ret){
		HandleSolutionErrorCode(L, "SetExcitationDipole", ret);
		return 0;
	}
	
	return 0;
}

/* Expected stack:
 *   1 Simulation
 *   2 incidence angles (numeric table of length 2; phi and theta in degrees)
 *   3 TE amplitude and phase (numeric table of length 2; phase in degrees)
 *   4 TM amplitude and phase (numeric table of length 2; phase in degrees)
 */
static int S4L_Simulation_SetExcitationPlanewave(lua_State *L){
	int j, ret;
	int order = 0;
	double angle[2];
	double pol_s[2]; /* s polarization; E out of plane */
	double pol_p[2]; /* p polarization; E in plane */
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "SetExcitationPlanewave: 'Simulation' object expected.");

	Simulation_DestroySolution(S);

	for(j = 0; j < 2; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 2);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "SetExcitationPlanewave: Incidence angle must be numeric.");
			return 0;
		}else{
			angle[j] = (M_PI/180.)*(double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	for(j = 0; j < 2; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 3);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "SetExcitationPlanewave: s-polarization amplitude and phase must be numeric.");
			return 0;
		}else{
			pol_s[j] = (double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	for(j = 0; j < 2; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 4);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "SetExcitationPlanewave: p-polarization amplitude and phase must be numeric.");
			return 0;
		}else{
			pol_p[j] = (double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	pol_s[1] *= (M_PI/180.);
	pol_p[1] *= (M_PI/180.);

	order = luaL_optunsigned(L, 5, 0);
	if(order > 0){ order--; }
	ret = Simulation_MakeExcitationPlanewave(S, angle, pol_s, pol_p, order);
	if(0 != ret){
		HandleSolutionErrorCode(L, "SetExcitationPlanewave", ret);
		return 0;
	}
	
	return 0;
}
static int S4L_Simulation_SetFrequency(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "SetFrequency: 'Simulation' object expected");

	Simulation_DestroySolution(S);

	S->omega[0] = 2*M_PI*luaL_checknumber(L, 2);
	S->omega[1] = 0;
	if(lua_gettop(L) > 2){
		S->omega[1] = 2*M_PI*lua_tonumber(L, 3);
	}
	if(S->omega[0] <= 0){
		S4L_error(L, "SetFrequency: Frequency must be positive.");
		return 0;
	}
	if(S->omega[1] > 0){
		S4L_error(L, "SetFrequency: Imaginary component of frequency must be negative.");
		return 0;
	}
	return 0;
}


static int S4L_Simulation_GetGList(lua_State *L){
	int *G;
	int n, i, j, ret;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "GetGList: 'Simulation' object expected.");
	
	ret = Simulation_InitSolution(S);
	if(0 != ret){
		HandleSolutionErrorCode(L, "GetGList", ret);
		return 0;
	}
	
	n = Simulation_GetNumG(S, &G);
	if(NULL == G){
		return 0;
	}
	
	lua_createtable(L, n, 0);
	for(i = 0; i < n; ++i){
		lua_pushinteger(L, i+1);
		lua_createtable(L, 2, 0);
		for(j = 0; j < 2; ++j){
			lua_pushinteger(L, j+1);
			lua_pushinteger(L, G[2*i+j]);
			lua_settable(L, -3);
		}
		lua_settable(L, -3);
	}
	
	return 1;
}


static int S4L_Simulation_GetDiffractionOrder(lua_State *L){
	int *G;
	int u, v = 0, n, i;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "GetDiffractionOrder: 'Simulation' object expected.");
	
	u = luaL_checkint(L, 2);
	if(lua_gettop(L) > 2){
		v = luaL_checkint(L, 3);
	}
	
	i = Simulation_InitSolution(S);
	if(0 != i){
		HandleSolutionErrorCode(L, "GetDiffractionOrder", i);
		return 0;
	}
	
	n = Simulation_GetNumG(S, &G);
	if(NULL == G){
		return 0;
	}
	
	for(i = 0; i < n; ++i){
		if(u == G[2*i+0] && v == G[2*i+1]){
			lua_pushinteger(L, i+1);
			return 1;
		}
	}
	
	return 0;
}


/* Expected stack:
 *   1 Simulation
 *   2 layer name string
 *   3 layer z-offset
 */
static int S4L_Simulation_GetPoyntingFlux(lua_State *L){
	double power[4];
	int ret;
	const char *layer_name;
	Layer *layer;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "GetPoyntingFlux: 'Simulation' object expected.");

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "GetPoyntingFlux: Layer named '%s' not found.", layer_name);
		return 0;
	}
	ret = Simulation_GetPoyntingFlux(S, 
		layer,
		luaL_checknumber(L, 3),
		power);

	lua_pushnumber(L, power[0]); /* real forw (time averaged) */
	lua_pushnumber(L, power[1]); /* real back (time averaged) */
	lua_pushnumber(L, power[2]); /* imag forw */
	lua_pushnumber(L, power[3]); /* imag back */
	
	if(0 != ret){
		HandleSolutionErrorCode(L, "GetPoyntingFlux", ret);
	}
	return 4;
}
static int S4L_Simulation_GetPoyntingFluxByOrder(lua_State *L){
	double *power;
	int *G;
	int n, i, j, ret;
	const char *layer_name;
	Layer *layer;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "GetPoyntingFluxByOrder: 'Simulation' object expected.");

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "GetPoyntingFluxByOrder: Layer named '%s' not found.", layer_name);
		return 0;
	}
	
	ret = Simulation_SolveLayer(S, layer);
	if(0 != ret){
		HandleSolutionErrorCode(L, "GetPoyntingFluxByOrder", ret);
		return 0;
	}
	
	n = Simulation_GetNumG(S, &G);
	if(NULL == G){
		return 0;
	}
	
	power = (double*)malloc(sizeof(double)*4*n);
	Simulation_GetPoyntingFluxByG(S, 
		layer,
		luaL_checknumber(L, 3),
		power);
	
	lua_createtable(L, n, 0);
	for(i = 0; i < n; ++i){
		lua_pushinteger(L, i+1);
		/* push the diffracted powers */
		lua_createtable(L, 4, 0);
		for(j = 0; j < 4; ++j){
			lua_pushinteger(L, j+1);
			lua_pushnumber(L, power[4*i+j]);
			lua_settable(L, -3);
		}
		lua_settable(L, -3);
	}
	free(power);
	return 1;
}

static int S4L_Simulation_GetAmplitudes(lua_State *L){
	double *amp;
	int *G;
	int n, i, k, ret;
	const char *layer_name;
	Layer *layer;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "GetAmplitudes: 'Simulation' object expected.");

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "GetAmplitudes: Layer named '%s' not found.", layer_name);
		return 0;
	}
	
	ret = Simulation_SolveLayer(S, layer);
	if(0 != ret){
		HandleSolutionErrorCode(L, "GetAmplitudes", ret);
		return 0;
	}
	
	n = Simulation_GetNumG(S, &G);
	if(NULL == G){
		return 0;
	}
	
	amp = (double*)malloc(sizeof(double)*8*n);
	Simulation_GetAmplitudes(S, 
		layer,
		luaL_checknumber(L, 3),
		amp, &amp[4*n]);
	
	lua_createtable(L, 2*n, 0);
	for(i = 0; i < 2*n; ++i){
		lua_pushinteger(L, i+1);
		/* push a complex number */
		lua_createtable(L, 2, 0);
		for(k = 0; k < 2; ++k){
			lua_pushinteger(L, k+1);
			lua_pushnumber(L, amp[4*n*0+2*i+k]);
			lua_settable(L, -3);
		}
		lua_settable(L, -3);
	}
	
	lua_createtable(L, 2*n, 0);
	for(i = 0; i < 2*n; ++i){
		lua_pushinteger(L, i+1);
		/* push a complex number */
		lua_createtable(L, 2, 0);
		for(k = 0; k < 2; ++k){
			lua_pushinteger(L, k+1);
			lua_pushnumber(L, amp[4*n*1+2*i+k]);
			lua_settable(L, -3);
		}
		lua_settable(L, -3);
	}
	
	free(amp);
	return 2;
}


/*
static int S4L_Simulation_EnableBasisFieldDump(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "EnableBasisFieldDump: 'Simulation' object expected.");
	
	S->options.do_vector_field_dump = 1;
	return 0;
}
*/
static int S4L_Simulation_UseDiscretizedEpsilon(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "UseDiscretizedEpsilon: 'Simulation' object expected.");
	
	if(lua_gettop(L) > 1){
		S->options.use_discretized_epsilon = lua_toboolean(L, 2);
	}else{
		S->options.use_discretized_epsilon = 1;
	}
	return 0;
}
static int S4L_Simulation_UseSubpixelSmoothing(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "UseSubpixelSmoothing: 'Simulation' object expected.");
	
	if(lua_gettop(L) > 1){
		S->options.use_subpixel_smoothing = lua_toboolean(L, 2);
	}else{
		S->options.use_subpixel_smoothing = 1;
	}
	return 0;
}
static int S4L_Simulation_UseLanczosSmoothing(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "UseLanczosSmoothing: 'Simulation' object expected.");
	
	if(lua_gettop(L) > 1){
		S->options.use_Lanczos_smoothing = lua_toboolean(L, 2);
	}else{
		S->options.use_Lanczos_smoothing = 1;
	}
	return 0;
}
static int S4L_Simulation_SetLanczosSmoothingWidth(lua_State *L){
	double prevwidth, width;
	int prevpower, power = 1;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "S4L_Simulation_SetLanczosSmoothingWidth: 'Simulation' object expected.");
	width = luaL_checknumber(L, 2);
	if(lua_gettop(L) > 2){
		power = luaL_checkint(L, 3);
	}
	
	prevwidth = S->options.lanczos_smoothing_width;
	prevpower = S->options.lanczos_smoothing_power;
	
	if(lua_gettop(L) > 2){
		S->options.lanczos_smoothing_width = width;
		S->options.lanczos_smoothing_power = power;
		lua_pushnumber(L, prevwidth);
		lua_pushinteger(L, prevpower);
		return 2;
	}else{
		S->options.lanczos_smoothing_width = width;
		lua_pushnumber(L, prevwidth);
		return 1;
	}
}
static int S4L_Simulation_UsePolarizationDecomposition(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "UsePolarizationDecomposition: 'Simulation' object expected.");
	
	if(lua_gettop(L) > 1){
		S->options.use_polarization_basis = lua_toboolean(L, 2);
	}else{
		S->options.use_polarization_basis = 1;
	}
	return 0;
}
static int S4L_Simulation_UseJonesVectorBasis(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "UseJonesVectorBasis: 'Simulation' object expected.");
	
	if(lua_gettop(L) > 1){
		S->options.use_jones_vector_basis = lua_toboolean(L, 2);
	}else{
		S->options.use_jones_vector_basis = 1;
	}
	return 0;
}
static int S4L_Simulation_UseNormalVectorBasis(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "UseNormalVectorBasis: 'Simulation' object expected.");
	
	if(lua_gettop(L) > 1){
		S->options.use_normal_vector_basis = lua_toboolean(L, 2);
	}else{
		S->options.use_normal_vector_basis = 1;
	}
	return 0;
}
static int S4L_Simulation_UseExperimentalFMM(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "UseExperimentalFMM: 'Simulation' object expected.");
	
	if(lua_gettop(L) > 1){
		S->options.use_experimental_fmm = lua_toboolean(L, 2);
	}else{
		S->options.use_experimental_fmm = 1;
	}
	return 0;
}

static int S4L_Simulation_SetBasisFieldDumpPrefix(lua_State *L){
	const char *prefix = NULL;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "SetBasisFieldDumpPrefix: 'Simulation' object expected.");
	
	if(lua_gettop(L) < 2 || LUA_TBOOLEAN == lua_type(L, 2)){
		if(NULL != S->options.vector_field_dump_filename_prefix){
			free(S->options.vector_field_dump_filename_prefix);
			S->options.vector_field_dump_filename_prefix = NULL;
		}
	}else{
		prefix = luaL_checklstring(L, 2, NULL);
		if(NULL != S->options.vector_field_dump_filename_prefix){
			free(S->options.vector_field_dump_filename_prefix);
		}
		S->options.vector_field_dump_filename_prefix = strdup(prefix);
	}
	
	return 0;
}

static int S4L_Simulation_SetLatticeTruncation(lua_State *L){
	const char *truncation_type;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "SetLatticeTruncation: 'Simulation' object expected.");
	truncation_type = luaL_checklstring(L, 2, NULL);
	
	if(0 == strcmp("Circular", truncation_type)){
		S->options.lattice_truncation = 0;
	}else if(0 == strcmp("Parallelogramic", truncation_type)){
		S->options.lattice_truncation = 1;
	}else{
		S4L_error(L, "SetLatticeTruncation: Truncation type must be Circular or Parallelogramic.");
	}
	
	return 0;
}

static int S4L_Simulation_SetVerbosity(lua_State *L){
	int verbosity;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "SetVerbosity: 'Simulation' object expected.");
	verbosity = luaL_checkint(L, 2);
	
	if(verbosity < 0){ verbosity = 0; }
	if(verbosity > 9){ verbosity = 9; }
	S->options.verbosity = verbosity;
	
	return 0;
}

static int S4L_Simulation_UseLessMemory(lua_State *L){
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "UseLessMemory: 'Simulation' object expected.");
	
	if(lua_gettop(L) > 1){
		S->options.use_less_memory = lua_toboolean(L, 2);
	}else{
		S->options.use_less_memory = 1;
	}
	return 0;
}

static int S4L_Simulation_OutputStructurePOVRay(lua_State *L){
	int ret;
	const char *filename;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "OutputStructurePOVRay: 'Simulation' object expected.");
	filename = luaL_optstring(L, 2, NULL);
	
	{
		FILE *fp = stdout;
		if(NULL != filename){
			fp = fopen(filename, "wb");
		}
		
		ret = Simulation_OutputStructurePOVRay(S, fp);
		if(0 != ret){
			HandleSolutionErrorCode(L, "OutputStructurePOVRay", ret);
		}
		
		if(NULL != filename){
			fclose(fp);
		}
	}
	return 0;
}

static int S4L_Simulation_OutputLayerPatternDescription(lua_State *L){
	int ret;
	const char *layer_name, *filename;
	Layer *layer;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "OutputLayerPatternDescription: 'Simulation' object expected.");
	
	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "OutputLayerPatternDescription: Layer named '%s' not found.", layer_name);
		return 0;
	}
	filename = luaL_optstring(L, 3, NULL);
	
	{
		FILE *fp = stdout;
		if(NULL != filename){
			fp = fopen(filename, "wb");
		}
		
		ret = Simulation_OutputLayerPatternDescription(S, layer, fp);
		if(0 != ret){
			HandleSolutionErrorCode(L, "OutputLayerPatternDescription", ret);
		}
		
		if(NULL != filename){
			fclose(fp);
		}
	}
	return 0;
}
static int S4L_Simulation_OutputLayerPatternRealization(lua_State *L){
	int ret;
	const char *layer_name, *filename;
	Layer *layer;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "OutputLayerPatternRealization: 'Simulation' object expected.");
	
	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "OutputLayerPatternRealization: Layer named '%s' not found.", layer_name);
		return 0;
	}
	filename = luaL_optstring(L, 5, NULL);
	
	{
		FILE *fp = stdout;
		if(NULL != filename){
			fp = fopen(filename, "wb");
		}
		
		ret = Simulation_OutputLayerPatternRealization(S,
			layer,
			luaL_checkint(L, 3),
			luaL_checkint(L, 4),
			fp);
		if(0 != ret){
			HandleSolutionErrorCode(L, "OutputLayerPatternRealization", ret);
		}
		
		if(NULL != filename){
			fclose(fp);
		}
	}
	return 0;
}
static int S4L_Simulation_GetEField(lua_State *L){
	int j, ret;
	double r[3], fE[6];
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "GetEField: 'Simulation' object expected.");
	
	for(j = 0; j < 3; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 2);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "GetEField: Point coordinate must be numeric.");
		}else{
			r[j] = (double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	
	ret = Simulation_GetField(S, r, fE, NULL);
	if(0 != ret){
		HandleSolutionErrorCode(L, "GetEField", ret);
	}
	
	for(j = 0; j < 6; ++j){
		lua_pushnumber(L, fE[j]);
	}
	return 6;
}
static int S4L_Simulation_GetHField(lua_State *L){
	int j, ret;
	double r[3], fH[6];
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "GetHField: 'Simulation' object expected.");
	
	for(j = 0; j < 3; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 2);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "GetHField: Point coordinate must be numeric.");
		}else{
			r[j] = (double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	
	ret = Simulation_GetField(S, r, NULL, fH);
	if(0 != ret){
		HandleSolutionErrorCode(L, "GetHField", ret);
	}
	
	for(j = 0; j < 6; ++j){
		lua_pushnumber(L, fH[j]);
	}
	return 6;
}
static int S4L_Simulation_GetFields(lua_State *L){
	int j, ret;
	double r[3], fE[6],fH[6];
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "GetFields: 'Simulation' object expected.");
	
	for(j = 0; j < 3; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 2);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "GetFields: Point coordinate must be numeric.");
		}else{
			r[j] = (double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	
	ret = Simulation_GetField(S, r, fE, fH);
	if(0 != ret){
		HandleSolutionErrorCode(L, "GetFields", ret);
	}
	
	for(j = 0; j < 3; ++j){
		lua_pushnumber(L, fE[j]);
	}
	for(j = 0; j < 3; ++j){
		lua_pushnumber(L, fH[j]);
	}
	for(j = 0; j < 3; ++j){
		lua_pushnumber(L, fE[3+j]);
	}
	for(j = 0; j < 3; ++j){
		lua_pushnumber(L, fH[3+j]);
	}
	return 12;
}
static int S4L_Simulation_GetSMatrixDeterminant(lua_State *L){
	int ret;
	double mant[2], base;
	int expo;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "GetSMatrixDeterminant: 'Simulation' object expected.");
	
	ret = Simulation_GetSMatrixDeterminant(S, mant, &base, &expo);
	if(0 != ret){
		HandleSolutionErrorCode(L, "GetSMatrixDeterminant", ret);
	}
	
	lua_pushnumber(L, mant[0]);
	lua_pushnumber(L, mant[1]);
	lua_pushnumber(L, base);
	lua_pushnumber(L, expo);

	return 4;
}
static int S4L_Simulation_GetEpsilon(lua_State *L){
	int j, ret;
	double r[3], feps[2];
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "GetEpsilon: 'Simulation' object expected.");
	
	for(j = 0; j < 3; ++j){
		lua_pushinteger(L, 1+j);
		lua_gettable(L, 2);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "GetEpsilon: Point coordinate must be numeric.");
		}else{
			r[j] = (double)lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
	}
	
	ret = Simulation_GetEpsilon(S, r, feps);
	if(0 != ret){
		HandleSolutionErrorCode(L, "GetEpsilon", ret);
	}
	
	lua_pushnumber(L, feps[0]);
	lua_pushnumber(L, feps[1]);
	return 2;
}


/* Expected stack:
 *   1 Simulation
 *   2 layer name string
 *   3 layer offset
 */
static int S4L_Simulation_GetStressTensorIntegral(lua_State *L){
	int ret;
	double Tint[6];
	const char *layer_name;
	Layer *layer;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "GetStressTensorIntegral: 'Simulation' object expected.");
	
	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "GetStressTensorIntegral: Layer named '%s' not found.", layer_name);
		return 0;
	}
	
	ret = Simulation_GetStressTensorIntegral(S, layer, luaL_checknumber(L, 3), Tint);
	if(0 != ret){
		HandleSolutionErrorCode(L, "GetStressTensorIntegral", ret);
	}
	
	lua_pushnumber(L, Tint[0]);
	lua_pushnumber(L, Tint[1]);
	lua_pushnumber(L, Tint[2]);
	lua_pushnumber(L, Tint[3]);
	lua_pushnumber(L, Tint[4]);
	lua_pushnumber(L, Tint[5]);
	return 6;
}

/* Expected stack:
 *   1 Simulation
 *   2 layer name string
 */
static int S4L_Simulation_GetLayerVolumeIntegral(lua_State *L, char which, const char *name){
	int ret;
	double integral[2];
	const char *layer_name;
	Layer *layer;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "Get*LayerIntegral: 'Simulation' object expected.");
	
	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "%s: Layer named '%s' not found.", name, layer_name);
		return 0;
	}
	
	ret = Simulation_GetLayerVolumeIntegral(S, layer, which, integral);
	if(0 != ret){
		HandleSolutionErrorCode(L, name, ret);
	}
	
	lua_pushnumber(L, integral[0]);
	lua_pushnumber(L, integral[1]);
	return 2;
}
static int S4L_Simulation_GetLayerEnergyDensityIntegral(lua_State *L){
	return S4L_Simulation_GetLayerVolumeIntegral(L, 'U', "GetLayerEnergyDensityIntegral");
}
static int S4L_Simulation_GetLayerElectricEnergyDensityIntegral(lua_State *L){
	return S4L_Simulation_GetLayerVolumeIntegral(L, 'E', "GetLayerElectricEnergyDensityIntegral");
}
static int S4L_Simulation_GetLayerMagneticEnergyDensityIntegral(lua_State *L){
	return S4L_Simulation_GetLayerVolumeIntegral(L, 'H', "GetLayerMagneticEnergyDensityIntegral");
}
static int S4L_Simulation_GetLayerElectricFieldIntensityIntegral(lua_State *L){
	return S4L_Simulation_GetLayerVolumeIntegral(L, 'e', "GetLayerElectricFieldIntensityIntegral");
}

static int S4L_Simulation_GetLayerZIntegral(lua_State *L){
	int ret, i;
	double integral[6], r[2];
	const char *layer_name;
	Layer *layer;
	Simulation *S = (Simulation *)luaL_checkudata(L, 1, "S4.Simulation");
	luaL_argcheck(L, S != NULL, 1, "GetLayerZIntegral: 'Simulation' object expected.");
	
	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "GetLayerZIntegral: Layer named '%s' not found.", layer_name);
		return 0;
	}
	
	luaL_checktype(L, 3, LUA_TTABLE);
	luaL_argcheck(L, lua_rawlen(L, 3) == 2, 3, "Position must be a pair of coordinates.");
	for(i = 0; i < 2; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, 3);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "Position must be a pair of coordinates.");
		}
		r[i] = lua_tonumber(L, -1);
		lua_pop(L, 1);
	}
	
	ret = Simulation_GetLayerZIntegral(S, layer, r, integral);
	if(0 != ret){
		HandleSolutionErrorCode(L, "GetLayerZIntegral", ret);
	}
	
	lua_pushnumber(L, integral[0]);
	lua_pushnumber(L, integral[1]);
	lua_pushnumber(L, integral[2]);
	lua_pushnumber(L, integral[3]);
	lua_pushnumber(L, integral[4]);
	lua_pushnumber(L, integral[5]);
	return 6;
}

static int S4_openlib(lua_State *L){
	static const struct luaL_Reg S4_lib[] = {
		{"NewSimulation", S4L_NewSimulation},
		{"NewSpectrumSampler", S4L_NewSpectrumSampler},
		{"NewInterpolator", S4L_NewInterpolator},
		{"SolveInParallel", S4L_SolveInParallel},
		{"ConvertUnits", S4L_ConvertUnits},
		{"Integrate", S4L_Integrate},
		{NULL, NULL}
	};
	static const struct luaL_Reg SpectrumSamplerObj[] = {
		{"IsDone", S4L_SpectrumSampler_IsDone},
		{"GetFrequency", S4L_SpectrumSampler_GetFrequency},
		{"SubmitResult", S4L_SpectrumSampler_SubmitResult},
		{"GetFrequencies", S4L_SpectrumSampler_GetFrequencies},
		{"SubmitResults", S4L_SpectrumSampler_SubmitResults},
		{"GetSpectrum", S4L_SpectrumSampler_GetSpectrum},
		{NULL, NULL}
	};
	static const struct luaL_Reg InterpolatorObj[] = {
		{"Get", S4L_Interpolator_Get},
		{NULL, NULL}
	};
	static const struct luaL_Reg SimulationObj[] = {
		{"Clone", S4L_Simulation_Clone},
		{"SetLattice", S4L_Simulation_SetLattice},
		{"SetNumG", S4L_Simulation_SetNumG},
		{"GetNumG", S4L_Simulation_GetNumG},
		{"SetResolution", S4L_Simulation_SetResolution},
		{"AddMaterial", S4L_Simulation_AddMaterial},
		{"SetMaterial", S4L_Simulation_SetMaterial},
		{"AddLayer", S4L_Simulation_AddLayer},
		{"SetLayer", S4L_Simulation_SetLayer},
		{"SetLayerThickness", S4L_Simulation_SetLayerThickness},
		{"AddLayerCopy", S4L_Simulation_AddLayerCopy},
		{"SetLayerPatternCircle", S4L_Simulation_SetLayerPatternCircle},
		{"SetLayerPatternEllipse", S4L_Simulation_SetLayerPatternEllipse},
		{"SetLayerPatternRectangle", S4L_Simulation_SetLayerPatternRectangle},
		{"SetLayerPatternPolygon", S4L_Simulation_SetLayerPatternPolygon},
		{"SetExcitationPlanewave", S4L_Simulation_SetExcitationPlanewave},
		{"SetExcitationDipoleK", S4L_Simulation_SetExcitationDipole},
		{"SetExcitationExterior", S4L_Simulation_SetExcitationExterior},
		{"SetExcitationInterior", S4L_Simulation_SetExcitationInterior},
		{"SetFrequency", S4L_Simulation_SetFrequency},
		/* Outputs requiring solutions */
		{"GetGList", S4L_Simulation_GetGList},
		{"GetDiffractionOrder", S4L_Simulation_GetDiffractionOrder},
		{"GetPoyntingFlux", S4L_Simulation_GetPoyntingFlux},
		{"GetPoyntingFluxByOrder", S4L_Simulation_GetPoyntingFluxByOrder},
		{"GetPowerFlux", S4L_Simulation_GetPoyntingFlux}, /* alias */
		{"GetPowerFluxByOrder", S4L_Simulation_GetPoyntingFluxByOrder}, /* alias */
		{"GetAmplitudes", S4L_Simulation_GetAmplitudes},
		{"GetStressTensorIntegral", S4L_Simulation_GetStressTensorIntegral},
		{"GetLayerEnergyDensityIntegral", S4L_Simulation_GetLayerEnergyDensityIntegral},
		{"GetLayerElectricEnergyDensityIntegral", S4L_Simulation_GetLayerElectricEnergyDensityIntegral},
		{"GetLayerMagneticEnergyDensityIntegral", S4L_Simulation_GetLayerMagneticEnergyDensityIntegral},
		{"GetLayerElectricFieldIntensityIntegral", S4L_Simulation_GetLayerElectricFieldIntensityIntegral},
		{"GetLayerZIntegral", S4L_Simulation_GetLayerZIntegral},
		{"GetEField", S4L_Simulation_GetEField},
		{"GetHField", S4L_Simulation_GetHField},
		{"GetFields", S4L_Simulation_GetFields},
		{"GetSMatrixDeterminant", S4L_Simulation_GetSMatrixDeterminant},
		/* Outputs not requiring solutions */
		{"GetReciprocalLattice", S4L_Simulation_GetReciprocalLattice},
		{"GetEpsilon", S4L_Simulation_GetEpsilon},
		{"OutputStructurePOVRay", S4L_Simulation_OutputStructurePOVRay},
		{"OutputLayerPatternDescription", S4L_Simulation_OutputLayerPatternDescription},
		{"OutputLayerPatternRealization", S4L_Simulation_OutputLayerPatternRealization},
		/* Options */
		{"UseDiscretizedEpsilon", S4L_Simulation_UseDiscretizedEpsilon},
		{"UseSubpixelSmoothing", S4L_Simulation_UseSubpixelSmoothing},
		{"UseLanczosSmoothing", S4L_Simulation_UseLanczosSmoothing},
		{"SetLanczosSmoothingWidth", S4L_Simulation_SetLanczosSmoothingWidth},
		{"UsePolarizationDecomposition", S4L_Simulation_UsePolarizationDecomposition},
		{"UseJonesVectorBasis", S4L_Simulation_UseJonesVectorBasis},
		{"UseNormalVectorBasis", S4L_Simulation_UseNormalVectorBasis},
		{"SetBasisFieldDumpPrefix", S4L_Simulation_SetBasisFieldDumpPrefix},
		{"SetLatticeTruncation", S4L_Simulation_SetLatticeTruncation},
		{"SetVerbosity", S4L_Simulation_SetVerbosity},
		{"UseExperimentalFMM", S4L_Simulation_UseExperimentalFMM},
		{"UseLessMemory", S4L_Simulation_UseLessMemory},
		{NULL, NULL}
	};
	
	luaL_newlib(L, S4_lib);

	luaL_newmetatable(L, "S4.Simulation");
	luaL_newlib(L, SimulationObj);
	lua_setfield(L, -2, "__index");
	lua_pushstring(L, "__gc");
	lua_pushcfunction(L, &S4L_Simulation__gc);
	lua_settable(L, -3);
	lua_pop(L, 1);
	
	luaL_newmetatable(L, "S4.SpectrumSampler");
	luaL_newlib(L, SpectrumSamplerObj);
	lua_setfield(L, -2, "__index");
	lua_pushstring(L, "__gc");
	lua_pushcfunction(L, S4L_SpectrumSampler__gc);
	lua_settable(L, -3);
	lua_pop(L, 1);
	
	luaL_newmetatable(L, "S4.Interpolator");
	luaL_newlib(L, InterpolatorObj);
	lua_setfield(L, -2, "__index");
	lua_pushstring(L, "__gc");
	lua_pushcfunction(L, S4L_Interpolator__gc);
	lua_settable(L, -3);
	lua_pop(L, 1);
	
	return 1;
}

void threadsafe_init(){
	extern void exactinit();
	extern float slamch_(const char *);
	extern double dlamch_(const char *);
	exactinit();
	slamch_("E");
	dlamch_("E");
	fft_init();
}
void threadsafe_destroy(){
	fft_destroy();
}

void usage(){
	printf("S4 [-a arg] [-h] [-t thread-count] [-v] [input-file]\n");
}
void version(){
	printf("Stanford Stratified Structure Solver (S4)\n");
	printf("Version %s\n", PACKAGE_VERSION);
#ifdef HAVE_MPI
	printf("  With MPI support\n");
#endif
}

/* from unistd.h */
extern char   *optarg;
extern int    optind, opterr, optopt;
int getopt(int argc, char * const argv[], const char *optstring);
#include <ctype.h>

int main(int argc, char *argv[]){
	char buff[256];
	int c;
	int index, error;
	unsigned int max_threads = get_max_threads();
	char *arg = NULL;
	int mpi_size = 1, mpi_rank = 0;
	lua_State *L;
	
#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
	
	opterr = 0;
	while((c = getopt(argc, argv, "a:ht:v")) != -1){
		switch(c){
		case 'a':
			arg = strdup(optarg);
			break;
		case 'h':
			usage();
			return EXIT_SUCCESS;
		case 't':
			max_threads = atoi(optarg);
			break;
		case 'v':
			version();
			return EXIT_SUCCESS;
		case '?':
			if('t' == optopt){
				fprintf(stderr, "Option -%c requires an argument.\n", optopt);
			}else if(isprint(optopt)){
				fprintf(stderr, "Unknown option -%c.\n", optopt);
			}else{
				fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
			}
			usage();
			return EXIT_FAILURE;
		default:
			abort();
		}
	}
	
	L = luaL_newstate(); /* opens Lua */
	
	threadsafe_init();
	
	luaL_requiref(L, "S4", &S4_openlib, 1);
	lua_pop(L, 1);
	
	luaL_openlibs(L); /* opens the standard libraries */

	S4_threads_init(L, max_threads);
	
	/* Set the argument if there is one */
	if(NULL != arg){
		lua_getglobal(L, "S4");
		lua_pushstring(L, "arg");
		lua_pushstring(L, arg);
		lua_settable(L, -3);
		lua_pop(L, 1);
	}
	
	lua_getglobal(L, "S4");
	lua_pushstring(L, "MPIRank");
	lua_pushinteger(L, mpi_rank);
	lua_settable(L, -3);
	lua_pushstring(L, "MPISize");
	lua_pushinteger(L, mpi_size);
	lua_settable(L, -3);
	lua_pop(L, 1);
	
	if(optind < argc){ /* has at least 1 argument */
		S4L_set_interactive(L, 0);
		
		for(index = optind; index < argc; ++index){
			error = luaL_loadfile(L, argv[index]);
			if(error){
				fprintf(stderr, "%s\n", lua_tostring(L, -1));
				lua_pop(L, 1); /* pop error message from the stack */
			}else{
				error = lua_pcall(L, 0, 0, 0);
				if(error){
					fprintf(stderr, "%s\n", lua_tostring(L, -1));
					lua_pop(L, 1); /* pop error message from the stack */
				}
			}
		}
	}else{ /* run in REPL mode */
		fprintf(stdout, "No input file given, running in interactive mode\n");
		S4L_set_interactive(L, 1);
		while(fgets(buff, sizeof(buff), stdin) != NULL){
			error = luaL_loadbuffer(L, buff, strlen(buff), "line") || lua_pcall(L, 0, 0, 0);
			if(error){
				fprintf(stderr, "%s", lua_tostring(L, -1));
				lua_pop(L, 1); /* pop error message from the stack */
			}
		}
	}
	
	S4_threads_destroy(L);
	
	lua_close(L);
	
	threadsafe_destroy();
	
	if(NULL != arg){ free(arg); }
#ifdef HAVE_MPI
	MPI_Finalize();
#endif
	return EXIT_SUCCESS;
}
