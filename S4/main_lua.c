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


#ifndef LUA_OK
#define LUA_OK 0
#endif

#if LUA_VERSION_NUM > 501
#define LUA_SETFUNCS(L, R) luaL_setfuncs(L, R, 0)
#else
#define LUA_SETFUNCS(L, R) luaL_register(L, NULL, R)
#define lua_rawlen lua_objlen
#endif

#ifndef luaL_newlibtable
#define luaL_newlibtable(L,l)	\
  lua_createtable(L, 0, sizeof(l)/sizeof((l)[0]) - 1)
#endif

#ifndef luaL_newlib
#define luaL_newlib(L,l)  \
  (luaL_newlibtable(L,l), LUA_SETFUNCS(L,l))
#endif

lua_State *new_S4_lua_state();

void fft_init();
void fft_destroy();

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
pthread_mutex_t g_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif
#ifdef S4_DEBUG
# include "debug.h"
# ifdef HAVE_LIBPTHREAD
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
/*
	if(!pthread_win32_process_attach_np()){
		fprintf(stderr, "Pthread initialization failed.\n");
		return;
	}
*/
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
//	pthread_win32_process_detach_np();
#endif
}
void S4_threads_join(lua_State *L, int n){
	int i, nthreads = n;
	pthread_t *thread;
pthread_mutex_lock(&g_mutex);
	lua_pushlightuserdata(L, (void *)&thread_count_key);
	lua_gettable(L, LUA_REGISTRYINDEX);
	nthreads = lua_tointeger(L, -1);
	lua_pop(L, 1);
	lua_pushlightuserdata(L, (void *)&thread_key);
	lua_gettable(L, LUA_REGISTRYINDEX);
	thread = (pthread_t*)lua_topointer(L, -1);
	lua_pop(L, 1);
	if(n < nthreads){ nthreads = n; }
pthread_mutex_unlock(&g_mutex);
	for(i = 0; i < nthreads; ++i){
//printf("In S4_threads_join, joining %d\n", i);
		pthread_join(thread[i], NULL);
	}
}
void S4_threads_run(lua_State *L, int i, void* (*func)(void*), void *data){
	int nthreads;
	pthread_t *thread;
pthread_mutex_lock(&g_mutex);
	lua_pushlightuserdata(L, (void *)&thread_count_key);
	lua_gettable(L, LUA_REGISTRYINDEX);
	nthreads = lua_tointeger(L, -1);
	lua_pop(L, 1);
	lua_pushlightuserdata(L, (void *)&thread_key);
	lua_gettable(L, LUA_REGISTRYINDEX);
	thread = (pthread_t*)lua_topointer(L, -1);
	lua_pop(L, 1);
pthread_mutex_unlock(&g_mutex);
//printf("In S4_threads_run, running on thread %d/%d\n", i, nthreads);
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

const char S4_simulation_typename[] = "S4.S4_Simulation";
S4_Simulation *S4L_get_simulation(lua_State *L, int index){
	S4_Simulation **pS = (S4_Simulation **)luaL_checkudata(L, index, S4_simulation_typename);
	if(NULL == pS){ return NULL; }
	return *pS;
}



typedef struct S4_solve_in_parallel_data_{
	S4_Simulation *S;
	S4_LayerID layer;
} S4_solve_in_parallel_data;

void* S4_solve_in_parallel(void *p){
	S4_solve_in_parallel_data *data = (S4_solve_in_parallel_data*)p;
	S4_Simulation_SolveLayer(data->S, data->layer);
	return NULL;
}
/* Expected stack:
 *   1 layer name string
 *   2 S4_Simulation object
 *   : S4_Simulation objects
 */
static int S4L_SolveInParallel(lua_State *L){
	int i, n;
	S4_solve_in_parallel_data *data;
	const char *layer_name = luaL_checklstring(L, 1, NULL);

	n = lua_gettop(L);
	data = (S4_solve_in_parallel_data*)malloc(sizeof(S4_solve_in_parallel_data)*(n-1));
	for(i = 2; i <= n; ++i){
		S4_LayerID layer;
		S4_Simulation *S = S4L_get_simulation(L, i);
		luaL_argcheck(L, S != NULL, i, "SolveInParallel: 'S4_Simulation' object expected.");
		layer = S4_Simulation_GetLayerByName(S, layer_name);
		if(layer < 0){
			S4L_error(L, "SolveInParallel: S4_Layer named '%s' not found.", layer_name);
		}
		data[i-2].S = S;
		data[i-2].layer = layer;

		S4_threads_run(L, i-2, &S4_solve_in_parallel, (void*)&data[i-2]);
	}
	S4_threads_join(L, n-1);
	free(data);

	return 0;
}


static void deep_copy_1(lua_State *Lfrom, int index, lua_State *Lto, int level){
	if(level > 8){
		luaL_error(Lfrom, "Too many levels of recursion in deep copy");
	}
	switch(lua_type(Lfrom, index)){
	case LUA_TNUMBER:
		lua_pushnumber(Lto, lua_tonumber(Lfrom, index));
		break;
	case LUA_TBOOLEAN:
		lua_pushboolean(Lto, lua_toboolean(Lfrom, index));
		break;
	case LUA_TSTRING: {
		size_t length;
		const char *string = lua_tolstring(Lfrom, index, &length);
		lua_pushlstring(Lto, string, length);
		break;
	}
	case LUA_TLIGHTUSERDATA: {
		lua_pushlightuserdata(Lto, lua_touserdata(Lfrom, index));
		break;
	}
	case LUA_TNIL:
		lua_pushnil(Lto);
		break;
	case LUA_TTABLE:
		/* make sure there is room on the new state for 3 values
		 * (table,key,value) */
		if (!lua_checkstack(Lto, 3)) {
			luaL_error(Lfrom, "To stack overflow");
		}
		/* make room on from stack for key/value pairs */
		luaL_checkstack(Lfrom, 2, "From stack overflow");
		lua_newtable(Lto);
		lua_pushnil(Lfrom);
		while(lua_next(Lfrom, index) != 0){
			/* key is at (top - 1), value at (top), but we need to normalize
			 * these to positive indices */
			int kv_pos = lua_gettop(Lfrom);
			deep_copy_1(Lfrom, kv_pos - 1, Lto, level+1);
			deep_copy_1(Lfrom, kv_pos    , Lto, level+1);
			/* Copied key and value are now at -2 and -1 in dest */
			lua_settable(Lto, -3);
			/* Pop value for next iteration */
			lua_pop(Lfrom, 1);
		}
		break;
	case LUA_TFUNCTION:
	case LUA_TUSERDATA:
	case LUA_TTHREAD:
	default:
		lua_pushfstring(Lto, "Unsupported value: %s: %p",
			lua_typename(Lfrom, lua_type(Lfrom, index)),
			lua_topointer(Lfrom, index)
		);
	}
}
static void deep_copy(int unpack, lua_State *Lfrom, int index, lua_State *Lto){
	if(unpack){
		int i, n;
		n = lua_rawlen(Lfrom, index);
		for(i = 0; i < n; ++i){
			lua_pushinteger(Lfrom, i+1);
			lua_gettable(Lfrom, index);
			deep_copy_1(Lfrom, lua_gettop(Lfrom), Lto, 0);
			lua_pop(Lfrom, 1);
		}
	}else{
		deep_copy_1(Lfrom, index, Lto, 0);
	}
}

/* ParallelInvoke(S_list, func_name, params_list)
 *   S_list:      List of S4_Simulation objects
 *   func_name:   String of function to invoke
 *   params_list: List of parameter lists to each S4_Simulation object
 * If len(params_list) < len(S_list), then
 *    S_list[i] is called with param_list[i%len(S_list)] (for 0 based indexing)
 *
 * Second attempt: Spin up a new lua_State for each thread!
 */

typedef struct{
	int i; // index of S in S_list, 0 based
	int iargs; // index of arguments
	lua_State *L; // master state
	S4_Simulation *S;
	const char *func_name;
} ParallelInvokeData;
static void stackDump (lua_State *L) {
      int i;
      int top = lua_gettop(L);
      for (i = 1; i <= top; i++) {  /* repeat for each level */
        int t = lua_type(L, i);
        switch (t) {

          case LUA_TSTRING:  /* strings */
            printf("`%s'", lua_tostring(L, i));
            break;

          case LUA_TBOOLEAN:  /* booleans */
            printf(lua_toboolean(L, i) ? "true" : "false");
            break;

          case LUA_TNUMBER:  /* numbers */
            printf("%g", lua_tonumber(L, i));
            break;

          default:  /* other values */
            printf("%s", lua_typename(L, t));
            break;

        }
        printf("  ");  /* put a separator */
      }
      printf("\n");  /* end the listing */
    }
void* ParallelInvoke_func(void *data){
	int j;
	ParallelInvokeData *pid = (ParallelInvokeData*)data;
	lua_State *Lm = pid->L;
	int iS, nargs = 0;

//printf("thread %d enter, master Lua stack top = %d\n", pid->i, lua_gettop(Lm));

	lua_State *L = new_S4_lua_state();
	// Add the existing simulation object to the new state
	S4_Simulation **pS = (S4_Simulation **)lua_newuserdata(L, sizeof(S4_Simulation*));
	luaL_getmetatable(L, "S4.S4_Simulation");
	lua_setmetatable(L, -2);
	*pS = pid->S;
	//S4_Simulation *S = *pS;

	// This is where things are tricky...
	{ // First set up the function call; each thread needs exclusive access to the Lua state
#ifdef HAVE_LIBPTHREAD
		pthread_mutex_lock(&g_mutex);
#endif
		// Current stack: ... S
		iS = lua_gettop(L);

		// First get the function
		lua_getmetatable(L, -1);
		lua_pushstring(L, "__index");
		lua_gettable(L, -2);
		/*{
			lua_pushnil(L);
			 while (lua_next(L, -2) != 0) {
				printf("%s - %s\n",
				lua_typename(L, lua_type(L, -2)),
				lua_typename(L, lua_type(L, -1)));
				printf("  key = %s\n", lua_tostring(L, -2));
				// removes 'value'; keeps 'key' for next iteration
				lua_pop(L, 1);
			}
		}*/

		lua_pushstring(L, pid->func_name);
		lua_gettable(L, -2);
		if(lua_isnil(L, -1)){
			S4L_error(Lm, "ParallelInvoke: function '%s' not found.", pid->func_name);
			lua_pop(L, 2);
#ifdef HAVE_LIBPTHREAD
			pthread_mutex_unlock(&g_mutex);
#endif
			return NULL;
		}
		// Current stack: ... S meta __index func
		lua_insert(L, -4);
		lua_pop(L, 2);
		// Current stack: ... func S

		// Now pull the args, assumed to be in a table at stack index 3
		if(pid->iargs >= 0){
			if(lua_istable(Lm, 3)){
				lua_pushinteger(Lm, pid->iargs + 1);
				lua_gettable(Lm, 3);

				// Current stack: ... func S [table of args]
				if(lua_istable(Lm, -1)){
					nargs = lua_rawlen(Lm, -1);

					deep_copy(1, Lm, lua_gettop(Lm), L);
				}else{
					nargs = 1;
					deep_copy(0, Lm, lua_gettop(Lm), L);
				}
				lua_pop(Lm, 1);
			}else{
				nargs = 1;
				deep_copy(0, Lm, 3, L);
			}
		}
		// Current stack: ... func S arg1 arg2 ...
#ifdef HAVE_LIBPTHREAD
		pthread_mutex_unlock(&g_mutex);
#endif
	}
	int ret = lua_pcall(L, nargs+1, LUA_MULTRET, 0);
	if(LUA_OK != ret){
#ifdef HAVE_LIBPTHREAD
		pthread_mutex_lock(&g_mutex);
#endif
		S4L_error(Lm, "ParallelInvoke: error in function: %s.", lua_tostring(L, -1));
		lua_pop(L, 1);
#ifdef HAVE_LIBPTHREAD
		pthread_mutex_unlock(&g_mutex);
#endif
	}

	{ // Now copy all results back
#ifdef HAVE_LIBPTHREAD
		pthread_mutex_lock(&g_mutex);
#endif
		int nret = lua_gettop(L) - iS;
		lua_createtable(Lm, nret, 0);
		for(j = 0; j < nret; ++j){
			deep_copy(0, L, j+1, Lm);
			lua_rawseti(Lm, -2, j+1);
		}
		lua_rawseti(Lm, 4, pid->i+1);

#ifdef HAVE_LIBPTHREAD
		pthread_mutex_unlock(&g_mutex);
#endif
	}

//printf("thread %d exit, master Lua stack top = %d\n", pid->i, lua_gettop(Lm));
	return NULL;
}
int S4L_ParallelInvoke(lua_State *L){
	int i;
	ParallelInvokeData *data = NULL;

//printf("Entering ParallelInvoke\n"); fflush(stdout);
	//// See if we have arguments
	int nargs = 0;
	{
		if(lua_gettop(L) < 3){
			lua_pushnil(L);
			nargs = 0;
		}else if(lua_istable(L, 3)){
			nargs = lua_rawlen(L, 3);
		}else{
			nargs = 1;
		}
	}

	//// Get vector of S4_Simulation pointers
	int nS = 0; // number of S objects
	{
		luaL_checktype(L, 1, LUA_TTABLE);
		nS = lua_rawlen(L, 1);
		if(0 == nS){ return 0; }
		data = (ParallelInvokeData*)malloc(sizeof(ParallelInvokeData) * nS);
		for(i = 0; i < nS; ++i){
			lua_pushinteger(L, i+1);
			lua_gettable(L, 1);
			data[i].i = i;
			data[i].L = L;
			data[i].S = S4L_get_simulation(L, -1);
			if(nargs > 0){
				data[i].iargs = i % nargs;
			}else{
				data[i].iargs = -1;
			}
			luaL_argcheck(L, data[i].S != NULL, 1, "ParallelInvoke: 'S4_Simulation' object expected.");
			lua_pop(L, 1);
		}
	}

	//// Pull in the function name
	char *func_name = strdup(luaL_checkstring(L, 2));
	for(i = 0; i < nS; ++i){
		data[i].func_name = func_name;
	}

	// also need to init the mutex

	// Create the return table
	lua_createtable(L, nS, 0);

	//// Spawn threads (this section is just an outline)
	for(i = 0; i < nS; ++i){
		S4_threads_run(L, i, &ParallelInvoke_func, &data[i]);
	}
//printf("called S4_threads_run\n"); fflush(stdout);
	S4_threads_join(L, nS);
//printf("called S4_threads_join\n"); fflush(stdout);

	free(func_name);
	free(data);

//printf("Exiting ParallelInvoke\n"); fflush(stdout);
	return 1;
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

struct IntegrateFunction_data{
	lua_State *L;
	int VectorArgument;
};
static void IntegrateFunction(
	unsigned ndim, const double *x, void *fdata,
	unsigned fdim, double *fval
){
	unsigned i;
	struct IntegrateFunction_data *ifd = (struct IntegrateFunction_data*)fdata;
	lua_State *L = ifd->L;

	fval[0] = 0;

	lua_pushvalue(L, -1); /* push a copy of the function */

	if(ifd->VectorArgument){
		unsigned j;
		lua_createtable(L, 1, 0);
		lua_pushinteger(L, 1);
		lua_createtable(L, ndim, 0);
		for(j = 0; j < ndim; ++j){
			lua_pushinteger(L, j+1);
			lua_pushnumber(L, x[j]);
			lua_settable(L, -3);
		}
		lua_settable(L, -3);
	}else{
		for(i = 0; i < ndim; ++i){
			lua_pushnumber(L, x[i]);
		}
	}
	if(0 != lua_pcall(L, ndim, fdim, 0)){
		S4L_error(L, "Error running Integrate function: %s", lua_tostring(L, -1));
        lua_pop(L, 1);
        return;
	}

	if(lua_isnumber(L, -1)){
		fval[0] = lua_tonumber(L, -1);
	}else if(lua_istable(L, -1)){
		lua_pushinteger(L, 1);
		lua_gettable(L, -2);
		if(lua_isnumber(L, -1)){
			fval[0] = lua_tonumber(L, -1);
		}else if(lua_istable(L, -1)){
			lua_pushinteger(L, 1);
			lua_gettable(L, -2);
			if(lua_isnumber(L, -1)){
				fval[0] = lua_tonumber(L, -1);
			}else{
				lua_pop(L, 1);
				S4L_error(L, "Integrate function must return a number or a nested table containing numbers");
				return;
			}
			lua_pop(L, 1);
		}else{
			lua_pop(L, 1);
			S4L_error(L, "Integrate function must return a number or a nested table containing numbers");
			return;
		}
		lua_pop(L, 1);
	}else{
		lua_pop(L, 1);
		S4L_error(L, "Integrate function must return a number or a nested table containing numbers");
		return;
	}
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
	int VectorArg = 0;

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
					}else if(0 == strcmp(key, "VectorArgument")){
						VectorArg = lua_toboolean(L, -1);
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
		struct IntegrateFunction_data data;
		data.L = L;
		data.VectorArgument = VectorArg;
		adapt_integrate(
			1, &IntegrateFunction, (void*)&data,
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
	S4_Simulation **pS;
	double Lr[4] = { 1, 0, 0, 1 };
	unsigned int nG = 1;
	unsigned int argflag = 0;

	lua_gc(L, LUA_GCCOLLECT, 0);
	if(lua_istable(L, 1)){
		lua_pushnil(L);
		while(0 != lua_next(L, 1)){
			const char *key;
			if(!lua_isstring(L, -2)){
				return luaL_argerror(L, 1, "Expected only named arguments");
			}
			key = lua_tostring(L, -2);
			if(0 == strcmp("Lattice", key)){
				if(lua_isnumber(L, -1)){
					Lr[0] = lua_tonumber(L, -1);
					Lr[1] = 0;
					Lr[2] = 0;
					Lr[3] = 0;
				}else if(lua_istable(L, -1)){
					int i, j;
					for(i = 0; i < 2; ++i){
						lua_pushinteger(L, i+1);
						lua_gettable(L, -2);
						if(!lua_istable(L, -1)){
							return luaL_argerror(L, 1, "Lattice must be a number or a pair of vectors");
						}
						for(j = 0; j < 2; ++j){
							lua_pushinteger(L, j+1);
							lua_gettable(L, -2);
							if(!lua_isnumber(L, -1)){
								return luaL_argerror(L, 1, "Encountered non-numeric lattice vector");
							}
							Lr[j+i*2] = lua_tonumber(L, -1);
							lua_pop(L, 1);
						}
						lua_pop(L, 1);
					}
				}
				argflag |= 0x1;
			}else if(0 == strcmp("BasisSize", key)){
				lua_Number n;
				int i;
				if(!lua_isnumber(L, -1) || ((n = lua_tonumber(L, -1)), (i = lua_tointeger(L, -1)), ((int)n != i || i <= 0))){
					return luaL_argerror(L, 1, "BasisSize must be a positive integer");
				}
				nG = (unsigned int)i;
				argflag |= 0x2;
			}else{
				return luaL_error(L, "Unrecognized named argument: %s", key);
			}
			lua_pop(L, 1);
		}
		if(0x3 != argflag){
			return luaL_argerror(L, 1, "Must specify Lattice and BasisSize");
		}
	}else{
		// deprecated path
	}
	pS = (S4_Simulation**)lua_newuserdata(L, sizeof(S4_Simulation*));
	luaL_getmetatable(L, S4_simulation_typename);
	lua_setmetatable(L, -2);
	*pS = S4_Simulation_New(Lr, nG, NULL);
	return 1;
}
static int S4L_Simulation__gc(lua_State *L){
	S4_Simulation *S = S4L_get_simulation(L, 1);
	S4_Simulation_Destroy(S);
	return 0;
}
static int S4L_Simulation_Clone(lua_State *L){
	S4_Simulation **pT;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "Clone: 'S4_Simulation' object expected.");

	pT = (S4_Simulation **)lua_newuserdata(L, sizeof(S4_Simulation*));
	luaL_getmetatable(L, S4_simulation_typename);
	lua_setmetatable(L, -2);

	*pT = S4_Simulation_Clone(S);
	return 1;
}

static int S4L_Simulation_SetLattice(lua_State *L){
	int i, j;
	S4_real Lr[4] = { 0, 0, 0, 0 };
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "SetLattice: 'S4_Simulation' object expected.");

	if(lua_isnumber(L, 2)){
		Lr[0] = lua_tonumber(L, 2);
		luaL_argcheck(L, Lr[0] > 0, 2, "1D lattice constant must be positive");
		S4_Simulation_SetLattice(S, Lr);
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
					Lr[2*i+j] = (S4_real)lua_tonumber(L, -1);
				}
				lua_pop(L, 1);
			}
		}
		i = S4_Simulation_SetLattice(S, Lr);
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
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetReciprocalLattice: 'S4_Simulation' object expected.");

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
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "SetNumG: 'S4_Simulation' object expected.");
	n = luaL_checkinteger(L, 2);
	if(n < 1){
		S4L_error(L, "SetNumG: Must have at least 1 G-vector.");
	}
	Simulation_SetNumG(S, n);
	return 0;
}
static int S4L_Simulation_GetNumG(lua_State *L){
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetNumG: 'S4_Simulation' object expected.");
	lua_pushinteger(L, Simulation_GetNumG(S, NULL));
	return 1;
}
static int S4L_Simulation_SetResolution(lua_State *L){
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "SetResolution: 'S4_Simulation' object expected.");
	S->options.resolution = luaL_checkinteger(L, 2);
	if(S->options.resolution < 2){
		S4L_error(L, "SetResolution: Resolution must be at least 2.");
	}
	return 0;
}
static int S4L_Simulation_AddMaterial(lua_State *L){
	int i, j;
	S4_real eps[18];
	int type = 0;
	S4_MaterialID M;
	const char *name;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "AddMaterial: 'S4_Simulation' object expected.");

	name = luaL_checklstring(L, 2, NULL);

	if(2 == lua_rawlen(L, 3)){
		for(j = 0; j < 2; ++j){
			lua_pushinteger(L, 1+j);
			lua_gettable(L, 3);
			if(!lua_isnumber(L, -1)){
				S4L_error(L, "AddMaterial: S4_Material epsilon must be numeric.");
			}else{
				eps[j] = (double)lua_tonumber(L, -1);
			}
			lua_pop(L, 1);
		}
		type = S4_MATERIAL_TYPE_SCALAR_COMPLEX;
	}else if(9 == lua_rawlen(L, 3)){
		for(i = 0; i < 9; ++i){
			lua_pushinteger(L, 1+i);
			lua_gettable(L, 3);
			for(j = 0; j < 2; ++j){
				lua_pushinteger(L, 1+j);
				lua_gettable(L, 4);
				if(!lua_isnumber(L, -1)){
					S4L_error(L, "AddMaterial: S4_Material epsilon must be numeric.");
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
		type = S4_MATERIAL_TYPE_XYTENSOR_COMPLEX;
	}else{
		S4L_error(L, "AddMaterial: Expected either a scalar or tensor value for material %s.", name);
	}

	M = S4_Simulation_SetMaterial(S, -1, name, type, eps);
	if(M < 0){
		S4L_error(L, "AddMaterial: There was a problem allocating the material named '%s'.", name);
		return 0;
	}

	return 0;
}
static int S4L_Simulation_SetMaterial(lua_State *L){
	int i, j;
	double eps[18];
	const char *name;
	S4_MaterialID M;
	int type = 0;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "SetMaterial: 'S4_Simulation' object expected.");

	name = luaL_checklstring(L, 2, NULL);
	M = S4_Simulation_GetMaterialByName(S, name);

	if(2 == lua_rawlen(L, 3)){
		for(j = 0; j < 2; ++j){
			lua_pushinteger(L, 1+j);
			lua_gettable(L, 3);
			if(!lua_isnumber(L, -1)){
				S4L_error(L, "AddMaterial: S4_Material epsilon must be numeric.");
			}else{
				eps[j] = (double)lua_tonumber(L, -1);
			}
			lua_pop(L, 1);
		}
		type = S4_MATERIAL_TYPE_SCALAR_COMPLEX;
	}else if(9 == lua_rawlen(L, 3)){
		for(i = 0; i < 9; ++i){
			lua_pushinteger(L, 1+i);
			lua_gettable(L, 3);
			for(j = 0; j < 2; ++j){
				lua_pushinteger(L, 1+j);
				lua_gettable(L, 4);
				if(!lua_isnumber(L, -1)){
					S4L_error(L, "AddMaterial: S4_Material epsilon must be numeric.");
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
		type = S4_MATERIAL_TYPE_XYTENSOR_COMPLEX;
	}else{
		S4L_error(L, "AddMaterial: Expected either a scalar or tensor value for material %s.", name);
	}

	M = S4_Simulation_SetMaterial(S, M, name, type, eps);
	if(M < 0){
		S4L_error(L, "AddMaterial: There was a problem allocating the material named '%s'.", name);
		return 0;
	}

	return 0;
}
static int S4L_Simulation_AddLayer(lua_State *L){
	S4_LayerID layer;
	const char *name;
	S4_real thickness;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "AddLayer: 'S4_Simulation' object expected.");
	name = luaL_checklstring(L, 2, NULL);
	thickness = luaL_checknumber(L, 3);
	const char *matname = luaL_checklstring(L, 4, NULL);

	S4_MaterialID M = S4_Simulation_GetMaterialByName(S, matname);
	if(M < 0){
		S4L_error(L, "AddLayer: Unknown material '%s'.", matname);
	}
	layer = S4_Simulation_SetLayer(S, -1, name, &thickness, -1, M);
	if(layer < 0){
		S4L_error(L, "AddLayer: There was a problem allocating the layer named '%s'.", name);
		return 0;
	}
	return 0;
}
static int S4L_Simulation_SetLayer(lua_State *L){
	S4_LayerID layer;
	const char *name;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "SetLayer: 'S4_Simulation' object expected.");
	S4_real thickness = luaL_checknumber(L, 3);

	name = luaL_checklstring(L, 2, NULL);
	layer = S4_Simulation_GetLayerByName(S, name);
	if(layer < 0){
		const char *matname = luaL_checklstring(L, 4, NULL);
		S4_MaterialID M = S4_Simulation_GetMaterialByName(S, matname);
		if(M < 0){
			S4L_error(L, "SetLayer: Unknown material '%s'.", matname);
		}
		layer = S4_Simulation_SetLayer(S, -1, name, &thickness, -1, M);
		if(layer < 0){
			S4L_error(L, "SetLayer: There was a problem allocating the layer named '%s'.", name);
			return 0;
		}
	}else{
		S4_Simulation_SetLayer(S, layer, NULL, &thickness, -1, -1);
	}
	return 0;
}
static int S4L_Simulation_SetLayerThickness(lua_State *L){
	S4_LayerID layer;
	const char *name;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "SetLayerThickness: 'S4_Simulation' object expected.");

	name = luaL_checklstring(L, 2, NULL);
	layer = S4_Simulation_GetLayerByName(S, name);
	if(layer < 0){
		S4L_error(L, "SetLayerThickness: S4_Layer named '%s' not found.", name);
	}else{
		S4_real thick = luaL_checknumber(L, 3);
		if(thick < 0){
			S4L_error(L, "SetLayerThickness: Thickness must be non-negative.");
		}
		S4_Simulation_SetLayer(S, layer, NULL, &thick, -1, -1);
	}
	return 0;
}

static int S4L_Simulation_AddLayerCopy(lua_State *L){
	S4_LayerID layer;
	const char *name;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "AddLayerCopy: 'S4_Simulation' object expected.");
	name = luaL_checklstring(L, 2, NULL);
	S4_real thickness = luaL_checknumber(L, 3);
	const char *copyname = luaL_checklstring(L, 4, NULL);
	S4_LayerID Lcopy = S4_Simulation_GetLayerByName(S, copyname);
	if(Lcopy < 0){
		S4L_error(L, "AddLayerCopy: Layer not found: '%s'.", copyname);
	}
	layer = S4_Simulation_SetLayer(S, -1, name, &thickness, Lcopy, -1);

	if(layer < 0){
		S4L_error(L, "AddLayerCopy: There was a problem allocating the layer named '%s'.", name);
		return 0;
	}
	return 0;
}

/* Expected stack:
 *   1 S4_Simulation
 *   2 layer name string
 *   3 circle material string
 *   4 circle center (numeric table of length 2)
 *   5 circle radius
 */
static int S4L_Simulation_SetLayerPatternCircle(lua_State *L){
	int j, ret;
	const char *layer_name;
	const char *material_name;
	S4_LayerID layer;
	S4_MaterialID M;
	double center[2];
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "SetLayerPatternCircle: 'Simulation' object expected.");
	S4_real radius = luaL_checknumber(L, 5);

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = S4_Simulation_GetLayerByName(S, layer_name);
	if(layer < 0){
		S4L_error(L, "SetLayerPatternCircle: Layer named '%s' not found.", layer_name);
		return 0;
	}
	if(S4_Layer_IsCopy(S, layer) > 0){
		S4L_error(L, "SetLayerPatternCircle: Cannot pattern a layer copy.");
		return 0;
	}
	material_name = luaL_checklstring(L, 3, NULL);
	M = S4_Simulation_GetMaterialByName(S, material_name);
	if(M < 0){
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
	S4_real hw[2] = { radius, radius };
	ret = S4_Layer_SetRegionHalfwidths(
		S, layer, M, S4_REGION_TYPE_ELLIPSE, hw, center, NULL
	);
	if(0 != ret){
		S4L_error(L, "SetLayerPatternCircle: There was a problem allocating the pattern.");
		return 0;
	}
	return 0;
}

/* Expected stack:
 *   1 S4_Simulation
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
	S4_LayerID layer;
	S4_MaterialID M;
	double center[2], halfwidths[2];
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "SetLayerPatternEllipse: 'S4_Simulation' object expected");

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = S4_Simulation_GetLayerByName(S, layer_name);
	if(layer < 0){
		S4L_error(L, "SetLayerPatternEllipse: S4_Layer named '%s' not found.", layer_name);
		return 0;
	}
	if(S4_Layer_IsCopy(S, layer) > 0){
		S4L_error(L, "SetLayerPatternEllipse: Cannot pattern a layer copy.");
		return 0;
	}
	material_name = luaL_checklstring(L, 3, NULL);
	M = S4_Simulation_GetMaterialByName(S, material_name);
	if(M < 0){
		S4L_error(L, "SetLayerPatternEllipse: S4_Material named '%s' not found.", material_name);
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
	S4_real angle = luaL_checknumber(L, 5) / 360.;
	ret = S4_Layer_SetRegionHalfwidths(
		S, layer, M, S4_REGION_TYPE_ELLIPSE, halfwidths, center, &angle
	);
	if(0 != ret){
		S4L_error(L, "SetLayerPatternEllipse: There was a problem allocating the pattern.");
		return 0;
	}
	return 0;
}

/* Expected stack:
 *   1 S4_Simulation
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
	S4_LayerID layer;
	S4_MaterialID M;
	double center[2], halfwidths[2];
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "SetLayerPatternRectangle: 'S4_Simulation' object expected.");

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = S4_Simulation_GetLayerByName(S, layer_name);
	if(layer < 0){
		S4L_error(L, "SetLayerPatternRectangle: S4_Layer named '%s' not found.", layer_name);
		return 0;
	}
	if(S4_Layer_IsCopy(S, layer) > 0){
		S4L_error(L, "SetLayerPatternRectangle: Cannot pattern a layer copy.");
		return 0;
	}
	material_name = luaL_checklstring(L, 3, NULL);

	M = S4_Simulation_GetMaterialByName(S, material_name);
	if(M < 0){
		S4L_error(L, "SetLayerPatternRectangle: S4_Material named '%s' not found.", material_name);
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
	S4_real angle = luaL_checknumber(L, 5) / 360.;
	ret = S4_Layer_SetRegionHalfwidths(
		S, layer, M, S4_REGION_TYPE_RECTANGLE, halfwidths, center, &angle
	);
	if(0 != ret){
		S4L_error(L, "SetLayerPatternRectangle: There was a problem allocating the pattern.");
		return 0;
	}
	return 0;
}

/* Expected stack:
 *   1 S4_Simulation
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
	S4_LayerID layer;
	S4_MaterialID M;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "SetLayerPatternPolygon: 'S4_Simulation' object expected.");

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = S4_Simulation_GetLayerByName(S, layer_name);
	if(layer < 0){
		S4L_error(L, "SetLayerPatternPolygon: S4_Layer named '%s' not found.", layer_name);
		return 0;
	}
	if(S4_Layer_IsCopy(S, layer) > 0){
		S4L_error(L, "SetLayerPatternPolygon: Cannot pattern a layer copy.");
		return 0;
	}
	material_name = luaL_checklstring(L, 3, NULL);
	M = S4_Simulation_GetMaterialByName(S, material_name);
	if(M < 0){
		S4L_error(L, "SetLayerPatternPolygon: S4_Material named '%s' not found.", material_name);
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
	S4_real angle = luaL_checknumber(L, 5) / 360.;
	ret = S4_Layer_SetRegionVertices(
		S, layer, M, S4_REGION_TYPE_POLYGON,
		nvert, vert, center, &angle
	);
	free(vert);
	if(0 != ret){
		S4L_error(L, "SetLayerPatternPolygon: There was a problem allocating the pattern.");
		return 0;
	}
	return 0;
}


/* Expected stack:
 *   1 S4_Simulation
 *   2 table
 * First table:
 *   Each entry a table of three entries, first is integer, second is 'x' or 'y', third is table of two entries re and im
 *   If integer is negative, indicates a backward propagating mode.
 */
static int S4L_Simulation_SetExcitationExterior(lua_State *L){
	S4_Simulation *S = S4L_get_simulation(L, 1);
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
	S4_Simulation_ExcitationExterior(S, n, exg, ex);
	return 0;
}

static int S4L_Simulation_SetExcitationInterior(lua_State *L){
	//S4_Simulation *S = S4L_get_simulation(L, 1);
	S4L_error(L, "SetExcitationInterior: Not implemented.");
	return 0;
}

/* Expected stack:
 *   1 S4_Simulation
 *   2 layer name string
 *   3 k (numeric table of length 2)
 *   4 position in plane (numeric table of length 2)
 *   5 dipole moment (numeric table of length 3)
 *   6 amplitude and phase (numeric table of length 2; phase in degrees)
 */
static int S4L_Simulation_SetExcitationDipole(lua_State *L){
	S4_Simulation *S = S4L_get_simulation(L, 1);
	const char *layer;
	int j, ret;
	double k[2], pos[2], moment[6], ampphase[2];
	luaL_argcheck(L, S != NULL, 1, "SetExcitationDipole: 'S4_Simulation' object expected.");
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
	ret = S4_Simulation_ExcitationDipole(S, k, layer, pos, moment);
	if(0 != ret){
		HandleSolutionErrorCode(L, "SetExcitationDipole", ret);
		return 0;
	}

	return 0;
}

/* Expected stack:
 *   1 S4_Simulation
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
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "SetExcitationPlanewave: 'S4_Simulation' object expected.");

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

	order = luaL_optinteger(L, 5, 0);
	if(order > 0){ order--; }
	ret = Simulation_MakeExcitationPlanewave(S, angle, pol_s, pol_p, order);
	if(0 != ret){
		HandleSolutionErrorCode(L, "SetExcitationPlanewave", ret);
		return 0;
	}

	return 0;
}
static int S4L_Simulation_SetFrequency(lua_State *L){
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "SetFrequency: 'S4_Simulation' object expected");
	S4_real freq[2] = {
		luaL_checknumber(L, 2),
		luaL_optnumber(L, 3, 0)
	};
	if(freq[0] <= 0){
		S4L_error(L, "SetFrequency: Frequency must be positive.");
		return 0;
	}
	if(freq[1] > 0){
		S4L_error(L, "SetFrequency: Imaginary component of frequency must be negative.");
		return 0;
	}
	S4_Simulation_SetFrequency(S, freq);
	return 0;
}


static int S4L_Simulation_GetGList(lua_State *L){
	int *G;
	int n, i, j;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetGList: 'S4_Simulation' object expected.");

	n = S4_Simulation_GetBases(S, NULL);
	G = (int*)malloc(sizeof(int)*2*n);
	if(NULL == G){
		return 0;
	}
	S4_Simulation_GetBases(S, G);

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
	free(G);
	return 1;
}


static int S4L_Simulation_GetDiffractionOrder(lua_State *L){
	int *G;
	int u, v = 0, n, i;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetDiffractionOrder: 'S4_Simulation' object expected.");

	u = luaL_checkinteger(L, 2);
	if(lua_gettop(L) > 2){
		v = luaL_checkinteger(L, 3);
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
 *   1 S4_Simulation
 *   2 layer name string
 *   3 layer z-offset
 */
static int S4L_Simulation_GetPoyntingFlux(lua_State *L){
	S4_real power[4];
	int ret;
	const char *layer_name;
	S4_LayerID layer;
	S4_real offset = 0;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetPoyntingFlux: 'S4_Simulation' object expected.");

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = S4_Simulation_GetLayerByName(S, layer_name);
	if(layer < 0){
		S4L_error(L, "GetPoyntingFlux: S4_Layer named '%s' not found.", layer_name);
		return 0;
	}
	offset = luaL_optnumber(L, 3, 0);
	ret = S4_Simulation_GetPowerFlux(S,
		layer,
		&offset,
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
	S4_LayerID layer;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetPoyntingFluxByOrder: 'S4_Simulation' object expected.");

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = S4_Simulation_GetLayerByName(S, layer_name);
	if(layer < 0){
		S4L_error(L, "GetPoyntingFluxByOrder: S4_Layer named '%s' not found.", layer_name);
		return 0;
	}

	ret = S4_Simulation_SolveLayer(S, layer);
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
		&S->layer[layer],
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
	S4_LayerID layer;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetAmplitudes: 'S4_Simulation' object expected.");

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = S4_Simulation_GetLayerByName(S, layer_name);
	if(NULL == layer){
		S4L_error(L, "GetAmplitudes: S4_Layer named '%s' not found.", layer_name);
		return 0;
	}

	ret = S4_Simulation_SolveLayer(S, layer);
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

static int S4L_Simulation_GetWaves(lua_State *L){
	double *waves;
	int n, n2, i, k, ret;
	const char *layer_name;
	S4_LayerID layer;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetWaves: 'S4_Simulation' object expected.");

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = S4_Simulation_GetLayerByName(S, layer_name);
	if(layer < 0){
		S4L_error(L, "GetWaves: S4_Layer named '%s' not found.", layer_name);
		return 0;
	}

	ret = S4_Simulation_SolveLayer(S, layer);
	if(0 != ret){
		HandleSolutionErrorCode(L, "GetWaves", ret);
		return 0;
	}

	n = S4_Simulation_GetBases(S, NULL);
	n2 = 2*n;

	waves = (double*)malloc(sizeof(double)*11*n2);
	S4_Simulation_GetWaves(S, layer, waves);

	lua_createtable(L, n2, 0);
	for(i = 0; i < n2; ++i){
		const double *wave = &waves[11*i];
		lua_pushinteger(L, i+1);
		/* a wave object is:
		 *   direction = {kx,ky,kzr,kzi}
		 *   polarization = {x,y,z}
		 *   cu = {re,im}
		 *   cv = {re,im}
		 */
		lua_createtable(L, 0, 4);
		{
			lua_createtable(L, 4, 0);
			for(k = 0; k < 4; ++k){
				lua_pushinteger(L, k+1);
				lua_pushnumber(L, wave[k]);
				lua_settable(L, -3);
			}
			lua_setfield(L, -2, "k");

			lua_createtable(L, 3, 0);
			for(k = 0; k < 3; ++k){
				lua_pushinteger(L, k+1);
				lua_pushnumber(L, wave[4+k]);
				lua_settable(L, -3);
			}
			lua_setfield(L, -2, "u");

			lua_createtable(L, 2, 0);
			for(k = 0; k < 2; ++k){
				lua_pushinteger(L, k+1);
				lua_pushnumber(L, wave[7+k]);
				lua_settable(L, -3);
			}
			lua_setfield(L, -2, "cu");

			lua_createtable(L, 2, 0);
			for(k = 0; k < 2; ++k){
				lua_pushinteger(L, k+1);
				lua_pushnumber(L, wave[9+k]);
				lua_settable(L, -3);
			}
			lua_setfield(L, -2, "cv");
		}
		lua_settable(L, -3);
	}

	free(waves);
	return 1;
}


/*
static int S4L_Simulation_EnableBasisFieldDump(lua_State *L){
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "EnableBasisFieldDump: 'S4_Simulation' object expected.");

	S->options.do_vector_field_dump = 1;
	return 0;
}
*/
static int S4L_Simulation_UseDiscretizedEpsilon(lua_State *L){
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "UseDiscretizedEpsilon: 'S4_Simulation' object expected.");

	if(lua_gettop(L) > 1){
		S->options.use_discretized_epsilon = lua_toboolean(L, 2);
	}else{
		S->options.use_discretized_epsilon = 1;
	}
	return 0;
}
static int S4L_Simulation_UseSubpixelSmoothing(lua_State *L){
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "UseSubpixelSmoothing: 'S4_Simulation' object expected.");

	if(lua_gettop(L) > 1){
		S->options.use_subpixel_smoothing = lua_toboolean(L, 2);
	}else{
		S->options.use_subpixel_smoothing = 1;
	}
	return 0;
}
static int S4L_Simulation_UseLanczosSmoothing(lua_State *L){
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "UseLanczosSmoothing: 'S4_Simulation' object expected.");

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
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "S4L_Simulation_SetLanczosSmoothingWidth: 'S4_Simulation' object expected.");
	width = luaL_checknumber(L, 2);
	if(lua_gettop(L) > 2){
		power = luaL_checkinteger(L, 3);
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
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "UsePolarizationDecomposition: 'S4_Simulation' object expected.");

	if(lua_gettop(L) > 1){
		S->options.use_polarization_basis = lua_toboolean(L, 2);
	}else{
		S->options.use_polarization_basis = 1;
	}
	return 0;
}
static int S4L_Simulation_UseJonesVectorBasis(lua_State *L){
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "UseJonesVectorBasis: 'S4_Simulation' object expected.");

	if(lua_gettop(L) > 1){
		S->options.use_jones_vector_basis = lua_toboolean(L, 2);
	}else{
		S->options.use_jones_vector_basis = 1;
	}
	return 0;
}
static int S4L_Simulation_UseNormalVectorBasis(lua_State *L){
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "UseNormalVectorBasis: 'S4_Simulation' object expected.");

	if(lua_gettop(L) > 1){
		S->options.use_normal_vector_basis = lua_toboolean(L, 2);
	}else{
		S->options.use_normal_vector_basis = 1;
	}
	return 0;
}
static int S4L_Simulation_UseExperimentalFMM(lua_State *L){
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "UseExperimentalFMM: 'S4_Simulation' object expected.");

	if(lua_gettop(L) > 1){
		S->options.use_experimental_fmm = lua_toboolean(L, 2);
	}else{
		S->options.use_experimental_fmm = 1;
	}
	return 0;
}

static int S4L_Simulation_SetBasisFieldDumpPrefix(lua_State *L){
	const char *prefix = NULL;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "SetBasisFieldDumpPrefix: 'S4_Simulation' object expected.");

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
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "SetLatticeTruncation: 'S4_Simulation' object expected.");
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
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "SetVerbosity: 'S4_Simulation' object expected.");
	verbosity = luaL_checkinteger(L, 2);

	if(verbosity < 0){ verbosity = 0; }
	if(verbosity > 9){ verbosity = 9; }
	S->options.verbosity = verbosity;

	return 0;
}

static int S4L_Simulation_UseLessMemory(lua_State *L){
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "UseLessMemory: 'S4_Simulation' object expected.");

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
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "OutputStructurePOVRay: 'S4_Simulation' object expected.");
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
	S4_Layer *layer;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "OutputLayerPatternDescription: 'S4_Simulation' object expected.");

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "OutputLayerPatternDescription: S4_Layer named '%s' not found.", layer_name);
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
	S4_Layer *layer;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "OutputLayerPatternRealization: 'S4_Simulation' object expected.");

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "OutputLayerPatternRealization: S4_Layer named '%s' not found.", layer_name);
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
			luaL_checkinteger(L, 3),
			luaL_checkinteger(L, 4),
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
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetEField: 'S4_Simulation' object expected.");

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
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetHField: 'S4_Simulation' object expected.");

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
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetFields: 'S4_Simulation' object expected.");

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

/* S:GetFieldPlane(z, {nu, nv}, format, filename)
 *   z: number specify global z coordinates
 *   format: string specifying return format
 *     'Array': returns a pair (E and H) of 2D arrays of 3D complex vectors
 *     'FileWrite': dumps to filename.E or filename.H (x, y, components)
 *     'FileAppend': appends to filename.E or filename.H (x, y, z, components)
 */
static int S4L_Simulation_GetFieldPlane(lua_State *L){
	int i, j, ret;
	int nxy[2];
	double z;
	double *Efields, *Hfields;
	const char *fmt;
	const char *fbasename; size_t len;
	char *filename;
	FILE *fp;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetFieldPlane: 'S4_Simulation' object expected.");

	z = luaL_checknumber(L, 2);
	luaL_argcheck(L, lua_istable(L, 3) && (lua_rawlen(L, 3) == 2), 3, "GetFieldPlane: pair of grid sample counts expected.");
	for(i = 0; i < 2; ++i){
		lua_rawgeti(L, 3, i+1);
		if(!lua_isnumber(L, -1)){
			S4L_error(L, "GetFieldPlane: pair of grid sample counts expected.");
		}
		nxy[i] = (int)lua_tointeger(L, -1);
		if(nxy[i] <= 0){
			S4L_error(L, "GetFieldPlane: grid sample counts must be positive.");
		}
		lua_pop(L, 1);
	}
	fmt = luaL_checkstring(L, 4);
	len = 5;
	fbasename = luaL_optlstring(L, 5, "field", &len);
	filename = (char*)malloc(sizeof(char) * (len+3));
	strcpy(filename, fbasename);
	filename[len+0] = '.';
	filename[len+2] = '\0';

	Efields = (double*)malloc(sizeof(double) * 2*3 * nxy[0] * nxy[1]);
	Hfields = (double*)malloc(sizeof(double) * 2*3 * nxy[0] * nxy[1]);

	ret = Simulation_GetFieldPlane(S, nxy, z, Efields, Hfields);
	if(0 != ret){
		HandleSolutionErrorCode(L, "GetFieldPlane", ret);
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
	}else{ /* Array */
		unsigned k, i3, i2;
		double *F[2] = { Efields, Hfields };

		for(k = 0; k < 2; ++k){
			lua_createtable(L, nxy[0], 0);
			for(i = 0; i < nxy[0]; ++i){
				lua_createtable(L, nxy[1], 0);
				for(j = 0; j < nxy[1]; ++j){
					lua_createtable(L, 3, 0);
					for(i3 = 0; i3 < 3; ++i3){
						lua_createtable(L, 2, 0);
						for(i2 = 0; i2 < 2; ++i2){
							lua_pushnumber(L, F[k][2*(3*(i+j*nxy[0])+i3)+i2]);
							lua_rawseti(L, -2, i2+1);
						}
						lua_rawseti(L, -2, i3+1);
					}
					lua_rawseti(L, -2, j+1);
				}
				lua_rawseti(L, -2, i+1);
			}
		}

		ret = 2;
	}

	free(Hfields);
	free(Efields);
	free(filename);
	return ret;
}

static int S4L_Simulation_GetSMatrixDeterminant(lua_State *L){
	int ret;
	double mant[2], base;
	int expo;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetSMatrixDeterminant: 'S4_Simulation' object expected.");

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
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetEpsilon: 'S4_Simulation' object expected.");

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
 *   1 S4_Simulation
 *   2 layer name string
 *   3 layer offset
 */
static int S4L_Simulation_GetStressTensorIntegral(lua_State *L){
	int ret;
	double Tint[6];
	const char *layer_name;
	S4_Layer *layer;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetStressTensorIntegral: 'S4_Simulation' object expected.");

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "GetStressTensorIntegral: S4_Layer named '%s' not found.", layer_name);
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
 *   1 S4_Simulation
 *   2 layer name string
 */
static int S4L_Simulation_GetLayerVolumeIntegral(lua_State *L, char which, const char *name){
	int ret;
	double integral[2];
	const char *layer_name;
	S4_Layer *layer;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "Get*LayerIntegral: 'S4_Simulation' object expected.");

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "%s: S4_Layer named '%s' not found.", name, layer_name);
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
	S4_Layer *layer;
	S4_Simulation *S = S4L_get_simulation(L, 1);
	luaL_argcheck(L, S != NULL, 1, "GetLayerZIntegral: 'S4_Simulation' object expected.");

	layer_name = luaL_checklstring(L, 2, NULL);
	layer = Simulation_GetLayerByName(S, layer_name, NULL);
	if(NULL == layer){
		S4L_error(L, "GetLayerZIntegral: S4_Layer named '%s' not found.", layer_name);
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

#ifdef _WIN32
extern __declspec(dllexport)
#else
LUALIB_API
#endif
int luaopen_RCWA(lua_State *L){
	static const struct luaL_Reg S4_lib[] = {
		{"NewSimulation", S4L_NewSimulation},
		{"NewSpectrumSampler", S4L_NewSpectrumSampler},
		{"NewInterpolator", S4L_NewInterpolator},
		{"SolveInParallel", S4L_SolveInParallel},
		{"ConvertUnits", S4L_ConvertUnits},
		{"Integrate", S4L_Integrate},
		{"ParallelInvoke", S4L_ParallelInvoke},
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
		{"GetWaves", S4L_Simulation_GetWaves},
		{"GetStressTensorIntegral", S4L_Simulation_GetStressTensorIntegral},
		{"GetLayerEnergyDensityIntegral", S4L_Simulation_GetLayerEnergyDensityIntegral},
		{"GetLayerElectricEnergyDensityIntegral", S4L_Simulation_GetLayerElectricEnergyDensityIntegral},
		{"GetLayerMagneticEnergyDensityIntegral", S4L_Simulation_GetLayerMagneticEnergyDensityIntegral},
		{"GetLayerElectricFieldIntensityIntegral", S4L_Simulation_GetLayerElectricFieldIntensityIntegral},
		{"GetLayerZIntegral", S4L_Simulation_GetLayerZIntegral},
		{"GetEField", S4L_Simulation_GetEField},
		{"GetHField", S4L_Simulation_GetHField},
		{"GetFields", S4L_Simulation_GetFields},
		{"GetFieldPlane", S4L_Simulation_GetFieldPlane},
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

	luaL_newmetatable(L, S4_simulation_typename);
	luaL_newlib(L, SimulationObj);
	lua_setfield(L, -2, "__index");
	lua_pushstring(L, "__gc");
	lua_pushcfunction(L, &S4L_Simulation__gc);
	lua_settable(L, -3);
	lua_pop(L, 1);
	/*
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
	*/
	return 1;
}

lua_State *new_S4_lua_state(){
	lua_State *L = luaL_newstate(); /* opens Lua */

#if LUA_VERSION_NUM < 502
	luaopen_RCWA(L);
	lua_setglobal(L, "S4");
#else
	luaL_requiref(L, "S4", &luaopen_RCWA, 1);
	lua_pop(L, 1);
#endif
	luaL_openlibs(L); /* opens the standard libraries */
	return L;
}

#ifndef LUA_BUILD_AS_DLL

void threadsafe_init(){
	extern void exactinit();
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



#include <signal.h>
#if LUA_VERSION_NUM < 502
static int traceback (lua_State *L) {
  if (!lua_isstring(L, 1))  /* 'message' not a string? */
    return 1;  /* keep it intact */
  lua_getfield(L, LUA_GLOBALSINDEX, "debug");
  if (!lua_istable(L, -1)) {
    lua_pop(L, 1);
    return 1;
  }
  lua_getfield(L, -1, "traceback");
  if (!lua_isfunction(L, -1)) {
    lua_pop(L, 2);
    return 1;
  }
  lua_pushvalue(L, 1);  /* pass error message */
  lua_pushinteger(L, 2);  /* skip this function and traceback */
  lua_call(L, 2, 1);  /* call debug.traceback */
  return 1;
}
#else
static int traceback (lua_State *L) {
  const char *msg = lua_tostring(L, 1);
  if (msg)
    luaL_traceback(L, L, msg, 1);
  else if (!lua_isnoneornil(L, 1)) {  /* is there an error object? */
    if (!luaL_callmeta(L, 1, "__tostring"))  /* try its 'tostring' metamethod */
      lua_pushliteral(L, "(no error message)");
  }
  return 1;
}
#endif
static lua_State *globalL = NULL;
static void lstop (lua_State *L, lua_Debug *ar) {
  (void)ar;  /* unused arg. */
  lua_sethook(L, NULL, 0, 0);
  luaL_error(L, "interrupted!");
}
static void laction (int i) {
  signal(i, SIG_DFL); /* if another SIGINT happens before lstop,
                              terminate process (default action) */
  lua_sethook(globalL, lstop, LUA_MASKCALL | LUA_MASKRET | LUA_MASKCOUNT, 1);
}
static int docall (lua_State *L, int narg, int nres) {
  int status;
  int base = lua_gettop(L) - narg;  /* function index */
  lua_pushcfunction(L, traceback);  /* push traceback function */
  lua_insert(L, base);  /* put it under chunk and args */
  globalL = L;  /* to be available to 'laction' */
  signal(SIGINT, laction);
  status = lua_pcall(L, narg, nres, base);
  signal(SIGINT, SIG_DFL);
  lua_remove(L, base);  /* remove traceback function */
  return status;
}

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

	L = new_S4_lua_state(); /* opens Lua */

	threadsafe_init();
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
				error = docall(L, 0, LUA_MULTRET);
				if(error){
					fprintf(stderr, "%s\n", lua_tostring(L, -1));
					lua_pop(L, 1); /* pop error message from the stack */
				}
			}
		}
	}else{ /* run in REPL mode */
		fprintf(stdout, "No input file given, running in interactive mode\n"); fflush(stdout);

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

#endif // LUA_BUILD_AS_DLL
