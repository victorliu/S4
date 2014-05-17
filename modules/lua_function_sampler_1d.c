/*
sampler = FunctionSampler1D.New{
	MaxCurvature = 10,
	YToleranceAbsolute = 0,
	YToleranceRelative = 1e-3,
	XTolerance = 1e-6,
	YBias = 'min', -- or 'max' or 'none' or ''
}

sampler:Add(x, y, id) -- id is optional
sampler:IsDone()
xvals = sampler:GetNext()
sampler:SetParameters{
	MaxCurvature = 10
}


sampler = FunctionSampler2D.New{
	MaxGaussianCurvature = 360,
	MaxPrincipalCurvature = 10,
	ZToleranceAbsolute = 0,
	ZToleranceRelative = 1e-3,
	XYTolerance = 1e-5,
}

*/
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "function_sampler_1d.h"

#include <stdio.h>

const static char *sampler_typename = "FunctionSampler1D";

typedef struct{
	function_sampler_1d sampler;
	int nx;
	double *x;
} sampler_obj;

static sampler_obj* getobj(lua_State *L, int i){
	return luaL_checkudata(L, i, sampler_typename);
}

static int FunctionSampler1D_new(lua_State *L){
	function_sampler_1d_options options;
	sampler_obj* obj;

	function_sampler_1d_options_defaults(&options);
	obj = (sampler_obj*)lua_newuserdata(L, sizeof(sampler_obj));
	luaL_setmetatable(L, sampler_typename);

	obj->sampler = function_sampler_1d_new(&options);
	obj->nx = 0;
	obj->x = NULL;
	return 1;
}
static int sampler__gc(lua_State *L){
	sampler_obj *obj = getobj(L, 1);
	function_sampler_1d_destroy(obj->sampler);
	free(obj->x);
	return 0;
}
static int sampler_add(lua_State *L){
	double x, y;
	int id;
	sampler_obj *obj = getobj(L, 1);

	x = luaL_checknumber(L, 2);
	y = luaL_checknumber(L, 3);
	id = luaL_optinteger(L, 4, 0);

	function_sampler_1d_add(obj->sampler, x, y, id);

	return 0;
}
static int sampler_clear(lua_State *L){
	sampler_obj *obj = getobj(L, 1);
	function_sampler_1d_clear(obj->sampler);
	return 0;
}
static int sampler_is_done(lua_State *L){
	sampler_obj *obj = getobj(L, 1);
	lua_pushboolean(L, function_sampler_1d_is_done(obj->sampler));
	return 1;
}
static int sampler_get_next(lua_State *L){
	int num_req, i;
	sampler_obj *obj = getobj(L, 1);

	num_req = luaL_optinteger(L, 2, 0);
	if(0 == num_req){
		num_req = function_sampler_1d_num_refine(obj->sampler);
	}
	if(num_req > obj->nx){
		obj->nx = num_req;
		obj->x = (double*)realloc(obj->x, sizeof(double) * obj->nx);
	}
	num_req = function_sampler_1d_get_refine(obj->sampler, num_req, obj->x);
	lua_createtable(L, num_req, 0);
	for(i = 0; i < num_req; ++i){
		lua_pushnumber(L, obj->x[i]);
		lua_rawseti(L, -2, i+1);
	}

	return 1;
}

static double get_number(lua_State *L, double min, double max, int iarg, const char *arg){
	double val;
	if(!lua_isnumber(L, -1)){
		luaL_error(L, "`%s' must be a number", arg);
	}
	val = lua_tonumber(L, -1);
	if(val < min || val > max){
		luaL_error(L, "`%s' is out of range", arg);
	}
	return val;
}

static int sampler_set_params(lua_State *L){
	function_sampler_1d_options *opts;
	sampler_obj *obj = getobj(L, 1);
	luaL_argcheck(L, lua_istable(L, 2), 2, "Expected table of parameters");

	opts = function_sampler_1d_get_options(obj->sampler);

	lua_pushnil(L);
	while(lua_next(L, 2)){
		const char *key;
		if(!lua_isstring(L, -2)){
			return luaL_argerror(L, 2, "Expected string keys");
		}
		key = lua_tostring(L, -2);
		if(strcmp("MaxCurvature", key) == 0){
			double val = get_number(L, 0, 360, 1, "MaxCurvature");
			opts->max_curvature = sin((M_PI/180.) * val);
		}else if(strcmp("XTolerance", key) == 0){
			double val = get_number(L, 0, DBL_MAX, 1, "XTolerance");
			opts->min_dx = val;
		}else if(strcmp("YToleranceRelative", key) == 0){
			double val = get_number(L, 0, DBL_MAX, 1, "YToleranceRelative");
			opts->min_dy_rel = val;
		}else if(strcmp("YToleranceAbsolute", key) == 0){
			double val = get_number(L, 0, DBL_MAX, 1, "YToleranceAbsolute");
			opts->min_dy_abs = val;
		}else if(strcmp("YBias", key) == 0){
			if(lua_isnil(L, -1)){ opts->range_bias = 0; }
			else{
				const char *str = lua_tostring(L, -1);
				if(!lua_isstring(L, -1)){
					luaL_error(L, "YBias must be one of `min', `max', `', `none', or nil");
				}
				if(strcasecmp("min", str) == 0){
					opts->range_bias = 1;
				}else if(strcasecmp("max", str) == 0){
					opts->range_bias = 2;
				}else if(strcasecmp("", str) == 0 || strcasecmp("none", str) == 0){
					opts->range_bias = 0;
				}
			}
		}
		lua_pop(L, 1);
	}
	return 0;
}

static int sampler_get_samples(lua_State *L){
	int i, n;
	sampler_obj *obj = getobj(L, 1);

	n = function_sampler_1d_num_samples(obj->sampler);

	lua_createtable(L, n, 0);
	for(i = 0; i < n; ++i){
		double x, y;
		int id;
		function_sampler_1d_get(obj->sampler, i, &x, &y, &id);
		lua_createtable(L, 3, 0);
		lua_pushnumber(L, x);
		lua_rawseti(L, -2, 1);
		lua_pushnumber(L, y);
		lua_rawseti(L, -2, 2);
		lua_pushinteger(L, id);
		lua_rawseti(L, -2, 3);

		lua_rawseti(L, -2, i+1);
	}

	return 1;
}
static int sampler__len(lua_State *L){
	sampler_obj *obj = getobj(L, 1);
	
	lua_pushinteger(L, function_sampler_1d_num_samples(obj->sampler));
	
	return 1;
}

#ifdef WIN32
__declspec(dllexport)
#endif
int luaopen_FunctionSampler1D(lua_State *L){
	static const luaL_Reg sampler_lib[] = {
		{ "Add", &sampler_add },
		{ "Clear", &sampler_clear },
		{ "IsDone", &sampler_is_done },
		{ "GetNext", &sampler_get_next },
		{ "SetParameters", &sampler_set_params },
		{ "GetSamples", &sampler_get_samples },
		{ NULL, NULL }
	};

	static const luaL_Reg FunctionSampler1D_lib[] = {
		{ "New", &FunctionSampler1D_new },
		{ NULL, NULL }
	};

	luaL_newlib(L, FunctionSampler1D_lib);

	// Stack: MyLib
	luaL_newmetatable(L, sampler_typename); // Stack: MyLib meta
	luaL_newlib(L, sampler_lib);
	lua_setfield(L, -2, "__index"); // Stack: MyLib meta

	lua_pushstring(L, "__gc");
	lua_pushcfunction(L, sampler__gc); // Stack: MyLib meta "__gc" fptr
	lua_settable(L, -3); // Stack: MyLib meta
	
	lua_pushstring(L, "__len");
	lua_pushcfunction(L, sampler__len); // Stack: MyLib meta "__gc" fptr
	lua_settable(L, -3); // Stack: MyLib meta
	
	lua_pop(L, 1); // Stack: MyLib

	return 1;
}



