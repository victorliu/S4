#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <stdlib.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include "S4.h"

#define LUA_S4_SIMULATION_TYPENAME "S4v2.Simulation"
#define LUA_S4_MATERIAL_TYPENAME   "S4v2.Material"
#define LUA_S4_LAYER_TYPENAME      "S4v2.Layer"

typedef struct{
	S4_Simulation *S;
	void *userdata;
	int externally_owned; // nonzero if we should not call __gc
} lua_S4_Simulation;

typedef struct{
	S4_Simulation *S;
	S4_LayerID L;
	int Sref;
} lua_S4_Layer;

typedef struct{
	S4_Simulation *S;
	S4_MaterialID M;
	int Sref;
} lua_S4_Material;

int lua_S4_len(lua_State *L, int arg);
int lua_S4_isinteger(lua_State *L, int arg);
int lua_S4_getcomplex(lua_State *L, int arg, S4_real *v);
int lua_S4_getvec3(lua_State *L, int arg, S4_real *v);
int lua_S4_gettensor(lua_State *L, int arg, S4_real *eps);

/* Forward declarations for object metamethods and helpers */
lua_S4_Material* lua_S4_Material_check(lua_State *L, int narg);
lua_S4_Layer* lua_S4_Layer_check(lua_State *L, int narg);
S4_Simulation* lua_S4_Simulation_check(lua_State *L, int narg);
int lua_S4_Simulation_push(lua_State *L, S4_Simulation *S, void *userdata, int ext_owned);
int lua_S4_Material_push(lua_State *L, S4_Simulation *S, S4_MaterialID M);

int lua_S4_Material__gc(lua_State *L);
int lua_S4_Layer__gc(lua_State *L);
int lua_S4_Simulation__gc(lua_State *L);

/* Forward declarations for object methods */
int lua_S4_NewSimulation(lua_State *L);

int lua_S4_Material_GetName(lua_State *L);
int lua_S4_Material_SetName(lua_State *L);
int lua_S4_Material_GetEpsilon(lua_State *L);
int lua_S4_Material_SetEpsilon(lua_State *L);

int lua_S4_Layer_GetName(lua_State *L);
int lua_S4_Layer_SetName(lua_State *L);
int lua_S4_Layer_GetThickness(lua_State *L);
int lua_S4_Layer_SetThickness(lua_State *L);
int lua_S4_Layer_ClearRegions(lua_State *L);
int lua_S4_Layer_SetRegion(lua_State *L);

int lua_S4_Simulation_Clone(lua_State *L);
int lua_S4_Simulation_SetLattice(lua_State *L);
int lua_S4_Simulation_GetLattice(lua_State *L);
int lua_S4_Simulation_SetBases(lua_State *L);
int lua_S4_Simulation_GetBases(lua_State *L);
int lua_S4_Simulation_SetFrequency(lua_State *L);
int lua_S4_Simulation_GetFrequency(lua_State *L);
int lua_S4_Simulation_AddMaterial(lua_State *L);
int lua_S4_Simulation_AddLayer(lua_State *L);
int lua_S4_Simulation_GetMaterial(lua_State *L);
int lua_S4_Simulation_GetLayer(lua_State *L);
int lua_S4_Simulation_ExcitationPlanewave(lua_State *L);
int lua_S4_Simulation_GetEpsilon(lua_State *L);
int lua_S4_Simulation_GetFields(lua_State *L);
int lua_S4_Layer_GetPowerFlux(lua_State *L);
//int lua_S4_Layer_GetPowerFluxes(lua_State *L);
int lua_S4_Layer_GetWaves(lua_State *L);

#ifdef _WIN32
extern __declspec(dllexport)
#endif
int luaopen_S4v2(lua_State *L){
	static const struct luaL_Reg S4v2_lib[] = {
		{"NewSimulation", lua_S4_NewSimulation},
		{NULL, NULL}
	};

	static const struct luaL_Reg lua_S4_Material_obj[] = {
		{"GetName", lua_S4_Material_GetName},
		{"SetName", lua_S4_Material_SetName},
		{"GetEpsilon", lua_S4_Material_GetEpsilon},
		{"SetEpsilon", lua_S4_Material_SetEpsilon},
		{NULL, NULL}
	};
	static const struct luaL_Reg lua_S4_Layer_obj[] = {
		{"GetName", lua_S4_Layer_GetName},
		{"SetName", lua_S4_Layer_SetName},
		{"GetThickness", lua_S4_Layer_GetThickness},
		{"SetThickness", lua_S4_Layer_SetThickness},
		{"ClearRegions", lua_S4_Layer_ClearRegions},
		{"SetRegion", lua_S4_Layer_SetRegion},
		/* Outputs */
		{"GetPowerFlux", lua_S4_Layer_GetPowerFlux},
		//{"GetPowerFluxes", lua_S4_Layer_GetPowerFluxes},
		//{"GetVolumeIntegral", lua_S4_Layer_GetVolumeIntegral},
		//{"GetZIntegral", lua_S4_Layer_GetZIntegral},
		//{"GetStressTensorIntegral", lua_S4_Layer_GetStressTensorIntegral},
		{"GetWaves", lua_S4_Layer_GetWaves},
		{NULL, NULL}
	};
	static const struct luaL_Reg lua_S4_Simulation_obj[] = {
		{"Clone", lua_S4_Simulation_Clone},
		{"SetLattice", lua_S4_Simulation_SetLattice},
		{"GetLattice", lua_S4_Simulation_GetLattice},
		{"SetBases", lua_S4_Simulation_SetBases},
		{"GetBases", lua_S4_Simulation_GetBases},
		{"SetFrequency", lua_S4_Simulation_SetFrequency},
		{"GetFrequency", lua_S4_Simulation_GetFrequency},
		{"AddMaterial", lua_S4_Simulation_AddMaterial},
		{"AddLayer", lua_S4_Simulation_AddLayer},
		{"GetMaterial", lua_S4_Simulation_GetMaterial},
		{"GetLayer", lua_S4_Simulation_GetLayer},
		{"ExcitationPlanewave", lua_S4_Simulation_ExcitationPlanewave},
		{"GetEpsilon", lua_S4_Simulation_GetEpsilon},
		{"GetFields", lua_S4_Simulation_GetFields},
		{NULL, NULL}
	};

	luaL_newmetatable(L, LUA_S4_MATERIAL_TYPENAME); /* stack: meta */
	luaL_newlib(L, lua_S4_Material_obj); /* stack: meta, lib */
	lua_setfield(L, -2, "__index"); /* stack: meta */
	lua_pushstring(L, "__gc"); /* stack: meta, "__gc" */
	lua_pushcfunction(L, &lua_S4_Material__gc); /* stack: meta, "__gc", __gc */
	lua_settable(L, -3); /* stack: meta */
	lua_pop(L, 1);

	luaL_newmetatable(L, LUA_S4_LAYER_TYPENAME); /* stack: meta */
	luaL_newlib(L, lua_S4_Layer_obj); /* stack: meta, lib */
	lua_setfield(L, -2, "__index"); /* stack: meta */
	lua_pushstring(L, "__gc"); /* stack: meta, "__gc" */
	lua_pushcfunction(L, &lua_S4_Layer__gc); /* stack: meta, "__gc", __gc */
	lua_settable(L, -3); /* stack: meta */
	lua_pop(L, 1);

	luaL_newmetatable(L, LUA_S4_SIMULATION_TYPENAME); /* stack: meta */
	luaL_newlib(L, lua_S4_Simulation_obj); /* stack: meta, lib */
	lua_setfield(L, -2, "__index"); /* stack: meta */
	lua_pushstring(L, "__gc"); /* stack: meta, "__gc" */
	lua_pushcfunction(L, &lua_S4_Simulation__gc); /* stack: meta, "__gc", __gc */
	lua_settable(L, -3); /* stack: meta */
	lua_pop(L, 1);

	luaL_newlib(L, S4v2_lib);
	return 1;
}

int lua_S4_len(lua_State *L, int arg){
#if LUA_VERSION_NUM < 502
	return lua_objlen(L, arg);
#else
	int n;
	if(arg < 0){ arg += lua_gettop(L)+1; }
	lua_len(L, -1);
	n = lua_tointeger(L, -1);
	lua_pop(L, 1);
	return n;
#endif
}
int lua_S4_isinteger(lua_State *L, int arg){
	lua_Number x;
	int i;
	if(!lua_isnumber(L, arg)){ return 0; }
	x = lua_tonumber(L, arg);
	i = (int)x;
	return ((lua_Number)i == x);
}
int lua_S4_getcomplex(lua_State *L, int arg, S4_real *v){
	int i;
	if(arg < 0){ arg += lua_gettop(L)+1; }
	if(!lua_istable(L, arg) || 2 != lua_S4_len(L, arg)){
		return luaL_error(L, "Expected complex number");
	}
	for(i = 0; i < 2; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, arg);
		if(!lua_isnumber(L, -1)){
			return luaL_error(L, "Expected complex number");
		}
		v[i] = lua_tonumber(L, -1);
		lua_pop(L, 1);
	}
	return 0;
}
int lua_S4_getvec3(lua_State *L, int arg, S4_real *v){
	int i;
	if(arg < 0){ arg += lua_gettop(L)+1; }
	if(!lua_istable(L, arg) || 3 != lua_S4_len(L, arg)){
		return luaL_error(L, "Expected 3-vector");
	}
	for(i = 0; i < 3; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, arg);
		if(!lua_isnumber(L, -1)){
			return luaL_error(L, "Expected 3-vector");
		}
		v[i] = lua_tonumber(L, -1);
		lua_pop(L, 1);
	}
	return 0;
}
int lua_S4_gettensor(lua_State *L, int arg, S4_real *eps){
	int i, j, k;
	if(arg < 0){ arg += lua_gettop(L)+1; }
	memset(eps, 0, sizeof(S4_real)*18);
	if(lua_isnumber(L, arg)){
		eps[0] = lua_tonumber(L, arg);
		return S4_MATERIAL_TYPE_SCALAR_REAL;
	}
	luaL_checktype(L, arg, LUA_TTABLE);
	int n1 = lua_S4_len(L, arg);
	if(2 == n1){
		for(i = 0; i < 2; ++i){
			lua_pushinteger(L, i+1);
			lua_gettable(L, arg);
			if(!lua_isnumber(L, -1)){
				return luaL_error(L, "Invalid format for scalar, complex, or tensor value (looks like a complex number)");
			}
			eps[i] = lua_tonumber(L, -1);
			lua_pop(L, 1);
		}
		return S4_MATERIAL_TYPE_SCALAR_COMPLEX;
	}
	if(3 == n1){
		for(i = 0; i < 3; ++i){
			lua_pushinteger(L, i+1);
			lua_gettable(L, arg);
			if(!lua_istable(L, -1)){
				return luaL_error(L, "Invalid format for scalar, complex, or tensor value (looks like a 3x3 matrix)");
			}
			for(j = 0; j < 3; ++j){
				lua_pushinteger(L, j+1);
				lua_gettable(L, -2);
				if(lua_isnumber(L, -1)){
					eps[(i+j*3)*2+0] = lua_tonumber(L, -1);
				}else if(lua_istable(L, -1)){
					if(2 != lua_S4_len(L, -1)){
						return luaL_error(L, "Invalid format for scalar, complex, or tensor value (looks like a 3x3 complex matrix)");
					}
					for(k = 0; k < 2; ++k){
						lua_pushinteger(L, k+1);
						lua_gettable(L, -2);
						if(!lua_isnumber(L, -1)){
							return luaL_error(L, "Invalid format for scalar, complex, or tensor value (looks like a 3x3 complex matrix)");
						}
						eps[(i+j*3)*2+k] = lua_tonumber(L, -1);
						lua_pop(L, 1);
					}
				}
				lua_pop(L, 1);
			}
			lua_pop(L, 1);
		}
		{
			S4_real abcde[10] = {
				eps[0], eps[1], eps[6], eps[7],
				eps[2], eps[3], eps[8], eps[9],
				eps[16], eps[17]
			};
			for(i = 0; i < 10; ++i){
				eps[i] = abcde[i];
			}
		}
		return S4_MATERIAL_TYPE_XYTENSOR_COMPLEX;
	}
	return 0;
}

lua_S4_Material* lua_S4_Material_check(lua_State *L, int narg){
	return (lua_S4_Material*)luaL_checkudata(L, narg, LUA_S4_MATERIAL_TYPENAME);
}
lua_S4_Layer* lua_S4_Layer_check(lua_State *L, int narg){
	return (lua_S4_Layer*)luaL_checkudata(L, narg, LUA_S4_LAYER_TYPENAME);
}
S4_Simulation* lua_S4_Simulation_check(lua_State *L, int narg){
	lua_S4_Simulation *pS = (lua_S4_Simulation*)luaL_checkudata(L, narg, LUA_S4_SIMULATION_TYPENAME);
	return pS->S;
}

int lua_S4_Material__gc(lua_State *L){
	lua_S4_Material *SM = lua_S4_Material_check(L, 1);
	luaL_unref(L, LUA_REGISTRYINDEX, SM->Sref);
	return 0;
}
int lua_S4_Layer__gc(lua_State *L){
	lua_S4_Layer *SL = lua_S4_Layer_check(L, 1);
	luaL_unref(L, LUA_REGISTRYINDEX, SL->Sref);
	return 0;
}
int lua_S4_Simulation__gc(lua_State *L){
	lua_S4_Simulation *pS = (lua_S4_Simulation*)luaL_checkudata(L, 1, LUA_S4_SIMULATION_TYPENAME);
	if(!pS->externally_owned){
		S4_Simulation_Destroy(pS->S);
	}
	if(NULL != pS->userdata){
		free(pS->userdata);
	}
	return 0;
}
int lua_S4_Simulation_push(lua_State *L, S4_Simulation *S, void *userdata, int ext_owned){
	lua_S4_Simulation *pS;
	pS = (lua_S4_Simulation*)lua_newuserdata(L, sizeof(lua_S4_Simulation));
	luaL_getmetatable(L, LUA_S4_SIMULATION_TYPENAME);
	lua_setmetatable(L, -2);
	pS->S = S;
	pS->userdata = userdata;
	pS->externally_owned = ext_owned;
	return 1;
}
int lua_S4_Material_push(lua_State *L, S4_Simulation *S, S4_MaterialID M){
	int ref = luaL_ref(L, LUA_REGISTRYINDEX);
	lua_S4_Material *SM = (lua_S4_Material*)lua_newuserdata(L, sizeof(lua_S4_Material));
	SM->S = S;
	SM->M = M;
	SM->Sref = ref;
	luaL_getmetatable(L, LUA_S4_MATERIAL_TYPENAME);
	lua_setmetatable(L, -2);
	return 1;
}
int lua_S4_Layer_push(lua_State *L, S4_Simulation *S, S4_LayerID pL){
	int ref = luaL_ref(L, LUA_REGISTRYINDEX);
	lua_S4_Layer *SL = (lua_S4_Layer*)lua_newuserdata(L, sizeof(lua_S4_Layer));
	SL->S = S;
	SL->L = pL;
	SL->Sref = ref;
	luaL_getmetatable(L, LUA_S4_LAYER_TYPENAME);
	lua_setmetatable(L, -2);
	return 1;
}
int lua_S4_Simulation_Material_check(lua_State *L, int index, S4_Simulation **SMS, S4_MaterialID *SMM){
	lua_S4_Material* SM = (lua_S4_Material*)luaL_checkudata(L, index, LUA_S4_MATERIAL_TYPENAME);
	*SMS = SM->S;
	*SMM = SM->M;
	return 0;
}
int lua_S4_Simulation_Layer_check(lua_State *L, int index, S4_Simulation **SLS, S4_LayerID *SLL){
	lua_S4_Layer* SL = (lua_S4_Layer*)luaL_checkudata(L, index, LUA_S4_LAYER_TYPENAME);
	*SLS = SL->S;
	*SLL = SL->L;
	return 0;
}
int lua_S4_Simulation_getmetatable(lua_State *L){
	luaL_getmetatable(L, LUA_S4_SIMULATION_TYPENAME);
	return 1;
}
int lua_S4_Layer_getmetatable(lua_State *L){
	luaL_getmetatable(L, LUA_S4_LAYER_TYPENAME);
	return 1;
}


int lua_S4_NewSimulation(lua_State *L){
	luaL_checktype(L, 1, LUA_TTABLE);
	S4_real Lr[4] = { 1, 0, 0, 1 };
	int nG = 1;
	int *G = NULL;
	unsigned int argflags = 0;

	lua_gc(L, LUA_GCCOLLECT, 0);
	lua_pushnil(L);
	while(0 != lua_next(L, 1)){
		const char *key;
		if(!lua_isstring(L, -2)){
			return luaL_argerror(L, 1, "Expected only named arguments");
		}
		key = lua_tostring(L, -2);
		if(0 == strcmp("lattice", key)){
			if(lua_isnumber(L, -1)){
				Lr[0] = lua_tonumber(L, -1);
				Lr[1] = 0;
				Lr[2] = 0;
				Lr[3] = 0;
			}else if(lua_istable(L, -1)){
				int i, j;
				if(2 != lua_S4_len(L, -1)){
					return luaL_argerror(L, 1, "lattice must be a number or a pair of vectors");
				}
				for(i = 0; i < 2; ++i){
					lua_pushinteger(L, i+1);
					lua_gettable(L, -2);
					if(!lua_istable(L, -1)){
						return luaL_argerror(L, 1, "lattice must be a number or a pair of vectors");
					}
					if(2 != lua_S4_len(L, -1)){
						return luaL_argerror(L, 1, "lattice must be a number or a pair of vectors");
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
			argflags |= 0x1;
		}else if(0 == strcmp("bases", key)){
			if(lua_S4_isinteger(L, -1)){
				nG = lua_tointeger(L, -1);
				if(nG <= 0){
					return luaL_argerror(L, 1, "bases must be a positive integer or an explicit basis set");
				}
				argflags |= 0x2;
			}else if(lua_istable(L, -1)){
				int i;
				nG = lua_S4_len(L, -1);
				if(nG <= 0){
					return luaL_argerror(L, 1, "bases must be a positive integer or an explicit basis set");
				}
				G = (int*)malloc(sizeof(int) * 2 * nG);
				for(i = 0; i < nG; ++i){
					G[2*i+0] = 0;
					G[2*i+1] = 0;

					lua_pushinteger(L, i+1);
					lua_gettable(L, -2);

					if(lua_S4_isinteger(L, -1)){
						G[2*i+0] = lua_tointeger(L, -1);
					}else if(lua_istable(L, -1)){
						int j;
						if(2 != lua_S4_len(L, -1)){
							return luaL_argerror(L, 1, "bases must be a positive integer or an explicit basis set");
						}
						for(j = 0; j < 2; ++j){
							lua_pushinteger(L, j+1);
							lua_gettable(L, -2);
							if(!lua_S4_isinteger(L, -1)){
								return luaL_argerror(L, 1, "bases must be a positive integer or an explicit basis set");
							}
							G[2*i+j] = lua_tointeger(L, -1);
							lua_pop(L, 1);
						}
					}

					lua_pop(L, 1);
				}
				argflags |= 0x2;
			}else{
			}
		}else{
			return luaL_error(L, "Unrecognized named argument: %s", key);
		}
		lua_pop(L, 1);
	}
	if(0x3 != argflags){
		return luaL_argerror(L, 1, "Must specify Lattice and BasisSize");
	}
	lua_S4_Simulation_push(L, S4_Simulation_New(Lr, nG, G), NULL, 0);
	free(G);
	return 1;
}

int lua_S4_Material_GetName(lua_State *L){
	const char *name;
	lua_S4_Material *SM = lua_S4_Material_check(L, 1);
	S4_Material_GetName(SM->S, SM->M, &name);
	lua_pushstring(L, name);
	return 1;
}
int lua_S4_Material_SetName(lua_State *L){
	lua_S4_Material *SM = lua_S4_Material_check(L, 1);
	const char *name = luaL_checkstring(L, 2);
	const char *prevname;
	S4_Material_GetName(SM->S, SM->M, &prevname);
	lua_pushstring(L, prevname);
	S4_Simulation_SetMaterial(SM->S, SM->M, name, 0, NULL);
	return 1;
}
int lua_S4_Material_GetEpsilon(lua_State *L){
	int i, j, k;
	S4_real eps[18];
	lua_S4_Material *SM = lua_S4_Material_check(L, 1);
	S4_Material_GetEpsilon(SM->S, SM->M, eps);
	lua_createtable(L, 3, 0);
	for(i = 0; i < 3; ++i){
		lua_pushinteger(L, i+1);
		lua_createtable(L, 3, 0);
		for(j = 0; j < 3; ++j){
			lua_pushinteger(L, j+1);
			lua_createtable(L, 2, 0);
			for(k = 0; k < 2; ++k){
				lua_pushinteger(L, k+1);
				lua_pushnumber(L, eps[k+(i+j*3)*2]);
				lua_settable(L, -3);
			}
			lua_settable(L, -3);
		}
		lua_settable(L, -3);
	}
	return 1;
}
int lua_S4_Material_SetEpsilon(lua_State *L){
	int i, j, k;
	S4_real eps[18] = { 0 };
	lua_S4_Material *SM = lua_S4_Material_check(L, 1);
	int type = lua_S4_gettensor(L, 2, eps);
	switch(type){
	case S4_MATERIAL_TYPE_SCALAR_REAL:
		S4_Simulation_SetMaterial(SM->S, SM->M, NULL, type, eps);
		break;
	case S4_MATERIAL_TYPE_SCALAR_COMPLEX:
		S4_Simulation_SetMaterial(SM->S, SM->M, NULL, type, eps);
		break;
	case S4_MATERIAL_TYPE_XYTENSOR_COMPLEX:
		S4_Simulation_SetMaterial(SM->S, SM->M, NULL, type, eps);
		break;
	default:
		return luaL_error(L, "Invalid format for epsilon");
	}
	return 0;
}

int lua_S4_Layer_GetName(lua_State *L){
	const char *name;
	lua_S4_Layer *SL = lua_S4_Layer_check(L, 1);
	S4_Layer_GetName(SL->S, SL->L, &name);
	lua_pushstring(L, name);
	return 1;
}
int lua_S4_Layer_SetName(lua_State *L){
	lua_S4_Layer *SL = lua_S4_Layer_check(L, 1);
	const char *name = luaL_checkstring(L, 2);
	const char *prevname;
	S4_Layer_GetName(SL->S, SL->L, &prevname);
	lua_pushstring(L, prevname);
	S4_Simulation_SetLayer(SL->S, SL->L, name, NULL, -1, -1);
	return 1;
}
int lua_S4_Layer_GetThickness(lua_State *L){
	lua_S4_Layer *SL = lua_S4_Layer_check(L, 1);
	S4_real thickness;
	S4_Layer_GetThickness(SL->S, SL->L, &thickness);
	lua_pushnumber(L, thickness);
	return 1;
}
int lua_S4_Layer_SetThickness(lua_State *L){
	lua_S4_Layer *SL = lua_S4_Layer_check(L, 1);
	S4_real thickness = luaL_checknumber(L, 2);
	S4_real prevthickness;
	S4_Layer_GetThickness(SL->S, SL->L, &prevthickness);
	S4_Simulation_SetLayer(SL->S, SL->L, NULL, &thickness, -1, -1);
	lua_pushnumber(L, prevthickness);
	return 1;
}
int lua_S4_Layer_SetRegion(lua_State *L){
	int type = S4_REGION_TYPE_INTERVAL;
	int i, j, nv;
	S4_real *v = NULL;
	S4_real center[2] = {0,0}, hw[2] = {0,0};
	S4_real angle_frac = 0;
	unsigned int argflags = 0;
	S4_MaterialID M = -1;

	lua_S4_Layer *SL = lua_S4_Layer_check(L, 1);
	luaL_checktype(L, 2, LUA_TTABLE);

	lua_pushnil(L);
	while(0 != lua_next(L, 2)){
		const char *key;
		if(!lua_isstring(L, -2)){
			return luaL_argerror(L, 2, "Expected only named arguments");
		}
		key = lua_tostring(L, -2);

		if(0 == strcmp("shape", key)){
			const char *shapestr;
			if(!lua_isstring(L, -1)){
				return luaL_argerror(L, 2, "Invalid format for shape");
			}
			shapestr = lua_tostring(L, -1);
			if(0 == strcmp("circle", shapestr)){
				type = S4_REGION_TYPE_CIRCLE;
			}else if(0 == strcmp("ellipse", shapestr)){
				type = S4_REGION_TYPE_ELLIPSE;
			}else if(0 == strcmp("rectangle", shapestr)){
				type = S4_REGION_TYPE_RECTANGLE;
			}else if(0 == strcmp("polygon", shapestr)){
				type = S4_REGION_TYPE_POLYGON;
			}
			argflags |= 0x01;
		}else if(0 == strcmp("center", key)){
			if(lua_isnumber(L, -1)){
				center[0] = lua_tonumber(L, -1);
			}else{
				if(!lua_istable(L, -1)){
					return luaL_argerror(L, 2, "Invalid format for center");
				}
				for(i = 0; i < 2; ++i){
					lua_pushinteger(L, i+1);
					lua_gettable(L, -2);
					if(!lua_isnumber(L, -1)){
						return luaL_argerror(L, 2, "Invalid format for center");
					}
					center[i] = lua_tonumber(L, -1);
					lua_pop(L, 1);
				}
			}
			argflags |= 0x02;
		}else if(0 == strcmp("angle_degrees", key)){
			if(!lua_isnumber(L, -1)){
				return luaL_argerror(L, 2, "Invalid format for angle_degrees");
			}
			angle_frac += lua_tonumber(L, -1) / 360.;
		}else if(0 == strcmp("angle_radians", key)){
			if(!lua_isnumber(L, -1)){
				return luaL_argerror(L, 2, "Invalid format for angle_radians");
			}
			angle_frac += lua_tonumber(L, -1) / (2*M_PI);
		}else if(0 == strcmp("halfwidths", key)){
			if(lua_isnumber(L, -1)){
				hw[0] = lua_tonumber(L, -1);
			}else{
				if(!lua_istable(L, -1)){
					return luaL_argerror(L, 2, "Invalid format for halfwidths");
				}
				for(i = 0; i < 2; ++i){
					lua_pushinteger(L, i+1);
					lua_gettable(L, -2);
					if(!lua_isnumber(L, -1)){
						return luaL_argerror(L, 2, "Invalid format for halfwidths");
					}
					hw[i] = lua_tonumber(L, -1);
					lua_pop(L, 1);
				}
			}
			argflags |= 0x04;
		}else if(0 == strcmp("radius", key)){
			if(!lua_isnumber(L, -1)){
				return luaL_argerror(L, 2, "Invalid format for radius");
			}
			hw[0] = lua_tonumber(L, -1);
			hw[1] = hw[0];
			argflags |= 0x08;
		}else if(0 == strcmp("vertices", key)){
			if(!lua_istable(L, -1)){
				return luaL_argerror(L, 2, "Invalid format for vertices");
			}
			nv = lua_S4_len(L, -1);
			if(nv < 2){
				return luaL_argerror(L, 2, "Invalid format for vertices");
			}
			v = (S4_real*)malloc(sizeof(S4_real) * 2 * nv);
			for(i = 0; i < nv; ++i){
				lua_pushinteger(L, i+1);
				lua_gettable(L, -2);
				if(!lua_istable(L, -1) || 2 != lua_S4_len(L, -1)){
					return luaL_argerror(L, 2, "Invalid format for vertices");
				}
				for(j = 0; j < 2; ++j){
					lua_pushinteger(L, j+1);
					lua_gettable(L, -2);
					if(!lua_isnumber(L, -1)){
						return luaL_argerror(L, 2, "Invalid format for vertices");
					}
					v[2*i+j] = lua_tonumber(L, -1);
					lua_pop(L, 1);
				}
				lua_pop(L, 1);
			}
			argflags |= 0x10;
		}else if(0 == strcmp("material", key)){
			if(lua_isstring(L, -1)){
				const char *matname;
				if(!lua_isstring(L, -1)){
					return luaL_argerror(L, 2, "Invalid format for material");
				}
				matname = lua_tostring(L, -1);
				M = S4_Simulation_GetMaterialByName(SL->S, matname);
				if(M < 0){
					return luaL_error(L, "Unknown material '%s'", matname);
				}
			}else{
				lua_S4_Material *SM = lua_S4_Material_check(L, -1);
				M = SM->M;
				if(M < 0){
					return luaL_error(L, "Invalid material");
				}
			}

			argflags |= 0x20;
		}else{
			return luaL_error(L, "Unrecognized argument: %s", key);
		}

		lua_pop(L, 1);
	}

	if(0 == (argflags & 0x20)){
		return luaL_argerror(L, 1, "Must specify material");
	}

	S4_real Lr[4];
	S4_Simulation_GetLattice(SL->S, Lr);
	int lattice1d = 0;
	if(0 == Lr[1] && 0 == Lr[2] && 0 == Lr[3]){
		lattice1d = 1;
	}

	if(S4_REGION_TYPE_INTERVAL == type){
		if(!lattice1d || 0 != angle_frac || 0 != center[1] || 0 != hw[1]){
			return luaL_argerror(L, 2, "Interval region must be used with 1d lattice and have zero y-components for center and halfwidths");
		}
		if(0 == (argflags & 0x02) || 0 == (argflags & 0x04)){
			return luaL_argerror(L, 2, "Interval region requires specifying center and halfwidths");
		}
		S4_Layer_SetRegionHalfwidths(SL->S, SL->L, M, type, hw, center, &angle_frac);
		return 0;
	}
	if(S4_REGION_TYPE_ELLIPSE == type || S4_REGION_TYPE_RECTANGLE == type){
		if(lattice1d){
			return luaL_argerror(L, 2, "Region type must be used with 2d lattice");
		}
		if(0 == (argflags & 0x02) || 0 == (argflags & 0x04)){
			return luaL_argerror(L, 2, "Region type requires specifying center and halfwidths");
		}
		S4_Layer_SetRegionHalfwidths(SL->S, SL->L, M, type, hw, center, &angle_frac);
		return 0;
	}
	if(S4_REGION_TYPE_CIRCLE == type){
		if(lattice1d){
			return luaL_argerror(L, 2, "Circle region must be used with 2d lattice");
		}
		if(0 == (argflags & 0x02) || 0 == (argflags & 0x08)){
			return luaL_argerror(L, 2, "Circle region requires specifying center and radius");
		}
		S4_Layer_SetRegionHalfwidths(SL->S, SL->L, M, type, hw, center, &angle_frac);
		return 0;
	}

	if(S4_REGION_TYPE_POLYGON == type){
		if(lattice1d){
			return luaL_argerror(L, 1, "Polygon region must be used with 2d lattice");
		}
		if(0 == (argflags & 0x10)){
			return luaL_argerror(L, 1, "Polygon region requires specifying vertices");
		}
		S4_Layer_SetRegionVertices(SL->S, SL->L, M, type, nv, v, center, &angle_frac);
		free(v);
		return 0;
	}
	free(v);
	return luaL_argerror(L, 1, "Unsupported region type");
}

int lua_S4_Layer_ClearRegions(lua_State *L){
	lua_S4_Layer *SL = lua_S4_Layer_check(L, 1);
	S4_Layer_ClearRegions(SL->S, SL->L);
	return 0;
}


int lua_S4_Simulation_Clone(lua_State *L){
	S4_Simulation *S = lua_S4_Simulation_check(L, 1);
	S4_Simulation **pS = (S4_Simulation**)lua_newuserdata(L, sizeof(S4_Simulation*));
	luaL_getmetatable(L, LUA_S4_SIMULATION_TYPENAME);
	lua_setmetatable(L, -2);
	*pS = S4_Simulation_Clone(S);
	return 1;
}
int lua_S4_Simulation_SetLattice(lua_State *L){
	S4_real Lr[4];
	S4_Simulation *S = lua_S4_Simulation_check(L, 1);
	if(lua_isnumber(L, 2)){
		S4_real period = lua_tonumber(L, 2);
		Lr[0] = period;
		Lr[1] = 0;
		Lr[2] = 0;
		Lr[3] = 0;
		S4_Simulation_SetLattice(S, Lr);
		return 0;
	}
	if(lua_istable(L, 2) && 2 == lua_S4_len(L, 2)){
		int i, j;
		for(i = 0; i < 2; ++i){
			lua_pushnumber(L, i+1);
			lua_gettable(L, 2);
			if(!lua_istable(L, -1) || 2 != lua_S4_len(L, -1)){
				return luaL_argerror(L, 2, "Invalid format for lattice");
			}
			for(j = 0; j < 2; ++j){
				lua_pushnumber(L, i+1);
				lua_gettable(L, 2);
				if(!lua_isnumber(L, -1)){
					return luaL_argerror(L, 2, "Invalid format for lattice");
				}
				Lr[i+j*2] = lua_tonumber(L, -1);
				lua_pop(L, 1);
			}
			lua_pop(L, 1);
		}
		S4_Simulation_SetLattice(S, Lr);
		return 0;
	}
	return luaL_argerror(L, 2, "Invalid format for lattice");
}
int lua_S4_Simulation_GetLattice(lua_State *L){
	S4_real Lr[4];
	int i, j;
	S4_Simulation *S = lua_S4_Simulation_check(L, 1);
	S4_Simulation_GetLattice(S, Lr);
	lua_createtable(L, 2, 0);
	for(i = 0; i < 2; ++i){
		lua_pushinteger(L, i+1);
		lua_createtable(L, 2, 0);
		for(j = 0; j < 2; ++j){
			lua_pushinteger(L, j+1);
			lua_pushnumber(L, Lr[i+j*2]);
			lua_settable(L, -3);
		}
		lua_settable(L, -3);
	}
	return 1;
}
int lua_S4_Simulation_SetBases(lua_State *L){
	int n;
	S4_Simulation *S = lua_S4_Simulation_check(L, 1);
	if(lua_S4_isinteger(L, 2)){
		n = lua_tointeger(L, 2);
		S4_Simulation_SetBases(S, n, NULL);
		return 0;
	}else if(lua_istable(L, 2)){
		int *G;
		int i, j;
		n = lua_S4_len(L, 2);
		G = (int*)malloc(sizeof(int)*2*n);
		for(i = 0; i < n; ++i){
			lua_pushinteger(L, i+1);
			lua_gettable(L, 2);
			if(lua_istable(L, -1)){
				for(j = 0; j < 2; ++j){
					lua_pushinteger(L, i+1);
					lua_gettable(L, 2);
					if(!lua_S4_isinteger(L, -1)){
						return luaL_argerror(L, 2, "Invalid format for basis");
					}
					G[2*i+j] = lua_tointeger(L, -1);
					lua_pop(L, 1);
				}
			}else if(lua_S4_isinteger(L, -1)){
				G[2*i+0] = lua_tointeger(L, -1);
				G[2*i+1] = 0;
			}else{
				return luaL_argerror(L, 2, "Invalid format for basis");
			}
			lua_pop(L, 1);
		}
		free(G);
		return 0;
	}else{
		return luaL_argerror(L, 2, "Invalid format for basis");
	}
}
int lua_S4_Simulation_GetBases(lua_State *L){
	S4_Simulation *S = lua_S4_Simulation_check(L, 1);
	int i, j, n;
	int *G = NULL;
	S4_real Lr[4];
	n = S4_Simulation_GetBases(S, G);
	G = (int*)malloc(sizeof(int)*2*n);
	S4_Simulation_GetBases(S, G);
	S4_Simulation_GetLattice(S, Lr);
	lua_createtable(L, n, 0);
	if(0 == Lr[1] && 0 == Lr[2] && 0 == Lr[3]){
		for(i = 0; i < n; ++i){
			lua_pushinteger(L, i+1);
			lua_pushinteger(L, G[2*i+0]);
			lua_settable(L, -3);
		}
	}else{
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
	}
	free(G);
	return 1;
}
int lua_S4_Simulation_SetFrequency(lua_State *L){
	S4_real freq[2] = { 0, 0 };
	S4_Simulation *S = lua_S4_Simulation_check(L, 1);
	if(lua_isnumber(L, 2)){
		freq[0] = lua_tonumber(L, 2);
		if(lua_isnumber(L, 3)){
			freq[1] = lua_tonumber(L, 3);
		}else if(!lua_isnoneornil(L, 3)){
			return luaL_argerror(L, 2, "Invalid format for frequency");
		}
		S4_real freq_prev[2];
		S4_Simulation_GetFrequency(S, freq_prev);
		S4_Simulation_SetFrequency(S, freq);
		lua_pushnumber(L, freq_prev[0]);
		lua_pushnumber(L, freq_prev[1]);
		return 2;
	}else{
		return luaL_argerror(L, 2, "Invalid format for frequency");
	}
}
int lua_S4_Simulation_GetFrequency(lua_State *L){
	S4_Simulation *S = lua_S4_Simulation_check(L, 1);
	S4_real freq[2];
	S4_Simulation_GetFrequency(S, freq);
	lua_pushnumber(L, freq[0]);
	lua_pushnumber(L, freq[1]);
	return 2;
}
int lua_S4_Simulation_AddMaterial(lua_State *L){
	S4_Simulation *S = lua_S4_Simulation_check(L, 1);
	luaL_checktype(L, 2, LUA_TTABLE);
	const char *name = NULL;
	S4_real eps[18];
	unsigned int argflags = 0;
	lua_S4_Material *SM;
	S4_MaterialID M;
	int type = 0;

	lua_pushnil(L);
	while(0 != lua_next(L, 2)){
		const char *key;
		if(!lua_isstring(L, -2)){
			return luaL_argerror(L, 1, "Expected only named arguments");
		}
		key = lua_tostring(L, -2);

		if(0 == strcmp("name", key)){
			if(!lua_isstring(L, -1)){
				return luaL_argerror(L, 1, "Invalid format for name");
			}
			name = lua_tostring(L, -1);
		}else if(0 == strcmp("epsilon", key)){
			type = lua_S4_gettensor(L, -1, eps);
			argflags |= 0x01;
		}else{
			return luaL_error(L, "Unrecognized argument: %s", key);
		}

		lua_pop(L, 1);
	}
	if(0 == (argflags & 0x01)){
		return luaL_argerror(L, 1, "Must specify epsilon");
	}
	M = S4_Simulation_SetMaterial(S, -1, name, type, eps);
	if(M < 0){
		return luaL_argerror(L, 1, "Failed to add material");
	}
	lua_pushvalue(L, 1);
	lua_S4_Material_push(L, S, M);
	return 1;
}
int lua_S4_Simulation_AddLayer(lua_State *L){
	S4_Simulation *S = lua_S4_Simulation_check(L, 1);
	luaL_checktype(L, 2, LUA_TTABLE);
	const char *name = NULL;
	S4_MaterialID M = -1;
	unsigned int argflags = 0;
	lua_S4_Layer *SL;
	S4_LayerID layer = -1;
	S4_LayerID copy = -1;
	S4_real thickness = 0;

	lua_pushnil(L);
	while(0 != lua_next(L, 2)){
		const char *key;
		if(!lua_isstring(L, -2)){
			return luaL_argerror(L, 2, "Expected only named arguments");
		}
		key = lua_tostring(L, -2);

		if(0 == strcmp("name", key)){
			if(!lua_isstring(L, -1)){
				return luaL_argerror(L, 2, "Invalid format for name");
			}
			name = lua_tostring(L, -1);
		}else if(0 == strcmp("material", key)){
			if(lua_isstring(L, -1)){
				const char *matname = lua_tostring(L, -1);
				M = S4_Simulation_GetMaterialByName(S, matname);
				if(M < 0){
					return luaL_error(L, "Unknown material '%s'", matname);
				}
			}else{
				lua_S4_Material *SM = lua_S4_Material_check(L, -1);
				M = SM->M;
				if(M < 0){
					return luaL_error(L, "Invalid material");
				}
			}
			argflags |= 0x01;
		}else if(0 == strcmp("thickness", key)){
			if(!lua_isnumber(L, -1)){
				return luaL_argerror(L, 2, "Invalid format for thickness");
			}
			thickness = lua_tonumber(L, -1);
			if(thickness < 0){
				return luaL_argerror(L, 2, "Thickness must be positive");
			}
			argflags |= 0x02;
		}else if(0 == strcmp("copy", key)){
			if(lua_isstring(L, -1)){
				const char *lname;
				if(!lua_isstring(L, -1)){
					return luaL_argerror(L, 2, "Invalid format for copy");
				}
				lname = lua_tostring(L, -1);
				copy = S4_Simulation_GetLayerByName(S, lname);
				if(copy < 0){
					return luaL_error(L, "Unknown layer '%s'", lname);
				}
			}else{
				lua_S4_Layer *SL = lua_S4_Layer_check(L, -1);
				copy = SL->L;
				if(copy < 0){
					return luaL_error(L, "Invalid layer");
				}
			}
			argflags |= 0x04;
		}else{
			return luaL_error(L, "Unrecognized argument: %s", key);
		}

		lua_pop(L, 1);
	}
	if(copy < 0){
		if(0x3 != (argflags & 0x03)){
			return luaL_argerror(L, 1, "Must specify material and thickness");
		}
		layer = S4_Simulation_SetLayer(S, -1, name, &thickness, -1, M);
	}else{
		if(0x0 != (argflags & 0x01)){
			return luaL_argerror(L, 1, "Must not specify material for a layer copy");
		}
		layer = S4_Simulation_SetLayer(S, -1, name, &thickness, copy, -1);
	}
	if(layer < 0){
		return luaL_argerror(L, 1, "Failed to add layer");
	}
	lua_pushvalue(L, 1);
	return lua_S4_Layer_push(L, S, layer);
}

int lua_S4_Simulation_GetMaterial(lua_State *L){
	S4_Simulation *S = lua_S4_Simulation_check(L, 1);
	const char *name = luaL_checkstring(L, 2);
	return lua_S4_Material_push(L, S, S4_Simulation_GetMaterialByName(S, name));
}

int lua_S4_Simulation_GetLayer(lua_State *L){
	S4_Simulation *S = lua_S4_Simulation_check(L, 1);
	const char *name = luaL_checkstring(L, 2);
	return lua_S4_Layer_push(L, S, S4_Simulation_GetLayerByName(S, name));
}
int lua_S4_Simulation_ExcitationPlanewave(lua_State *L){
	S4_Simulation *S = lua_S4_Simulation_check(L, 1);
	luaL_checktype(L, 2, LUA_TTABLE);
	unsigned int argflags = 0;
	S4_real k[3], u[3], cu[2], cv[2];

	lua_pushnil(L);
	while(0 != lua_next(L, 2)){
		const char *key;
		if(!lua_isstring(L, -2)){
			return luaL_argerror(L, 2, "Expected only named arguments");
		}
		key = lua_tostring(L, -2);

		if(0 == strcmp("k", key)){
			lua_S4_getvec3(L, -1, k);
			argflags |= 0x01;
		}else if(0 == strcmp("u", key)){
			lua_S4_getvec3(L, -1, u);
			argflags |= 0x02;
		}else if(0 == strcmp("cu", key)){
			lua_S4_getcomplex(L, -1, cu);
			argflags |= 0x04;
		}else if(0 == strcmp("cv", key)){
			lua_S4_getcomplex(L, -1, cv);
			argflags |= 0x08;
		}else{
			return luaL_error(L, "Unrecognized argument: %s", key);
		}
		lua_pop(L, 1);
	}
	if(0x0F != argflags){
		return luaL_argerror(L, 2, "Must specify k, u, cu, cv");
	}
	S4_Simulation_ExcitationPlanewave(S, k, u, cu, cv);
	return 0;
}
int lua_S4_Simulation_GetEpsilon(lua_State *L){
	int n;
	S4_Simulation *S = lua_S4_Simulation_check(L, 1);
	luaL_checktype(L, 2, LUA_TTABLE);

	n = lua_S4_len(L, 2);
	if(3 == n){
		int i;
		double r[3], eps[2];
		for(i = 0; i < 3; ++i){
			lua_pushinteger(L, i+1);
			lua_gettable(L, 2);
			r[i] = lua_tonumber(L, -1);
			lua_pop(L, 1);
		}
		Simulation_GetEpsilon(S, r, eps);
		lua_pushnumber(L, eps[0]);
		lua_pushnumber(L, eps[1]);
		return 2;
	}else{
		return luaL_argerror(L, 2, "Must specify query point");
	}
}
int lua_S4_Simulation_GetFields(lua_State *L){
	int n = 0;
	
	S4_Simulation *S = lua_S4_Simulation_check(L, 1);
	luaL_checktype(L, 2, LUA_TTABLE);
	
	n = lua_S4_len(L, 2);
	
	if(3 == n){ // evaluation at a single point
		static const int n11[2] = { 1, 1 };
		int i;
		S4_real r[3];
		S4_real EH[12];
		for(i = 0; i < 3; ++i){
			lua_pushinteger(L, i+1);
			lua_gettable(L, 2);
			if(!lua_isnumber(L, -1)){
				return luaL_argerror(L, 2, "Evaluation point must be numeric");
			}
			r[i] = lua_tonumber(L, -1);
			lua_pop(L, 1);
		}
		S4_Simulation_GetFieldPlane(S, n11, r, &EH[0], &EH[6]);
		
		for(i = 0; i < 2; ++i){
			int j;
			lua_createtable(L, 3, 0);
			for(j = 0; j < 3; ++j){
				int k;
				lua_pushinteger(L, j+1);
				lua_createtable(L, 2, 0);
				for(k = 0; k < 2; ++k){
					lua_pushinteger(L, k+1);
					lua_pushnumber(L, EH[(3*i+j)*2+k]);
					lua_settable(L, -3);
				}
				lua_settable(L, -3);
			}
		}
		return 2;
	}
	return 0;
}

int lua_S4_Layer_GetPowerFlux(lua_State *L){
	S4_real power[4];
	int i;
	lua_S4_Layer *SL = lua_S4_Layer_check(L, 1);
	S4_Simulation_GetPowerFlux(SL->S, SL->L, NULL, power);
	for(i = 0; i < 4; ++i){
		lua_pushnumber(L, power[i]);
	}
	return 4;
}
int lua_S4_Layer_GetWaves(lua_State *L){
	lua_S4_Layer *SL = lua_S4_Layer_check(L, 1);
	S4_real *waves = NULL;
	int n, i, j, k, ret;

	n = S4_Simulation_GetBases(SL->S, NULL);
	const int n2 = 2*n;

	waves = (S4_real*)malloc(sizeof(S4_real)*11*n2);
	S4_Simulation_GetWaves(SL->S, SL->L, waves);

	for(j = 0; j < 2; ++j){
		lua_createtable(L, n, 0);
		for(i = 0; i < n; ++i){
			const double *wave = &waves[11*(i+n*j)];
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
	}
	free(waves);
	return 2;
}
