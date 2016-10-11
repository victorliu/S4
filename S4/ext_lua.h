#ifndef S4_EXT_LUA_H_INCLUDED
#define S4_EXT_LUA_H_INCLUDED

#include "S4.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

int luaopen_S4v2(lua_State *L);

S4_Simulation* lua_S4_Simulation_check(lua_State *L, int narg);
int lua_S4_Simulation_push(lua_State *L, S4_Simulation *S, void *userdata, int ext_owned);
int lua_S4_Material_push(lua_State *L, S4_Simulation *S, S4_MaterialID M);
int lua_S4_Layer_push(lua_State *L, S4_Simulation *S, S4_LayerID pL);
int lua_S4_Simulation_Material_check(lua_State *L, int index, S4_Simulation **pS, S4_MaterialID *pM);
int lua_S4_Simulation_Layer_check(lua_State *L, int index, S4_Simulation **pS, S4_LayerID *pL);
int lua_S4_Simulation_getmetatable(lua_State *L);
int lua_S4_Layer_getmetatable(lua_State *L);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* S4_EXT_LUA_H_INCLUDED */
