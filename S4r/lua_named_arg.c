#include "lua_named_arg.h"
#include <lualib.h>
#include <lauxlib.h>
#include <string.h>
#include <stdio.h>

int lua_named_arg_error(lua_State *L, int narg, const char *name, const char *extramsg){
	lua_Debug ar;
	if(NULL == name){ return luaL_argerror(L, narg, extramsg); }
	if(!lua_getstack(L, 0, &ar)){ /* no stack frame? */
		return luaL_error(L, "bad argument " LUA_QS " (%s)", name, extramsg);
	}
	lua_getinfo(L, "n", &ar);
	if(0 == strcmp(ar.namewhat, "method")){
		narg--;  /* do not count `self' */
		if(0 == narg){ /* error is in the self argument itself? */
			return luaL_error(L, "calling " LUA_QS " on bad self (%s)", ar.name, extramsg);
		}
	}
	if(NULL == ar.name){
		ar.name = "?";
	}
	return luaL_error(L, "bad argument " LUA_QS " to " LUA_QS " (%s)", name, ar.name, extramsg);
}

// Parses the object at stack index `index' as a 2-vector, placing the values into `v'
// If the object was a named argument, then an optional `argname' may be provided for
// more helpful error messages.
int lua_arg_parse_vec2(lua_State *L, int index, const char *argname, double *v){
	int i;
	if(!lua_istable(L, index) || 2 != lua_rawlen(L, index)){
		return lua_named_arg_error(L, index, argname, "expected argument type: vec2");
	}
	for(i = 0; i < 2; ++i){
		lua_rawgeti(L, index, i+1);
		if(!lua_isnumber(L, -1)){
			return lua_named_arg_error(L, index, argname, "argument of type vec2 does not contain numbers");
		}
		v[i] = lua_tonumber(L, -1);
		lua_pop(L, 1);
	}
	return 0;
}
int lua_arg_parse_vec3(lua_State *L, int index, const char *argname, double *v){
	int i;
	if(!lua_istable(L, index) || 3 != lua_rawlen(L, index)){
		return lua_named_arg_error(L, index, argname, "expected type: vec2");
	}
	for(i = 0; i < 3; ++i){
		lua_rawgeti(L, index, i+1);
		if(!lua_isnumber(L, -1)){
			return lua_named_arg_error(L, index, argname, "argument of type vec2 does not contain numbers");
		}
		v[i] = lua_tonumber(L, -1);
		lua_pop(L, 1);
	}
	return 0;
}

#define lua_arg_parse_ctensor3_err() \
	lua_named_arg_error(L, index, argname, "expected type: scalar, complex, or ctensor3")
	
int lua_arg_parse_ctensor3(lua_State *L, int index, const char *argname, double *m){
	int i, j, k, l, ikind = 0;
	for(i = 0; i < 18; ++i){ m[i] = 0.; }
	if(lua_isnumber(L, index)){
		m[0] = lua_tonumber(L, index);
		m[8] = m[0];
		m[16] = m[0];
		return 0;
	}else if(!lua_istable(L, index)){
		return lua_arg_parse_ctensor3_err();
	}
	
	if(index < 0){ index += 1+lua_gettop(L); }
	lua_len(L, index); l = lua_tointeger(L, -1); lua_pop(L, 1);
	if(2 == l){ /* complex number */
		for(k = 0; k < 2; ++k){
			lua_pushinteger(L, k+1);
			lua_gettable(L, index);
			if(!lua_isnumber(L, -1)){
				return lua_named_arg_error(L, index, argname, "non-numeric value in what looks like a complex number");
			}
			m[k] = lua_tonumber(L, -1);
			lua_pop(L, 1);
		}
		m[8] = m[0];
		m[9] = m[1];
		m[16] = m[0];
		m[17] = m[1];
		return 0;
	}else if(3 != l){
		return lua_arg_parse_ctensor3_err();
	}
	/* Now we should expect a 3 vector or a 3x3 matrix */
	/* Iterate over the rows */
	for(i = 0; i < 3; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, index);
		
		if(0 == ikind){
			if(lua_isnumber(L, -1)){
				ikind = 1;
			}else if(!lua_istable(L, -1)){
				return lua_arg_parse_ctensor3_err();
			}else{
				lua_len(L, -1); l = lua_tointeger(L, -1); lua_pop(L, 1);
				if(2 == l){
					ikind = 2;
				}else if(3 == l){
					ikind = 3;
				}else{
					return lua_arg_parse_ctensor3_err();
				}
			}
		}
		
		if(1 == ikind){ /* real diagonal (3-vector) */
			if(!lua_isnumber(L, -1)){ return lua_arg_parse_ctensor3_err(); }
			m[8*i] = lua_tonumber(L, -1);
		}else if(2 == ikind){ /* complex diagonal (3-cvector) */
			if(!lua_istable(L, -1)){ return lua_arg_parse_ctensor3_err(); }
			lua_len(L, -1); l = lua_tointeger(L, -1); lua_pop(L, 1);
			if(2 != l){ return lua_arg_parse_ctensor3_err(); }
			for(k = 0; k < 2; ++k){
				lua_pushinteger(L, k+1);
				lua_gettable(L, -2);
				m[8*i+k] = lua_tonumber(L, -1);
				lua_pop(L, 1);
			}
		}else if(3 == ikind){ /* 3x3 matrix */
			int jkind = 0;
			if(!lua_istable(L, -1)){ return lua_arg_parse_ctensor3_err(); }
			lua_len(L, -1); l = lua_tointeger(L, -1); lua_pop(L, 1);
			if(3 != l){ return lua_arg_parse_ctensor3_err(); }
			for(j = 0; j < 3; ++j){
				lua_pushinteger(L, j+1);
				lua_gettable(L, -2);
				
				if(0 == jkind){
					if(lua_isnumber(L, -1)){
						jkind = 1;
					}else if(!lua_istable(L, -1)){
						return lua_arg_parse_ctensor3_err();
					}else{
						lua_len(L, -1); l = lua_tointeger(L, -1); lua_pop(L, 1);
						if(2 == l){
							jkind = 2;
						}else{
							return lua_arg_parse_ctensor3_err();
						}
					}
				}
				
				if(1 == jkind){ /* real diagonal (3-vector) */
					if(!lua_isnumber(L, -1)){ return lua_arg_parse_ctensor3_err(); }
					m[(i+j*3)*2] = lua_tonumber(L, -1);
				}else if(2 == jkind){ /* complex diagonal (3-cvector) */
					if(!lua_istable(L, -1)){ return lua_arg_parse_ctensor3_err(); }
					lua_len(L, -1); l = lua_tointeger(L, -1); lua_pop(L, 1);
					if(2 != l){ return lua_arg_parse_ctensor3_err(); }
					for(k = 0; k < 2; ++k){
						lua_pushinteger(L, k+1);
						lua_gettable(L, -2);
						m[k+(i+j*3)*2] = lua_tonumber(L, -1);
						lua_pop(L, 1);
					}
				}
				lua_pop(L, 1);
			}
		}
		
		lua_pop(L, 1);
		i++;
	}
	
	return 0;
}

/* Checks to make sure there is a table of named arguments at stack index `index'. */
int lua_named_arg_check(lua_State *L, int index, int nreq, int nopt, const char *args[]){
	int i;
	if(!lua_istable(L, index)){
		return luaL_argerror(L, index, "expected named arguments");
	}
	
	lua_pushnil(L);
	while(lua_next(L, index) != 0){
		/* key at index -2 and value at index -1 */
		const char *key;
		int found = 0;
		if(!lua_isstring(L, -2)){
			return luaL_argerror(L, index, "non-string key in named argument table");
		}
		key = lua_tostring(L, -2);
		for(i = 0; i < nreq+nopt; ++i){
			if(NULL != args[i] && 0 == strcmp(key, args[i])){
				args[i] = NULL;
				found = 1;
				break;
			}
		}
		if(!found){
			return luaL_error(L, "unexpected argument: " LUA_QS, key);
		}
		lua_pop(L, 1); /* remove value, keep key for next iteration */
	}
	for(i = 0; i < nreq; ++i){
		if(NULL != args[i]){
			return luaL_error(L, "missing argument: " LUA_QS, args[i]);
		}
	}
	return 0;
}

/* For a table of named arguments at stack index `index', returns a string with
 * argument name `argname'.
 */
int lua_named_arg_get_string(lua_State *L, int index, const char *argname, const char **value){
	lua_getfield(L, index, argname);
	if(!lua_isstring(L, -1)){
		return lua_named_arg_error(L, index, argname, "expected type: string");
	}
	*value = lua_tostring(L, -1);
	lua_pop(L, 1);
	return 0;
}
int lua_named_arg_get_number(lua_State *L, int index, const char *argname, double *value){
	lua_getfield(L, index, argname);
	if(!lua_isnumber(L, -1)){
		return lua_named_arg_error(L, index, argname, "expected type: number");
	}
	*value = lua_tonumber(L, -1);
	lua_pop(L, 1);
	return 0;
}
int lua_named_arg_exists(lua_State *L, int index, const char *argname){
	int ret;
	lua_getfield(L, index, argname);
	ret = !lua_isnil(L, -1);
	lua_pop(L, 1);
	return ret;
}
int lua_named_arg_get(lua_State *L, int index, const char *argname){
	lua_getfield(L, index, argname);
	if(lua_isnil(L, -1)){
		lua_pop(L, 1);
		return 0;
	}
	return 1;
}
