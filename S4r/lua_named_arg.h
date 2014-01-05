#ifndef LUA_NAMED_ARG_H_INCLUDED
#define LUA_NAMED_ARG_H_INCLUDED

#include <lua.h>

/* Checks to make sure there is a table of named arguments
 * at stack index `index'.
 */
int lua_named_arg_check(
	lua_State *L, int index, int nreq, int nopt, const char *args[]
);

/* Produces an error message for the named argument `name'
 * with additional message `extramsg'. If name is NULL, then
 * the error indicates the argument number specified by `index'.
 */
int lua_named_arg_error(
	lua_State *L, int index, const char *name, const char *extramsg
);

/* Parses the object at stack index `i' as a 2-vector,
 * placing the values into `v'. If the object was a named
 * argument, then an optional `argname' may be provided for
 * more helpful error messages.
 */
int lua_arg_parse_vec2(lua_State *L, int i, const char *argname, double *v);
int lua_arg_parse_vec3(lua_State *L, int i, const char *argname, double *v);
int lua_arg_parse_ctensor3(lua_State *L, int index, const char *argname, double *m);

/* For a table of named arguments at stack index `index',
 * returns a string with argument name `argname'.
 */
int lua_named_arg_get_string(
	lua_State *L, int index, const char *argname,
	const char **value
);
int lua_named_arg_get_number(
	lua_State *L, int index, const char *argname,
	double *value
);

/* Check to see if a named argument exists
 * (for optional arguments)
 */
int lua_named_arg_exists(lua_State *L, int index, const char *argname);

/* Pushes the named argument onto the stack if it exists and returns 1,
 * otherwise returns 0.
 */
int lua_named_arg_get(lua_State *L, int index, const char *argname);

#endif /* LUA_NAMED_ARG_H_INCLUDED */
