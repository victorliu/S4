#include <S4r/S4r.hpp>
#include <cstdio>
using namespace S4r;

extern "C" {
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
#include "lua_named_arg.h"
}

int lua_pushcomplex(lua_State *L, const doublecomplex &z){
	lua_createtable(L, 2, 5);
	
	lua_pushinteger(L, 1);
	lua_pushnumber(L, z.real());
	lua_settable(L, -3);
	
	lua_pushinteger(L, 2);
	lua_pushnumber(L, z.imag());
	lua_settable(L, -3);
	
	lua_pushnumber(L, z.real());
	lua_setfield(L, -2, "real");
	lua_pushnumber(L, z.imag());
	lua_setfield(L, -2, "imag");
	
	lua_pushnumber(L, std::abs(z));
	lua_setfield(L, -2, "abs");
	lua_pushnumber(L, std::norm(z));
	lua_setfield(L, -2, "abs2");
	lua_pushnumber(L, std::arg(z));
	lua_setfield(L, -2, "arg");
	
	return 1;
}

class LuaSimulation{
	static const char className[];

    static Simulation *checkarg(lua_State *L, int narg){
		luaL_checktype(L, narg, LUA_TUSERDATA);
		void *ud = luaL_checkudata(L, narg, className);
		if(!ud){
			luaL_argerror(L, narg, "expected Simulation object");
			return NULL;
		}
		return *(Simulation**)ud;
    }
    static int gc(lua_State *L) {
		Simulation *a = *(Simulation**)(lua_touserdata(L, 1));
		delete a;
		return 0;
    }
public:
    static int create(lua_State *L){
		const char *args[] = { "Lattice", "Mesh", "NumModes" };
		lua_named_arg_check(L, 1, 2, 1, args);

		// Get the lattice
		double Lr[4];
		lua_named_arg_get(L, 1, "Lattice");
		if(!lua_istable(L, -1) || 2 != lua_rawlen(L, -1)){
			return lua_named_arg_error(L, 1, "Lattice", "must be a pair of vec2's");
		}
		for(unsigned int i = 0; i < 2; ++i){
			lua_pushinteger(L, i+1);
			lua_gettable(L, -2);
			if(!lua_istable(L, -1) || 2 != lua_rawlen(L, -1)){
				return lua_named_arg_error(L, 1, "Lattice", "must be a pair of vec2's");
			}
			for(unsigned j = 0; j < 2; ++j){
				lua_pushinteger(L, j+1);
				lua_gettable(L, -2);
				if(!lua_isnumber(L, -1)){
					return lua_named_arg_error(L, 1, "Lattice", "non-numeric value in a vec2");
				}
				Lr[j+i*2] = lua_tonumber(L, -1);
				lua_pop(L, 1);
			}
			lua_pop(L, 1);
		}
		lua_pop(L, 1);
		
		// Get the mesh specifications
		size_t ngrid[2];
		lua_named_arg_get(L, 1, "Mesh");
		if(lua_isnumber(L, -1)){
			int ni = lua_tointeger(L, -1);
			if(ni < 1){
				return lua_named_arg_error(L, 1, "Mesh", "number of grid cells must be positive");
			}
			ngrid[0] = ni;
			ngrid[1] = ni;
		}else if(lua_istable(L, -1)){
			lua_len(L, -1);
			int tablen = lua_tointeger(L, -1); lua_pop(L, 1);
			if(2 != tablen){
				return lua_named_arg_error(L, 1, "Mesh", "must be an integer or pair of integers");
			}
			for(unsigned int i = 0; i < 2; ++i){
				int ni;
				lua_rawgeti(L, -1, i+1);
				if(!lua_isnumber(L, -1)){
					return lua_named_arg_error(L, 1, "Mesh", "must be a pair of integers");
				}
				ni = lua_tointeger(L, -1);
				if(ni < 1){
					return lua_named_arg_error(L, 1, "Mesh", "number of grid cells must be positive");
				}
				ngrid[i] = ni;
				lua_pop(L, 1);
			}
		}else{
			return lua_named_arg_error(L, 1, "Mesh", "must be an integer or pair of integers");
		}
		lua_pop(L, 1);
		
		// Get the default number of modes
		size_t nmodes = 0;
		if(lua_named_arg_get(L, 1, "NumModes")){
			if(!lua_isnumber(L, -1)){
				return lua_named_arg_error(L, 1, "NumModes", "must be a non-negative integer");
			}
			nmodes = lua_tointeger(L, -1);
			lua_pop(L, 1);
		}
		
		Simulation *S = new Simulation(Lr, ngrid, nmodes);
		*(Simulation**)lua_newuserdata(L, sizeof(Simulation*)) = S;
		luaL_getmetatable(L, className);
		lua_setmetatable(L, -2);
		return 1;
    }
    static int SetMaterial(lua_State *L){
		Simulation *S = checkarg(L, 1);
		const char *args[] = { "Name", "Epsilon", "Mu" };
		lua_named_arg_check(L, 2, 1, 2, args);
		
		const char *name;
		lua_named_arg_get_string(L, 2, "Name", &name);
		
		static const doublecomplex zero(0.);
		static const doublecomplex one(1.);
		
		double eps[18] = {
			1,0, 0,0, 0,0,
			0,0, 1,0, 0,0,
			0,0, 0,0, 1,0
		};
		double mu[18] = {
			1,0, 0,0, 0,0,
			0,0, 1,0, 0,0,
			0,0, 0,0, 1,0
		};
		
		if(lua_named_arg_get(L, 2, "Epsilon")){
			lua_arg_parse_ctensor3(L, -1, "Epsilon", eps);
			lua_pop(L, 1);
		}
		if(lua_named_arg_get(L, 2, "Mu")){
			lua_arg_parse_ctensor3(L, -1, "Mu", mu);
			lua_pop(L, 1);
		}
		
		Material *mat = new Material(eps, mu);
		
		S->SetMaterial(name, mat);
		
		return 0;
	}
    static int AddLayer(lua_State *L){
		Simulation *S = checkarg(L, 1);
		const char *args[] = { "Name", "Thickness", "Material", "NumModes" };
		lua_named_arg_check(L, 2, 3, 1, args);
		
		const char *name, *matname;
		double thickness;
		lua_named_arg_get_string(L, 2, "Name", &name);
		lua_named_arg_get_number(L, 2, "Thickness", &thickness);
		lua_named_arg_get_string(L, 2, "Material", &matname);
		
		Layer *layer = S->AddLayer(name, thickness, matname);
		
		if(lua_named_arg_get(L, 2, "NumModes")){
			if(!lua_isnumber(L, -1)){
				return lua_named_arg_error(L, 2, "NumModes", "must be a non-negative integer");
			}
			layer->description.num_modes = lua_tointeger(L, -1);
			lua_pop(L, 1);
		}

		return 0;
	}
    static int SetRegionCircle(lua_State *L){
		Simulation *S = checkarg(L, 1);
		const char *args[] = { "Layer", "Material", "Center", "Radius" };
		lua_named_arg_check(L, 2, 4, 0, args);
	
		const char *layer, *matname;
		double c[2], r;
		lua_named_arg_get_string(L, 2, "Layer", &layer);
		lua_named_arg_get_string(L, 2, "Material", &matname);
		lua_named_arg_get_number(L, 2, "Radius", &r);
		lua_named_arg_get(L, 2, "Center");
		lua_arg_parse_vec2(L, -1, "Center", c); lua_pop(L, 1);
		
		S->SetRegionCircle(layer, matname, Vec2(c[0],c[1]), r);
		return 0;
	}
    static int SetRegionEllipse(lua_State *L){
		Simulation *S = checkarg(L, 1);
		const char *args[] = { "Layer", "Material", "Center", "Angle", "Semiaxes" };
		lua_named_arg_check(L, 2, 5, 0, args);
	
		const char *layer, *matname;
		double c[2], h[2], angle;
		lua_named_arg_get_string(L, 2, "Layer", &layer);
		lua_named_arg_get_string(L, 2, "Material", &matname);
		lua_named_arg_get_number(L, 2, "Angle", &angle);
		lua_named_arg_get(L, 2, "Center");
		lua_arg_parse_vec2(L, -1, "Center", c); lua_pop(L, 1);
		lua_named_arg_get(L, 2, "Semiaxes");
		lua_arg_parse_vec2(L, -1, "Semiaxes", h); lua_pop(L, 1);
		
		angle *= 180./M_PI;
		
		S->SetRegionEllipse(layer, matname, Vec2(c[0],c[1]), angle, Vec2(h[0],h[1]));
		return 0;
	}
    static int SetRegionRectangle(lua_State *L){
		Simulation *S = checkarg(L, 1);
		const char *args[] = { "Layer", "Material", "Center", "Angle", "Semiaxes" };
		lua_named_arg_check(L, 2, 5, 0, args);
	
		const char *layer, *matname;
		double c[2], h[2], angle;
		lua_named_arg_get_string(L, 2, "Layer", &layer);
		lua_named_arg_get_string(L, 2, "Material", &matname);
		lua_named_arg_get_number(L, 2, "Angle", &angle);
		lua_named_arg_get(L, 2, "Center");
		lua_arg_parse_vec2(L, -1, "Center", c); lua_pop(L, 1);
		lua_named_arg_get(L, 2, "Semiaxes");
		lua_arg_parse_vec2(L, -1, "Semiaxes", h); lua_pop(L, 1);
		
		angle *= 180./M_PI;
		
		S->SetRegionRectangle(layer, matname, Vec2(c[0],c[1]), angle, Vec2(h[0],h[1]));
		return 0;
	}
    static int SetRegionPolygon(lua_State *L){
		Simulation *S = checkarg(L, 1);
		const char *args[] = { "Layer", "Material", "Center", "Angle", "Vertices" };
		lua_named_arg_check(L, 2, 5, 0, args);
	
		const char *layer, *matname;
		double c[2], angle;
		lua_named_arg_get_string(L, 2, "Layer", &layer);
		lua_named_arg_get_string(L, 2, "Material", &matname);
		lua_named_arg_get_number(L, 2, "Angle", &angle);
		lua_named_arg_get(L, 2, "Center");
		lua_arg_parse_vec2(L, -1, "Center", c); lua_pop(L, 1);
		
		lua_named_arg_get(L, 2, "Vertices");
		if(!lua_istable(L, -1)){
			return lua_named_arg_error(L, 1, "Vertices", "must be a list of Vec2's");
		}
		std::vector<Vec2> vert;
		size_t len = 0;
		lua_len(L, -1); len = lua_tointeger(L, -1); lua_pop(L, 1);
		vert.reserve(len);
		for(size_t i = 0; i < len; ++i){
			lua_pushinteger(L, i+1);
			lua_gettable(L, -2);
			double v[2];
			lua_arg_parse_vec2(L, 2, "Vertices", v);
			vert.push_back(Vec2(v[0],v[1]));
			lua_pop(L, 1);
		}
		lua_pop(L, 1);
		
		angle *= 180./M_PI;
		
		S->SetRegionPolygon(layer, matname, Vec2(c[0],c[1]), angle, vert);
		return 0;
	}
    static int RemoveLayerRegions(lua_State *L){
		Simulation *S = checkarg(L, 1);
		const char *args[] = { "Layer" };
		lua_named_arg_check(L, 2, 1, 0, args);
		
		const char *layer;
		lua_named_arg_get_string(L, 2, "Layer", &layer);
		
		S->RemoveLayerRegions(layer);
		return 0;
	}
    static int SetLayerThickness(lua_State *L){
		Simulation *S = checkarg(L, 1);
		const char *args[] = { "Layer", "Thickness" };
		lua_named_arg_check(L, 2, 2, 0, args);
		
		const char *layer;
		double thickness;
		lua_named_arg_get_string(L, 2, "Layer", &layer);
		lua_named_arg_get_number(L, 2, "Thickness", &thickness);
		
		S->SetLayerThickness(layer, thickness);
		return 0;
	}
	
    static int SetK(lua_State *L){
		Simulation *S = checkarg(L, 1);
		double k[2];
		lua_arg_parse_vec2(L, 2, NULL, k);
		S->SetK(k);
		return 0;
	}
    static int SetFrequency(lua_State *L){
		Simulation *S = checkarg(L, 1);
		double f = luaL_checknumber(L, 2);
		S->SetFrequency(f);
		return 0;
	}
    static int SolveAllLayers(lua_State *L){
		Simulation *S = checkarg(L, 1);
		S->SolveAllLayers();
		return 0;
	}
	
	static int GetReciprocalLattice(lua_State *L){
		return 0;
	}
	static int GetEpsilon(lua_State *L){
		return 0;
	}
	static int OutputLayerDescription(lua_State *L){
		Simulation *S = checkarg(L, 1);
		const char *args[] = { "Layer", "Filename" };
		lua_named_arg_check(L, 2, 1, 1, args);
		
		const char *layer, *filename = NULL;
		lua_named_arg_get_string(L, 2, "Layer", &layer);
		if(lua_named_arg_exists(L, 2, "Filename")){
			lua_named_arg_get_string(L, 2, "Filename", &filename);
		}
		
		FILE *fp = stdout;
		if(NULL != filename){
			fp = fopen(filename, "wb");
		}
		S->OutputLayerDescription(layer, fp);
		if(NULL != filename){
			fclose(fp);
		}
		
		return 0;
	}
	
	static int GetPowerFlux(lua_State *L){
		Simulation *S = checkarg(L, 1);
		const char *args[] = { "Layer", "ZOffset" };
		lua_named_arg_check(L, 2, 1, 1, args);
		
		const char *layer;
		double z = 0;
		lua_named_arg_get_string(L, 2, "Layer", &layer);
		if(lua_named_arg_exists(L, 2, "ZOffset")){
			lua_named_arg_get_number(L, 2, "ZOffset", &z);
		}
		
		doublecomplex forw, back;
		S->GetPowerFlux(layer, z, forw, back);
		lua_pushnumber(L, forw.real());
		lua_pushnumber(L, back.real());
		lua_pushnumber(L, forw.imag());
		lua_pushnumber(L, back.imag());
		return 4;
	}
	static int GetFields(lua_State *L){
		Simulation *S = checkarg(L, 1);
		luaL_argcheck(L, 2, lua_istable(L, 2) && (lua_rawlen(L, 2) == 3), "must be a 3-vector");
		
		Vec3 r;
		for(unsigned i = 0; i < 3; ++i){
			lua_pushinteger(L, i+1);
			lua_gettable(L, 2);
			luaL_argcheck(L, 2, lua_isnumber(L, -1), "must be a 3-vector");
			r[i] = lua_tonumber(L, -1);
			lua_pop(L, 1);
		}
		
		CVec3 eh[2];
		S->GetFields(r, eh[0], eh[1]);
		
		for(unsigned j = 0; j < 2; ++j){
			lua_createtable(L, 3, 0);
			for(unsigned i = 0; i < 3; ++i){
				lua_pushinteger(L, i+1);
				lua_pushcomplex(L, eh[j][i]);
				lua_settable(L, -3);
			}
		}
		
		return 2;
	}

	static void Register(lua_State* L){
		static const luaL_Reg SimulationLib[] = {
			{"__gc", &LuaSimulation::gc},
			{"SetMaterial", &LuaSimulation::SetMaterial},
			{"AddLayer"          , &LuaSimulation::AddLayer},
			{"SetRegionCircle"   , &LuaSimulation::SetRegionCircle},
			{"SetRegionEllipse"  , &LuaSimulation::SetRegionEllipse},
			{"SetRegionRectangle", &LuaSimulation::SetRegionRectangle},
			{"SetRegionPolygon"  , &LuaSimulation::SetRegionPolygon},
			{"RemoveLayerRegions", &LuaSimulation::RemoveLayerRegions},
			{"SetLayerThickness" , &LuaSimulation::SetLayerThickness},
			{"SetK"         , &LuaSimulation::SetK},
			{"SetFrequency" , &LuaSimulation::SetFrequency},
			{"OutputLayerDescription", &LuaSimulation::OutputLayerDescription},
			{"SolveAllLayers", &LuaSimulation::SolveAllLayers},
			{"GetPowerFlux", &LuaSimulation::GetPowerFlux},
			{"GetFields", &LuaSimulation::GetFields},
			{NULL, NULL}
		};
		
		luaL_newmetatable(L, className);
		lua_pushvalue(L, -1);
		lua_setfield(L, -2, "__index");
		luaL_setfuncs(L, SimulationLib, 0);
		lua_pop(L, 1);  /* pop new metatable */
    }
};

const char LuaSimulation::className[] = "S4r::Simulation";

extern "C" int luaopen_S4r(lua_State *L){
	static const luaL_Reg S4r_lib[] = {
		{"NewSimulation", &LuaSimulation::create},
		{NULL, NULL}
	};
	luaL_newlib(L, S4r_lib);
	
	LuaSimulation::Register(L);
	return 1;
}
