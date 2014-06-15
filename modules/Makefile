OPTFLAGS = -O3
LUA_INCLUDE = -I../lua-5.2.3/install/include
LUA_LIB = -L../lua-5.2.3/install/lib -llua
OBJ_DIR = build

all: FunctionSampler1D.so FunctionSampler2D.so
FunctionSampler1D.so: function_sampler_1d.c function_sampler_1d.h lua_function_sampler_1d.c
	gcc -c $(OPTFLAGS) -fpic -Wall -I. function_sampler_1d.c -o ../$(OBJ_DIR)/function_sampler_1d.o
	gcc $(OPTFLAGS) -shared -fpic -Wall $(LUA_INCLUDE) -o FunctionSampler1D.so ../$(OBJ_DIR)/function_sampler_1d.o lua_function_sampler_1d.c $(LUA_LIB)
FunctionSampler2D.so: function_sampler_2d.c function_sampler_2d.h lua_function_sampler_2d.c
	gcc -c $(OPTFLAGS) -fpic -Wall -I. function_sampler_2d.c -o ../$(OBJ_DIR)/function_sampler_2d.o
	gcc -c -O2 -fpic -Wall -I. predicates.c -o ../$(OBJ_DIR)/mod_predicates.o
	gcc $(OPTFLAGS) -shared -fpic -Wall $(LUA_INCLUDE) -o FunctionSampler2D.so ../$(OBJ_DIR)/function_sampler_2d.o ../$(OBJ_DIR)/mod_predicates.o lua_function_sampler_2d.c $(LUA_LIB)
clean:
	rm -f FunctionSampler2D.so FunctionSampler1D.so
