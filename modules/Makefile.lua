OPTFLAGS = -O3
LUA_INCLUDE = -I/usr/include/lua5.2
OBJDIR = ../obj

all: FunctionSampler1D.so FunctionSampler2D.so
FunctionSampler1D.so: function_sampler_1d.c function_sampler_1d.h lua_function_sampler_1d.c
	gcc -c $(OPTFLAGS) -fpic -Wall -I. function_sampler_1d.c -o $(OBJDIR)/function_sampler_1d.o
	gcc $(OPTFLAGS) -shared -fpic -Wall $(LUA_INCLUDE) -o FunctionSampler1D.so $(OBJDIR)/function_sampler_1d.o lua_function_sampler_1d.c
FunctionSampler2D.so: function_sampler_2d.c function_sampler_2d.h lua_function_sampler_2d.c
	gcc -c $(OPTFLAGS) -fpic -Wall -I. function_sampler_2d.c -o $(OBJDIR)/function_sampler_2d.o
	gcc -c -O2 -fpic -Wall -I. predicates.c -o $(OBJDIR)/mod_predicates.o
	gcc $(OPTFLAGS) -shared -fpic -Wall $(LUA_INCLUDE) -o FunctionSampler2D.so $(OBJDIR)/function_sampler_2d.o $(OBJDIR)/mod_predicates.o lua_function_sampler_2d.c
clean:
	rm -f *.o FunctionSampler2D.so FunctionSampler1D.so