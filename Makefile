# Set defaults
ALIB_EXT = a
SHLIB_EXT = so

# Set platform specific variables
include Makefile.plat

# The following must be defined already:
#   CC
#   CXX
#   OPTFLAGS
#   LUA_INC
#   LUA_LIB
#   LUA_MODULE_LIB
#   LA_LIBS
#   OBJDIR
#   SHLIB_EXT (so, dll)
#   SHLIB_FLAGS

CPPFLAGS += -IS4 -IS4/RNP -IS4/kiss_fft

#### Compilation targets

S4_LIBNAME = libS4.$(ALIB_EXT)
S4_LUA_LIBNAME = libS4_lua.$(ALIB_EXT)
S4_LUA_MODNAME = S4v2.$(SHLIB_EXT)

all: $(OBJDIR)/$(S4_LIBNAME) $(OBJDIR)/$(S4_LUA_MODNAME)

objdir:
	mkdir -p $(OBJDIR)
	mkdir -p $(OBJDIR)/S4k
	
S4_LIBOBJS = \
	$(OBJDIR)/S4k/S4.o \
	$(OBJDIR)/S4k/rcwa.o \
	$(OBJDIR)/S4k/fmm_common.o \
	$(OBJDIR)/S4k/fmm_FFT.o \
	$(OBJDIR)/S4k/fmm_kottke.o \
	$(OBJDIR)/S4k/fmm_closed.o \
	$(OBJDIR)/S4k/fmm_PolBasisNV.o \
	$(OBJDIR)/S4k/fmm_PolBasisVL.o \
	$(OBJDIR)/S4k/fmm_PolBasisJones.o \
	$(OBJDIR)/S4k/fmm_experimental.o \
	$(OBJDIR)/S4k/fft_iface.o \
	$(OBJDIR)/S4k/pattern.o \
	$(OBJDIR)/S4k/intersection.o \
	$(OBJDIR)/S4k/predicates.o \
	$(OBJDIR)/S4k/numalloc.o \
	$(OBJDIR)/S4k/gsel.o \
	$(OBJDIR)/S4k/sort.o \
	$(OBJDIR)/S4k/kiss_fft.o \
	$(OBJDIR)/S4k/kiss_fftnd.o
	
ifndef LAPACK_LIB
  S4_LIBOBJS += $(OBJDIR)/S4k/Eigensystems.o
endif

$(OBJDIR)/$(S4_LIBNAME): objdir $(S4_LIBOBJS)
	$(AR) crvs $@ $(S4_LIBOBJS)

$(OBJDIR)/S4k/S4.o: S4/S4.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/rcwa.o: S4/rcwa.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/fmm_common.o: S4/fmm/fmm_common.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/fmm_FFT.o: S4/fmm/fmm_FFT.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/fmm_kottke.o: S4/fmm/fmm_kottke.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/fmm_closed.o: S4/fmm/fmm_closed.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/fmm_PolBasisNV.o: S4/fmm/fmm_PolBasisNV.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/fmm_PolBasisVL.o: S4/fmm/fmm_PolBasisVL.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/fmm_PolBasisJones.o: S4/fmm/fmm_PolBasisJones.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/fmm_experimental.o: S4/fmm/fmm_experimental.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/fft_iface.o: S4/fmm/fft_iface.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/pattern.o: S4/pattern/pattern.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/intersection.o: S4/pattern/intersection.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/predicates.o: S4/pattern/predicates.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/numalloc.o: S4/numalloc.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/gsel.o: S4/gsel.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/sort.o: S4/sort.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/kiss_fft.o: S4/kiss_fft/kiss_fft.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/kiss_fftnd.o: S4/kiss_fft/tools/kiss_fftnd.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@
$(OBJDIR)/S4k/Eigensystems.o: S4/RNP/Eigensystems.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@

#### Python extension

S4_pyext: objdir $(OBJDIR)/libS4.a
	echo "$(LIBS)" > $(OBJDIR)/tmp.txt
	sh gensetup.py.sh $(OBJDIR) $(OBJDIR)/$(S4_LIBNAME)
	python setup.py build

clean:
	rm -rf $(OBJDIR)

#### Lua extension
$(OBJDIR)/$(S4_LUA_LIBNAME): $(OBJDIR)/$(S4_LIBNAME) S4/ext_lua.c
	$(CC) -c $(LUA_INC) S4/ext_lua.c -o $(OBJDIR)/ext_lua.o
	$(AR) crvs $@ $(OBJDIR)/ext_lua.o
$(OBJDIR)/$(S4_LUA_MODNAME): $(OBJDIR)/$(S4_LIBNAME) $(OBJDIR)/$(S4_LUA_LIBNAME) S4/ext_lua.c
	$(CC) $(SHLIB_FLAGS) $(LUA_INC) S4/ext_lua.c -o $@ $(LUA_MODULE_LIB) -L$(OBJDIR) -lS4 $(LA_LIBS) -lstdc++
