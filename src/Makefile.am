# S4: Stanford Stratified Structure Solver
# Copyright (C) 2011, Victor Liu
#
# This file may be used without restriction.  It is in the public
# domain, and is NOT restricted by the terms of any GNU license.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# Victor Liu can be contacted at vkl@stanford.edu.

##############################################################################

bin_PROGRAMS = S4@S4_SUFFIX@

S4@S4_SUFFIX@_SOURCES = \
  S4.cpp \
  rcwa.cpp \
  fmm/fmm_common.cpp \
  fmm/fmm_FFT.cpp \
  fmm/fmm_kottke.cpp \
  fmm/fmm_closed.cpp \
  fmm/fmm_PolBasisNV.cpp \
  fmm/fmm_PolBasisVL.cpp \
  fmm/fmm_PolBasisJones.cpp \
  fmm/fmm_experimental.cpp \
  fmm/fft_iface.cpp \
  pattern/pattern.c \
  pattern/intersection.c \
  pattern/predicates.c \
  numalloc.c \
  gsel.c \
  sort.c \
  kiss_fft/kiss_fft.c \
  kiss_fft/tools/kiss_fftnd.c \
  SpectrumSampler.c \
  cubature.c \
  Interpolator.c \
  convert.c \
  main_lua.c

noinst_HEADERS = \
	RNP \
	Interpolator.h \
	S4.h \
	SpectrumSampler.h \
	convert.h \
	cubature.h \
	debug.h \
	gsel.h \
	numalloc.h \
	rcwa.h \
	sort.h \
	pattern/pattern.h \
	pattern/intersection.h \
	fmm/fmm.h \
	fmm/fft_iface.h \
	kiss_fft/kiss_fft.h \
	kiss_fft/_kiss_fft_guts.h \
	kiss_fft/tools/kiss_fftnd.h


S4@S4_SUFFIX@_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/RNP -I$(top_srcdir)/src/kiss_fft $(LUA_INCLUDE)
S4@S4_SUFFIX@_LDADD = $(LUA_LIB)

# what is printed out when invoking your program with --version:
VERSION_STRING = "S4@S4_SUFFIX@ @VERSION@, Copyright (C) 2011 by Stanford University"
