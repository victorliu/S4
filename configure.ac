#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT(S4, 1.0.0, vkl@stanford.edu, S4)
AM_INIT_AUTOMAKE(S4, 1.0.0)
AC_CONFIG_SRCDIR([src/rcwa.cpp])
AC_CONFIG_HEADERS([src/config.h])

S4_VERSION=$PACKAGE_VERSION
AC_SUBST(S4_VERSION)

# Get C compiler.
AC_MSG_CHECKING([for vendor cc to be used instead of gcc])
AC_MSG_RESULT()
AC_CHECK_PROGS(CC, cc xlc) # prefer vendor cc, to stay in sync with Fortran
AC_PROG_CC
AM_PROG_CC_C_O

AC_MSG_CHECKING([for vendor cxx to be used instead of g++])
AC_MSG_RESULT()
AC_CHECK_PROGS(CXX, CC xlC) # prefer vendor cc, to stay in sync with Fortran
AC_PROG_CXX

# Checks for programs.
AC_PROG_AWK
AC_PROG_MAKE_SET

# check how to transform the name of the installed program:
AC_ARG_PROGRAM

# Add lots of compiler warnings to check for if we are using gcc:
# (The variable $GCC is set to "yes" by AC_PROG_CC if we are using gcc.)
if test "$GCC" = "yes"; then
	CFLAGS="$CFLAGS -Wall -W -Wbad-function-cast -Wcast-qual -Wpointer-arith -Wcast-align -pedantic"
fi


##############################################################################
#                       Checks for libraries:
##############################################################################

##############################################################################
# Check for Lua library
sinclude(m4/ax_lua.m4)
AX_PROG_LUA([5.1], [5.2], [], [])
AX_LUA_HEADERS([], [
	AC_MSG_ERROR([Couldn't find the required Lua headers.])])
AX_LUA_LIBS([], [
	AC_MSG_ERROR([Couldn't find the required Lua libraries.])])
#AC_CHECK_LIB(lua, luaL_newstate, [], [
#	AC_MSG_ERROR([Couldn't find the required liblua library.])], [-lm])

##############################################################################
# Check for CHOLMOD library
AC_CHECK_LIB(cholmod, cholmod_start, [], [], [-lamd -lcolamd -lcamd -lccolamd])

##############################################################################
# Check for FFTW3 library
AC_CHECK_LIB(fftw3, fftw_malloc)

##############################################################################
# Checks for BLAS and LAPACK libraries:
sinclude(m4/ax_blas.m4)
sinclude(m4/ax_lapack.m4)
AX_BLAS([], [AC_MSG_ERROR([BLAS was not found!])])
AX_LAPACK([], [AC_MSG_ERROR([LAPACK was not found!])])
LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

#AC_ARG_ENABLE(blacs, AS_HELP_STRING([--with-blacs], [Blacs support]))
#AC_ARG_ENABLE(scalapack, AS_HELP_STRING([--with-scalapack], [ScaLAPACK support, add before --with-blacs]))
#ACX_BLACS
#ACX_SCALAPACK

##############################################################################
# Checks for math library:
AC_CHECK_LIB([m], [sqrt])
 If we have a Fortran compiler, link to its libraries, since these
 are probably needed to link with the BLAS/LAPACK libraries.
if test x != x"$F77"; then
	AC_F77_LIBRARY_LDFLAGS
	AC_F77_WRAPPERS
fi

##############################################################################
# Check for MPI library

AC_ARG_WITH(mpi, [AC_HELP_STRING([--with-mpi],[enable MPI parallelization])], with_mpi=$withval, with_mpi=no)

if test "x$with_mpi" = "xyes"; then
	sinclude(m4/ax_mpi.m4)
	AX_MPI([], [AC_MSG_ERROR(could not find mpi library for --with-mpi)])
        CC="$MPICC"
	LIBS="$MPILIBS $LIBS"

	# Make installed program mpb-mpi instead of mpb
	MPB_SUFFIX="${MPB_SUFFIX}_mpi"

	AC_DEFINE(HAVE_MPI,1,[Define if you have & link an MPI library.])
fi

AC_CHECK_LIB([gfortran], [_gfortran_abort])
AC_CHECK_LIB([pthread], [pthread_create])

##############################################################################
# Checks for header files.
AC_FUNC_ALLOCA
AC_CHECK_HEADERS([float.h inttypes.h limits.h malloc.h stddef.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
# 2012-01-23 vkl Fix ::malloc not found error
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memset pow sqrt strcasecmp strchr strdup])

AC_SUBST(S4_SUFFIX)

AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 examples/Makefile
                 src/Makefile])
AC_OUTPUT
