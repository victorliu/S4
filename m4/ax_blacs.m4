dnl NOT available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_blacs.html
dnl
AC_DEFUN([ACX_BLACS], [
acx_blacs_ok=no

dnl We cannot use BLACS if MPI is not found
if test "x$acx_mpi_ok" != xyes; then
  acx_blacs_ok=nompi
else

dnl Get fortran linker name of BLACS function to check for.
dnl if not compiling with fortran, convert the names
m4_if(_AC_LANG, Fortran, [blacs_pinfo=blacs_pinfo], [AC_FC_FUNC(blacs_pinfo)])

dnl Check if the library was given in the command line
if test $acx_blacs_ok = no; then
  AC_ARG_WITH(blacs, [AS_HELP_STRING([--with-blacs=<lib>], [use BLACS library <lib>])])
  case $with_blacs in
    yes | "") ;;
    no) acx_blacs_ok=disable ;;
    -* | */* | *.a | *.so | *.so.* | *.o) LIBS_BLACS="$with_blacs" ;;
    *) LIBS_BLACS="-l$with_blacs" ;;
  esac
fi

dnl Backup LIBS 
acx_blacs_save_LIBS="$LIBS"
LIBS="$LIBS_BLACS $LIBS_LAPACK $LIBS_BLAS $LIBS $FLIBS"

dnl First, check LIBS_BLACS environment variable
if test $acx_blacs_ok = no; then
  AC_MSG_CHECKING([for $blacs_pinfo in $LIBS_BLACS])
  AC_TRY_LINK_FUNC($blacs_pinfo, [acx_blacs_ok=yes], [])
  if test $acx_blacs_ok = no; then
    AC_MSG_RESULT([$acx_blacs_ok ($LIBS_BLACS)])
  else
    AC_MSG_RESULT([$acx_blacs_ok ($LIBS_BLACS)])
  fi
fi

dnl Generic BLACS library?
for blacs in blacs blacs-openmpi; do
  if test x"$blacs" = xblacs-openmpi; then       
    blacsinit="blacsF77init-openmpi"
  else
    blacsinit="blacsF77init"
  fi
  if test $acx_blacs_ok = no; then
    AC_CHECK_LIB($blacs -l$blacsinit -l$blacs, $blacs_pinfo,
      [acx_blacs_ok=yes; LIBS_BLACS="$LIBS_BLACS -l$blacs -l$blacsinit -l$blacs"], [], [$FLIBS])
  fi
done

AC_SUBST(LIBS_BLACS)
LIBS="$acx_blacs_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blacs_ok" = xyes; then
  AC_DEFINE(HAVE_BLACS,1,[Defined if you have BLACS library.])
  $1
else
  AC_MSG_WARN([Could not find Blacs library (required for Scalapack). 
               *** Will compile without Scalapack support])
  $2
fi
fi
])dnl ACX_BLACS
