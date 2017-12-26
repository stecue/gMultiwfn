AC_INIT([gMultiwfn], [3.4.1.99.20171225])
: ${FFLAGS=""}
: ${FCFLAGS=""}
AC_PROG_F77([ifort gfortran])
AC_PROG_FC([ifort gfortran])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([src/function.f90])
AC_CONFIG_MACRO_DIR([config])
AM_CONDITIONAL([IFORT],[test x$FC = xifort ])
AC_FC_SRCEXT(F)
AC_FC_SRCEXT(f90)
dnl AC_OPENMP
AX_LAPACK([],AC_MSG_ERROR([lapack-devel not found!]))
AM_INIT_AUTOMAKE(gMultiwfn, 3.4.1.99.20171225)
AC_CONFIG_FILES([
  Makefile
  src/Makefile
])
AC_OUTPUT
