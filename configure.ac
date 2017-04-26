AC_INIT([gMultiwfn], [3.3.9])
AC_CONFIG_SRCDIR([src/function.f90])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config])
AM_INIT_AUTOMAKE(gMultiwfn, 3.3.9)
AC_PROG_F77([gfortran])
AC_PROG_FC([gfortran])
AX_LAPACK([],AC_MSG_ERROR([lapack-devel not found!]))
AC_CONFIG_FILES([
  Makefile
  src/Makefile
])
AC_OUTPUT
