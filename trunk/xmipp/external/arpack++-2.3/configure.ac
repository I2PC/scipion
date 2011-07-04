#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.58)
AC_INIT(arpack++, 2.3, christophe.prudhomme@ujf-grenoble.fr)
AM_INIT_AUTOMAKE(arpack++, 2.3)

AC_CONFIG_SRCDIR([configure.ac])
AM_CONFIG_HEADER([config.h])


# Checks for programs.
AC_PROG_CXX


# Checks for libraries.
AC_CHECK_LIB([m], [log], [all_libs="-lm"],[ echo "you need libm.";exit 1; ])
dnl AC_CHECK_LIB([gfortran], [s_copy], [all_libs="-lg2c"],[ echo "you need g2c.";exit 1; ])
AC_CHECK_LIB([blas], [daxpy_], [all_libs="-lblas ${all_libs}"],[echo "blas might be missing.";all_libs="-lblas ${all_libs}"])
AC_CHECK_LIB([lapack], [dgemm_],[all_libs="-llapack ${all_libs}"],[echo "lapack might be missing.";all_libs="-llapack ${all_libs}"], [-lblas])
AC_CHECK_LIB([arpack], [dnaupd_], [all_libs="-larpack ${all_libs}"],[ echo "arpack might be missing.";all_libs="-larpack ${all_libs}"], [-llapack -lblas])

AC_SUBST(all_libs)
LIBS="${all_libs} ${LIBS}"

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE

AC_PROG_LIBTOOL

# Checks for library functions.
AC_CHECK_FUNCS([sqrt strchr])

AC_CONFIG_FILES([ Makefile include/Makefile src/Makefile
                  examples/Makefile
                  examples/matrices/Makefile
                  examples/matrices/sym/Makefile
                  examples/matrices/nonsym/Makefile
                  examples/matrices/complex/Makefile
                  examples/matprod/Makefile
                  examples/matprod/sym/Makefile
                  examples/matprod/nonsym/Makefile
                  examples/matprod/complex/Makefile
                  examples/reverse/Makefile
                  examples/reverse/sym/Makefile
                  examples/reverse/nonsym/Makefile
                  examples/reverse/complex/Makefile
                  examples/band/Makefile
                  examples/band/sym/Makefile
                  examples/band/nonsym/Makefile
                  examples/band/complex/Makefile
                  examples/dense/Makefile
                  examples/dense/sym/Makefile
                  examples/dense/nonsym/Makefile
                  examples/dense/complex/Makefile
                  examples/product/Makefile
                  examples/product/simple/Makefile
                  examples/product/sym/Makefile
                  examples/product/nonsym/Makefile
                  examples/product/complex/Makefile
])
# need superlu
#                  examples/areig/Makefile
#                  examples/harwell/Makefile
#                  examples/harwell/sym/Makefile
#                  examples/harwell/nonsym/Makefile
#                  examples/harwell/complex/Makefile
#])

AC_OUTPUT
