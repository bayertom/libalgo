#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([libalgo], [1.0.0], [bayertom@natur.cuni.cz])
AM_INIT_AUTOMAKE([-Wall foreign subdir-objects])
AC_CONFIG_SRCDIR([src/libalgo/source/comparators/sortSamplesByVoronoiCellTARRatio.h])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([floor pow sqrt strchr])

AC_CONFIG_FILES([
        Makefile
        src/Makefile
])

AC_OUTPUT
