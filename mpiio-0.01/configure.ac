#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(mpiio, 0.01, g.marras@cineca.it)
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR([partition.c])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AC_ENABLE_DISABLE
AC_PROG_RANLIB

# Checks for libraries.

# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

pwd=$(pwd)
AC_SUBST(pwd)
AC_PREFIX_DEFAULT([../..])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
