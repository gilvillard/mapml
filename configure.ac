# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([mapml],[0.0.1])
AC_CONFIG_SRCDIR([maple/conversion.c])
AM_MAINTAINER_MODE([enable])
AC_CANONICAL_TARGET


AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_HEADERS([maple/config.h])

DEFAULT_CXXFLAGS="-O3 -Wall -I/usr/local/include"

CXXFLAGS=${CXXFLAGS:-$DEFAULT_CXXFLAGS}

LT_INIT


AC_PROG_INSTALL

# Checks for programs.
AC_LANG([C++])
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AX_CXX_COMPILE_STDCXX([11],[noext],[mandatory])


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_C_VOLATILE

# Checks for library functions.
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor pow rint sqrt strtol])


CHECK_ATOMIC()

AC_ARG_WITH(gmp,
      AS_HELP_STRING([--with-gmp=@<:@=DIR@:>@], [GMP install directory]), [
         CPPFLAGS="$CXXFLAGS -I$withval/include"
         LDFLAGS="$LDFLAGS -L$withval/lib"
         mpir_lib_path="$withval/lib"
      ])
AC_CHECK_LIB(gmp, __gmp_get_memory_functions, , [AC_MSG_ERROR(
   [GMP version >= 4.2.0 needed, see http://gmplib.org])])





AC_CONFIG_FILES([Makefile maple/Makefile])


AC_OUTPUT

