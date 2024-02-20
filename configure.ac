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

DEFAULT_CPPFLAGS="-O3 -Wall -I/usr/local/include"

CPPFLAGS=${CPPFLAGS:-$DEFAULT_CPPFLAGS}

#Defaults 


AC_SUBST(MAPLE_INCDIR,"/Library/Frameworks/Maple.framework/Versions/2023/extern/include")
AC_SUBST(MAPLE_LIBDIR,"/Library/Frameworks/Maple.framework/Versions/2023/bin.APPLE_ARM64_MACOS")

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

#################################################
AS_BOX([          mapml config               ])
#################################################

AC_ARG_WITH(gmp,
      AS_HELP_STRING([--with-gmp=@<:@=DIR@:>@], [GMP install directory]), [
         CPPFLAGS="$CPPFLAGS -I$withval/include"
         LDFLAGS="$LDFLAGS -L$withval/lib"
      ])

AC_CHECK_LIB(gmp, __gmp_get_memory_functions, , [AC_MSG_ERROR(
   [GMP version >= 4.2.0 needed, see http://gmplib.org])])


AC_ARG_WITH(maplec,
      AS_HELP_STRING([--with-maplec=@<:@=DIR@:>@], [Maple include]), [
         MAPLE_INCDIR="$withval"
      ])


CPPFLAGS="$CPPFLAGS -I${MAPLE_INCDIR}"


AC_CHECK_HEADER([maplec.h],
                [AC_DEFINE([HAVE], [1],
                   [with maplec.h])],
                [AC_MSG_ERROR([Sorry, maplec.h location missing])])

AC_SUBST(MAPLE_INCDIR)

AS_ECHO(["##--------------------------------------##"])

AC_CONFIG_FILES([Makefile maple/Makefile])


AC_OUTPUT

