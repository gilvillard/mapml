
# 
# Copyright (C) 2024 Gilles Villard
#
# Created from https://github.com/fplll/fplll/blob/master/configure.ac
#
#   This file is part of mapml. mapml is free software: you
#   can redistribute it and/or modify it under the terms of the GNU Lesser
#   General Public License as published by the Free Software Foundation,
#   either version 2.1 of the License, or (at your option) any later version.
#
#   mapml is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with mapml. If not, see <http://www.gnu.org/licenses/>. 


# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([mapml],[0.0.1])
AC_CONFIG_SRCDIR([src/conversion.c])
AM_MAINTAINER_MODE([enable])
AC_CANONICAL_TARGET


AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_HEADERS([src/config.h])

DEFAULT_CPPFLAGS="-O3 -Wall -I/usr/local/include"

CPPFLAGS=${CPPFLAGS:-$DEFAULT_CPPFLAGS}

DEFAULT_LDFLAGS="-L/usr/local/lib"

LDFLAGS=${LDFLAGS:-$DEFAULT_LDFLAGS}

#Defaults 


AC_SUBST(MAPLE_INCDIR,"/Library/Frameworks/Maple.framework/Versions/2024/extern/include")
AC_SUBST(MAPLE_LIBDIR,"/Library/Frameworks/Maple.framework/Versions/2024/bin.APPLE_ARM64_MACOS")

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


# Check gmp
#----------

AC_ARG_WITH(gmp,
      AS_HELP_STRING([--with-gmp=@<:@=DIR@:>@], [GMP install directory]), [
         CPPFLAGS="$CPPFLAGS -I$withval/include"
         LDFLAGS="$LDFLAGS -L$withval/lib"
      ])

AC_CHECK_LIB(gmp, __gmp_get_memory_functions, , [AC_MSG_ERROR(
   [GMP version >= 4.2.0 needed, see http://gmplib.org])])


# Check flint
#------------

AC_ARG_WITH(flint,
      AS_HELP_STRING([--with-flint=@<:@=DIR@:>@], [Flint install directory]), [
         CPPFLAGS="$CPPFLAGS -I$withval/include"
         LDFLAGS="$LDFLAGS -L$withval/lib"
      ])

AC_CHECK_LIB(flint, __flint_get_memory_functions, , [AC_MSG_ERROR(
   [ERROR FLINT])])

LIBS="$LIBS -lm -lpthread"

# Check pml
#----------

#GV bien voir pour cette bibliothèque, à reprendre 

PML_INCDIR="/usr/local/include"
PML_LIBDIR="/usr/local/lib"

AC_ARG_WITH(pml,
      AS_HELP_STRING([--with-pml=@<:@=DIR@:>@], [PML install directory]), [
         PML_INCDIR="$withval/include"
         PML_LIBDIR="$withval"
      ])

CPPFLAGS="$CPPFLAGS -I${PML_INCDIR}"
LDFLAGS="$LDFLAGS -L${PML_LIBDIR}"

LIBS="$LIBS -lpml"


AC_CHECK_FILE([$PML_INCDIR/pml/nmod_extra.h],
                [AC_DEFINE([HAVE_PML], [1],
                [with nmod_extra.h])],
                [AC_MSG_ERROR([Sorry, PML nmod_extra.h location missing])],
                [])

AC_SUBST(PML_INCDIR)
AC_SUBST(PML_LIBDIR)


# Check maple
#------------

AC_ARG_WITH(maplec,
      AS_HELP_STRING([--with-maplec=@<:@=DIR@:>@], [Maple include]), [
         MAPLE_INCDIR="$withval/include"
      ])


CPPFLAGS="$CPPFLAGS -I${MAPLE_INCDIR}"

AC_CHECK_HEADER([maplec.h],
                [AC_DEFINE([HAVE_MAPLEC], [1],
                   [with maplec.h])],
                [AC_MSG_ERROR([Sorry, maplec.h location missing])])

# Note: no lib directory here 

AC_SUBST(MAPLE_INCDIR)

AC_ARG_WITH(maplelib,
      AS_HELP_STRING([--with-maplelib=@<:@=DIR@:>@], [Maple lib]), [
         MAPLE_LIBDIR="$withval"
      ])


LDFLAGS="$LDFLAGS -L${MAPLE_LIBDIR}"
LIBS="$LIBS -lmaplec"

#GV Does not seem to work, howto do ? 

#AC_CHECK_LIB(maplec, __EvalMapleProc, , [AC_MSG_ERROR(
 #  [ERROR MAPLELIB])])


AS_ECHO(["##--------------------------------------##"])

AC_CONFIG_FILES([Makefile src/Makefile maple/Makefile])


AC_OUTPUT

