# Copyright (c) XXX
# This file is part of XXX (and probably stolen here and there)

#  ========LICENCE========
# This file is part of the library LinBox.
#
# LinBox is free software: you can redistribute it and/or modify
# it under the terms of the  GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========
#/

#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.


AC_PREREQ([2.72])
AC_CONFIG_MACRO_DIR([macros])
AC_INIT([mapml],[0.0.1])

AC_CONFIG_SRCDIR([maple/dum.c])

#Recompile for any modification in autotool chain 
#AM_MAINTAINER_MODE(enable)

AC_CONFIG_HEADERS([maple/config.h])

AM_INIT_AUTOMAKE([1.10 -Wall -Wno-portability foreign])
#AX_PREFIX_CONFIG_H(config.h, __MAPML)

DEFAULT_CFLAGS="-O3 -Wall"

CXXFLAGS=${CXXFLAGS:-$DEFAULT_CFLAGS}

LT_INIT 

AC_PROG_INSTALL

AC_LANG([C++])
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET

#AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

AC_CHECK_INCLUDES_DEFAULT

AC_PROG_EGREP

AC_CHECK_HEADERS([float.h limits.h stddef.h stdlib.h string.h sys/time.h stdint.h pthread.h])

#
AC_PROG_EGREP
AC_PROG_SED


AC_ARG_WITH(gmp,
   [  --with-gmp=DIR          GMP install directory ], [
      if test "$withval" != "" ; then
      CXXFLAGS="$CXXFLAGS -I$withval/include"
       LDFLAGS="$LDFLAGS -L$withval/lib -Wl,-rpath,$withval/lib"
       gmp_lib_path="$withval/lib"
       fi
   ])


AC_ARG_WITH(mpfr,
   [  --with-mpfr=DIR          MPFR install directory ], [
     if test "$withval" != "" ; then
      CXXFLAGS="$CXXFLAGS -I$withval/include"
       LDFLAGS="$LDFLAGS -L$withval/lib -Wl,-rpath,$withval/lib"
       mpfr_lib_path="$withval/lib"
      fi
   ])


#################################################
AS_BOX([          START MAPML CONFIG          ])
#################################################



# Set OPTIM_FLAGS, DEBUG_FLAGS depending on compiler and command line arguments
SET_FLAGS

# Append -march=native or -mcpu=native (if recognized by the compiler) to
# OPTIM_FLAGS if not present in CXXFLAGS and not cross-compiling and
# --without-archnative is not set
ARCH_FLAGS

# Append -mfpmath=sse to OPTIM_FLAGS on i386 and i686 architecture with SSE
FPMATH_FLAGS


AS_ECHO([---------------------------------------])

AC_CHECK_LIB(gmp, __gmpz_init, , [AC_MSG_ERROR([GNU MP not found, see http://gmplib.org])])


AC_CHECK_LIB(mpfr, mpfr_add, , [AC_MSG_ERROR([MPFR not found, see http://www.mpfr.org])])

#GIV_CHECK_GMP(40)


AS_ECHO([---------------------------------------])

#LB_CHECK_MPFR


AC_CONFIG_FILES([Makefile maple/Makefile])

AC_OUTPUT












