dnl linbox miscellaneous functonnnalities
dnl Copyright (c) the LinBox group
dnl This file is part of LinBox

 dnl ========LICENCE========
 dnl This file is part of the library LinBox.
 dnl 
 dnl LinBox is free software: you can redistribute it and/or modify
 dnl it under the terms of the  GNU Lesser General Public
 dnl License as published by the Free Software Foundation; either
 dnl version 2.1 of the License, or (at your option) any later version.
 dnl 
 dnl This library is distributed in the hope that it will be useful,
 dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
 dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 dnl Lesser General Public License for more details.
 dnl 
 dnl You should have received a copy of the GNU Lesser General Public
 dnl License along with this library; if not, write to the Free Software
 dnl Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 dnl ========LICENCE========
 dnl



AC_DEFUN([LB_MISC],
[

AC_ARG_WITH(default,
[AC_HELP_STRING([--with-default=<path>], [Add <path> to the default path for external package
  		        checking. Set as default with /usr and /usr/local.
])],
	    [if test "$withval" = yes ; then
			echo "Default path = /usr /usr/local "
			DEFAULT_CHECKING_PATH="/usr /usr/local "
	      else
			echo "Default path = $withval /usr /usr/local "
			DEFAULT_CHECKING_PATH="$withval /usr /usr/local "
	     fi
	     ],
	     [
		echo "Default path = /usr /usr/local "
		DEFAULT_CHECKING_PATH="/usr /usr/local "
             ])


AC_ARG_WITH(all,
[AC_HELP_STRING([--with-all=<path>|yes|no], [Use all external packages. If the argument is no,
  	      		   you not sure that all libraries are reachable with
			   the default path. If the argument is yes or <empty>,
			   that means that all libraries are reachable with the
			   default path. Otherwise add <path> to default path
			   and enable all external packages.
])],
	    [if test "$withval" = yes ; then
			check_all="yes"
			echo "Checking all external packages in ${DEFAULT_CHECKING_PATH}"

	      elif test "$withval" != no ; then
			check_all="yes"
			DEFAULT_CHECKING_PATH="$withval ${DEFAULT_CHECKING_PATH}"
			echo "Checking all external packages in ${DEFAULT_CHECKING_PATH}"
	     fi
	     ],
	     [])

if test -n "$check_all"; then

	GMP_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	MPFR_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	FLINT_HOME_PATH="${DEFAULT_CHECKING_PATH}"
fi


])
