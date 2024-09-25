#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# $Id: configure.ac,v 1.92 2014/07/19 01:12:41 mp Exp $
AC_PREREQ(2.62)

# for release replace x with appropriate version (4.2.4 was 2010-Jul-7)
#     release 4.2.5 was 2012-Nov-2
#     release 4.2.5.1 was 2012-Dec-13
#     release 4.2.6 was 2014-Jun
AC_INIT([autodock], [4.2.6], [autodock@scripps.edu])
# the m4_esyscmd following creates a string such as VERSION='4.2.5.x.2014-12-25'
# indicating when the "autoreconf" was most recently run.
# This date is also inserted into PACKAGE_VERSION and PACKAGE_STRING
#AC_INIT([autodock], m4_esyscmd([date +4.2.6.%Y-%m-%d|tr -d '\n']), [autodock@scripps.edu])
AC_CONFIG_SRCDIR([analysis.cc])

AM_INIT_AUTOMAKE()
#AM_INIT_AUTOMAKE([foreign])
#AM_INIT_AUTOMAKE([-Wall -Werror foreign])


# Checks for programs.
AC_PROG_AWK

# don't let AC_PROG_CXX force "-g -O2"
saved_CXXFLAGS=${CXXFLAGS}
AC_PROG_CXX
CXXFLAGS=${saved_CXXFLAGS}

AC_PROG_INSTALL
AC_PROG_RANLIB

AC_CANONICAL_HOST

# Checks for header files.
AC_HEADER_STDC
#AC_CHECK_HEADERS([fcntl.h limits.h stdlib.h string.h sys/param.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_MALLOC
# gethostname does not seem available in WIN32 MinGW so we include a local equiv in mingw_sysconf.h
# getpid has a WIN32 equivalent GetProcessId but WIN32 MinGW supplies getpid
# fchmod and constants S_IRGRP, S_IROTH not available in WIN32 MinGW so block using (output_state.cc)
# setlinebuf not available in WIN32 MinGW so we block using it (setflags.cc)
# sysconf not available in WIN32 MinGW so we include a local subset in mingw_sysconf.h
# times   not available in WIN32 MinGW so we include a local equiv in mingw_sys_times.h
# getrusage not available in WIN32 MinGW so we include a local equiv in mingw_getrusage.cc
AC_CHECK_FUNCS([gethostname getpid GetProcessId fchmod setlinebuf sysconf times getrusage pow sqrt strncasecmp])

case $host_os in
  *cygwin* ) CYGWIN=yes;;
         * ) CYGWIN=no;;
esac
AM_CONDITIONAL([WANT_BIGSTACK], [test "$CYGWIN" = yes])

# M Pique: Open MP detection taken from Joachim Wuttke at
# http://stackoverflow.com/questions/9037434/enable-openmp-checks-and-flags-in-autoconf-and-automake
# who suggests also:   AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS")
AC_LANG(C++)
AC_PROG_CXX
AC_OPENMP  # requires autoconf 2.62 or newer, sets OPENMP_CXXFLAGS if successful

AM_CONDITIONAL([WANT_OMP], [test "x$OPENMP_CXXFLAGS" != "x"])


# See "automake" manual, version 1.10, section 10.5 (Python)
#AM_PATH_PYTHON([2.1],, [:])
#AM_CONDITIONAL([HAVE_PYTHON], [test "$PYTHON" != :])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
