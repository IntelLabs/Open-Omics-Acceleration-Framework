## Process this file with automake to produce Makefile.in

# $Id: Makefile.am,v 1.31 2014/07/04 01:27:43 mp Exp $
# 
# AutoGrid 
# 
# Copyright (C) 1989-2009, TSRI
# All Rights Reserved.
# 
# AutoGrid is a Trade Mark of The Scripps Research Institute.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#
# Makefile.am to build AutoGrid 4
#

# Autotools

EXTRA_DIST = Tests

dist-hook:
	rm -rf `find $(distdir)/Tests -name CVS -o -name '*.glg'` 

bin_PROGRAMS = autogrid4 

autogrid4_SOURCES = main.cpp \
       check_size.cpp \
       setflags.cpp \
       timesys.cpp \
       timesyshms.cpp \
       printhms.cpp \
       prHMSfixed.cpp \
       printdate.cpp \
       strindex.cpp \
       banner.cpp \
       gpfparser.cpp \
       parsetypes.cpp \
       atom_parameter_manager.cpp \
       read_parameter_library.cpp \
       parse_param_line.cpp \
       distdepdiel.cpp \
       mingw_getrusage.cpp \
       calc_vina_potential.cpp \
       stop.cpp \
            atom_parameter_manager.h \
            autoglobal.h \
            autogrid.h \
            constants.h \
            default_parameters.h \
            distdepdiel.h \
            gpftoken.h \
            grid.h \
	    mingw_sysconf.h \
	    mingw_sys_times.h \
            parameters.h \
            parm.h \
            parse_param_line.h \
            partokens.h \
            printdate.h \
            printhms.h \
            prototypes.h \
            ranlib.h \
            read_parameter_library.h \
            times.h \
            timesys.h \
            timesyshms.h \
            typedefs.h \
            util.h \
            calc_vina_potential.h \
            gridmap.h \
            version.h 


# Define the AutoGrid 4 source code files:
# See GNU Automake documentation, section 10.4 re BUILT_SOURCES;
# "default_parameters.h" is generated using the rules below
# which invoke "paramdat2h.csh" to convert the AD4.1_bound.dat
# parameter library file into a C++ header file.
# This same step runs "copyfilesfromAD" to copy common source files 
# from the autodock source directory to the autogrid build directory.
BUILT_SOURCES = default_parameters.h


AM_CFLAGS = -Wall -O3 # All warnings, Agressive optimization

AM_CPPFLAGS = -Wall -O3 -I $(srcdir) -I$(srcdir)/../autodock # All warnings, Agressive optimization
AM_CPPFLAGS += -Wno-write-strings # increase compatibility with pre 4.x compilers for string handling
#  M Pique - add USE_8A_NBCUTOFF
AM_CPPFLAGS += -DUSE_8A_NBCUTOFF
#  M Pique - add USE_DOUBLE to define Real as double instead of float, see typedefs.h
AM_CPPFLAGS += -DUSE_DOUBLE

# autogrid4_CXXFLAGS = -Wall -O3 -ftree-vectorize # All warnings, Agressive optimization, GCC vectorization in Mac OS X

autogrid4_CPPFLAGS=$(AM_CPPFLAGS)

# LINKOPT = $(CSTD) $(OPT) -Wl,--stack=0x2000000 # Cygwin, 32MB stacksize

if WANT_BIGSTACK
    autogrid4_LDFLAGS = -Wl,--stack,8388608 # Cygwin, 32MB stacksize
endif


# DBUG = -DNDEBUG # No debugging and no assert code.
# DBUG = # Use assert code.
# DBUG = -g # dbx.
# DBUG = -g -DDEBUG # dbx + DEBUG-specific code.
# DBUG = -g3 # dbx + optimization.
# DBUG = -g3 -DDEBUG # dbx + optimization, + DEBUG-specific code.
# DBUG = -DDEBUG # Just DEBUG-specific code.
# DBUG = -DDEBUG2 # Just DEBUG2-specific code for tracking prop.selection.
# DBUG = -DDEBUG3 # Just DEBUG3-specific code for print age of individuals.

check : test

test : autogrid4 $(srcdir)/Tests/test_autogrid4.py
#if HAVE_PYTHON
	@echo " "
	@echo "Running the Python test_autogrid4.py unit-tests on autogrid4"
	@echo " "
# for technical reasons, autogrid must be run in the directory where the
# test files reside because it writes its maps into its current directory
# and expects to find AD4.1_bound.dat in that same directory.
# We create a run-time Testruns directory for this purpose.
# The name is different to allow "configure;make;make check" in the source dir.
	outputdir=`pwd` ;\
	umask 0 ;\
	test -d $$outputdir/Testruns && chmod -fR a+w $$outputdir/Testruns;\
	test -d $$outputdir/Testruns || mkdir $$outputdir/Testruns ;\
	cp $(srcdir)/Tests/AD4.1_bound.dat   $$outputdir/Testruns;\
	cp $(srcdir)/Tests/*.py   $$outputdir/Testruns;\
	cp $(srcdir)/Tests/*.gpf  $$outputdir/Testruns;\
	cp $(srcdir)/Tests/*.pdbqt  $$outputdir/Testruns;\
	cp $(srcdir)/Tests/*.map*  $$outputdir/Testruns;\
	chmod -fR a+w $$outputdir/Testruns;\
	cd $$outputdir/Testruns;\
	python ./test_autogrid4.py -d . -e ../autogrid4$(EXEEXT) -o .
#else
#	@echo " "
#	@echo "Python version 2.1 or higher needed for 'make test', but not found."
#	@echo "Skipping AutoGrid4 unit tests."
#	@echo " "
#endif


main.o : $(srcdir)/../autodock/autocomm.h

default_parameters.h : $(srcdir)/../autodock/paramdat2h.csh  $(srcdir)/../autodock/AD4_parameters.dat $(srcdir)/../autodock/AD4.1_bound.dat
	rm -f $@
	csh $(srcdir)/../autodock/paramdat2h.csh $(srcdir)/../autodock/AD4_parameters.dat  $(srcdir)/../autodock/AD4.1_bound.dat > $@
	sh $(srcdir)/copyfilesfromAD

# 2007-04-30 MP & GMM: this line is necessary to be able to compile read_parameter_library.o
read_parameter_library.o : default_parameters.h

# see Autoconf manual 4.7.3 (p. 25) support for VPATH
# enabling multi-platform builds
srcdir = @srcdir@
VPATH  = @srcdir@

# see Autoconf manual for how to extend 'make clean' to directories:
clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	sh $(srcdir)/copyfilesfromAD -r
	-rm -f default_parameters.h
	-rm -rf Testruns
#
# EOF
#
