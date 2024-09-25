## Process this file with automake to produce Makefile.in

# $Id: Makefile.am,v 1.50 2014/07/08 19:36:27 mp Exp $
# 
# AutoDock  
# 
# Copyright (C) 1989-2007,  Garrett M. Morris,   TSRI
# All Rights Reserved.
# 
# AutoDock is a Trade Mark of The Scripps Research Institute.
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
# Makefile.am to build AutoDock 4
#

# Autotools

EXTRA_DIST = AD4_PARM99.dat AD4_parameters.dat AD4.1_bound.dat AD4.1_extended.dat AD4.1_compact.dat paramdat2h.csh Tests

dist-hook:
	rm -rf `find $(distdir)/Tests -name CVS -o -name '*.dlg'` 

if WANT_OMP
bin_PROGRAMS = autodock4   autodock4.omp
else
bin_PROGRAMS = autodock4   
endif


CLEANFILES = default_parameters.h

noinst_LIBRARIES = 

# Define the AutoDock 4 source code files:

autodock4_SOURCES = \
    main.cc

autodock4_SOURCES += \
    alea.cc \
    analysis.cc \
	atom_parameter_manager.cc \
    banner.cc \
    bestpdb.cc \
	calculateEnergies.cc \
    call_glss.cc \
    call_ls.cc \
    center_ligand.cc \
    changeState.cc \
    check_header_float.cc \
    check_header_int.cc \
    check_header_line.cc \
    cluster_analysis.cc \
    clmode.cc \
    cnv_state_to_coords.cc \
    com.cc \
    conformation_sampler.cc \
    dimLibrary.cc\
	distdepdiel.cc \
    stateLibrary.cc \
    readfield.cc \
    readmap.cc \
    readPDBQT.cc \
	read_parameter_library.cc \
    eval.cc \
    gencau.cc \
    getrms.cc \
    get_atom_type.cc \
    getInitialState.cc \
    getpdbcrds.cc \
    gs.cc \
    initautodock.cc \
    input_state.cc \
    investigate.cc \
    linpack.cc \
    ls.cc \
    mapping.cc \
    mingw_getrusage.cc \
    minmeanmax.cc \
    mkNewState.cc \
    mkTorTree.cc \
    mkRandomState.cc \
    nonbonds.cc \
    openfile.cc \
    output_state.cc \
    parse_dpf_line.cc \
    parse_param_line.cc \
    parse_PDBQT_line.cc \
    parse_trj_line.cc \
	parsetypes.cc \
    print_2x.cc \
    print_atomic_energies.cc \
    print_avsfld.cc \
    writePDBQT.cc \
    print_rem.cc \
    printdate.cc \
    printEnergies.cc \
    printhms.cc \
    prClusterHist.cc \
    prInitialState.cc \
    prTorConList.cc \
    pso.cc \
    qmultiply.cc \
    qtransform.cc \
    quicksort.cc \
    ranlib.cc \
    rep.cc \
    scauchy.cc \
    setflags.cc \
    simanneal.cc \
    sort_enrg.cc \
    stack.cc \
    stop.cc \
    strindex.cc \
    success.cc \
    support.cc \
    swap.cc \
    threadlog.cc \
    timesys.cc \
    timesyshms.cc \
    torNorVec.cc \
    torsion.cc \
    trilinterp.cc \
    usage.cc \
    weedbonds.cc \
    warn_bad_file.cc \
    coliny.cc \
    eintcal.cc \
    eintcalPrint.cc \
    intnbtable.cc \
    nbe.cc \
		PDBQT_tokens.h \
		analysis.h \
		atom_parameter_manager.h \
		alea.h \
		autocomm.h \
		autoglobal.h \
		banner.h \
		bestpdb.h \
		calculateEnergies.h \
		call_glss.h \
		call_ls.h \
		center_ligand.h \
		changeState.h \
		check_header_float.h \
		check_header_int.h \
		check_header_line.h \
		clmode.h \
		cluster_analysis.h \
		cnv_state_to_coords.h \
		coliny.h \
		conformation_sampler.h \
		constants.h \
		default_parameters.h \
		dimLibrary.h \
		distdepdiel.h \
		dpftoken.h \
		eintcal.h \
		eintcalPrint.h \
		energy.h \
		eval.h \
		gencau.h \
		getInitialState.h \
		get_atom_type.h \
		getpdbcrds.h \
		getrms.h \
		globals.h \
		grid.h \
		gs.h \
		hybrids.h \
		initautodock.h \
		input_state.h \
		intnbtable.h \
		investigate.h \
		ls.h \
		main.h \
		mdist.h \
		mingw_sysconf.h \
		mingw_sys_times.h \
		mkNewState.h \
		mkRandomState.h \
		mkTorTree.h \
		nbe.h \
		nonbonds.h \
		openfile.h \
		output_state.h \
		parameters.h \
		parse_PDBQT_line.h \
		parse_dpf_line.h \
		parse_param_line.h \
		parse_trj_line.h \
		parsetypes.h \
		partokens.h \
		prClusterHist.h \
		prInitialState.h \
		prTorConList.h \
		printEnergies.h \
		print_2x.h \
		print_atomic_energies.h \
		print_avsfld.h \
		print_rem.h \
		printdate.h \
		printhms.h \
        pso.h \
		qmultiply.h \
		qtransform.h \
		quicksort.h \
		ranlib.h \
		readPDBQT.h \
		read_parameter_library.h \
		readfield.h \
		readmap.h \
		rep.h \
		rep_constants.h \
		setflags.h \
		simanneal.h \
		sort_enrg.h \
		stateLibrary.h \
		stack.h \
		stop.h \
		strindex.h \
		structs.h \
		success.h \
		support.h \
		swap.h \
		threadlog.h \
		timesys.h \
		timesyshms.h \
		torNorVec.h \
		torsion.h \
		trilinterp.h \
		trjtokens.h \
		typedefs.h \
		usage.h \
		version.h \
		warn_bad_file.h \
		weedbonds.h \
		writePDBQT.h
if WANT_OMP
autodock4_omp_SOURCES = \
     main.cc 
autodock4_omp_SOURCES += \
   alea.cc \
    analysis.cc \
	atom_parameter_manager.cc \
    banner.cc \
    bestpdb.cc \
	calculateEnergies.cc \
    call_glss.cc \
    call_ls.cc \
    center_ligand.cc \
    changeState.cc \
    check_header_float.cc \
    check_header_int.cc \
    check_header_line.cc \
    cluster_analysis.cc \
    clmode.cc \
    cnv_state_to_coords.cc \
    com.cc \
    conformation_sampler.cc \
    dimLibrary.cc\
	distdepdiel.cc \
    stateLibrary.cc \
    readfield.cc \
    readmap.cc \
    readPDBQT.cc \
	read_parameter_library.cc \
    eval.cc \
    gencau.cc \
    getrms.cc \
    get_atom_type.cc \
    getInitialState.cc \
    getpdbcrds.cc \
    gs.cc \
    initautodock.cc \
    input_state.cc \
    investigate.cc \
    linpack.cc \
    ls.cc \
    mapping.cc \
    mingw_getrusage.cc \
    minmeanmax.cc \
    mkNewState.cc \
    mkTorTree.cc \
    mkRandomState.cc \
    nonbonds.cc \
    openfile.cc \
    output_state.cc \
    parse_dpf_line.cc \
    parse_param_line.cc \
    parse_PDBQT_line.cc \
    parse_trj_line.cc \
	parsetypes.cc \
    print_2x.cc \
    print_atomic_energies.cc \
    print_avsfld.cc \
    writePDBQT.cc \
    print_rem.cc \
    printdate.cc \
    printEnergies.cc \
    printhms.cc \
    prClusterHist.cc \
    prInitialState.cc \
    prTorConList.cc \
    pso.cc \
    qmultiply.cc \
    qtransform.cc \
    quicksort.cc \
    ranlib.cc \
    rep.cc \
    scauchy.cc \
    setflags.cc \
    simanneal.cc \
    sort_enrg.cc \
    stack.cc \
    stop.cc \
    strindex.cc \
    success.cc \
    support.cc \
    swap.cc \
    threadlog.cc \
    timesys.cc \
    timesyshms.cc \
    torNorVec.cc \
    torsion.cc \
    trilinterp.cc \
    usage.cc \
    weedbonds.cc \
    warn_bad_file.cc \
    coliny.cc \
    eintcal.cc \
    eintcalPrint.cc \
    intnbtable.cc \
    nbe.cc \
		PDBQT_tokens.h \
		analysis.h \
		atom_parameter_manager.h \
		alea.h \
		autocomm.h \
		autoglobal.h \
		banner.h \
		bestpdb.h \
		calculateEnergies.h \
		call_glss.h \
		call_ls.h \
		center_ligand.h \
		changeState.h \
		check_header_float.h \
		check_header_int.h \
		check_header_line.h \
		clmode.h \
		cluster_analysis.h \
		cnv_state_to_coords.h \
		coliny.h \
		conformation_sampler.h \
		constants.h \
		default_parameters.h \
		dimLibrary.h \
		distdepdiel.h \
		dpftoken.h \
		eintcal.h \
		eintcalPrint.h \
		energy.h \
		eval.h \
		gencau.h \
		getInitialState.h \
		get_atom_type.h \
		getpdbcrds.h \
		getrms.h \
		globals.h \
		grid.h \
		gs.h \
		hybrids.h \
		initautodock.h \
		input_state.h \
		intnbtable.h \
		investigate.h \
		ls.h \
		main.h \
		mdist.h \
		mingw_sysconf.h \
		mingw_sys_times.h \
		mkNewState.h \
		mkRandomState.h \
		mkTorTree.h \
		nbe.h \
		nonbonds.h \
		openfile.h \
		output_state.h \
		parameters.h \
		parse_PDBQT_line.h \
		parse_dpf_line.h \
		parse_param_line.h \
		parse_trj_line.h \
		parsetypes.h \
		partokens.h \
		prClusterHist.h \
		prInitialState.h \
		prTorConList.h \
		printEnergies.h \
		print_2x.h \
		print_atomic_energies.h \
		print_avsfld.h \
		print_rem.h \
		printdate.h \
		printhms.h \
        pso.h \
		qmultiply.h \
		qtransform.h \
		quicksort.h \
		ranlib.h \
		readPDBQT.h \
		read_parameter_library.h \
		readfield.h \
		readmap.h \
		rep.h \
		rep_constants.h \
		setflags.h \
		simanneal.h \
		sort_enrg.h \
		stateLibrary.h \
		stack.h \
		stop.h \
		strindex.h \
		structs.h \
		success.h \
		support.h \
		swap.h \
		threadlog.h \
		timesys.h \
		timesyshms.h \
		torNorVec.h \
		torsion.h \
		trilinterp.h \
		trjtokens.h \
		typedefs.h \
		usage.h \
		version.h \
		warn_bad_file.h \
		weedbonds.h \
		writePDBQT.h

endif


BUILT_SOURCES = $(autodock4_SOURCES) default_parameters.h

AM_CXXFLAGS = -O3 # Agressive optimization
#AM_CXXFLAGS = -g # debug

    # if gnu c++ compiler
AM_CXXFLAGS += -Wall # All warnings
AM_CXXFLAGS += -Wno-write-strings # increase compatibility with pre 4.x compilers for string handling

if WANT_OMP  
# M Pique - I'm not sure which of these are really necessary
# could this be replaced by :  AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS")
autodock4_omp_CFLAGS = $(AM_CFLAGS)  -fopenmp  $(OPENMP_CXXFLAGS)
autodock4_omp_CPPFLAGS = $(AM_CPPFLAGS)  $(OPENMP_CXXFLAGS)
autodock4_omp_LDFLAGS = $(AM_LDFLAGS) -fopenmp  $(OPENMP_CXXFLAGS)

autodock4_omp_CXXFLAGS = $(AM_CXXFLAGS) $(autodock4_CXXFLAGS) -fopenmp  $(OPENMP_CXXFLAGS)
autodock4_omp_LDADD = $(autodock4_LDADD)
endif



#  M Pique - add ASSERTQUATOK if you want error checking in qmultiply.cc
#  M Pique - add USE_DOUBLE to define Real as double instead of float, see typedefs.h
#  M Pique - USE_INT_AS_FOURBYTELONG is for random number seeds only and is appropriate for
#     all post-2005 computers, which treat "int" as four bytes.
#  M Pique - add -g for debugging
#  M Pique - add -O3 for optimization
AM_CPPFLAGS = -DNOSQRT \
	      -DUSE_DOUBLE \
              -DUSE_8A_NBCUTOFF \
	      -DUSE_INT_AS_FOURBYTELONG \
	      -O3
# AM_LDFLAGS = -g # -g is MP debug


#autodock4_CXXFLAGS = -Wall -O3 # All warnings, Agressive optimization
#autodock4_CXXFLAGS = -Wall -O3 -ftree-vectorize # All warnings, Agressive optimization + GCC vectorization # MacOSX only?

autodock4_CPPFLAGS = $(AM_CPPFLAGS)
autodock4_LDFLAGS = $(AM_LDFLAGS)


autodock4_LDADD = $(COLINY_LDFLAGS) -lm

if WANT_BIGSTACK
    autodock4_LDFLAGS += -Wl,--stack,8388608 # Cygwin, big stacksize
endif


#
# If you need to use debugging or profiling, these should also be 
# modified appropriately:
# 	DBUG & PROF
#
# If you want to use the Coliny solver library, uncomment the following:
# COLINY = yes
#

# Abbreviations:
#
# Alpha   = Compaq/Digital Equipment Corp., Alpha
# Convex  = Convex, c2
# Cygwin  = Cygwin running on Microsoft Windows
# Darwin  = Darwin
# HP      = Hewlett Packard Precision Architecture, hppa
# Linux   = Any platform that runs GNU/Linux, Linux
# MacOSX  = Apple Mac OS X 10.0 & higher, MacOSX
# SGI     = Silicon Graphics Inc., sgi4D
# Sun     = Sun Microsystems, sun4

# CFLAGS = $(CSTD) $(OPT) -DUSE_8A_NBCUTOFF # SGI, HP, Alpha, Sun, Convex, Cygwin, Linux, MacOSX: Standard accuracy, but faster; no crossover in quaternion
# CFLAGS = $(CSTD) $(OPT) -DUSE_8A_NBCUTOFF # SGI, HP, Alpha, Sun, Convex, Linux, MacOSX: Standard accuracy, but faster
# CFLAGS = $(CSTD) $(OPT) -DUSE_8A_NBCUTOFF -DUSE_DOUBLE # SGI, HP, Alpha, Sun, Convex, Linux, MacOSX: Standard accuracy, but faster; also use Double precision throughout
# CFLAGS = $(CSTD) $(OPT) # SGI, HP, Alpha, Sun, Convex, Cygwin, Linux, MacOSX

# CSTD = $(DBUG) $(PROF) $(WARN) # SGI, Sun, Linux, MacOSX
# CSTD = $(DBUG) $(PROF) $(WARN) -DUSE_XCODE # Smaller memory footprint, good for Xcode
# CSTD = $(DBUG) $(PROF) $(WARN) # SGI, Sun, Linux, MacOSX
# CSTD = -std arm -verbose $(PROF) $(DBUG) $(WARN) # Alpha. sarah

# OPTLEVEL = -O3 -ftree-vectorize # Agressive optimization + GCC vectorization
# OPTLEVEL = -O3 -ffast-math # Agressive optimization, for Intel Itanium, ia64Linux2
# OPTLEVEL = -O3 # Agressive optimization
# OPTLEVEL = -O2 # High optimization
# OPTLEVEL = -O1 # Do optimizations that can be done quickly; default.  Recommended for unit testing
# OPTLEVEL = -O0 # Do not optimize



# LINKOPT = $(CSTD) $(OPT) # All platforms except Cygwin and Sun
# LINKOPT = $(CSTD) $(OPT) -fno-stack-limit # Cygwin, 32MB stacksize
# LINKOPT = $(CSTD) $(OPT) -Wl,--stack=0x4000000 # Cygwin, 64MB stacksize
# LINKOPT = $(CSTD) $(OPT) -L/opt/sfw/lib # Sun

# LINK = $(LINKOPT) # Linking flags

# LIB = -lm # for all platforms

# DBUG = # Use assert code
# DBUG = -DNDEBUG # No debugging and no assert code
# DBUG = -g # dbx, or Gnu gdb
# DBUG = -g -DDEBUG # dbx + DEBUG-specific code
# DBUG = -g3 # dbx + optimization
# DBUG = -g3 -DDEBUG # dbx + optimization, + DEBUG-specific code
# DBUG = -DDEBUG # Just DEBUG-specific code
# DBUG = -DDEBUG2 # Just DEBUG2-specific code for tracking prop.selection
# DBUG = -DDEBUG3 # Just DEBUG3-specific code for print age of individuals
# DBUG = -g -DDEBUG -DDEBUG2 -DDEBUG3 # Debug everything
# DBUG = -g3 -DDEBUG_QUAT # Use assert code, & assert quaternions are 4D-normalised
# DBUG = -g3 -DDEBUG_QUAT -DDEBUG_QUAT_PRINT # Use assert code, & assert quaternions are 4D-normalised, and print out
# DBUG = -g3 -DDEBUG_QUAT -DDEBUG_QUAT_PRINT -DDEBUG # Use assert code, & assert quaternions are 4D-normalised, and use DEBUG-specific code, and print out
# DBUG = -g3 -DDEBUG_MUTATION # dbx + optimization; Use assert code, & print out information anywhere involving mutation
# DBUG = -g3 -DDEBUG_MUTATION -DDEBUG # dbx + optimization; Use assert code, & print out information anywhere involving mutation + all debug statements
# DBUG = -g3 -DDEBUG_MUTATION -DDEBUG -DDEBUG3 # dbx + optimization; Use assert code, & print out information anywhere involving mutation + all debug statements incl. 3
# DBUG = -DASSERTQUATOK # Uses the macro "assertQuatOK()" to ensure the quaternion magnitude is near unity

# PROF = # No profiling
# PROF = -p # CC Profiling
# PROF = -pg # For GNU gcc & GNU g++ Profiling with GNU gprof

# WARN = # Default warning level
# WARN = -woff all # For no warnings
# WARN = -Wall # All warnings, gcc -- Recommended for developers
# WARN = -fullwarn -ansiE -ansiW # For full warnings during compilation

### ifeq ($(COLINY),yes)
###   COLINYLIB = libcoliny.a # Using Coliny
### else
#COLINY_LDADD = # Not using Coliny
### endif

#
# Makefile targets
#

autodock4.html :
	@echo " "
	@echo "Creating a colour-coded HTML report of the AutoDock source code"
	@echo " "
	enscript -E --color -Whtml --toc -p$@ *.h *.cc

check : test

if WANT_OMP
test : autodock4 autodock4.omp $(srcdir)/Tests/test_autodock4.py
else
test : autodock4 $(srcdir)/Tests/test_autodock4.py
endif
#if HAVE_PYTHON
	@echo " "
	@echo "Running the Python test_autodock4.py unit-tests on autodock4"
	@echo " "
# for technical reasons, autodock must be run in the directory where the
# test files reside because it writes its logs into its current directory
# and expects to find AD4.1_bound.dat in that same directory.
# We create a run-time Testruns directory for this purpose.
# The name is different from Tests to allow "configure;make;make check" in the source dir.
	outputdir=`pwd` ;\
	umask 0 ;\
	test -d $$outputdir/Testruns && chmod -fR a+w $$outputdir/Testruns;\
	test -d $$outputdir/Testruns || mkdir $$outputdir/Testruns ;\
	cp $(srcdir)/AD4.1_bound.dat   $$outputdir/Testruns;\
	cp $(srcdir)/AD4_parameters.dat   $$outputdir/Testruns;\
	cp $(srcdir)/Tests/*.py   $$outputdir/Testruns;\
	cp $(srcdir)/Tests/*.dpf  $$outputdir/Testruns;\
	cp $(srcdir)/Tests/*.pdbqt  $$outputdir/Testruns;\
	cp $(srcdir)/Tests/*.map*  $$outputdir/Testruns;\
	chmod -fR a+w $$outputdir/Testruns;
	-cd Testruns; rm -f test_*.dlg ; \
	python ./test_autodock4.py -d . -e ../autodock4$(EXEEXT) -o . ;
if WANT_OMP
	@echo " "
	@echo "Running the Python test_autodock4.py unit-tests on autodock4.omp"
	@echo " "
	-cd Testruns; rm -f test_*.dlg ; \
	python ./test_autodock4.py -d . -e ../autodock4.omp$(EXEEXT) -o . ;
endif

#else
#	@echo " "
#	@echo "Python version 2.1 or higher needed for 'make check', but not found."
#	@echo "Skipping AutoDock4 unit tests."
#	@echo " "
#endif

#$(COLINYLIB) : $(OBJS) $(OBJNOSQRT) $(OBJNOMINPT)
#	@echo " "
#	@echo Making the Coliny library
#	@echo " "
#	$(AR) $(ARFLAGS) $(COLINYLIB) $(?:.cc=.o)
#	$(RANLIB) $(RANLIBFLAGS) $(COLINYLIB)
###
### alternative from another CVS user
###
# coliny: $(ACRO)

### acro:
### 	@echo " "
### 	@echo Making the Coliny library
### 	@echo " "
### 	@(cd ../acro;\
### 	if [ ! -e Makefile ]; then\
### 	   ./setup configure;\
### 	fi;\
### 	$(MAKE))


#
# Object dependencies:
#

gs.o : gs.cc gs.h ranlib.h eval.h rep.h support.h writePDBQT.h
	$(CXXCOMPILE) -DCHECK_ISNAN -c $(srcdir)/gs.cc

# 2007-04-30 MP & GMM: this line is necessary to be able to compile read_parameter_library.o
read_parameter_library.o : default_parameters.h

default_parameters.h : AD4_parameters.dat AD4.1_bound.dat paramdat2h.csh
	rm -f $@
	csh $(srcdir)/paramdat2h.csh $(srcdir)/AD4_parameters.dat $(srcdir)/AD4.1_bound.dat > $@

#coliny.o : coliny.h
#	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(ACRO_INCLUDES) -c coliny.cc -o coliny.o

# see Autoconf manual 4.7.3 (p. 25) support for VPATH
# enabling multi-platform builds
srcdir = @srcdir@
VPATH  = @srcdir@

# see Autoconf manual for how to extend 'make clean' to directories:
clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-rm -rf Testruns
#
# EOF
#
