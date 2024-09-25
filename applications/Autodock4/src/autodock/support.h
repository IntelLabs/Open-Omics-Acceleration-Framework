/*

 $Id: support.h,v 1.24 2014/06/23 23:41:28 mp Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
 All Rights Reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

//  These are the class used to support the Representation classes.
//  Genotypes are the representations that the Global_Search class 
//  and its derivations acts on.  Local_Search and its children act
//  on the Phenotype classes.  Phenotypes are basically what results
//  from mapping the Genotype to the solution domain.  It has the 
//  fundamental characteristic of fitness.  We need to make sure that 
//  the const pointers are never used to indirectly change values!
//  We can factor Genotype and Phenotype into 
//  a common base class Chromosome
//  rsh 07/08/95

/*
** $Log: support.h,v $
** Revision 1.24  2014/06/23 23:41:28  mp
** if-def'ing out some constructor debugging code
**
** Revision 1.23  2014/06/12 01:44:08  mp
** General updates to accommodate OpenMP parallelization, mostly to pass
** logFile references and to make the timesys and random-number functions
** thread-safe.  New "threadlog" utility functions added to code.
**
** Revision 1.22  2011/03/08 04:18:37  mp
** Incorporation of Steffan Moeller patch set 20101104 and 2010114.
** Virtually all changes proposed are incorporated except for declaring "static"
** the AutoDock-unused functions in ranlib.cc/ranlib.h.
**
**  Modified Files:
** 	Makefile.am alea.cc alea.h call_cpso.cc check_header_int.cc
** 	configure.ac conformation_sampler.cc conformation_sampler.h
** 	distdepdiel.cc distdepdiel.h gs.cc hybrids.h intnbtable.cc
** 	intnbtable.h ls.cc mkTorTree.cc paramdat2h.csh
** 	parse_dpf_line.cc parse_param_line.cc qmultiply.cc qmultiply.h
** 	ranlib.cc ranlib.h readmap.cc rep.cc rep.h simanneal.cc
** 	stateLibrary.cc support.cc support.h torNorVec.cc typedefs.h
**
** Revision 1.21  2010/10/01 22:51:40  mp
** Applied patches 2010-09-29 from Steffan Moeller "patches_introducing_references"
**
** Revision 1.20  2010/08/27 00:05:08  mp
** Integration of Steffan Möller <steffen_moeller@gmx.de> "const" contribution
** from his patch files 4 August 2010, as adapted slightly by Michael Pique
**
** Note: patches applied so far only to AutoDock, not AutoGrid.
**
**  Modified Files:
**  	alea.cc alea.h analysis.cc analysis.h banner.cc banner.h
**  	calculateEnergies.cc calculateEnergies.h call_cpso.cc
**  	call_glss.cc changeState.cc changeState.h
**  	check_header_float.cc check_header_float.h check_header_int.cc
**  	check_header_int.h check_header_line.cc check_header_line.h
**  	clmode.cc clmode.h cluster_analysis.cc cluster_analysis.h
**  	cnv_state_to_coords.cc cnv_state_to_coords.h coliny.cc
**  	coliny.h com.cc configure.ac dimLibrary.cc dimLibrary.h
**  	distdepdiel.cc distdepdiel.h eintcal.cc eintcal.h
**  	eintcalPrint.h eval.cc eval.h gencau.cc getInitialState.cc
**  	getInitialState.h get_atom_type.cc getpdbcrds.cc getpdbcrds.h
**  	getrms.cc getrms.h gs.cc gs.h hybrids.h initautodock.cc
**  	initautodock.h input_state.cc input_state.h intnbtable.cc
**  	intnbtable.h investigate.cc investigate.h linpack.cc ls.cc
**  	ls.h main.cc main.h minmeanmax.cc mkNewState.cc mkNewState.h
**  	mkTorTree.cc mkTorTree.h nbe.cc nbe.h nonbonds.cc nonbonds.h
**  	openfile.cc openfile.h output_state.cc output_state.h
**  	parse_PDBQT_line.cc parse_PDBQT_line.h parse_dpf_line.cc
**  	parse_dpf_line.h parse_param_line.cc parse_param_line.h
**  	parse_trj_line.cc parse_trj_line.h parsetypes.cc parsetypes.h
**  	prClusterHist.cc prClusterHist.h prInitialState.cc
**  	prInitialState.h prTorConList.cc prTorConList.h
**  	printEnergies.cc printEnergies.h print_2x.cc print_2x.h
**  	print_atomic_energies.cc print_atomic_energies.h
**  	print_avsfld.cc print_avsfld.h print_rem.cc print_rem.h
**  	printdate.cc printdate.h printhms.cc printhms.h qmultiply.cc
**  	qmultiply.h qtransform.cc qtransform.h quicksort.cc
**  	quicksort.h ranlib.cc ranlib.h readGridMap.cc readPDBQT.cc
**  	readPDBQT.h read_parameter_library.cc read_parameter_library.h
**  	readfield.cc readfield.h readmap.cc readmap.h rep.cc rep.h
**  	setflags.cc setflags.h simanneal.cc simanneal.h sort_enrg.cc
**  	sort_enrg.h stack.cc stack.h stateLibrary.cc stateLibrary.h
**  	stop.cc stop.h success.cc success.h support.cc support.h
**  	swap.cc swap.h timesys.cc timesys.h torNorVec.cc torNorVec.h
**  	torsion.cc torsion.h trilinterp.cc trilinterp.h usage.cc
**  	usage.h warn_bad_file.cc warn_bad_file.h weedbonds.cc
**  	weedbonds.h writePDBQT.cc writePDBQT.h
**  ----------------------------------------------------------------------
**
** Revision 1.19  2010/05/19 19:47:13  mp
** Implemented number-of-evaluations trigger, set by "output_population_statistics"
** DPF keyword, most of the logic is in call_glss.cc
**  Modified Files:	call_glss.cc configure.ac ls.cc ls.h main.cc support.h
**
** Revision 1.18  2010/05/14 21:25:51  mp
** Added printing of "Detailed state:" and "QState:".
** Added DPF keyword "output_population_statistics" to control the detailed
**   "Population at Generation:" logging, with partial implementation,
**   triggered by number of generations only, not number of evaluations (yet).
**  Modified Files:
**  	call_glss.cc call_gs.cc configure.ac dpftoken.h gs.cc gs.h
**  	hybrids.h main.cc parse_dpf_line.cc stateLibrary.cc structs.h
**  	support.cc support.h writePDBQT.cc
**
** Revision 1.17  2010/03/22 20:40:56  mp
** Added reporting state vector for best individual to "Population at Generation:"
** line as underscore-separated string
**  Modified Files: call_glss.cc configure.ac stateLibrary.cc support.cc support.h
**
** Revision 1.16  2010/01/08 20:13:47  mp
** Extended population generations statistics, turned on if outlev>1.
**  Modified Files: RELEASENOTES call_glss.cc support.cc support.h
**
** Revision 1.15  2009/12/12 18:44:20  mp
** Added "Population at Generation:" statistics.
**  Modified Files: call_glss.cc configure.ac gs.cc support.cc support.h
**
** Revision 1.14  2009/05/08 23:02:18  rhuey
** Updated copyright notice in 188 source files
**
** Revision 1.13  2009/05/08 21:46:11  rhuey
** removed debugging comments and print-out
**
** Revision 1.12  2009/04/28 21:12:19  rhuey
** Changed so now Individual does mapping of its genotype into its phenotype and inverse_mapping of its phenotype into its genotype; in both cases returns a reference to itself; added a check for self-assignment
**
** Revision 1.11  2008/06/09 22:27:51  garrett
** Added "end_of_branch[MAX_TORS]" to the Population class, the logic being that every Individual in the Population is the same, so rather than add the overhead to all the Individuals, we added it to the Population.  Also introducted two new methods, set_eob() and get_eob(), to set the end_of_branch[] array, and get values given a key torsion number.  These changes are to support the new "Branch Crossover mode".
**
** Revision 1.10  2007/04/27 06:01:51  garrett
** Added the files necessary for GNU Autotools and the "dot-slash-configure dance"...
**
** Revision 1.9  2007/03/21 06:30:56  garrett
** Created a branch of AutoDock 4 with internal representation of orientations changed from axis-angle nx,ny,nz,ang to quaternion-components qx,qy,qz,qw.  This is intended to avoid rotation singularities of axis-angles near ((1,0,0),0 radians), and to avoid orientational bias in dockings.
**
** Revision 1.8  2006/11/03 02:10:48  garrett
** Significant change.  The initial population is now generated in a different way; previously, the axis that defined the rotation was created by generating uniformly-distributed random numbers in the range REALV_LOW to REALV_HIGH.  The same for the rotation angle (and torsion angles).  Now, we use a range of +/- 1 for the initial unit vector, and +/- PI for the rotation angle (and torsion angles).\
**
** Revision 1.7  2005/10/14 03:10:01  garrett
** Completed the "printPopulationAsCoordsEnergies" member function of the "Population" class, so that it now prints the nonbonded energy and the electrostatics energy, in addition to the translation and total energy for each member of the population.  These numbers are written to the population file at the end of each generation.  The DPF keyword "output_pop_file" expects the name of this population file; if this keyword is not given before a given "ga_run" command, then no population file will be written.
**
** Revision 1.6  2005/09/29 03:34:42  garrett
** Added a new method to the Population class, called "printPopulationAsCoordsEnergies", which is used to print out the translation of the centre of each individual and its total interaction energy.
**
** Revision 1.5  2004/12/07 02:07:53  gillet
** -- fix problem of compilation with g++ 3.3.2 on Linux:
** 	added Genotype(Genotype const &); in support.h
** 	it s definition in support.cc
**
** 	added Individual(Individual const &)
**
** Use the following message to resolve problem:
** http://gcc.gnu.org/ml/gcc-help/2003-10/msg00121.html
**
** You should put a copy constructor in your Anton class.
**
** "Anton(Anton& a)" isn't a copy constructor.  You need an "Anton(Anton const& a)".
**
** If Anton objects cannot be used in a copy constructor, then there are certain operations which Anton objects cannot perform.
**
** It appears that you hit upon one of them.
**
** If the "some code which _needs_ to modify a" is doing so in such a way that the LOGICAL state of the object is not affected, then those data members which are being modified should be marked as "mutable" in the class itself.  For example, certain reference counting schemes.  Another example is the std::string's c_str() method, which is a const method.
**
** If the "some code which _needs_ to modify a" does modify the LOGICAL state of the Anton object being copied from, then that's not kosher for use in a copy constructor.  C'est la vie.
**
** The Standard C++ Library auto_ptr<> is an example of a template class which modifies the state of the copied-from object.  That's one of the reasons that auto_ptr<>'s and STL don't mix (by-and-large).
**
** Or to say it another way, auto_ptr<> doesn't satisfy the contract requirements of STL containers.  Generally speaking.  If someone is REALLY careful, they may be able to use auto_ptr<> in a STL container... but I tend to recommend against it.   The BOOST <www.boost.org> folks have some smart pointer classes that are STL friendly.
**
** Revision 1.4  2004/11/16 23:42:56  garrett
** This is the result of merging the existing CVS respository with the AutoDock 4.0 code.  We have tested the code with a variety of problems: rigid ligand-rigid protein, rigid ligand-flexible protein, flexible ligand-rigid protein and flexible ligand-flexible protein: all four tests passed.  There was a bug fix regarding the flexible ligand-rigid protein case, to do with the absence of a BEGIN_RES record in the PDBQ file. -- GMM & RH
**
** Revision 1.3  2004/02/12 04:32:16  garrett
**
** After a first round of compilation using Apple's IDE, Xcode, various
** warnings have been eliminated (mainly unsigned ints and ints being
** interchanged); and
**
** After using Apple's Shark tool for profiling source code, the
** internal energy calculation has been optimized.
**
** The non-bonded cutoff used in the internal energy calculation has been
** reduced from 64 Angstroms to 8 Angstroms.  Most contributions beyond
** 8 Angstroms are very small, less than -0.001 kcal/mol, even for the
** largest atoms. Also, the conversion from double to int used to
** be done before the if to decide if we were within the cutoff; now
** the square of the distance is used in the comparison, and only if
** we are within the cutoff, do we convert from the double to int.
**
** The version checked in here still uses the type array to lookup
** the energy of interaction for a nonbond; this level of indirection
** can be pre-computed, and this should appear in my next round of checkins
**
** -- Garrett
**
** Revision 1.2  2002/10/30 01:49:15  garrett
** Commented out the #include <iostream.h> lines, since these appeared
** to conflict with <stdio.h>.  Also, added -lsupc++ to the linker
** options for Mac OS X 10.2, which now uses GCC 3.1; this may be
** necessary on GNU/Linux systems that use GCC 3.1.
**
** -- Lindy and Garrett
**
** Revision 1.1.1.1  2001/08/13 22:05:53  gillet
**  import initial of autodock sources
**
*/

#ifndef _SUPPORT_H
#define _SUPPORT_H

#include <stdio.h>
#include "rep.h"
#include "eval.h"
#include "structs.h"

extern FILE *logFile; // for DEBUG only

enum EvalMode { Reset, Always_Eval, Normal_Eval, Always_Eval_Nonbond, Always_Eval_Elec };

typedef struct
{
   unsigned int vector;
   unsigned int index;
} Lookup;

//  For class Genotype, right now assume the user implements the
//  default constructor.
class Genotype
{
   //friend void debug(Genotype &);
   protected:
      //  Could some of these be made static?
      unsigned int number_of_genes;
      unsigned int number_of_vectors; // #vectors in rep_vector
      Lookup *lookup;		      // a table that helps in looking up a gene
      Representation **rep_vector; /* the actual representation of the genotype
				      like arrays of reals, bits, ints */
      unsigned modified : 1; /* used in caching for genotype operators, 
				e.g. crossover */

   public:
      Genotype(void);
      Genotype(const Genotype &);
      Genotype(unsigned int, Representation **const); /* creates a genotype from the
					     representation & total # vectors */
					/* Steffen's comment - representation is apparently unused */
      ~Genotype(void); /* destructor */
      Genotype &operator=(const Genotype &);
      unsigned int num_vectors(void) const; /* e.g. "real,bit,bit,int" would = 4 */
      unsigned int num_genes(void) const; /* returns number_of_genes (see above) */
      RepType gtype(int) const; /* returns the type (real,bit,int) for 
							    a particular gene */
      const Element gread(const int) const;
      const Representation *vread(int) const;
      void write(const Element&, const int); /* not const */
      void write(const unsigned char&, const int); /* not const */
      void write(const int&, const int); /* not const */
      void write(ConstDouble, const int); /* not const */
      void write(const Representation &, const int); /* not const */
      Quat readQuat() const;
      void writeQuat( const Quat& q ); /* not const */
};

//  Should Phenotype automatically evaluate itself upon construction?
class Phenotype
{
   //friend void debug(Phenotype &);
   protected:
      unsigned int number_of_dimensions, number_of_points;
      Lookup *lookup;
      Representation **value_vector;
      double value;
      unsigned evalflag : 1;  //  =1 means that this is the current evaluation
      unsigned modified : 1;  //  =1 means that this has been modified

   public:
      Phenotype(void);
      Phenotype(Eval *);
      Phenotype(const Phenotype &);
      //Phenotype(const Genotype &);//to do
      Phenotype(const unsigned int, Representation **const);
      Phenotype(const unsigned int, Representation **const, Eval *);
      ~Phenotype(void);
      Phenotype &operator=(const Phenotype &);
      Eval *pevaluate; // MPique 2014  could be friend of Individual I think
      RepType gtype(int) const;
      const Element gread(const int) const;
      const Representation *vread(int) const;
      void write(const Element&, const int); /* not const */
      void write(const unsigned char&, const int); /* not const */
      void write(const int& value, const int gene_number); /* not const */
      void write(ConstDouble value, const int gene_number); /* not const */
      void write(const Representation &, const int);
      double evaluate(const EvalMode&) /* not const */ ;  //  This should return evaluation if that's the right answer, and it should evaluate otherwise.
      State make_state(int) const;
      unsigned int num_dimensions(void) const; // Steffen : implementation not found
      unsigned int num_pts(void) const;
      Quat readQuat() const;
      void writeQuat( const Quat& q ); /* not  const */
};

//  This should be an encapsulated class within Population
class Individual
{
   //friend void debug(Individual &);
   public:
      Genotype genotyp;   /* Genotype  is operated upon by *global search* operators */
      Phenotype phenotyp; /* Phenotype  "     "      "   " *local search*  operators, eg SW */
      Molecule *mol;		/* molecule */
      Eval *ievaluate; /* evaluation object from Phenotype  MPique 2014 */
      unsigned long age;	/* age of this individual; gmm, 1998-07-10 */

      Individual(void);
      Individual(Eval *);
      Individual(Individual &); /* copy constructor */
      Individual(Individual const &);
      Individual(Genotype &, Phenotype &);
      ~Individual(void); /* destructor */
      Individual &operator=(const Individual &); /* assignment function for
						    individuals */
      Individual &mapping(void);         //updates phenotype from current genotype values 
      Individual &inverse_mapping(void); //updates genotype from current phenotype values 
      //Phenotype mapping(void); /* takes the genotype and converts it into a phenotype.  */
      //Genotype inverse_mapping(void);  // Scott should do: Also copy Phenotype's value
      double value(EvalMode); /* not const */ /* evaluation of the individual gives its value */ /* not const */
      State state(const int) const; /* state variables in AutoDock */
      void  getMol(Molecule * /* not const */) const; /* converts phenotype to mol's state and returns this individual's mol data */
      void printIndividualsState(FILE *const, const int, const int) const; /* print out the state of this individual */
      void incrementAge(); /* not const */ /* make individual grow 1 generation older */
      int serial; // serial number of this individual
};

class Population
{
   //friend void debug(Population &);
   protected:
      int lhb;  //  These keep track of the lower & upper heap bounds
      int size; /* the number of individuals in the population */
      Individual *heap; /* a heap of individuals -- special binary tree */
      void swap(Individual &, Individual &) const; /* for maintaining the heap order*/
      void SiftUp(void); /* not const */ /* for maintaining the heap order*/
      void SiftDown(void); /* not const */ /* for maintaining the heap order*/
      int end_of_branch[MAX_TORS]; // For Branch Crossover Mode

   public:
      Eval *evaluate; /* evaluation object for individuals within population */
      Population(void);
      Population(int); /* create a pop. with this many individuals */
      Population(int, Eval *); /* create a pop. with this many individuals and eval fcn */
      Population(int, Eval *, Individual *); /* takes an array of ind's and turns into pop. */
      Population(int, Eval *, Individual &); /* takes a prototype ind and turns into pop. */
      Population(const Population &); /* copy constructor */
      ~Population(void); /* destructor */
      Individual &operator[](const int) const;  /* for accessing a particular indiv.in pop*/
      Population &operator=(const Population &);
      unsigned int num_individuals(void) const; /* returns the size of the pop. */
      void msort(const int); /* not const */ /* sorts the first m individuals using heap properties */
      // void print(ostream &, int); /* prints top int energies */
      void print(FILE * const,  const int) const; /* like above */
      // best_e and best_i added M Pique 2010-03 strictly for statistics in log
      //  see printPopulationStatistics (TODO - put in better place)
      double best_e; // best energy
      int best_i; // index in heap[] of indiv with best energy
      int printPopulationStatistics(FILE * const, const int, const char [])  /* not const changed in support.cc hack TODO */; /* prints best, worse, mean, etc energies */
      int printPopulationStatisticsVerbose(FILE *const, const unsigned int, const long int, const int, const char [])  /* not const changed in support.cc hack TODO */; /* print with generations & #evals */
      unsigned long nevals_last_pop_stats; // when pop stats were last printed, see call_glss.cc
      void printPopulationAsStates(FILE *, int, int) const; /*prints energies,states of top energies */
      void printPopulationAsCoordsEnergies(FILE *const, const int, const int) const; /*prints energies,states of top energies */
      void set_eob(int init_end_of_branch[MAX_TORS]); // For Branch Crossover Mode
      int get_eob(int init_tor) const; // For Branch Crossover Mode
};

/**************************************************************************
      Inline Functions
**************************************************************************/

//  The following should be the user's default constructor.  For now,
//  we'll deal with just RealVectors
inline Genotype::Genotype(void)
{
   number_of_genes = number_of_vectors = 0;
   modified = 0;
   rep_vector = (Representation **)NULL;
   lookup = (Lookup *)NULL;
}

inline unsigned int Genotype::num_genes(void) const
{
   return(number_of_genes);
}

inline unsigned int Genotype::num_vectors(void) const
{
   return(number_of_vectors);
}

inline RepType Genotype::gtype(int gene_number) const
{
   return(rep_vector[lookup[gene_number].vector]->type());
}

inline const Element Genotype::gread(const int gene_number) const
{
   return(rep_vector[lookup[gene_number].vector]->gene(lookup[gene_number].index));
}

inline const Representation *Genotype::vread(int vector_number) const
{
   return(rep_vector[vector_number]);
}

//  More user definable stuff
inline Phenotype::Phenotype(void)
{
   value_vector = (Representation **)NULL;
   lookup = (Lookup *)NULL;
   pevaluate = (Eval *)NULL;
   number_of_dimensions = 0;
   number_of_points = 0;
   value = 0;
   evalflag = 0;
}
inline Phenotype::Phenotype(Eval *init_pevaluate)
: pevaluate(init_pevaluate)
{
   value_vector = (Representation **)NULL;
   lookup = (Lookup *)NULL;
   pevaluate = init_pevaluate;
   number_of_dimensions = 0;
   number_of_points = 0;
   value = 0;
   evalflag = 0;
}
inline RepType Phenotype::gtype(int gene_number) const
{
   return(value_vector[lookup[gene_number].vector]->type());
}

inline const Element Phenotype::gread(const int gene_number) const
{
   return(value_vector[lookup[gene_number].vector]->gene(lookup[gene_number].index));
}

inline const Representation *Phenotype::vread(int vector_number) const
{
   return(value_vector[vector_number]);
}

inline unsigned int Phenotype::num_pts(void) const
{
   return(number_of_points);
}

//  Constructs an Individual using the default constructors
inline Individual::Individual(void)
{
   ievaluate = (Eval *) NULL;
   age = 0;
}
inline Individual::Individual(Eval *init_evaluate)
{
   age = 0L;
   ievaluate = init_evaluate;
   phenotyp.pevaluate = init_evaluate; // copy Eval reference into phenotype
}


inline Individual::Individual(Individual &original) // 'copy' constructor Individual i = oldi;
: genotyp(original.genotyp), phenotyp(original.phenotyp), ievaluate(original.ievaluate)
{
}

// caution, does not do mapping
inline Individual::Individual(Genotype &init_genotyp, Phenotype &init_phenotyp)
: genotyp(init_genotyp), phenotyp(init_phenotyp), ievaluate(init_phenotyp.pevaluate)
{
}

inline Individual::~Individual(void)
{
}

// MPique TODO how is this different from copy constructor above??
inline Individual &Individual::operator=(const Individual &original)
{
   if (this == &original) {//Prevent self assignment
      return *this ;
   }
   genotyp = original.genotyp;
   phenotyp = original.phenotyp;
   ievaluate = original.ievaluate;
   mol = original.mol;
   age = original.age;
   return(*this);
}

inline double Individual::value(EvalMode mode) /* not const */
{ // TODO: check if mapping from genotyp to phenotyp is up-to-date
  // note that phenotyp.evaluate only does evaluation if evalflag is false
   return(phenotyp.evaluate(mode));
}

inline Population::Population(void)
:lhb(-1), size(0), heap((Individual *)NULL)
{
#ifdef DEBUGPOPCONSTRUCTORS
fprintf(logFile, "support.h %d Population::Population(void)\n",  __LINE__); 
#endif
    evaluate = (Eval *) NULL;
    for (int i=0; i<MAX_TORS; i++) {
        end_of_branch[i] = -1;
    }
}

inline Population::Population(const int num_inds, Eval *init_evaluate)
: lhb(num_inds-1), size(num_inds), evaluate(init_evaluate)
{
    heap = new Individual[num_inds];
    for (int i=0; i<num_inds; i++) heap[i] = Individual(evaluate);
    for (int i=0; i<MAX_TORS; i++) {
        end_of_branch[i] = -1;
    }
}

inline Population::Population(const int num_inds, Eval *init_evaluate, 
 Individual & prototype)
: size(num_inds), evaluate(init_evaluate)
{
#ifdef DEBUGPOPCONSTRUCTORS
fprintf(logFile, "support.h %d Population::Population(%d/%d,Eval*,Indiv&)\n",
 __LINE__, size, num_inds); // DEBUG MPique
#endif
    heap = new Individual[num_inds];
    for (int i=0; i<num_inds; i++) {
	heap[i] = Individual(prototype);
	heap[i].ievaluate = evaluate;
	heap[i].phenotyp.pevaluate = evaluate;
    }
 //  Do initialization stuff
    for (int i=0; i<MAX_TORS; i++) {
        end_of_branch[i] = -1;
    }
}

inline Population::Population(const int newpopsize, Eval *init_evaluate, 
 Individual *const newpop)
: size(newpopsize), heap(newpop), evaluate(init_evaluate)
{
#ifdef DEBUGPOPCONSTRUCTORS
fprintf(logFile, "support.h %d Population::Population(%d/%d,Eval*,Indiv*)\n",__LINE__,size,newpopsize); // DEBUG MPique
#endif
   //  Do initialization stuff
    for (int i=0; i<MAX_TORS; i++) {
        end_of_branch[i] = -1;
    }
}

inline Population::~Population(void)
{
   if(heap != (Individual *)NULL)
   {
      delete [] heap;
   }
}

inline unsigned int Population::num_individuals(void) const
{
   return(size);
}

#endif
