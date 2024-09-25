/*

 $Id: gs.h,v 1.26 2014/06/12 01:44:07 mp Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

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

//  These are the classes in the global search hierarchy.  Notice that
//  Global_Search is an abstract base class and as such it can never
//  be instantiated.  For now, the only global search operator is the
//  Genetic Algorithm.  rsh 07/08/95
//  ParticleSwarmGS added as another global search spring-summer 2011 rsh and mp

#ifndef _GLOBAL_SEARCH_H
#define _GLOBAL_SEARCH_H

#include "support.h"
#include "ls.h"

enum M_mode { ERR = -1, BitFlip, CauchyDev, IUniformSub };
enum Selection_Mode { Proportional=0, LinearRanking=1, Tournament=2, Boltzmann=3 };
enum Xover_Mode { TwoPt=0, OnePt=1, Uniform=2, Arithmetic=3, Branch=4 };
enum Worst_Mode { AverageOfN, OfN, Ever };

class Global_Search
{
   public:
      Global_Search(unsigned int init_max_evals, unsigned int init_max_generations);
      virtual ~Global_Search(void);
      virtual int search(Population &, Eval *, int, FILE *) = 0;
      virtual int localsearch(Population &, Local_Search *, Eval *, int, FILE *) = 0;
      virtual int terminate(void) = 0;
      virtual void reset(void) = 0;
      virtual void reset(const Output_pop_stats&) = 0;
      virtual char * shortname(void) = 0;
      virtual char * longname(void) = 0;
      unsigned int max_evals ;
      unsigned int max_generations ;
      // the next four are only applicable to Genetic_Algorithm
      unsigned int cg_count; // statistics - crossover gene-by-gene count
      unsigned int ci_count; // statistics - crossover indiv-by-indiv count
      unsigned int mg_count; // statistics - mutation gene-by-gene count
      unsigned int mi_count; // statistics - mutation indiv-by-indiv count
};
inline Global_Search::Global_Search(
        const unsigned int init_max_evals, 
        const unsigned int init_max_generations
    ) : max_evals(init_max_evals), max_generations(init_max_generations),
    cg_count(0), ci_count(0), mg_count(0), mi_count(0)
{
}


// The class Genetic_Algorithm is a Global_Search method,
// 
class Genetic_Algorithm : public Global_Search
{
//   friend void debug(Genetic_Algorithm &, Population &);
   private:
      EvalMode e_mode;
      Selection_Mode s_mode;
      Xover_Mode c_mode;
      Worst_Mode w_mode;
      unsigned int elitism;
	  Real c_rate;
	  Real m_rate;
      Real localsearch_freq;
      unsigned int window_size;
      Real alpha;
	  Real beta;
      Real tranStep, quatStep, torsStep;
      int low, high; // should these be int or Real?
      unsigned int generations; 
      Output_pop_stats output_pop_stats;// gmm 2000.11.1,2003.08.18, MPique 2010.05 
      unsigned int converged; // gmm 7-jan-98
 	  Real *alloc;
      Real *mutation_table;
      unsigned int *ordering;	  
      unsigned int m_table_size;
      double worst, avg;
      double *worst_window;
      Real linear_ranking_selection_probability_ratio;

      double worst_this_generation(const Population &, int outlev, FILE *logFile);
      void set_worst(const Population &, int outlev, FILE *logFile);
      void make_table(int, ConstReal, int outlev, FILE *logFile);
      int check_table(ConstReal, int outlev, FILE *logFile);
      M_mode m_type(const RepType) const;
      void mutate(Genotype &, const int, int outlev, FILE *logfile);
      void mutation(Population &, int outlev, FILE *logfile);
      void crossover(Population &, int outlev, FILE *logfile);
      void crossover_2pt(Genotype &, Genotype &, const unsigned int, const unsigned int, int outlev, FILE *logFile);
      void crossover_uniform(Genotype &, Genotype &, const unsigned int, int outlev, FILE *logFile);
      void crossover_arithmetic(Genotype &, Genotype &, ConstReal, int outlev, FILE *logFile );
      void selection_proportional(Population &, Individual* const, int outlev, FILE *logFile);
      void selection_linear_ranking(/* sorted */ Population &, /* not const */ Individual *const, int outlev, FILE *logFile);
      void selection_tournament(Population &, Individual* const, int outlev, FILE *logFile);
      Individual *selection(Population &, int outlev, FILE *logfile);

   public:
      Genetic_Algorithm(void);
      // Genetic_Algorithm(EvalMode, Selection_Mode, Xover_Mode, Worst_Mode, int, Real, Real, int, unsigned int); // before 2000.11.1
      //Genetic_Algorithm(EvalMode, Selection_Mode, Xover_Mode, Worst_Mode, int, Real, Real, int, unsigned int, unsigned int); // after 2000.11.1
      Genetic_Algorithm(const EvalMode init_e_mode, const Selection_Mode init_s_mode,
			const Xover_Mode init_c_mode, const Worst_Mode init_w_mode, const int init_elitism,
                        ConstReal  init_c_rate, ConstReal  init_m_rate,
			ConstReal init_localsearch_freq,
                        const int init_window_size, 
			const unsigned int init_max_evals,
			const unsigned int init_max_generations,
                        const Output_pop_stats&); // after 2010.05
      ~Genetic_Algorithm(void);
      void initialize(unsigned int, unsigned int, int, FILE *);
      void mutation_values(int, int, ConstReal, ConstReal, ConstReal, ConstReal, ConstReal);
      unsigned int num_generations(void) const;
      void reset(void);
      void reset(const Output_pop_stats&);
      char * shortname(void);
      char * longname(void);
      int terminate(void);
      int search(Population &, Eval *, int, FILE *);
      int localsearch(Population &, Local_Search *, Eval *, int, FILE *);
      int set_linear_ranking_selection_probability_ratio(ConstReal );
      
};

//  Inline Functions

inline Global_Search::~Global_Search(void)
{
}

#ifdef VOIDCONSTRUCTOR
// Default values set in this constructor.
inline Genetic_Algorithm::Genetic_Algorithm(void)
: alloc(NULL), mutation_table(NULL), ordering(NULL), m_table_size(0), worst_window(NULL)
{
   generations = 0;
   elitism = window_size = low = high = 0;
   m_rate = 0.02;
   c_rate = 0.80;
   localsearch_freq = 0.06;
   alpha = beta = 0.0;
   tranStep = 2.0;
   quatStep = torsStep = DegreesToRadians( 30.0 );
   worst = avg = 0.0L;
   converged = 0; // gmm 7-jan-98
   output_pop_stats.level = 0;
   output_pop_stats.everyNgens = OUTLEV1_GENS; // gmm 2000-nov-1
   output_pop_stats.everyNevals = 0;
   linear_ranking_selection_probability_ratio = 2.0; //mp+rh 10/2009
}
#endif

inline Genetic_Algorithm::~Genetic_Algorithm(void)
{
   if (worst_window!=NULL) {
      delete [] worst_window;
   }

   if (alloc!=NULL) {
      delete [] alloc;
   }

   if (ordering!=NULL) {
      delete [] ordering;
   }

   if (mutation_table!=NULL) {
      delete [] mutation_table;
   }
}

inline void Genetic_Algorithm::mutation_values(const int init_low, const int init_high, 
        ConstReal init_alpha, ConstReal  init_beta, 
        ConstReal init_tranStep, ConstReal init_quatStep, ConstReal init_torStep )
{
   low = init_low;
   high = init_high;
   alpha = init_alpha;
   beta = init_beta;
   tranStep = init_tranStep;
   quatStep = init_quatStep;
   torsStep = init_torStep;
}

inline char * Genetic_Algorithm::shortname(void)
{
        return "GA";
}

inline char * Genetic_Algorithm::longname(void)
{
        return "GENETIC ALGORITHM";
}
inline unsigned int Genetic_Algorithm::num_generations(void) const
{
   return(generations);
}

inline int Genetic_Algorithm::terminate(void)
{
   if (max_generations>0) {
      // before 7-jan-98, was: return(generations>=max_generations);
      return((generations>=max_generations)||(converged==1)); // gmm 7-jan-98
   } else {
      return(0);  //  Don't terminate
   }
}

inline void Genetic_Algorithm::reset(void)
{
   generations = 0;
   converged = 0; // gmm 7-jan-98
   cg_count = ci_count = mg_count = mi_count = 0; // restart statistics
}

inline void Genetic_Algorithm::reset(const Output_pop_stats& extOutput_pop_stats)
{
   output_pop_stats = extOutput_pop_stats; // gmm 2000.11.1   MPique 2010.05
   generations = 0;
   converged = 0; // gmm 7-jan-98
   cg_count = ci_count = mg_count = mi_count = 0; // restart statistics
}

#endif
