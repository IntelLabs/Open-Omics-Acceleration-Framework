/*

 $Id: gs.cc,v 1.55 2014/06/12 01:44:07 mp Exp $

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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/********************************************************************
     These are the methods for Global_Search and its derivations.

                                rsh 9/95
*********************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include "gs.h"
#include "ranlib.h"
#include "eval.h"
#include "rep.h"
#include "rep_constants.h"
#include "assert.h"

#ifdef sgi
    #include <ieeefp.h>
#endif
#ifdef sun
    #include <ieeefp.h>
#endif

#include "constants.h"
#include "autocomm.h"
#include "writePDBQT.h"
#include "qmultiply.h"


extern int sel_prop_count;//debug
extern int debug;//debug
//#define DEBUG
//#define DEBUG2
//#define DEBUG3
//#define DEBUG_MUTATION


static double
worst_in_window(const double *const window, const int size, int outlev, FILE *logFile) 
{
   register int i;

   double worst = window[0];

#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/double worst_in_window(double *window, int size)_________________________\n");//debug
#endif

   for (i=1; i<size; i++) {

#ifdef DEBUG2
      (void)fprintf(logFile, "gs.cc/window[%d]= %.3f\tworst= %.3f\n", i, window[i], worst);//debug
#endif

      if (window[i]>worst) {
         worst = window[i];

#ifdef DEBUG2
         (void)fprintf(logFile, "gs.cc/i= %d\t(window[i]>worst)\tUpdating: worst= %.3f\n", i, worst);//debug
#endif

      }
   }// for i

#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/Returning: worst= %.3f\n\n", worst);//debug
#endif

   return(worst);
}

static double
avg_in_window(const double *const window, const int size, FILE *logFile) 
{
   register int i;
   double mysum = 0.0;

#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/avg_in_window(double *window, int size)_________________________\n");//debug
#endif
   for (i=0; i<size; i++) {
      mysum += window[i];
#ifdef DEBUG2
      (void)fprintf(logFile, "gs.cc/mysum= %.3f\twindow[%d]= %.3f\n",mysum, i, window[i]);//debug
#endif
   }
   const double myavg = mysum / size;
#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/Returning: myavg= %.3f\n\n",myavg);//debug
#endif

   return(myavg);
}

//  Also set avg -- and because of avg this is not a const function
double Genetic_Algorithm::worst_this_generation(const Population &pop, int outlev, FILE *logFile)
{
   register unsigned int i;
   double worstval, avgval;

#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/worst_this_generation(Population &pop, logFile)_________________________\n");
#endif

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/double Genetic_Algorithm::worst_this_generation(Population &pop)\n");
#endif /* DEBUG */

   avgval = worstval = pop[0].value(Normal_Eval);
#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/avgval= %.3f\tworstval= %.3f\n", avgval, worstval);
#endif
   for (i=1; i<pop.num_individuals(); i++) {
      avgval += pop[i].value(Normal_Eval);
#ifdef DEBUG2
      (void)fprintf(logFile, "gs.cc/avgval= %.3f\tpop[%d].value(Normal_Eval)= %.3f\n", avgval, i, pop[i].value(Normal_Eval));
#endif
      if (pop[i].value(Normal_Eval)>worstval) {
         worstval = pop[i].value(Normal_Eval);
#ifdef DEBUG2
         (void)fprintf(logFile, "gs.cc/(pop[i].value(Normal_Eval)>worstval): Updating: worstval= %.3f\n", worstval);
#endif
      }
   }

   avg = avgval/pop.num_individuals();
#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/Returning: avg= %.3f, worstval= %.3f\n\n", avg, worstval);
#endif
   return(worstval);
}

//  This could be made inline

Genetic_Algorithm::Genetic_Algorithm( const EvalMode init_e_mode, 
                                      const Selection_Mode init_s_mode, 
                                      const Xover_Mode init_c_mode,
                                      const Worst_Mode init_w_mode, 
                                      const int init_elitism, 
                                      ConstReal  init_c_rate, 
                                      ConstReal  init_m_rate, 
                                      ConstReal  init_localsearch_freq, 
                                      const int init_window_size, 
                                      const unsigned int init_max_evals,
                                      const unsigned int init_max_generations,
                                      const Output_pop_stats& init_output_pop_stats)

: Global_Search(init_max_evals, init_max_generations),
e_mode(init_e_mode),
s_mode(init_s_mode),
c_mode(init_c_mode),
w_mode(init_w_mode),
elitism(init_elitism),
c_rate(init_c_rate),
m_rate(init_m_rate),
localsearch_freq(init_localsearch_freq),
window_size(init_window_size),
alpha(1.0),
beta(0.0),
tranStep(2.0),                         //  2 Angstroms
quatStep( DegreesToRadians( 30.0 ) ),  // 30 degrees
torsStep( DegreesToRadians( 30.0 ) ),  // 30 degrees
low(-100),
high(100),
generations(0),
output_pop_stats(init_output_pop_stats),
converged(0),
alloc(NULL),
mutation_table(NULL),
ordering(NULL),
m_table_size(0),
worst(0.0L),
avg(0.0L),
linear_ranking_selection_probability_ratio(2.0)

{
#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/Genetic_Algorithm::Genetic_Algorithm(EvalMode init_e_mode,...\n");
#endif /* DEBUG */

   worst_window = new double[window_size];
}

int Genetic_Algorithm::set_linear_ranking_selection_probability_ratio(ConstReal  r)
{
    if (r<0.) return -1;  //ERROR!
    linear_ranking_selection_probability_ratio = r;
    return 1;
}



void Genetic_Algorithm::set_worst(const Population &currentPop, int outlev, FILE *logFile)
{
   double temp = 0.0;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::set_worst(Population &currentPop)\n");
#endif /* DEBUG */

   worst_window[generations%window_size] = worst_this_generation(currentPop, outlev, logFile);
   switch(w_mode)
   {
      //  Assume for this case that there's a window_size of 1
      case Ever:
         if (generations!=0) {
            if (temp>worst)
               worst = worst_window[0];
         } else {
            worst = worst_window[0];
         }
         break;
      case OfN:
         if (generations>=window_size) {
            worst = worst_in_window(worst_window, window_size, outlev, logFile);
         } else {
            worst = worst_in_window(worst_window, generations+1, outlev, logFile);
         }
         break;
      case AverageOfN:
         if (generations>=window_size) {
            worst = avg_in_window(worst_window, window_size, logFile);
         } else {
            worst = avg_in_window(worst_window, generations+1, logFile);
         }
         break;
      default:
         (void)fprintf(logFile,"gs.cc/Unable to set the individual with the worst fitness!\n");
   }
}

M_mode Genetic_Algorithm::m_type(const RepType type) const
{

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/M_mode Genetic_Algorithm::m_type(RepType type)\n");
#endif /* DEBUG */
   switch(type)
   {
      case T_BitV:
         return(BitFlip);
      case T_RealV:
      case T_CRealV:
         return(CauchyDev);
      case T_IntV:
         return(IUniformSub);
      default:
         stop("gs.cc/Unrecognized Type (The allowed types are:  T_BitV, T_RealV, T_CRealV and T_IntV)!\n");
         return(ERR); // NOTREACHED
   }
}

void Genetic_Algorithm::make_table(const int size, ConstReal  prob, int outlev, FILE *logFile)
{
   register int i, j;
   double L = 0.0L;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::make_table(int size=%d, Real prob=%f)\n",size, prob);
#endif /* DEBUG */

#ifdef DEBUG_MUTATION
   (void)fprintf(logFile, "\ngs.cc/void Genetic_Algorithm::make_table(int size=%d, Real prob=%f)\n",size, prob);
#endif /* DEBUG */

   m_table_size = size;
   mutation_table = new Real[size+1];

   mutation_table[0] = pow(1-prob, size);
   mutation_table[size] = 1;

#ifdef DEBUG_MUTATION
   fprintf(logFile,"mutation_table[0] = %.3f\n",  mutation_table[0]);
   fprintf(logFile,"mutation_table[%d] = %.3f\n",  size,mutation_table[size]);
#endif
   i = 1;
   while (i<=(int)size*prob) {
      L = 0.0;
      for (j=1; j<=i; j++) {
         L += log(size+1-j) - log(j);
#ifdef DEBUG_MUTATION
         fprintf(logFile,"j= %d (size+1-j)= %d log(_)= %.4f log(j)=%.4f L= %.4f\n", j, (size+1-j), log(size+1-j), log(j), L);
#endif
      }
      L += i*log(prob) + (size-i)*log(1-prob);
#ifdef DEBUG_MUTATION
      fprintf(logFile,"j= %d (size+1-j)= %d log(_)= %.4f log(j)=%.4f L= %.4f\n", j, (size+1-j), log(size+1-j), log(j), L);
#endif

      assert( i > 0 && i<=size); // M Pique 2009-12 TODO - suspect problem if m_rate>1
      mutation_table[i] = mutation_table[i-1]+exp(L);
#ifdef DEBUG_MUTATION
      fprintf(logFile,"mutation_table[%d] = %.3f\n", i, mutation_table[i]);
#endif
      i++;
   } // end while

#ifdef DEBUG_MUTATION
      fprintf(logFile,"L= %.3f  exp(L)= %.3f\n", L, exp(L) );
#endif
   L = exp(L);
   for (; i<size; i++) {
#ifdef DEBUG_MUTATION
      fprintf(logFile,"i= %d  L= %.3f  prob= %.3f size= %d  (size+1-i)= %d  i*(1-prob)= %.3f\n", i, L, prob, size, (size+1-i), i*(1-prob) );
#endif
      L = (L*prob*(size+1-i))/(i*(1-prob));
#ifdef DEBUG_MUTATION
      fprintf(logFile,"L= %.3f\n", L );
#endif
      mutation_table[i] = mutation_table[i-1]+L;
#ifdef DEBUG_MUTATION
      fprintf(logFile,"mutation_table[%d] = %.3f\n", i, mutation_table[i]);
#endif
   }
#ifdef DEBUG_MUTATION
      fflush(logFile);
#endif
}

int Genetic_Algorithm::check_table(ConstReal  prob, int outlev, FILE *logFile)
{
   int low, high;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/int Genetic_Algorithm::check_table(Real prob=%f)\n",prob);
#endif /* DEBUG */

#ifdef DEBUG_MUTATION
   (void)fprintf(logFile, "\ngs.cc/int Genetic_Algorithm::check_table(Real prob=%f)\n",prob);
#endif /* DEBUG */

   low = 0; high = m_table_size;
#ifdef DEBUG_MUTATION
      fprintf(logFile,"prob= %.3f  low= %d  high= %d\n", prob, low, high );
#endif

   while (high-low>1) {
#ifdef DEBUG_MUTATION
      fprintf(logFile,"high-low= %d\n", high-low );
      fprintf(logFile,"(high+low)/2= %d  mutation_table[(high+low)/2]= %.3f  prob= %.3f\n", (high+low)/2,  mutation_table[(high+low)/2], prob );
#endif
      if (mutation_table[(high+low)/2]<prob) {
         low = (high+low)/2;
#ifdef DEBUG_MUTATION
         fprintf(logFile,"low= %d\n", low );
#endif
      } else if (mutation_table[(high+low)/2]>prob) {
         high = (high+low)/2;
#ifdef DEBUG_MUTATION
         fprintf(logFile,"high= %d\n", high );
#endif
      } else {
         high = low = (high+low)/2;
#ifdef DEBUG_MUTATION
         fprintf(logFile,"low= %d  high= %d\n", low, high );
#endif
      }
   }
   return(low);
}

void Genetic_Algorithm::initialize(const unsigned int pop_size, const unsigned int num_poss_mutations, int outlev, FILE *logFile)
{
   register unsigned int i;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::initialize(unsigned int pop_size=%d, ",pop_size);
   (void)fprintf(logFile, "unsigned int num_poss_mutations=%d)\n",num_poss_mutations);
#endif /* DEBUG */

   if (alloc!=NULL) {
      delete [] alloc;
   }

   if (ordering!=NULL) {
      delete [] ordering;
   }

   if (mutation_table!=NULL) {
      delete [] mutation_table;
   }

   alloc = new Real[pop_size];

   ordering = new unsigned int[pop_size];
   for (i=0; i<pop_size; i++) {
      ordering[i] = i;
      assert(ordering[i] < pop_size);//debug
      alloc[i] = 1.0; // changed by gmm, 12-sep-1997.
   }

   make_table(pop_size*num_poss_mutations, m_rate, outlev, logFile);
}

void Genetic_Algorithm::mutate(Genotype &mutant, const int gene_number, int outlev, FILE *logFile)
{
   Element tempvar;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::mutate(Genotype &mutant, int gene_number=%d)\n",gene_number);
#endif /* DEBUG */

#ifdef DEBUG_MUTATION
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::mutate(Genotype &mutant, int gene_number=%d)\n",gene_number);
#endif /* DEBUG */

   switch(m_type(mutant.gtype(gene_number))) {

      case BitFlip:
#ifdef DEBUG_MUTATION
         fprintf( logFile, "case BitFlip:  gene_number=%d\n", gene_number );
#endif
         //((unsigned char *)gene)[point] = 1 - ((unsigned char *)gene)[point];
         //  Read the bit
         tempvar = mutant.gread(gene_number);
         //  Flip it
         tempvar.bit = 1 - tempvar.bit;
         //  write it
         mutant.write(tempvar, gene_number);
         break;

      case CauchyDev:
#ifdef DEBUG_MUTATION
         fprintf( logFile, "case CauchyDev:  gene_number=%d\n", gene_number );
#endif
         // gene_numbers 3, 4, 5 and 6 correspond to the 
         // four components of the quaternion, (x,y,z,w)
         if ( is_rotation_index( gene_number ) ) {
#ifdef DEBUG_MUTATION
            fprintf( logFile, "is_rotation_index( gene_number=%d )\n", gene_number );
#endif
            // Mutate all four comopnents of the quaternion, (x,y,z,w) simultaneously:
            // Generate a uniformly-distributed quaternion
            Quat q_change;
	    // MP TODO here we could use quat by amount for better mutations
            q_change = randomQuat();
#ifdef DEBUG_MUTATION
            fprintf( logFile, "q_change -- after randomQuat\n" );
            printQuat_q( logFile, q_change );
#endif
            Quat q_current = mutant.readQuat();
#ifdef DEBUG_MUTATION
            fprintf( logFile, "q_current -- after mutant.readQuat()\n" );
            printQuat_q( logFile, q_current );
#endif
            Quat q_new;
#ifdef DEBUG_MUTATION
            fprintf( logFile, "q_current -- about to call qmultiply\n" );
#endif
#ifdef DEBUG_QUAT
#ifdef DEBUG_QUAT_PRINT
            pr( logFile, "DEBUG_QUAT: mutate() -- q_current\n" );
            (void) fflush(logFile);
#endif // DEBUG_QUAT_PRINT
            //  Make sure the quaternion is suitable for 3D rotation
            assertQuatOK( q_current );
#ifdef DEBUG_QUAT_PRINT
            pr( logFile, "DEBUG_QUAT: mutate() -- q_change\n" );
            (void) fflush(logFile);
#endif // DEBUG_QUAT_PRINT
            //  Make sure the quaternion is suitable for 3D rotation
            assertQuatOK( q_change );
#endif // DEBUG_QUAT
            // Multiply the quaternions, applying the rotation to the current orientation
            qmultiply( &q_new, &q_current, &q_change );
#ifdef DEBUG_QUAT
#ifdef DEBUG_QUAT_PRINT
            pr( logFile, "DEBUG_QUAT: mutate() -- q_current\n" );
            (void) fflush(logFile);
#endif
            //  Make sure the quaternion is suitable for 3D rotation
            assertQuatOK( q_new );
#endif // DEBUG_QUAT
#ifdef DEBUG_MUTATION
            fprintf( logFile, "q_new - after qmultiply\n" );
            printQuat_q( logFile, q_new );
#endif
            mutant.writeQuat( q_new );
         } else {
             //  Read the real
             tempvar = mutant.gread(gene_number);
#ifdef DEBUG_MUTATION
             (void)fprintf(logFile, "   ---CauchyDev---\n" );
             (void)fprintf(logFile, "   Before mutating:        tempvar= %.3f\n", tempvar.real );
             (void)fprintf(logFile, "   gene_number= %d\n", gene_number );
             (void)fprintf(logFile, "   tempvar.real += scauchy2()\n" );
#endif
             //  Add deviate
             //  We never vary alpha and beta, so just use the faster "scauchy2()" function:
             tempvar.real += scauchy2();
#ifdef DEBUG_MUTATION
             (void)fprintf(logFile, "   Add Cauchy deviate:     tempvar= %.3f\n", tempvar.real );
             (void)fflush(logFile );
#endif
             //  Write it
             mutant.write(tempvar, gene_number);
         }
         break;

      case IUniformSub:
#ifdef DEBUG_MUTATION
         fprintf( logFile, "case IUniformSub:  gene_number=%d\n", gene_number );
#endif
         //((int *)gene)[point] = ignuin(low, high);
         //  Generate the new integer
         tempvar.integer = ignuin(low, high);
         //  Write it
         //mutant.write((void *)(&tempvar.integer), gene_number);
         mutant.write(tempvar, gene_number);
         break;

      default:
         (void)fprintf(logFile,"gs.cc/Unrecognized mutation Mode!\n");
         break; 
   }
}

void Genetic_Algorithm::mutation(Population &pure, int outlev, FILE *logFile)
{
   int num_mutations, individual, gene_number;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::mutation(Population &pure)\n");
#endif /* DEBUG */

#ifdef DEBUG_MUTATION
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::mutation(Population &pure)\n");
#endif /* DEBUG */

   num_mutations = check_table(ranf(), outlev, logFile);
#ifdef DEBUG_MUTATION
   (void)fprintf(logFile, "num_mutations= %d\n", num_mutations );
#endif /* DEBUG */

   if(num_mutations<=0) return;

   Boole *individual_mutated = new Boole[pure.num_individuals()]; // for statistics only
   for(unsigned int i=0;i<pure.num_individuals();i++) individual_mutated[i]=FALSE;

   //  Note we don't check to see if we mutate the same gene twice.
   //  So, effectively we are lowering the mutation rate, etc...
   //  Might want to check out Bentley's chapter on selection.
   for (; num_mutations>0; num_mutations--) {
      individual = ignlgi()%pure.num_individuals();
      gene_number = ignlgi()%pure[individual].genotyp.num_genes();
#ifdef DEBUG_MUTATION
      (void)fprintf(logFile, "  @__@  mutate(pure[individual=%d].genotyp, gene_number=%d);\n\n", individual, gene_number );
#endif /* DEBUG */
      mutate(pure[individual].genotyp, gene_number, outlev, logFile);
      pure[individual].age = 0L;
      pure[individual].mapping();//map genotype to phenotype and evaluate

      mg_count++; // update statistics: count of mutated genes per run
      individual_mutated[individual] = TRUE; // for statistics only
   }
   // update statistics: count of mutated individuals per run
   for(unsigned int i=0;i<pure.num_individuals();i++) \
     if(individual_mutated[i]) mi_count++; 
   delete [] individual_mutated;
}

void Genetic_Algorithm::crossover(Population &original_population, int outlev, FILE *logFile)
{
   register unsigned int i;
   int first_point, second_point, temp_index, temp_ordering;
   Real alpha = 0.5;
   
#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover(Population &original_population)\n");
#endif /* DEBUG */

   //  Shuffle the population
   //  Permute the ordering of the population, "original_population"
   for (i=0; i<original_population.num_individuals(); i++) {
      temp_ordering = ordering[i];
      // ignlgi is GeNerate LarGe Integer and is in com.cc
      temp_index = ignlgi()%original_population.num_individuals();
      ordering[i] = ordering[temp_index];
      assert(ordering[i] < original_population.num_individuals());//debug
      assert(ordering[i] >= 0);//debug
      ordering[temp_index] = temp_ordering;
      assert(ordering[temp_index] < original_population.num_individuals());//debug
      assert(ordering[temp_index] >= 0);//debug
   }
                                                                    
   //  How does Goldberg implement crossover?

   // Loop over individuals in population
   for (i=0; i<original_population.num_individuals()-1; i=i+2) {
      // The two individuals undergoing crossover are original_population[ordering[i]] and original_population[ordering[i+1]]

      if (ranf() < c_rate) {
         // Perform crossover with a probability of c_rate
         int fi = ordering[i]; //index of father
         int mi = ordering[i+1]; //index of mother 
         // Assert the quaternions of the mother and father are okay:
         Genotype father = original_population[fi].genotyp;
         Genotype mother = original_population[mi].genotyp;
         Quat q_father, q_mother;
#ifdef DEBUG_QUAT_PRINT
         pr( logFile, "DEBUG_QUAT: crossover()  q_father (individual=%d)\n", fi );
#endif // endif DEBUG_QUAT_PRINT
         //  Make sure the quaternion is suitable for 3D rotation
         q_father = father.readQuat();
#ifdef DEBUG_QUAT_PRINT
         printQuat( logFile, q_father );
         (void) fflush(logFile);
#endif // endif DEBUG_QUAT_PRINT
         assertQuatOK( q_father );
#ifdef DEBUG_QUAT_PRINT
         pr( logFile, "DEBUG_QUAT: crossover()  q_mother (individual=%d)\n", mi );
         (void) fflush(logFile);
#endif // endif DEBUG_QUAT_PRINT
         //  Make sure the quaternion is suitable for 3D rotation
         q_mother = mother.readQuat();
#ifdef DEBUG_QUAT_PRINT
         printQuat( logFile, q_mother );
         (void) fflush(logFile);
#endif // endif DEBUG_QUAT_PRINT
         assertQuatOK( q_mother );

                // Pos.  Orient. Conf.
                // 0 1 2 3 4 5 6 7 8 9
                // X Y Z x y z w t t t
                //                     num_genes() = 10
         switch(c_mode) {
            case OnePt:
                // To guarantee 1-point always creates 2 non-empty partitions,
                // the crossover point must lie in range from 1 to num_genes-1.
                // Choose one random integer in this range
                // Make sure we do not crossover inside a quaternion...
                do {
                    // first_point = ignlgi()%original_population[i].genotyp.num_genes();
                    first_point = ignuin(1, original_population[i].genotyp.num_genes()-1);
                } while ( is_within_rotation_index( first_point ) ) ;
                //  We can accomplish one point crossover by using the 2pt crossover operator
                crossover_2pt( original_population[fi].genotyp, 
                               original_population[mi].genotyp,
                               first_point, 
                               original_population[fi].genotyp.num_genes(), outlev, logFile);//either one works
                break;
            case TwoPt:
                // To guarantee 2-point always creates 3 non-empty partitions,
                // each crossover point must lie in range from 1 to num_genes-1.
                // Choose two different random integers in this range
                // Make sure we do not crossover inside a quaternion...
                do {
                    first_point = ignuin(1, original_population[i].genotyp.num_genes()-1);
                } while ( is_within_rotation_index( first_point ) ) ;
                do {
                    second_point = ignuin(1, original_population[i].genotyp.num_genes()-1);
                } while ( is_within_rotation_index( second_point ) || second_point == first_point );
                // Do two-point crossover, with the crossed-over offspring replacing the parents in situ:
                crossover_2pt( original_population[fi].genotyp, 
                               original_population[mi].genotyp, 
                               min( first_point, second_point ),
                               max( first_point, second_point), outlev, logFile );
                break;
            case Branch:
                // New crossover mode, designed to exchange just one corresponding sub-trees (or "branch")
                // between two individuals.
                // If there are only position and orientation genes, there will be only
                // 7 genes; this mode would not change anything in such rigid-body dockings.
                if (original_population[i].genotyp.num_genes() <= 7) {
                    // Rigid body docking, so no torsion genes to crossover.
                    continue; //TODO raise a fatal error if branch crossover with no torsions
                 } else {
                    // Pick a random torsion gene
                    first_point = ignuin(7, original_population[i].genotyp.num_genes()-1);
                    second_point = original_population.get_eob( first_point - 7 );
                    // Do two-point crossover, with the crossed-over offspring replacing the parents in situ:
                    crossover_2pt( original_population[fi].genotyp, 
                                   original_population[mi].genotyp, 
                                   min( first_point, second_point ),
                                   max( first_point, second_point ), outlev, logFile );
                    break;
                 }
            case Uniform:
                crossover_uniform( original_population[fi].genotyp, 
                                   original_population[mi].genotyp,
                                   original_population[mi].genotyp.num_genes(), outlev, logFile);
                break;
            case Arithmetic:
               // select the parents A and B
               // create new offspring, a and b, where
               // a = x*A + (1-x)*B, and b = (1-x)*A + x*B    -- note: x is alpha in the code
               alpha = (Real) ranf();
#ifdef DEBUG
               (void)fprintf(logFile, "gs.cc/  alpha = " FDFMT "\n", alpha);
               (void)fprintf(logFile, "gs.cc/ About to call crossover_arithmetic with original_population[%d] & [%d]\n", fi, mi);
#endif /* DEBUG */
               crossover_arithmetic( original_population[fi].genotyp, 
                                     original_population[mi].genotyp, 
                                     alpha, outlev, logFile );
               break;
            default:
                (void)fprintf(logFile,"gs.cc/ Unrecognized crossover mode!\n");
         }
         original_population[fi].age = 0L;
         original_population[mi].age = 0L;
         //map genotype to phenotype and evaluate energy
         original_population[fi].mapping();
         original_population[mi].mapping();
	 // update statistics
         ci_count++; // count of crossovers, individ-by-individ
      }
   }
}


void Genetic_Algorithm::crossover_2pt(Genotype &father, Genotype &mother, const unsigned int pt1, const unsigned int pt2, int outlev, FILE *logFile)
{
    // Assumes that 0<=pt1<=pt2<=number_of_pts

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover_2pt(Genotype");
   (void)fprintf(logFile, "&father, Genotype &mother, unsigned int pt1, unsigned int pt2)\n");
   (void)fprintf(logFile,"gs.cc/Performing crossover from %d to %d \n", pt1,pt2);
#endif /* DEBUG */

   // loop over genes to be crossed over
   for (unsigned int i=pt1; i<pt2; i++) {
#ifdef DEBUG
      //(void)fprintf(logFile,"gs.cc/1::At pt %d   father: %.3lf   mother: %.3lf\n",
      //i, *((double *)father.gread(i)), *((double *)mother.gread(i)) );
#endif /* DEBUG */
      Element temp = father.gread(i);
      father.write(mother.gread(i), i);
      mother.write(temp, i);
      cg_count++; // count of crossovers, gene-by-gene
#ifdef DEBUG
      //(void)fprintf(logFile,"gs.cc/1::At pt %d   father: %.3lf   mother: %.3lf\n",
      //i, *((double *)father.gread(i)), *((double *)mother.gread(i)) );
#endif /* DEBUG */
   }

#ifdef DEBUG_QUAT
   Quat q_father, q_mother;
#ifdef DEBUG_QUAT_PRINT
   pr( logFile, "DEBUG_QUAT: crossover_2pt()  q_father\n" );
   pr( logFile, "DEBUG_QUAT: crossover_2pt()  pt1=%d, pt2=%d\n", pt1, pt2 );
#endif // endif DEBUG_QUAT_PRINT
   //  Make sure the quaternion is suitable for 3D rotation
   q_father = father.readQuat();
#ifdef DEBUG_QUAT_PRINT
   printQuat( logFile, q_father );
   (void) fflush(logFile);
#endif // endif DEBUG_QUAT_PRINT
   assertQuatOK( q_father );
#ifdef DEBUG_QUAT_PRINT
   pr( logFile, "DEBUG_QUAT: crossover_2pt()  q_mother\n" );
   (void) fflush(logFile);
#endif // endif DEBUG_QUAT_PRINT
   //  Make sure the quaternion is suitable for 3D rotation
   q_mother = mother.readQuat();
#ifdef DEBUG_QUAT_PRINT
   printQuat( logFile, q_mother );
   (void) fflush(logFile);
#endif // endif DEBUG_QUAT_PRINT
   assertQuatOK( q_mother );
#endif // endif DEBUG_QUAT

}


void Genetic_Algorithm::crossover_uniform(Genotype &father, Genotype &mother, const unsigned int num_genes, int outlev, FILE *logFile)
{
#ifdef DEBUG
    (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover_uniform(Genotype");
    (void)fprintf(logFile, "&father, Genotype &mother, unsigned int num_genes)\n");
#endif /* DEBUG */

    for (unsigned int i=0; i<num_genes; i++) {
        // Choose either father's or mother's gene/rotation gene set, with a 50/50 probability
        if (ranf() > 0.5 ) continue; // 50% probability of crossing a gene or gene-set; next i
        if ( ! is_rotation_index( i ) ) {
            // Exchange parent's genes
            Element temp = father.gread(i);
            father.write(mother.gread(i), i);
            mother.write(temp, i);
        } else if ( is_first_rotation_index(i) ) {
            // don't crossover within a quaternion, only at the start; next i
#ifdef DEBUG_QUAT
            Quat q_father, q_mother;
#ifdef DEBUG_QUAT_PRINT
            pr( logFile, "DEBUG_QUAT: pre-crossover_uniform()  q_father\n" );
#endif // endif DEBUG_QUAT_PRINT
            //  Make sure the quaternion is suitable for 3D rotation
            q_father = father.readQuat();
#ifdef DEBUG_QUAT_PRINT
            printQuat( logFile, q_father );
            (void) fflush(logFile);
#endif // endif DEBUG_QUAT_PRINT
            assertQuatOK( q_father );
#ifdef DEBUG_QUAT_PRINT
            pr( logFile, "DEBUG_QUAT: pre-crossover_uniform()  q_mother\n" );
            (void) fflush(logFile);
#endif // endif DEBUG_QUAT_PRINT
            //  Make sure the quaternion is suitable for 3D rotation
            q_mother = mother.readQuat();
#ifdef DEBUG_QUAT_PRINT
            printQuat( logFile, q_mother );
            (void) fflush(logFile);
#endif // endif DEBUG_QUAT_PRINT
            assertQuatOK( q_mother );
#endif // endif DEBUG_QUAT
            // Exchange father's or mother's set of rotation genes
            for (unsigned int j=i; j<i+4; j++) {
                Element temp = father.gread(j);
                father.write(mother.gread(j), j);
                mother.write(temp, j);
            }
            // Increment gene counter, i, by 3, to skip the 3 remaining rotation genes
            i=i+3;
#ifdef DEBUG_QUAT
#ifdef DEBUG_QUAT_PRINT
            pr( logFile, "DEBUG_QUAT: post-crossover_uniform()  q_father\n" );
#endif // endif DEBUG_QUAT_PRINT
            //  Make sure the quaternion is suitable for 3D rotation
            q_father = father.readQuat();
#ifdef DEBUG_QUAT_PRINT
            printQuat( logFile, q_father );
            (void) fflush(logFile);
#endif // endif DEBUG_QUAT_PRINT
            assertQuatOK( q_father );
#ifdef DEBUG_QUAT_PRINT
            pr( logFile, "DEBUG_QUAT: post-crossover_uniform()  q_mother\n" );
            (void) fflush(logFile);
#endif // endif DEBUG_QUAT_PRINT
            //  Make sure the quaternion is suitable for 3D rotation
            q_mother = mother.readQuat();
#ifdef DEBUG_QUAT_PRINT
            printQuat( logFile, q_mother );
            (void) fflush(logFile);
#endif // endif DEBUG_QUAT_PRINT
            assertQuatOK( q_mother );
#endif // endif DEBUG_QUAT
        } // is_rotation_index( i )
    cg_count++; // count of crossovers, gene-by-gene
    } // next i
}

void Genetic_Algorithm::crossover_arithmetic(Genotype &A, Genotype &B, ConstReal  alpha, int outlev, FILE *logFile)
{
   register unsigned int i;
   Element temp_A, temp_B;
   Real one_minus_alpha;

   one_minus_alpha = 1.0 - alpha;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover_arithmetic(Genotype");
   (void)fprintf(logFile, "&A, Genotype &B, Real alpha)\n");
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover_arithmetic");
   (void)fprintf(logFile, "/Trying to perform arithmetic crossover using alpha = " FDFMT "\n", alpha);
   (void)fflush(logFile);
#endif /* DEBUG */

   // loop over genes to be crossed over
   // a = alpha*A + (1-alpha)*B
   // b = (1-alpha)*A + alpha*B
   for (i=0; i<A.num_genes(); i++) {
       if ( is_translation_index(i) ) {
#ifdef DEBUG
           (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover_arithmetic");
           (void)fprintf(logFile, "/looping over genes to be crossed over, i = %d\n", i);
           (void)fflush(logFile);
#endif /* DEBUG */
           temp_A = A.gread(i);
           temp_B = B.gread(i);
#ifdef DEBUG
           (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover_arithmetic");
           (void)fprintf(logFile, "/temp_A = %.3f  &  temp_B = %.3f\n", temp_A.real, temp_B.real);
           (void)fflush(logFile);
#endif
           A.write( (one_minus_alpha * temp_A.real  +  alpha * temp_B.real), i);
           B.write( (alpha * temp_A.real  +  one_minus_alpha * temp_B.real), i);
#ifdef DEBUG
           (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover_arithmetic");
           (void)fprintf(logFile, "/A = %.3f  &  B = %.3f\n", A.gread(i).real, B.gread(i).real);
           (void)fflush(logFile);
#endif
       } else if ( is_rotation_index(i) ) {
           // Interpolate the father's and mother's set of 4 quaternion genes
           Quat q_A;
           q_A = slerp( A.readQuat(), B.readQuat(), alpha );
           B.writeQuat( slerp( A.readQuat(), B.readQuat(), one_minus_alpha ) );
           A.writeQuat( q_A );
           // Increment gene counter, i, by 3, to skip the 3 remaining quaternion genes
           i=i+3;
       } else if ( is_conformation_index(i) ) {
           // Use anglular interpolation, alerp(a,b,fract), to generate properly interpolated torsion angles
           temp_A = A.gread(i);
           temp_B = B.gread(i);
           A.write( alerp(temp_A.real, temp_B.real, alpha), i);
           B.write( alerp(temp_A.real, temp_B.real, one_minus_alpha), i);
       } else {
           // MP: BUG CHECK!
	   char msg[100];
           (void)sprintf(msg, "BUG: Invalid gene type at i=%d\n", i);
	   stop(msg); // exits
       }
   cg_count++; // count of crossovers, gene-by-gene
   }
}

/* * */

/*
 * Proportional Selection
 *
 *
 *  We want to perform minimization on a function.  To do so, we
 *      take fitness values in a given generation and normalize them s.t.
 *      they are non-negative, and such that smaller f's are given a higher
 *      weighting.
 *  If f in [a,b], then f' in [A,B] s.t. f(a) => f'(B) and f(b) => f'(A) is
 *
 *      f' = (B-A) * (1 - (f-a)/(b-a)) + A
 *
 *         = (B-A) * (b-f)/(b-a) + A
 *
 *  Note that it suffices to map f into [0,1], giving
 *
 *      f' = (b-f)/(b-a)
 *
 *  Now the expected number of samples generated from a given point is
 *
 *      N * f'_i / \sum f'_i  =  N * ((b-f_i)/(b-a)) / \sum ((b-f_i)/(b-a))
 *
 *                            =  N * (b-f_i) / (N*b - \sum f_i)
 *
 *  Within a given generation, let 'b' be the maximal (worst) individual and 
 *      let 'a' be the minimal individual. 
 *
 *  Note:
 *      (1) This calculation is invariant to the value of B, but _not_ 
 *              invariant to the value of A.
 *      (2) This selection strategy works fine for functions bounded above,
 *              but it won't necessarily work for functions which are unbounded
 *              above.
 *
 *  The 'b' parameter is represented as 'Worst' and is selected by the 
 *    scale_fitness method.  If a value is greater than Worst, it is given a 
 *    value of zero for it's expectation.
 */

/* not static */
void Genetic_Algorithm::selection_proportional(Population &original_population, Individual *const new_pop, int outlev, FILE *logFile)
{
   register unsigned int i=0, start_index = 0;
   int temp_ordering, temp_index;
#ifdef DEBUG2
   int allzero = 1;//debug
   Molecule *individualMol;//debug
#endif

#undef CHECK_ISNAN
#ifdef CHECK_ISNAN
   int allEnergiesEqual = 1;
   double diffwa = 0.0, invdiffwa = 0.0, firstEnergy = 0.0;
#endif

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::");
   (void)fprintf(logFile, "selection_proportional(Population &original_population, Individual *new_pop)\n");
#endif /* DEBUG */

#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/At the start of sel_prop:  sel_prop_count= %d, start_index= %d\n\n",sel_prop_count, start_index); //debug

   //original_population.printPopulationAsStates(logFile, original_population.num_individuals(), global_ntor);//old debug
#endif

#ifdef CHECK_ISNAN

   /*
    * This is the new code to check for the NaN case that could
    * arise; this will preempt any fatal errors.
    */

   // Calculate expected number of children for each individual

   assert(finite(worst));
   assert(finite(avg));
   assert(!ISNAN(worst));
   assert(!ISNAN(avg));

   diffwa = worst - avg;

   assert(finite(diffwa));
   assert(!ISNAN(diffwa));

   if (diffwa != 0.0) { // added by gmm, 4-JUN-1997

      invdiffwa = 1.0 / diffwa;

      if (ISNAN(invdiffwa)) {
          (void)fprintf(logFile,"WARNING!  While doing proportional selection, not-a-number was detected (NaN).\n");
          (void)fprintf(logFile,"All members of the population will be arbitrarily allocated 1 child each.\n\n");
          for (i=0;  i < original_population.num_individuals();  i++) {
             alloc[i] = 1.0;  // arbitrary
          }
          // Assume run has converged:
          converged = 1;
          (void)fprintf(stderr,"WARNING!  The population appears to have converged, so this run will shortly terminate.\n");
      } else {
          assert(finite(invdiffwa));
          assert(finite(worst));
          assert(finite(original_population.num_individuals()));
          assert(!ISNAN(invdiffwa));
          assert(!ISNAN(worst));
          assert(!ISNAN(original_population.num_individuals()));

          for (i=0;  i < original_population.num_individuals();  i++) {
             alloc[i] = (worst - original_population[i].value(e_mode)) * invdiffwa;

             assert(finite(original_population[i].value(e_mode)));
             assert(finite(alloc[i]));
             assert(!ISNAN(original_population[i].value(e_mode)));
             assert(!ISNAN(alloc[i]));
          }// for i
      }// endif (ISNAN(invdiffwa))
   } else if (original_population.num_individuals() > 1) {
      // diffwa = 0.0,  so worst = avg
      // This implies the population may have converged.
      converged = 1;

      // Write warning message to both stderr and logFile...
      (void)fprintf(stderr,"WARNING!  While doing proportional selection, worst (%6.2le) = avg (%6.2le).\n", worst, avg);
      (void)fprintf(stderr,"          This would cause a division-by-zero error.\n");
      (void)fprintf(stderr,"          All members of the population will be arbitrarily allocated 1 child each.\n");
      (void)fprintf(stderr,"WARNING!  The population appears to have converged, so this run will shortly terminate.\n\n");

      (void)fprintf(logFile,"WARNING!  While doing proportional selection, worst (%6.2le) = avg (%6.2le).\n", worst, avg);
      (void)fprintf(logFile,"          This would cause a division-by-zero error.\n");
      (void)fprintf(logFile,"          All members of the population will be arbitrarily allocated 1 child each.\n");
      (void)fprintf(logFile,"WARNING!  The population appears to have converged, so this run will shortly terminate.\n\n");

      alloc[0] = 1.0; // Added by gmm, 2-APR-1997
      firstEnergy = original_population[0].value(e_mode);
      allEnergiesEqual = 1;
      for (i=1;  i < original_population.num_individuals();  i++) {
         alloc[i] = 1.0; // Added by gmm, 2-APR-1997
         allEnergiesEqual = allEnergiesEqual && (firstEnergy == original_population[i].value(e_mode));
      }
      if (allEnergiesEqual) {
          (void)fprintf(logFile,"          All individuals in the population have the same fitness (%6.2le)\n", firstEnergy);
      } else {
          (void)fprintf(logFile,"          Here are the fitness values of the population:\n\n");
          for (i=0;  i < original_population.num_individuals();  i++) {
              (void)fprintf(logFile,"%3d = %6.2le,  ", i+1, original_population[i].value(e_mode) );
          }
      }
      (void)fprintf(logFile,"\n");
   } // diffwa = 0.0,  so worst = avg

   /*
    * Commented out because this writes too many
    * errors out -- (worst-avg) is constant inside loop,
    * so can be brought outside for-loop.
    * for (i=0; i<original_population.num_individuals(); i++) {
    * if ((worst - avg) != 0.0) { // added by gmm, 4-JUN-1997
    * //  In function minimization, the max energy is the worst
    * alloc[i] = (worst - original_population[i].value(e_mode))/(worst - avg);
    * } else {
    * (void)fprintf(logFile,"gs.cc/WARNING!  While doing proportional selection,
    * worst (%6.2le) and avg (%6.2le) were found equal, which would cause a 
    * division-by-zero error; this was just prevented, for individual %d\n", 
    * worst, avg, i); // Added by gmm, 4-JUN-1997
    * alloc[i] = 1.0; // Added by gmm, 2-APR-1997
    * }
    * }
    */

#else 

   /*
    * This is how the code used to be, before the worst=avg problem
    * was observed.
    */

   // Calculate expected number of children for each individual
   for (i=0;  i < original_population.num_individuals();  i++)
   {
      //  In our case of function minimization, the max individual is the worst
      if(avg==worst) alloc[i] = 1./(original_population.num_individuals()+1); // HACK TODO  investigate 2008-11
       else alloc[i] = (worst - original_population[i].value(e_mode))/(worst - avg);
   }

#endif /* not CHECK_ISNAN */

#ifdef DEBUG2
   allzero = 1; //debug
   unsigned int J;//debug
   (void)fprintf(logFile, "gs.cc: checking that all alloc[] variables are not all zero...\n"); //debug
   for (J=0;  J < original_population.num_individuals();  J++) {//debug
       allzero = allzero & (alloc[J] == (Real)0.0);//debug
   }//debug
   if (allzero) {//debug
       (void)fprintf(logFile, "gs.cc:  W A R N I N G !  all alloc variables are zero!\n"); //debug
   }//debug
     (void) fprintf(logFile, "gs.cc worst = %f   avg = %f   worst-avg = %f\n", worst, avg, worst-avg);
     (void) fprintf(logFile, "gs.cc score[..] = ");
   for (J=0;  J < original_population.num_individuals() && J<20;  J++) //debug
     fprintf(logFile, "%.3f ", (worst - original_population[J].value(e_mode))*0.0000001);
     (void) fprintf(logFile, "\n");

     (void) fprintf(logFile, "gs.cc alloc[..] = ");
     double ta=0; // alloc total
   for (J=0;  J < original_population.num_individuals() && J<20;  J++) {//debug
     fprintf(logFile, "%.3f ", alloc[J]);
     ta += alloc[J];
     }
     (void) fprintf(logFile, "  total=%.3f\n", ta);
#endif

   //  Permute the individuals

#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc:  Permuting beginning: original_population.num_individuals()= %d\n", original_population.num_individuals()); //debug
#endif

   for (i=0;  i < original_population.num_individuals();  i++) {
      temp_ordering = ordering[i];

      assert(ordering[i] < original_population.num_individuals());//debug

      temp_index = ignlgi()%(original_population.num_individuals());

      assert(ordering[temp_index] < original_population.num_individuals());//debug

      ordering[i] = ordering[temp_index];
      ordering[temp_index] = temp_ordering;

#ifdef DEBUG2
      //(void)fprintf(logFile, "gs.cc: permuting in sel_prop: ordering check: ordering[i=%d]= %d, ordering[temp_index=%d]= %d\n", ordering[i],i, ordering[temp_index], temp_index);//debug
#endif
   }// endfor i

   //  We might get some savings here if we sorted the individuals before calling
   //  this routine
   for (i=0; (i < original_population.num_individuals()) && (start_index < original_population.num_individuals()); i++) {
      for (; (alloc[i] >= 1.0) && (start_index < original_population.num_individuals());  alloc[i]-= 1.0) {
         new_pop[start_index] = original_population[i];
         //new_pop[start_index].incrementAge();
         ++start_index;
      }
   }

#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::"); //debug
   (void)fprintf(logFile, "selection_proportional(Population &original_population, Individual *new_pop)\n"); //debug
#endif

   i = 0;

#ifdef DEBUG2
   int count = 0;//debug
   (void)fprintf(logFile, "gs.cc/beginning \"while(start_index < original_population.num_individuals()) {\" loop\n"); //debug
#endif

   // ??? start_index = 0; // gmm, 1998-07-13 ???

   while (start_index < original_population.num_individuals()) {
       Real r; // local 
#ifdef DEBUG2
      (void)fprintf(logFile, "gs.cc:596/inside \"while(start_index(=%d) < original_population.num_individuals()(=%d)) \" loop:  count= %d\n", start_index, original_population.num_individuals(), ++count); //debug
#endif
      assert(ordering[i] < original_population.num_individuals());//debug
      r = ranf();
#ifdef DEBUG2
      (void)fprintf(logFile, "gs.cc:599/inside debug_ranf= %.3f, alloc[ordering[i]]= %.3e, ordering[i]= %d,  i= %d\n", r, alloc[ordering[i]], ordering[i], i); // debug
#endif //  DEBUG2
      if (r < alloc[ordering[i]]) {
#ifdef DEBUG2
         (void)fprintf(logFile, "gs.cc:603/inside (debug_ranf < alloc[ordering[i]]) is true!\n"); //debug
         (void)fprintf(logFile, "gs.cc:604/inside about to increment start_index in:  \"new_pop[start_index++] = original_population[ordering[i]];\"; right now, start_index= %d\n", start_index); //debug
#endif
         new_pop[start_index] = original_population[ordering[i]];
         //new_pop[start_index].incrementAge();
         start_index++;
#ifdef DEBUG2
         (void)fprintf(logFile, "gs.cc:605/inside just incremented start_index, now start_index= %d\n", start_index); //debug
#endif
      }// endif (ranf() < alloc[ordering[i]])

#ifdef DEBUG2
      (void)fprintf(logFile, "gs.cc:608/inside i= %d, original_population.num_individuals()= %d\n", i, original_population.num_individuals()); //debug
      (void)fprintf(logFile, "gs.cc:609/inside about to \"i = (i+1)%%original_population.num_individuals();\"\n"); //debug
#endif
      i = (i+1)%original_population.num_individuals();
#ifdef DEBUG2
      (void)fprintf(logFile, "gs.cc:611/inside just done \"i = (i+1)%%original_population.num_individuals();\"\n"); //debug
      (void)fprintf(logFile, "gs.cc:612/inside i= %d  _____________________________________________________\n\n", i); //debug

       allzero = 1;//debug
       for (J=0;  J < original_population.num_individuals();  J++) {//debug
           allzero = allzero & (alloc[J] == (Real)0.0);//debug
       }//debug
       if (allzero) {//debug
           (void)fprintf(logFile, "gs.cc:  W A R N I N G !  all alloc variables are zero!\n"); //debug
       }//debug
#endif

   }// endwhile (start_index < original_population.num_individuals())

#ifdef DEBUG2
  (void)fprintf(logFile, "gs.cc/finished \"while(start_index < original_population.num_individuals()) \" loop\n"); //debug
#endif

}

/* Probabilistic Binary Tournament Selection
 *
 * This type of tournament selection was described in Goldberg and Deb (1991),
 * and performs the same type of selection as is seen in linear rank
 * selection.  Two individuals are chosen at random, and the better individual
 * is selected with probability P, 0.5 <= P <= 1.  If we let 2P = C, then
 * we see that the expected number of samples generated by rank r_i is
 *
 *      N * [2P - 2(2P-1) r_i]            , 1 >= P >= 0.5
 *
 * We can again parameterize this ranking with K, the relative probability
 * between the best and worst individual.  Since 2P = C,
 * P = K/(1+K).
 */
void Genetic_Algorithm::selection_tournament(Population &original_population, Individual *const new_pop, int outlev, FILE *logFile)
{
   register unsigned int i = 0, start_index = 0;
   int temp_ordering, temp_index;
//#define DEBUG
#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::");
   (void)fprintf(logFile, "selection_tournament(Population &original_population, Individual *new_pop)\n");
#endif /* DEBUG */
   Real tournament_prob =  0.; // Dummy value - this code seems unfinished M Pique Oct 2009
   original_population.msort(original_population.num_individuals());
   for (i=0; i<original_population.num_individuals(); i++) {
      alloc[i] = original_population.num_individuals()*(2*tournament_prob - i*(4*tournament_prob - 2));
   }

   for (i=0;  (i < original_population.num_individuals()) && (start_index < original_population.num_individuals());  i++) {
      for (; (alloc[i] >= 1.0) && (start_index < original_population.num_individuals());  alloc[i] -= 1.0) {
         new_pop[start_index++] = original_population[i];
      }
   }

   for (i=0; i < original_population.num_individuals(); i++) {
      temp_ordering = ordering[i];
      temp_index = ignlgi()%original_population.num_individuals();
      ordering[i] = ordering[temp_index];
      ordering[temp_index] = temp_ordering;
   }

   i = 0;
   while (start_index < original_population.num_individuals()) {
      if (ranf() < alloc[ordering[i]]) {
         new_pop[start_index++] = original_population[ordering[i]];
      }
      i = (i+1)%original_population.num_individuals();
   }
}
/* Linear Ranking Selection
 
  This type of selection was described in Goldberg and Deb (1991),
  Section 4 and subsection 4.2. "Sort the population from best to
  worst, assign the number of copies [alpha, or 'alloc[]' here] that
  each individual should receive ..., and then perform proportionate
  selection according to that assignment."
 
 
  We can again parameterize this ranking with K, the relative probability
  between the best and worst individual.  

  M Pique: set c0 and c1 (selection coeffient for first individual, slope)
   givens: k: ratio  k=c0/(c0-c1)   and that c1=2*(c0-1) (see text)
  
   k=c0/(c0-2*(c0-1))
   c0= k*(c0-2*(c0-1))
   c0= k*c0 -2*k*c0 +2*k
   c0 = -k*c0 + 2*k
   c0+(k*c0) =  2*k
   (1+k)*c0= 2*k
   c0 = 2*k/(1+k)
  
   ratio k is   c0 / (  c0 - c1), so k*c0 - k*c1 = c0, (k-1)*c0 = k*c1
   Goldberg & Deb : c1=2*(c0-1) so 2*c0 = c1+2

 ***/
void Genetic_Algorithm::selection_linear_ranking(/* not const (msort) */ Population &original_population,
						 /* not const */ Individual *const new_pop,
						int outlev, FILE *logFile) 
{
   register unsigned int i = 0, start_index = 0;
//#define DEBUG
#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::");
   (void)fprintf(logFile, "selection_linear_ranking(Population &original_population, Individual *new_pop)\n");
   (void)fprintf(logFile, "gs.cc/ linear_ranking_selection_probability_ratio=%f\n", linear_ranking_selection_probability_ratio);
#endif /* DEBUG */
   Real c0 =  2*linear_ranking_selection_probability_ratio/(1+linear_ranking_selection_probability_ratio); // K
   Real c1 = 2*(c0-1);
   unsigned int num_indiv = original_population.num_individuals(); // abbreviation
   original_population.msort(original_population.num_individuals());

   for (i=0; i<num_indiv; i++) {
      Real x; // range 0..1 as i ranges 0..num_indiv-1
      if(num_indiv<2) x=0;
      else x=i/(num_indiv-1);
      alloc[i] = (c0-c1*x); // alpha
   }

   // allocate whole number parts of alloc (alpha) distribution
   for (i=0;  i<num_indiv && start_index<num_indiv;  i++) {
      while (alloc[i]>=1.0 && start_index<num_indiv) {
         new_pop[start_index++] = original_population[i];
	 alloc[i] -= 1.0;
      }
   }

   // allocate fractional (residual) parts of alloc distribution
   // The method here is not completly "fair" but is quick 
   //  choose a random indiv (without replacement) and decide selection
   // from biased coin flip
   // Again see Goldberg & Deb 1991    -= M Pique October 2009
   // WHILE loop has emergency breakout for "cant happen" cases... MP
   int emergencydoor=25+num_indiv/10; // just in case
   unsigned int fracalloc_start = start_index; // for statistics logging
   while ( start_index < num_indiv && --emergencydoor>0 ) {
	i = ignlgi()%num_indiv;
	if(alloc[i]<=0) continue;
	if(alloc[i]>ranf()) {
	  new_pop[start_index++] = original_population[i];
	  }
	  alloc[i]=0; // remove from further consideration
   }
   if(start_index!=num_indiv) {
     // put in log
   (void)fprintf(logFile, "WARNING: gs.cc: linear_ranking_selection found no indiv for %d slot%s of %d.\n", num_indiv-start_index, (num_indiv-start_index)==1?"":"s", num_indiv-fracalloc_start);
     // just copy most fit indiv into rest of array
   while(start_index < num_indiv) new_pop[start_index++] = original_population[0];
   }
}

Individual *Genetic_Algorithm::selection(Population &solutions, int outlev, FILE *logFile)
{
   Individual *next_generation;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/Individual *Genetic_Algorithm::selection(Population &solutions)\n");
#endif /* DEBUG */

   next_generation = new Individual[solutions.num_individuals()];
   
   set_worst(solutions, outlev, logFile);
   switch(s_mode)
   {
      case Proportional:
         selection_proportional(solutions, next_generation, outlev, logFile);
         break;
      case LinearRanking:
         selection_linear_ranking(solutions, next_generation, outlev, logFile);
         break;
      case Tournament:
         // M Pique October 2009 - does not work so disallowing for now
	 // selection_tournament(solutions, next_generation);
         (void)fprintf(logFile,"gs.cc/Unimplemented Selection Method - using proportional\n");
         selection_proportional(solutions, next_generation, outlev, logFile);
         break;
      case Boltzmann:
         (void)fprintf(logFile,"gs.cc/Unimplemented Selection Method - using proportional\n");
         selection_proportional(solutions, next_generation, outlev, logFile);
         break;
      default:
         (void)fprintf(logFile,"gs.cc/Unknown Selection Mode!\n");
   }

   return(next_generation);
}

//  Global search is a GA or a PSO
//
//  This is where the action is... SEARCH!
//
int Genetic_Algorithm::search(Population &solutions, Eval *evaluate, int outlev, FILE *logFile)
{
   register unsigned int i;

#ifdef DEBUG /* DEBUG { */
   (void)fprintf(logFile, "gs.cc/int Genetic_Algorithm::search(Population &solutions)\n");
#endif /* } DEBUG */

#ifdef DEBUG3 /* DEBUG3 { */
   for (i=0; i<solutions.num_individuals(); i++) {
       (void)fprintf(logFile,"%ld ", solutions[i].age);
   }
   (void)fprintf(logFile,"\n");
#endif /* } DEBUG3 */

   // search no longer doing Mapping on the incoming population
   // assumed already mapped earlier, may 2009 mp+rh
   
#ifdef DEBUG3 /* DEBUG3 { */
   (void)fprintf(logFile,"About to perform Selection on the solutions.\n");
   for (i=0; i<solutions.num_individuals(); i++) {
       (void)fprintf(logFile,"%ld ", solutions[i].age);
   }
   (void)fprintf(logFile,"\n");
#endif /* } DEBUG3 */

   //
   // Perform selection
   //
   Population newPop(solutions.num_individuals(), solutions.evaluate, selection(solutions, outlev, logFile));

#ifdef DEBUG3 /* DEBUG3 { */
   (void)fprintf(logFile,"About to perform Crossover on the population, newPop.\n");
   for (i=0; i<solutions.num_individuals(); i++) {
       (void)fprintf(logFile,"%ld ", newPop[i].age);
   }
   (void)fprintf(logFile,"\n");
#endif /* } DEBUG3 */

   //
   // Perform crossover
   // 
   crossover(newPop, outlev, logFile);

#ifdef DEBUG3 /* DEBUG3 { */
   (void)fprintf(logFile,"About to perform mutation on the population, newPop.\n");
   for (i=0; i<solutions.num_individuals(); i++) {
       (void)fprintf(logFile,"%ld ", newPop[i].age);
   }
   (void)fprintf(logFile,"\n");
#endif /* } DEBUG3 */

   //
   // Perform mutation
   // 
   mutation(newPop, outlev, logFile);

#ifdef DEBUG3 /* DEBUG3 { */
   (void)fprintf(logFile,"About to perform elitism, newPop.\n");
   for (i=0; i<solutions.num_individuals(); i++) {
       (void)fprintf(logFile,"%ld ", newPop[i].age);
   }
   (void)fprintf(logFile,"\n");
#endif /* } DEBUG3 */

   //
   // Copy the n best individuals to the next population, if the elitist flag is set, where n is the value of elitism.
   //
   if (elitism > 0) {
      solutions.msort(elitism);
      for (i=0; i<elitism; i++) {
         newPop[solutions.num_individuals()-1-i] = solutions[i];
      }
   }

#ifdef DEBUG3 /* DEBUG3 { */
   (void)fprintf(logFile,"About to Update the Current Generation, newPop.\n");
   for (i=0; i<solutions.num_individuals(); i++) {
       (void)fprintf(logFile,"%ld ", newPop[i].age);
   }
   (void)fprintf(logFile,"\n");
#endif /* } DEBUG3 */

   //
   // Update current generation 
   //  
   solutions = newPop;

   //
   // Increase the number of generations
   //
   generations++;

   //
   // Increment the age of surviving individuals...
   // 
   for (i=0; i<solutions.num_individuals(); i++) {
       solutions[i].incrementAge();
   }
   return(0);
}
int Genetic_Algorithm::localsearch(Population &thisPop, Local_Search *local_method, 
 Eval *evaluate, int outlev, FILE *logFile)
{
	// apply local search (if not NULL) to population with 
	// frequency search_freq
	// (this code moved from call_glss.cc into methods of
	//  the GA/LGA and PSO since their needs differed so - M Pique Jun 2011)
	if(local_method != NULL) for (unsigned int i=0; i<thisPop.num_individuals(); i++) {
#ifdef LOCALSEARCHDEBUG
// MP June 2011 - disabled after move of code TODO since values not handy
            if (outlev >= LOGRUNV) {
                (void)fprintf( logFile, "LS: %d",generations-1); 
                (void)fprintf( logFile, " %d",i+1); 
                (void)fprintf( logFile, " %f",thisPop[i].value(localEvalMode)); 
           }
#endif
           if(ranf() < localsearch_freq ) local_method->search(thisPop[i], 
	    evaluate, outlev, logFile);
#ifdef LOCALSEARCHDEBUG
           if (outlev >= LOGRUNV) {
                (void)fprintf( logFile, " %f",thisPop[i].value(localEvalMode)); 
                (void)fprintf( logFile, " \n"); 
            }
#endif
	} // if a local_method is passed
   return(0);
}
