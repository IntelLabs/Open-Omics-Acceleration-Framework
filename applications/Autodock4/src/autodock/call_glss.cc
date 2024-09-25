/*

 $Id: call_glss.cc,v 1.75 2014/06/12 01:44:07 mp Exp $ 
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/********************************************************************
     Call_glss:  Invokes either GA-only or a GA-LS hybrid to try and solve the
                 docking problem.

                                rsh 9/95
********************************************************************/

#include <string.h>
#include "gs.h"
#include "ls.h"
#include "support.h"
#include "eval.h"
#include "hybrids.h"
#include "constants.h"
#include "structs.h"
#include "openfile.h"
#include "qmultiply.h"

#include <sys/types.h>      /*time_t time(time_t *tloc); */
#include <time.h>           /*time_t time(time_t *tloc); */
#include "timesyshms.h"

extern char *programname;
extern int debug;


Representation **generate_R(const int num_torsions, GridMapSetInfo *const info, FILE *logFile)
{
   Representation **retval;
   Quat q;

#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R()  about to create a new Representation with 5 elements, retval...\n");
#endif
   retval = new Representation *[5];
#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R()  done creating   a new Representation with 5 elements, retval...\n");
#endif
   // Set the x-translation
   retval[0] = new RealVector( 1, info->lo[X], info->hi[X] );
   // Set the y-translation
   retval[1] = new RealVector( 1, info->lo[Y], info->hi[Y] );
   // Set the z-translation
   retval[2] = new RealVector( 1, info->lo[Z], info->hi[Z] );

   // Generate a uniformly-distributed random quaternion for a random rotation (UDQ)
   q = randomQuat();
#ifdef DEBUG
   printQuat( logFile, q );
#endif

   // Set the unit vector components (the "axis"), for the rotation about axis
   retval[3] = new RealVector( 4, -1., 1., q.x, q.y, q.z, q.w ); // uniformly-distributed quaternion (UDQ)

   // Set the angle for any torsion angles
   retval[4] = new RealVector( num_torsions, -PI, PI );  // torsion angle is uniformly distributed, -PI to PI

#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R()  done assigning each of the retval[0-4] elements...\n");
#endif

   return(retval);
}

Representation **generate_R_quaternion(const int num_torsions, const GridMapSetInfo *const info, FILE *logFile)
{
   Representation **retval;
   Quat q;

#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/Representation **generate_R_quaternion()  about to create a new Representation with 5 elements, retval...\n");
#endif
   retval = new Representation *[5];
#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R_quaternion()  done creating   a new Representation with 5 elements, retval...\n");
#endif
   // Set the x-translation
   retval[0] = new RealVector( 1, info->lo[X], info->hi[X] );
   // Set the y-translation
   retval[1] = new RealVector( 1, info->lo[Y], info->hi[Y] );
   // Set the z-translation
   retval[2] = new RealVector( 1, info->lo[Z], info->hi[Z] );

   // Generate a uniformly-distributed random quaternion for a random rotation (UDQ)
   q = randomQuat();
#ifdef DEBUG
   printQuat( logFile, q );
#endif

#ifdef DEBUG_QUAT
#ifdef DEBUG_QUAT_PRINT
    pr( logFile, "DEBUG_QUAT: generate_R_quaternion()\n" );
    (void) fflush(logFile);
#endif
    //  Make sure the quaternion is suitable for 3D rotation
    assertQuatOK( q );
#endif

   // Set the quaternion (x,y,z,w) genes
   retval[3] = new RealVector( 4, -1., 1., q.x, q.y, q.z, q.w ); // uniformly-distributed quaternion (UDQ)
   // TODO retval[3] = new ConstrainedRealVector( 4, -1., 1., q.x, q.y, q.z, q.w ); // uniformly-distributed quaternion (UDQ)

   // Set the torsion angles
   retval[4] = new RealVector( num_torsions, -PI, PI );

#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R_quaternion()  done assigning each of the retval[0-5] elements...\n\n");
#endif

   return(retval);
}

Genotype generate_Gtype(const int num_torsions, const GridMapSetInfo *const info, int outlev, FILE *logFile)
{
#ifdef DEBUG
    // (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Gtype() about to call Genotype temp(5, generate_R())...\n");
    (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Gtype() about to call Genotype temp(5, generate_R_quaternion())...\n");
#endif
   // Genotype temp((unsigned int)5, generate_R(num_torsions, info));
   Genotype temp((unsigned int)5, generate_R_quaternion(num_torsions, info, logFile));
#ifdef DEBUG
   // (void)fprintf(logFile,"call_glss.cc/Genotype generate_Gtype() done calling  Genotype temp(5, generate_R())...\n\n");
   (void)fprintf(logFile,"call_glss.cc/Genotype generate_Gtype() done calling  Genotype temp(5, generate_R_quaternion())...\n\n");
#endif

   return(temp);
}

Phenotype generate_Ptype(const int num_torsions, const GridMapSetInfo *const info, Eval *evaluate, int outlev, FILE *logFile) 
{
#ifdef DEBUG
    // (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Ptype() about to call Phenotype temp(5, generate_R())...\n");
    (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Ptype() about to call Phenotype temp(5, generate_R_quaternion())...\n");
#endif
   // Phenotype temp((unsigned int)5, generate_R(num_torsions, info));
   Phenotype temp((unsigned int)5, generate_R_quaternion(num_torsions, info, logFile), evaluate);
#ifdef DEBUG
   // (void)fprintf(logFile,"call_glss.cc/Genotype generate_Ptype() done calling  Phenotype temp(5, generate_R())...\n\n");
   (void)fprintf(logFile,"call_glss.cc/Genotype generate_Ptype() done calling  Phenotype temp(5, generate_R_quaternion())...\n\n");
#endif

   return(temp);
}

Individual random_ind(const int num_torsions,  const GridMapSetInfo *const info, Eval *evaluate, int outlev, FILE *logFile) 
{

#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/Individual random_ind()  About to generate_Gtype()...\n");
#endif
   Genotype temp_Gtype = generate_Gtype(num_torsions, info, outlev, logFile);
#ifdef DEBUG
   (void)fprintf(logFile,"call_glss.cc/Individual random_ind()  About to generate_Ptype()...\n");
#endif
   Phenotype temp_Ptype = generate_Ptype(num_torsions, info, evaluate, outlev, logFile); 

#ifdef DEBUG
   (void)fprintf(logFile,"call_glss.cc/Individual random_ind()  About to Individual temp(temp_Gtype, temp_Ptype)...\n");
#endif
   //shotgun wedding: does not map genotype to phenotype
   Individual temp(temp_Gtype, temp_Ptype);
   temp.mapping();

#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Individual random_ind()  Done     Individual temp(temp_Gtype, temp_Ptype)...\n\n");
#endif

   return(temp);
}

#ifdef MOVEDTOCONFORMATIONSAMPLER
// this block moved to conformation_sampler.cc  
Individual set_ind(const int num_torsions,  const GridMapSetInfo *const info, const State state, FILE *logFile)
{
   Genotype temp_Gtype;
   Phenotype temp_Ptype;
   Quat q;
   int i;

   temp_Gtype = generate_Gtype(num_torsions, info, logFile);
   temp_Ptype = generate_Ptype(num_torsions, info, logFile);

   // use the state to generate a Genotype
   temp_Gtype.write(state.T.x, 0);
   temp_Gtype.write(state.T.y, 1);
   temp_Gtype.write(state.T.z, 2);

#ifdef DEBUG_QUAT
#ifdef DEBUG_QUAT_PRINT
    pr( logFile, "DEBUG_QUAT: set_ind()\n" );
    (void) fflush(logFile);
#endif
    //  Make sure the quaternion is suitable for 3D rotation
    assertQuatOK( q );
#endif

   temp_Gtype.write( q.x, 3);
   temp_Gtype.write( q.y, 4);
   temp_Gtype.write( q.z, 5);
   temp_Gtype.write( q.w, 6);

   for (i=0;i<state.ntor; i++) {
       temp_Gtype.write(state.tor[i], 7+i);
   };

   Individual temp(temp_Gtype, temp_Ptype);   

   // use mapping to generate a Phenotype
   //temp.phenotyp =  temp.mapping();
   temp.mapping();

   return(temp);
}
#endif

State call_glss(/* not const */ Global_Search *global_method, 
                /* not const */ Local_Search *local_method, 
                const State& sInit, 
                const unsigned int num_evals, const unsigned int pop_size, 
                const int outlev,  FILE *logFile,
		const Output_pop_stats& output_pop_stats,
                Molecule * const mol, 
		Eval *evaluate,
                const Boole B_RandomTran0, const Boole B_RandomQuat0, const Boole B_RandomDihe0,
                const GridMapSetInfo *info, const char *const FN_pop_file,
                /* not const */ int end_of_branch[MAX_TORS])
{
    register unsigned int i;
    register int j;
    int allEnergiesEqual = 1, numTries = 0;
    int indiv = 0; // Number of Individual in Population to set initial state variables for.
    int max_numTries = 1000;
    double firstEnergy = 0.0;
    EvalMode localEvalMode = Normal_Eval;
    static FILE *pop_fileptr = NULL;


    if(outlev>=LOGRUNV) \
    (void)fprintf(logFile, "call_glss:  global_method %s   local_method %s\n",
      global_method?global_method->shortname():"NULL",
      local_method?local_method->shortname():"NULL");

    if (global_method) global_method->reset(output_pop_stats);
    if (local_method) local_method->reset();
    evaluate->reset();

    if(outlev>=LOGRUNV) \
    (void)fprintf( logFile, "Creating an initial population of %u individuals.\n", pop_size);
    Population thisPop(pop_size, evaluate);
    //  Pass in the end_of_branch tree for Branch Crossover Mode.
    thisPop.set_eob( end_of_branch );
    thisPop.nevals_last_pop_stats = 0;  // reset last time stats were printed

#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/State call_glss():  {\n");
#endif
    do {
        ++numTries;
        // Create a population of pop_size random individuals...
        for (i=0; i<pop_size; i++) {
#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/State call_glss():  Creating individual thisPop[i=%d] using random_ind(%d,info)...\n", i, sInit.ntor);
#endif
            thisPop[i] = random_ind( sInit.ntor, info, evaluate, outlev, logFile);
#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/State call_glss(): Created  individual i= %d in thisPop[i]\n\n", i);
#endif
            thisPop[i].mol = mol;
            thisPop[i].age = 0L;
        }

        // If initial values were supplied, put them in thisPop[0] and remap
        if (!B_RandomTran0) {
            if (outlev >= LOGRUNVV)  (void)fprintf(logFile, "Setting the initial translation (tran0) for individual number %d to %.2lf %.2lf %.2lf\n\n", indiv+1, sInit.T.x, sInit.T.y, sInit.T.z); 
            thisPop[indiv].genotyp.write( sInit.T.x, 0 );
            thisPop[indiv].genotyp.write( sInit.T.y, 1 );
            thisPop[indiv].genotyp.write( sInit.T.z, 2 );
            // Remember to keep the phenotype up-to-date
            thisPop[indiv].mapping();
        }
        if (!B_RandomQuat0) {
	        if (outlev >= LOGRUNVV) {
		AxisAngle aa = QuatToAxisAngle(sInit.Q);
                (void)fprintf(logFile,  
		 "Setting the initial orientation using quaterion values (x,y,z,w) for individual number %d to %.6lf %.6lf %.6lf %.6lf\n\n", 
		 indiv+1, sInit.Q.x, sInit.Q.y, sInit.Q.z, sInit.Q.w);
                (void) fprintf(logFile, 
		"which corresponds to the axis-angle (x,y,z,degree) values:  %.2lf %.2lf %.2lf %.2lf\n\n",
		aa.nx, aa.ny, aa.nz, aa.ang); 
            }
            thisPop[indiv].genotyp.write( sInit.Q.x, 3 );
            thisPop[indiv].genotyp.write( sInit.Q.y, 4 );
            thisPop[indiv].genotyp.write( sInit.Q.z, 5 );
            thisPop[indiv].genotyp.write( sInit.Q.w, 6 );
            // Remember to keep the phenotype up-to-date
            thisPop[indiv].mapping();
        }
        if (sInit.ntor > 0) {
            if (!B_RandomDihe0) {
                if (outlev >= LOGRUNVV) 
                (void)fprintf(logFile, "Setting the initial torsions (dihe0) for individual number %d to ", indiv+1); 
                for (j=0; j<sInit.ntor; j++) {
                    thisPop[indiv].genotyp.write( sInit.tor[j], 7+j );
                    if (outlev >=LOGRUNVV)  (void)fprintf(logFile, "%.2lf ", RadiansToDegrees(sInit.tor[j])); 
                };
                if (outlev >= LOGRUNVV) (void)fprintf(logFile, " deg\n\n");
                // Remember to keep the phenotype up-to-date
                thisPop[indiv].mapping();
            }
        }

#ifdef DEBUG
    (void)fprintf(logFile,"\n\ncall_glss.cc  // ensuring there is variation in the energies\n\n");
#endif
        // Now ensure that there is some variation in the energies...
        firstEnergy = thisPop[0].value(localEvalMode);
#ifdef DEBUG
    (void)fprintf(logFile,"\n\ncall_glss.cc  // ensuring there is variation in the energies, firstEnergy=%lf\n\n", firstEnergy);
#endif
        for (i=1; i<pop_size; i++) {
#ifdef DEBUG
    (void)fprintf(logFile,"\n\ncall_glss.cc  // ensuring there is variation in the energies, i=%d, thisPop[i].value=%lf\n\n", i, thisPop[i].value(localEvalMode));
#endif
             allEnergiesEqual = allEnergiesEqual && (thisPop[i].value(localEvalMode) == firstEnergy);
        }
        if ( pop_size>1 && allEnergiesEqual) {
            (void)fprintf(logFile,"NOTE: All energies are equal in population; re-initializing. (Try Number %d)\n", numTries);
        }
        if (numTries > max_numTries) {
            (void)fprintf(logFile,"WARNING: the number of tries has exceeded the maximum number of tries permitted.\nWARNING: AutoDock will attempt continue with the currently-generated random population.\n");
            break;
        }
    } while (pop_size>1 && allEnergiesEqual);


#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/State call_glss():  }\n");
    if (outlev > 2) { 
    thisPop.printPopulationAsCoordsEnergies( logFile, pop_size, sInit.ntor );
    }
#endif

    if (outlev >= LOGRUNVV ) { 
        (void)fprintf( logFile, "The initial population consists of the following %d individuals:\n\n", pop_size);
        (void)fprintf( logFile, "<generation t=\"%d\" after_performing=\"initialisation of population\">\n", 0);
        (void)fprintf( logFile, "</generation>\n\n\n");
    }


// We now have a mapped and evaluated population suitable for search

    // next line will resemble "Beginning LAMARCKIAN GENETIC ALGORITHM (LGA), with .."
    if(outlev >= LOGFORADT )
    (void)fprintf( logFile, "Beginning %s%s (%s%s), with a maximum of %u energy evaluations.\n\n", 
        local_method?"LAMARCKIAN ":"",
      global_method?global_method->longname():"NULL",
        local_method?"L":"",
      global_method?global_method->shortname():"NULL",
    num_evals);


     // major loop over generations - terminated by logic within
     // generation 0 is searchless and is used to print initial population statistics
     //
     // skeleton is:
     //   while ( ! terminate ) increment generation 0 to ...
     //      if ( generation > 0 )  {  search population globally and locally }
     //      if ( search is terminating  or   generation is "due" ) print stats
     //      if ( search is terminating ) break
     //      end while

     Boole terminate = FALSE;
     for(int generation=0; ! terminate ; generation++) {

	struct tms tms_genStart;
	struct tms tms_genEnd;
	Clock genEnd;
	Clock genStart = times( &tms_genStart );

        if (outlev >= LOGRUNVV ) (void)fprintf( logFile, "Global-Local Search Iteration: %d\n", generation);

      if(generation>0) {
        if (outlev >= LOGRUNVV)  (void)fprintf( logFile, "Performing Global Search.\n"); 
        global_method->search(thisPop, evaluate, outlev, logFile);

        if (outlev >= LOGRUNVVV) {
            (void)fprintf( logFile, "<generation t=\"%d\" after_performing=\"global search\">\n", generation);
            thisPop.printPopulationAsStates( logFile, pop_size, sInit.ntor );
            (void)fprintf( logFile, "</generation>\n\n\n");
        }

        if (pop_size > 1 && outlev >= LOGRUNVVV) minmeanmax( logFile, thisPop, generation, info );

        // call the global method's local search method if any
	if(local_method != NULL) {
	   if (outlev >= LOGRUNVV ) (void)fprintf( logFile, "Performing Local Search.\n");
	   if (outlev >= LOGRUNVVV ) for (i=0; i<pop_size; i++) {
                (void)fprintf( logFile, "LS: %d",generation); 
                (void)fprintf( logFile, " %d",i+1); 
                (void)fprintf( logFile, " %f",thisPop[i].value(localEvalMode)); 
           }
           global_method->localsearch(thisPop, local_method, evaluate, outlev, logFile);
	   if (outlev >= LOGRUNVVV ) for (i=0; i<pop_size; i++) {
                (void)fprintf( logFile, " %f",thisPop[i].value(localEvalMode)); 
                (void)fprintf( logFile, " \n"); 
            }

          if (outlev >= LOGRUNVV ) {
            (void)fprintf( logFile, "<generation t=\"%d\" after_performing=\"local search\">\n", generation);
            thisPop.printPopulationAsStates( logFile, pop_size, sInit.ntor );
            (void)fprintf( logFile, "</generation>\n\n\n");
          }
	} // if a local_method is active

        if (pop_size > 1 && outlev >= LOGRUNVVV ) minmeanmax( logFile, thisPop, generation, info );
	} // end if generation > 0

    // note we terminate without searching if num_evals is 0
    // current global methods' "terminate()" are generation count or convergence based
    terminate = global_method->terminate()  || evaluate->evals() >= num_evals;

    genEnd = times( &tms_genEnd );


    if(pop_size>0 && generation==0 && outlev >= LOGRUNV) {
	// print initial best energy for statistics
	// M Pique 28 Oct 2009 - adding a (harmless...) thisPop.msort(1) here
	// broke the unbound-extended test. Investigating why
	// Conclusion: the test is sensitive to population order since it does
	// only one generation of glss.  So I'm removing the msort(1) and doing
	// the search for lowest-energy individual by hand.
       double bestenergy = thisPop[0].value(Normal_Eval);
       for (i=1; i<pop_size; i++) if(bestenergy>thisPop[i].value(Normal_Eval)) \
         bestenergy=thisPop[i].value(Normal_Eval);
       (void)fprintf(logFile,"Initial-Value: %.3f\n", bestenergy);
       }


       // Print basic generation statistics 
       // (code moved here from gs.cc search() May 2011 MP)

      if (debug > 0) {
       (void)fprintf(logFile,
        "DEBUG:  Generation: %3u, output_pop_stats.everyNgens = %3d, generation%%output_pop_stats.everyNgens = %u\n",
       generation, output_pop_stats.everyNgens,
       output_pop_stats.everyNgens>0?generation%output_pop_stats.everyNgens:0);
     }
       /* Only output statistics if the output level is not 0. */
       if (output_pop_stats.everyNgens != 0 
         && generation%output_pop_stats.everyNgens == 0) {


         // print "Generation:" line (basic info, no mean/median/stddev...
         (void)fprintf(logFile,"Generation: %3u   ", generation);
#ifdef DEBUG3
         // medium (with age/pop info) output level, no newline at end
         (void) thisPop.printPopulationStatistics(logFile, 2, "");
#else
         // lowest output level, no newline at end
         (void) thisPop.printPopulationStatistics(logFile, 1, ""); 
#endif /* DEBUG3 */
         (void)fprintf(logFile,"    Num.evals.: %d   Timing: ", 
               evaluate->evals() );
         timesyshms( genEnd - genStart, &tms_genStart, &tms_genEnd, logFile);
       }

       // Print extended generational population statistics, when "due" :
       //
       //  at generation 0
       //  or..
       //  upon search termination
       //  or..
       //  every generation 0 to 20
       //  every tenth generation 20 to output_pop_stats.everyNgens (typically 100)
       //  every "output_pop_stats.everyNgens" (typically 100) if greater than that
       //
       //  or..
       //  every "output_pop_stats.everyNevals" (typically 1,000,000)
       // The number of evals when last printed is kept in 'thisPop'.
       //
       // This is purely for studying population convergence and is
       // controlled by the "output_population_statistics" DPF keyword. - M Pique 2010
       if (outlev>= LOGRUNV && output_pop_stats.level > 0 &&  (
	  generation==0
	  ||
	  terminate 
	  ||
            (output_pop_stats.everyNgens != 0 && (
	    (generation <= 20) ||
	    (generation <= (int)output_pop_stats.everyNgens 
	       && generation%10 == 0) ||
            (generation%output_pop_stats.everyNgens == 0 ) ) )
	  || (output_pop_stats.everyNevals != 0 &&
	     evaluate->evals() - thisPop.nevals_last_pop_stats 
	       >= output_pop_stats.everyNevals)
             )
	 ) {
	       // print "Population at Generation:" line 
	       //   with low/high/mean/median/stddev/state_of_best_indiv...
	       // followed by GA global search stats (crossover count, mutation count)
	       // and local search stats (invocation count)
	       (void) thisPop.printPopulationStatisticsVerbose(logFile, 
		 generation, evaluate->evals(), sInit.ntor, "");
		 if (global_method) {
		     fprintf(logFile, " cg_count: %u", global_method->cg_count);
		     fprintf(logFile, " ci_count: %u", global_method->ci_count);
		     fprintf(logFile, " mg_count: %u", global_method->mg_count);
		     fprintf(logFile, " mi_count: %u", global_method->mi_count);
		    }
		 if (local_method) fprintf(logFile, " ls_count: %u", local_method->ls_count);
		 fprintf(logFile, "\n");
		 thisPop.nevals_last_pop_stats = evaluate->evals();
		 }

        if (strlen(FN_pop_file) > 0) { // YES, do print!
            if (pop_fileptr==NULL) {
	       // attempt open 
	       if ((pop_fileptr = ad_fopen( FN_pop_file, "w", logFile)) == NULL) {
	        char msg[PATH_MAX+200];
                sprintf(msg, "%s: ERROR:  I'm sorry, I cannot create\"%s\".", programname, FN_pop_file);
		stop(msg);
		} else {
		fprintf(logFile, "Opened population file \"%s\" for writing.\n",FN_pop_file);
		fprintf(pop_fileptr, "Opened from call_glss\n");
		fflush(pop_fileptr);
		}

            } else {
		fprintf(pop_fileptr, "Generation %d\n", generation);
                fflush( pop_fileptr ); // debug
                // MP TODO FAILS 2012-08-17 thisPop.printPopulationAsCoordsEnergies( pop_fileptr, pop_size, sInit.ntor); 
                fflush( pop_fileptr );
            }
        }

        (void)fflush(logFile);
    } // end while not terminating loop over generation

    if (pop_size>1) thisPop.msort(1);
    (void)fprintf(logFile,"Final-Value: %.3f\n", thisPop[0].value(Normal_Eval));

    return( thisPop[0].state(sInit.ntor) );
}
