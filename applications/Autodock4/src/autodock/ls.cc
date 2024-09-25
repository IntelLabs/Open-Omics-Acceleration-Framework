/*

 $Id: ls.cc,v 1.23 2014/06/12 01:44:07 mp Exp $

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
//#define DEBUG0  // mild
//#define DEBUG
//#define DEBUG2
/********************************************************************
      These are the methods of the local searchers

				rsh 9/95

      Modifications to the class heirarchy made 2/21/96 rsh
********************************************************************/
#include "ls.h"

//  This function adds sign * (deviates + bias) to all the reals in the representation
Phenotype genPh(const Phenotype &original, ConstReal sign, const Real *const deviates, Real *const bias, int outlev, FILE *logFile)
{
   RepType genetype;
   register unsigned int i, index = 0;
   Phenotype retval(original);

#ifdef DEBUG2
//   (void)fprintf(logFile, "ls.cc/Phenotype genPh(const Phenotype &original, Real *deviates, Real *bias)\n");
#endif /* DEBUG */

   for (i=0; i < retval.num_pts(); i++) {
      genetype = retval.gtype(i);
      if ((genetype == T_RealV)||(genetype == T_CRealV)) {
         // 4/2009 experiment with gene-by-gene scaling,mp+rh
         //if(index>=0 && index<=2) scale = 1;//hack translation
         //else if(index>=3 && index<=6) scale = QSCALE;//hack quaternion
         //else scale = 1;//hack torsion

         retval.write(retval.gread(i).real + sign * (deviates[index] + bias[index]), i);
         index++;
      }
   }

   Quat q;
   q = retval.readQuat();

#ifdef DEBUG_QUAT
   pr( logFile, "DEBUG_QUAT: genPh()  q\n" );
   printQuat( logFile, q );
   assertQuatOK( q );
#endif // endif DEBUG_QUAT

   retval.writeQuat( normQuat( q ) );

   return(retval);
}

//  What Solis & Wets does is add random deviates to every
//  real number in the Phenotype.
//  
//  This has only one value of rho, for all genes.
// SW mapped individual
// SW returns TRUE if it modifies the individual
Boole Solis_Wets::SW(/* not const */ Phenotype &vector, int outlev, FILE *logFile)
{
   register unsigned int i, j, num_successes = 0, num_failures = 0;
   register Real temp_rho = rho;
   Phenotype newPh;
   Boole modified = FALSE;
   Quat thisQuat, prevQuat, startQuat;
   
#ifdef DEBUG
   (void)fprintf(logFile, "ls.cc/void Solis_Wets::SW(Phenotype &vector)\n");

#endif /* DEBUG */

   //  Reset bias
   for (i=0; i < size; i++) {
      bias[i] = 0.0;
   }

   for (i=0; i < max_its; i++) {
#ifdef TRACELS
//convenience function for debugging
#define traceState(msg,vector) printDState(logFile,msg,vector,i,prevxyz,startxyz,prevQuat,startQuat,num_successes,num_failures,temp_rho,bias,deviates)
void printDState(FILE *const logFile, const char *const msg, Phenotype &newPh, const int i, const Real prevxyz[3],
                 const Real startxyz[3], const Quat& prevQuat, const Quat& startQuat, 
                 const unsigned int num_successes, const unsigned int num_failures,
                 ConstReal temp_rho, Real *const  bias, Real *const  deviates);
   Real xyz[3];
   Real prevxyz[3];
   Real startxyz[3];
   //save previous values of state
   if (i>0) {
      prevQuat = thisQuat;
      for ( int d=0;d<3;d++)prevxyz[d] = xyz[d];
   }

   //update current values of state
   thisQuat = vector.readQuat();//the individual 'vector' holds where we are now
   for ( int d=0;d<3;d++) xyz[d] = vector.gread(d).real;

   // special case for first iteration
   if (i==0) {
      startQuat = thisQuat;
      for ( int d=0;d<3;d++) startxyz[d] = xyz[d];
      prevQuat = thisQuat;
      for ( int d=0;d<3;d++) prevxyz[d] = xyz[d];
      traceState("initial", vector);
   } 
   
#endif /* TRACELS */
#ifdef DEBUG
#ifdef INTERNAL
   Real dt; // translation step scalar
   Real ct; // cumulative translation from step 0
   static Quat thisQuat, prevQuat, startQuat; // for DEBUG only
   dt=0;
   for ( int d=0;d<3;d++) dt += (prevxyz[d]- xyz[d])*(prevxyz[d]- xyz[d]);
   dt = sqrt(dt);
   ct=0;
   for ( int d=0;d<3;d++) ct += (startxyz[d] - xyz[d])*(startxyz[d]- xyz[d]);
   ct = sqrt(ct);
   (void)fprintf(logFile, "\nLS::    %3d #S=%3d #F=%3d %+8.4f p=%4.2f b=(%5.2f %5.2f %5.2f) dev=(%5.2f %5.2f %5.2f)", 
                                  i, num_successes, num_failures, vector.evaluate(Normal_Eval), temp_rho, bias[0], bias[1], bias[2], 
                                  deviates[0], deviates[1], deviates[2]);
   (void)fprintf(logFile, " xyz=(");
    fprintf(logFile, "%5.2f %5.2f %5.2f", newPh.gread(0).real, newPh.gread(1).real,newPh.gread(2).real);
   (void)fprintf(logFile, ")");
    fprintf(logFile, "dT=%5.2f ", dt);
    fprintf(logFile, "cT=%5.2f ", ct);
    // assuming x,y,z,w ignoring structs.h order 
    //fprintf(logFile, "%5.2f %5.2f %5.2f %5.2f",
    thisQuat.x = newPh.gread(3).real; 
    thisQuat.y = newPh.gread(4).real;
    thisQuat.z = newPh.gread(5).real;
    thisQuat.w = newPh.gread(6).real;
    
   (void)fprintf(logFile, " quat=(");
    fprintf(logFile, "%5.2f %5.2f %5.2f %5.2f", newPh.gread(3).real, newPh.gread(4).real,newPh.gread(5).real,newPh.gread(6).real);
   //?? newPh.printIndividualsState(logFile, 7, 3);
   (void)fprintf(logFile, ")");
   if (i>0){
    fprintf(logFile, " dQ=%5.2lf", quatDifferenceToAngleDeg(prevQuat,thisQuat)); 
    fprintf(logFile, " cQ=%5.2lf", quatDifferenceToAngleDeg(startQuat,thisQuat)); 
   } else {
    startQuat = thisQuat;
   };
   prevQuat = thisQuat;
#endif /* INTERNAL */
#endif /* DEBUG */

      // Generate deviates
      for (j=0; j < size; j++) {
         deviates[j] = gen_deviates(temp_rho);
      }

      // zeta = x + bias + deviates
      newPh = genPh(vector, +1., deviates, bias, outlev, logFile); // zeta; +1 means 'forward' step
      // Evaluate
      if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
         num_successes++;
         num_failures = 0;
         vector = newPh;
         modified = TRUE;
#ifdef TRACELS
         traceState("accept+", newPh);
#endif /* TRACELS */

         for (j=0; j < size; j++) {
            // bias[j] = 0.20*bias[j] + 0.40*deviates[j];  // original & questionable
            bias[j] = 0.60*bias[j] + 0.40*deviates[j]; // strict Solis+Wets
         }
      } else {
         // We need to check if the opposite move does any good (move = bias[j] + deviates[j])
         newPh = genPh(vector, -1., deviates, bias, outlev, logFile); // 2x - zeta = x - move
         if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
            num_successes++;
            num_failures = 0;
            vector = newPh;
            modified = TRUE;
#ifdef TRACELS
            traceState("accept-", newPh);
#endif /* TRACELS */
            for (j=0; j < size; j++) {
               // bias[j] -= 0.40*deviates[j]; // incorrect
               bias[j] = 0.60*bias[j] - 0.40*deviates[j]; // correct if deviates is not changed
            }
         } else {
            num_failures++;
            num_successes = 0;
            // vector is unchanged  // x
#ifdef TRACELS
            traceState("reject ", newPh);
#endif /* TRACELS */
            for (j=0; j < size; j++) {
               bias[j] *= 0.50;
            }
         }
      }

      // Check to see if we need to expand or contract
      if (num_successes >= max_successes) {
         temp_rho *= expansion;
         num_successes = num_failures = 0;
      } else if (num_failures >= max_failures) {
         temp_rho *= contraction;
         num_successes = num_failures = 0;
      }
         
      if (temp_rho < lower_bound_on_rho)
         break;  // GMM - this breaks out of the i loop...
   } // i-loop
   return modified;
} // void Solis_Wets::SW(Phenotype &vector)


//  This is pseudo-Solis & Wets in that it adds random deviates to every dimension
//  of the current solution, but the variances vary across dimensions.
//
//  This has a different value of rho for each gene.
// PSW mapped individual
// PSW returns TRUE if it modifies the individual
Boole Pseudo_Solis_Wets::SW(/* not const */ Phenotype &vector, int outlev, FILE *logFile)
{
   register unsigned int i, j, num_successes = 0, num_failures = 0,  all_rho_stepsizes_too_small = 1;
    
   Phenotype newPh;
   Boole modified = FALSE;

#ifdef DEBUG
   (void)fprintf(logFile, "ls.cc/void Pseudo_Solis_Wets::SW(Phenotype &vector)\n");
#endif /* DEBUG */

   //  Initialize the temp_rho's
   for (i=0; i < size; i++) {
      temp_rho[i] = rho[i];
   }
   //  Reset bias or 'momentum'
   for (i=0; i < size; i++) {
      bias[i] = 0.0;
   }

   for (i=0; i < max_its; i++) {
#ifdef DEBUG0
   Real xyz[3];
   for ( int d=0;d<3;d++) xyz[d] = vector.gread(d).real;
   (void)fprintf(logFile, "LS:: %3d #S=%3d #F=%3d %+8.4f p0=%.4f b0=%+7.4f xyz=(%.6f %.4f %.4f)\n", 
                                  i, num_successes, num_failures, vector.evaluate(Normal_Eval), temp_rho[0], bias[0],
                                  xyz[0],xyz[1],xyz[2]);
#endif /* DEBUG0 */
      // Generate deviates
      for (j=0; j < size; j++) {
         deviates[j] = gen_deviates(temp_rho[j]);
      }

      newPh = genPh(vector, +1., deviates, bias, outlev, logFile);
      // Evaluate
      if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
         num_successes++;
         num_failures = 0;
         vector = newPh;
         modified = TRUE;
         for (j=0; j < size; j++) {
            // bias[j] = 0.20*bias[j] + 0.40*deviates[j]; 
            bias[j] = 0.60*bias[j] + 0.40*deviates[j]; // strict Solis+Wets
         }
      } else  {
         //  We need to check if the opposite move does any good (move = bias[j] - deviates[j])
         newPh = genPh(vector, -1., deviates, bias, outlev, logFile);
         if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
            num_successes++;
            num_failures = 0;
            vector = newPh;
            modified = TRUE;
            for (j=0; j < size; j++) {
               // bias[j] -= 0.40*deviates[j];
               bias[j] = 0.60*bias[j] - 0.40*deviates[j]; // correct if deviates is not changed
            }
         } else {
            num_failures++;
            num_successes = 0;
            // vector is unchanged  // x
            for (j=0; j < size; j++) {
               bias[j] *= 0.50;
            }
         }
      }

// DEBUG TRACE used to be here, mp+rh 4/09
      // Check to see if we need to expand or contract
      if (num_successes >= max_successes) {
         for(j=0; j < size; j++) {
            temp_rho[j] *= expansion;
         }
         num_successes = num_failures = 0;
      } else if (num_failures >= max_failures) {
         for(j=0; j < size; j++) {
            temp_rho[j] *= contraction;
         }
         num_successes = num_failures = 0;
      }
      
      //  WEH - Scott's code doesn't do anything!!! no stopping based upon step scale!!!
      //  GMM - corrected Scott's code; this does now stop correctly, based upon step scale.
      //  GMM - This version only exits if all the step sizes are too small...
      all_rho_stepsizes_too_small = 1;
      for(j=0; j < size; j++) {   
         if (temp_rho[j]>= lower_bound_on_rho[j]){
            all_rho_stepsizes_too_small = FALSE;
            break;
         }
      } //  j-loop
      if (all_rho_stepsizes_too_small) {
         break; //  GMM - THIS breaks out of i loop, which IS what we want...
      }
   } //  i-loop
   return modified;
} // Boole Pseudo_Solis_Wets::SW(Phenotype &vector)


int Solis_Wets_Base::search(Individual &solution, Eval *evaluate, int outlev, FILE *logFile)
{

#ifdef DEBUG2
   (void)fprintf(logFile, "ls.cc/int Solis_Wets_Base::search(Individual &solution)\n");
#endif /* DEBUG */

      // Do inverse mapping if SW changed phenotyp 
   if (SW(solution.phenotyp, outlev, logFile)) {
      solution.inverse_mapping();
      }
   ls_count++;      
   return(0);
}

Pattern_Search::Pattern_Search(void)
{
}

Pattern_Search::Pattern_Search(const unsigned int init_size, const unsigned int init_max_success, ConstReal init_step_size, ConstReal init_step_threshold, ConstReal init_expansion, ConstReal init_contraction, ConstReal init_search_frequency)
: size(init_size), max_success(init_max_success), step_size(init_step_size), step_threshold(init_step_threshold), expansion(init_expansion), contraction(init_contraction), localsearch_frequency(init_search_frequency)
{
  current_step_size = step_size;
  pattern = new Real[size];
  index = new unsigned int[size];
  reset_pattern();
  reset_indexes();
  successes = 0;
}

Pattern_Search::~Pattern_Search(void)
{
	delete []pattern;
	delete []index;
}

void Pattern_Search::reset()
{
  current_step_size = step_size;
  reset_pattern();
  reset_indexes();
  successes = 0;
}

void Pattern_Search::reset_pattern() {
  for (unsigned int i=0; i < size; i++) {
    pattern[i] = 0.0;
  }
}

void Pattern_Search::reset_indexes() {
	for (unsigned int i=0; i < size; i++) {
		index[i] = i;
	}
}

void Pattern_Search::shuffle_indexes() {
	int select;
	unsigned int temp;
	for (unsigned int i=size; i > 1; i--) {
		select = rand() % i;
		temp = index[select];
		index[select] = index[i-1];
		index[i-1] = temp;
	}
}

int Pattern_Search::terminate(void) const
{
   return (0);
}

int Pattern_Search::search(Individual &solution, Eval *evaluate, int outlev, FILE * logFile)
{
  // TODO: implement scaling?
  // MP TODO I dont think the Eval parm is used.

  if (ranf() >= localsearch_frequency) {
    return(0);
  }

	reset();
  Phenotype base = solution.phenotyp;
  Phenotype newPoint;
  // evaluate function at base point
  while (current_step_size > step_threshold) {
    // do exploratory moves
    //fprintf(stderr, "base point energy: %f\n", base.evaluate(Normal_Eval));
    newPoint = exploratory_move(base);
    //fprintf(stderr, "newPoint energy: %f\n", newPoint.evaluate(Normal_Eval));
    if (newPoint.evaluate(Normal_Eval) < base.evaluate(Normal_Eval)) {
      // new point is more favorable than base point
      // set new point as base point
      base = newPoint;

      while (true) {
        newPoint = pattern_explore(base);
        if (newPoint.evaluate(Normal_Eval) < base.evaluate(Normal_Eval)) {
					successes++;
          base = newPoint;
        }
        else {
					break;
					successes = 0;
				}

				if (successes > max_success) {
					//fprintf(stderr, "Expanding step size\n");
					successes = 0;
					current_step_size *= expansion;
				}
      }
    }

    else {
      current_step_size *= contraction;
			successes = 0;
      reset_pattern();
      //fprintf(stderr, "Contracted to %f after %ld evaluations.\n", current_step_size, evaluate.evals());
    }
  }
  
  solution.phenotyp = base;
  solution.inverse_mapping();
  return (0);
}

Phenotype Pattern_Search::exploratory_move(const Phenotype& base) /* not const */ {
  Phenotype newBase(base);
	shuffle_indexes(); 
	unsigned int current_index;
	int direction;

  for (unsigned int i=0; i < size; i++) {
    Phenotype trialPoint(newBase);

		current_index = index[i];
		// pick a random direction
		if (rand()%2 == 0) {
			direction = 1;
		}
		else {
			direction = -1;
		}
    // try first coordinate direction
    trialPoint.write(trialPoint.gread(current_index).real+current_step_size*direction, current_index);
    // if successful, keep new point
    if (trialPoint.evaluate(Normal_Eval) < newBase.evaluate(Normal_Eval)) {
      newBase = trialPoint;
      pattern[current_index] += current_step_size*direction;
    }
    // otherwise, try opposite coordinate and test again
    else {
      trialPoint.write(trialPoint.gread(current_index).real-2.0*current_step_size*direction, current_index);
      if (trialPoint.evaluate(Normal_Eval) < newBase.evaluate(Normal_Eval)) {
        newBase = trialPoint;
        pattern[current_index] -= current_step_size*direction;
      }
    }
  }
  return newBase;
}

Phenotype Pattern_Search::pattern_explore(const Phenotype& base) /* not const */ {
  Phenotype newPoint = pattern_move(base);
  reset_pattern();
  Phenotype newBase = exploratory_move(newPoint);
  return newBase;
}

Phenotype Pattern_Search::pattern_move(const Phenotype& base) const {
  Phenotype newPoint(base);
  for (unsigned int i=0; i < size; i++) {
    newPoint.write(newPoint.gread(i).real + pattern[i] , i);
  }
  return newPoint;
}
//void printDState(logFile,msg,vector,i,prevxyz,startxyz, prevQuat,startQuat,num_successes,num_failures,temp_rho,bias,deviates); 
void printDState(FILE *logFile, const char *const msg,Phenotype &newPh, const int i, const Real prevxyz[3],
                 const Real startxyz[3], const Quat& prevQuat, const Quat& startQuat, 
                 const unsigned int num_successes, const unsigned int num_failures,
                 ConstReal temp_rho, const Real *const  bias, const Real *const  deviates)
{
#ifdef DEBUG
   Real dt; // translation step scalar
   Real ct; // cumulative translation from step 0
   static Quat thisQuat; // for DEBUG only
   Real xyz[3];
   for ( int d=0;d<3;d++) xyz[d] = newPh.gread(d).real;
   dt=0;
   for ( int d=0;d<3;d++) dt += (prevxyz[d]- xyz[d])*(prevxyz[d]- xyz[d]);
   dt = sqrt(dt);
   ct=0;
   for ( int d=0;d<3;d++) ct += (startxyz[d] - xyz[d])*(startxyz[d]- xyz[d]);
   ct = sqrt(ct);
   (void)fprintf(logFile, "\nLS::    %3d %s #S=%d #F=%d %+12.4f p=%4.2f b=(%5.2f %5.2f %5.2f) dev=(%5.2f %5.2f %5.2f)", 
                                  i, msg, num_successes, num_failures, newPh.evaluate(Normal_Eval), temp_rho, bias[0], bias[1], bias[2], 
                                  deviates[0], deviates[1], deviates[2]);
   (void)fprintf(logFile, " xyz=(");
    fprintf(logFile, "%5.2f %5.2f %5.2f", xyz[0], xyz[1], xyz[2]);
   (void)fprintf(logFile, ") ");
    fprintf(logFile, "dT=%5.2f ", dt);
    fprintf(logFile, "cT=%5.2f ", ct);
    // assuming x,y,z,w ignoring structs.h order 
    //fprintf(logFile, "%5.2f %5.2f %5.2f %5.2f",
    thisQuat = newPh.readQuat();
    
   (void)fprintf(logFile, " quat=(");
   fprintf(logFile, "%5.2f %5.2f %5.2f %5.2f", newPh.gread(3).real, newPh.gread(4).real,newPh.gread(5).real,newPh.gread(6).real);
   (void)fprintf(logFile, ")");
   //if (i>0){
    fprintf(logFile, " dQ=%6.1lf", quatDifferenceToAngleDeg(prevQuat,thisQuat)); 
    fprintf(logFile, " cQ=%6.1lf", quatDifferenceToAngleDeg(startQuat,thisQuat)); 
    fprintf(logFile, " \n"); 
   //}
#endif /* DEBUG */

}
