/*

 $Id: call_ls.cc,v 1.11 2014/06/12 01:44:07 mp Exp $

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
     Call_ls:  Invokes a local searcher on a docking to try and 
               find the locally optimal solution.  So, the docking
               must be specified BEFORE calling this routine.

				rsh 2/5/96
********************************************************************/
// possibly unnecessary // #include <iostream.h>
#include "ls.h"
#include "support.h"
#include "eval.h"

   #include "constants.h"
   #include "structs.h"
   #include "qmultiply.h"

Representation **cnv_state_to_rep(const State &state)
{
   register int i;
   Representation **retval;

   retval = new Representation *[5];
   retval[0] = new RealVector(1);
   retval[0]->write(state.T.x, 0);
   retval[1] = new RealVector(1);
   retval[1]->write(state.T.y, 0);
   retval[2] = new RealVector(1);
   retval[2]->write(state.T.z, 0);
   retval[3] = new RealVector(4);
   retval[3]->write(state.Q.x, 0);
   retval[3]->write(state.Q.y, 1);
   retval[3]->write(state.Q.z, 2);
   retval[3]->write(state.Q.w, 3);
   retval[4] = new RealVector(state.ntor);
   for(i=0; i<state.ntor; i++) {
      retval[4]->write(state.tor[i], i);
   }

   return(retval);
}

Individual cnv_state_to_ind(const State &original, Eval *evaluate)
{
   // BEGIN DELETION
   // return(Individual(Genotype(5, cnv_state_to_rep(original)), Phenotype(5, cnv_state_to_rep(original))));
   // END DELETION

   // BEGIN ADDITION
   // Added by gmm, 27-MAR-97, to solve these compiler warnings:
   //
   // call_ls.cc:59: warning: In this statement, the initialization of a non-const reference requires a temporary for "Genotype(5,cnv_state_to_rep(original))". (reftemporary)
   // call_ls.cc:59: warning: In this statement, the initialization of a non-const reference requires a temporary for "Phenotype(5,cnv_state_to_rep(original))". (reftemporary)
   // call_ls.cc:59: warning: In this statement, the initialization of a non-const reference requires a temporary for "(Individual(Genotype(5,cnv_state_to_rep(original)),Phenotype(5,cnv_state_to_rep(original))))". (reftemporary)

   Genotype temp_Gtype;
   Phenotype temp_Ptype;

   temp_Gtype = Genotype(5, cnv_state_to_rep(original));
   temp_Ptype = Phenotype(5, cnv_state_to_rep(original), evaluate);

   Individual temp(temp_Gtype, temp_Ptype);

   return(temp);
   // END ADDITION

}

State
call_ls(Local_Search *const local_method, const State& now, 
 const unsigned int pop_size, Molecule *const mol, 
  Eval *evaluate, const int outlev, FILE *logFile) 
{
   unsigned int i;

   evaluate->reset();
   local_method->reset();
   //MP Individual prototype;

   //MP prototype = cnv_state_to_ind(now, evaluate);
   //MPPopulation thisPop(pop_size, evaluate, &prototype); // MPique TODO this evaluate is redundant
   Population thisPop(pop_size, evaluate);

   // MPique 2013:
   // for thread-safety, each population must have its own local search workspaces
   // MP TODO how know if SW or PSW? Ask Ruth


/*
   Real (*ibias)[pop_size], (*irho)[pop_size], (*idev)[pop_size]; 
   local_method (*iSearch)[pop_size];
   ibias = new *Real[pop_size];
   irho= new *Real[pop_size];
   idev= new *Real[pop_size];
   iSearch= new *local_method[pop_size];
*/



      for(i=0;i<pop_size;i++)
      {
	thisPop[i] = cnv_state_to_ind(now, evaluate);  // MP 2014 need this ??
	thisPop[i].mol = mol;
      }


#ifdef DEBUG0
// MP Sept 2011
double e_start;
#endif

for(i=0; i<pop_size; i++)
{
#ifdef DEBUG0
// MP Sept 2011
if(i==0)fprintf(logFile, "LS:: ->call_ls(thisPop[%3d]) %10.6f\n", i,
e_start=thisPop[i].value(Normal_Eval));
#endif
      local_method->search( thisPop[i], evaluate, outlev, logFile );
#ifdef DEBUG0
fprintf(logFile, "LS:: <-call_ls(thisPop[%3d]) %10.6f\n", i,
thisPop[i].value(Normal_Eval));
#endif
   }

   if (pop_size > 1)
	 	thisPop.msort(1);
    (void)fprintf(logFile,"Final-Value: %.3f\n", thisPop[0].value(Normal_Eval));
#ifdef DEBUG0
double e_final = thisPop[0].value(Normal_Eval);
fprintf(logFile, "LS:: orig %10.6f best %10.6f (improvement %10.6f) \n", 
 e_start, e_final, e_start-e_final);
#endif
   return(thisPop[0].state(now.ntor));
}
