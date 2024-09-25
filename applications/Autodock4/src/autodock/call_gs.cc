/*

 $Id: call_gs.cc,v 1.16 2014/06/12 01:44:07 mp Exp $

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
     Call_gs:  Invokes a Global Searcher object on a randomly
               generated population of solution to the docking 
               problem.

				rsh 3/12/96
********************************************************************/
// possibly unnecessary // #include <iostream.h>
#include "gs.h"
#include "support.h"
#include "eval.h"
#include "hybrids.h"

   #include "constants.h"
   #include "structs.h"


State call_gs(Global_Search *global_method, State& now, unsigned int num_evals, unsigned int pop_size,
              Molecule *mol,
	      Eval *evaluate,
              Output_pop_stats& extOutput_pop_stats,
              GridMapSetInfo *info,
              int end_of_branch[MAX_TORS], int outlev, FILE * logFile)
{
   register unsigned int i;

   global_method->reset(extOutput_pop_stats);
   evaluate->reset();

   Population thisPop(pop_size, evaluate);
   thisPop.set_eob(end_of_branch);

   for (i=0; i<pop_size; i++) {
      thisPop[i] = random_ind(now.ntor, info); //random_ind does mapping because global search no longer does 2009/04
      thisPop[i].mol = mol;
   }

    while ((evaluate->evals() < num_evals) && (!global_method->terminate()))
      global_method->search(thisPop, outlev, logFile);

   if (pop_size>1) thisPop.msort(1);

   return( thisPop[0].state(now.ntor) );
}
