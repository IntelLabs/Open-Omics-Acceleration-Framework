/*

 $Id: hybrids.h,v 1.23 2014/06/12 01:44:07 mp Exp $

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

/*******************************************************************************
	Due to the fact that *.c files are compiled differently than *.cc
	(because of the mangled C++ names), this is a separate header file
	that serves much the same purpose as other prototypes headers.

			2/6/96  rsh
    
        Can get rid of the individual #if-#endif brackets around the fcns.
*******************************************************************************/
#include "constants.h"
#ifndef _STRUCTS_H
#include "structs.h"
#endif
#include "gs.h"
#include "ls.h"
#include "support.h"
#include <stdio.h>

#ifndef CALL_GLSS
#define CALL_GLSS

State call_glss(Global_Search *const global_method, Local_Search *const local_method, 
		const State& now, const unsigned int num_evals, const unsigned int pop_size, 
		const int outlev, FILE *logFile, const Output_pop_stats& extOutput_pop_stats, Molecule *const mol,
		Eval *const evaluate,
		const Boole B_RandomTran0, const Boole B_RandomQuat0, const Boole B_RandomDihe0,
        const GridMapSetInfo *const info, const char *const FN_pop_file,
        /* not const */ int end_of_branch[MAX_TORS]);


Representation **generate_R(int num_torsions, GridMapSetInfo *info );

Representation **generate_R_quaternion(int num_torsions, GridMapSetInfo *info );

Genotype generate_Gtype(const int num_torsions, const GridMapSetInfo *const info, int outlev, FILE *logFile);

Phenotype generate_Ptype(const int num_torsions, const GridMapSetInfo *const info, Eval *evaluate, int outlev, FILE *logFile);

Individual random_ind(const int num_torsions, const GridMapSetInfo *const info, Eval *evaluate, int outlev, FILE *logFile );

#endif


#ifndef CALL_GLSS_TORS
#define CALL_GLSS_TORS

/* currently unused in AutoDock 2014: ***
State call_glss_tors(Global_Search *const global_method, Local_Search *const local_method, 
		const State& now, const unsigned int num_evals, const unsigned int pop_size, 
		const int outlev, FILE *logFile, const Output_pop_stats& extOutput_pop_stats, Molecule *const mol,
		const Boole B_RandomDihe0,
        const GridMapSetInfo *const info, const char *const FN_pop_file);

Representation **generate_R_tors(int num_torsions, GridMapSetInfo *info );
Genotype generate_Gtype_tors(int num_torsions, GridMapSetInfo *info );
Phenotype generate_Ptype_tors(int num_torsions, GridMapSetInfo *info );
Individual random_ind_tors(const int num_torsions, const GridMapSetInfo *const info );
 ***/

#endif

#ifndef CALL_LS
#define CALL_LS

State call_ls(Local_Search *local_method, const State& now, unsigned int pop_size, Molecule *mol, Eval *evaluate, int outlev, FILE *logFile);

#endif


#ifndef CALL_GS
#define CALL_GS

State call_gs(Global_Search *global_method, State& now, unsigned int num_evals, unsigned int pop_size,
              Molecule *mol,
              Output_pop_stats& extOutput_pop_stats,
              GridMapSetInfo *info,
              int end_of_branch[MAX_TORS],
	      int outlev,
	      FILE * logFile);

#endif


#ifndef CALL_PSO
#define CALL_PSO

//State call_cpso(Local_Search *const local_method, const State& sInit, const int n_exec, const int S, const int D, 
//                double *const xmin, double *const xmax, const unsigned int num_evals, const int K,
//		ConstDouble c1, ConstDouble c2, const int outlev);

#endif

#ifndef MMM
#define MMM

void minmeanmax( FILE *const fp, const Population &pop, const int num_its, const GridMapSetInfo *const info );

#endif
