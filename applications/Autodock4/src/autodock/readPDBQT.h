/*

 $Id: readPDBQT.h,v 1.21 2014/06/12 01:44:08 mp Exp $

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

#ifndef READPDBQT
#define READPDBQT

#include "structs.h"
#include "constants.h"
#include "stop.h"
#include "get_atom_type.h"
#include "print_2x.h"
#include "mkTorTree.h"
#include "nonbonds.h"
#include "weedbonds.h"
#include "torNorVec.h"
#include "success.h"
#include "openfile.h"
#include "stack.h"

void  readPDBQTLine( char line[LINE_LEN],
                     int  *const ptr_pdbatomnumber,
                     Real crd[SPACE], 
                     Real *const P_q,
                     ParameterEntry *thisparm,
		     const int outlev,
		     FILE *logFile);

Molecule readPDBQT( char  line[LINE_LEN],

              int   num_atm_maps,

              int   *const P_natom,
              Real crdpdb[MAX_ATOMS][NTRN],
              Real charge[MAX_ATOMS],
              Boole *const P_B_haveCharges,
              int   type[MAX_ATOMS],
              int   bondtype[MAX_ATOMS],
              char  pdbaname[MAX_ATOMS][5],

              char  *const pdbqFileName,
              char  *const FN_flexres,
              Boole B_have_flexible_residues,

              char  atomstuff[MAX_ATOMS][MAX_CHARS],
              int   *const P_n_heavy_atoms_in_ligand,
	      int   *const P_true_ligand_atoms,

              Boole *const P_B_constrain,
              int   *const P_atomC1,
              int   *const P_atomC2,
              Real  *const P_sqlower,
              Real  *const P_squpper,

              int   *const P_ntor,
              int   *const P_ntor_ligand,
              int   tortree[MAX_TORS][MAX_ATOMS],
              Real vt[MAX_TORS][NTRN],

              int   *const P_Nnb,
	      int  Nnb_array[3],
              NonbondParam *const nonbondlist,

              const Clock& jobStart,
              const struct tms& tms_jobStart,
              const char  hostnm[MAX_CHARS],

              int   *const P_ntorsdof,

              int   ignore_inter[MAX_ATOMS],
              
              const int   B_include_1_4_interactions,
              
              /* not const */ char  PDBQT_record[MAX_RECORDS][LINE_LEN],

              /* not const */ int end_of_branch[MAX_TORS],
              const int   debug,
              const int   outlev,
	      FILE *logFile
              );
#endif
