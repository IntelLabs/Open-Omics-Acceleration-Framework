/*

 $Id: analysis.h,v 1.28 2012/04/13 06:22:10 mp Exp $

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

#include "constants.h"
#include "getpdbcrds.h"
#include "stateLibrary.h"
#include "cnv_state_to_coords.h"
#include "sort_enrg.h"
#include "cluster_analysis.h"
#include "prClusterHist.h"
#include "getrms.h"
#include "eintcal.h"
#include "trilinterp.h"
#include "print_rem.h"
#include "strindex.h"
#include "print_avsfld.h"

#ifndef ANALYSIS
#define ANALYSIS

void  analysis( const int   Nnb, 
	        int Nnb_array[3],
	        GroupEnergy *group_energy,
		const int true_ligand_atoms,
                const char  atomstuff[MAX_ATOMS][MAX_CHARS], 
                const Real charge[MAX_ATOMS], 
                const Real abs_charge[MAX_ATOMS], 
                const Real qsp_abs_charge[MAX_ATOMS], 
                const Boole B_calcIntElec,
                ConstReal clus_rms_tol, 
                const Real crdpdb[MAX_ATOMS][SPACE], 
                
                const EnergyTables *ptr_ad_energy_tables,
                #include "map_declare.h"
                const Real econf[MAX_RUNS], 
                const int   irunmax, 
                const int   natom, 
                const NonbondParam *nonbondlist, 
                const int   nconf, 
                const int   ntor, 
                State hist[MAX_RUNS],  // modified in analysis.cc (hack)
                const char  *smFileName, 
                const Real sml_center[SPACE], 
                const Boole B_symmetry_flag, 
                const Boole B_unique_pair_flag, 
                const int   tlist[MAX_TORS+1][MAX_ATOMS], 
                const int   type[MAX_ATOMS], 
                const Real vt[MAX_TORS][SPACE],
                const char  *rms_ref_crds,
                ConstReal   torsFreeEnergy,
                const Boole B_write_all_clusmem,
                const int ligand_is_inhibitor,
                const int   ignore_inter[MAX_ATOMS],
                const Boole   B_include_1_4_interactions,
                ConstReal   scale_1_4,

                ConstReal   unbound_internal_FE,
                
                const GridMapSetInfo *const info,
                const Boole B_use_non_bond_cutoff,
                const Boole B_have_flexible_residues,
                const Boole B_rms_atoms_ligand_only,
                const Unbound_Model ad4_unbound_model,
                const Boole B_rms_heavy_atoms_only,
                const int h_index,
                const int   outlev,
		FILE *logFile
               );
#endif
