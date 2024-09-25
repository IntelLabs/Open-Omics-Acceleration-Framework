/*

 $Id: investigate.h,v 1.25 2012/04/13 06:22:10 mp Exp $

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

#ifndef INVESTIGATE
#define INVESTIGATE

#include "constants.h"
#include "getpdbcrds.h"
#include "mkRandomState.h"
#include "cnv_state_to_coords.h"
#include "getrms.h"
#include "trilinterp.h"
#include "eintcal.h"
#include "changeState.h"
#include "stateLibrary.h"

#define NUMRMSBINS 80 /* int   NumRmsBins = 40; // NumRmsBins = MaxRms / RmsBinSize; */

extern FILE *logFile;
extern char *programname;


void investigate( const int   Nnb, int Nnb_array[3], GroupEnergy *group_energy,
                    const Real charge[MAX_ATOMS],
                    const Real abs_charge[MAX_ATOMS],
                    const Real qsp_abs_charge[MAX_ATOMS],
                    const Boole B_calcIntElec,
                    /* not const */ Real crd[MAX_ATOMS][SPACE], // modified in cnv_state_to_coords
                    const Real crdpdb[MAX_ATOMS][SPACE],

                    const EnergyTables *const ptr_ad_energy_tables,

                    const int   maxTests,
                #include "map_declare.h"
                    const int   natom,
                    const NonbondParam *const nonbondlist,
                    const int   ntor,
                    const int   tlist[MAX_TORS+1][MAX_ATOMS],
                    const int   type[MAX_ATOMS],
                    const Real vt[MAX_TORS][SPACE],
                    const Boole B_isGaussTorCon,
                    const unsigned short US_torProfile[MAX_TORS][NTORDIVS],
                    const Boole B_isTorConstrained[MAX_TORS],
                    const Boole B_ShowTorE,
                    /* not const */ unsigned short US_TorE[MAX_TORS],
                    /* not const */ Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
                    /* not const */ int   N_con[MAX_TORS], // modifed in mkRandomState
                    const Boole B_symmetry_flag,
                    const Boole B_unique_pair_flag,
                    const char  *const FN_rms_ref_crds,
                    const int   OutputEveryNTests,
                    const int   NumLocalTests,
                    ConstReal trnStep,
                    ConstReal torStep,
                    
                    const int   ignore_inter[MAX_ATOMS],
                    
                    const Boole         B_include_1_4_interactions,
                    ConstReal scale_1_4,
                    ConstReal scale_eintermol,


                    ConstReal unbound_internal_FE,
                    /* not const */ GridMapSetInfo *const info, // modified in mkRandomState
                    const Boole B_use_non_bond_cutoff,
                    const Boole B_have_flexible_residues,
                    const Boole B_rms_heavy_atoms_only,
                    const int h_index,
		    const int true_ligand_atoms,
                    const int   outlev,
		    FILE *logFile);
#endif
#if MPIQUE
void investigate(
                int   Nnb,
                Real charge[MAX_ATOMS],
                Real abs_charge[MAX_ATOMS],
                Real qsp_abs_charge[MAX_ATOMS],
                Boole B_calcIntElec,
                Real crd[MAX_ATOMS][SPACE],
                Real crdpdb[MAX_ATOMS][SPACE],

                EnergyTables *ptr_ad_energy_tables,

                int   maxTests,
                #include "map_declare.h"
                int   natom,
                NonbondParam *nonbondlist,
                int   ntor,
                int   outlev,
                int   tlist[MAX_TORS+1][MAX_ATOMS],
                int   type[MAX_ATOMS],
                Real vt[MAX_TORS][SPACE],
                Boole B_isGaussTorCon,
               unsigned short US_torProfile[MAX_TORS][NTORDIVS],
                Boole B_isTorConstrained[MAX_TORS],
                Boole B_ShowTorE,
               unsigned short US_TorE[MAX_TORS],
                Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
                int   N_con[MAX_TORS],
                Boole B_symmetry_flag,
                Boole B_unique_pair_flag,
                char  *FN_rms_ref_crds,
                int   OutputEveryNTests,
                int   NumLocalTests,
                ConstReal trnStep,
                ConstReal torStep,
                
                int   ignore_inter[MAX_ATOMS],
                
                const Boole         B_include_1_4_interactions,
                ConstReal scale_1_4,
                ConstReal scale_eintermol,
                

                ConstReal unbound_internal_FE,
                GridMapSetInfo *info,
                Boole B_use_non_bond_cutoff,
                Boole B_have_flexible_residues,
                Boole B_rms_heavy_atoms_only,
                int h_index);
#endif
