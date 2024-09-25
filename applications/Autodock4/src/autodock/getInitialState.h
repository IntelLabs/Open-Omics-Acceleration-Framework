/*

 $Id: getInitialState.h,v 1.23 2012/08/18 00:00:29 mp Exp $

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

#ifndef GETINITIALSTATE
#define GETINITIALSTATE

#include "constants.h"
#include "qmultiply.h"
#include "stateLibrary.h"
#include "initautodock.h"
#include "trilinterp.h"
#include "eintcal.h"
#include "cnv_state_to_coords.h"
#include "prInitialState.h"
#include "timesys.h"

void getInitialState(  
            /* not const */ Real *const Addr_e0,
            ConstReal  e0max,

	    /* not const */ State *const sInit,
	    /* not const */ State *const sMin,
	    /* not const */ State *const sLast,

            const Boole B_RandomTran0,
            const Boole B_RandomQuat0,
            const Boole B_RandomDihe0,

            const Real charge[MAX_ATOMS],
            const Real abs_charge[MAX_ATOMS],
            const Real qsp_abs_charge[MAX_ATOMS],
            /* not const */ Real crd[MAX_ATOMS][SPACE],
            const Real crdpdb[MAX_ATOMS][SPACE],
            const char  atomstuff[MAX_ATOMS][MAX_CHARS],
            EnergyComponent	*peratomE,        // output if not NULL - intermolecular energies

            const EnergyTables *const ptr_ad_energy_tables,

            const Boole B_calcIntElec,
                #include "map_declare.h"
            const int   natom,
            const int   Nnb,
	    int   Nnb_array[3],
	    GroupEnergy   *group_energy,
	    const int   true_ligand_atoms,
            const NonbondParam *const nonbondlist,
            const int   ntor,
            const int   tlist[MAX_TORS+1][MAX_ATOMS],
            const int   type[MAX_ATOMS],
            const Real vt[MAX_TORS][SPACE],
            const int   irun1,
	    const int   MaxRetries,

	    ConstReal  torsFreeEnergy,

            const int   ligand_is_inhibitor,

            const int ignore_inter[MAX_ATOMS],

            const Boole         B_include_1_4_interactions,
            ConstReal  scale_1_4,
	    ConstReal  scale_eintermol,


            ConstReal  unbound_internal_FE,

            const GridMapSetInfo *const info,
            const Boole B_use_non_bond_cutoff,
            const Boole B_have_flexible_residues,
            const Unbound_Model ad4_unbound_model,
            const int   outlev,
	    FILE *logFile
           );

#endif
