/*

 $Id: trilinterp.h,v 1.17 2012/08/18 00:00:29 mp Exp $

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

#ifndef TRILINTERP
#define TRILINTERP

#define ENERGYPENALTY 500.0	/* Energy factor which is multiplied by distance
				   squared from centre of grid, to penalize atoms
				   outside grid */

#define NULL_EVDW ((Real *)NULL)
#define NULL_ELEC ((Real *)NULL)
#define NULL_EVDW_TOTAL ((Real *)NULL)
#define NULL_ELEC_TOTAL ((Real *)NULL)
#define NULL_IGNORE_INTERMOL ((int *)NULL)
#define NULL_ENERGY_BREAKDOWN ((EnergyComponent *)NULL)  /* also used in eintcal calls */

#include "constants.h"
#include "structs.h"

Real  trilinterp( CONST_INT first_atom,
 CONST_INT last_atom,
 CONST_FLOAT tcoord[MAX_ATOMS][SPACE], // temporary coordinates
 CONST_FLOAT charge[MAX_ATOMS], // partial atomic charges
 CONST_FLOAT abs_charge[MAX_ATOMS], 
 CONST_INT   type[MAX_ATOMS], // atom type of each atom
 #include "map_declare.h"
 const GridMapSetInfo *const info, // info->lo[X],info->lo[Y],info->lo[Z],    minimum coordinates in x,y,z
 const int ignore_inter[MAX_ATOMS], // array of booleans, says to ignore computation intermolecular energies per atom
    EnergyComponent	peratomE[MAX_ATOMS],        // output if not NULL - intermolecular energies
    EnergyComponent	*p_totalE,        // output if not NULL - total energy components
 /* not const */ EnergyComponent *energy_component // set if not NULL - breakdown by elec/vdW_Hb/desolv
 );

#endif

