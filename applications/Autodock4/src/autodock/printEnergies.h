/*

 $Id: printEnergies.h,v 1.15 2014/06/12 01:44:07 mp Exp $

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

#ifndef PRINTENERGIES
#define PRINTENERGIES

#include "autocomm.h"
#include "constants.h"
#include "structs.h"

void printEnergies( const EnergyBreakdown *const eb,
                    const char  *prefixString, 
                    const int ligand_is_inhibitor,
                    ConstReal emap_total,
                    ConstReal elec_total,
                    const Boole B_have_flexible_residues,
                    ConstReal emap_flexres_total,
                    ConstReal elec_flexres_total,
                    const Unbound_Model ad4_unbound_model,
		    int outlev,
		    FILE *logFile);

void printStateEnergies( const EnergyBreakdown *eb,
			 const char *const prefixString, 
			 const int ligand_is_inhibitor,
			 int outlev,
			 FILE *stateFile);
#endif
