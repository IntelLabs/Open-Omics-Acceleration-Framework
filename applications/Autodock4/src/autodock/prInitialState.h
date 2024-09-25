/*

 $Id: prInitialState.h,v 1.12 2014/06/12 01:44:07 mp Exp $

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

#ifndef PRINITIALSTATE
#define PRINITIALSTATE

#include "constants.h"
#include "print_atomic_energies.h"
#include "printEnergies.h"

void  prInitialState( 
    const EnergyBreakdown *eb,
    const int natom,
    const int true_ligand_atoms,
    const Real crd[MAX_ATOMS][SPACE],
    const char atomstuff[MAX_ATOMS][MAX_CHARS],
    const int type[MAX_ATOMS],
    EnergyComponent peratomE[MAX_ATOMS],
    const Real charge[MAX_ATOMS],
    const int ligand_is_inhibitor,
    const Boole B_have_flexible_residues,
    const Unbound_Model ad4_unbound_model,
    int outlev,
    FILE *logFile
    );

#endif
