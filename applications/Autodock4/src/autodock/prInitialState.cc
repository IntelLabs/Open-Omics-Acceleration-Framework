/*

 $Id: prInitialState.cc,v 1.18 2014/06/12 01:44:07 mp Exp $

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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "prInitialState.h"
#include "writePDBQT.h"


extern int keepresnum;

void prInitialState(
    const EnergyBreakdown *p_eb,
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
    const int outlev,
    FILE *logFile
    )

{
    char *descriptor = "INITIAL STATE:  ";
    register int i = 0;
    int a = 0;
    Real emap_total = 0.0;
    Real elec_total = 0.0;
    Real emap_flexres_total = 0;
    Real elec_flexres_total = 0;


    pr( logFile, "\n\t\t%s\n\t\t______________\n\n\n", descriptor );

    pr( logFile, "%sUSER    Transformed Initial Coordinates\n", descriptor );
    for (i = 0;  i < natom;  i++) {
        if (keepresnum > 0) print_PDBQT_atom_resstr( logFile, descriptor, i, atomstuff[i],  crd, 
	     1.0, 0.0, charge[i], "", "\n");
        else  print_PDBQT_atom_resnum( logFile, descriptor, i, atomstuff[i],  0, crd, 
	     1.0, 0.0, charge[i], "", "\n");
    } /* i */
    pr( logFile, "%sTER\n\n\n", descriptor );

    pr( logFile, "\t\tINITIAL ENERGY BREAKDOWN\n" );
    pr( logFile, "\t\t________________________\n" );
    pr( logFile, "\n\nEnergy of starting position of Small Molecule by atom: \n\n" );

    print_atomic_energies( natom, atomstuff, type, peratomE, charge, outlev, logFile);

    for (a=0; a<true_ligand_atoms; a++) {
        emap_total += peratomE[a].vdW_Hb+peratomE[a].desolv;
        elec_total += peratomE[a].elec;
    }
    for (a=true_ligand_atoms; a<natom; a++) {
        emap_flexres_total += peratomE[a].vdW_Hb+peratomE[a].desolv;
        elec_flexres_total += peratomE[a].elec;
    }
    
	pr( logFile, "\n\n" );
    printEnergies( p_eb, "Initial ", ligand_is_inhibitor, emap_total, elec_total,
     B_have_flexible_residues, emap_flexres_total, elec_flexres_total, ad4_unbound_model, outlev, logFile );
    pr( logFile, "\n\n" );

    flushLog;
}
/* EOF */
