/*

 $Id: calculateEnergies.cc,v 1.21 2012/08/18 00:00:29 mp Exp $

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

/* calculateEnergies.cc */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "calculateEnergies.h"
#include "trilinterp.h"
#include "eintcal.h"


EnergyBreakdown calculateEnergies(
    const int            natom,                     // input  number of atoms
    const int            ntor,                      // input  number of torsions
    ConstReal            unbound_internal_FE,       // input  pre-calculated internal energy of unbound state
    ConstReal            torsFreeEnergy,            // input  constant times number of freely-rotatable bonds
    const Boole          B_have_flexible_residues,  // input  boolean whether we have flexible residues in protein

    // trilinterp
    const Real           tcoord[MAX_ATOMS][SPACE],  // input  coordinates of atoms to be trilinearly-interpolated
    CONST_FLOAT          charge[MAX_ATOMS],         // input  partial atomic charges
    CONST_FLOAT          abs_charge[MAX_ATOMS],     // input  absolute magnitude of partial charges
    CONST_INT            type[MAX_ATOMS],           // input  atom type of each atom
    #include "map_declare.h"
    const GridMapSetInfo *const info,               // input  info->lo[X],info->lo[Y],info->lo[Z],    minimum coordinates in x,y,z
    const int            ignore_inter[MAX_ATOMS],   // input  array of booleans, says to ignore computation intermolecular energies per atom
    EnergyComponent	peratomE[MAX_ATOMS],        // output if not NULL - intermolecular energies
    EnergyComponent	*p_totalE,        // output if not NULL - total energy components

    // eintcal
    const NonbondParam *const nonbondlist,          // input  list of nonbonds
    const EnergyTables *const ptr_ad_energy_tables, // input  pointer to AutoDock intermolecular, dielectric, solvation lookup tables
    const int            Nnb,                       // input  total number of nonbonds
    int Nnb_array[3],
    GroupEnergy *group_energy,
    const int true_ligand_atoms,
    const Boole          B_calcIntElec,             // input  boolean whether we must calculate internal electrostatics
    const Boole          B_include_1_4_interactions,// input  boolean whether to include 1,4 interactions as non-bonds
    const Real           scale_1_4,                 // input  scaling factor for 1,4 interactions, if included
    const Real           scale_eintermol,                 // input  scaling factor for intermolecular energies
    const Real           qsp_abs_charge[MAX_ATOMS], // input  q-solvation parameters
    const Boole          B_use_non_bond_cutoff,     // input  boolean whether to use a nonbond distance cutoff
    const int outlev,
    FILE *logFile

)

{
    EnergyBreakdown eb;

    initialise_energy_breakdown( &eb, torsFreeEnergy, unbound_internal_FE );

    // computing trilinear-interpolated energies from atom = 0 to atom < true_ligand_atoms
    // gives the intermolecular energy between the ligand and the protein
    eb.e_inter_moving_fixed = scale_eintermol * trilinterp( 0, true_ligand_atoms, tcoord, charge, abs_charge, type, map, 
	 info, ignore_inter, peratomE, p_totalE,
	 &group_energy->inter_moving_fixed);

    if (B_have_flexible_residues) {
        // computing trilinear-interpolated energies from atom = true_ligand_atoms to atom < natom
        // gives the intramolecular energy within the protein
        // we can ignore the elec_total and emap_total breakdown here
        eb.e_intra_moving_fixed_rec = trilinterp( true_ligand_atoms, natom, tcoord, charge, abs_charge, type, map, 
	     info, ignore_inter, peratomE, NULL, &group_energy->intra_moving_fixed_rec);
    }

    if (ntor > 0) {
        // computing all the nonbond interaction energies fills group_energy structure
        // with intramolecular energy of ligand, intermolecular energy, and intramolecular energy of receptor
        (void) eintcal( nonbondlist, ptr_ad_energy_tables, tcoord, Nnb, 
	Nnb_array, group_energy,
	B_calcIntElec, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, 
	B_use_non_bond_cutoff, B_have_flexible_residues, outlev, logFile) ;
        
        eb.e_intra_moving_moving_lig = group_energy->intra_moving_moving_lig.total;
        eb.e_inter_moving_moving = group_energy->inter_moving_moving.total;
        eb.e_intra_moving_moving_rec = group_energy->intra_moving_moving_rec.total;
    }

    // update the totals in the energy breakdown structure
    update_energy_breakdown( &eb );

    return eb;
} // calculateEnergies()

void update_energy_breakdown( /* not const */ EnergyBreakdown *const eb )
{
    // total intermolecular energy = (1) + (4)
    eb->e_inter     = eb->e_inter_moving_fixed + eb->e_inter_moving_moving;

    // total intramolecular energy = (3)  +  (2) + (5)
    eb->e_intra     = eb->e_intra_moving_moving_lig + eb->e_intra_moving_fixed_rec + eb->e_intra_moving_moving_rec;

    // ligand intramolecular energy = (3)
    eb->e_intra_lig = eb->e_intra_moving_moving_lig;

    // receptor intramolecular energy = (2) + (5)
    eb->e_intra_rec = eb->e_intra_moving_fixed_rec + eb->e_intra_moving_moving_rec;

    // estimated free energy upon binding
    eb->deltaG = eb->e_inter + eb->e_intra + eb->e_torsFreeEnergy - eb->e_unbound_internal_FE;
}

void initialise_energy_breakdown ( /* not const */ EnergyBreakdown *const eb,
                                   ConstReal   torsFreeEnergy, 
                                   ConstReal   unbound_internal_FE )
{
    eb->e_inter_moving_fixed = 0.0;      // (1)  // trilinterp( 0, true_ligand_atoms, ...)
    eb->e_intra_moving_fixed_rec = 0.0;  // (2)  // trilinterp( true_ligand_atoms, natom, ...)
    eb->e_intra_moving_moving_lig = 0.0; // (3)  // eintcal( 0, nb_array[0], ...)            // group_energy[INTRA_LIGAND]
    eb->e_inter_moving_moving = 0.0;     // (4)  // eintcal( nb_array[0], nb_array[1], ...)  // group_energy[INTER]
    eb->e_intra_moving_moving_rec = 0.0; // (5)  // eintcal( nb_array[1], nb_array[2], ...)  // group_energy[INTRA_RECEPTOR]

    eb->e_inter = 0.0;                   // total    intermolecular energy = (1) + (4)
    eb->e_intra = 0.0;                   // total    intramolecular energy = (3)  +  (2) + (5)
    eb->e_intra_lig = 0.0;               // ligand   intramolecular energy = (3)
    eb->e_intra_rec = 0.0;               // receptor intramolecular energy = (2) + (5)

    eb->e_torsFreeEnergy = torsFreeEnergy;            // empirical torsional free energy penalty
    eb->e_unbound_internal_FE = unbound_internal_FE;  // computed internal free energy of the unbound state
    eb->deltaG = 0.0;                    // estimated change in free energy upon binding
}


// calculateBindingEnergies() should only be used at the end of a docking, not during a docking:
// this is because the internal energy will be ignored in the Unbound_same_as_bound case,
// which might create conformations with high internal energies.

EnergyBreakdown calculateBindingEnergies(
    const int            natom,                     // input  number of atoms
    const int            ntor,                      // input  number of torsions
    ConstReal            unbound_internal_FE,       // input  pre-calculated internal energy of unbound state
    ConstReal            torsFreeEnergy,            // input  constant times number of freely-rotatable bonds
    const Boole          B_have_flexible_residues,  // input  boolean whether we have flexible residues in protein

    // trilinterp
    const Real           tcoord[MAX_ATOMS][SPACE],  // input  coordinates of atoms to be trilinearly-interpolated
    CONST_FLOAT          charge[MAX_ATOMS],         // input  partial atomic charges
    CONST_FLOAT          abs_charge[MAX_ATOMS],     // input  absolute magnitude of partial charges
    CONST_INT            type[MAX_ATOMS],           // input  atom type of each atom
    #include "map_declare.h"
    const GridMapSetInfo *const info,               // input  info->lo[X],info->lo[Y],info->lo[Z],    minimum coordinates in x,y,z
    const int            ignore_inter[MAX_ATOMS],   // input  array of booleans, says to ignore computation intermolecular energies per atom
    EnergyComponent	peratomE[MAX_ATOMS],        // output if not NULL - intermolecular energies
    EnergyComponent	*p_totalE,        // output if not NULL - total energy components

    // eintcal
    const NonbondParam *const nonbondlist,          // input  list of nonbonds
    const EnergyTables *const ptr_ad_energy_tables, // input  pointer to AutoDock intermolecular, dielectric, solvation lookup tables
    const int            Nnb,                       // input  total number of nonbonds
    int Nnb_array[3],
    GroupEnergy *group_energy,
    const int true_ligand_atoms,
    const Boole          B_calcIntElec,             // input  boolean whether we must calculate internal electrostatics
    const Boole          B_include_1_4_interactions,// input  boolean whether to include 1,4 interactions as non-bonds
    ConstReal            scale_1_4,                 // input  scaling factor for 1,4 interactions, if included
    const Real           qsp_abs_charge[MAX_ATOMS], // input  q-solvation parameters
    const Boole          B_use_non_bond_cutoff,     // input  boolean whether to use a nonbond distance cutoff
    const Unbound_Model  ad4_unbound_model,
    const int outlev,
    FILE *logFile

)

{
    EnergyBreakdown eb;

    initialise_binding_energy_breakdown( &eb, torsFreeEnergy, unbound_internal_FE , ad4_unbound_model);

    // computing trilinear-interpolated energies from atom = 0 to atom < true_ligand_atoms
    // gives the intermolecular energy between the ligand and the protein
    eb.e_inter_moving_fixed = trilinterp( 0, true_ligand_atoms, tcoord, charge, abs_charge, type, map, 
	 info, ignore_inter, peratomE, p_totalE,
	 &group_energy->inter_moving_fixed);

    if (B_have_flexible_residues) {
        // computing trilinear-interpolated energies from atom = true_ligand_atoms to atom < true_ligand_atoms
        // gives the intramolecular energy within the protein
        // we can ignore the elec_total and emap_total breakdown here
        eb.e_intra_moving_fixed_rec = trilinterp( true_ligand_atoms, natom, tcoord, charge, abs_charge, type, map, 
	     info, ignore_inter, peratomE, NULL, &group_energy->intra_moving_fixed_rec);
    }

    if (ntor > 0) {
        // computing all the nonbond interaction energies fills group_energy structure
        // with intramolecular energy of ligand, intermolecular energy, and intramolecular energy of receptor
        (void) eintcal( nonbondlist, ptr_ad_energy_tables, tcoord, Nnb, 
	Nnb_array, group_energy,
	B_calcIntElec, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, B_use_non_bond_cutoff, B_have_flexible_residues, outlev, logFile) ;
        
        eb.e_intra_moving_moving_lig = group_energy->intra_moving_moving_lig.total;
        eb.e_inter_moving_moving = group_energy->inter_moving_moving.total;
        eb.e_intra_moving_moving_rec = group_energy->intra_moving_moving_rec.total;
    }

    // update the totals in the energy breakdown structure
    update_binding_energy_breakdown( &eb, ad4_unbound_model );

    return eb;
} // calculateBindingEnergies()

void update_binding_energy_breakdown( /* not const */ EnergyBreakdown *const eb, const Unbound_Model ad4_unbound_model )
{
    // total intermolecular energy = (1) + (4)
    eb->e_inter     = eb->e_inter_moving_fixed + eb->e_inter_moving_moving;

    // total intramolecular energy = (3) + (2) + (5)
    eb->e_intra     = eb->e_intra_moving_moving_lig + eb->e_intra_moving_fixed_rec + eb->e_intra_moving_moving_rec;

    // ligand intramolecular energy = (3)
    eb->e_intra_lig = eb->e_intra_moving_moving_lig;

    // receptor intramolecular energy = (2) + (5)
    eb->e_intra_rec = eb->e_intra_moving_fixed_rec + eb->e_intra_moving_moving_rec;

    // Set the internal energy of the unbound state
    switch (ad4_unbound_model) {
        // in AutoDock 4.2, the default unbound model is "unbound is same as bound"
        case Unbound_Default:
        case Unbound_Same_As_Bound:
        default:
            // Update the unbound internal energy, setting it to the current internal energy
            eb->e_unbound_internal_FE = eb->e_intra;  // current internal energy of the ligand + flexres unbound state (2)+(3)+(5)
            eb->deltaG = eb->e_inter +  eb->e_torsFreeEnergy;
            break;
        case User:
        case Extended:
        case Compact:
            // The unbound internal energy has already been set in 
            // initialise_energy_breakdown() to the value passed in unbound_internal_FE
            // There is no need to update here.
            eb->deltaG = eb->e_inter + eb->e_intra + eb->e_torsFreeEnergy - eb->e_unbound_internal_FE;
            break;
    }

    // estimated free energy upon binding
    eb->deltaG = eb->e_inter + eb->e_intra + eb->e_torsFreeEnergy - eb->e_unbound_internal_FE;
}

void initialise_binding_energy_breakdown( /* not const */ EnergyBreakdown *const eb,
                                          ConstReal   torsFreeEnergy, 
                                          ConstReal   unbound_internal_FE ,
                                          const Unbound_Model ad4_unbound_model)
        {
    eb->e_inter_moving_fixed = 0.0;      // (1)  // trilinterp( 0, true_ligand_atoms, ...)
    eb->e_intra_moving_fixed_rec = 0.0;  // (2)  // trilinterp( true_ligand_atoms, natom, ...)
    eb->e_intra_moving_moving_lig = 0.0; // (3)  // eintcal( 0, Nnb_array[0], ...)             // group_energy[INTRA_LIGAND]
    eb->e_inter_moving_moving = 0.0;     // (4)  // eintcal( Nnb_array[0], Nnb_array[1], ...)  // group_energy[INTER]
    eb->e_intra_moving_moving_rec = 0.0; // (5)  // eintcal( Nnb_array[1], Nnb_array[2], ...)  // group_energy[INTRA_RECEPTOR]

    eb->e_inter = 0.0;                   // total    intermolecular energy = (1) + (4)
    eb->e_intra = 0.0;                   // total    intramolecular energy = (3) + (2) + (5)
    eb->e_intra_lig = 0.0;               // ligand   intramolecular energy = (3)
    eb->e_intra_rec = 0.0;               // receptor intramolecular energy = (2) + (5)

    eb->e_torsFreeEnergy = torsFreeEnergy; // empirical torsional free energy penalty

    // Set the internal energy of the unbound state
    switch (ad4_unbound_model) {
        // in AutoDock 4.2, the default unbound model is "unbound is same as bound"
        case Unbound_Same_As_Bound:
        default:
            // Update the unbound internal energy, setting it to the current internal energy
            eb->e_unbound_internal_FE = eb->e_intra_lig;  // current internal energy of the ligand unbound state
            break;
        case User:
        case Extended:
        case Compact:
            // Set the unbound internal energy to the value passed in unbound_internal_FE
            eb->e_unbound_internal_FE = unbound_internal_FE;  // supplied internal energy of the ligand unbound state
            break;
    }

    eb->deltaG = 0.0;                    // estimated change in free energy upon binding
}

// EOF
