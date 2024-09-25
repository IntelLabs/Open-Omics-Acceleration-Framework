/*

 $Id: calculateEnergies.h,v 1.15 2012/08/18 00:00:29 mp Exp $

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

#ifndef CALCULATEENERGIES
#define CALCULATEENERGIES
#include <stdio.h>
#include "autocomm.h"
#include "constants.h"
#include "structs.h"

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
    GridMapSetInfo       *info,                     // input  info->lo[X],info->lo[Y],info->lo[Z],    minimum coordinates in x,y,z
    const int                  ignore_inter[MAX_ATOMS],   // input  array of booleans, says to ignore computation intermolecular energies per atom
    EnergyComponent	peratomE[MAX_ATOMS],        // output if not NULL - intermolecular energies
    EnergyComponent	*p_totalE,        // output if not NULL - total energy components

    // eintcal
    const NonbondParam * const         nonbondlist,       // input  list of nonbonds
    const EnergyTables   *ptr_ad_energy_tables,     // input  pointer to AutoDock intermolecular, dielectric, solvation lookup tables
    const int            Nnb,                       // input  total number of nonbonds
    const Boole          B_calcIntElec,             // input  boolean whether we must calculate internal electrostatics
    const Boole          B_include_1_4_interactions,// input  boolean whether to include 1,4 interactions as non-bonds
    ConstReal            scale_1_4,                 // input  scaling factor for 1,4 interactions, if included
    ConstReal            scale_eintermol,                 // input  scaling factor for intermolecular energies
    const Real           qsp_abs_charge[MAX_ATOMS], // input  q-solvation parameters
    const Boole          B_use_non_bond_cutoff,     // input  boolean whether to use a nonbond distance cutoff
    const Unbound_Model ad4_unbound_model,
    const int outlev,
    FILE *logFile

);

void update_energy_breakdown( /* not const */ EnergyBreakdown * eb ); //FIXME: Steffen:  this misses the extra parameter const Unbound_Model ad4_unbound_model

void initialise_energy_breakdown ( /* not const */ EnergyBreakdown *const eb,
                                   ConstReal   torsFreeEnergy, 
                                   ConstReal   unbound_internal_FE );

EnergyBreakdown calculateBindingEnergies(
    int                  natom,                     // input  number of atoms
    int                  ntor,                      // input  number of torsions
    ConstReal                  unbound_internal_FE,       // input  pre-calculated internal energy of unbound state
    ConstReal                  torsFreeEnergy,            // input  constant times number of freely-rotatable bonds
    Boole                B_have_flexible_residues,  // input  boolean whether we have flexible residues in protein

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

);

void update_binding_energy_breakdown( /* not const */ EnergyBreakdown *const eb, const Unbound_Model ad4_unbound_model);

void initialise_binding_energy_breakdown ( EnergyBreakdown *const eb,
                                           ConstReal   torsFreeEnergy, 
                                           ConstReal   unbound_internal_FE,
                                           const Unbound_Model ad4_unbound_model);
#endif
