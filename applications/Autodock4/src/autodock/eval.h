/*

 $Id: eval.h,v 1.32 2014/06/12 01:44:07 mp Exp $

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

/********************************************************************
    The header file for the eval class

                                rsh 09/06/95
********************************************************************/
#ifndef _EVAL_H
#define _EVAL_H

#include <stdio.h>
#include "structs.h"
#include "rep.h"
#include "qmultiply.h"
#include "cnv_state_to_coords.h"
#include "trilinterp.h"
#include "eintcal.h"
#include "energy.h"

#ifdef DEBUG
extern FILE *logFile;
#endif

#if defined(USING_COLINY)
void make_state_from_rep(const double *const x, const int n, State *const now);
#endif

void make_state_from_rep(const Representation *const *const rep, State *const stateNow, int outlev, FILE *logFile);

class Eval
{
   private:
	// the following are constant for a particular receptor-ligand docking
	//  and so can be shared read-only by concurrent evaluations:
      unsigned int num_evals;
      int natom, Nnb;
      int *Nnb_array;
      GridMapSetInfo *info;
      MapType *map;
      Boole B_calcIntElec, B_isGaussTorCon, B_ShowTorE;
      unsigned short *US_TorE, (*US_torProfile)[NTORDIVS]; // unused? M Pique
      int *type, (*tlist)[MAX_ATOMS];
      NonbondParam *nonbondlist;
      Real *charge, *abs_charge, *qsp_abs_charge;
      Real (*crd)[SPACE], (*vt)[SPACE];
      EnergyTables *ptr_ad_energy_tables;
      Boole *B_isTorConstrained;
      Molecule mol;
      int *ignore_inter;  // [MAX_ATOMS] gmm 2002-05-21, for CA, CB in flexible sidechains
      Boole         B_include_1_4_interactions; // gmm 2005-01-8, for scaling 1-4 nonbonds
      Real scale_1_4;                  // gmm 2005-01-8, for scaling 1-4 nonbonds
      Real scale_eintermol;  // for scaling intermolecular energy term
      //ParameterEntry *parameterArray;
      Real  unbound_internal_FE;
      Boole B_compute_intermol_energy; // use for computing unbound state
      Boole B_use_non_bond_cutoff;  // set this to FALSE if we are computing unbound extended conformations
      Boole B_have_flexible_residues;
      int true_ligand_atoms;
      int outlev;
	// the following are not constant for a particular ligand docking
	//  and so must be allocated in cases of concurrent evaluations:
      State stateNow;
      Real (*crdpdb)[SPACE];
      EnergyComponent	*peratomE;        // output if not NULL - intermolecular energies
      GroupEnergy  *group_energy;
      FILE *logFile;

   public:
      Eval(void);
      void setup( /* not const */ Real init_crd[MAX_ATOMS][SPACE], // not const since pointers are copied, not contents
                  /* not const */ Real  init_charge[MAX_ATOMS],
                  /* not const */ Real  init_abs_charge[MAX_ATOMS],
                  /* not const */ Real  init_qsp_abs_charge[MAX_ATOMS],
                  /* not const */ int            init_type[MAX_ATOMS], int init_natom,
                  GridMapSetInfo *init_info,
                  MapType  *init_map,

		    EnergyComponent	*peratomE,        // output if not NULL - intermolecular energies

                  /* not const */ NonbondParam *init_nonbondlist,
                  /* not const */ EnergyTables   *init_ptr_ad_energy_tables,
                  const int init_Nnb,
		  int *init_Nnb_array,
		  GroupEnergy *init_group_energy,
                  const Boole          init_B_calcIntElec,
                  const Boole          init_B_isGaussTorCon,
		  /* not const */ Boole init_B_isTorConstrained[MAX_TORS],
                  const Boole          init_B_ShowTorE,
		  /* not const */ unsigned short init_US_TorE[MAX_TORS],
                  /* not const */ unsigned short init_US_torProfile[MAX_TORS][NTORDIVS],
                  /* not const */ Real  init_vt[MAX_TORS][SPACE],
		  /* not const */ int init_tlist[MAX_TORS+1][MAX_ATOMS],
                  /* not const */ Real  init_crdpdb[MAX_ATOMS][SPACE], 
                  const State& stateInit, 
		  const Molecule& molInit,
                  /* not const */ int            init_ignore_inter[MAX_ATOMS], // values are copied
                  const Boole          init_B_include_1_4_interactions, // gmm 2005-01-8, for scaling 1-4 nonbonds
                  ConstReal   init_scale_1_4,                   // gmm 2005-01-8, for scaling 1-4 nonbonds
                  ConstReal   init_scale_eintermol,
                  ConstReal   init_unbound_internal_FE,
                  const Boole  init_B_use_non_bond_cutoff,  // set this to FALSE if we are computing unbound extended conformations
                  const Boole  init_B_have_flexible_residues,
		  int init_true_ligand_atoms,
		  int init_outlev,
		  FILE *init_logFile
                  );
      void update_crdpdb( Real init_crdpdb[MAX_ATOMS][SPACE], 
                        Real init_vt[MAX_TORS][SPACE] );

      double operator()(const Representation *const *const );
      double operator()(const Representation *const *const , int); // GMM - allows calculation of a particular term of the total energy
#if defined(USING_COLINY)
      double operator()(double*, int);
#endif
      double eval();    // WEH - a basic change that facilitates the use of Coliny
      double eval(const int); // GMM - allows calculation of a particular term of the total energy
      unsigned int evals(void);
      void reset(void);
      int write(const Representation *const *const rep, const int true_ligand_atoms, const int outlev, FILE *logFile);
      void compute_intermol_energy(const Boole init_B_compute_intermol_energy); // for computing unbound state
};

inline Eval::Eval(void)
: num_evals(0)
{
}

inline void Eval::setup(/* not const */ Real init_crd[MAX_ATOMS][SPACE], // not const since pointers are copied, not values
                        /* not const */ Real init_charge[MAX_ATOMS],
                        /* not const */ Real init_abs_charge[MAX_ATOMS],
                        /* not const */ Real init_qsp_abs_charge[MAX_ATOMS],
                        /* not const */ int init_type[MAX_ATOMS],
                        const int init_natom,
                        GridMapSetInfo *const init_info,
                        MapType *init_map,

			    EnergyComponent	*init_peratomE,        // output if not NULL - intermolecular energies
                        /* not const */ NonbondParam *const init_nonbondlist,
                        /* not const */ EnergyTables   *const init_ptr_ad_energy_tables,
                        const int init_Nnb,
			int *init_Nnb_array,
			GroupEnergy *init_group_energy,
                        const Boole init_B_calcIntElec, 
                        const Boole init_B_isGaussTorCon,
                        /* not const */ Boole init_B_isTorConstrained[MAX_TORS], // values are not copied but pointers
                        const Boole init_B_ShowTorE,
                        /* not const */ unsigned short init_US_TorE[MAX_TORS],
                        /* not const */ unsigned short init_US_torProfile[MAX_TORS][NTORDIVS],
                        /* not const */ Real init_vt[MAX_TORS][SPACE],
                        /* not const */ int init_tlist[MAX_TORS+1][MAX_ATOMS],
                        /* not const */ Real init_crdpdb[MAX_ATOMS][SPACE],
                        const State& stateInit,
                        const Molecule& molInit,

                        /* not const */ int init_ignore_inter[MAX_ATOMS],

                        const Boole init_B_include_1_4_interactions,
                        ConstReal  init_scale_1_4,
                        ConstReal  init_scale_eintermol, 


                        ConstReal  init_unbound_internal_FE,
                        const Boole init_B_use_non_bond_cutoff,  // set this to FALSE if we are computing unbound extended conformations
                        const Boole init_B_have_flexible_residues,
		       int init_true_ligand_atoms,
		       int init_outlev,
		       FILE *init_logFile
                       )

{
    register int i;

    crd = init_crd;
    charge = init_charge;
    abs_charge = init_abs_charge;
    qsp_abs_charge = init_qsp_abs_charge;
    type = init_type;
    natom = init_natom;
    map = init_map;

    nonbondlist = init_nonbondlist;
    ptr_ad_energy_tables = init_ptr_ad_energy_tables;
    Nnb = init_Nnb;
    Nnb_array= init_Nnb_array;
    group_energy= init_group_energy;
    B_calcIntElec = init_B_calcIntElec;
    B_isGaussTorCon = init_B_isGaussTorCon;
    B_isTorConstrained = init_B_isTorConstrained;
    B_ShowTorE = init_B_ShowTorE;
    US_TorE = init_US_TorE;
    US_torProfile = init_US_torProfile;
    vt = init_vt;
    tlist = init_tlist;
    crdpdb = init_crdpdb;
    // set all of the components of the State, one at a time...
    copyState( &stateNow, stateInit );
#ifdef DEBUG
    pr(logFile, "\n\nstateNow:\n");
    printState( logFile, stateNow, 2 );
#endif
    num_evals = 0;
    ignore_inter = init_ignore_inter;
    peratomE = init_peratomE;
    static EnergyComponent zeroEC;
    for (i=0; i<MAX_ATOMS; i++) {
       peratomE[i] = zeroEC;
    }
    mol = molInit;

    B_include_1_4_interactions = init_B_include_1_4_interactions;
    scale_1_4 = init_scale_1_4;
    scale_eintermol = init_scale_eintermol;

    //parameterArray = init_parameterArray;

    unbound_internal_FE = init_unbound_internal_FE;

    info = init_info;
    map = init_map;
    B_compute_intermol_energy = TRUE; // default is "Yes, calculate the intermolecular energy".

    B_use_non_bond_cutoff = init_B_use_non_bond_cutoff;  // set this to FALSE if we are computing unbound extended conformations
    B_have_flexible_residues = init_B_have_flexible_residues;
    true_ligand_atoms = init_true_ligand_atoms;
    outlev = init_outlev;
    logFile = init_logFile;
}

inline void Eval::update_crdpdb( Real init_crdpdb[MAX_ATOMS][SPACE], 
                               Real init_vt[MAX_TORS][SPACE])
{
    crdpdb = init_crdpdb;
    vt = init_vt;
}

inline void Eval::compute_intermol_energy(const Boole init_B_compute_intermol_energy)
    // For computing the conformation and the internal energy of the unbound state.
{
    B_compute_intermol_energy = init_B_compute_intermol_energy;
}


inline unsigned int Eval::evals(void)
{
   return(num_evals);
}

inline void Eval::reset(void)
{
   num_evals = 0;
}

#endif
