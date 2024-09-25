/*

 $Id: structs.h,v 1.36 2012/10/15 17:48:28 mp Exp $

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

#ifndef _STRUCTS_H
#define _STRUCTS_H

#include "constants.h"
#include "typedefs.h"

/* *****************************************************************************
 *      Name: structs.h                                                       *
 *  Function: Defines structures used in Molecular Applications.              *
 *Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
 *----------------------------------------------------------------------------*
 *    Author: Garrett Matthew Morris, The Scripps Research Institute          *
 *      Date: SEP/07/1995                                                     *
 *----------------------------------------------------------------------------*
 *    Inputs: none                                                            *
 *   Returns: nothing                                                         *
 *   Globals: none                                                            *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 02/28/95 GMM     This header added                                         *
 ***************************************************************************** */

/* ____________________________________________________________________________ */

typedef double Vector [3];  // for vectors and points

/* ____________________________________________________________________________ */

typedef struct coord
{
  double x;			/* Cartesian x-coordinate */
  double y;			/* Cartesian y-coordinate */
  double z;			/* Cartesian z-coordinate */
} Coord;

/* ____________________________________________________________________________ */

typedef struct quat
{
  double x;			/* quaternion's x-component */
  double y;			/* quaternion's y-component */
  double z;			/* quaternion's z-component */
  double w;			/* quaternion's w-component */
} Quat;

typedef struct axisangle  /* used in older AutoDock code as "Rot" */
{
  double nx;			/* unit vector's x-component */
  double ny;			/* unit vector's y-component */
  double nz;			/* unit vector's z-component */
  double ang;			/* angle of rotation about unit-vector, radians */
} AxisAngle;

#ifdef AXANG
typedef struct axis
{
    double v[3];
} Axis;

// MP in progress - which is easier to use
typedef struct axisang
{
    union axis {
    Axis axis;
    struct xyz {
	  double nx;			/* unit vector's x-component */
	  double ny;			/* unit vector's y-component */
	  double nz;			/* unit vector's z-component */
	  } xyz;
	};
    double angle;
} AxisAngle;
#endif


/* ____________________________________________________________________________ */

typedef struct energy
{
  double total;			/* total energy */
  double intra;			/* intramolecular energy, a.k.a. "internal" energy */
  double inter;			/* intermolecular energy */
  double FE;			/* estimated Free Energy of binding */
} Energy;

/* ____________________________________________________________________________ */

typedef struct state
{
  Coord T;			/* coordinates of center of molecule */
  Quat Q;			/* rigid-body orientation */
  double tor[MAX_TORS];		/* torsion angles in radians */
  int ntor;			/* number of torsions in molecule */
  int hasEnergy;		/* if 0, this state has an undefined energy */
  Energy e;			/* energy structure */
  Coord Center;			/* original input ligand center (pivot point)
     					as set by 'about' DPF keyword */
} State;

/* ____________________________________________________________________________ */

typedef struct molecule
{
  Real crdpdb[MAX_ATOMS][SPACE];	    /* original coordinates of atoms */
  Real crd[MAX_ATOMS][SPACE];      	/* current coordinates of atoms */
  char atomstr[MAX_ATOMS][MAX_CHARS];	/* strings describing atoms, from PDB file, cols,1-30. */
  int natom;			                /* number of atoms in molecule */
  Real vt[MAX_TORS][SPACE];        	/* vectors  of torsions */
  int tlist[MAX_TORS+1][MAX_ATOMS];	    /* torsion list of movable atoms+root atoms at end */
  State S;		                    	/* state of molecule */
} Molecule;

/* ____________________________________________________________________________ */

typedef struct rotamer
{
  double tor[MAX_TORS_IN_ROTAMER];	/* torsion angles in radians */
  int ntor;			/* number of torsions */
} Rotamer;

/* ____________________________________________________________________________ */

typedef struct chargestruct
{
    double charge;
    double abs_charge;
    double qsp_abs_charge;
} Charge;

/* ____________________________________________________________________________ */

typedef struct atom
{
  double    coords[3];			    /* transformed coordinates */
  double    crdpdb[3];			    /* input PDB coordintates */
  Boole     has_charge;			    /* TRUE if the atom has a charge */

  double    charge;
  double    abs_charge;
  double    qsp_abs_charge;

  int       type;			        /* atom type as integer */
  char      type_string[MAX_CHARS]; /* atom type as string */
  Boole     is_hydrogen;		    /* TRUE if atom is a hydrogen */

  int       serial;			        /* serial ID */
  char      name[5];			    /* PDB atom name; formerly "pdbaname" */
  char      stuff[MAX_CHARS];       /* PDB atom string; formerly "atomstuff" */

  int       nnb;			        /* number of non-bonds for this atom */
} Atom;
/* ____________________________________________________________________________ */

typedef struct bond
{
  Atom *atom1;
  Atom *atom2;
  double bondLength;
  Coord bondVector;
} Bond;
/* ____________________________________________________________________________ */

typedef struct pair_id
{
  Atom *atom1;			/* pointer to one atom in pair */
  Atom *atom2;			/* pointer to other atom */
} PairID;
/* ____________________________________________________________________________ */

typedef struct dist_constraint
{
  PairID bond;			/* two atoms defining distance constraint */
  double length;		/* current bond length */
  double lower;			/* lower bound on distance */
  double upper;			/* upper bound on distance */
} DisCon;

/* ____________________________________________________________________________ */

typedef struct torsion
{
  PairID rotbnd;		/* atom serial-IDs of rotatable bond */
  int nmoved;			/* number of atoms moved by this */
  int IDmove[MAX_ATOMS];	/* atom serial-IDs of atoms moved by this */
  Coord vt;			/* bond-vector of rotatable bond */
} Torsion;

/* ______________________________________________________________________________
** Molecular fragments, side-chains, even entire ligands */

typedef struct group
{
  int natom;			/* Number of atoms in fragment */
  Atom atm[MAX_ATOMS];		/* Atom data */
  int ntor;			/* Number of torsions in fragment */
  Torsion tors[MAX_TORS];	/* Torsion data */
  int nnb;			/* Number of non-bonds in fragment */
  PairID nbs[MAX_NONBONDS];	/* Non-bond data */
  Boole B_constrain;		/* TRUE if any distance constraints */
  DisCon distcon;		/* Distance constraint data */
  char pdbqfilnam[PATH_MAX];	/* PDBQ filename holding these data */
} Group;

#include "grid.h"

/* ______________________________________________________________________________
** Parameter Dictionary */

#include "parameters.h"
/* ______________________________________________________________________________ */

typedef struct linear_FE_model
{
    double coeff_vdW;                 // Free energy coefficient for van der Waals term
    double coeff_hbond;               // Free energy coefficient for H-bonding term
    double coeff_estat;               // Free energy coefficient for electrostatics term
    double coeff_desolv;              // Free energy coefficient for desolvation term
    double coeff_tors;                // Free energy coefficient for torsional term

    double stderr_vdW;                // Free energy standard error for van der Waals term
    double stderr_hbond;              // Free energy standard error for H-bonding term
    double stderr_estat;              // Free energy standard error for electrostatics term
    double stderr_desolv;             // Free energy standard error for desolvation term
    double stderr_tors;               // Free energy standard error for torsional term
} Linear_FE_Model;

/* ______________________________________________________________________________ */
/* Energy Lookup Tables */

typedef struct energy_tables
{
    // e_vdW_Hb is sized for distances only up to NBC, the non-bond-cutoff
    Real e_vdW_Hb[NEINT][MAX_ATOM_TYPES][MAX_ATOM_TYPES];  // vdW & Hb energies
    // the other tables are sized for "unlimited" distances
    Real sol_fn[NDIEL];                            // distance-dependent desolvation function
    Real epsilon_fn[NDIEL];                        // distance-dependent dielectric function
    Real r_epsilon_fn[NDIEL];                      // r * distance-dependent dielectric function
    Boole is_hbond[MAX_ATOM_TYPES][MAX_ATOM_TYPES]; // for eintcalprint use
} EnergyTables;

/* ______________________________________________________________________________ */
/* Nonbonded pair parameters */
typedef struct nonbond_param
{
    int a1;           // ATM1
    int a2;           // ATM2
    int nonbond_type; // NBTYPE  0 = not 1_4     4 = is 1_4
    double desolv;
    double q1q2;      // product of atom partial charges
    int t1;           // TYPE1
    int t2;           // TYPE2
    Boole is_hbond;

    nonbond_param() : a1(0), a2(0) {}
} NonbondParam;

/* ______________________________________________________________________________ */
/* Statistics */

typedef struct statistics 
{
    int number;
    Real minimum;
    Real mean;
    Real maximum;
    /* Real standard_deviation; */
} Statistics;

/* ______________________________________________________________________________ */
/* Output_pop_stats (Output Population Statistics)   added 2010-05 M Pique */

typedef struct output_pop_stats
{
    unsigned int level;  // 0 means minimal, 1 means "basic" controlled by ngens, nevals:
    unsigned int everyNgens;
    unsigned int everyNevals;
} Output_pop_stats;


/* ______________________________________________________________________________ */
/* EnergyBreakdown */

typedef struct energy_breakdown
{
    Real e_inter_moving_fixed;      // (1)  // trilinterp( 0, true_ligand_atoms, ...)
    Real e_intra_moving_fixed_rec;  // (2)  // trilinterp( true_ligand_atoms, natom, ...)
    Real e_intra_moving_moving_lig; // (3)  // eintcal( 0, nb_array[0], ...)            // nb_group_energy[INTRA_LIGAND]
    Real e_inter_moving_moving;     // (4)  // eintcal( nb_array[0], nb_array[1], ...)  // nb_group_energy[INTER]
    Real e_intra_moving_moving_rec; // (5)  // eintcal( nb_array[1], nb_array[2], ...)  // nb_group_energy[INTRA_RECEPTOR]

    Real e_inter;                   // total    intermolecular energy = (1) + (4)
    Real e_intra;                   // total    intramolecular energy = (3)  +  (2) + (5)
    Real e_intra_lig;               // ligand   intramolecular energy = (3)
    Real e_intra_rec;               // receptor intramolecular energy = (2) + (5)

    Real e_torsFreeEnergy;          // empirical torsional free energy penalty
    Real e_unbound_internal_FE;     // computed internal free energy of the unbound state
    Real deltaG;                    // estimated change in free energy upon binding

} EnergyBreakdown;

/* component-level energy breakdowns  -  can be per-atom or for totals */
typedef struct energy_component
{
   Real total;
   Real elec;
   Real vdW_Hb;
   Real desolv;
   Real vdW;  // not available for grid-based values
   Real Hb;   // not available for grid-based values
   } EnergyComponent;
/* component-level energy breakdowns for each of the five atomic interaction groups in AD4 */
typedef struct group_energy {
   EnergyComponent inter_moving_fixed;      // true lig  .   rec maps  : set in trilinterp
   EnergyComponent intra_moving_fixed_rec;  // flex res  .   rec maps  : set in trilinterp
   EnergyComponent intra_moving_moving_lig; // true lig  x   true lig  : set in eintcal
   EnergyComponent inter_moving_moving;     // true lig  x   flex res  : set in eintcal
   EnergyComponent intra_moving_moving_rec; // flex res  x   flex res  : set in eintcal
   } GroupEnergy;

#endif
/* EOF */
