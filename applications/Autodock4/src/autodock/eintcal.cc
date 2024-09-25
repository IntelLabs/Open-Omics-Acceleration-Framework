/*

 $Id: eintcal.cc,v 1.33 2014/06/12 01:44:07 mp Exp $

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

#include <stdio.h>
#include <math.h>
#include "eintcal.h"
#include "constants.h"
#include "distdepdiel.h"
#include "stop.h"


extern Linear_FE_Model AD4;


#ifndef EINTCALPRINT

// Calculate internal energy
Real eintcal( const NonbondParam * const nonbondlist,
              const EnergyTables  *ptr_ad_energy_tables,
              const Real tcoord[MAX_ATOMS][SPACE],
              const int           Nnb,
	      int Nnb_array[3],
  	      GroupEnergy * group_energy, // sets structure internals unless NULL
              const Boole         B_calcIntElec,
              const Boole         B_include_1_4_interactions,
              ConstReal  scale_1_4,
              const Real qsp_abs_charge[MAX_ATOMS],
              const Boole B_use_non_bond_cutoff,
              const Boole B_have_flexible_residues,  // if the receptor has flexibile residues, this will be set to TRUE
	      const int outlev,
	      FILE *logFile
             )

#else 

// eintcalPrint [


// Calculate internal energy and print out a detailed report
Real eintcalPrint( const NonbondParam * const nonbondlist,
                   const EnergyTables  *ptr_ad_energy_tables,
                   const Real tcoord[MAX_ATOMS][SPACE],
                   const int           Nnb,
	           int Nnb_array[3],
  	           GroupEnergy *group_energy, // sets structure components if not NULL
                   const Boole         B_calcIntElec,
                   const Boole         B_include_1_4_interactions,
                   ConstReal  scale_1_4,
                   const Real qsp_abs_charge[MAX_ATOMS],
                   const Boole B_use_non_bond_cutoff,
                   const Boole B_have_flexible_residues, // if the receptor has flexibile residues, this will be set to TRUE
		   const int natom,
		   const int type[],
		   char const atom_type_name[MAX_MAPS][3],
		   const int outlev,
		   FILE *logFile
                  )
// eintcalPrint ]

#endif 

/* *****************************************************************************/
/*       Name: eintcal                                                         */
/*   Function: Calculate the Internal Energy of the Small Molecule.            */
/*             Accelerated non-square-rooting, dx,dy,dz version.               */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/* ____________________________________________________________________________*/
/*    Authors: Garrett M. Morris, TSRI                                         */
/*             David Goodsell, UCLA                                            */
/*       Date: 16/03/94                                                        */
/* ____________________________________________________________________________*/
/*     Inputs: nonbondlist, ptr_ad_energy_tables, tcoord, type, Nnb            */
/*    Returns: total_e_total                                                */
/*    Globals: NEINT, MAX_ATOMS, SPACE                                         */
/* ____________________________________________________________________________*/
/*  Modification Record                                                        */
/*  Date     Inits   Comments                                                  */
/*  07/05/92 DSG     Original FORTRAN                                          */
/*  15/05/92 GMM     Translated into C                                         */
/*  15/05/92 GMM     hypotenuse macro                                          */
/*  19/11/93 GMM     Accelerated non-square-rooting version.                   */
/*  16/03/94 GMM     Accelerated dx,dy,dz version.                             */
/*  10/02/04 GMM     Reduced NBC from 64.0 to 8.0                              */
/*  04/03/05 GMM     Added the new internal desolvation term                   */
/* *****************************************************************************/

{

  // if r is less than the non-bond-cutoff, 
  //  -OR-
  // If we are computing the unbound conformation then we ignore the non bond cutoff, NBC
        // if we have defined USE_8A_CUTOFF, then NBC = 8
        const double nbc2 = B_use_non_bond_cutoff ? NBC2 : 999 * 999;


    // strutures for tallying and reporting energy components by atom group
    static EnergyComponent zero_components; // always all zero

    // totals over entire system - all groups
    double total_e_total=0.0L; // total_e_total = eint
#ifdef EINTCALPRINT
    double total_e_elec=0.0L;
    double total_e_vdW_Hb=0.0L;
    double total_e_vdW=0.0L;
    double total_e_Hb=0.0L;
    double total_e_desolv=0.0L;

    double peratom_e_elec[MAX_ATOMS];
    double peratom_e_vdW[MAX_ATOMS];
    double peratom_e_Hb[MAX_ATOMS];
    double peratom_e_desolv[MAX_ATOMS];
    for(int a=0;a<MAX_ATOMS;a++) peratom_e_elec[a]=peratom_e_vdW[a]=peratom_e_Hb[a]=peratom_e_desolv[a]=0;
#endif

    int nb_group_max;

    // By default, we have one nonbond group, (1) intramolecular in the ligand
    // If we have flexible residues, we need to consider three groups of nonbonds:
    // (1) intramolecular in the ligand, (2) intermolecular and (3) intramolecular in the receptor
    nb_group_max = (B_have_flexible_residues) ? 3: 1;

    // Loop over the nonbonding groups --
    // Either (intramolecular ligand nonbonds)
    // or (intramolecular ligand nonbonds, intermolecular nonbonds, and intramolecular receptor nonbonds)
    for (int nb_group = 0;  nb_group < nb_group_max;  nb_group++) {
        int inb_from, inb_to;
        EnergyComponent grouptotal=zero_components;  // component totals for current group (0, 1, or 2)

#ifdef EINTCALPRINT
        if (nb_group == 0) {
            pr(logFile, "\n\n\t\tLigand Intramolecular Energy Analysis\n");
            pr(logFile,     "\t\t=====================================\n\n");
        }
        if (nb_group == 1) {
            pr(logFile, "\n\n\t\tLigand-Receptor Moving-Atom Intermolecular Energy Analysis\n");
            pr(logFile,     "\t\t==========================================================\n\n");
        }
        if (nb_group == 2) {
            pr(logFile, "\n\n\t\tReceptor Moving-Atom Intramolecular Energy Analysis\n");
            pr(logFile,     "\t\t===================================================\n\n");
        }
#define H1 "Non-bond  Atom1-Atom2  Distance   Total  "
#define U1 "________  ___________  ________   ______ "
#define H2 "     vdW        Hb     Desolv     Sol_fn   Type Dielectric"
#define U2 "  ________ ________  ________   ________   ____ __________"

        pr( logFile, "%s", H1);
        if (B_calcIntElec) pr( logFile, "      Elec");
        pr( logFile, "%s\n", H2);

        pr( logFile, "%s", U1);
        if (B_calcIntElec) pr( logFile, " _________");
        pr( logFile, "%s\n", U2);

#endif

        if (nb_group == 0) inb_from = 0;
        else inb_from = Nnb_array[nb_group-1];

        inb_to   = Nnb_array[nb_group];

        // Loop over the non-bonds in this nonbond "group", "inb",
        for (int inb = inb_from;  inb < inb_to;  inb++) {

	    int a1, a2;
	    int index_lt_NDIEL;

	    // energy components for single non-bond interaction:
            double e_total;  // e_total = epair (total)
	    double e_elec;
            double e_desolv;    // e_desolv = dpair
	    double e_vdW_Hb=0;
#ifdef EINTCALPRINT
	    double e_vdW=0;
	    double e_Hb=0;
#endif

	    double dx, dy, dz;
	    register double r2;

	    int nonbond_type; // if = 4, it is a 1_4;  otherwise it is another kind of nonbond

            a1 = nonbondlist[inb].a1;
            a2 = nonbondlist[inb].a2;

            dx = tcoord[a1][X] - tcoord[a2][X];
            dy = tcoord[a1][Y] - tcoord[a2][Y];
            dz = tcoord[a1][Z] - tcoord[a2][Z];

            // Calculate the van der Waals and/or H-bonding energy & the desolvation energy.
            //|
            //| desolvation energy = sol_fn[dist] * ( rec.vol * (lig.solpar + qsolpar * |lig.charge|)
            //|                                     + lig.vol * (rec.solpar + qsolpar * |rec.charge|) );
            //|

            r2 = sqhypotenuse(dx,dy,dz); // r2, the square of the separation between the atoms a1 and a2 in this non-bond, inb, 
            r2 = clamp(r2, (RMIN_ELEC*RMIN_ELEC));

	    // convert real-valued distance to an index for energy lookup tables
#ifndef NOSQRT 
            // Use square-root, slower...
            const int index = Ang_to_index(sqrt(r2)); 
#else   
            //  Non-square-rooting version, faster...
            const int index = SqAng_to_index(r2);
#endif  // NOSQRT

            index_lt_NDIEL = BoundedNdiel(index);  // guarantees that index_lt_NDIEL is never greater than (NDIEL - 1)

	    nonbond_type = nonbondlist[inb].nonbond_type;
            double nb_desolv = nonbondlist[inb].desolv;

            if (B_calcIntElec) {
                //  Calculate  Electrostatic  Energy
                double r_dielectric = ptr_ad_energy_tables->r_epsilon_fn[index_lt_NDIEL];
                e_elec = nonbondlist[inb].q1q2 * r_dielectric;
            }
	    else e_elec = 0;

	    e_total = e_elec;

            e_desolv = ptr_ad_energy_tables->sol_fn[index_lt_NDIEL] * nb_desolv;
            if  ( r2 < nbc2 ) {   
		int t1, t2; 
		t1 = nonbondlist[inb].t1; // t1 is a map_index
		t2 = nonbondlist[inb].t2; // t2 is a map_index
		 int index_lt_NEINT = BoundedNeint(index);  // guarantees that index_lt_NEINT is never greater than (NEINT - 1) (scaled NBC, non-bond cutoff)
		e_vdW_Hb= ptr_ad_energy_tables->e_vdW_Hb[index_lt_NEINT][t2][t1];
#ifdef EINTCALPRINT
		if( nonbondlist[inb].is_hbond ) e_Hb = e_vdW_Hb;
		else e_vdW = e_vdW_Hb;
#endif
                if (B_include_1_4_interactions && nonbond_type==4 ) {
                    //| Compute a scaled 1-4 interaction,
		    e_vdW_Hb *= scale_1_4;
		    e_desolv *= scale_1_4;
#ifdef EINTCALPRINT
		    e_vdW *=  scale_1_4;
		    e_Hb  *=  scale_1_4;
#endif
		}
                e_total += e_vdW_Hb + e_desolv;
	   }
	   else e_total += e_desolv; // no NBC-based cutoff for desolvation


          total_e_total += e_total;
#ifdef EINTCALPRINT // eintcalPrint [
          total_e_vdW_Hb   += e_vdW_Hb;
          total_e_vdW      += e_vdW;
          total_e_Hb       += e_Hb;
          total_e_desolv   += e_desolv;
          total_e_elec     += e_elec;
          double dielectric = ptr_ad_energy_tables->epsilon_fn[index_lt_NDIEL];

	  peratom_e_vdW[a1] += e_vdW/2;
	  peratom_e_vdW[a2] += e_vdW/2;
	  peratom_e_Hb[a1]  += e_Hb/2;
	  peratom_e_Hb[a2]  += e_Hb/2;
	  peratom_e_desolv[a1] += e_desolv/2;
	  peratom_e_desolv[a2] += e_desolv/2;
	  peratom_e_elec[a1] += e_elec/2;
	  peratom_e_elec[a2] += e_elec/2;


          pr( logFile, " %6d   %5d-%-5d %8.4lf %+9.4lf ",
                    (int)(inb+1), (int)(a1+1), (int)(a2+1), (double)sqrt(r2), (double)e_total);

          if (B_calcIntElec) pr( logFile, " %+9.4lf", (double)e_elec);

          pr( logFile, " %+9.4lf %+9.4lf %+9.4lf  %+9.4lf   %d  %8.3lf\n", 
                    (double)e_vdW, (double) e_Hb, (double)e_desolv, 
                    (double)ptr_ad_energy_tables->sol_fn[index_lt_NDIEL], 
		    (int)nonbond_type, (double)dielectric);

#endif // eintcalPrint ]
	if(group_energy!=NULL) {
		grouptotal.total += e_total;
		grouptotal.elec += e_elec;
		grouptotal.vdW_Hb += e_vdW_Hb;
#ifdef EINTCALPRINT
		grouptotal.vdW += e_vdW;
		grouptotal.Hb += e_Hb;
#endif
		grouptotal.desolv += e_desolv;
		}

        } //  inb -- next non-bond interaction

	// note: next operations always occur in specified order
	if(group_energy!=NULL)  switch ( nb_group ) {
	case INTRA_LIGAND:  // [0] Intramolecular energy of ligand
            group_energy->intra_moving_moving_lig = grouptotal;
	    break;
	case INTER: // [1]  intermolecular energy
            group_energy->inter_moving_moving = grouptotal;
	    break;
        case INTRA_RECEPTOR:  // [2]  intramolecular energy of receptor
            group_energy->intra_moving_moving_rec = grouptotal;
	    break;
	default: stop("bug: bad group index in eintcal");
        }

    } // nb_group -- intra lig, inter, intra rec


#ifdef EINTCALPRINT
    if (B_calcIntElec) {
        pr( logFile, "                                ________   ________  ________  ________  ________\n");
        pr( logFile, "Total                          %+9.4lf  %+9.4lf %+9.4lf  %+9.4lf %+9.4lf\n", total_e_total, total_e_elec, total_e_vdW, total_e_Hb, total_e_desolv);
        pr( logFile, "                                ________  ________  ________  ________ ________\n");
        pr( logFile, "                                   Total       Elec       vdW        Hb    Desolv\n");
    } else {
        pr( logFile, "                                ________   ________  ________ ________\n");
        pr( logFile, "Total                          %+9.4lf  %+9.4lf %+9.4lf %+9.4lf\n", total_e_total, total_e_vdW, total_e_Hb, total_e_desolv);
        pr( logFile, "                                ________   ________  ________ ________\n");
        pr( logFile, "                                   Total        vdW        Hb    Desolv\n");
    }
#endif

#ifdef EINTCALPRINT
            pr(logFile, "\n\n\tPer-atom Intramolecular Energy Analysis\n");
            pr(logFile,     "\t=======================================\n\n");
#define PAH1 "Atom Type   Total  "
#define PAU1 "____ ___  ________ "
#define PAH2 "    vdW        Hb     Desolv "
#define PAU2 " ________  ________  ________"

        pr( logFile, "%s", PAH1);
        if (B_calcIntElec) pr( logFile, "   Elec   ");
        pr( logFile, "%s\n", PAH2);

        pr( logFile, "%s", PAU1);
        if (B_calcIntElec) pr( logFile, " ________ ");
        pr( logFile, "%s\n", PAU2);

	for(int a=0;a<natom;a++) {
	 // print atom number, atom type (chars), Total E, [optional] e_elec, vdW, Hb, Desolv
	pr( logFile, "%4d  %-2s %+9.4lf ", a+1, atom_type_name[type[a]], 
	  peratom_e_vdW[a]+peratom_e_Hb[a]+peratom_e_elec[a]+peratom_e_desolv[a]);
	if(B_calcIntElec) pr( logFile, "%+9.4lf ", peratom_e_elec[a]);
	pr( logFile, "%+9.4lf %+9.4lf %+9.4lf\n", 
	  peratom_e_vdW[a],peratom_e_Hb[a],peratom_e_desolv[a]);
	}
   pr(logFile,"\n");

#endif

#ifdef EINTCALPRINT
    pr( logFile, "\nTotal Intramolecular Interaction Energy   = %+.3lf kcal/mol\n", (double)total_e_total); // eintcalPrint
#endif

    return (Real) total_e_total;
}
/*  EOF */
