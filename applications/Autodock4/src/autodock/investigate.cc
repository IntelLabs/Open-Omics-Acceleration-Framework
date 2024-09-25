/*

 $Id: investigate.cc,v 1.36 2014/02/01 05:14:53 mp Exp $

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
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/param.h>
#include <time.h>
#include "structs.h"
#include "investigate.h"
#include "stop.h"

#define RANDOM_MODE 1
#define CHANGE_MODE 2

extern char *programname;


void investigate( const int   Nnb, int Nnb_array[3], GroupEnergy *group_energy,
                    const Real charge[MAX_ATOMS],
                    const Real abs_charge[MAX_ATOMS],
                    const Real qsp_abs_charge[MAX_ATOMS],
                    const Boole B_calcIntElec,
                    /* not const */ Real crd[MAX_ATOMS][SPACE], // modified in cnv_state_to_coords
                    const Real crdpdb[MAX_ATOMS][SPACE],

                    const EnergyTables *const ptr_ad_energy_tables,

                    const int   maxTests,
                #include "map_declare.h"
                    const int   natom,
                    const NonbondParam *const nonbondlist,
                    const int   ntor,
                    const int   tlist[MAX_TORS+1][MAX_ATOMS],
                    const int   type[MAX_ATOMS],
                    const Real vt[MAX_TORS][SPACE],
                    const Boole B_isGaussTorCon,
                    const unsigned short US_torProfile[MAX_TORS][NTORDIVS],
                    const Boole B_isTorConstrained[MAX_TORS],
                    const Boole B_ShowTorE,
                    /* not const */ unsigned short US_TorE[MAX_TORS],
                    /* not const */ Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
                    /* not const */ int   N_con[MAX_TORS], // modifed in mkRandomState
                    const Boole B_symmetry_flag,
                    const Boole B_unique_pair_flag,
                    const char  *const FN_rms_ref_crds,
                    const int   OutputEveryNTests,
                    const int   NumLocalTests,
                    ConstReal trnStep,
                    ConstReal torStep,
                    
                    const int   ignore_inter[MAX_ATOMS],
                    
                    const Boole         B_include_1_4_interactions,
                    ConstReal scale_1_4,
                    ConstReal scale_eintermol,


                    ConstReal unbound_internal_FE,
                    /* not const */ GridMapSetInfo *const info, // modified in mkRandomState
                    const Boole B_use_non_bond_cutoff,
                    const Boole B_have_flexible_residues, 
                    const Boole B_rms_heavy_atoms_only, 
                    const int h_index,
		    const int true_ligand_atoms,
		    const int outlev,
		    FILE *logFile)

{
    Boole B_outside = FALSE;

    int Itor = 0;
    register int Test = -1;
    int indx;
    int ref_natoms = -1;
    register int i = 0;
    //register int XYZ = 0;

    Real e = 0.;
    Real ref_crds[MAX_ATOMS][SPACE];
    Real rms;
    Real MaxRms = 20.0;
    Real RmsBinSize = 0.25;
    Real MinEnergyInRmsBin[NUMRMSBINS];
    int   NumberInRmsBin[NUMRMSBINS];
    int   NumberRandomInRmsBin[NUMRMSBINS];
    int   NumberChangeInRmsBin[NUMRMSBINS];
    int   RmsBinNum = 0;
    //int   NumMaxRms = 0;
    register int NumOutside = 0;
    register int LocalTest = 0;
    int   mode;

    State sNow; /* qtnNow, torNow */


/*  Initialize
*/
    for (i=0; i<NUMRMSBINS; i++) {
        MinEnergyInRmsBin[i] = BIG;
        NumberInRmsBin[i] = 0;
        NumberRandomInRmsBin[i] = 0;
        NumberChangeInRmsBin[i] = 0;
    }
    sNow.ntor = ntor;

/*  Read in reference coordinates
*/
    if (strncmp(FN_rms_ref_crds,"unspecified filename",20) != 0) {
        if ((ref_natoms = getpdbcrds( FN_rms_ref_crds, ref_crds, logFile)) == -1) {
            fprintf( logFile, "%s: Problems while reading \"%s\".\n", programname, FN_rms_ref_crds);
            fprintf( logFile, "Will attempt to use the input PDBQ file coordinates as reference.\n");
        } else if (ref_natoms != natom) {
	    char msg[200];
	    sprintf(msg, 
	    "%s: ERROR!  Wrong number of atoms in reference structure.\n\
            Input PDBQ structure has %d atoms, but reference structure has %d atoms.\n",
	    programname, natom, ref_natoms);
            stop(msg); // exits
        }
    }

/* Begin investigating the force field,
   by recording the lowest energy for this rmsd
   from the crystal structure
   and from the minimized crystal structure.
*/

    pr( logFile, "\n\n\t\tBEGINNING INVESTIGATION OF FORCE FIELD\n");
    pr( logFile,     "\t\t______________________________________\n\n\n\n" );

    for ( Test = 0;  Test < maxTests;  Test++ ) {

        for (LocalTest = 0; LocalTest < NumLocalTests; LocalTest++, Test++ ) {

            if (LocalTest == 0) {
                mode = RANDOM_MODE;
            } else {
                mode = CHANGE_MODE;
            }

            do { /* while (rms > MaxRms); */
                do { /* while (B_outside); */
                    if (mode == RANDOM_MODE) {
                        sNow = mkRandomState( ntor, F_TorConRange, N_con, info );
                        if (outlev >= LOGRUNV) {
                            fprintf(logFile, "mkRandomState:  ");
                            writeState(logFile, sNow);
                            fflush(logFile);
                        }
                    } else {
                        sNow = changeState( sNow, trnStep, torStep,
                                              ntor, F_TorConRange, N_con);
                        if (outlev >= LOGRUNV) {
                            fprintf(logFile, "changeState:  ");
                            writeState(logFile, sNow);
                            fflush(logFile);
                        }
                    }
                    cnv_state_to_coords( sNow, vt, tlist, ntor, crdpdb, crd, natom,
		     true_ligand_atoms, outlev, logFile);
     
                    /* Check to see if any atom is outside the grid...  */
                    for (i = 0;  i < natom;  i++) {
                        B_outside= is_out_grid_info(crd[i][X], crd[i][Y], crd[i][Z]);
                        if ( B_outside ) {  /* Outside grid! */
                            ++NumOutside;
                            if (mode == CHANGE_MODE) {
                                /* changing pushes out of grid, so switch mode*/
                                mode = RANDOM_MODE;
                            }
                            break;/*...out of i*/
                        }/*outside*/
                    }/*for atoms i*/
                    /* If an atom is outside, do again */
                } while (B_outside);
                /* Now, ligand is inside grid */
                /* Calculate RMSD from reference structure */
                rms = getrms( crd, ref_crds, B_symmetry_flag, B_unique_pair_flag, natom, type, B_rms_heavy_atoms_only, h_index);
            } while (rms > MaxRms);
            /* Calculate Energy of System, */
            e = scale_eintermol * trilinterp( 0, natom, crd, charge, abs_charge, type, map, info, 
                ignore_inter, NULL, NULL, NULL_ENERGY_BREAKDOWN)
                 + eintcal( nonbondlist, ptr_ad_energy_tables, crd,
		     Nnb, Nnb_array, NULL_GROUP_ENERGY,
                     B_calcIntElec, B_include_1_4_interactions,
                     scale_1_4, qsp_abs_charge, 
                     B_use_non_bond_cutoff, B_have_flexible_residues, 
		     outlev, logFile);
            if (B_isGaussTorCon) {
                for (Itor = 0; Itor < ntor; Itor++) {
                    if (B_isTorConstrained[Itor] == 1) {
                        indx = RadiansToDivs( sNow.tor[Itor] );
                        if (B_ShowTorE) {
                            e += (Real)( US_TorE[Itor] 
                                          = US_torProfile[Itor][indx] );
                        } else {
                            e += (Real)US_torProfile[Itor][indx];
                        }
                    }
                }
            }
            /* Update minimum energy for this RMSD bin */
            RmsBinNum = (int)(rms / RmsBinSize);
            ++NumberInRmsBin[RmsBinNum];
            if (mode == RANDOM_MODE) {
                ++NumberRandomInRmsBin[RmsBinNum];
            } else if (mode == CHANGE_MODE) {
                ++NumberChangeInRmsBin[RmsBinNum];
            }
            if (e <= MinEnergyInRmsBin[RmsBinNum]) {
                MinEnergyInRmsBin[RmsBinNum] = e;
            }
            /* Output if it is time, */
            if (outlev >= LOGBASIC ) {
                if ((Test+1)%OutputEveryNTests == 0) {
                    fprintf(logFile, "NumberOfTests= %d\n", Test+1);
                    fprintf(logFile, "-------------\n");
                    for (i=0; i<NUMRMSBINS; i++) {
                        fprintf(logFile, "%2d %5.2f-%5.2f:  %9.2f\t%7d\t%7d\t%7d\n", i+1, i*RmsBinSize, (i+1)*RmsBinSize, MinEnergyInRmsBin[i], NumberInRmsBin[i], NumberRandomInRmsBin[i], NumberChangeInRmsBin[i]);
                    }
                    fprintf(logFile, "\n");
                    fprintf(logFile, "NumOutside= %d\n", NumOutside);
                    fprintf(logFile, "\n");
                    fprintf(logFile, "\n");
                    fflush(logFile);
                }
            }

        } /*LocalTest*/
    } /* Loop over Test */
}
/* EOF */
