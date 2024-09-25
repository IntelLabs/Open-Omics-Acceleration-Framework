/*

 $Id: getInitialState.cc,v 1.42 2014/07/02 00:01:21 mp Exp $

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
#include <stdio.h>
#include <sys/types.h>
#include <sys/param.h>
#include <time.h>
#include "timesys.h"
#include "getInitialState.h"
#include "trilinterp.h"
#include "calculateEnergies.h"


extern char *programname;


void getInitialState(

            /* not const */ Real *const Addr_e0total,
            ConstReal  e0max,

            /* not const */ State *const sInit, /* was qtn0[QUAT] and tor0[MAX_TORS] */
            /* not const */ State *const sMinm, /* was qtnMin[QUAT] and torMin[MAX_TORS] */
            /* not const */ State *const sLast, /* was qtnLast[QUAT] and torLast[MAX_TORS] */

            const Boole B_RandomTran0,
            const Boole B_RandomQuat0,
            const Boole B_RandomDihe0,

            const Real charge[MAX_ATOMS],
            const Real abs_charge[MAX_ATOMS],
            const Real qsp_abs_charge[MAX_ATOMS],
            /* not const */ Real crd[MAX_ATOMS][SPACE],
            const Real crdpdb[MAX_ATOMS][SPACE],
            const char  atomstuff[MAX_ATOMS][MAX_CHARS],
            /* not const */ EnergyComponent	*peratomE,        // output if not NULL - intermolecular energies

            const EnergyTables *const ptr_ad_energy_tables,

            const Boole B_calcIntElec,
                #include "map_declare.h"
            const int   natom,
            const int   Nnb,
	    int Nnb_array[3],
	    GroupEnergy *group_energy,
	    const int true_ligand_atoms,
            const NonbondParam *const nonbondlist,
            const int   ntor,
            const int   tlist[MAX_TORS+1][MAX_ATOMS],
            const int   type[MAX_ATOMS],
            const Real vt[MAX_TORS][SPACE],
            const int   irun, // 0-origin, includes nconf
            const int   MaxRetries,

            ConstReal  torsFreeEnergy,

            const int   ligand_is_inhibitor,

            const int   ignore_inter[MAX_ATOMS],

            const Boole         B_include_1_4_interactions,
            ConstReal  scale_1_4,
	    ConstReal            scale_eintermol,  // input  scaling factor for intermolecular energies


            ConstReal  unbound_internal_FE,

            const GridMapSetInfo *const info,
            const Boole B_use_non_bond_cutoff,
            const Boole B_have_flexible_residues,
            const Unbound_Model ad4_unbound_model,
            const int   outlev,
	    FILE *logFile
           )

{
    Real e0total = 0.;
    Real e0inter = 0.;
    Real e0intra = 0.;
    Real e0min = BIG;
    int   retries = 0;
    register int i = 0;
    Clock  initStart;
    Clock  initEnd;
    struct tms tms_initStart;
    struct tms tms_initEnd;
    EnergyBreakdown eb;


    /* Store time at start of initialization process... */

    initStart = times( &tms_initStart );

    retries = 0;
    if ((B_RandomTran0) || (B_RandomQuat0) || (B_RandomDihe0)) {
        do {
            /*
            ** while e0total, initial energy, is too high,
            */

            /*
            ** Initialize all state variables...
            */
            if (B_RandomTran0) {
                sInit->T.x = random_range( info->lo[X], info->hi[X] );
                sInit->T.y = random_range( info->lo[Y], info->hi[Y] );
                sInit->T.z = random_range( info->lo[Z], info->hi[Z] );
                if (outlev >= LOGRUNV) {
                    pr( logFile, "Random initial translation,  tran0 %.3f %.3f %.3f\n", sInit->T.x, sInit->T.y, sInit->T.z);
                }
            }/*if*/
            if (B_RandomQuat0) {
                sInit->Q = randomQuat(); // generate a uniformly-distributed quaternion

                if (outlev >= LOGRUNV ) {
                    pr( logFile, "Random initial quaternion,  quaternion0 %.3f %.3f %.3f %.3f\n", sInit->Q.x, sInit->Q.y, sInit->Q.z,  sInit->Q.w );
                     // convert from qx,qy,qz,qw to axis-angle (nx,ny,nz,ang)
		    AxisAngle aa = QuatToAxisAngle(sInit->Q);
                    pr( logFile, "Random initial axis-angle,  axisangle0 %.3f %.3f %.3f %.1f\n", aa.nx, aa.ny, aa.nz, RadiansToDegrees( aa.ang ) );
                }
            }/*if*/
            if ( B_RandomDihe0 && (ntor > 0) ) {
                if (outlev >= LOGRUNV ) {
                    pr( logFile, "Random initial torsions, ndihe = %d\ndihe0 = ", ntor);
                }
                sInit->ntor = ntor;
                for (i=0; i<ntor; i++) {
                    sInit->tor[i] = random_range(-180.,180.);
                    if (outlev >= LOGRUNV ) {
                        pr( logFile, "%7.2f ", sInit->tor[i] ); /*in degrees*/
                    }
                    sInit->tor[i] = DegreesToRadians( sInit->tor[i] ); /*now in radians*/
                }
                if (outlev >= LOGRUNV ) {
                    pr( logFile, "\n");
                }
            }/*if*/

            copyState( sMinm, *sInit );
            copyState( sLast, *sInit );

    /* _________________________________________________________________________
    **
    ** Initialize the automated docking simulation,
    ** _________________________________________________________________________
    */
            initautodock( atomstuff, crd, crdpdb, 
                natom, ntor, sInit, tlist, vt, true_ligand_atoms,outlev, info);
            
            e0inter = scale_eintermol * trilinterp( 0, natom, crd, charge, abs_charge, type, map, 
                        info, ignore_inter, peratomE, NULL,
                        NULL_ENERGY_BREAKDOWN);
            e0intra = eintcal( nonbondlist, ptr_ad_energy_tables, crd,
			  Nnb, Nnb_array, group_energy,
                          B_calcIntElec, B_include_1_4_interactions,
                          scale_1_4, qsp_abs_charge,
                          B_use_non_bond_cutoff, B_have_flexible_residues,
			  outlev, logFile);
            e0total = e0inter + e0intra;

            if (e0total < e0min) {
                /*
                ** This energy is lower so update
                ** the initialization minimum-energy state variables,
                */
                e0min = e0total;
                copyState( sMinm, *sInit );
            }

            /*
            if ( (e0total > e0max) && ((!B_RandomTran0)||(!B_RandomQuat0)) ) {
                B_RandomTran0 = TRUE;
                B_RandomQuat0 = TRUE;
            }
            */

            ++retries;
            if ((retries > 0) && (retries < MaxRetries) &&
                (e0total > e0max) && (outlev >= LOGRUNVV)) {

                pr(logFile, "Initial total energy, e0total = %.3f, too high!\n", e0total);
                pr(logFile, "Number of attempts = %d (run %d)\n\n", retries, irun+1);
                pr(logFile, "Will try again...\n");

            } else if (retries >= MaxRetries) {

		if(outlev >= LOGRUNVV)
                pr( logFile, "Sorry, too many retries (%d).  Continuing...\n\nWill use the state with the lowest energy found, %.2f\n\n", MaxRetries, e0min);
                e0total = e0min;
                copyState( sInit, *sMinm );

                break;
            }
            fflush( logFile );
        } while ( e0total > e0max );
    } // endif

    cnv_state_to_coords( *sInit, vt, tlist, ntor, crdpdb, crd, natom,
     true_ligand_atoms, outlev, logFile);

    eb = calculateBindingEnergies( natom, ntor, unbound_internal_FE, torsFreeEnergy, B_have_flexible_residues,
         crd, charge, abs_charge, type, map, info,
         ignore_inter, peratomE, NULL,
         nonbondlist, ptr_ad_energy_tables,
	 Nnb, Nnb_array, group_energy, true_ligand_atoms,
         B_calcIntElec, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, B_use_non_bond_cutoff, ad4_unbound_model,
	 outlev, logFile);

    copyState( sMinm, *sInit );
    copyState( sLast, *sInit );

    if(outlev>=LOGRUNVV) {
      prInitialState( &eb, natom, true_ligand_atoms, crd, atomstuff, type, peratomE, charge, 
       ligand_is_inhibitor, B_have_flexible_residues, ad4_unbound_model,
	outlev, logFile);
      initEnd = times( &tms_initEnd );
      pr(logFile, "Number of initialization attempts = %d (run %d)\n", retries, irun+1);
      pr(logFile, "Time spent initializing: (Real, CPU, System): ");
      timesys( initEnd - initStart, &tms_initStart, &tms_initEnd, logFile);
    }

    *Addr_e0total = eb.e_inter + eb.e_intra;

}
/* EOF */
