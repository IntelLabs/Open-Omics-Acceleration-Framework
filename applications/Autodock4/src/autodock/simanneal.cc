/*

 $Id: simanneal.cc,v 1.50 2014/07/10 23:27:37 mp Exp $

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
#include <sys/types.h>
#include <sys/param.h>
#include <time.h>
#include "simanneal.h"
#include "energy.h"
#include "timesys.h"


extern char *programname;


void simanneal ( int   *const Addr_nconf, /* absolute serial number of "run", modified here */
                const int   Nnb,
		int Nnb_array[3],
		GroupEnergy *group_energy,
		const int true_ligand_atoms,
                ConstReal WallEnergy,
                const char  atomstuff[MAX_ATOMS][MAX_CHARS],
                Real charge[MAX_ATOMS],
                Real abs_charge[MAX_ATOMS],
                Real qsp_abs_charge[MAX_ATOMS],
                const Boole B_calcIntElec,
                Real crd[MAX_ATOMS][SPACE],
                Real crdpdb[MAX_ATOMS][SPACE],
                const char  *const FN_dpf,
        
                EnergyTables *const ptr_ad_energy_tables,

                /* not const */ Real econf[MAX_RUNS],
                const Boole B_either,
		EnergyComponent	*peratomE,
                const int   NcycMax,
                const int   nruns,
		FourByteLong runseed[][2], /* [MAX_RUNS] */
                const Clock& jobStart,
                #include "map_declare.h"
                const int   naccmax,
                const int   natom,
                NonbondParam *const nonbondlist,
                const int   nrejmax,
                const int   ntor,

                /* not const */ State sInit, /* tor0, qtn0 */
                /* not const */ State sHist[MAX_RUNS], /* was qtnHist, torHist */

                ConstReal qtwFac,
                const Boole B_qtwReduc,
                ConstReal qtwStep0,
                const Boole B_selectmin,
                const char  *const FN_ligand,
                const Real lig_center[SPACE],
                ConstReal RT0,
                const Boole B_RTChange,
                ConstReal RTFac,
                const struct tms& tms_jobStart,
                const int   tlist[MAX_TORS][MAX_ATOMS],
                ConstReal torFac,
                const Boole B_torReduc,
                ConstReal torStep0,
                const char  *const FN_trj,
                const int   trj_cyc_max,
                const int   trj_cyc_min,
                const int   trj_freq,
                ConstReal trnFac,
                const Boole B_trnReduc,
                ConstReal trnStep0,
                const int   type[MAX_ATOMS],
                Real vt[MAX_TORS][SPACE],
                const Boole B_writeTrj,
                const Boole B_constrain,
                const int   atomC1,
                const int   atomC2,
                ConstReal sqlower,
                ConstReal squpper,
                const Boole B_linear_schedule,
                ConstReal RTreduc,
                /*Real maxrad,*/
                const Boole B_watch,
                const char  *const FN_watch,
                const Boole B_isGaussTorCon,
                /* not const */ unsigned short US_torProfile[MAX_TORS][NTORDIVS],
                const Boole B_isTorConstrained[MAX_TORS],
                const Boole B_ShowTorE,
                /* not const */ unsigned short US_TorE[MAX_TORS],
                const Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
                const int N_con[MAX_TORS],
                const Boole B_RandomTran0,
                const Boole B_RandomQuat0,
                const Boole B_RandomDihe0,
                ConstReal e0max,
        
                ConstReal torsFreeEnergy,
                
        const int   MaxRetries,
                
        const int   ligand_is_inhibitor,
                
        const int   ignore_inter[MAX_ATOMS],
        
        const Boole         B_include_1_4_interactions,
        ConstReal scale_1_4,
        ConstReal scale_eintermol,
        
        const ParameterEntry parameterArray[MAX_ATOM_TYPES],

        ConstReal unbound_internal_FE,

        GridMapSetInfo *const info,
        const Boole B_use_non_bond_cutoff,
        const Boole B_have_flexible_residues, 
        const char PDBQT_record[MAX_RECORDS][LINE_LEN],
        const Unbound_Model ad4_unbound_model,
	const int outlev,
	FILE *logFile
        )

{
    Real local_unbound_internal_FE = unbound_internal_FE;

    FILE *FP_trj;
        FP_trj = NULL;

    State sNow; /* qtnNow, torNow */
    State sChange; /* qtnChange, torChange */
    State sLast; /* qtnLast, torLast */
    State sMin; /* qtnMin, torMin */
    State sSave; /* qtnSave, torSave */

    Real d[SPACE];
    Real e = 0.;
    Real e0 = 0.;
    Real einter = 0.;
    Real eintra = 0.;
    Real eLast = 0.;
    Real eMin = BIG_ENERGY;
    Real etot = 0.0;
    Real inv_RT = 0.;
    Real qtwStep;
    Real RT = 616.;
    Real torTmp;
    Real torStep;
    Real trnStep;
    /* ** Real xloTrn; ** Real xhiTrn; ** Real yloTrn; ** Real yhiTrn; ** Real zloTrn; ** Real zhiTrn; ** Real lo[SPACE]; ** Real trnStepHi; ** Real qtwStepHi; ** Real torStepHi; */

    Boole B_inRange = FALSE;
    Boole B_outside = FALSE;
    Boole B_within_constraint = TRUE;

    int count = 0;
    int Itor = 0;
    int irun;
    int indx;

    register int i = 0;
    register int nacc = 0;
    register int nAcc = 0;
    register int nAccProb = 0;
    register int nedge = 0;
    register int nrej = 0;

    struct tms tms_cycStart;
    struct tms tms_cycEnd;
    struct tms tms_jobEnd;

    Clock cycStart;
    Clock cycEnd;
    Clock jobEnd;

    /* trnStepHi = HI_NRG_JUMP_FACTOR * trnStep; ** qtwStepHi = HI_NRG_JUMP_FACTOR * qtwStep; ** torStepHi = HI_NRG_JUMP_FACTOR * torStep; */
    /* lo[X] = xlo;  lo[Y] = ylo;  lo[Z] = zlo;*/
    /* xloTrn = xlo + maxrad; ** xhiTrn = xhi - maxrad; ** yloTrn = ylo + maxrad; ** yhiTrn = yhi - maxrad; ** zloTrn = zlo + maxrad; ** zhiTrn = zhi - maxrad; */
    
    // Initialise these State variables:
    initialiseState( &sNow );
    initialiseState( &sChange );
    initialiseState( &sLast );
    initialiseState( &sMin );
    initialiseState( &sSave );

/* Open the trajectory file for writing, =====================================*/

    if ( B_writeTrj ) {
                FP_trj = fopen(FN_trj, "w");
        if ( FP_trj == NULL ) {
	    char error_message[LINE_LEN];
            prStr( error_message, "\n%s: can't create trajectory file %s\n", programname, FN_trj);
            pr_2x( stderr, logFile, error_message );
            jobEnd = times( &tms_jobEnd );
            timesys( jobEnd - jobStart, &tms_jobStart, &tms_jobEnd, logFile);
            pr_2x( logFile, stderr, UnderLine );
            stop(error_message);
        }/*END PROGRAM*/
    }/*endif*/

/* Begin the automated docking simulation, ===================================*/

    if(outlev>LOGBASIC) {
    pr( logFile, "\n\n\t\tBEGINNING MONTE CARLO SIMULATED ANNEALING\n");
    pr( logFile, "     \t\t_________________________________________\n\n\n\n" );
    }


    int nconf = *(Addr_nconf); /* absolute serial number of this run */
    for ( irun = 0;  irun < nruns;  irun++ ) { /*===========================*/

	/* set RNG seed using global run number, see main.cc and com.cc */
	if(nconf==0&&irun==0) getsd(&runseed[nconf][0], &runseed[nconf][1]);
	else setsd(runseed[nconf+irun][0], runseed[nconf+irun][1]); 


        if (outlev >= LOGRUNV) 
            pr(logFile, "\n\tINITIALIZING AUTOMATED DOCKING SIMULATION\n" );
        pr( logFile, "Run: %d Seed: %ld %ld [ Run %d of %d SA ]\n", nconf+irun+1,
		(long)runseed[nconf+irun][0], (long)runseed[nconf+irun][1],
		irun+1, nruns);

        if ( B_writeTrj ) {
            pr( FP_trj, "ntorsions %d\nrun %d\n", ntor, nconf+irun+1 );
            fflush( FP_trj );
        }

        getInitialState( &e0, e0max,
                     &sInit, &sMin, &sLast, 
                     B_RandomTran0, B_RandomQuat0, B_RandomDihe0, 
                     charge, abs_charge, qsp_abs_charge, crd, crdpdb, atomstuff,
                     peratomE, ptr_ad_energy_tables, B_calcIntElec,
                     map, natom, Nnb, Nnb_array, group_energy, 
		     true_ligand_atoms, nonbondlist,
                     ntor, tlist, type, vt, nconf+irun, MaxRetries,
                     torsFreeEnergy, ligand_is_inhibitor,
                     ignore_inter,
                     B_include_1_4_interactions, scale_1_4, scale_eintermol,
                     unbound_internal_FE, info, 
                     B_use_non_bond_cutoff, B_have_flexible_residues,
                     ad4_unbound_model, outlev, logFile);

        /* Initialize the "annealing" temperature */
        RT = RT0;                

        if (RT <= APPROX_ZERO) {
            RT = 616.;
        }
        inv_RT = 1. / RT;

        eMin    = min(BIG_ENERGY, e0);
        eLast   = e0;
        trnStep = trnStep0;        /* translation*/
        qtwStep = qtwStep0;        /* quaternion angle, w*/
        torStep = torStep0;        /* torsion angles*/

        if (outlev >= LOGRUNV) {
            pr( logFile, "\n\n\t\tBEGINNING SIMULATED ANNEALING");
            pr( logFile, "\n\t\t_____________________________\n\n");
            pr( logFile, "\n      \t      \tMinimum     Average     | Acc/    Accepted:    Rejected:     |          |  xyz-Translation  |        Time:        \n");
              pr( logFile, "Run#  \tCycle:\tEnergy:     Energy:     |   /Rej: Down:  Up:   Total: Edge:  |   RT:    |   of Min.Energy   |  Real, CPU, System  \n" );
              pr( logFile, "______\t______\t___________ ___________ | ______ ______ ______ ______ ______ | ________ | _________________ | ____________________\n" );
/*                                          12345678901 12345678901   123456 123456 123456 123456 123456   12345678   12345 12345 12345
**                         "%D /%D\T%D /%D\T%+11.2F     %+11.2F       %6.2F  %6D    %6D    %6D    %6D      %8.1F     %5.2F %5.2F %5.2F   ",
**                         IRUN1, IRUNMAX, ICYCLE1, ICYCLEMAX, EmIN, ETOT/NTOT, QTNmIN[x], QTNmIN[y], QTNmIN[z], (NREJ!=0) ? (FLOAT)NACC/NREJ : -999., NACC, NREJ, NEDGE, rt
*/
        }
        fflush(logFile);
/*____________________________________________________________________________*/

        for ( int icycle = 0;  icycle < NcycMax;  icycle++ ) {

	    int ntot;
            cycStart = times( &tms_cycStart );

	    ntot = 1;
            B_inRange = (icycle >= trj_cyc_min) && (icycle <= trj_cyc_max);
            if ( B_writeTrj && B_inRange ) {
                pr( FP_trj, "cycle %d\ntemp %f\n", icycle+1, RT );
                fflush(  FP_trj  );
            }
            nAcc = nAccProb = nacc = nrej = nedge = 0;
            etot = eLast;

/*____________________________________________________________________________*/

            do {
                /*
                ** Do one Monte Carlo step,
                ** while the number of accepted steps < naccmax 
                ** and number of rejected steps < nrejmax...
                */

                mkNewState( &sNow, &sLast, &sChange,
                    vt, tlist, ntor, crd, crdpdb, natom,
                    trnStep,
                    qtwStep,
                    torStep,
                    F_TorConRange, N_con,
		    true_ligand_atoms,
		    outlev, logFile);

                if (B_constrain) {
		    Real rsqC1C2;
                    for (int xyz = 0;  xyz < SPACE;  xyz++) {
                        d[xyz] = crd[atomC1][xyz] - crd[atomC2][xyz];
                    }
                    rsqC1C2 = sqhypotenuse(d[X],  d[Y], d[Z]);
                    if (! (B_within_constraint = (rsqC1C2 > sqlower) && 
                                                 (rsqC1C2 < squpper))) {
                        copyState( &sNow, sLast );
                    }
                }/*if B_constrain*/

                if (B_within_constraint) {
                    /* 
                    ** Normally, this is true.
                    ** If the distance-constraint was set,
                    ** and was just violated, this is false.
                    */

                    for (i = 0;  i < natom;  i++) {
                        B_outside= is_out_grid_info(crd[i][X], crd[i][Y], crd[i][Z]);
                        if ( B_outside ) {
                            /*
                            ** Outside grid!
                            */
                            ++nedge;
                            ++nrej;
                            /*pr(logFile, "e"); / *###*/
                            /* etot += WallEnergy; ++ntot;*/
                            /*
                            ** Undo this move,
                            */
                            copyState( &sNow, sLast );
                            /*
                            ** Subtract just the translation;
                            ** hence "doubling back" from the edge.
                            */
                            sNow.T.x -= sChange.T.x;
                            sNow.T.y -= sChange.T.y;
                            sNow.T.z -= sChange.T.z;

                            if ( B_writeTrj && B_inRange && B_either && (++count == trj_freq)) {
                                count = 0;
                                e = eintra = WallEnergy;
                                output_state( FP_trj, sNow, ntor, nacc+nrej, e,
                                    eintra, (char)'e', B_watch, FN_watch, 
                                    atomstuff, natom, crd);
                            }/*writeTrj*/

                            break;/*...out of i*/

                        }/*outside*/
                    }/*for atoms i*/

                    if ( !B_outside ){ /*inside grid maps*/

                         /* Calculate Energy of System, =======================*/

                        /*
                        ** MORE ACCURATE METHOD, (SLOWER):
                        */
                        e = scale_eintermol * trilinterp( 0, natom, crd, charge, abs_charge, type, map, 
                                        info, ignore_inter, NULL, NULL,
                                        NULL_ENERGY_BREAKDOWN)
                           + (eintra = eintcal(nonbondlist, ptr_ad_energy_tables, crd, Nnb,
				   Nnb_array, group_energy,  // perhaps no breakdown needed
                                   B_calcIntElec, B_include_1_4_interactions,
                                   scale_1_4, qsp_abs_charge, 
                                   B_use_non_bond_cutoff, B_have_flexible_residues,
				   outlev, logFile)
                               );

                        if (B_isGaussTorCon) {
                            /*** This looks wrong... for (Itor = 0; Itor <= ntor; Itor++)   MP ***/
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

                        /* Apply the Metropolis energy test... ===============*/

                        if (e <= eLast) {
                            /*
                            **  Accept this move immediately.
                            */
                            ++nacc;
                            ++nAcc;
                            etot += (eLast = e);
                            ++ntot;
                            copyState( &sLast, sNow );

                            /* pr(logFile, "A"); / *###*/
                            if (e < eMin) {
                                /*
                                ** Update minimum-energy state variables,
                                */
                                eMin = e;
                                copyState( &sMin,  sNow );
                            }
                            if ( B_writeTrj && B_inRange && (++count == trj_freq)){
                                count = 0;
                                output_state( FP_trj, sNow, ntor, nacc+nrej, e,
                                    eintra, (char)'A', B_watch, FN_watch, 
                                    atomstuff, natom, crd);
                            }/*write trajectory*/

                        } else {

                            /* Probabilistic move. ===========================*/

                            if (exp((double)((eLast-e)*inv_RT))<local_random()){
                                /*
                                ** Failed the probability test. 
                                ** Reject this move.
                                */
                                ++nrej;
                                /* pr(logFile, "R"); / *###*/
                                copyState( &sNow, sLast );

                                if ( B_writeTrj && B_inRange && B_either && 
                                    (++count == trj_freq)) {
                                    count = 0;
                                    output_state( FP_trj, sNow, ntor, nacc+nrej,
                                        e, eintra, (char)'R', B_watch, FN_watch,
                                        atomstuff, natom, crd);
                                }/*write trajectory*/

                            } else {
                                /*
                                ** Passed the probability test.  
                                ** Accept this move.
                                ** A chance to escape a local minimum...
                                */
                                ++nacc;
                                ++nAccProb;
                                etot += e;
                                ++ntot;
                                /*pr(logFile, "a"); / *###*/

                                eLast = e;
                                copyState( &sLast, sNow );
                                if ( B_writeTrj && B_inRange && (++count == trj_freq))  {
                                    count = 0;
                                    output_state(FP_trj, sNow, ntor, nacc+nrej,
                                        e, eintra, (char)'a', B_watch, FN_watch,
                                        atomstuff, natom, crd);
                                }/*write trajectory*/
                            }/*passed Monte Carlo probablility test*/
                        }/*e > eLast, Probabilistic move...*/
                    }/*inside grid maps*/
                }/*within_constraint*/

            } while ( nacc < naccmax  &&  nrej < nrejmax );
/*____________________________________________________________________________*/

            if ((nacc == 0) && (nedge == nrejmax)) {
                /*
                **  Clear indication that ligand got stuck on an edge...
                */
                pr(logFile, "\n\n>>> Ligand appears to be stuck on an edge: forced to re-initialize. <<<\n\n");
                --icycle;
                getInitialState( &e0, e0max,
                         &sInit, &sMin, &sLast, 
                         TRUE, TRUE, TRUE, 
                         charge, abs_charge, qsp_abs_charge, crd, crdpdb, atomstuff,
                         peratomE, ptr_ad_energy_tables, B_calcIntElec,
                         map, natom, 
			 Nnb, Nnb_array, group_energy, true_ligand_atoms,
			 nonbondlist, ntor, tlist, type, vt, nconf+irun, MaxRetries,
                         torsFreeEnergy, ligand_is_inhibitor,
                         ignore_inter,
                         B_include_1_4_interactions, scale_1_4, scale_eintermol,
                         unbound_internal_FE,
                         info, B_use_non_bond_cutoff,
                         B_have_flexible_residues,
                         ad4_unbound_model,
			 outlev, logFile);

            } else {

                if ( B_trnReduc )  trnStep *= trnFac;       
                if ( B_qtwReduc )  qtwStep *= qtwFac;
                if ( B_torReduc )  torStep *= torFac;
                /*
                **  Output-level dependent diagnostics...
                */
                if (outlev >= LOGRUNV) {
                    /*pr(logFile, "\n"); / *###*/
                    pr( logFile, "%d \t%d /%d\t%+11.2f %+11.2f   %6.2f %6d %6d %6d %6d   %8.1f   %5.2f %5.2f %5.2f   ", nconf+irun+1, icycle+1, NcycMax, eMin, etot/ntot, (nrej!=0) ? (Real)nacc/nrej : 999.99, nAcc, nAccProb, nrej, nedge, RT, sMin.T.x, sMin.T.y, sMin.T.z );
                    cycEnd = times( &tms_cycEnd );
                    timesys( cycEnd - cycStart, &tms_cycStart, &tms_cycEnd, logFile);
                    if (outlev > LOGRUNV ) {
                        pr( logFile, "\tEnergy:   \tState:\n\t__________\t____________________________________________________________\nMinimum\t%+6.2f\t(%+.2f,%+.2f,%+.2f), q = [x,y,z,w] = [%5.1f deg, (%+.2f,%+.2f,%+.2f)],\n", eMin, sMin.T.x, sMin.T.y, sMin.T.z, sMin.Q.x, sMin.Q.y, sMin.Q.z, sMin.Q.w );
                        pr( logFile, "\nLast\t%+6.2f\t(%+.2f,%+.2f,%+.2f), q = [x,y,z,w] = [%5.1f deg, (%+.2f,%+.2f,%+.2f)],\n", eLast, sLast.T.x, sLast.T.y, sLast.T.z, sLast.Q.x, sLast.Q.y, sLast.Q.z, sLast.Q.w );
                        if (ntor > 0) {
                            pr( logFile, "Minimum:\t(" );
                            for (i=0; i<ntor; i++) {
                                pr( logFile, "%.1f%s ", RadiansToDegrees(sMin.tor[i]), (i < ntor-1)?",":" deg)" );
                            }
                            pr( logFile, "\nLast:\t(" );
                            for (i=0; i<ntor; i++) {
                                pr( logFile, "%.1f%s ", RadiansToDegrees(sLast.tor[i]), (i < ntor-1)?",":" deg)" );
                            }
                            pr( logFile, "\n" );
                        }
                        if ( B_trnReduc )
                            pr( logFile, "\nTranslation step size reduced; now =\t\t +/- %.2f A\n", trnFac);
                        if ( B_qtwReduc )
                            pr( logFile, "\nQuaternion Rotation step size reduced; now =\t +/- %.2f deg\n", qtwFac);
                        if ( B_torReduc )
                            pr( logFile, "\nTorsion step size reduced; now =\t\t +/- %.2f deg\n", torFac);
                    }/*outlev > LOGRUNV */
                    flushLog;
                }/*outlev >= LOGRUNV*/
                /*
                ** Reduce temperature,
                */
                if ( B_linear_schedule ) {
                    RT -= RTreduc;

                    if (RT <= APPROX_ZERO) {
                        inv_RT = 1. / APPROX_ZERO;
                    } else {
                        inv_RT = 1. / RT;
                    }
                } else if ( B_RTChange ) {
                    inv_RT = 1./( RT *= RTFac );
                }
                /*
                ** Start next cycle at minimum state?
                */
                if ( B_selectmin ) {
                    eLast = eMin;
                    copyState( &sLast, sMin );
                }
            } /* Prepares for next cycle */

        } /* icycle */
/*____________________________________________________________________________*/

        if ( B_selectmin ) {
            copyState( &sSave, sMin );
        } else {
            copyState( &sSave, sLast );
        }

        for (i=0; i<ntor; i++) {
            sSave.tor[i] = WrpRad( ModRad( sSave.tor[i] ) );
        }
   
	if(outlev>=LOGRUNV) pr( logFile, "\tFINAL DOCKED STATE\n" );
        if(outlev>=LOGRUNVV)pr( logFile, "Run Number %d\n", nconf+irun+1);
        if(outlev>=LOGRUNVV)pr( logFile, "Final Energy = %+.2f\n", eLast);
	pr(logFile,"Final-Value: %.3f\n", eLast);

	if(outlev>=LOGRUNV) {
        pr( logFile, "Final Translation = %.2f, %.2f, %.2f\n", sSave.T.x, sSave.T.y, sSave.T.z );
        pr( logFile, "Final Quaternion = ( %+.2f, %+.2f, %+.2f, %+.2f )\n", sSave.Q.x, sSave.Q.y, sSave.Q.z, sSave.Q.w );
	AxisAngle aa = QuatToAxisAngle( sSave.Q );
        pr( logFile, "Final Rotation Axis = ( %+.2f, %+.2f, %+.2f )\n", aa.nx, aa.ny, aa.nz );
        pr( logFile, "Final Rotation Angle = %5.1f deg\n", RadiansToDegrees(WrpRad( ModRad(aa.ang))) );
	}

        copyState( &sHist[ nconf+irun ], sSave );

        for (i=0; i<ntor; i++) sHist[ nconf+irun ].tor[i] = sSave.tor[i];

        if (ntor > 0 && outlev>=LOGRUNVV) {
	    pr( logFile, "Final Torsions:\n" );
            for (i=0; i<ntor; i++) {
                torTmp = RadiansToDegrees( sSave.tor[i] );
                torTmp = ModDeg( torTmp );
                torTmp = WrpDeg( torTmp );
                sHist[ nconf+irun ].tor[i] = sSave.tor[i];
                pr( logFile, "          %2d = %7.2f deg", i+1, torTmp);
                if (B_isTorConstrained[i] && B_ShowTorE) {
                    pr(logFile, ", Energetic penalty = %uhd\n", US_TorE[i]);
                } else {
                    pr(logFile, "\n");
                }
            }
        }

        cnv_state_to_coords( sSave, vt, tlist, ntor, crdpdb, crd, natom,
	 true_ligand_atoms, outlev, logFile);

        if (ntor > 0) {
            eintra = eintcal( nonbondlist, ptr_ad_energy_tables, crd,
	       Nnb, Nnb_array, group_energy,
               B_calcIntElec, B_include_1_4_interactions,
               scale_1_4, qsp_abs_charge, 
               B_use_non_bond_cutoff, B_have_flexible_residues,
	       outlev, logFile);
        } else {
            eintra = 0.0 ;
        }
        einter = scale_eintermol * trilinterp( 0, natom, crd, charge, abs_charge, type, map, 
                    info, ignore_inter, peratomE, NULL,
                    NULL_ENERGY_BREAKDOWN);

	if(outlev>=LOGMIN)
        writePDBQT( nconf+irun, runseed[nconf], FN_ligand, FN_dpf, lig_center, sSave, ntor,
                &eintra, &einter, natom, atomstuff, crd, peratomE,
                charge, abs_charge, qsp_abs_charge,
                ligand_is_inhibitor, torsFreeEnergy, 
                vt, tlist, crdpdb, nonbondlist, 
                ptr_ad_energy_tables,
                type, Nnb, Nnb_array, group_energy, true_ligand_atoms,
		B_calcIntElec,
                map,
                ignore_inter,
                B_include_1_4_interactions, scale_1_4, parameterArray, unbound_internal_FE,
                info, 1 /* = DOCKED */, PDBQT_record, 
                B_use_non_bond_cutoff, B_have_flexible_residues,
                ad4_unbound_model, outlev, logFile);

        // See also "calculateEnergies.cc", switch(ad4_unbound_model)
        if (ad4_unbound_model == Unbound_Same_As_Bound) {
            // Update the unbound internal energy, setting it to the current internal energy
            local_unbound_internal_FE = eintra;
        } else {
            local_unbound_internal_FE = unbound_internal_FE;
        }
        // originally: econf[nconf+irun] = eLast;
        econf[nconf+irun] = eintra + einter + torsFreeEnergy - local_unbound_internal_FE;

    } /* Loop over runs ======================================================*/
    *Addr_nconf += nruns;

    if ( B_writeTrj ) {
            fclose( FP_trj );
    }
}
/* EOF */
