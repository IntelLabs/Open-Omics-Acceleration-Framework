/*

 $Id: clmode.cc,v 1.24 2014/06/12 05:04:57 mp Exp $

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
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/types.h>
#include "timesyshms.h"
#include "clmode.h"


extern char *programname;

void  clmode( const int   num_atm_maps,
              ConstReal   clus_rms_tol,
              const char  *const hostnm,
              const Clock& jobStart,
              const struct tms& tms_jobStart,
              const Boole write_all_clusmem,
              const char  *const clusFN,
              const Real crdpdb[MAX_ATOMS][SPACE],
              const Real sml_center[SPACE],
              const Boole B_symmetry_flag,
              const Boole B_unique_pair_flag,
              const char  *const rms_ref_crds,
              const Boole B_rms_heavy_atoms_only,
              const int h_index,
	      const int outlev,
	      FILE *logFile)

{
    FILE *clusFile;
    register int xyz = 0;
    int   anum = 0;
    char  atomstuff[MAX_ATOMS][MAX_CHARS];
    Real econf[MAX_RUNS];
    Real eSave[2];
    Boole haveAtoms = FALSE;
    Boole haveTypes = FALSE;
    int   ii = 0;
    int   lastanum = -1;
    char  line[LINE_LEN];
    int   atomCounter = 0;
    int   natom = 0;
    int   natom_1 = -1;
    int   nconf = 0;
    int   confCounter = 0;
    int   ntype[MAX_ATOMS];
    char  pdbaname[MAX_ATOMS][5];
    Real q = 0.;
    char  rec5[5];
    int   nsaved = 0;
    char  anumStr[5];
    int   type[MAX_ATOMS];
    // the next three arrays are "static" solely to prevent stack size
    // failures on some computers when MAX_RUNS is large (e.g. 2000)
    static Real crdSave[MAX_RUNS][MAX_ATOMS][SPACE];
    static Real clu_rms[MAX_RUNS][MAX_RUNS];
    static int   cluster[MAX_RUNS][MAX_RUNS];
    register int i = 0;
    register int j = 0;
    int   irunmax = -1;
    int   isort[MAX_RUNS];
    int   ncluster = 0;
    int   num_in_clu[MAX_RUNS];
    Real ref_crds[MAX_ATOMS][SPACE];
    int   ref_natoms = -1;
    Real ref_rms[MAX_RUNS];
    Boole haveEnergy = FALSE;
    ParameterEntry thisparm;

    for (j = 0; j < MAX_RUNS; j++) {
        num_in_clu[j] = 0;
        isort[j] = j;
        econf[j] = 0.;
    }
    /*
     * Open file containing coordinates to be clustered...
     */
    if ( openFile( clusFN , "r", &clusFile, jobStart, tms_jobStart, TRUE, logFile) ) {
        pr( logFile, "Conformations to be clustered are in this file: \"%s\"\n\n", clusFN );
    }
    /*
     * Read in the conformations
     * All we need are the xyz's of each conformation,
     * and their Energies, plus the Run number/parent dlg file.
     */
    while ( fgets( line, LINE_LEN, clusFile) != NULL ) {

        pr( logFile, "INPUT-PDBQT: %s", line);

        for (ii = 0; ii < 4; ii++) { rec5[ii] = tolower( (int)line[ii] ); };

        if (( strindex( line, "USER    Total Interaction Energy of Complex") >= 0 )
         || ( strindex( line, "REMARK  Total Interaction Energy of Complex") >= 0 )
         || ( strindex( line, "USER    Estimated Free Energy of Binding") >= 0 )  // Added for 4.2.x compatibility SF
        ){
            /*
             * Read in the energy of this conformation;
             * This is preferred over "Final Docked Energy" because this is never
             * printed out rounded up as +4.25e+03, but always as +4246.45, e.g.:
             */
            if ( haveAtoms ) {
                econf[confCounter] = 0.;
                sscanf( line, "%*s %*s %*s %*s %*s %*s %*s " FDFMT, &econf[confCounter]);
                haveEnergy = TRUE;
            } else {
                /* ! haveAtoms
                 * We have not seen any atoms yet, so save this energy. 
                 */
                eSave[nsaved]=0.;
                sscanf( line, "%*s %*s %*s %*s %*s %*s %*s " FDFMT, &eSave[nsaved] );
                ++nsaved;
            }

        } else if ( (( strindex( line, "USER    Final Docked Energy") >= 0 ) 
                  || ( strindex( line, "REMARK  Final Docked Energy") >= 0 )) && ( ! haveEnergy ) ) {
            /*
             * Read in the energy of this conformation if we don't already
             * have an energy:
             */
            if ( haveAtoms ) {
                econf[confCounter] = 0.;
                sscanf( line, "%*s %*s %*s %*s %*s " FDFMT, &econf[confCounter]);
                haveEnergy = TRUE;
            } else {
                /* ! haveAtoms
                 * We have not seen any atoms yet, so save this energy. 
                 */
                eSave[nsaved]=0.;
                sscanf( line, "%*s %*s %*s %*s %*s " FDFMT, &eSave[nsaved] );
                ++nsaved;
            }

        } else if (equal( rec5,"atom", 4) || equal( rec5,"heta", 4)) {

            int serial;

            /* 
             * This line should contain coordinates, partial charge & 
             * atom type for one atom.
             * Let's save the coordinates for this atom, atomCounter.
             */
            readPDBQTLine( line, &serial, crdSave[confCounter][atomCounter], &q, &thisparm, outlev, logFile );

            if ( ! haveAtoms ) {
                /*
                 * We do not have any atoms for this conformation,
                 */
                sscanf( &line[6], "%s", anumStr );

                if ((anum = atoi(anumStr)) < lastanum) { /* initially, lastanum is -1, while anum is probably never -1, so this is false initially */
                    /*
                     * haveAtoms is FALSE, but this line begins with "atom" or
                     * "heta", so this must be the...
                     *
                     * Start of next conformation,
                     */
                    /* This is an atom line, so haveAtoms must be set to true: */
                    haveAtoms = TRUE;
                    /* We must also have read in the atom types: */
                    haveTypes = TRUE;
                    /* Transfer the saved energies to the econf arrays: */
                    econf[0] = eSave[0];
                    econf[1] = eSave[1];
                    /* Now we have the energy: */
                    haveEnergy = TRUE;
                    /* Increment the number of conformations */
                    ++confCounter;
                    for (xyz = 0;  xyz < SPACE;  xyz++) {
                        crdSave[confCounter][0][xyz] = crdSave[confCounter-1][atomCounter][xyz]; 
                    }
                    natom = atomCounter; /* number of atoms is set to the atom counter*/
                    natom_1 = natom - 1;  /* number of atoms minus 1, for 0-based counting */
                    atomCounter = 0; /* reset the counter "atomCounter" */
                } else {
                    /* 
                     * First of all, determine the atom types for all the atoms in the
                     * molecule: 
                     */
                    strncpy( atomstuff[atomCounter], line, (size_t)30 );
                    atomstuff[atomCounter][30] = '\0';
                    if ( ! haveTypes ) {
                        //sscanf( &line[12], "%s", pdbaname[atomCounter] );
                        if(strlen(line)<77+1) {
                            // TODO MPique 2010-06 write test for this 77+1
                            pr( logFile, "\nNOTE: Atom number %d, line too short, cannot determine atom type\n", atomCounter+1);
                            strcpy(pdbaname[atomCounter], "?");
                        }
				
                        else sscanf( &line[77], "%s", pdbaname[atomCounter] ); //  Added for 4.2.x compatibility SF
                        /*
                         * Determine this atom's atom type:
                         */
                        type[atomCounter] = get_atom_type(pdbaname[atomCounter]);

                        if (type[atomCounter] == -1) {
                            pr( logFile, "\nNOTE: Atom number %d, unknown atom type \n\n", atomCounter+1);
                            stop("Unknown atom type");
                        }
                        /* 
                         * Increment the number of atoms with this atom type:
                         */
                        ++ntype[ type[atomCounter] ];
                    }
                }
                /*
                 * Update the value of the last atom's serial number:
                 */
                lastanum = anum;
            } else if ( atomCounter == natom_1 ) {  /* initially, atomCounter=0, and natom_1= -1, so this is not true initially */
                /*
                 * We have all the atoms we expect for one molecule:
                 * Increment total number of conformations,
                 */
                ++confCounter;
                atomCounter = -1; /*  Pre-zero out the "atomCounter" counter... */
                haveEnergy = FALSE; /* we don't have energy yet for next conf. */
            }
            /* 
             * Just increment the number of atoms, atomCounter:
             */
            ++atomCounter;
            
        } /* This was an "atom" or "heta" line */
    } /* end while there is a new line. */

    irunmax = confCounter;
    nconf = confCounter;

    pr( logFile, "\nNumber of conformations found = %d\n", nconf );

    pr( logFile, "\nNumber of atoms per conformation = %d\n\n", natom );

    for (i=0; i<num_atm_maps; i++) {
        pr( logFile, "Number of atoms with type %d = %d\n", i+1, ntype[i]);
    }

    if (strncmp(rms_ref_crds,"unspecified filename",20) != 0) {
        /*
         * Read in reference structure, specified by the "rmsref" command...
         */
        if ((ref_natoms = getpdbcrds( rms_ref_crds, ref_crds, logFile)) == -1) {
     
            fprintf( logFile, "%s: Problems while reading \"%s\".\n", programname, rms_ref_crds);
            fprintf( logFile, "Will attempt to use the input PDBQ file coordinates as reference.\n");
     
        } else if (ref_natoms != natom) {
     
            pr( logFile, "%s: ERROR!  Wrong number of atoms in reference structure.\n", programname);
            pr( logFile, "Input PDBQ structure has %d atoms, but reference structure has %d atoms.\n\n", natom, ref_natoms);
            ref_natoms = -1;
     
        }
    }

    if (nconf <= 1) {

        pr( logFile, "\nSorry!  Unable to perform cluster analysis, because not enough structures were read in.\n");

    } else {

#ifdef DEBUG
     for (i = 0;  i < nconf;  i++) { 
        pr( logFile, "i=%-3d\tisort[i]=%-3d\teconf[isort[i]]=%+7.2f\n", 
                      i, isort[i], econf[isort[i]] );
    }
#endif /* DEBUG */

        pr( logFile, "\nSorting %d conformations by their energy.\n", irunmax);
        flushLog;

        sort_enrg( econf, isort, nconf );

#ifdef DEBUG
     for (i = 0;  i < nconf;  i++) { pr( logFile, "i=%-3d\tisort[i]=%-3d\teconf[isort[i]]=%+7.2f\n", i, isort[i], econf[isort[i]] ); }
#endif /* DEBUG */

        pr( logFile, "\nPerforming cluster analysis, using a cluster RMS tolerance of %.1f\n", clus_rms_tol );
        flushLog;

        ncluster = cluster_analysis( clus_rms_tol, cluster, num_in_clu, isort, 
                                     nconf, natom, type, crdSave, crdpdb, 
                                     sml_center, clu_rms, 
                                     B_symmetry_flag, B_unique_pair_flag,
                                     ref_crds, ref_natoms, ref_rms,
                                     B_rms_heavy_atoms_only, h_index);

        pr( logFile, "\nOutputting structurally similar clusters, ranked in order of increasing energy.\n" );
        flushLog;

        prClusterHist( ncluster, irunmax, clus_rms_tol, num_in_clu, 
                       cluster, econf, clu_rms, ref_rms, outlev, logFile);

        bestpdb( ncluster, num_in_clu, cluster, econf, crdSave, 
                 atomstuff, natom, write_all_clusmem, ref_rms, outlev, logFile);

    }/*if we have more than 1 conformation... */

/*
 *  End cluster_mode and END PROGRAM...
 */
    success( hostnm, jobStart, tms_jobStart, logFile);

    exit(EXIT_SUCCESS); // POSIX, defined in stdlib.h
}
/* EOF */
