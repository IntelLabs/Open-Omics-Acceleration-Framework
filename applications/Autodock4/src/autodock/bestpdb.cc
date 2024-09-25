/*

 $Id: bestpdb.cc,v 1.11 2014/06/12 01:44:07 mp Exp $

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

/* bestpdb.cc */

#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "print_rem.h"
#include "strindex.h"
#include "print_avsfld.h"
#include "bestpdb.h"


extern int keepresnum;
extern char dock_param_fn[];

void bestpdb( const int ncluster, 
	      const int num_in_clu[MAX_RUNS],
	      const int cluster[MAX_RUNS][MAX_RUNS],
	      const Real econf[MAX_RUNS],
	      const Real crd[MAX_RUNS][MAX_ATOMS][SPACE],
	      const char atomstuff[MAX_ATOMS][MAX_CHARS],
	      const int natom,
	      const Boole B_write_all_clusmem,
	      const Real ref_rms[MAX_RUNS],
	      const int outlev,
	      FILE *const logFile)

{
    register int  i=0,
	          j=0,
	          k=0,
		  confnum=0;

    int           c = 0,
		  kmax = 0,
		  /* imol = 0, */
		  indpf = 0,
		  off[7],
		  nframes = 0,
	          stride = 0,
		  c1 = 1,
		  i1 = 1;
		  

    char  filnm[PATH_MAX],
		  label[MAX_CHARS];
    char AtmNamResNamNumInsCode[20]; /* PDB record 0-origin indices 11-29 (from blank after serial_number to just before xcrd */

    pr( logFile, "\n\tLOWEST ENERGY DOCKED CONFORMATION from EACH CLUSTER");
    pr( logFile, "\n\t___________________________________________________\n\n\n" );

    if (keepresnum > 0 ) {
	pr( logFile, "\nKeeping original residue number (specified in the input PDBQ file) for outputting.\n\n");
    } else {
	pr( logFile, "\nResidue number will be the conformation's rank.\n\n");
    }

    for (i = 0;  i < ncluster;  i++) {
        i1 = i + 1;

        if (B_write_all_clusmem) {
            kmax = num_in_clu[i];
        } else {
            kmax = 1;	/* write lowest-energy only */
        }

        for ( k = 0; k < kmax; k++ ) {
            c = cluster[i][k];
            c1 = c + 1;

            fprintf( logFile, "USER    DPF = %s\n", dock_param_fn);
            fprintf( logFile, "USER    Conformation Number = %d\n", ++confnum);

            print_rem(logFile, i1, num_in_clu[i], c1, ref_rms[c]);


            if (keepresnum > 0) {
                fprintf( logFile, "USER                              x       y       z   Rank Run  Energy    RMS\n");
                for (j = 0;  j < natom;  j++) {
                    sprintf(AtmNamResNamNumInsCode, "%-19.19s", &atomstuff[j][11]);
                    // retain original residue number (in fact, all fields
                    // from blank after atom serial number to start of coords)
                    // replace occupancy by cluster index, 
                    // tempfactor by conformation index within cluster,
                    // add two non-standard fields with energy and RMSD from reference
                    #define FORMAT_PDBQT_ATOM_RANKRUN_STR     "ATOM  %5d%-19.19s%8.3f%8.3f%8.3f%6d%6d    %+6.2f %8.3f\n"
                    fprintf(logFile, FORMAT_PDBQT_ATOM_RANKRUN_STR, j+1, AtmNamResNamNumInsCode, 
                           crd[c][j][X], crd[c][j][Y], crd[c][j][Z],  i1, c1, econf[c], ref_rms[c] );

                } /* j */
            } else {
                fprintf( logFile, "USER                   Rank       x       y       z    Run   Energy    RMS\n");
                for (j = 0;  j < natom;  j++) {
                    sprintf(AtmNamResNamNumInsCode, "%-11.11s%4d%-4.4s", &atomstuff[j][11], 
                           i1, &atomstuff[j][26]);
                    // replace original residue number by cluster index
                    // replace occupancy by conformation index within cluster
                    // tempfactor by energy
                    // add one non-standard field with RMSD from reference
                    #define FORMAT_PDBQT_ATOM_RANKRUN_NUM          "ATOM  %5d%-19.19s%8.3f%8.3f%8.3f%6d%+6.2f    %6.3f\n"
                    fprintf(logFile, FORMAT_PDBQT_ATOM_RANKRUN_NUM, j+1, AtmNamResNamNumInsCode, 
                           crd[c][j][X], crd[c][j][Y], crd[c][j][Z],  c1, econf[c], ref_rms[c] );
                } /* j */
            }
            fprintf( logFile, "TER\n" );
            fprintf( logFile, "ENDMDL\n" );
            fflush( logFile );

            nframes++;
        } /* for k */

    } /* for i */

    fprintf( logFile, "\n" );

    strcpy(label, "x y z Rank Run Energy RMS\0" );
    if (keepresnum > 0) {
        off[0]=5; off[1]=6; off[2]=7; off[3]=8; off[4]=9; off[5]=10; off[6]=11;
        stride=12;
    } else {
        off[0]=5; off[1]=6; off[2]=7; off[3]=4; off[4]=8; off[5]=9; off[6]=10;
        stride=11;
    } /* if */
     
    indpf = strindex( dock_param_fn, ".dpf" );
    strncpy( filnm, dock_param_fn, (size_t)indpf );
    filnm[ indpf ] = '\0';
    strcat( filnm, ".dlg.pdb\0" );
     
    print_avsfld( logFile, 7, natom, nframes, off, stride, label, filnm );
}
/* EOF */
