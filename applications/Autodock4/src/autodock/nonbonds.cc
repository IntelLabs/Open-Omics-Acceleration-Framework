/*

 $Id: nonbonds.cc,v 1.18 2012/10/15 17:48:28 mp Exp $

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
#include "nonbonds.h"
#include "quicksort.h"
#include "mdist.h"  /* mindist and maxdist are here */


extern char *programname;

using namespace std;

// nonbonds - returns 0 if OK, else non-zero for error
int
nonbonds(const Real  crdpdb[MAX_ATOMS][SPACE],
		      /* not const */ int         nbmatrix[MAX_ATOMS][MAX_ATOMS],
	      const int   natom, 
              const int   bond_index[MAX_ATOMS],
              const int         B_include_1_4_interactions,
	      const int	nbonds[MAX_ATOMS],
              const int         bonded[MAX_ATOMS][MAX_NBONDS],
	      const int debug,
	      const int outlev,
	      FILE *logFile)
{
	int i,j,k,l;
    int nonbond_type;
    int errorcode = 0;

    if (debug>0) {
        printbonds(natom, nbonds, bonded, "\nDEBUG:  3. INSIDE nonbonds, bonded[][] array is:\n\n", 0, outlev, logFile);
    }

    //
    // in "nbmatrix", the values 1 (and 4) mean this pair of atoms will be included in the internal, non-bonded list
    //                           0                                         ignored
 
    // set all nonbonds in nbmatrix to 1, except "1-1 interactions" (self interaction)
    for (i = 0; i<MAX_ATOMS; i++) {
        for (j = 0; j<MAX_ATOMS; j++){
            nbmatrix[i][j] = 1;
        } // j
        nbmatrix[i][i] = 0; /* 2005-01-10 RH & GMM */
    } // i

	for (i = 0; i<natom; i++) {  // loop over atoms, "i", from 1 to natom
        for (j=0; j<nbonds[i]; j++) {  // loop over this atom's "j" bonds
            // Ignore 1-2 Interactions
            nbmatrix[ i            ][ bonded[i][j] ] = 0;
            nbmatrix[ bonded[i][j] ][ i            ] = 0;
        } // j
    } // i

	for (i=0; i<natom; i++) {  // loop over all the atoms, "i"; there are "natom" atoms
		for (j=0; j<nbonds[i]; j++) {  // loop over all the atoms "j" bonded to the current atom, "i"
			for (k=0; k<nbonds[bonded[i][j]]; k++) {  // loop over each atom "k" bonded to the current atom "j"

				// Ignore "1-3 Interactions"
				nbmatrix[bonded[bonded[i][j]][k]][i] = 0;
				nbmatrix[i][bonded[bonded[i][j]][k]] = 0;

                if (B_include_1_4_interactions == FALSE) {
                    // include is FALSE,  i.e.: ignore
                    // Ignore "1-4 Interactions"
                    nonbond_type = 0;
                } else {
                    // remember that this is a 1-4 interaction so we can scale it later...
                    // Remember "1-4 Interactions"
                    nonbond_type = 4;
                } //  end if we are ignoring "1-4 interactions"

                // Loop over each atom "l" bonded to the atom "k" bonded to the atom "j" bonded to the atom "i"...
                if (debug>0) {
                    (void)fprintf(logFile, "bonded[ %d ][ %d ] = %d\n", bonded[i][j], k, bonded[bonded[i][j]][k] );
                    (void)fprintf(logFile, "nbonds[ %d ] = %d\n", bonded[bonded[i][j]][k], nbonds[bonded[bonded[i][j]][k]] );
                }
                // validate the bonded[][] values before they used as indices in the for l-loop.
                if (bonded[i][j] < 0) {
                    (void)fprintf(logFile, "%s:  WARNING! The count of bonds, %d, to atom %d does not match the list of bonded atoms.\n\n", programname, nbonds[i], i);
                    (void)fprintf(logFile, "%s:  WARNING! bonded[%d][] = %d, %d, %d, %d, %d, nbonds=%d.\n\n", programname, i, bonded[i][0], bonded[i][1], bonded[i][2], bonded[i][3], bonded[i][4], nbonds[i]);
                    errorcode |= 1;
                }
                if (bonded[ bonded[i][j] ][k] < 0) {
                    (void)fprintf(logFile, "%s:  WARNING! The count of bonds, %d, to atom %d does not match the list of bonded atoms.\n\n", programname, nbonds[bonded[i][j]], bonded[i][j]);
                    (void)fprintf(logFile, "%s:  WARNING! bonded[%d][] = %d, %d, %d, %d, %d, nbonds=%d.\n\n", programname, bonded[i][j], bonded[bonded[i][j]][0], bonded[bonded[i][j]][1], bonded[bonded[i][j]][2], bonded[bonded[i][j]][3], bonded[bonded[i][j]][4], nbonds[bonded[i][j]]); // gadzooks
                    errorcode |= 2;
                }
                for (l=0; l<nbonds[bonded[bonded[i][j]][k]]; l++) { 
                    nbmatrix[i][bonded[bonded[bonded[i][j]][k]][l]] = nonbond_type;
                    nbmatrix[bonded[bonded[bonded[i][j]][k]][l]][i] = nonbond_type;
                } //  l

			} //  k
		} //  j
	} //  i
	return errorcode;
} // end of nonbonds

/*----------------------------------------------------------------------------*/

// getbonds - returns 0 if OK, else non-zero for error
int
getbonds(const Real crdpdb[MAX_ATOMS][SPACE],
              const int from_atom,
              const int to_atom,
              const int bond_index[MAX_ATOMS],
	      /* not const */ int nbonds[MAX_ATOMS],
              /* not const */ int bonded[MAX_ATOMS][MAX_NBONDS],
	      const int debug,
	      const int outlev,
	      FILE *logFile
	      )
{
	int i,j;
	double dist,dx,dy,dz;
    Real distance[MAX_ATOMS];
    int isort[MAX_ATOMS];
    int errorcode = 0;

    // set up all the minimum and maximum possible distances for bonds
	mdist();

    // determine the bonded atoms or "1-2 interactions"

    // loop over atoms, "i", from "from_atom" to "to_atom"
	for (i = from_atom;  i < to_atom;  i++) {

        // compute the distance to all the atoms
        for (j = i+1; j < to_atom; j++) {
            isort[j] = j;
        }
        for (j = i+1; j < to_atom; j++) {
			dx = crdpdb[i][X] - crdpdb[j][X];
			dy = crdpdb[i][Y] - crdpdb[j][Y];
			dz = crdpdb[i][Z] - crdpdb[j][Z];
			distance[j] = sqrt(dx*dx + dy*dy + dz*dz);  // calculate the distance from "i" to "j"
        }
        // sort by distance
        quicksort( distance, isort, i+1, to_atom-1 );
        if (debug>0) {
            // print out distances
            for (j = i+1; j < to_atom; j++) {
                (void)fprintf(logFile, "i=%d, j=%d, distance[j]=%.3lf, distance[%2d]=%.3lf, isort[j]=%d\n", i, j, distance[j], isort[j], distance[isort[j]], isort[j]);
            }
        }

        // loop over the remaining atoms, from closest to furthest...
        // while on atom "i", loop over atoms "j_index", from "i+1" to "to_atom"
		for (int j_index = i+1;  j_index < to_atom;  j_index++) {
	
            // the isort array points to the indices sorted by distance
            j = isort[ j_index ];

            // find the shortest distance
			dist = distance[ j ];

            if (debug>0) {
                // print out distances
                (void)fprintf(logFile,"j_index=%d, j=%d, dist=%.3lf\n", j_index, j, dist);
            }

            // if distance from "i" to "j" is in range for their atom types, 
            // set one of the atoms to which "i" is bonded to "j", and vice-versa.
			if (dist >= mindist[bond_index[i]][bond_index[j]] && 
                dist <= maxdist[bond_index[i]][bond_index[j]]) {  

	            // Remember:   int bonded[MAX_ATOMS][MAX_NBONDS];	

                    if (nbonds[i] >= MAX_NBONDS || nbonds[j] >= MAX_NBONDS) {
                    if (nbonds[i] >= MAX_NBONDS) {
                        // bonded array for the i-th atom is full; return failure code here.
                        fprintf( logFile, "%s: WARNING!  Atom %d has too many bonds %d (%d is the maximum): ", programname, i+1, nbonds[i], MAX_NBONDS );
                        for(int k=0;k<MAX_NBONDS;k++) fprintf( logFile, "%d ", bonded[i][k]+1);
                        fprintf( logFile, "[%d...]\n", j+1);
			errorcode |= 1;
			}
                    if (nbonds[j] >= MAX_NBONDS) {
                        // bonded array for the j-th atom is full; return failure code here.
                        fprintf( logFile, "%s: WARNING!  Atom %d has too many bonds %d (%d is the maximum): ", programname, j+1, nbonds[j], MAX_NBONDS );
                        for(int k=0;k<MAX_NBONDS;k++) fprintf( logFile, "%d ", bonded[j][k]+1);
                        fprintf( logFile, "[%d...]\n", i+1);
			errorcode |= 2;
			}

		    } else {
                // add a bond between i and j
				bonded[i][ nbonds[i] ] = j;
				bonded[j][ nbonds[j] ] = i;
                // increment the number of bonds to i and j
				nbonds[i] += 1;
				nbonds[j] += 1;
		   }

                if (debug>0) {
                    // print out distances
                    (void)fprintf(logFile,"Adding a bond between %d and %d (dist=%.3lf, mindist=%.3lf, maxdist=%.3lf)\n", i, j, 
                                  dist, mindist[bond_index[i]][bond_index[j]], maxdist[bond_index[i]][bond_index[j]]);
                    (void)fprintf(logFile,"  bonded[%d][ nbonds[%d]=%d ] = %d\n", i, i, nbonds[i], j);
                    (void)fprintf(logFile,"  bonded[%d][ nbonds[%d]=%d ] = %d\n", j, j, nbonds[j], i);
                    (void)fprintf(logFile,"  nbonds[%d]=%d\n", i, nbonds[i]);
                    (void)fprintf(logFile,"  nbonds[%d]=%d\n", j, nbonds[j]);
                }
				
			} //  dist is in bonding range for these atom types
		} //  j
 	} //  i

	return errorcode;
} // end of get bonds

/*----------------------------------------------------------------------------*/

void printbonds(const int natom, const int nbonds[MAX_ATOMS], const int bonded[MAX_ATOMS][MAX_NBONDS], const char *message, const int B_print_all_bonds,
const int outlev, FILE *logFile)
{
    pr(logFile, "%s", message);
    for (int i = 0; i<natom; i++) {  // loop over atoms, "i", from 1 to natom
        pr(logFile, "DEBUG:  atom %d  bonded to ", i+1);
        for (int j=0; j< 
	    B_print_all_bonds 
	    ? MAX_NBONDS // loop over all the slots for this atom's "j" bonds
	    : nbonds[i]; // loop over this atom's "j" bonds
	    j++) pr(logFile, "  %d", bonded[i][j]+1);
        pr(logFile, "\n");
    } // i
} // end of printbonds

/*----------------------------------------------------------------------------*/

void print_1_4_message(const Boole B_include_1_4_interactions, const Real scale_1_4, const int outlev, FILE *logFile)
{
    if (B_include_1_4_interactions == FALSE) {
	if(outlev >= LOGBASIC)
        pr(logFile, "1,4-interactions will be _ignored_ in the non-bonded internal energy calculation.\n\n");
    } else if( outlev >= LOGMIN) {
        pr(logFile, "1,4-interactions will be _included_ in the non-bonded internal energy calculation.\n\n");
        pr(logFile, "1,4-interaction energies will be will be scaled by a factor of %.2lf .\n\n", (double)scale_1_4);
        pr(logFile, "NOTE:  Computed internal energies will differ from the standard AutoDock free energy function.\n\n");
    }
    (void) fflush(logFile);
} // end of print_1_4_message

/* EOF */
