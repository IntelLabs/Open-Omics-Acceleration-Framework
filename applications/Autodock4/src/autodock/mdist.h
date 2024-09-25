/*

 $Id: mdist.h,v 1.12 2014/06/12 01:44:07 mp Exp $

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

#include "autocomm.h"

// expand the allowed bond length ranges by the BOND_LENGTH_TOLERANCE
#define BOND_LENGTH_TOLERANCE 0.1
#define set_minmax( a1, a2, min, max)  \
    mindist[(a1)][(a2)] = (mindist[(a2)][(a1)] = (min)-BOND_LENGTH_TOLERANCE);\
    maxdist[(a1)][(a2)] = (maxdist[(a2)][(a1)] = (max)+BOND_LENGTH_TOLERANCE)

void mdist();

enum {C=0,N=1,O=2,H=3,XX=4,P=5,S=6};  // see "bond_index" in the "AD4.1_bound.dat" or "AD4_parameters.dat" file.
#define NUM_ENUM_ATOMTYPES 7 // this should be the length of the enumerated atom types above
	
double mindist[NUM_ENUM_ATOMTYPES][NUM_ENUM_ATOMTYPES];
double maxdist[NUM_ENUM_ATOMTYPES][NUM_ENUM_ATOMTYPES];

void mdist() {

	register int i,j;

    // Set all the mindist and maxdist elements to the defaults for AutoDock versions 1 - 3...
	for (i=0; i<   NUM_ENUM_ATOMTYPES; i++) {
		for (j=0; j<   NUM_ENUM_ATOMTYPES; j++) {
			mindist[i][j] = 0.9 - BOND_LENGTH_TOLERANCE;
			maxdist[i][j] = 2.1 + BOND_LENGTH_TOLERANCE;
		}
	}

    /*
     * These values were contributed by Peter Reilly et al.:
     *
    mindist[C][C] = 1.32;
    maxdist[C][C] = 1.545;
    mindist[C][N] = 1.32;
    maxdist[C][N] = 1.39;
    mindist[C][O] = 1.20;
    maxdist[C][O] = 1.43;
    mindist[C][XX] = 1.07;
    maxdist[C][XX] = 1.15;
    mindist[C][S] = 1.80;
    maxdist[C][S] = 1.84;
    mindist[N][H] = 0.99;
    maxdist[N][H] = 1.03;
    mindist[N][XX] = 0.99;
    maxdist[N][XX] = 1.03;
    mindist[O][H] = 0.94;
    maxdist[O][H] = 0.98;
    mindist[O][P] = 1.47;
    maxdist[O][P] = 1.63;
    mindist[H][S] = 1.316;
    maxdist[H][S] = 1.356;
    mindist[XX][S] = 1.316;
    maxdist[XX][S] = 1.356;
    mindist[S][S] = 2.018;
    maxdist[S][S] = 2.058;
     */

    /*
     * These values, unless otherwise stated,
     * are taken from "handbook of Chemistry and Physics"
     * 44th edition(!)
     */
    set_minmax(C, C, 1.20, 1.545); // mindist[C][C] = 1.20, p. 3510 ; maxdist[C][C] = 1.545, p. 3511
    set_minmax(C, N, 1.1, 1.479); // mindist[C][N] = 1.1, p. 3510 ; maxdist[C][N] = 1.479, p. 3511
    set_minmax(C, O, 1.15, 1.47); // mindist[C][O] = 1.15, p. 3510 ; maxdist[C][O] = 1.47, p. 3512
    set_minmax(C, H, 1.022, 1.12);  // p. 3518, p. 3517
    set_minmax(C, XX, 0.9, 1.545); // mindist[C][XX] = 0.9, AutoDock 3 defaults ; maxdist[C][XX] = 1.545, p. 3511
    set_minmax(C, P, 1.85, 1.89); // mindist[C][P] = 1.85, p. 3510 ; maxdist[C][P] = 1.89, p. 3510
    set_minmax(C, S, 1.55, 1.835); // mindist[C][S] = 1.55, p. 3510 ; maxdist[C][S] = 1.835, p. 3512
    set_minmax(N, N, 1.0974, 1.128); // mindist[N][N] = 1.0974, p. 3513 ; maxdist[N][N] = 1.128, p. 3515
    set_minmax(N, O, 1.0619, 1.25); // mindist[N][O] = 1.0975, p. 3515 ; maxdist[N][O] = 1.128, p. 3515
    set_minmax(N, H, 1.004, 1.041); // mindist[N][H] = 1.004, p. 3516 ; maxdist[N][H] = 1.041, p. 3515
    set_minmax(N, XX, 0.9, 1.041); // mindist[N][XX] = 0.9, AutoDock 3 defaults ; maxdist[N][XX] = 1.041, p. 3515
    set_minmax(N, P, 1.4910, 1.4910); // mindist[N][P] = 1.4910, p. 3515 ; maxdist[N][P] = 1.4910, p. 3515
    set_minmax(N, S, 1.58, 1.672); // mindist[N][S] = 1.58, 1czm.pdb sulfonamide ; maxdist[N][S] = 1.672, J. Chem. SOC., Dalton Trans., 1996, Pages 4063-4069 
    set_minmax(O, O, 1.208, 1.51); // p.3513, p.3515
    set_minmax(O, H, 0.955, 1.0289); // mindist[O][H] = 0.955, p. 3515 ; maxdist[O][H] = 1.0289, p. 3515
    set_minmax(O, XX, 0.955, 2.1); // AutoDock 3 defaults
    set_minmax(O, P, 1.36, 1.67); // mindist[O][P] = 1.36, p. 3516 ; maxdist[O][P] = 1.67, p. 3517
    set_minmax(O, S, 1.41, 1.47); // p. 3517, p. 3515
    set_minmax(H, H, 100.,-100.); // impossible values to prevent such bonds from forming.
    set_minmax(H, XX, 0.9, 1.5); // AutoDock 4 defaults
    set_minmax(H, P, 1.40, 1.44); // mindist[H][P] = 1.40, p. 3515 ; maxdist[H][P] = 1.44, p. 3515
    set_minmax(H, S, 1.325, 1.3455); // mindist[H][S] = 1.325, p. 3518 ; maxdist[H][S] = 1.3455, p. 3516
    set_minmax(XX, XX, 0.9, 2.1); // AutoDock 3 defaults
    set_minmax(XX, P, 0.9, 2.1); // AutoDock 3 defaults
    set_minmax(XX, S, 1.325, 2.1); // mindist[XX][S] = 1.325, p. 3518 ; maxdist[XX][S] = 2.1, AutoDock 3 defaults
    set_minmax(P, P, 2.18, 2.23); // mindist[P][P] = 2.18, p. 3513 ; maxdist[P][P] = 2.23, p. 3513
    set_minmax(P, S, 1.83, 1.88); // mindist[P][S] = 1.83, p. 3516 ; maxdist[P][S] = 1.88, p. 3515
    set_minmax(S, S, 2.03, 2.05); // mindist[S][S] = 2.03, p. 3515 ; maxdist[S][S] = 2.05, p. 3515
    /* end values from Handbook of Chemistry and Physics */
}
