/*

 $Id: getrms.cc,v 1.10 2011/10/10 17:42:24 rhuey Exp $

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
#include "getrms.h"
#include <stdio.h>


Real getrms ( const Real Crd[MAX_ATOMS][SPACE], 
               const Real CrdRef[MAX_ATOMS][SPACE], 
               const Boole B_symmetry_flag, 
               const Boole B_unique_pair_flag,
               const int natom, 
               const int type[MAX_ATOMS],
               const Boole B_rms_heavy_atoms_only,
               const int h_index)

{
    double sqrSum, sqrMin, dc[SPACE];
    Boole is_available[MAX_ATOMS]; // available to be chosen as atom i's symmetry mate 
    register int i, j, xyz;
    int nhydrogen = 0;
    int num_atoms = natom;

    sqrSum = 0.;

    if (B_symmetry_flag) {
        for (i = 0;  i < natom;  i++) {
            is_available[i] = TRUE;
            if (type[i]==h_index) {
              nhydrogen++;
            }
        }
        for (i = 0;  i < natom;  i++) {
            if (B_rms_heavy_atoms_only&&type[i]==h_index) {
                continue; 
            } else {//skip this one
            sqrMin = BIG;
            int nearest_j_to_i=0;
            for (j = 0;  j < natom;  j++) {                
                if (type[i] == type[j] && is_available[j] ) { //include this atom if (1) h_index not set so there ARE no hydrogens; if (2) rms is to use all atoms even hydrogensx OR (3) this ith atom is not a hydrogen  
                    for (xyz = 0;  xyz < SPACE;  xyz++) {
                        dc[xyz]= Crd[i][xyz] - CrdRef[j][xyz];
                    } /* xyz */
                    double dist2= sqhypotenuse(dc[X], dc[Y], dc[Z]);
                    if(dist2<sqrMin) {
                        sqrMin=dist2;
                        nearest_j_to_i = j;
                    }
                }
            } /*  next j  */
            if(B_unique_pair_flag) is_available[nearest_j_to_i] = FALSE;
            sqrSum += sqrMin;
            };//include in rms
        } /*  next i  */
    } else {
        for (i = 0;  i < natom;  i++) {
            if (B_rms_heavy_atoms_only&&type[i]==h_index) continue; 
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                dc[xyz]= Crd[i][xyz] - CrdRef[i][xyz];
            } /* xyz */
            sqrSum += sqhypotenuse( dc[X], dc[Y], dc[Z] );
        } /*  next i  */
    } //else
    if (h_index >=0 && B_rms_heavy_atoms_only) num_atoms = natom - nhydrogen;
    return ( sqrt( sqrSum / (double)num_atoms )  );
}
/* EOF */
