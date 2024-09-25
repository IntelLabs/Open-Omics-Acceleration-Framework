/*

 $Id: cluster_analysis.cc,v 1.10 2011/10/10 17:42:24 rhuey Exp $

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
#include "cluster_analysis.h"


int cluster_analysis( ConstReal   clus_rms_tol, 
                      /* not const */ int cluster[MAX_RUNS][MAX_RUNS], 
                      /* not const */ int num_in_clus[MAX_RUNS], 
                      const int isort[MAX_RUNS], 
                      const int nconf, 
                      const int natom, 
                      const int type[MAX_ATOMS],
                      const Real crd[MAX_RUNS][MAX_ATOMS][SPACE], 
                      const Real crdpdb[MAX_ATOMS][SPACE], 
                      const Real sml_center[SPACE], 
                      /* not const */ Real clu_rms[MAX_RUNS][MAX_RUNS], 
                      const Boole B_symmetry_flag,
                      const Boole B_unique_pair_flag,
                      /* not const */ Real ref_crds[MAX_ATOMS][SPACE],
                      const int ref_natoms,
                      /* not const */ Real ref_rms[MAX_RUNS],
                      const Boole B_rms_heavy_atoms_only,
                      const int h_index)
{
/* __________________________________________________________________________
  | Cluster Analysis                                                         |
  |__________________________________________________________________________|
  |  If conformations are within clus_rms_tol, scored as the same.           |
  |  Compares atoms with same type, not atom name, to find propellers.       |
  |__________________________________________________________________________|*/

    register int compare = 0,
                 i = 0;

    int          nClusters = 1,
          thisconf = 0,
          new_conf = FALSE;

    Real rms = 0.;

    if (ref_natoms == -1) {
        // No reference coordinates were defined, 
        //
        // Assume the Input Ligand PDBQT file is the reference structure
        // we must un-center the original PDBQT coordinates, 

        for (i = 0;  i < natom;  i++) {
                ref_crds[i][0] = sml_center[0] + crdpdb[i][0];
                ref_crds[i][1] = sml_center[1] + crdpdb[i][1];
                ref_crds[i][2] = sml_center[2] + crdpdb[i][2];
        }/*i*/
    }

/* Assign the index of the lowest energy to 0,0 in "cluster" */

    thisconf = cluster[0][0] = isort[0]; 

/* Set number in 0-th cluster to 1 */
/* Also initialize total number of clusters to 1 */

    num_in_clus[0] = nClusters = 1;         
    clu_rms[0][0]  = getrms(crd[thisconf], ref_crds,
                            B_symmetry_flag, B_unique_pair_flag, natom, type,
                            B_rms_heavy_atoms_only, h_index);

/* Go through *all* conformations... */
    for ( i=0; i<nconf; i++) {

/* Calculate the RMSD to the reference structure: */
        ref_rms[i] = getrms(crd[i], ref_crds, B_symmetry_flag, B_unique_pair_flag, natom, type, 
                            B_rms_heavy_atoms_only, h_index);
    }

/* Go through all conformations *except* 0-th... */
    for ( i=1; i<nconf; i++) {        

/* "thisconf" is the index of the energy-sorted i-th conformation: */
        thisconf = isort[i];

/* Assume this is a new conformation until proven otherwise... */
        new_conf = TRUE;

        for ( compare=0; compare<nClusters; compare++ ) {

            rms = getrms(crd[thisconf], crd[cluster[compare][0]],
                         B_symmetry_flag, B_unique_pair_flag, natom, type,
                         B_rms_heavy_atoms_only, h_index);

/* Check rms; if greater than tolerance, */
            if ( rms > clus_rms_tol ) {

                continue; /* to compare next conformation, */

            } else {

/* Otherwise, add this conformation to the current cluster, */
                cluster[compare][ num_in_clus[compare] ] = thisconf;
                clu_rms[compare][ num_in_clus[compare] ] = rms;
                num_in_clus[compare]++;
                new_conf = FALSE;
                break; /* Go on to next i iteration... */
            }
        } /* next comparison... */

        if (new_conf) {

/* Start a new cluster...  */
            cluster[ nClusters ][0] = thisconf;
            clu_rms[ nClusters ][0] = getrms(crd[thisconf], ref_crds,
                                             B_symmetry_flag, B_unique_pair_flag, natom, type,
                                             B_rms_heavy_atoms_only, h_index);
            num_in_clus[ nClusters ] = 1;

/* Increment the number of clusters... */
            nClusters++;

        } /* endif */

    } /* next i... */

    return nClusters;
}
/* EOF */
