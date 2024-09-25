/*

 $Id: cnv_state_to_coords.cc,v 1.13 2012/04/05 01:39:31 mp Exp $

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
#include <string.h>
#include <stdio.h>
#include "constants.h"
#include "structs.h"
#include "cnv_state_to_coords.h"
#include "stateLibrary.h"

#ifdef LSQFIT
#include "rotlsqfit.h"
#endif

void cnv_state_to_coords( const State& now,
                          const Real vt[MAX_TORS][SPACE],
                          const int tlist[MAX_TORS+1][MAX_ATOMS],
                          const int ntor,
                          const Real crdpdb[MAX_ATOMS][SPACE],
                          /* not const */ Real crd[MAX_ATOMS][SPACE],
                          const int natom,
			  const int true_ligand_atoms,
			  const int outlev,
			  FILE *logFile)

{
#ifdef DEBUG
    if (checkState(&now)) { // only if the state is valid can we proceed...
#endif

    //  Setting back to the original PDB
    //  coordinates ensures that cumulative
    //  rounding errors do not occur.

    (void) memcpy( crd, crdpdb,  natom * 3 * sizeof(Real));

    //  memcpy is about 100x faster than these nested for-loops...
    //  for (i = 0;  i < natom;  i++) for (XYZ = 0;  XYZ < SPACE;  XYZ++) crd[i][XYZ] = crdpdb[i][XYZ];

    //  Apply torsions, if any
    if (ntor > 0) {
        torsion( now, crd, vt, tlist, ntor ); /* all const except for crd */
    }

#ifdef LSQFIT
    //  Optionally, least-squares fit the current conformation on the original input conformation
    static double crdfit[MAX_ATOMS][SPACE];
    // double weights[MAX_ATOMS];
    // double vector[3];
    // double matrix[3][3];
    // rotlsqfit(crd, crdpdb, crdfit, weights, natom, matrix, vector)
    (void) rotlsqfit( crd, crdpdb, crdfit, NULL, natom, NULL, NULL);
    (void) memcpy( crd, crdfit,  natom * 3 * sizeof(Real));
#endif

#ifdef DEBUG_QUAT
#ifdef DEBUG_QUAT_PRINT
    pr( logFile, "DEBUG_QUAT: cnv_state_to_coords()\n" );
    (void) fflush(logFile);
#endif
    //  Make sure the quaternion is suitable for 3D rotation
    assertQuatOK( now.Q );
#endif

    //  Apply quaternion rigid-body rotation and translation...
    qtransform( now.T, now.Q, crd, true_ligand_atoms ); // all const except crd

#ifdef DEBUG
    } else {
        (void)fprintf(logFile,"WARNING: About to convert state to coordinates, but there is a problem with state! (The function 'checkState' failed).\n");
        // exit(-1);
    } // endif checkState
#endif
}

/* EOF */
