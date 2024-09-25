/*

 $Id: initautodock.cc,v 1.18 2012/04/05 01:39:32 mp Exp $

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
#include "initautodock.h"

#define TINYDELTA 0.001         /* To nudge ligand into grid... */
#define AddNewHardCon(iCon,low,upp)  F_TorConRange[i][iCon][LOWER]=low;F_TorConRange[i][iCon][UPPER]=upp

extern int   keepresnum;
extern FILE *logFile;
extern char  *programname;

void initautodock(  const char  atomstuff[MAX_ATOMS][MAX_CHARS],
                    /* not const */ Real crd[MAX_ATOMS][SPACE],
                    const Real crdpdb[MAX_ATOMS][SPACE],
                    const int   natom,
                    const int   ntor,
                    /* not const */ State *const s0,
                    const int   tlist[MAX_TORS+1][MAX_ATOMS],
                    const Real vt[MAX_TORS][SPACE],
		    const int true_ligand_atoms,
                    const int   outlev,
                    const GridMapSetInfo *const info )

{
    Boole B_change = FALSE;
    Boole B_eq_and_opp = TRUE;
    Boole B_move_inside = FALSE;
    Boole B_outside;

    char  note[LINE_LEN];
    char  rec8[10];
    static char  axis[] = {"xyz"};

    Real delta[MAX_ATOMS][SPACE];
    Real delta_max[SPACE];
    Real delta_min[SPACE];
    Real last_delta[SPACE];

    int   ip[SPACE];
    int   ip_max[SPACE];
    int   ip_min[SPACE];

    register int   xyz = 0;
    register int   i = 0;

/*
**  Initialize the automated docking simulation,
*/
    /* Initialize the delta arrays... */
    for (xyz = 0;  xyz < SPACE;  xyz++) {
        delta_min[xyz] =  BIG;
        delta_max[xyz] = -BIG;
        ip[xyz] = ip_min[xyz] = ip_max[xyz] = 0;
    }
    for (i = 0;  i < natom;  i++) {
        for (xyz = 0;  xyz < SPACE;  xyz++) {
            delta[i][xyz] = 0.;
        }
    }

    if (outlev > LOGLIGREAD) {
        pr(logFile, "Allowable atom-coordinates are within these grid extents:\n\n");
        pr(logFile, "\t%7.3f < x <%7.3f\n\t%7.3f < y <%7.3f\n\t%7.3f < z <%7.3f\n\n", (double)info->lo[X], (double)info->hi[X], (double)info->lo[Y], (double)info->hi[Y], (double)info->lo[Z], (double)info->hi[Z]);
    }

    do {
/*
** Re-position the Small Molecule *until* it is inside grid...
*/
        /* Assume inside, until proved otherwise */
        B_outside = FALSE;

        /* if inside, there's no need to move inside. */
        B_move_inside = FALSE;        

        for (xyz = 0;  xyz < SPACE;  xyz++) {
            last_delta[xyz] =  delta[ip[xyz]][xyz];
        }

        if (outlev >= LOGLIGREAD  ) {
            pr( logFile, "Initializing the correction vectors and pointers.\n");
        }

        /* Initialize the delta arrays... */

        for (xyz = 0;  xyz < SPACE;  xyz++) {
            delta_min[xyz] =  BIG;
            delta_max[xyz] = -BIG;
            ip[xyz] = ip_min[xyz] = ip_max[xyz] = 0;
        }
        for (i = 0;  i < natom;  i++) {
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                delta[i][xyz] = 0.;
            }
        }

        if (outlev >= LOGLIGREAD ) {
            mkUnitQuat( &(s0->Q) );  // assure is normalized
	    AxisAngle aa = QuatToAxisAngle(s0->Q);
            pr( logFile, "Undoing all previous ligand transformations.\n");
            pr( logFile, "Resetting ligand to input PDBQ coordinates centred on \"about\" coordinates.\n" );
            if (ntor > 0) {
                pr( logFile, "Applying initial (relative) torsions... (\"0.0\" implies unchanged)\n" );
                pr( logFile,"\ndihe0 ");
                for (i=0; i<ntor; i++) {
                    pr( logFile," %+.1f", s0->tor[i]);
                }
            }
            pr( logFile, "\nApplying initial translation...\n");
            pr( logFile, "Ligand translated to:  %+.3f  %+.3f  %+.3f\n\n", s0->T.x, s0->T.y, s0->T.z );
            pr( logFile, "Applying initial quaternion...\n");
            pr( logFile, "Ligand rigid-body-rotated by: %+.1f degrees,   about unit vector: %+.3f %+.3f %+.3f.\n\n", RadiansToDegrees(aa.ang), aa.nx, aa.ny, aa.nz);
            flushLog;
        }

        cnv_state_to_coords( *s0,  vt, tlist, ntor,  crdpdb, crd, natom,
	 true_ligand_atoms, outlev, logFile); // all const except crd

        for (i = 0;  i < natom;  i++) {
            B_outside = is_out_grid_info( crd[i][X], crd[i][Y], crd[i][Z] );
            if ( B_outside ) {

                strncpy( rec8, &atomstuff[i][13], (size_t)8);
                rec8[8]='\0';
                if (outlev >= LOGLIGREAD) {
                    pr( logFile, "WARNING: Atom %s is outside grid!\n", rec8);
                }
/*
** Remember to move ligand inside grid...
*/
                B_move_inside = TRUE;
/*
** Figure out the deltas needed to move back into grid,
*/
                for (xyz = 0;  xyz < SPACE;  xyz++) {
                    if ( crd[i][xyz] < info->lo[xyz] )  {
                        delta[i][xyz] = info->lo[xyz] - crd[i][xyz] + TINYDELTA;
                        if (outlev >= LOGLIGREAD ) {
                            pr( logFile,"%s: atom %d, %c=%.3f, is too low, since grid-min=%.3f; ", programname, i+1, axis[xyz], crd[i][xyz], info->lo[xyz] );
                            pr( logFile,"increase %c by at least %.3f A\n", axis[xyz], (double)delta[i][xyz] );
                        }
                    } else if ( crd[i][xyz] ==  info->lo[xyz] )  {
                        if (outlev >= LOGLIGREAD ) {
                            pr( logFile,"%s: \"is_out_grid\"//lo macro failure on atom %d, %c-axis. Overriding move_inside instruction.\n", programname,i+1,axis[xyz]);
                        }
                        B_move_inside = FALSE;
                    }
                    if ( crd[i][xyz] > info->hi[xyz] )  {
                        delta[i][xyz] = info->hi[xyz] - crd[i][xyz] - TINYDELTA;
                        if (outlev >= LOGLIGREAD ) {
                            pr( logFile,"%s: atom %d, %c=%.3f, is too high, since grid-max=%.3f;", programname, i+1, axis[xyz], crd[i][xyz], info->hi[xyz] );
                            pr( logFile,"decrease %c by at least %.3f A\n", axis[xyz], -(double)delta[i][xyz] );
                        }
                    } else if ( crd[i][xyz] ==  info->hi[xyz] )  {
                        if (outlev >= LOGLIGREAD )  {
                            pr( logFile,"%s: \"is_out_grid\"//hi macro failure on atom %d, %c-axis. Overriding move_inside instruction.\n", programname,i+1,axis[xyz]);
                        }
                        B_move_inside = FALSE;
                    }
                }/*xyz*/
                if (outlev >= LOGLIGREAD ) {
                    pr( logFile,"\n" );
                }
            }/*if B_outside*/
        }/*i*/
        flushLog;

        if ( B_move_inside ) {
            if (outlev >= LOGLIGREAD ) {
                pr(logFile,"Axis Atom delta   delta_max ip_max delta_min ip_min axis\n");
                pr(logFile,"____ ____ _______ _________ ______ _________ ______ ____\n");
            }
            for (i = 0;  i < natom;  i++) {
                for (xyz = 0;  xyz < SPACE;  xyz++) {
                    B_change = FALSE;
                    if (delta[i][xyz] > delta_max[xyz]) {
                        delta_max[xyz] = delta[i][xyz]; 
                        ip_max[xyz] = i;
                        B_change = TRUE;
                    }
                    if (delta[i][xyz] < delta_min[xyz]) {
                        delta_min[xyz] = delta[i][xyz]; 
                        ip_min[xyz] = i;
                        B_change = TRUE;
                    }
                    if ((i>0) && B_change && (outlev >= LOGLIGREAD )) {
                        pr(logFile," %c   %3d %7.3f %7.3f    %3d   %7.3f    %3d\n", axis[xyz], (i+1), delta[i][xyz], delta_max[xyz], ip_max[xyz], delta_min[xyz], ip_min[xyz]);
                    }
                } /*xyz*/
            } /*i*/

            for (xyz = 0;  xyz < SPACE;  xyz++) {
                if ( fabs((double)delta_min[xyz]) > fabs((double)delta_max[xyz]) ) {
                    ip[xyz] = ip_min[xyz];
                } else {
                    ip[xyz] = ip_max[xyz];
                }
            } /*xyz*/

            if (outlev >= LOGLIGREAD ) {
                pr( logFile, "\n%s:\nLast tran0 correction vector was:\t(%.3f, %.3f, %.3f) \n", programname, (double)last_delta[X], (double)last_delta[Y], (double)last_delta[Z] );
                pr( logFile, "tran0 correction vector is now:\t(%.3f, %.3f, %.3f) \n", (double)delta[ip[X]][X], (double)delta[ip[Y]][Y], (double)delta[ip[Z]][Z] );
            }

            /*Assume*/ B_eq_and_opp = TRUE;
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                B_eq_and_opp = (delta[ip[xyz]][xyz] == -last_delta[xyz]) 
                               && B_eq_and_opp;
            }
            if (B_eq_and_opp) {

                prStr( note, "\n>>> Ligand does not fit within grid, while in this orientation!\n");
                pr_2x( stderr, logFile, note );
                prStr( note, ">>> Trying a new, randomly-generated rigid body rotation. (overriding user-specified initial rotation values)\n");
                pr_2x( stderr, logFile, note );

                s0->Q = randomQuat();

            }/*if (B_eq_and_opp)*/

            if (outlev > 1) {
                pr( logFile, "\n%s:\n*** Changing initial-translation (tran0 override):\n", programname );
                pr( logFile, "from:\ttran0  %.3f %.3f %.3f \n", s0->T.x, s0->T.y, s0->T.z );
            }

            s0->T.x +=  delta[ip[X]][X];
            s0->T.y +=  delta[ip[Y]][Y];
            s0->T.z +=  delta[ip[Z]][Z];

            if (outlev > 1) {
                pr( logFile, "to:\ttran0  %.3f %.3f %.3f \n\n", s0->T.x, s0->T.y, s0->T.z );
            }

        }/*if B_move_inside*/

    } while ( B_move_inside );

    flushLog;
}
/* EOF */
