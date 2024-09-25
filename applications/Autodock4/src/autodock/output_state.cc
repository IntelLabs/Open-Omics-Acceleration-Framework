/*

 $Id: output_state.cc,v 1.14 2012/04/17 23:08:42 mp Exp $

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
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "structs.h"
#include "output_state.h"

/* note: removed the not-so-portable file locking in favor of this read-noread method.
 *  I do not believe this facility is used widely.  M Pique, spring 2012
 */
// WIN32 MinGW seems to lack fchmod, S_IRGRP, and S_IROTH
#ifdef HAVE_FCHMOD
#define PERM_BUSY 0        /* watch-file when not to be readable by monitoring process */
#define PERM_DONE  S_IRUSR|S_IRGRP|S_IROTH        /* watch-file when readable by monitoring process */
#else
#define PERM_BUSY 0444    /* dont even try to fiddle with the permissions */
#endif

/*----------------------------------------------------------------------------*/
void output_state( FILE *const fp,
		   const State& S,
                   const int ntor,
                   const int istep,
                   ConstReal energy,
                   ConstReal eint,
                   const char lastmove,
                   const Boole B_watch,
                   const char *const FN_watch,
                   const char atomstuff[MAX_ATOMS][MAX_CHARS],
                   const int natom,
                   const Real crd[MAX_ATOMS][SPACE])
/*----------------------------------------------------------------------------*/
{
 // Writes state to fp (logFile) and optionally coords to "watch file" as pseudo-PDB.
 // This appears to be used only by simanneal.  M Pique spring 2012
    int i;
	
    int FD_watch;
    FILE *FP_watch;

    // note: state is printed as quaternion not axis-angle
    fprintf(fp, "state %d %c %f %f  %lf %lf %lf  %lf %lf %lf %lf\n",
        istep, lastmove, energy, eint, S.T.x, S.T.y, S.T.z,
        S.Q.x, S.Q.y, S.Q.z, S.Q.w );

    for (i=0; i<ntor; i++) {
        fprintf(fp, "%f\n", RadiansToDegrees( S.tor[i]) );
    }


    if (B_watch) {
        if ((FD_watch = creat( FN_watch, PERM_BUSY )) != -1) {;
            /* creates new file, or re-write old one */

            if ((FP_watch = fdopen( FD_watch, "w")) != NULL ) {
#ifdef HAVE_FCHMOD
		fchmod( FD_watch, PERM_BUSY);
#endif

                for (i = 0;  i < natom;  i++) {
		    fprintf( FP_watch, "%30s%8.3f%8.3f%8.3f\n",
		    atomstuff[i], crd[i][X], crd[i][Y], crd[i][Z]);
                }
                fflush( FP_watch );
#ifdef HAVE_FCHMOD
		fchmod( FD_watch, PERM_DONE);
#endif
                fclose( FP_watch );
            }
        }
    }

    return;
}
/* EOF */
