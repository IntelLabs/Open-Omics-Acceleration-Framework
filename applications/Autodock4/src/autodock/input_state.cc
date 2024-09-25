/*

 $Id: input_state.cc,v 1.9 2011/06/03 05:31:36 mp Exp $

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
#include <string.h>
#include "input_state.h"

#define LINELEN 1000

/* TODO MPique 2011-06 appears unused in AutoDock or AutoGrid, untested 
 *   after adapting to read state as quaternion x,y,z,w rather than axis-angle
 */
int input_state( State *const S,
		 FILE  *const fp,
		 const char  *line,
		 const int   ntor,
		 int   *const p_istep,
		 Real  *const p_energy,
		 Real  *const p_eint,
		 char  *const p_lastmove )
{
    int i, istep, status;
    Real energy, eint;
    char lastmove;
    char myline[LINELEN];

#ifdef DEBUG
    fprintf(stderr, "line=|%s|\n", line);
#endif /* DEBUG */

    status = sscanf(line, "%*s %d %1s " FDFMT " " FDFMT " %lf %lf %lf %lf %lf %lf %lf", &istep, &lastmove, &energy, &eint,  &(S->T.x), &(S->T.y), &(S->T.z),  &(S->Q.x), &(S->Q.y), &(S->Q.z),  &(S->Q.w) );

    if (status != 0) {
	mkUnitQuat( &(S->Q) );

    *p_istep = istep;
	*p_energy = energy;
	*p_eint = eint;
	*p_lastmove = lastmove;

        for (i=0; i<ntor; i++) {
	    (void) fgets(myline, LINELEN, fp);
            sscanf(myline, "%lf", &(S->tor[i]) ); /* input torsions are in degrees */
            S->tor[i] = DegreesToRadians( S->tor[i] );    /* now in radians */
        }
    }
    return( status );
}
/* EOF */
