/*

 $Id: torNorVec.cc,v 1.11 2012/02/07 20:47:30 mp Exp $

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
#include <stdlib.h>
#include "torNorVec.h"
#ifdef DEBUG
#include <ctype.h>
#endif /* DEBUG */

void torNorVec( const Real crdpdb[MAX_ATOMS][SPACE],
		const int ntor,
		const int tlist[MAX_TORS][MAX_ATOMS],
		/* not const */ Real vt[MAX_TORS][SPACE],
		const int outlev,
		FILE *logFile)
{

    register int xyz = 0;
    register int j = 0;

    Real magVec = 0.;
    Real imagVec = 0.;
    Real v[SPACE];

    /*_____________________________________________________________
      | Calculate normal vectors of torsion bonds,                 |
      |____________________________________________________________|
      |      Note: torsions must be rotated from leaves to root    |
      |       or this pre-calculation of normal vectors will fail. |
      |____________________________________________________________| */

    if(outlev>=LOGLIGREAD)
    pr( logFile, "\nCalculating normalized torsion vectors.\n\n" );

    for (j=0; j<ntor; j++) {
	for (xyz = 0;  xyz < SPACE;  xyz++) {
	    v[xyz] = crdpdb[ tlist[j][ATM2] ][xyz] - crdpdb[ tlist[j][ATM1] ][xyz];

#ifdef DEBUG
		pr( logFile, "\n__norm__ Torsion %d, crdpdb[ %d,ATM2 ][%c] = %.3f, crdpdb[ %d,ATM1 ][%c] = %.3f\n",
		    j, tlist[j][ATM2], toascii(88+xyz), crdpdb[tlist[j][ATM2]][xyz], 
		    tlist[j][ATM1], toascii(88+xyz), crdpdb[tlist[j][ATM1]][xyz] );
		flushLog; 
#endif /* DEBUG */

	} /* xyz */

	magVec = hypotenuse( v[X], v[Y], v[Z] );  /* Magnitude of vector v[xyz] */

	if (magVec == 0.) {
	    char error_message[LINE_LEN];
	    prStr( error_message, "Torsion %d, normal vector, magVec, is 0; imminent division by zero caught.", j );
	    stop( error_message );
	}

	imagVec = 1. / magVec;

	for (xyz = 0;  xyz < SPACE;  xyz++) {
	    vt[j][xyz] = v[xyz] * imagVec;  /* Normalize. */
	} /* xyz */

    } /* j */
    return;
}

void update_torsion_vectors( const Real crdpdb[MAX_ATOMS][SPACE],
                             const int ntor,
                             const int  tlist[MAX_TORS][MAX_ATOMS],
                             /* not const */ Real vt[MAX_TORS][SPACE],
                             /* not const */ Molecule *ligand,
                             const int debug,
			     const int outlev,
			     FILE *logFile)
{ // Update the unit vectors for the torsion rotations
    register int i=0, j=0;
    if (debug > 0) {
        pr(logFile, "Calculating unit vectors for each torsion.\n\n");
    }
    torNorVec(crdpdb, ntor, tlist, vt,outlev,logFile);
    for (i = 0; i < MAX_TORS; i++) {
        ligand->vt[i][X] = vt[i][X];
        ligand->vt[i][Y] = vt[i][Y];
        ligand->vt[i][Z] = vt[i][Z];
        for (j = 0; j < MAX_ATOMS; j++) {
            ligand->tlist[i][j] = tlist[i][j];
        }
    }
    return;
} // Update the unit vectors for the torsion rotations

/* EOF */
