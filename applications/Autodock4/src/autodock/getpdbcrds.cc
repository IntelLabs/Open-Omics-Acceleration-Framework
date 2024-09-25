/*

 $Id: getpdbcrds.cc,v 1.8 2014/02/01 05:14:53 mp Exp $

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
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "getpdbcrds.h"


extern char *programname;


int getpdbcrds( const char *const rms_ref_crds_FN,
		/* not const */ Real ref_crds[MAX_ATOMS][SPACE], FILE *logFile)
{
    int ii=0;
    int natoms=0;
    char line[LINE_LEN];
    char str[4][WORDLEN];
    char rec5[5];
    FILE *rms_ref_FilePtr;

    if ( !openfile( rms_ref_crds_FN, "r", &rms_ref_FilePtr, logFile)) {
	fprintf( logFile, "%s: ERROR!  Sorry, could not open file \"%s\" for reading.\n", programname,  rms_ref_crds_FN );
	return -1;
    }

    pr (logFile, "\nRMS Reference Coordinates from \"%s\":-\n\n", rms_ref_crds_FN);

    while ( fgets(line, LINE_LEN, rms_ref_FilePtr) != NULL ) {

	for (ii = 0; ii < 4; ii++) {
	    rec5[ii] = (char)tolower( (int)line[ii] );
	}
	if (equal(rec5,"atom",4) || equal(rec5,"heta",4)) {

	    if (natoms < MAX_ATOMS) {
		sscanf( &line[30], "%s %s %s", str[X], str[Y], str[Z] );
		ref_crds[natoms][X] = atof( str[X] );
		ref_crds[natoms][Y] = atof( str[Y] );
		ref_crds[natoms][Z] = atof( str[Z] );
		pr (logFile, "Atom %5d,  x,y,z = %8.3f %8.3f %8.3f\n", natoms+1, ref_crds[natoms][X], ref_crds[natoms][Y], ref_crds[natoms][Z]);
	    } else {
		fprintf( logFile, "%s: ERROR!  Sorry, too many atoms in file \"%s\"\n", programname,  rms_ref_crds_FN );
		return -1;
	    }
	    ++natoms;
	} /* End if "atom" or "heta" */

    } /* End while there's a line to read... */

    /* Close the reference coordinates file... */
    (void) fclose(rms_ref_FilePtr);

    return natoms;
}
/* EOF */
