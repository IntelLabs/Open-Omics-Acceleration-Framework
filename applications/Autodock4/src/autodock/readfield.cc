/*

 $Id: readfield.cc,v 1.11 2014/02/01 05:14:53 mp Exp $

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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#include "readfield.h"

//FIXME: tms_jobStart could be passed by reference
void readfield( /* not const */ GridMapSetInfo *const info,
                /* not const */ char line[LINE_LEN],
                const Clock& jobStart,
                const struct tms& tms_jobStart,
		const int outlev,
		FILE *logFile)
{
    FILE *fldFile;
    char rec9[9], inputline[LINE_LEN];
    double localSpacing;
    register int i = 0;

    /*
    ** GRID_DATA_FILE
    ** Read the (AVS-format) grid data file, .fld
    ** _____________________________________________________________
    **/
    (void) sscanf( line, "%*s %s", info->FN_gdfld);
     
    if ( openFile( info->FN_gdfld, "r", &fldFile, jobStart, tms_jobStart, TRUE, logFile)) {
	if(outlev >= LOGFORADT) 
        pr( logFile, "Opening Grid Map Dimensions file:\t\t%s\n", info->FN_gdfld);
    }

    /*
    ** Skip over the AVS-readable .fld  header comments, until 
    **
    ** #SPACING
    */
    while( fgets(line, LINE_LEN, fldFile) != NULL ) {
        (void) sscanf(line,"%s", rec9);
        if (equal(rec9,"#SPACING", 8)) {
        (void) sscanf(line,"%*s %lf", &localSpacing);
            info->spacing = localSpacing;
            break;
        }
    } /* endwhile */
    info->inv_spacing = 1. / (info->spacing);
    if(outlev>=LOGFORADT)
    pr( logFile, "Grid Point Spacing =\t\t\t\t%.3f Angstroms\n\n", info->spacing);

    /*
    ** #NELEMENTS 
    */
    (void) fgets(inputline, LINE_LEN, fldFile);
    (void) sscanf(inputline,"%*s %d %d %d", &info->num_points[X], &info->num_points[Y], &info->num_points[Z]);
    /* Check that these are all even numbers... */
    for (i=0; i<SPACE; i++) {
        if ( (info->num_points[i])%2 != 0 ) {
            stop("the number of user-specified grid points must be even in the \"#NELEMENTS\" line in the \".fld\" file.");
        }
    }

    if(outlev>=LOGFORADT)
    pr( logFile, "Even Number of User-specified Grid Points =\t%d x-points\n\t\t\t\t\t\t%d y-points\n\t\t\t\t\t\t%d z-points\n\n", info->num_points[X],info->num_points[Y],info->num_points[Z]);
    for (i = 0;  i < SPACE;  i++) {
        info->num_points1[i] = info->num_points[i] + 1;
	info->num_alloc[i] = info->num_points1[i]; // this is an odd number
    } /* i */
    //TODO useful? pr( logFile, "Adding the Central Grid Point makes:\t\t%d x-points\n\t\t\t\t\t\t%d y-points\n\t\t\t\t\t\t%d z-points\n\n", info->num_points1[X], info->num_points1[Y], info->num_points1[Z]);
    if ( (info->num_points[X] <= 0)||(info->num_points[Y] <= 0)||(info->num_points[Z] <= 0) ) {
        stop("insufficient grid points." );
    } else if ((info->num_points[X] > MAX_GRID_PTS)||(info->num_points[Y] > MAX_GRID_PTS)||(info->num_points[Z] > MAX_GRID_PTS)) {
        stop("too many grid points." );
    }

    /*
    ** #CENTER 
    */
    (void) fgets(inputline, LINE_LEN, fldFile);
    (void) sscanf(inputline,"%*s %lf %lf %lf", &info->center[X], &info->center[Y], &info->center[Z]);
    if(outlev >= LOGFORADT)
	    pr( logFile, "Coordinates of Central Grid Point of Maps =\t(%.3f, %.3f, %.3f)\n", info->center[X],  info->center[Y],  info->center[Z]);

    /*
    ** #MACROMOLECULE 
    */
    (void) fgets(inputline, LINE_LEN, fldFile);
    (void) sscanf(inputline,"%*s %s", info->FN_receptor);
    pr( logFile, "Macromolecule file used to create Grid Maps =\t%s\n", info->FN_receptor);

    /*
    ** #GRID_PARAMETER_FILE 
    */
    (void) fgets(inputline, LINE_LEN, fldFile);
    (void) sscanf(inputline,"%*s %s", info->FN_gpf);
    if(outlev>=LOGFORADT)
    pr( logFile, "Grid Parameter file used to create Grid Maps =\t%s\n", info->FN_gpf);

    /*
    ** Close Grid-dimensions data file 
    */
    fclose(fldFile);

    /*
    ** Determine the dimensions of the grids,
    */
    for (i = 0;  i < SPACE;  i++) {
        info->lo[i] = info->center[i] - (info->num_points[i]/2) * (info->spacing);
        info->hi[i] = info->center[i] + (info->num_points[i]/2) * (info->spacing);
    }

    if(outlev>=LOGFORADT) {
    pr( logFile, "Minimum coordinates in grid = (%.3f, %.3f, %.3f)\n",   info->lo[X], info->lo[Y], info->lo[Z]);
    pr( logFile, "Maximum coordinates in grid = (%.3f, %.3f, %.3f)\n\n", info->hi[X], info->hi[Y], info->hi[Z]);
    fflush( logFile );
    }
}
/*
** EOF
*/
