/*

 $Id: readmap.cc,v 1.24 2014/06/27 01:17:43 mp Exp $

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
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include "readmap.h"
#include "timesys.h" // for struct tms

extern char dock_param_fn[];
extern char *programname;
extern int ignore_errors;
extern int ElecMap;
extern int debug;

Statistics readmap( char           line[LINE_LEN],
                    const int      outlev,
                    const Clock&   jobStart,
                    const struct tms& tmsJobStart,
                    const Boole    B_charMap,
    /* not const */ Boole *const   P_B_HaveMap, 
                    const int      num_maps, 
                    const GridMapSetInfo *const info,
		    #include "map_declare.h"
                    const char     map_type,
		    FILE *logFile
                  )

{
    FILE *map_file;

    char FileName[PATH_MAX];
    char FldFileName[PATH_MAX];
    char GpfName[PATH_MAX];
    char mmFileName[PATH_MAX];
    const static char xyz_str[]="xyz";
    char C_mapValue;  // Caution: may be unsigned on some platforms. M Pique
    char map_line[LINE_LEN];
    char inputline[LINE_LEN];
    char atom_type_name[MAX_CHARS];

    Real cen[SPACE];
    Real spacing = 0.;
    double map_max;
    double map_min;
    double map_total;

    int nel[SPACE];
    int nv=0;
    int nvExpected;

    register int xyz = 0;
    register int i = 0;
    register int j = 0;
    register int k = 0;

    struct tms tms_jobEnd;
    struct tms tms_loadEnd;
    struct tms tms_loadStart;

    Clock jobEnd;
    Clock loadEnd;
    Clock loadStart;

    Statistics map_stats;


    //maps->atom_type = num_maps;

    /*
    \  ATOMIC AFFINITY or ELECTROSTATIC GRID MAP
     \  Read in active site grid map...
      \____________________________________________________________
     */

    strcpy(FileName, "");
    (void) sscanf( line, "%*s %s", FileName );
    if ( openFile( FileName, "r", &map_file, jobStart,tmsJobStart,TRUE, logFile)) {
        *P_B_HaveMap = TRUE;
        if (debug > 0) {
            for (i=0; i < info->num_atom_types; i++) {
                (void) fprintf(logFile, "info->atom_type_name[%d] = \"%s\"\n", i, info->atom_type_name[i] );
            }
        }

        if (map_type == 'e') {
            strcpy(atom_type_name, "e");
        } else if (map_type == 'd') {
            strcpy(atom_type_name, "d");
        } else {
            strcpy(atom_type_name, info->atom_type_name[num_maps]);
        }
	if(outlev>=LOGRECREAD)
        pr( logFile, "Opened Grid Map %d (%s):\t\t\t\t%s\n", num_maps+1, atom_type_name, FileName );

        if (outlev>=LOGRECREAD && !ignore_errors) {
            pr( logFile, "Checking header information.\n" );
        }
         /*
         \ Check header lines of grid map... 
         /
         \ :Line 1  GRID_PARAMETER_FILE 
        */
        if (fgets(inputline, LINE_LEN, map_file) == NULL) {
            warn_bad_file( FileName,"Could not read GRID_PARAMETER_FILE line." );
        } else {
            (void) sscanf(inputline, "%*s %s", GpfName);
        } /* endif */
         /*
         \ :Line 2  GRID_DATA_FILE 
        */
        if (fgets(inputline, LINE_LEN, map_file) == NULL) {
            warn_bad_file( FileName,"Could not read \".fld\" GRID_DATA_FILE line." );
        } else {
            (void) sscanf(inputline, "%*s %s", FldFileName);
            if (!ignore_errors) {
                check_header_line( FldFileName, info->FN_gdfld );
            } /* endif */
        } /* endif */
         /*
         \ :Line 3  MACROMOLECULE 
        */
        if (fgets(inputline, LINE_LEN, map_file) == NULL) {
            warn_bad_file( FileName,"Could not read MACROMOLECULE line." );
        } else {
            (void) sscanf(inputline,"%*s %s", mmFileName);
            check_header_line( mmFileName, info->FN_receptor );
        } /* endif */
         /*
         \ :Line 4  SPACING 
        */
        if (fgets(inputline, LINE_LEN, map_file) == NULL) {
            warn_bad_file( FileName,"Could not read SPACING line." );
        } else {
            (void) sscanf(inputline,"%*s " FDFMT, &spacing);
            check_header_float(spacing, info->spacing, "grid point spacing", FileName );
        } /* endif */
         /*
         \ :Line 5  NELEMENTS 
        */
        if (fgets(inputline, LINE_LEN, map_file) == NULL) {
            warn_bad_file( FileName,"Could not read NELEMENTS line." );
        } else {
            (void) sscanf(inputline,"%*s %d %d %d", &nel[X], &nel[Y], &nel[Z]);
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                //maps->num_points[xyz] = nel[xyz];
                //maps->num_points1[xyz] = nel[xyz] + 1;
                check_header_int( nel[xyz], info->num_points[xyz], xyz_str[xyz], FileName );
            } /* xyz */
        } /* endif */
         /* 
         \ :Line 6  CENTER
        */
        if (fgets(inputline, LINE_LEN, map_file) == NULL) {
            warn_bad_file( FileName,"Could not read CENTER line." );
        } else {
            (void) sscanf(inputline,"%*s " FDFMT3, &cen[X], &cen[Y], &cen[Z]);
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                //maps->center[xyz] = cen[xyz];
                check_header_float(cen[xyz], info->center[xyz], "grid-map center", FileName );
            } /* xyz */
        } /* endif */
    } /* endif */
    flushLog;

    /*
    \   Now find the extrema of the grid-map energies,
     \  While reading in the values...
      \____________________________________________________________
     */

    map_max = -BIG;
    map_min =  BIG;
    map_total = 0.;
    nvExpected = info->num_points1[X] * info->num_points1[Y] * info->num_points1[Z];
    nv = 0;

    if(outlev>=LOGRECREAD) {
    pr( logFile, "Number of grid points expected in  x-dimension:  %d\n", info->num_points1[X] );
    pr( logFile, "Number of grid points expected in  y-dimension:  %d\n", info->num_points1[Y] );
    pr( logFile, "Number of grid points expected in  z-dimension:  %d\n", info->num_points1[Z] );
    pr( logFile, "Looking for %d energies from Grid Map %d... \n", nvExpected, num_maps+1 );
    flushLog;
    }
    loadStart = times( &tms_loadStart );

    for ( k = 0;  k < info->num_points1[Z];  k++) {
        for ( j = 0;  j < info->num_points1[Y];  j++) {
            for ( i = 0;  i < info->num_points1[X];  i++) {
		float thisval;
	        if (fgets( map_line, LINE_LEN, map_file) == NULL) continue; // eof or error
                if (B_charMap) {
                    if (sscanf( map_line,  "%c",  &C_mapValue ) != 1) continue;
		    thisval=mapc2f(C_mapValue);
                } else {
                    if( sscanf( map_line,  "%f",  &thisval) != 1) continue;
                }
		SetMap(map,info,k,j,i,num_maps,thisval);
		map_max = max( map_max, thisval );
                map_min = min( map_min, thisval );
                map_total += thisval;
		nv++;
            }
        }
    }

    map_stats.number = nv;
    map_stats.minimum = map_min;
    map_stats.maximum = map_max;
    map_stats.mean = map_total / (double) nv;

    /*
    if (map_stats.number > 1) {
        double deviation = 0.;
        double sum_squares = 0.;
        for ( k = 0;  k < info->num_points1[Z];  k++) {
            for ( j = 0;  j < info->num_points1[Y];  j++) {
                for ( i = 0;  i < info->num_points1[X];  i++) {
#ifdef MAPSUBSCRIPT 
                    deviation = map[k][j][i][num_maps] - map_stats.mean;
#else
                    deviation = GetMap(map,info,k,j,i,num_maps) - map_stats.mean;
#endif
                    sum_squares += deviation * deviation;
                }
            }
        }
        map_stats.standard_deviation = sqrt(sum_squares / (map_stats.number - 1));
    } else {
        map_stats.standard_deviation = 0.;
    }
    */

    if(outlev >= LOGRECREAD) pr( logFile, "Closing file.\n" );
    fclose( map_file );
    if(outlev >= LOGRECREAD) {
    pr( logFile, "%d energies found for map %d\n", nv, num_maps+1 );
    if (map_type == 'e') {
        pr( logFile, "Minimum electrostatic potential = %.2f,  maximum electrostatic potential = %.2f\n\n", map_min, map_max );
    } else {
        pr( logFile, "Minimum energy = %.2f,  maximum energy = %.2f\n\n", map_min, map_max );
    }
    pr( logFile, "Time taken (s): " );

    loadEnd = times( &tms_loadEnd );
    timesys( loadEnd - loadStart, &tms_loadStart, &tms_loadEnd, logFile);

    pr( logFile, "\n" );
    } // if outlev

    if (nv != nvExpected ) {
        char message[LINE_LEN];
        prStr( message, "\n%s: too few values read in. Check grid map '%c' !\n\n",
		 programname, map_type  );
        pr_2x( stderr, logFile, message );

        jobEnd = times( &tms_jobEnd );
        timesys( jobEnd - jobStart, &tmsJobStart, &tms_jobEnd, logFile);
        pr_2x( logFile, stderr, UnderLine );

        /* END PROGRAM */
        stop("ERROR in readmap");
    } 

    flushLog;

    return map_stats;
}

Real mapc2f(const char numin)
{
    Real numout;
    if (numin == 0) {
        numout = 0.;
    } else if (numin > 0) {
        numout = numin * 10.;
    } else {
        numout = numin /10.;
    }
    return numout;
}
/* EOF */
