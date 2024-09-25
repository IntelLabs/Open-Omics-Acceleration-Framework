/*

 $Id: prTorConList.cc,v 1.9 2014/06/12 01:44:07 mp Exp $

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
#include "prTorConList.h"

/*
 * Number of columns for torsion profile output.
 * Keep NCOLS less than 72, to ensure whole line
 * fits into 80 columns.
 */

#define NCOLS 64

/* Number of rows for torsion profile output. */

#define NROWS 16


#define AddNewHardCon(iCon,low,upp)  F_TorConRange[i][iCon][LOWER]=low;F_TorConRange[i][iCon][UPPER]=upp

void prTorConList( const int    ntor,
		    const Boole B_isTorConstrained[MAX_TORS],
		    const unsigned short US_torProfile[MAX_TORS][NTORDIVS],
		    /* not const */ Real  F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		    /* not const */ int   N_con[MAX_TORS], 
		    const int outlev,
		    FILE *logFile)

{
    char graph[NROWS][LINE_LEN];

    int col = 0;
    register int i = 0;
    int iCon = 0;
    int inc;
    int j = 0;
    int jncol216;
    int jncol2;
    int N_ConChecks;

    Real lower;
    Real upper;

    unsigned short US_row = 0;
    unsigned short US_yscale;

    jncol2 = (int)((Real)NCOLS / 2. - 2.);
    jncol216 = jncol2 - 16;

    pr( logFile, UnderLine );
    pr( logFile, "\n\nTorsion Constraint List\n");
    pr( logFile,     "_______________________\n\n");
    for (i=0; i<ntor; i++) {
	switch ( B_isTorConstrained[i] ) {
	    case 0:
	    default:
		pr( logFile, "Torsion %d is Not Constrained.\n", i+1);
		break;
	    case 1:
		pr( logFile, "\nTorsion %d, Gaussian Constraint Profile: (gausstorcon)\n\n", i+1);
		for ( j = 0;  j <= NROWS;  j++ ) {
		    /*
		    ** Initially all blanks,
		    */
		    strcpy(graph[j],"                                                                       ");
		}/*j*/
		inc = NTORDIVS / NCOLS;
		inc = ((inc<=0)? 1 : inc);
		US_yscale = (unsigned short) (65536 / NROWS);
		col = 0;
		for (j=0; j<NTORDIVS; j += inc) {
		    US_row = US_torProfile[ i ][j] / US_yscale;
		    graph[US_row][col] = '*';
		    col++;
		}
		pr( logFile,"Energy|\n      |\n");
		for (j = (NROWS - 1); j >= 0; j--) {
		    pr(logFile, "%5uhd |%s\n", US_yscale*(unsigned short)j, graph[j]);
		}
		pr( logFile,"------|");
		for (j=0;j<NCOLS;j++) {
		    pr(logFile,"-");
		} 
		pr(logFile,"\n    -180");
		flushLog;
		for (j=0;j<jncol2;j++) {
		    pr(logFile," ");
		    /* DEBUG: pr(logFile,"A-<%d> [%d]",j,jncol2); */
		} 
		pr(logFile,"0");
		flushLog;
		for (j=0; j<jncol2; j++) {
		    pr(logFile," ");
		    /* DEBUG: pr(logFile,"B-<%d> [%d]",j,jncol2); */
		} 
		pr(logFile,"180\n       ");
		flushLog;
		for (j=0; j<jncol216; j++) {
		    pr(logFile," ");
		    /* DEBUG: pr(logFile,"C-<%d> [%d]",j,jncol2); */
		} 
		pr(logFile,"      Relative Torsion Angle / degrees\n\n");
		pr(logFile,"\n\nTors\tAngle\tDegs\tPenalty\n");
		pr(logFile,"____\t_____\t_____\t_______\n\n");
		for (j=0; j<NTORDIVS; j++) {
		    pr(logFile,"%3d\t%3d\t%5.1f\t%6uhd\n",i+1,j,WrpDeg( ModDeg( RadiansToDegrees( DivsToRadians(j) )) ), US_torProfile[ i ][j]);
		}
		flushLog;
		break;
	    case 2:
		pr( logFile, "\nTorsion %d is Hard Constrained: (hardtorcon)\n", i+1);
		N_ConChecks = N_con[i];
		for (j=0; j<N_ConChecks; j++) {
		    lower = F_TorConRange[i][j][LOWER];
		    upper = F_TorConRange[i][j][UPPER];
		    if (lower < -180.) {
			AddNewHardCon( j, -180., upper );
			iCon = N_con[i] + 1;
			if ( iCon < MAX_TOR_CON ) {
			    pr(logFile,"Splitting requested constraint into two, since lower value less than -180.deg.\n");
			    lower += 360.;
			    upper  = 180.;
			    pr(logFile,"Placing %.1f into lower and %.1f into upper limits of HardCon %d\n",lower,upper,N_con[i]+1);
			    F_TorConRange[i][N_con[i]][LOWER]=lower;
			    F_TorConRange[i][N_con[i]][UPPER]=upper;
			    N_con[i] = iCon;
			} else {
			    pr(logFile,"I attempted to split the requested constraint %d into two, since lower value is less than -180.deg.\n",j+1);
			    pr(logFile,"\n\nBut I'm sorry, you can only have %d (=MAX_TOR_CON) torsion constraints.\nIf you need more, change the \"#define MAX_TOR_CON\" line in \"constants.h\".\n\n",MAX_TOR_CON);
			}
		    } else if (upper > 180.) {
			AddNewHardCon( j, lower, 180. );
			iCon = N_con[i] + 1;
			if ( iCon < MAX_TOR_CON ) {
			    pr(logFile,"Splitting requested constraint into two, since upper value exceeds 180.deg.\n");
			    lower  = -180.;
			    upper -= -360.;
			    pr(logFile,"Placing %.1f into lower and %.1f into upper limits of HardCon %d\n",lower,upper,N_con[i]+1);
			    F_TorConRange[i][N_con[i]][LOWER]=lower;
			    F_TorConRange[i][N_con[i]][UPPER]=upper;
			    N_con[i] = iCon;
			} else {
			    pr(logFile,"I attempted to split the requested constraint %d into two, since the upper value exceeded 180.deg.\n",j+1);
			    pr(logFile,"\n\nBut I'm sorry, you can only have %d (=MAX_TOR_CON) torsion constraints.\nIf you need more, change the \"#define MAX_TOR_CON\" line in \"constants.h\".\n\n",MAX_TOR_CON);
			}
		    } else {
		       AddNewHardCon( j, lower, upper );
		    }
		}/*j*/
		for (j=0; j<N_con[i]; j++) {
		    pr(logFile,"HardCon %d, Allowed Torsion Angles are from %.1f to %.1f degrees.\n", j+1, F_TorConRange[i][j][LOWER], F_TorConRange[i][j][UPPER]);
		    F_TorConRange[i][j][LOWER] = DegreesToRadians( F_TorConRange[i][j][LOWER] );
		    F_TorConRange[i][j][UPPER] = DegreesToRadians( F_TorConRange[i][j][UPPER] );
		}/*j*/
		break;
	}/*switch*/
    }/** for (i=0; i<ntor; i++) **/
    pr( logFile, UnderLine );
    flushLog;

}
/* EOF */
