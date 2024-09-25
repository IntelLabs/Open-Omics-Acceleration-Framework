/*

 $Id: stateLibrary.cc,v 1.24 2014/07/02 20:25:40 mp Exp $

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
#include "stateLibrary.h"
#include "qmultiply.h"

static void writeStateAQ( FILE *const fp, /* not const */ State S, const char);

void initialiseState( /* not const */ State *const S )
{
    register int i;
    S->T.x = 0.0;
    S->T.y = 0.0;
    S->T.z = 0.0;
    S->Q = identityQuat();
    S->ntor = 0;
    for (i = 0; i  < MAX_TORS;  i++ ) {
        S->tor[i] = 0.0;
    }
    S->Center.x = 0.0;
    S->Center.y = 0.0;
    S->Center.z = 0.0;
}
void initialiseQuat ( Quat *Q )
{
	*Q = identityQuat();
}

void copyState( /* not const */ State *const D,  /* Destination -- copy to here */
                const State& S ) /* Source      -- copy this.   */
{
    register int i;
        
    D->T = S.T;
    
    D->Q = S.Q;
 
    D->ntor   = S.ntor;
 
    for ( i=0; i < S.ntor; i++ ) {
            D->tor[i] = S.tor[i];
    }

    D->hasEnergy = S.hasEnergy;

    D->e = S.e;
    D->Center = S.Center;
}

static double RadCanonicalDeg( const double x) {
   // convert torsion value from radians to degrees, 
   // wrapping angles at all points.
   // Probably could be simplified.
	return WrpDeg( ModDeg( RadiansToDegrees( WrpRad( ModRad( x )))));
}

void printState( FILE *const fp, 
                 /* not const */ State S, 
                 const int detail )
{
    register int i;
    AxisAngle aa;

    switch( detail ) {
        case 0:
        case 1:
            writeStateAQ(fp,S,'A'); // short format, axis-angle convention
            break;

        case 2:
        default:
	    // comprehensive format, no initial newline
            (void)fprintf( fp, "STATE VARIABLES:\n________________\n\n" );
            (void)fprintf( fp, "Translation x,y,z         = %.3f %.3f %.3f\n", S.T.x, S.T.y, S.T.z );
            (void)fprintf( fp, "Quaternion x,y,z,w        = %.6f %.6f %.6f %.6f\n", S.Q.x, S.Q.y, S.Q.z, S.Q.w );
	    aa = QuatToAxisAngle( S.Q );
            (void)fprintf( fp, "Axis-Angle nx,ny,nz,angle = %.3f %.3f %.3f %.3lf\n", aa.nx, aa.ny, aa.nz, RadCanonicalDeg(aa.ang) );
            (void)fprintf( fp, "Center x,y,z         = %.3f %.3f %.3f\n", S.Center.x, S.Center.y, S.Center.z );
            (void)fprintf( fp, "Number of Torsions        = %d\n", S.ntor );
            if (S.ntor > 0) {
                (void)fprintf( fp, "Torsions (degrees)        =");
                for (i=0; i<S.ntor; i++) {
                    S.tor[i] = WrpRad( ModRad( S.tor[i] ) );
                }
                for (i=0; i<S.ntor; i++) {
                    pr( fp, " %.2lf", RadCanonicalDeg( S.tor[i]));
                    //if ((B_isTorConstrained[i] == 1) && B_ShowTorE) {
                        //pr( fp, ", Energetic penalty = %uhd\n", US_TorE[i]);
                    //} else {
                        //pr( fp, "\n");
                    //}
                }
            }
            (void)fprintf( fp, "\n\n");
            break;

        case 3:
            // Writes only the translation component of the state
            (void)fprintf( fp, "%.3f %.3f %.3f", S.T.x, S.T.y, S.T.z );
            break;
        case 4:
	    // Writes in compact format with underscore separations, no newline
	    // Used by PrintPopulationStatisticsVerbose (M Pique 2010-03)
            (void)fprintf( fp, "%.3f_%.3f_%.3f", S.T.x, S.T.y, S.T.z );
            (void)fprintf( fp, "_%.6f_%.6f_%.6f_%.6f", S.Q.x, S.Q.y, S.Q.z, S.Q.w );
            for (i=0; i<S.ntor; i++) {
                    pr( fp, "_%.3lf", RadCanonicalDeg(S.tor[i]) );
            }
            break;
        case 5:
            writeStateAQ(fp,S,'Q'); // short format, quaternion convention
            break;
        case 6:
	    // verbose state including center and ntors
	    //  for writing "Detailed state:" line into DLG  (M Pique 2010-05) 
	    // Torsions are written out in degrees
	    // Does not append newline
            (void)fprintf( fp, " trans %.3f %.3f %.3f", S.T.x, S.T.y, S.T.z );
            (void)fprintf( fp, " quatxyzw %.6f %.6f %.6f %.6f", 
	      S.Q.x, S.Q.y, S.Q.z, S.Q.w );
            (void)fprintf( fp, " center %.3f %.3f %.3f", 
	      S.Center.x, S.Center.y, S.Center.z );
            (void)fprintf( fp, " ntor %d", S.ntor);
            for (i=0; i<S.ntor; i++) pr( fp, " %.4lf", RadCanonicalDeg(S.tor[i]));
            break;
    }
}

static void writeStateAQ( FILE *const fp, /* not const */ State S, const char convention )
{


    // Write translation.
    (void)fprintf( fp, "%7.3f %7.3f %7.3f  ", S.T.x, S.T.y, S.T.z );

    switch (convention) {
       case 'Q':
       case 'q':
	    // quaternion
	    (void)fprintf( fp, "%.5f %.5f %.5f %.5f  ", S.Q.x, S.Q.y, S.Q.z, S.Q.w);
	    break;

       case 'A':
       case 'a':
       default:
	    //  axis-angle.
	    AxisAngle aa = QuatToAxisAngle( S.Q );

	    float ang = RadCanonicalDeg(WrpRad( ModRad( aa.ang )));
	    (void)fprintf( fp, "%6.3f %6.3f %6.3f %6.3f  ",
	      aa.nx, aa.ny, aa.nz, ang);
	    break;
    }
    // Write torsion angles.
    if (S.ntor > 0) {
        for (int i=0; i<S.ntor; i++) {
            pr( fp, " %7.2lf", RadCanonicalDeg(S.tor[i]));
        }
    }
}
void writeState( FILE *const fp, /* not const */ State S )
{
	// simple wrapper to write state with axis-angle convention
	writeStateAQ( fp, S, 'A');
}

int checkState( FILE *const logFile, const State *const D)
{
// return 0 if failure   1 if OK
    register int i;
    int retval = 1;
    double magnitude_q;
        
    if (ISNAN(D->T.x)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in x translation\n");
        retval = 0;
    }
    if (ISNAN(D->T.y)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in y translation\n");
        retval = 0;
    }
    if (ISNAN(D->T.z)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in z translation\n");
        retval = 0;
    
    }

    if (ISNAN(D->Q.x)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in x quaternion\n");
        retval = 0;
    }
    if (ISNAN(D->Q.y)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in y quaternion\n");
        retval = 0;
    }
    if (ISNAN(D->Q.z)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in z quaternion\n");
        retval = 0;
    }
    if (ISNAN(D->Q.w)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in w quaternion\n");
        retval = 0;
    }
    magnitude_q = hypotenuse4(D->Q.x,  D->Q.y,  D->Q.z,  D->Q.w);
    if (ISNAN(magnitude_q)) {
        (void)fprintf(logFile,"checkState: After computing the magnitude of quaternion, (NaN) was detected\n");
        retval = 0;
    }
 
    for ( i=0; i < D->ntor; i++ ) {
            if (ISNAN(D->tor[i])) {
                (void)fprintf(logFile,"checkState: (NaN) detected in torsion %d\n",i+1);
                retval = 0;
            }
    }

    return(retval);
}

Molecule copyStateToMolecule(const State *const S, /* not const */ Molecule *const mol) /* S is the source */
{
    register int i;
    mol->S.T = S->T;
    mol->S.Q = S->Q;
    mol->S.ntor = S->ntor;
    for (i = 0; i  < MAX_TORS;  i++ ) {
        mol->S.tor[i] = S->tor[i];
    }
    return *mol;
}
/* EOF */
