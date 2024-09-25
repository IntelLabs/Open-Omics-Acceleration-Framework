/*

 $Id: autocomm.h,v 1.23 2012/04/17 04:06:10 mp Exp $

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

/* autocomm.h */

#ifndef _AUTOCOMM
#define _AUTOCOMM

#include <sys/types.h>
#include <time.h>
/* include stdio to pick up definition of FILENAME_MAX and possibly PATH_MAX */
#include <stdio.h>
#include "typedefs.h"

/*******************************************************************************
**      Name: autocomm.h                                                      **
**  Function: Defines Constants, common to both AUTOGRID & AUTODOCK...        **
**Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
**----------------------------------------------------------------------------**
**    Author: Garrett Matthew Morris, The Scripps Research Institute          **
**      Date: 02/28/1995                                                      **
**----------------------------------------------------------------------------**
**    Inputs: none                                                            **
**   Returns: nothing                                                         **
**   Globals: all defines                                                     **
**----------------------------------------------------------------------------**
** Modification Record                                                        **
** Date     Inits   Comments                                                  **
** 02/28/95 GMM     This header was added.                                    **
*******************************************************************************/

/*
** Constants,
*/

#define FALSE        0      /* Logical constant                               */
#define TRUE         1      /* Logical constant                               */

const Real PI=3.14159265358979323846;   /* Mathematical constant, pi */
const Real TWOPI=6.28318530717958647692;
const Real HALF_PI=1.57079632679489661923;

#define X            0      /* x-coordinate                                   */
#define Y            1      /* y-coordinate                                   */
#define Z            2      /* z-coordinate                                   */
#define XYZ          3      /* Dimensions of Cartesian Space                  */
#define SPACE        3      /* Dimensions of Cartesian Space                  */

const Real APPROX_ZERO=1.0E-6; /* To avoid division-by-zero errors...            */
const Real BIG=1.0E12; /* Very large constant                            */
const unsigned int MAX_CHARS=128;    /* Number of characters in atom data & filenames  */
const unsigned int MAX_LINES=256;    /* Number of lines in parameter file              */
#ifndef PATH_MAX
#define PATH_MAX     FILENAME_MAX
#endif

#ifdef USE_XCODE
const /* not unsigned */ int LINE_LEN=140;    /* Line length in characters                      */
#else
const /* not unsigned */ int LINE_LEN=256;    /* Line length in characters                      */
#endif 

// Note: MAX_GRID_PTS is a sanity check, and could be enlarged at will   M Pique 2012
const /* not unsigned */ int MAX_GRID_PTS=1025;	/* Maximum number of grid points in 1 dimension */

const Real EINTCLAMP=100000.; /* Clamp pairwise internal energies (kcal/mol )  */

#define MAX_MAPS_PAD 0       // Use this to pad MAX_MAPS to a power of 2, for presumably-faster memory access
#define NUM_NON_VDW_MAPS 2   // Number of electrostatic and desolvation maps
#define MAX_ATOM_TYPES (16 - NUM_NON_VDW_MAPS)    /* Maximum number of atom types set to keep MAX_MAPS a power of 2 */
#define MAX_MAPS (MAX_ATOM_TYPES + NUM_NON_VDW_MAPS + MAX_MAPS_PAD) /* Maximum number of energy maps        */
                            /* 0,1,2,... are for atomic interactions          */
                            /* last two are for electrostatics and desolvation */

#define VECLENMAX    16     /* For AVS fld files...                           */

#define UnderLine "________________________________________________________________________________\n\n"

/*
** Common Macros...
*/

#define pr              (void) fprintf
#define pr_2x           print_2x
#define prStr           (void) sprintf
#define flushLog        (void) fflush(logFile)

//FIXME: this should eventually become an inline function
#define dist(x1,y1,z1,x2,y2,z2,r) _dx=((x2)-(x1)),_dy=((y2)-(y1)),_dz=((z2)-(z1)),r=sqrt(_dx*_dx + _dy*_dy + _dz*_dz)

/*
** New types...
*/


#include "typedefs.h"


typedef int Boole; // this could become the C++ type bool  (MP 2011 TODO)


typedef struct AtomDesc {

	Real crd[XYZ];
	Real q;
	int   type;

	} AtomDesc;


/*
** Note the following differing definitions of "times" and "time":-
**
** Arch. times()				time()
** ----- ----------------------------------	--------------------------
** Sun	 clock_t times(struct tms *buffer);	time_t time(time_t *tloc);
**
** ----- -------------------------------									-----------------------------------------------
** Sun	 void srand48(long seedval);										struct tm *localtime(const time_t *clock);
** timesys and timesyshms used to use Clock, should use time_t
**
*/

#define Clock clock_t


/*
 * assert that quaternions are OK
 */
#include <assert.h> // for assert in assertQuatOK
#include <math.h> // for sqrt in assertQuatOK

const Real ONE_MINUS_EPSILON=0.999;
const Real  ONE_PLUS_EPSILON=1.001;

/*
 * void assertQuatOK( const Quat q )
 * {
 *     register double mag4 = hypotenuse4( q.x, q.y, q.z, q.w );
 *     assert((mag4 > ONE_MINUS_EPSILON) && (mag4 < ONE_PLUS_EPSILON));
 * }
 */
#define assertQuatOK( q ) {fflush(logFile);register double aQOK_mag4 = hypotenuse4( (q).x, (q).y, (q).z, (q).w ); assert((aQOK_mag4 > ONE_MINUS_EPSILON) && (aQOK_mag4 < ONE_PLUS_EPSILON)); }


#endif /*_AUTOCOMM*/

/*
** EOF
*/
