/*

 $Id: torsion.cc,v 1.9 2010/10/01 22:51:40 mp Exp $

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
#include "constants.h"
#include "structs.h"
#include "torsion.h"


void torsion( const State& now,
    /* not const */ Real crd[MAX_ATOMS][SPACE],
              const Real v[MAX_TORS][SPACE],
              const int tlist[MAX_TORS][MAX_ATOMS],
              const int ntor )

/******************************************************************************/
/*      Name: torsion                                                         */
/*  Function: Apply the Torsion rotation(s) to the Small Molecule.            */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 03/14/94                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: now, v, tlist, ntor                                             */
/*   Returns: crd                                                             */
/*   Globals: MAX_TORS, SPACE, MAX_ATOMS.                                     */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 05/14/92 GMM     Translated into C                                         */
/* 03/14/94 GMM     optimized equations, redundant operations removed...      */
/******************************************************************************/

{
    register double crdtemp[SPACE];
    register double  d[SPACE];
    register double sv[SPACE];
    register double ov[SPACE];
    register double  k[SPACE][SPACE];
    register double s, c, o, vni, this_tor;                /* "o" is: 1. - c, "One minus c". */
    register int n, a, mvatm, atmnum, numatmmoved;

    for (n = 0;  n < ntor;  n++) {                /* "n"-th Torsion */
          s = sin(this_tor = ModRad(now.tor[n]));
          o = 1. - (c = cos(this_tor));
          /*
          atmnum = Serial number of Atom 0 in Torsion n, (0)--(1) 
          This atom is at the origin of the current torsion
          vector, about which we are rotating...
          */
          atmnum = tlist[n][0];
          crdtemp[X] = (double)crd[atmnum][X];
          crdtemp[Y] = (double)crd[atmnum][Y];
          crdtemp[Z] = (double)crd[atmnum][Z];

          sv[X] = s * (vni = v[n][X]);
          k[X][X] = (ov[X] = o * vni) * vni + c;

          sv[Y] = s * (vni = v[n][Y]);
          k[Y][Y] = (ov[Y] = o * vni) * vni + c;

          sv[Z] = s * (vni = v[n][Z]);
          k[Z][Z] = (ov[Z] = o * vni) * vni + c;

          k[Y][Z]  =  v[n][Y] * ov[Z]  -  sv[X];
          k[Z][X]  =  v[n][Z] * ov[X]  -  sv[Y];
          k[X][Y]  =  v[n][X] * ov[Y]  -  sv[Z];

          k[Z][Y]  =  v[n][Z] * ov[Y]  +  sv[X];
          k[X][Z]  =  v[n][X] * ov[Z]  +  sv[Y];
          k[Y][X]  =  v[n][Y] * ov[X]  +  sv[Z];

          numatmmoved = tlist[n][NUM_ATM_MOVED] + 3; 
          for (a = 3;  a < numatmmoved;  a++ )  {        
              mvatm = tlist[n][a]; /* mvatm = Serial Num of Atom to be moved by this Torsion */
              d[X] = (double)crd[mvatm][X] - crdtemp[X];
              d[Y] = (double)crd[mvatm][Y] - crdtemp[Y];
              d[Z] = (double)crd[mvatm][Z] - crdtemp[Z];
              crd[mvatm][X] = (double)crdtemp[X] + d[X] * k[X][X] + d[Y] * k[X][Y] + d[Z] * k[X][Z]; 
              crd[mvatm][Y] = (double)crdtemp[Y] + d[X] * k[Y][X] + d[Y] * k[Y][Y] + d[Z] * k[Y][Z]; 
              crd[mvatm][Z] = (double)crdtemp[Z] + d[X] * k[Z][X] + d[Y] * k[Z][Y] + d[Z] * k[Z][Z]; 
          }/*a*/
    } /*n*/
}
/* EOF */
