/*

 $Id: qtransform.cc,v 1.15 2012/02/07 05:14:55 mp Exp $

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
#include "qtransform.h"
#include <stdio.h>
#include <string.h>


void qtransform( const Coord& T,
                 const Quat&  q,
                 /* not const */ Real tcoord[MAX_ATOMS][SPACE],
                 const int   natom)

/******************************************************************************/
/*      Name: qtransform                                                      */
/*  Function: Accelerated quaternion transformation                           */
/*            Performs both a rigid-body translation and rotation.            */
/*            Assumes quaternion is normalized outside this routine.          */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*   Authors: Garrett M. Morris, The Scripps Research Institute               */
/*            David Goodsell, UCLA                                            */
/*      Date: 11/23/94                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: q, tcoord, natom                                                */
/*            q[] contains (xt,yt,zt,x,y,z,w(radians))                        */
/*   Returns: tcoord                                                          */
/*   Globals: QUAT, MAX_ATOMS                                                 */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 05/15/92 GMM     Translated into C.                                        */
/* 11/23/92 GMM     Introduced (15+,9*) version, replacing (12+,36*) version. */
/* 11/23/92 GMM     Introduced h and inv_qmag, since * is faster than /.      */
/******************************************************************************/
{
    register int a;
    Coord  tmp;
    double w, x, y, z;
    double tx, ty, tz;
    double omtxx;
    double twx, txy, txz;
    double twy, tyy, tyz;
    double twz, tzz;
    double r11, r12, r13, r21, r22, r23, r31, r32, r33;

    w = q.w;
    x = q.x;
    y = q.y;
    z = q.z;


/*  12 adds, 36 multiplies...     */
/*    r11 = 1. - 2.*y*y - 2.*z*z; */
/*    r12 =      2.*x*y + 2.*w*z; */
/*    r13 =      2.*x*z - 2.*w*y; */
/*    r21 =      2.*x*y - 2.*w*z; */
/*    r22 = 1. - 2.*x*x - 2.*z*z; */
/*    r23 =      2.*y*z + 2.*w*x; */
/*    r31 =      2.*x*z + 2.*w*y; */
/*    r32 =      2.*y*z - 2.*w*x; */
/*    r33 = 1. - 2.*x*x - 2.*y*y; */

/*  15 adds, 9 multiplies...      */
    tx  = x+x;
    ty  = y+y;
    tz  = z+z;

    twx = w*tx;
    omtxx = 1. - x*tx;
    txy = y*tx;
    txz = z*tx;

    twy = w*ty;
    tyy = y*ty;
    tyz = z*ty;

    twz = w*tz;
    tzz = z*tz;

    r11 = 1. - tyy - tzz;
    r12 =      txy - twz;
    r13 =      txz + twy;
    r21 =      txy + twz;
    r22 = omtxx    - tzz;
    r23 =      tyz - twx;
    r31 =      txz - twy;
    r32 =      tyz + twx;
    r33 = omtxx    - tyy;

    for (a = 0;  a < natom;  a++) {
        tmp.x = ((double)tcoord[a][X])*r11 + ((double)tcoord[a][Y])*r21 + ((double)tcoord[a][Z])*r31 + T.x;
        tmp.y = ((double)tcoord[a][X])*r12 + ((double)tcoord[a][Y])*r22 + ((double)tcoord[a][Z])*r32 + T.y;
        tmp.z = ((double)tcoord[a][X])*r13 + ((double)tcoord[a][Y])*r23 + ((double)tcoord[a][Z])*r33 + T.z;
        tcoord[a][X] = tmp.x;
        tcoord[a][Y] = tmp.y;
        tcoord[a][Z] = tmp.z;
    }
}

void reorient( FILE *const logFile, 
               const int true_ligand_atoms, 
               const char atomstuff[MAX_ATOMS][MAX_CHARS],
               /* not const */ Real crdpdb[MAX_ATOMS][SPACE],  // original PDB coordinates from input
               const Real charge[MAX_ATOMS],
               const int type[MAX_ATOMS],
               const ParameterEntry  parameterArray[MAX_ATOM_TYPES],
               const Quat& q_reorient,
               const Coord& origin,
               const int ntor,
               const int tlist[MAX_TORS][MAX_ATOMS],
               /* not const */ Real vt[MAX_TORS][SPACE],
               /* not const */ Molecule *ptr_ligand,
               const int debug, const int outlev )
 {
    // Print out the un-reoriented coordinates
    pr( logFile, "\nUn-reoriented ligand's coordinates:\n" );
    pr( logFile, "-----------------------------------\n\n" );
    print_PDBQT( logFile, "UN-REORIENTED: ", true_ligand_atoms, atomstuff, crdpdb, charge, parameterArray, type, "\n" );

    // Print message about q_reorient
    print_q_reorient_message( logFile, q_reorient );

    // Apply the rotation defined by q_reorient to the input coordinates of the ligand, "crdpdb"
    qtransform( origin, q_reorient, crdpdb, true_ligand_atoms );

    // Update the unit vectors for the torsion rotations
    update_torsion_vectors( crdpdb, ntor, tlist, vt, ptr_ligand, 
     debug, outlev, logFile);
    
    // Print out the re-oriented coordinates
    pr( logFile, "Reoriented ligand's coordinates:\n" );
    pr( logFile, "--------------------------------\n\n" );
    print_PDBQT( logFile, "REORIENTED:  ", true_ligand_atoms, atomstuff, crdpdb, charge, parameterArray, type, "\n" );
}
/* EOF */
