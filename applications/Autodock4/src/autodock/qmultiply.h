/*

 $Id: qmultiply.h,v 1.17 2014/06/20 23:03:52 mp Exp $

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

#ifndef QMULTIPLY
#define QMULTIPLY

#include <stdio.h>
#include "constants.h"
#include "structs.h"

Quat randomQuat( void );
AxisAngle QuatToAxisAngle( const Quat& q );
Quat AxisAngleToQuat( const AxisAngle& aa );
Quat raaToQuat( const Real raa[3], ConstReal angle );
Quat raaDoubleToQuat( const double raa[3], const double angle );
int mkUnitQuat( Quat *const q ); // normalize in place, return 0 if error
Quat normQuat( Quat q ); // return normalized copy of quat q
AxisAngle normAxisAngle( const  AxisAngle& aa ); // return normalized axis, current angle
Real quatDifferenceToAngle( const Quat& ql, const Quat& qr );
Real quatDifferenceToAngleDeg( const Quat& ql, const Quat& qr );
Quat conjugate( const Quat& q );
Quat inverse( const Quat& q );
Quat slerp( const Quat& qa, const Quat& qb, ConstDouble t );
Quat slerp0( const Quat& qa, const Quat& qb, ConstDouble t );
Quat slerp1( const Quat& qa, const Quat& qb, ConstDouble t );
Quat axisRadianToQuat( ConstReal ax, ConstReal ay, ConstReal az, ConstReal angle );
Quat axisDegreeToQuat( ConstReal ax, ConstReal ay, ConstReal az, ConstReal angle );
Quat quatComponentsToQuat( ConstReal qx, ConstReal qy, ConstReal qz, ConstReal qw );

void qmultiply( Quat *const q, register const Quat *const ql, register const Quat *const qr );
void qconjmultiply( Quat *const q, register const Quat *const ql, register const Quat *const qr );
void printQuat_q( FILE *const fp, const Quat& q );
void printQuat_r( FILE *const fp, const Quat& q );
void printQuat( FILE *const fp, const Quat& q );
void debugQuat( FILE *const fp, const Quat& q, const unsigned int linenumber, const char *const message );
Quat randomQuatByAmount( ConstReal amount );
void print_q_reorient_message( FILE *const logFile, const Quat& q_reorient );
//void assertQuatOK( const Quat q );
Quat identityQuat();
AxisAngle identityAxisAngle();
Real a_range_reduction( ConstReal a );
Real alerp( ConstReal a, ConstReal b, ConstReal fract );
#endif
