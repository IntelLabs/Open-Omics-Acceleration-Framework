/*

 $Id: stateLibrary.h,v 1.8 2014/06/12 01:44:08 mp Exp $

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

#ifndef COPYSTATE
#define COPYSTATE

#include "structs.h"
#include "constants.h"

void initialiseState( /* not const */ State *const S );

void initialiseQuat( /* not const */ Quat *const Q );

void copyState( State *const destination,
		const State& source);

void printState( FILE *const fp,
		 /* not const */ State state, 
		 const int detail );

void writeState( /* not const */ FILE *const fp, 
		 /* not const */ State state );

int checkState( FILE *const fp, const State *const D );

Molecule copyStateToMolecule(const State *const source, /* not const */ Molecule *const mol);
#endif
