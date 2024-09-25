/*

 $Id: qtransform.h,v 1.12 2012/04/05 01:39:32 mp Exp $

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


#ifndef QTRANSFORM
#define QTRANSFORM

#include <stdio.h>
#include "constants.h"
#include "structs.h"
#include "writePDBQT.h"
#include "qmultiply.h"
#include "torNorVec.h"

void qtransform( const Coord& T,
	 	 const Quat&  q,
                 /* not const */ Real tcoord[MAX_ATOMS][SPACE],
		 const int   natom);

void reorient( FILE *const logFile, 
               const int true_ligand_atoms, 
               const char atomstuff[MAX_ATOMS][MAX_CHARS],
               /* not const */ Real crdpdb[MAX_ATOMS][SPACE],  // original PDB coordinates from input
               const Real charge[MAX_ATOMS],
               const int type[MAX_ATOMS],
               const ParameterEntry parameterArray[MAX_ATOM_TYPES], // input  nonbond and desolvation parameters
               const Quat& q_reorient,
               const Coord& origin,
               const int ntor,
               const int tlist[MAX_TORS+1][MAX_ATOMS],
               /* not const */ Real vt[MAX_TORS][SPACE],
               /* not const */ Molecule *ptr_ligand,
               const int debug,
	       const int outlev);
#endif
