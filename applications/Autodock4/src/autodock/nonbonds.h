/*

 $Id: nonbonds.h,v 1.12 2012/02/02 02:16:47 mp Exp $

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

#ifndef NONBONDS
#define NONBONDS
#include "constants.h"
int
nonbonds( const Real crdpdb[MAX_ATOMS][SPACE],  
	/* not const */  int nbmatrix[MAX_ATOMS][MAX_ATOMS],
	const int natom, 
	const int bond_index[MAX_ATOMS],
	const int B_include_1_4_interactions,
	const int nbonds[MAX_ATOMS], // per atom
	const int bonded[MAX_ATOMS][MAX_NBONDS],
	const int debug,
	const int outlev,
	FILE *logFile);
#endif

#ifndef GETBONDS
#define GETBONDS
#include "constants.h"
int
getbonds(const Real crdpdb[MAX_ATOMS][SPACE], 
              const int from_atom,
              const int to_atom,
	      const int bond_index[MAX_ATOMS],
	      /* not const */ int nbonds[MAX_ATOMS], // per atom
              /* not const */ int bonded[MAX_ATOMS][MAX_NBONDS],
	      const int debug,
	      const int outlev,
	      FILE *logFile);
#endif

#ifndef PRINTBONDS
#define PRINTBONDS
#include "constants.h"
void printbonds(const int natom, const int nbonds[MAX_ATOMS], const int bonded[MAX_ATOMS][MAX_NBONDS], const char *message, const int B_print_all_bonds, const int outlev, FILE *logFile);
#endif

#ifndef PRINT14
#define PRINT14
#include "constants.h"
#include <stdio.h>
void print_1_4_message(Boole B_include_1_4_interactions,  Real scale_1_4,
const int outlev, FILE *logFile);
#endif
