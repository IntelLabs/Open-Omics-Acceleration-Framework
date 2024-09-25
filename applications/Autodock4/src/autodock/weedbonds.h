/*

 $Id: weedbonds.h,v 1.13 2012/04/05 01:39:32 mp Exp $

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

#ifndef WEEDBONDS
#define WEEDBONDS

#include "constants.h"
#include "stop.h"
#include "structs.h"

void  weedbonds( const int   natom,
                 const char  pdbaname[MAX_ATOMS][5],
                 const int   piece[MAX_ATOMS],
                 const int   ntor,
                 const int   tlist[MAX_TORS+1][MAX_ATOMS],
       /* not const */ int   nbmatrix_binary[MAX_ATOMS][MAX_ATOMS],
       /* not const */ int   *P_Nnb,
       /* not const */ NonbondParam *const nonbondlist,
      /* not const */	int   Nnb_array[3],
		 const int   true_ligand_atoms,
		 const int   type[MAX_ATOMS],
		 const int   debug,
                 const int   outlev,
		 FILE *logFile);
#endif


#ifndef PRINT_NONBONDS
#define PRINT_NONBONDS

#include "constants.h"
#include "stop.h"
#include "structs.h"

void print_nonbonds(
                const int natom,
                const char pdbaname[MAX_ATOMS][5],
                const int piece[MAX_ATOMS],
                const int ntor,
                const int tlist[MAX_TORS][MAX_ATOMS],
                /* not const */ int nbmatrix[MAX_ATOMS][MAX_ATOMS],
                const int Nnb,
                const NonbondParam *nonbondlist,
                const int type[MAX_ATOMS],
                const int outlev,
		FILE *logFile);
#endif

