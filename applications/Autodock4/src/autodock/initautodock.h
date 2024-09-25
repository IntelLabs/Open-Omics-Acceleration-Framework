/*

 $Id: initautodock.h,v 1.9 2012/04/05 01:39:32 mp Exp $

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

#ifndef INITAUTODOCK
#define INITAUTODOCK

#include "constants.h"
#include "qmultiply.h"
#include "cnv_state_to_coords.h"
#include "print_2x.h"

void initautodock(
	    const char  atomstuff[MAX_ATOMS][MAX_CHARS],
	    /* not const */ Real crd[MAX_ATOMS][SPACE],
	    const Real crdpdb[MAX_ATOMS][SPACE],
	    const int   natom,
	    const int   ntor,
	    /* not const */ State *const s0,
	    const int   tlist[MAX_TORS+1][MAX_ATOMS],
	    const Real vt[MAX_TORS][SPACE],
	    const int   true_ligand_atoms,
	    const int   outlev,
            const GridMapSetInfo *const info );
#endif
