/*

 $Id: center_ligand.h,v 1.1 2012/04/05 04:58:42 mp Exp $

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

#include "autocomm.h"


void center_ligand( 
                 /* not const */ Real crdorig[MAX_ATOMS][SPACE],
	       const Boole B_auto_center_ligand, 
               const int natom,
               const int true_ligand_atoms, 
	       const int tlist[MAX_TORS+1][MAX_ATOMS], // used for auto-centering
	       const int ntor, // used for auto-centering
               /* not const */ Real crdpdb[MAX_ATOMS][SPACE],  // set
	       /* not const */ Real lig_center[SPACE], // set if "auto center", else used
	       /* not const */ Coord *sInit_Center, // set to final lig_center
               /* not const */ Coord *ligand_S_Center, // set to final lig_center
	       const Boole B_report_maxrad,
	       const int outlev,
		FILE *const logFile );
