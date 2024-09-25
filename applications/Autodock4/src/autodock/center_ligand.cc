/*

 $Id: center_ligand.cc,v 1.1 2012/04/05 04:58:42 mp Exp $

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

#include <stdio.h>
#include "structs.h"
#include "autocomm.h"
#include "constants.h"


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
		FILE *const logFile )
 {

 /* if "B_auto".. compute suitable ligand center (orientation pivot point) 
  *  else use lig_center as given
  *
  * set ligand center into lig_center, sInit, and ligand_S_center vectors

  * copy   crdorig [0..natom] to crdpdb, 
  *   subtracting ligand center from range [0..true_ligand_atoms]
  * report maxrad if B_report_maxrad
  *
  * this code moved from main.cc and adapted by M Pique spring 2012
  */
Real  maxrad;
Real origin[SPACE]= {0,0,0};

	if(B_auto_center_ligand) {
		// find mean of root atoms (if any)
		// This uses the last entry in tlist, which refers to the root group.
		// Logically, root should be tlist[0] but is tlist[ntor] for historical reasons
		// This root group includes only "true ligand" (not flexres) atoms
		Real sum[SPACE]={0,0,0};

		int nrootatoms=tlist[ntor][NUM_ATM_MOVED];
		for (int a=0;a<nrootatoms;a++) for(int xyz=0;xyz<SPACE;xyz++) 
		  sum[xyz] += crdorig[tlist[ntor][NUM_ATM_MOVED+1+a]][xyz];

		// use origin (0,0,0) in unlikely case of no root atoms - which shouldn't happen
		for(int d=0;d<SPACE;d++) lig_center[d] = nrootatoms>0 ? sum[d]/nrootatoms : origin[d];

		if(outlev>=LOGBASIC) fprintf(logFile, "auto-centering ligand on root atoms: %.3f %.3f %.3f\n",
		  lig_center[X], lig_center[Y], lig_center[Z]);
		}
	else if(outlev>=LOGBASIC) fprintf(logFile, "centering ligand on specified point: %.3f %.3f %.3f\n",
	  lig_center[X], lig_center[Y], lig_center[Z]);

	/* record center used as part of overall State */
	if(sInit_Center) {
		sInit_Center->x=lig_center[X];
		sInit_Center->y=lig_center[Y];
		sInit_Center->z=lig_center[Z];
	}
	if(ligand_S_Center) {
		ligand_S_Center->x = lig_center[X];
		ligand_S_Center->y = lig_center[Y];
		ligand_S_Center->z = lig_center[Z];
	}

	//fprintf(logFile, "@@ sInit_Center    =  %.3f %.3f %.3f \n", sInit_Center->x, sInit_Center->y, sInit_Center->z);
	//fprintf(logFile, "@@ ligand_S_Center =  %.3f %.3f %.3f \n",ligand_S_Center->x, ligand_S_Center->y, ligand_S_Center->z);
        /*
        **  Center the ligand's 'true ligand atoms', copying flexres atoms unchanged
        */
            if (outlev >= LOGLIGREAD) {
                pr( logFile, "Translating small molecule by:\t" );
                pr( logFile, "(%+.3f, %+.3f, %+.3f)\n\n", -lig_center[X], -lig_center[Y], -lig_center[Z]);
            }
            /*
            **  set:   crdpdb[0..true_ligand_atoms] = crdorig - lig_center
            **  set:   crdpdb[true_ligand_atoms .. natom] = crdorig 
            */
            maxrad = 0;
            for ( int i=0; i<natom; i++ ) { /*new, gmm, 6-23-1998*/
                Real r2sum=0.0;
                for (int xyz = 0;  xyz < SPACE;  xyz++) {
		    if(i<true_ligand_atoms) {
			    crdpdb[i][xyz] = crdorig[i][xyz] - lig_center[xyz];
			    r2sum += crdpdb[i][xyz] * crdpdb[i][xyz];
			    }
		    else crdpdb[i][xyz] = crdorig[i][xyz];
                } /* xyz */
                maxrad = max(maxrad,sqrt(r2sum));
            } /* i */
            if (B_report_maxrad && true_ligand_atoms>0) {
                pr( logFile, "Furthest true ligand atom from \"about\" center is %.3f Angstroms (maxrad).\n",maxrad);
            }
 }
