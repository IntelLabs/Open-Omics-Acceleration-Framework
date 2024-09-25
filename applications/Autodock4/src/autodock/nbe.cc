/*

 $Id: nbe.cc,v 1.11 2014/07/01 23:27:17 mp Exp $

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
#include <stdio.h>
#include "nbe.h"


#ifdef NOSQRT

/*  ACCELERATED NON SQUARE-ROOTING VERSION;
     *  Look-up internal non-bond energy based on square-of-the-distance,
     *  in square Angstroms. This saves a square-root operation for each
     *  non-bonded interaction.
     */

#define        LookUpProc(i)        sqrt( index_to_SqAng( i ) )

#else

/*  SQUARE-ROOTING VERSION;
     *  Look-up internal non-bond energy based on distance,
     *  in Angstroms.
     */

#define        LookUpProc(i)        index_to_Ang( i )

#endif

void nbe( const GridMapSetInfo *const info,
          const EnergyTables *const ptr_ad_energy_tables,
          const int num_atm_maps, int outlev, FILE *const logFile)

{
 
    static int NUMPTS = 640; // MPique - arbitrary limit for output table
    register int i = 0;
    register int j = 0;
    register int k = 0;
    Real r = 0.;

    pr( logFile,"SUMMARY OF PAIRWISE-ATOMIC NON-BONDED INTERNAL ENERGIES\n" );
    pr( logFile,"________________________________________________________\n\n");
    pr( logFile,"Clamp pairwise-atomic interaction energies at: %.2f kcal/mol\n", EINTCLAMP );
 
    pr( logFile, "     \n r   Look-up " );
    for ( i = 0; i < num_atm_maps; i++) {
        for ( j = i; j < num_atm_maps; j++) {
            pr( logFile, "  E    " );
        }
    }
    pr( logFile, "\n /Ang  Index" );
    for ( i = 0; i < num_atm_maps; i++) {
        for ( j = i; j < num_atm_maps; j++) {
            pr( logFile, "  %2s,%-2s", info->atom_type_name[i], info->atom_type_name[j] );
        }
    }
    pr( logFile, "\n______ _____ " );
    for ( i = 0; i < num_atm_maps; i++) {
        for ( j = i; j < num_atm_maps; j++) {
            pr( logFile, " ______" );
        }
    }
    pr( logFile, "\n" );
    for ( k = 10;  k <= NUMPTS;  k += 10 ) {
        r = LookUpProc( k );
        pr( logFile, "%6.3f %5d ", r, k );
        for ( i = 0;  i < num_atm_maps; i++) {
            for ( j = i;  j < num_atm_maps; j++) {
		// confine to 6 digits, plus leading space
		// We need worry only about positive values being of large magnitude
                double e= ptr_ad_energy_tables->e_vdW_Hb[k][j][i];
		char * fmt;
		if(e>=10000) fmt=" %6.0f";  
		else if(e>=1000) fmt=" %6.1f";
		else fmt=" %6.2f";
                pr( logFile, fmt, ptr_ad_energy_tables->e_vdW_Hb[k][j][i] );
            } /*  j  */
        } /*  i  */
        pr( logFile, "\n" );
    } /*  k  */
    flushLog;
}
/* EOF */
