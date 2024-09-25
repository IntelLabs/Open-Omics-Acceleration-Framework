/*

 $Id: gridmap.h,v 1.1 2011/09/28 17:15:15 rhuey Exp $

 AutoGrid 

Copyright (C) 2011 The Scripps Research Institute. All rights reserved.

 AutoGrid is a Trade Mark of The Scripps Research Institute.

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

#ifndef INITGRIDMAP
#define INITGRIDMAP
#include "typedefs.h"
#include "structs.h"
#include "autogrid.h"  // for NUM_RECEPTOR_TYPES


typedef struct mapObject {
    int    atom_type; /*corresponds to receptor numbers????*/
    int    map_index;
    int    is_covalent;
    int    is_hbonder;
    FILE   *map_fileptr;
    char   map_filename[MAX_CHARS];
    char   type[3]; /*eg HD or OA or NA or N*/
    double constant; /*this will become obsolete*/
    double energy_max;
    double energy_min;
    double energy;
    double vol_probe;
    double solpar_probe;
    /*new 6/28*/
    double Rij;
    double epsij;
    hbond_type hbond; /*hbonding character: */
    double Rij_hb;
    double epsij_hb;
    /*per receptor type parameters, ordered as in receptor_types*/
    double nbp_r[NUM_RECEPTOR_TYPES]; /*radius of energy-well minimum*/
    double nbp_eps[NUM_RECEPTOR_TYPES];/*depth of energy-well minimum*/
    int xA[NUM_RECEPTOR_TYPES]; /*generally 12*/
    int xB[NUM_RECEPTOR_TYPES]; /*6 for non-hbonders 10 for h-bonders*/
    int hbonder[NUM_RECEPTOR_TYPES];
} MapObject;


MapObject *  init_gridmap( const int num_ligand_types, 
                          const char ligand_types[][3], 
                          const int num_receptor_types, 
                          const char receptor_types[][3], 
                          const int num_maps, 
                          const int num_grid_points_per_map,
                          FILE * logFile
                          );


#endif
