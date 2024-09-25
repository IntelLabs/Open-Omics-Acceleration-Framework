/*

 $Id: grid.h,v 1.5 2010/06/12 05:54:04 mp Exp $

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

#ifndef GRID_MACROS
#define GRID_MACROS

/* ______________________________________________________________________________*/

/*
 * Grid Map 
 */

typedef struct      grid_map_set_info
{
    int             num_alloc_maps; // allocated number of maps, >= num_all_maps
    int             num_alloc[3]; // the dimensions allocated, >= num_points1
    int             num_atom_types; // number of atom types
    int             num_all_maps; // number of all maps, = num_atom_types + 2
    int             num_points[3]; // the actual dimensions of the grid minus 1; should be an even number
    int             num_points1[3]; // the actual dimensions of the grid; should be an odd number
    int		    map_alloc_size; // max legal (allocated) index
    double          spacing; // uniform grid spacing in Angstroms
    double          inv_spacing; // reciprocal of the uniform grid spacing in Angstroms^-1
    double          hi[3]; // maximum coordinates, in Angstroms
    double          lo[3]; // minimum coordinates, in Angstroms
    double          center[3]; // central coordinates, in Angstroms
    char            FN_gdfld[PATH_MAX]; // filename of the field file
    char            FN_gpf[PATH_MAX]; // filename of the AutoGrid parameter file
    char            FN_receptor[PATH_MAX]; // filename of the receptor used to calculate the grids
    char            atom_type_name[MAX_MAPS][3]; // array of atom type names, corresponding to the grids
    
}                   GridMapSetInfo;


/*
 * Macros for accessing dynamically-allocated multidimensional arrays.
 *
 * Thanks to Mike Pique for his helpful discussions of the Dot source code!
 */


#define SetMap(map,info,z,y,x,m,val) map[GridIndex(info,z,y,x,m)]=(val)
#define GetMap(map,info,z,y,x,m) (map[GridIndex(info,z,y,x,m)])

/* 
 * GridIndex(info,z,y,x,m)

 *
 * Convert from a grid map index, m, and spatial x, y, z- indices 
 * to the corresponding index into the multidimensional array, 
 * where m is the fastest-varying index, as in [z][y][x][m]. 
 *
 * info is a pointer to a GridMapSetInfo structure, typically "info"
 * m is the map index (the atom affinity, electrostatic potential or desolvation map)
 * x is the index in the X-dimension
 * y is the index in the Y-dimension
 * z is the index in the Z-dimension
 *
 * Note also: this macro could also be used for a 3D grid of vectors, where each vector is of length (info)->num_all_maps, 
 */

#define GridIndex(info,z,y,x,m)  (((info)->num_alloc_maps)*(((info)->num_alloc[X])*(((info)->num_alloc[Y])*(z) + (y)) + (x)) + (m))


#endif
