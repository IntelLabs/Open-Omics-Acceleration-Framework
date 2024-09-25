/*

 $Id: grid.h,v 1.3 2009/05/08 23:17:34 rhuey Exp $

 AutoGrid 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

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

#ifndef GRID_MACROS
#define GRID_MACROS

/* ______________________________________________________________________________*/

/*
 * Grid Map 
 */

typedef struct      grid_map_set_info
{
    double          spacing; // uniform grid spacing in Angstroms
    double          inv_spacing; // reciprocal of the uniform grid spacing in Angstroms^-1
    char            FN_gdfld[MAX_CHARS]; // filename of the field file
    char            FN_gpf[MAX_CHARS]; // filename of the AutoGrid parameter file
    int             num_points[3]; // the actual dimensions of the grid minus 1; should be an even number
    int             num_points1[3]; // the actual dimensions of the grid; should be an odd number
    int             num_alloc[3]; // the dimensions allocated, >= num_points1; if this is a power of 2, it should be faster
    double          hi[3]; // maximum coordinates, in Angstroms
    double          lo[3]; // minimum coordinates, in Angstroms
    double          center[3]; // central coordinates, in Angstroms
    char            FN_receptor[MAX_CHARS]; // filename of the receptor used to calculate the grids
    char            atom_type_name[MAX_MAPS][3]; // array of atom type names, corresponding to the grids
    int             num_atom_types; // number of atom types
    int             num_all_maps; // number of all maps, = num_atom_types + 2
    int             num_alloc_maps; // allocated number of maps, >= num_all_maps
    
}                   GridMapSetInfo;


/*
 * Macros for accessing dynamically-allocated multidimensional arrays.
 *
 * Thanks to Mike Pique for his helpful discussions of the Dot source code!
 */

// we need space for all the atom maps (info->num_atom_types), 
// plus the electrostatic potential and the desolvation map

#define GridMapSetSize(gp) (((gp)->num_alloc[Z]) * ((gp)->num_alloc[Y]) * ((gp)->num_alloc[X]) * ((gp)->num_alloc_maps))

#define NewGridMapSet(gp) (double *) malloc( sizeof(double) * GridMapSetSize(gp) )

/* 
 * GridIndex(gp,m,x,y,z)
 *
 * Convert from a grid map index, m, and spatial x, y, z- indices 
 * to the corresponding index into the multidimensional array, 
 * where m is the fastest-varying index, as in [z][y][x][m]. 
 *
 * gp is a pointer to a GridMapSetInfo structure, typically "info"
 * m is the map index (the atom affinity, electrostatic potential or desolvation map)
 * x is the index in the X-dimension
 * y is the index in the Y-dimension
 * z is the index in the Z-dimension
 *
 * Note also: this macro could also be used for a 3D grid of vectors, where each vector is of length (gp)->num_all_maps, 
 */

#define GridIndex(gp,m,x,y,z)  (((gp)->num_alloc_maps)*(((gp)->num_alloc[X])*(((gp)->num_alloc[Y])*(z) + (y)) + (x)) + (m))


#endif
