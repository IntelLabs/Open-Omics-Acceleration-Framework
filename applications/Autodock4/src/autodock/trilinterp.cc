/*

 $Id: trilinterp.cc,v 1.27 2014/06/12 01:44:08 mp Exp $

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
#include "trilinterp.h"

/* linear interpolation from l (when a=0) to h (when a=1)*/
/* (equal to (a*h)+((1-a)*l) )*/
#define LERP(a,l,h)	((l)+(((h)-(l))*(a)))

extern int ElecMap;
extern int DesolvMap;

#ifdef DEBUG
#include <stdio.h>
extern FILE *logFile; // DEBUG only
#endif

Real trilinterp( 

 const int first_atom, // loop begins at this atom  for (i=first_atom;
 const int last_atom, // loop ends at this atom - 1       i<last_atom; i++)
 const Real tcoord[MAX_ATOMS][SPACE], // temporary coordinates
 const Real charge[MAX_ATOMS], // partial atomic charges
 const Real abs_charge[MAX_ATOMS], // absolute magnitude of partial charges
 const int   type[MAX_ATOMS], // atom type of each atom
 #include "map_declare.h"
 const GridMapSetInfo *const info, // info->lo[X],info->lo[Y],info->lo[Z],    minimum coordinates in x,y,z
 const int ignore_inter[MAX_ATOMS], // array of booleans, says to ignore computation intermolecular energies per atom
    EnergyComponent	peratomE[MAX_ATOMS],        // output if not NULL - intermolecular energies
    EnergyComponent	*p_totalE,        // output if not NULL - total energy components
 /* not const */ EnergyComponent *energy_component // set if not NULL - breakdown by elec/vdW_Hb/desolv

 )

/******************************************************************************/
/*      Name: trilinterp                                                      */
/*  Function: Trilinear interpolation of interaction energies from map[]      */
/*            using the coordinates in tcoord[].                              */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*   Authors: Garrett M. Morris, TSRI, Accelerated C version 2.2              */
/*            David Goodsell, UCLA, Original FORTRAN version 1.0              */
/*      Date: 10/06/94                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: tcoord, charge, type, total_atoms, map, inv_spacing, lo         */
/*   Returns: total energy                                                    */
/*   Globals: MAX_ATOMS, SPACE, MAX_ATOMS, MAX_GRID_PTS, MAX_MAPS.            */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 05/05/91 GMM     Translated into C.                                        */
/* 01/03/94 GMM     Optimized code by examining 'cc -S trilinterp.c' output.  */
/* 10/06/94 GMM     Optional 10% gain in speed, using nearest point, not      */
/*                  trilinear interpolation. Compile with -DMINPOINT flag.    */
/******************************************************************************/

{
    double elec_total=0, emap_total=0, dmap_total=0, vdW_Hb_total=0, desolv_total=0;
    register int i;               /* i-th atom */
    static EnergyComponent ECzero; // const

    for (i=first_atom; i<last_atom; i++) {
        register double e, m, d; 
        register double u,   v,   w;
        float p0u, p0v, p0w;
        float p1u, p1v, p1w;
        register int AtomType;        /* atom type */
        register int u0,  v0,  w0;
        register int u1,  v1,  w1;
	register double x,y,z;

        if (ignore_inter[i]) {
            if (peratomE != NULL) peratomE[i] = ECzero;
            continue;
        }

	x = tcoord[i][X];
	y = tcoord[i][Y];
	z = tcoord[i][Z];
	if (is_out_grid_info(x,y,z)) {
                double epenalty;
                x -= info->center[X];
                y -= info->center[Y];
                z -= info->center[Z];
                // sqhypotenuse(x,y,z) is the square of the distance from grid's centre to atom
//#define GRIDEHACK
#ifdef GRIDEHACK
// TODO note GRIDEHACK
		epenalty = hypotenuse(x,y,z);
#else
                epenalty = sqhypotenuse(x,y,z) * ENERGYPENALTY;
#endif
                if (peratomE != NULL) peratomE[i].elec = peratomE[i].vdW_Hb = peratomE[i].total = epenalty;
                elec_total += epenalty;
                emap_total += epenalty;
                continue;
        }

        AtomType = type[i];

	/* MP: note u0 is < u1, weight for u==u0 value is p1u (v,w same)
	 *
	 *    u0 ............. u1
	 *    |                |
	 *  p0u=0 increases  p0u=1
	 *  p1u=1 decreases  p1u=0
	 *  
	 */
        u1  = (u0 = (int) (u = ((double)tcoord[i][X]-(double)info->lo[X]) * (double)info->inv_spacing)) + 1;
        p1u = 1.0L - (p0u = u - (double) u0);

        v1  = (v0 = (int) (v = ((double)tcoord[i][Y]-(double)info->lo[Y]) * (double)info->inv_spacing)) + 1;
        p1v = 1.0L - (p0v = v - (double) v0);

        w1  = (w0 = (int) (w = ((double)tcoord[i][Z]-(double)info->lo[Z]) * (double)info->inv_spacing)) + 1;
        p1w = 1.0L - (p0w = w - (double) w0);

#ifdef MINPOINT
        register int ix,iy,iz;                      /*MINPOINT*/
        ix = (p0u < p1u)? u0 : u1;				    /*MINPOINT*/
        iy = (p0v < p1v)? v0 : v1;				    /*MINPOINT*/
        iz = (p0w < p1w)? w0 : w1;				    /*MINPOINT*/

#ifdef MAPSUBSCRIPT
        e = map[iz][iy][ix][ElecMap];               /*MINPOINT*/
        m = map[iz][iy][ix][AtomType]; 	            /*MINPOINT*/
        d = map[iz][iy][ix][DesolvMap]; 	        /*MINPOINT*/
#else
	e = GetMap(map, info, iz, iy, ix, ElecMap); // MINPOINT
	m = GetMap(map, info, iz, iy, ix, AtomType); // MINPOINT
	d = GetMap(map, info, iz, iy, ix, DesolvMap); // MINPOINT

#endif
#else
        e = m = d = 0.0L;
#ifdef MAPSUBSCRIPT
        e += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][ElecMap];
        m += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][AtomType];
        d += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][DesolvMap];

        d += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][DesolvMap];
        m += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][AtomType];
        e += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][ElecMap];

        e += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][ElecMap];
        m += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][AtomType];
        d += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][DesolvMap];

        d += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][DesolvMap];
        m += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][AtomType];
        e += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][ElecMap];

        e += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][ElecMap];
        m += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][AtomType];
        d += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][DesolvMap];

        d += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][DesolvMap];
        m += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][AtomType];
        e += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][ElecMap];

        e += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][ElecMap];
        m += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][AtomType];
        d += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][DesolvMap];

        d += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][DesolvMap];
        m += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][AtomType];
        e += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][ElecMap];
#else
	float pu[2] = { p1u, p0u };
	float pv[2] = { p1v, p0v };
	float pw[2] = { p1w, p0w };
	for(int w=0;w<=1;w++) for(int v=0;v<=1;v++) for(int u=0;u<=1;u++) {
	 e += pu[u]*pv[v]*pw[w]*GetMap(map,info,w0+w, v0+v, u0+u, ElecMap);
	 m += pu[u]*pv[v]*pw[w]*GetMap(map,info,w0+w, v0+v, u0+u, AtomType);
	 d += pu[u]*pv[v]*pw[w]*GetMap(map,info,w0+w, v0+v, u0+u, DesolvMap);
	 }
#endif /* not MAPSUBSCRIPT */
#endif /* not MINPOINT */

//#define PRINTATOMCOMPONENTS
 /* Hack to print each atom's terms to stdout M Pique  2011-12 */
#ifdef PRINTATOMCOMPONENTS
  fprintf(stdout, "atom ligand %3d e= %8.5f m= %8.5f d= %8.5f xyz( %8.4f %8.4f %8.4f ) \n",
      i+1, e, m, d, x, y, z);
#endif
 

        elec_total += e * charge[i];
        emap_total += m;
        dmap_total += d * abs_charge[i]; 

	if (energy_component != NULL)  {
	   desolv_total += d * abs_charge[i];
	   vdW_Hb_total += m;
	   }

        if (peratomE != NULL) {
		peratomE[i].elec = e * charge[i];
		peratomE[i].vdW_Hb = m;
		peratomE[i].desolv = d * abs_charge[i];
		peratomE[i].total =  e * charge[i]  + m + d * abs_charge[i];
		}


    } // for (i=first_atom; i<last_atom; i++)

    if (p_totalE != NULL) {
    	p_totalE->elec = elec_total;
    	p_totalE->vdW_Hb = emap_total;
    	p_totalE->desolv = desolv_total;
	p_totalE->total = elec_total + emap_total + desolv_total;
	}
    if (energy_component != NULL)  {
    	energy_component->elec = elec_total;
    	energy_component->vdW_Hb = vdW_Hb_total;
    	energy_component->desolv = desolv_total;
	energy_component->total = elec_total + vdW_Hb_total + desolv_total;
	}

    return( (Real)elec_total+emap_total+dmap_total );
}

/*----------------------------------------------------------------------------*/

/* EOF */
