/*

 $Id: parameters.h,v 1.6 2010/10/23 00:05:02 rhuey Exp $

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

#ifndef _PARAMETERS_H
#define _PARAMETERS_H

/* ______________________________________________________________________________
** Parameter Dictionary */

#define MAX_NUM_AUTOGRID_TYPES 100
#define MAX_LEN_AUTOGRID_TYPE 7

enum hbond_type
{ NON, DS, D1, AS, A1, A2 };	/* hbonding character: */

typedef struct parameter_entry
{				// was "parm_info" in earlier AutoGrid 4 code
  char autogrid_type[MAX_LEN_AUTOGRID_TYPE + 1];	/* autogrid_type is a string based on babel_types assigned by PyBabel */
  double Rij;			/* Lennard-Jones equilibrium separation */
  double epsij;			/* Lennard-Jones energy well-depth */
  double xs_radius;     /* AutoDock vina xs_vdw_radius */
  double vol;			/* solvation volume */
  double solpar;		/* solvation parameter */
  hbond_type hbond;		/* hbonding character: 
				   NON: none, 
				   DS: spherical donor 
				   D1: directional donor
				   AS: spherical acceptor
				   A1: acceptor of 1 directional hbond
				   A2: acceptor of 2 directional hbonds */
  double Rij_hb;		/* 12-10 Lennard-Jones equilibrium separation */
  double epsij_hb;		/* 12-10 Lennard-Jones energy well-depth */
  int rec_index;		/* used to set up receptor atom_types */
  int map_index;		/* used to set up map atom_types */
  int bond_index;		/* used to set up bonds; corresponds to the enum in mdist.h */
} ParameterEntry;

#define is_hydrogen_type(t) ((strcmp(t,"H")==0) || (strcmp(t,"HD")==0) || (strcmp(t,"HS")==0))

#endif
/* EOF */
