/*

 $Id: eintcalPrint.h,v 1.19 2012/08/17 02:25:05 mp Exp $

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

#ifndef _EINTCALPRINT
#define _EINTCALPRINT

Real  eintcalPrint( const NonbondParam * const nonbondlist,
                     const EnergyTables  *ad_energy_tables,
                     const Real tcoord[MAX_ATOMS][SPACE],
                     const int   Nnb,
		     int Nnb_array[3],
  	             GroupEnergy *group_energy,
                     const Boole B_calcIntElec,
                     const Boole B_include_1_4_interactions,
                     ConstReal  scale_1_4,
                     const Real qsp_abs_charge[MAX_ATOMS],
                     const Boole B_use_non_bond_cutoff,
                     const Boole B_have_flexible_residues, // if the receptor has flexibile residues, this will be set to TRUE
		     int natom,
		     const int type[],
		     char const atom_type_name[MAX_MAPS][3],
		     const int outlev,
		     FILE *logFile);  

#endif

/* EOF */
