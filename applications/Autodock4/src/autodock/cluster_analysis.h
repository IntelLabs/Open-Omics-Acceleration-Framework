/*

 $Id: cluster_analysis.h,v 1.9 2011/10/10 17:42:24 rhuey Exp $

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

#ifndef CLUSTER_ANALYSIS
#define CLUSTER_ANALYSIS
#include "constants.h"
#include "getrms.h"


int  cluster_analysis( ConstReal   clus_rms_tol, 
                       /* not const */ int   cluster[MAX_RUNS][MAX_RUNS], 
                       /* not const */ int   num_in_clus[MAX_RUNS], 
                       const int   isort[MAX_RUNS], 
                       const int   nconf, 
                       const int   natom, 
                       const int   type[MAX_ATOMS], 
                       const Real crd[MAX_RUNS][MAX_ATOMS][SPACE], 
                       const Real crdpdb[MAX_ATOMS][SPACE], 
                       const Real sml_center[SPACE], 
                       /* not const */ Real clu_rms[MAX_RUNS][MAX_RUNS], 
                       const Boole B_symmetry_flag,
                       const Boole B_unique_pair_flag,
                       /* not const */ Real ref_crds[MAX_ATOMS][SPACE],
                       const int   ref_natoms,
                       /* not const */ Real ref_rms[MAX_RUNS],
                       const Boole B_rms_heavy_atoms_only,
                       const int h_index);
#endif
