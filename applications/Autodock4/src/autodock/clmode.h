/*

 $Id: clmode.h,v 1.12 2012/02/02 02:16:47 mp Exp $

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


#ifndef CLMODE
#define CLMODE
#include "constants.h"
#include "strindex.h"
#include "readPDBQT.h"
#include "get_atom_type.h"
#include "getpdbcrds.h"
#include "sort_enrg.h"
#include "cluster_analysis.h"
#include "prClusterHist.h"
#include "bestpdb.h"
#include "success.h"
#include "qmultiply.h"
#include "openfile.h"
void  clmode( const int   num_atm_maps, 
              ConstReal   clus_rms_tol, 
              const char *const hostnm, 
              const Clock& jobStart,
              const struct tms& tms_jobStart, 
              const Boole B_write_all_clusmem, 
              const char  *const clusFN, 
              const Real crdpdb[MAX_ATOMS][SPACE], 
              const Real sml_center[SPACE], 
              const Boole B_symmetry_flag,
              const Boole B_unique_pair_flag,
              const char  *const rms_ref_crds,
              const Boole B_rms_heavy_atoms_only,
              const int h_index,
	      const int outlev,
	      FILE *logFile
              );
#endif
