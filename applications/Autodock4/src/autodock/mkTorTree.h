/*

 $Id: mkTorTree.h,v 1.13 2014/06/12 01:44:07 mp Exp $

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

#ifndef MKTORTREE
#define MKTORTREE

#include "constants.h"
#include "parse_PDBQT_line.h"
#include "stop.h"

void  mkTorTree(const int   atomnumber[MAX_RECORDS],
		const int   pdbatomnumber[2][MAX_RECORDS],
                const char  record[MAX_RECORDS][LINE_LEN],
                const int   nrecord,
                int   tlist[MAX_TORS+1][MAX_ATOMS],
                int   *const P_ntor,
                int   *const P_ntor_ligand,
                const char  *const smFileName,
                const char  pdbaname[MAX_ATOMS][5],
                Boole *const P_B_constrain,
                int   *const P_atomC1,
                int   *const P_atomC2,
                Real *const P_sqlower,
                Real *const P_squpper,
                int   *const P_ntorsdof,
                int   ignore_inter[MAX_ATOMS],
		int true_ligand_atoms,
		int outlev,
		FILE *logFile);
#endif
