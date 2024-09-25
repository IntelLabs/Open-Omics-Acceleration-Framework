/*

 $Id: prClusterHist.h,v 1.8 2014/06/12 01:44:07 mp Exp $

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

#ifndef PRCLUSTERHIST
#define PRCLUSTERHIST

#include "constants.h"

void  prClusterHist(const int   ncluster,
                    const int   irunmax,
                    ConstReal clus_rms_tol,
                    const int   num_in_clu[MAX_RUNS],
                    const int   cluster[MAX_RUNS][MAX_RUNS],
                    const Real econf[MAX_RUNS],
                    const Real clu_rms[MAX_RUNS][MAX_RUNS],
		    const Real ref_rms[MAX_RUNS],
		    const int outlev,
		    FILE *const logFile);
#endif
