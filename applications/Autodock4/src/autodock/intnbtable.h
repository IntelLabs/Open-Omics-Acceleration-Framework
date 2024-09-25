/*

 $Id: intnbtable.h,v 1.17 2012/05/01 00:22:29 mp Exp $

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

#ifndef INTNBTABLE
#define INTNBTABLE
#include "constants.h"
#include "timesys.h"
#include "structs.h"

void intnbtable(/* not const */ Boole *const P_B_havenbp,
                const int   a1,
                const int   a2,
                const GridMapSetInfo *const info,
                ConstReal   cA,
                ConstReal   cB,
                const int   xA,
                const int   xB,
		const Boole is_hbond,
                ConstReal   r_smooth,
                const Linear_FE_Model AD4,
                ConstDouble sigma,
                /* not const */ EnergyTables *const ad_tables,
                const Boole B_is_unbound_calculation,
		FILE *logFile,
		const int outlev);
#endif
