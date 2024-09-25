/*

 $Id: changeState.h,v 1.7 2010/10/01 22:51:39 mp Exp $

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

#ifndef CHANGESTATE
#define CHANGESTATE

State  changeState( const State& last,      /* ...must be a normalized quaternion! */
		    ConstReal   trnStep,
		    ConstReal   torStep,
		    const int   ntor,
		    const Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		    const int   N_con[MAX_TORS]);
#endif
