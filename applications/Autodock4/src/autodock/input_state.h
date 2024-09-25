/*

 $Id: input_state.h,v 1.7 2011/06/03 05:31:36 mp Exp $

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

#ifndef INPUT_STATE
#define INPUT_STATE
#include "constants.h"
#include "qmultiply.h"

int input_state( State *S,
		 FILE  *const fp, 
                 const char  *const line, 
                 const int   ntor, 
		 int   *const P_istep, 
                 Real  *const P_energy, 
		 Real  *const P_eint, 
                 char  *const P_lastmove );
#endif
