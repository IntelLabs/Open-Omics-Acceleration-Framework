/*

 $Id: sort_enrg.h,v 1.6 2010/08/27 00:05:08 mp Exp $

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


#ifndef SORT_ENERGY
#define SORT_ENERGY

#include <stdio.h>
#include "constants.h"
#include "quicksort.h"

void sort_enrg( const Real econf[MAX_RUNS],
      /* not const */ int  isort[MAX_RUNS],
                const int  nconf );
#endif
