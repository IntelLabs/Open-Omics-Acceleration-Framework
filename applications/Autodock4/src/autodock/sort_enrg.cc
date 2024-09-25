/*

 $Id: sort_enrg.cc,v 1.8 2013/05/23 20:06:02 mp Exp $

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sort_enrg.h"


void sort_enrg( const Real econf[MAX_RUNS],
      /* not const */ int  isort[MAX_RUNS],
		const int  nconf )

{
/*__________________________________________________________________________
 | Sort conformations on energy                                             |
 |__________________________________________________________________________|
 | Searches through all conformations;  puts in isort[0] the index of the   |
 | lowest energy, in isort[1] the next lowest energy's index, and so on.    |
 |__________________________________________________________________________|
 | WARNING: Fails if any 2 or more econf[] energies are equal.              |
 |__________________________________________________________________________|*/

    quicksort( econf, isort, 0, nconf-1 );
}
/* EOF */
