/*

 $Id: quicksort.cc,v 1.7 2014/06/12 01:44:08 mp Exp $

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

#include "quicksort.h"

#ifdef DEBUG
#include <stdio.h>
#include <string.h>
#endif /* DEBUG */


/* 
\  Based on C.A.R.Hoare's original
 \ algorithm of 1962...
*/

#ifdef DEBUG
extern FILE *logFile; // DEBUG only
#endif /* DEBUG */

void quicksort( const Real e[], 
      /* not const */ int  isort[],
		const int  left,
		const int  right )

{
    int i, last;

#ifdef DEBUG
            char array[101];
#endif  /* DEBUG */

    if (left >= right) {
	return;
    }

#ifdef DEBUG
            strncpy( array, "----------------------------------------------------------------------------------------------------", (size_t)100 );
            array[100] = '\0';
            for (i=left+1; i<right; i++) array[i]=' ';
            array[left] = 'L';
            array[right] = 'R';
            fprintf( logFile, "%s\n", array );
#endif  /* DEBUG */

    swap( isort, left, (left+right)/2 );

    last = left;

    for (i = left+1; i <= right; i++) {
	if (e[isort[i]] < e[isort[left]]) {
	    swap( isort, ++last, i );
	}
    }

    swap( isort, (int)left, (int)last );

    quicksort( e, isort, left,  last-1 );
    quicksort( e, isort, last+1, right  );
}
/* EOF */
