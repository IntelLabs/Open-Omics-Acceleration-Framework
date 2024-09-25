/*

 $Id: swap.cc,v 1.6 2010/08/27 00:05:08 mp Exp $

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

#include "swap.h"

#ifdef DEBUG
#include <stdio.h>
#include <string.h>

extern FILE *logFile;
#endif /* DEBUG */

void swap ( int v[],
	    const int i, 
	    const int j )

{
    int temp;

#ifdef DEBUG
    int k;
    char array[101];
    strncpy( array, "----------------------------------------------------------------------------------------------------", (size_t)100 );
    array[100] = '\0';
    for (k=i+1; k<j; k++) array[k]=' ';
    array[i] = '<';
    array[j] = '>';
    fprintf( logFile, "%s", array );
    fprintf( logFile, " swapping %d & %d.\n", i, j);
#endif /* DEBUG */

    temp = v[i];
    v[i] = v[j];
    v[j] = temp;
}
/* EOF */
