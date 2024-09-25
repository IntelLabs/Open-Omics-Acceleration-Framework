/*

 $Id: check_header_int.cc,v 1.7 2011/03/08 04:18:36 mp Exp $

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

#include <stdio.h>
#include "check_header_int.h"


extern char *programname;
extern FILE *logFile;


void check_header_int( const int i1, 
		       const int i2, 
		       const char axis, 
		       const char *const filename )

{
    if ( i1 != i2 ) { 
        char message[LINE_LEN];

	sprintf( message, "%s: Wrong number of %c grid-points in grid-map file \"%s\".\n", programname, (char)axis, filename );

        print_2x( logFile, stderr, message );

	sprintf( message, "%s: Use either %d or %d throughout!\n", programname, i1, i2 );

        print_2x( logFile, stderr, message );

        stop("Using wrong grid-map file.\n");
    }
}
/* EOF */
