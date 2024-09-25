/*

 $Id: check_header_float.cc,v 1.7 2010/10/01 22:51:39 mp Exp $

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
#include "stop.h"
#include "check_header_float.h"


extern char *programname;
extern FILE *logFile;


void check_header_float( ConstReal   f1, ConstReal   f2, const char keyword[], const char filename[] )

{
    if ( f1 != f2 ) { 
        fprintf(logFile, "Wrong %s in grid-map file \"%s\".\n", keyword, filename);
        fprintf(stderr, "%s: Wrong %s in grid-map file \"%s\".\n", programname, keyword, filename);

        fprintf(logFile, "Use either %.3f or %.3f throughout!\n\n", f1,f2);
        fprintf(stderr, "%s: Use either %.3f or %.3f throughout!\n", programname, f1,f2);

        stop("Using wrong grid-map file.\n");
    }
}
/* EOF */
