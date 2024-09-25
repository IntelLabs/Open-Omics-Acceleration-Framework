/*

 $Id: readmap.h,v 1.14 2014/02/15 01:45:56 mp Exp $

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


#ifndef MAPC2F
#define MAPC2F
#include "constants.h"
#include "openfile.h"
#include "warn_bad_file.h"
#include "strindex.h"
#include "print_2x.h"
#include "check_header_line.h"
#include "warn_bad_file.h"
#include "check_header_float.h"
#include "check_header_int.h"
#include "timesys.h"
Real   mapc2f( const char C_mapValue ) ;
#endif

#ifndef READMAP
#define READMAP
#include "constants.h"
#include "openfile.h"
#include "warn_bad_file.h"
#include "strindex.h"
#include "print_2x.h"
#include "check_header_line.h"
#include "warn_bad_file.h"
#include "check_header_float.h"
#include "check_header_int.h"
#include "timesys.h"
#include "structs.h"

Statistics readmap( char line[LINE_LEN],
             const int outlev,

             const Clock& jobStart,
             const struct tms& tmsJobStart,

             const Boole B_charMap,

/* not const */ Boole *const P_B_HaveMap, 
             const int num_maps, 
             
             const GridMapSetInfo *const info,
             #include "map_declare.h"
             const char map_type, 
	     FILE *logFile
             );

#endif
