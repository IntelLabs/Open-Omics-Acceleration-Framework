/*

 $Id: readfield.h,v 1.9 2012/02/07 05:14:55 mp Exp $

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

#ifndef READFIELD
#define READFIELD

#include "constants.h"
#include "openfile.h"
#include "stop.h"
#include "structs.h"

void readfield( GridMapSetInfo *const info, // *ptr_map_set_info
                /* not const */ char line[LINE_LEN],
                const Clock& jobStart,
                const struct tms& tms_jobStart,
		const int outlev,
		FILE *logFile);


#endif
