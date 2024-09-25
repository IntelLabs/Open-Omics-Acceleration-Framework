/*

 $Id: timesys.h,v 1.7 2014/06/12 01:44:08 mp Exp $

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


#ifndef TIMESYS
#define TIMESYS

#include <sys/types.h>

// note: only timesys.h and timesyshms.h include the real sys/times.h 
#ifdef HAVE_TIMES
#include <sys/times.h>
#else
#include "mingw_sys_times.h"
#endif

#include <time.h>
#include "autocomm.h"

void  timesys( const Clock&            duration,
               const struct tms *const start,
               const struct tms *const end,
		FILE *logFile );
#endif
