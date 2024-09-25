/*

 $Id: timesys.cc,v 1.12 2014/06/12 01:44:08 mp Exp $

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

#include <unistd.h>

#ifdef HAVE_CONFIG_H
#   include <config.h>
#endif

#include <stdio.h>
#include "timesys.h"

extern	Real	idct;

/*----------------------------------------------------------------------------*/

void timesys( const Clock&            duration,
              const struct tms *const start,
              const struct tms *const end, FILE *logFile)

/*----------------------------------------------------------------------------*/

{
	fprintf( logFile, "Real= %.2f,  CPU= %.2f,  System= %.2f\n",     (Real)duration * idct,
                         (Real)(end->tms_utime  - start->tms_utime) * idct,
                         (Real)(end->tms_stime  - start->tms_stime) * idct );
}
/* EOF */
