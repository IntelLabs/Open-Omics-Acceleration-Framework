/*

 $Id: openfile.h,v 1.9 2014/02/01 05:14:53 mp Exp $

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

#ifndef OPENFILE
#define OPENFILE

#include "constants.h"
#include <sys/types.h>      /*time_t time(time_t *tloc); */
#include <time.h>           /*time_t time(time_t *tloc); */
#include "timesys.h"
#include "print_2x.h"

int  openfile( const char *const filename,
               const char mode[],
               FILE  **const fp, 
	       FILE *logFile);


int openFile( const char *const filename,
              const char        mode[],
              FILE      **const fp,
              const Clock&      start,
              const struct tms& tms_start,
	      const Boole       mayExit,
	      FILE *logFile);

FILE *ad_fopen(const char *const path, const char *const mode, FILE *logFile);

#endif
