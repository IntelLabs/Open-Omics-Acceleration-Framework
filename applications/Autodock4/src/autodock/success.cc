/*

 $Id: success.cc,v 1.11 2014/06/12 01:44:08 mp Exp $

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
#include <sys/types.h>
#include <time.h>
#include "success.h"
#include "timesyshms.h"


extern char *programname;

void success( const char *const hostnm,
	      const Clock& jobStart,
	      const struct tms& tms_jobStart, FILE *logFile)

{
    char message[LINE_LEN];
    Clock jobEnd;
    struct tms tms_jobEnd;

    pr( logFile, "\n" );
    pr( logFile, UnderLine );
    prStr( message, "%s: Successful Completion on \"%s\"\n\n", programname, hostnm );
    pr( logFile, "%s", message );

    jobEnd = times( &tms_jobEnd );

    timesyshms( jobEnd - jobStart, &tms_jobStart, &tms_jobEnd, logFile );

    pr( logFile, "%s", UnderLine );
}
/* EOF */
