/*

 $Id: openfile.cc,v 1.12 2014/06/12 01:44:07 mp Exp $

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

/*

 $Id: openfile.cc,v 1.12 2014/06/12 01:44:07 mp Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* openfile.cc */

#include <stdio.h>
#include <stdlib.h>
/* the BOINC API header file */
#ifdef BOINC
#include "diagnostics.h"
#include "boinc_api.h" 
#include "filesys.h" 		// boinc_fopen(), etc... */
#endif
#include "openfile.h"
#include "stop.h"


extern char *programname;
/*----------------------------------------------------------------------------*/
/* fopen rewrite to either use BOINC api or normal system call */
FILE *ad_fopen(const char *const path, const char *const mode, FILE *logFile)
{
  FILE *filep;
#ifdef BOINC
  char resolved_name[512];
  const int rc = boinc_resolve_filename(path, resolved_name, sizeof(resolved_name));
  if (rc){
      fprintf(stderr, "BOINC_ERROR: cannot open filename.%s\n",path);
      boinc_finish(rc);    /* back to BOINC core */
    }
    // Then open the file with boinc_fopen() not just fopen()
    filep = boinc_fopen(resolved_name, mode);
#else
    filep = fopen(path,mode);
#endif
    return filep;
}

/*----------------------------------------------------------------------------*/
int openfile( const char *const filename,
	      const char mode[],
	      FILE **const fp,
	      FILE *logFile)

{
	if ( (*fp = ad_fopen(filename, mode, logFile)) == NULL ) {
		fprintf(stderr, "\n%s: I'm sorry; I can't find or open \"%s\"\n", programname, filename);
		fprintf(logFile,"\n%s: I'm sorry; I can't find or open \"%s\"\n", programname, filename);
		return( FALSE );
	}
	return( TRUE );
}

/*----------------------------------------------------------------------------*/
int openFile( const char *const filename,
	      const char        mode[],
	      FILE      **const fp,
	      const Clock&      start,
	      const struct tms& tms_start,
	      const Boole       mayExit,
	      FILE *logFile)

{
    Clock  jobEnd;
    struct tms tms_jobEnd;

    if ( (*fp = ad_fopen(filename, mode, logFile)) == NULL ) {
	char error_message[400];
	sprintf(error_message, "%s: I'm sorry; I can't find or open \"%s\"\n", programname, filename);

	jobEnd = times( &tms_jobEnd );
	timesys( jobEnd - start, &tms_start, &tms_jobEnd, logFile);
	pr_2x( logFile, stderr, error_message );

	if (mayExit) {
	    stop(error_message); // exits
	} else {
	    return( FALSE );
	}
    }
    return( TRUE );
}
/*----------------------------------------------------------------------------*/
