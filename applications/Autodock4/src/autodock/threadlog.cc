/* thread-safe log file utility functions - M Pique, 2014 */

/*

 $Id: threadlog.cc,v 1.2 2014/06/12 01:44:08 mp Exp $

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

#include "threadlog.h"
#include "constants.h"
#include "stop.h"
/* include stdlib.h for "free" and unistd.h for "unlink"  */ 
/* tempnam is in <stdio.h>  */
#include <stdlib.h>
#include <unistd.h>

static char *tfilename[MAX_RUNS];
static FILE *tfileptr[MAX_RUNS];

FILE *
threadLogOpen(int j)
{
	// note that tempnam does its own malloc() and is thread-safe MPique
	tfilename[j] = tempnam(NULL, "autod");
	tfileptr[j] = fopen(tfilename[j], "w");

	if(NULL==tfileptr[j]) stop("cannot allocate or open temp log file");
	return tfileptr[j];
	}
void
threadLogClose(int j) {
	if(NULL==tfileptr[j]) stop("closing non-open temp log file");
	fclose(tfileptr[j]);
	}
void
threadLogConcat(FILE * logFile, int j) {
	int c;
	FILE * tmpfd = fopen(tfilename[j], "r");
	if(NULL==tmpfd) stop("cannot obtain new fd to concatenate log file");
	fflush(logFile);
	while( EOF != (c=getc(tmpfd)) ) putc(c, logFile);
	fflush(logFile);
	fclose(tmpfd);
	}
void
threadLogFree(int j) {
	if(NULL==tfileptr[j]) stop("freeing non-active temp log file");
#pragma omp critical
{
	unlink(tfilename[j]);
	free(tfilename[j]);
	tfileptr[j]=NULL;
}
	}
void
threadLogFreeAll(void) {
	// for emergency cleanup of the temporary files on abnormal exit
   for(int j=0;j<MAX_RUNS;j++) {
	if(NULL==tfileptr[j]) continue;
	unlink(tfilename[j]);
	//free(tfilename[j]);
	//tfileptr[j]=NULL;
   }
}
