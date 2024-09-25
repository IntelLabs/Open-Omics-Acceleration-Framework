/*

 $Id: test_times.cc,v 1.10 2014/06/12 01:44:08 mp Exp $

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

typedef float Real;

#include <stdio.h>
#include <sys/types.h> // time_t time(time_t *tloc);
#include <time.h>      // time_t time(time_t *tloc);
#include <sys/times.h>
#include <unistd.h> // sysconf

// next includes autocomm.h which currently #define Clock clock_t
// clock_t is generally a long
//#include "timesyshms.h"

//#include <sys/types.h>
//#include <sys/times.h>
//#include <time.h>
//#include "autocomm.h"
//#include "printhms.h"

#define Clock clock_t

// POSIX migration:
#include <sys/resource.h>

void  timesyshms( const Clock&  duration,
                  const struct tms *const start,
                  const struct tms *const end,
		  FILE *logFile);


Real idct;


#ifdef USE_INT_AS_FOURBYTELONG
    typedef int  FourByteLong;
    typedef unsigned int UnsignedFourByteLong;
#else
    typedef long FourByteLong;
    typedef unsigned long UnsignedFourByteLong;
#endif


//int main( int argc, char **argv, char **envp );

int main( int argc, char **argv)
{
    static FourByteLong clktck = 0;
    struct tms tms_jobStart;
    struct tms tms_jobEnd;
    Clock  jobStart;
    Clock  jobEnd;
    long i=0L, j=0L;

    if (clktck == 0) {        /* fetch clock ticks per second first time */
        if ( (clktck = sysconf(_SC_CLK_TCK)) < (FourByteLong)0L) {
            (void) printf("\"sysconf(_SC_CLK_TCK)\" command failed in \"main.c\"\n");
	    return(1);
        } else {
            idct = (Real)1. / (Real)clktck;
            (void) printf("\n\nFYI:  Number of clock ticks per second = %d\nFYI:  Elapsed time per clock tick = %.3e seconds\n\n\n\n", clktck, idct);
        }
    }

    jobStart = times( &tms_jobStart );

    for (i=0; i<1e8; i++) {
      j = i;
      }

/*
** Get the time at the start of the run...
*/
    jobEnd = times( &tms_jobEnd );
    (void) printf( "\nRun completed;  time taken for this run:\n");
    timesyshms( jobEnd - jobStart, &tms_jobStart, &tms_jobEnd, logFile);


return 0;
}

    // #include <sys/types.h>
    // #include <sys/times.h>
    // #include <time.h>
    // #include <unistd.h>
    // #include "timesyshms.h"


//extern	Real	idct;

/*----------------------------------------------------------------------------*/

void timesyshms( const Clock&  duration,
		 const struct tms *const start,
		 const struct tms *const end )

/*----------------------------------------------------------------------------*/

{
    int   h,
          m;
    Real t,
	  T,
	  s;
    const Real min = 60.,
                hrs = 3600.;
 

    (void)fprintf( stdout, "Real= " );
    t = (Real)duration * idct;
    h = (int)(t/hrs);
    T = t - ((Real)h)*hrs;
    m = (int)(T/min);
    s = T - ((Real)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(stdout,       "%.2fs",       s );
        else
            (void)fprintf(stdout,    "%dm %05.2fs",    m, s );
    } else {
            (void)fprintf(stdout, "%dh %02dm %05.2fs", h, m, s );
    }

    (void)fprintf( stdout, ",  CPU= " );
    t = (Real)((end->tms_utime  - start->tms_utime) * idct);
    h = (int)(t/hrs);
    T = t - ((Real)h)*hrs;
    m = (int)(T/min);
    s = T - ((Real)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(stdout,       "%.2fs",       s );
        else
            (void)fprintf(stdout,    "%dm %05.2fs",    m, s );
    } else {
            (void)fprintf(stdout, "%dh %02dm %05.2fs", h, m, s );
    }

    (void)fprintf( stdout, ",  System= " );
    t = (Real)((end->tms_stime  - start->tms_stime) * idct);
    h = (int)(t/hrs);
    T = t - ((Real)h)*hrs;
    m = (int)(T/min);
    s = T - ((Real)m)*min;
    if (h == 0) {
        if (m == 0)
            (void)fprintf(stdout,       "%.2fs",       s );
        else
            (void)fprintf(stdout,    "%dm %05.2fs",    m, s );
    } else {
            (void)fprintf(stdout, "%dh %02dm %05.2fs", h, m, s );
    }

    (void)fprintf( stdout, "\n" );
}
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
