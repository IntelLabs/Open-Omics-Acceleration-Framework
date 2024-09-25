// AutoDock WIN32 MinGW replacment for times()
// adapted from MinGW include file  M Pique 2012
// $Id: mingw_sys_times.h,v 1.2 2012/04/18 05:17:46 mp Exp $
#ifndef	_SYS_TIMES_H
#ifdef __cplusplus
extern "C" {
#endif
#define	_SYS_TIMES_H
#include <time.h>   // will set or define "clock_t CLK_TCK"
#include <sys/time.h>   // for MinGW struct timeval
/*  Get Process Times, P1003.1b-1993, p. 92 */
struct tms {
	clock_t	tms_utime;		/* user time */
	clock_t	tms_stime;		/* system time */
	clock_t	tms_cutime;		/* user time, children */
	clock_t	tms_cstime;		/* system time, children */
};

//#include <sys/resource.h> // defines struct rusage

struct rusage {
	struct timeval ru_utime;
	struct timeval ru_stime;
};
int getrusage(int, struct rusage *);
#define RUSAGE_SELF 0

//clock_t _EXFUN(times,(struct tms *));
 __CRT_INLINE clock_t __cdecl times( struct tms *buffer )
{
static struct rusage self;

(void) getrusage(RUSAGE_SELF, &self);


buffer->tms_utime  = (clock_t) (( (double)self.ru_utime.tv_sec + self.ru_utime.tv_usec/1000000.) * CLK_TCK);
buffer->tms_stime  =  (clock_t)(( (double)self.ru_stime.tv_sec + self.ru_stime.tv_usec/1000000.) * CLK_TCK);
buffer->tms_cutime =  (clock_t) 0; // dummy
buffer->tms_cstime =  (clock_t) 0;
return clock();
}
#ifdef __cplusplus
}
#endif
#endif	/* !_SYS_TIMES_H */
