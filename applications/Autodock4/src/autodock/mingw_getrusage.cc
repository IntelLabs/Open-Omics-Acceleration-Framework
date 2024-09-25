// AutoDock WIN32 MinGW replacment for getrusage()
// $Id: mingw_getrusage.cc,v 1.1 2012/04/18 05:17:41 mp Exp $

#ifndef HAVE_GETRUSAGE
#ifdef __cplusplus
extern "C" {
#endif
//#include <_ansi.h>
//#include <machine/types.h>
#undef X  // used in MinGW include files and in AutoDock
#undef Y  // used in MinGW include files and in AutoDock
#undef min // used in MinGW include files and in AutoDock
#undef max // used in MinGW include files and in AutoDock
#include <time.h>   // will set or define "clock_t CLK_TCK"
#include <sys/time.h>   // for MinGW struct timeval
//#include <sys/resource.h> // defines struct rusage

struct rusage {
	struct timeval ru_utime;
	struct timeval ru_stime;
};
#define RUSAGE_SELF 0

// getrusage() by Danny Smith  http://permalink.gmane.org/gmane.comp.gnu.mingw.announce/1240 (2008)
//   no copyright claimed
// Another source for getrusage could be 
//  www.mail-archive.com/bug-gnulib@gnu.org/msg26777.html (2012)
//   M Pique 2012
#include <errno.h>
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

int __cdecl getrusage(int who, struct rusage *r_usage)
{
	 FILETIME dummy;
	 ULARGE_INTEGER KernelTime;
	 ULARGE_INTEGER UserTime;

	 if(!r_usage) {
		 errno = EFAULT;
		 return -1;
	 }
	 if(who != RUSAGE_SELF) {
		 errno= EINVAL;
		 return -1;
	 }
	 GetProcessTimes (GetCurrentProcess(), &dummy, &dummy,
			  (LPFILETIME) &KernelTime, (LPFILETIME) &UserTime);
	 r_usage->ru_stime.tv_sec = (long)  (0x7fffffff & (unsigned long) KernelTime.QuadPart/10000000);
	 r_usage->ru_utime.tv_sec =  (long) (0x7fffffff & (unsigned long) UserTime.QuadPart/10000000);
	 r_usage->ru_stime.tv_usec = (long) (0x7fffffff & (unsigned long) ((KernelTime.QuadPart%10000000)/10));
	 r_usage->ru_utime.tv_usec = (long) (0x7fffffff & (unsigned long) ((UserTime.QuadPart%10000000)/10));
	 return 0;
}
#ifdef __cplusplus
}
#endif
#endif // HAVE_GETRUSAGE
