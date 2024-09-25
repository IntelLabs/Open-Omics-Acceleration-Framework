/* AutoDock subset from Mingw32 - for sysconf() and gethostname()
 * Adapted from MinGW  M Pique 2012
 * $Id: mingw_sysconf.h,v 1.3 2012/04/18 05:17:46 mp Exp $
 */
/*
 * This file is part of the Mingw32 package.
 *
 * unistd.h maps (roughly) to io.h
 */
#include <io.h>
#include <sys/locking.h>
 
#include <time.h>
#ifndef _SC_CLK_TCK
#define _SC_CLK_TCK 3
#endif
 
__CRT_INLINE long sysconf(int name)
{
if(name == _SC_CLK_TCK) {
return CLK_TCK;
}
else {
//printf("Only _SC_CLK_TCK can be used in simulated sysconf");
//exit(-1);
	return (long) -1;
}
}
 
#ifdef __MINGW32__
#define HAVE_GETHOSTNAME
#include <windows.h>
//int gethostname_mingw (char *, size_t);
 
__CRT_INLINE int  __cdecl gethostname_mingw (char *name, size_t len)
{
  DWORD dlen = len;
  return (GetComputerName (name, &dlen) ? 0 : -1);
}
#define gethostname gethostname_mingw
#endif
