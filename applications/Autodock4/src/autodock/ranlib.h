/*

 $Id: ranlib.h,v 1.14 2014/06/12 01:44:08 mp Exp $

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

/* Prototypes for all user accessible RANLIB routines: 
   note only a few are used in AutoDock 
   Most are functions in ranlib.cc, a few are from the underlying com.cc
 */

#ifndef _RANLIB_H
#define _RANLIB_H

#include "typedefs.h"

// set configuration of max number of random number generators, see com.cc
// In the OpenMP (OMP) version, this is also the max number of threads allowed
#define NUMG 8

//extern void advnst(const FourByteLong k);
extern Real genbet(ConstReal  aa,ConstReal  bb);
//extern Real genchi(ConstReal  df);
extern Real genexp(ConstReal  av);
extern Real genf(ConstReal  dfn, ConstReal  dfd);
//extern Real gengam(ConstReal  a,ConstReal  r);
extern void genmn(const Real *const parm, /* not const */ Real *const x, /* not const */ Real *const work);
extern void genmul(const FourByteLong n,const Real *const p,const FourByteLong ncat,/* not const */ FourByteLong *const ix);
//extern Real gennch(ConstReal  df,ConstReal  xnonc);
extern Real gennf(ConstReal  dfn, ConstReal  dfd, ConstReal  xnonc);
extern Real gennor(ConstReal  av,ConstReal  sd); // referenced by ls.h
extern void genprm(/* not const */FourByteLong *const iarray,const int larray);
extern Real genunf(ConstReal  low,ConstReal  high); // referenced by ls.h
//extern FourByteLong ignbin(const FourByteLong n,ConstReal  pp);
extern FourByteLong ignnbn(const FourByteLong n,ConstReal  p);
extern FourByteLong ignlgi_t(int); // in com.cc, referenced by gs.cc (get random long integer)

#ifdef _OPENMP
/* redefine random number functions to be thread-safe */
#include <omp.h>
#define ignlgi() ignlgi_t(omp_get_thread_num())
#define setsd(i,j) setsd_t((i),(j), omp_get_thread_num())
#define getsd(i,j) getsd_t((i),(j), omp_get_thread_num())
#else
#define ignlgi() ignlgi_t(0)
#define setsd(i,j) setsd_t((i),(j),0)
#define getsd(i,j) getsd_t((i),(j),0)
#endif
extern int gscgn(const int g); // in com.cc, referenced by main.cc (sets index of current generator)
extern void getsd_t(FourByteLong *const iseed1,FourByteLong *const iseed2, int thread); // in com.cc (get seeds)
extern void setsd_t(const FourByteLong iseed1,const FourByteLong iseed2, int thread); // in com.cc (set seeds)
//extern FourByteLong ignpoi(ConstReal  mu);
extern int ignuin(const int low,const int high); // referenced by gs.cc
extern FourByteLong mltmod(const FourByteLong a,const FourByteLong s,const FourByteLong m); // referenced by com.cc

/* ranf() changed to macro to avoid extra function call, M Pique 2013, see ranlib.cc */
//extern Real ranf(void); // referenced by gs.cc
/*
     4.656613057E-10 is 1/M1  M1 is set in a data statement in IGNLGI
      and is currently 2147483563. If M1 changes, change this also.
*/
#define ranf() ((Real)(ignlgi()*4.656613057E-10))

extern void setall(const FourByteLong iseed1,const FourByteLong iseed2); // referenced by main.cc
extern Real scauchy2(void); // referenced by gencau.cc

#endif
