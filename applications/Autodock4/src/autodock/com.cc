/*

 $Id: com.cc,v 1.12 2014/06/12 01:44:07 mp Exp $

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
#include <stdlib.h>
#include "ranlib.h"
#include "structs.h"
#include "stop.h"

/* this software can be configured to provide multiple concurrent
 random number generators (e.g., 32), however, AutoDock uses only one
 at present unless compiled for OpenMP (OMP).
 This number NUMG is defined in ranlib.h M Pique 2014-02
 */

/* this software can be configured to optionally return antithetic values,
 * see function Xqanti below. Not used in AutoDock.
 * to enable, #define ANTITHETIC
 */



/* former FORTRAN COMMON block, now statically initialized as constants: */
/*
     V=20;                            W=30;
     A1W = MOD(A1**(2**W),M1)         A2W = MOD(A2**(2**W),M2)
     A1VW = MOD(A1**(2**(V+W)),M1)    A2VW = MOD(A2**(2**(V+W)),M2)
   If V or W is changed A1W, A2W, A1VW, and A2VW need to be recomputed.
    An efficient way to precompute a**(2*j) MOD m is to start with
    a and square it j times modulo m using the function MLTMOD.
*/
#define Xm1 ((FourByteLong) 2147483563L)
#define Xm2 ((FourByteLong) 2147483399L)
#define Xa1 ((FourByteLong) 40014L)
#define Xa2 ((FourByteLong) 40692L)
#define Xa1w ((FourByteLong) 1033780774L)
#define Xa2w ((FourByteLong) 1494757890L)
#define Xa1vw ((FourByteLong) 2082007225L)
#define Xa2vw ((FourByteLong) 784306273L)

/* former FORTRAN EXTERN, now static to this source file: */
static FourByteLong Xcg1[NUMG],Xcg2[NUMG],Xig1[NUMG],Xig2[NUMG],Xlg1[NUMG],Xlg2[NUMG];

#ifdef ANTITHETIC
static int Xqanti[NUMG]; /* boolean, initially zero */
#endif

static int qqssd=false; /* have seeds been set (was gsssd() in code MP */
static int curntg=0; /* global: current generator (was 'g' from gscgn() in code) MP */

void advnst(const int k)  // not used in AutoDock code
/*
**********************************************************************
     void advnst(const int k)
               ADV-a-N-ce ST-ate
     Advances the state  of  the current  generator  by 2^K values  and
     resets the initial seed to that value.
     This is  a  transcription from   Pascal to  Fortran    of  routine
     Advance_State from the paper
     L'Ecuyer, P. and  Cote, S. "Implementing  a  Random Number Package
     with  Splitting   Facilities."  ACM  Transactions  on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     k -> The generator is advanced by 2^K values
**********************************************************************
*/
{
FourByteLong ib1,ib2;


    ib1 = Xa1;
    ib2 = Xa2;
    for (int i=0; i<k; i++) {
        ib1 = mltmod(ib1,ib1,Xm1);
        ib2 = mltmod(ib2,ib2,Xm2);
    }
    setsd_t(mltmod(ib1,Xcg1[curntg],Xm1),mltmod(ib2,Xcg2[curntg],Xm2),curntg);
/*
     NOW, IB1 = A1**K AND IB2 = A2**K
*/
}
void getsd_t(FourByteLong *const iseed1,FourByteLong *const iseed2, int curntg)
/* thread-safe if curntg is current thread number
**********************************************************************
     void getsd_t(FourByteLong *iseed1,FourByteLong *iseed2)
               GET SeeD
     Returns the value of two integer seeds of the current generator
     This  is   a  transcription from  Pascal   to  Fortran  of routine
     Get_State from the paper
     L'Ecuyer, P. and  Cote,  S. "Implementing a Random Number  Package
     with   Splitting Facilities."  ACM  Transactions   on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     iseed1 <- First integer seed of generator G
     iseed2 <- Second integer seed of generator G
**********************************************************************
*/
{

    *iseed1 = Xcg1[curntg];
    *iseed2 = Xcg2[curntg];
}
FourByteLong ignlgi_t(int curntg)
/* thread-safe if curntg is thread index (0 to NUMG-1)
**********************************************************************
     FourByteLong ignlgi(void)
               GeNerate LarGe Integer
     Returns a random integer following a uniform distribution over
     (1, 2147483562) using the current generator.
     This is a transcription from Pascal to Fortran of routine
     Random from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
**********************************************************************
*/
{
FourByteLong k,s1,s2,z;
/*
     IF THE RANDOM NUMBER PACKAGE HAS NOT BEEN INITIALIZED YET, DO SO.
     IT CAN BE INITIALIZED IN ONE OF TWO WAYS : 1) THE FIRST CALL TO
     THIS ROUTINE  2) A CALL TO SETALL.
*/
    if(!qqssd) setall(1234567890L,123456789L);
/*
     Get Current Generator
*/
    s1 = Xcg1[curntg];
    s2 = Xcg2[curntg];
    k = s1/53668L;
    s1 = Xa1*(s1-k*53668L)-k*12211;
    if(s1 < 0) s1 += Xm1;
    k = s2/52774L;
    s2 = Xa2*(s2-k*52774L)-k*3791;
    if(s2 < 0) s2 += Xm2;
    Xcg1[curntg] = s1;
    Xcg2[curntg] = s2;
    z = s1-s2;
    if(z < 1) z += (Xm1-1);
#ifdef ANTITHETIC
    if(Xqanti[curntg]) z = Xm1-z;
#endif
    return z;
}

void initgn(const int g, const int isdtyp)
/*
**********************************************************************
     void initgn(const int g, const int isdtyp)
          INIT-ialize current G-e-N-erator
     Reinitializes the state of the generator g
     This is a transcription from Pascal to Fortran of routine
     Init_Generator from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     g the generator index 0..NUMG-1
     isdtyp -> The state to which the generator is to be set
          isdtyp = -1  => sets the seeds to their initial value
          isdtyp =  0  => sets the seeds to the first value of
                          the current block
          isdtyp =  1  => sets the seeds to the first value of
                          the next block
**********************************************************************
*/
{

    switch (isdtyp) {
    case -1:
	Xlg1[g] = Xig1[g];
	Xlg2[g] = Xig2[g];
	break;
    case 0:
	break;
    case 1:
	Xlg1[g] = mltmod(Xa1w,Xlg1[g],Xm1);
	Xlg2[g] = mltmod(Xa2w,Xlg2[g],Xm2);
	break;
    default:
        stop("isdtyp not in range in INITGN");
    }

    Xcg1[g] = Xlg1[g];
    Xcg2[g] = Xlg2[g];
}
void setall(const FourByteLong iseed1,const FourByteLong iseed2)
/*
**********************************************************************
     void setall(FourByteLong iseed1,FourByteLong iseed2)
               SET ALL random number generators
     Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
     initial seeds of the other generators are set accordingly, and
     all generators states are set to these seeds.
     This is a transcription from Pascal to Fortran of routine
     Set_Initial_Seed from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     iseed1 -> First of two integer seeds
     iseed2 -> Second of two integer seeds
**********************************************************************
*/
{
/*
     TELL IGNLGI, THE ACTUAL NUMBER GENERATOR, THAT THIS ROUTINE
      HAS BEEN CALLED.
*/
// MPique - omp critical probably not needed but being on the safe side 2014
#pragma omp critical
{
    qqssd=true;

    Xig1[0] = iseed1;
    Xig2[0] = iseed2;
    initgn(0, -1);
    for (int g=1; g<NUMG; g++) {
        Xig1[g] = mltmod(Xa1vw,Xig1[g-1],Xm1);
        Xig2[g] = mltmod(Xa2vw,Xig2[g-1],Xm2);
        initgn(g, -1);
    }
    }
}
#ifdef ANTITHETIC
void setant(const FourByteLong qvalue)
/*
**********************************************************************
     void setant(FourByteLong qvalue)
               SET ANTithetic
     Sets whether the current generator produces antithetic values.  If
     X   is  the value  normally returned  from  a uniform [0,1] random
     number generator then 1  - X is the antithetic  value. If X is the
     value  normally  returned  from a   uniform  [0,N]  random  number
     generator then N - 1 - X is the antithetic value.
     All generators are initialized to NOT generate antithetic values.
     This is a transcription from Pascal to Fortran of routine
     Set_Antithetic from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     qvalue -> nonzero if generator G is to generating antithetic
                    values, otherwise zero
**********************************************************************
*/
{

    Xqanti[curntg] = qvalue;
}
#endif
void setsd_t(const FourByteLong iseed1,const FourByteLong iseed2, int curntg)
/* thread-safe if curntg is current thread number
**********************************************************************
     void setsd_t(FourByteLong iseed1,FourByteLong iseed2,curntg)
               SET S-ee-D of current generator
     Resets the initial  seed of  the current  generator to  ISEED1 and
     ISEED2. The seeds of the other generators remain unchanged.
     This is a transcription from Pascal to Fortran of routine
     Set_Seed from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     iseed1 -> First integer seed
     iseed2 -> Second integer seed
**********************************************************************
*/
{

    if(iseed1==0 || iseed2==0) {
	    char errmsg[200];
            sprintf(errmsg, 
	    "BUG: random number seed is zero in setsd_t(iseed1=%ld, iseed2=%ld, g=%d)",
 		(long) iseed1, (long) iseed2, curntg);
	    stop(errmsg);
	    }
    Xig1[curntg] = iseed1;
    Xig2[curntg] = iseed2;
    initgn(curntg, -1);
}
int gscgn(const int g)
/*
**********************************************************************
     int gscgn(const int g)
                         Get/Set GeNerator
     Sets the global number of the current generator curntg to g
     Returns previous value 
                              Arguments
     g <-- Number of the current random number generator (0..NUMG-1)
**********************************************************************
*/
{
int otg = curntg;
    if( g >= 0 && g < NUMG) curntg = g;
    else {
	    char errmsg[200];
            sprintf(errmsg, 
	    "BUG: Generator number %d out of range %d to %d in GSCGN",g,0,NUMG-1);
	    stop(errmsg);
        }
    return otg;
}
