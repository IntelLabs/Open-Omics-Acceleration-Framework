/*

 $Id: gencau.cc,v 1.11 2014/06/12 01:44:07 mp Exp $

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

#include <math.h>
#include "ranlib.h"


Real rcauchy(ConstReal  alpha, ConstReal  beta)
/*
**********************************************************************
     Real rcauchy(Real alpha, Real beta)
         GENerate random deviate from a CAUchy distribution
                              Function
     Generates a single random deviate from a Cauchy distribution
     with mean, alpha, and parameter, beta.
                              Arguments
**********************************************************************
*/
{
    // scauchy2 is faster than scauchy1
    // according to testing involving looping over the scauchy variants and timing each
    return beta*scauchy2()+alpha;
}

