/*

 $Id: coliny.h,v 1.7 2010/08/27 00:05:07 mp Exp $

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

//
// coliny.h
//
// Interface to Coliny optimizers
//

#ifndef __coliny_h
#define __coliny_h

#if defined(USING_COLINY)

#include <vector>

//
// Initialize the 'algname' coliny optimizer over 'domain'
//
void coliny_init(const char *const algname, const char *const domain, const int num_vars);

//
// Perform minimization with a given seed and initial point. Return
// summary statistics
//
void coliny_minimize(const int seed, const std::vector<double>& initpt,
				/* not const */ std::vector<double>& finalpt,
				/* unused */ const int& neval, /* unused */ const int& niters);

#endif

#endif
