/*

 $Id: distdepdiel.cc,v 1.6 2011/03/08 04:18:36 mp Exp $

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

#include <math.h>
#include "distdepdiel.h"


double calc_ddd_Mehler_Solmajer( ConstDouble distance, ConstDouble approx_zero ) {
    /*____________________________________________________________________________
     * Distance-dependent dielectric ewds: Mehler and Solmajer, Prot Eng 4, 903-910.
     *____________________________________________________________________________*/
    const double lambda = 0.003627L;
    const double epsilon0 = 78.4L;
    const double A = -8.5525L;
    const double B = epsilon0 - A;
    const double rk= 7.7839L;
    const double lambda_B = -lambda * B;
        
    const double epsilon = A + B / (1.0L + rk*exp(lambda_B * distance));
    
    if (epsilon < approx_zero) {
        return 1.0L;
    }
    else {
        return epsilon;
    }
}
 /* EOF */
