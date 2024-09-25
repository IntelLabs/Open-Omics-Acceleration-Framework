/*

 $Id: mapping.cc,v 1.6 2014/06/12 01:44:07 mp Exp $

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

/********************************************************************
     These are the user defined functions that perform the
     mapping between Genotype and Phenotype and its inverse

				rsh 9/95
********************************************************************/
#include "support.h"

extern FILE *logFile; // DEBUG only

//  This should be made more efficient.  As it is now, we (de facto) AlwaysEval!!!!
//Phenotype Individual::mapping(void)
Individual &Individual::mapping(void)
{

#ifdef DEBUG
   (void)fprintf(logFile, "mapping.cc/Phenotype Individual::mapping(void)\n");
#endif /* DEBUG */
   //phenotyp.write sets phenotyp.evalflag to FALSE (0)
   phenotyp.write(*genotyp.vread(0), 0);
   phenotyp.write(*genotyp.vread(1), 1);
   phenotyp.write(*genotyp.vread(2), 2);
   phenotyp.write(*genotyp.vread(3), 3);
   phenotyp.write(*genotyp.vread(4), 4);
   //Possible future improvement
   //for (int i=0;i<5;i++) {
   //   if (*phenotyp.vread(i)!=*genotyp.vread(i)) phenotyp.write(*genotyp.vread(i), i);
   //}
   
   value(Normal_Eval);

   return(*this);
   //return(phenotyp);
}

//Genotype Individual::inverse_mapping(void)
Individual &Individual::inverse_mapping(void)
{

#ifdef DEBUG
   (void)fprintf(logFile, "mapping.cc/Genotype Individual::inverse_mapping(void)\n");
#endif /* DEBUG */

   genotyp.write(*phenotyp.vread(0), 0);
   genotyp.write(*phenotyp.vread(1), 1);
   genotyp.write(*phenotyp.vread(2), 2);
   genotyp.write(*phenotyp.vread(3), 3);
   genotyp.write(*phenotyp.vread(4), 4);

   return(*this);
   //return(genotyp);
}
