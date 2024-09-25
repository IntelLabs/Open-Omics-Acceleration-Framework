/*

 $Id: print_rem.cc,v 1.8 2014/06/12 01:44:07 mp Exp $

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
#include "print_rem.h"


void print_rem( FILE *const logFile,
		const int Rank,
		const int NumMem,
		const int Run,
		ConstReal ref_rms)
{
    fprintf( logFile, "MODEL     %4d\n", Run );
    fprintf( logFile, "USER    Run = %d\n", Run );
    fprintf( logFile, "USER    Cluster Rank = %d\n", Rank );
    fprintf( logFile, "USER    Number of conformations in this cluster = %d\n", NumMem );
    fprintf( logFile, "USER  \n");
    fprintf( logFile, "USER    RMSD from reference structure       = %.3f A\n", ref_rms );
    fprintf( logFile, "USER  \n");
}
/* EOF */
