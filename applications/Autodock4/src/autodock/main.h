/*

 $Id: main.h,v 1.12 2010/08/27 00:05:07 mp Exp $

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

#ifndef MAIN
#define MAIN

#include "analysis.h"
#include "assert.h"
#include "atom_parameter_manager.h"
#include "autoglobal.h"
#include "banner.h"
#include "clmode.h"
#include "cnv_state_to_coords.h"
#include "constants.h"
#include "eintcalPrint.h"
#include "intnbtable.h"
#include "investigate.h"
#include "nbe.h"
#include "parse_dpf_line.h"
#include "parse_param_line.h"
#include "parsetypes.h"
#include "printEnergies.h"
#include "print_2x.h"
#include "printdate.h"
#include "qmultiply.h"
#include "readPDBQT.h"
#include "readfield.h"
#include "readmap.h"
#include "read_parameter_library.h"
#include "setflags.h"
#include "simanneal.h"
#include "stateLibrary.h"
#include "stop.h"
#include "strindex.h"
#include "structs.h"
#include "success.h"
#include "timesyshms.h"
#include "writePDBQT.h"

#define UNBOUND 0
#define DOCKED 1

int  main( int  argc, const char **argv);
//int main (int argc, char * const argv[], char * const envp[]);

#endif

