/*

 $Id: check_size.cpp,v 1.13 2012/04/24 23:33:30 mp Exp $

 AutoGrid 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoGrid is a Trade Mark of The Scripps Research Institute.

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

#include <iostream>
#include <math.h>
#include <cstdlib>
#include "autogrid.h" // for print_error 


extern char *programname;
extern FILE *logFile;

/*----------------------------------------------------------------------------*/
int check_size(int nelements, 
	       char axischar)

/*----------------------------------------------------------------------------*/

/******************************************************************************/
/*      Name: check_size                                                      */
/*  Function: Checks that number of grid elements is valid.                   */ 
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 13/07/92                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: nelements, axischar                                             */
/*   Returns: nelements                                                       */
/*   Globals: MAX_GRID_PTS                                                    */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 04/01/93 GMM     Created for use in makefile.                              */
/******************************************************************************/

{
    int oldnelements;
    char msg[1000];

    if (nelements < 0) {
        fprintf(stderr, "\n%s: Error! Negative number of %c-grid elements!  Aborting.\n\n", programname, axischar);
        sprintf(msg, "\n%s: Error! Negative number of %c-grid elements!  Aborting.\n\n", programname, axischar);
        print_error( logFile, FATAL_ERROR, msg );  // exits

    }
    if (nelements == 0) {
        fprintf(stderr, "\n%s: Error!! 0 %c-grid elements!\n\n", programname, axischar);
        sprintf(msg, "\n%s: Error!! 0 %c-grid elements!\n\n", programname, axischar);
        print_error( logFile, FATAL_ERROR, msg );  // exits
    }
    if (nelements>MAX_GRID_PTS) {
        fprintf(stderr, "\n%s: Error! Maximum number of %c-grid elements allowed is %d.\n", programname, axischar, MAX_GRID_PTS);
        sprintf(msg, "\n%s: Error! Maximum number of %c-grid elements allowed is %d.\n", programname, axischar, MAX_GRID_PTS);
        print_error( logFile, FATAL_ERROR, msg);  // exits
    }
    oldnelements = nelements;
    nelements = (int) ((nelements/2) * 2); // N.B.: integer divide truncates remainder.
    if (oldnelements != nelements)
        fprintf(logFile, "%s: Number of grid elements must be even; %c-elements changed to: %d\n", programname, axischar, nelements);

    return nelements;
}
 
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
