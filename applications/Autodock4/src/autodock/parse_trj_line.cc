/*

 $Id: parse_trj_line.cc,v 1.5 2010/08/27 00:05:08 mp Exp $

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
#include <string.h>
#include <ctype.h>
#include "parse_trj_line.h"
#include "trjtokens.h"


int parse_trj_line( const char line[LINE_LEN] )

/******************************************************************************/
/*      Name: parse_trj_line                                                  */
/*  Function: Parse the trajectory file line                                  */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 26/01/93                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: line                                                            */
/*   Returns: integer token describing the keyword found.                     */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 26/01/93 GMM     Entered code.                                             */
/******************************************************************************/

{
    int l, i, token = -1 ;	       /* return -1 if nothing is recognized. */
    char c[LINE_LEN];

    l = strindex(line, " ");

    if (l == -1)
        l = strlen(line);
    
    for (i=0; i<l; i++)
        c[i] = (char)tolower( (int)line[i] );

    if ((c[0]=='\n')||(c[0]=='\0')) {
        token = TRJ_NULL;
    } else if (strncmp(c,"ntorsions",9)==0) {
        token = TRJ_NTOR;
    } else if (strncmp(c,"run",3)==0) {
        token = TRJ_RUN;
    } else if (strncmp(c,"cycle",5)==0) {
        token = TRJ_CYCLE;
    } else if (strncmp(c,"state",5)==0) {
        token = TRJ_STATE;
    } else if (strncmp(c,"temp",4)==0) {
        token = TRJ_TEMP;
    }

    return(token);
}
/* EOF */
