/*

 $Id: parse_PDBQT_line.cc,v 1.8 2014/06/12 01:44:07 mp Exp $

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
#include "parse_PDBQT_line.h"


int parse_PDBQT_line( const char line[LINE_LEN])

/******************************************************************************/
/*      Name: parse_PDBQT_line                                                 */
/*  Function: Parse the PDBQ file line.                                       */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 11/06/93                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: line                                                            */
/*   Returns: integer token describing the keyword found.                     */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 11/06/93 GMM     Entered code.                                             */
/******************************************************************************/

{
    int l, i, token = PDBQ_UNRECOGNIZED;
    char c[LINE_LEN];

    if ((l = strindex(line, " ")) == -1) {
        l = strlen(line);
    }
    for (i=0; i<l && i<(LINE_LEN-1); i++) {
        c[i] = (char)tolower( (int)line[i] );
    }
    c[i] = '\0';

    if ((c[0]=='\n')||(c[0]=='\0')) {
		token = PDBQ_NULL;
    } else if (l >= 4) {
        if ((strncmp(c,"rema",4)==0) || (strncmp(c,"user",4)==0)) {
            token = PDBQ_REMARK;
        } else if (strncmp(c,"root",4)==0) {
            token = PDBQ_ROOT;
        } else if (strncmp(c,"endr",4)==0) {
            token = PDBQ_ENDROOT;
        } else if (strncmp(c,"atom",4)==0) {
            token = PDBQ_ATOM;
        } else if (strncmp(c,"heta",4)==0) {
            token = PDBQ_HETATM;
        } else if ((strncmp(c,"tors",4)==0) && (strncmp(c,"torsdof",7)!=0)) {
            token = PDBQ_TORS; // obsolete synonym for "branch"
        } else if ((strlen(c) >= 7) && (strncmp(c,"torsdof",7)==0)) {
            token = PDBQ_TORSDOF;
        } else if (strncmp(c,"tdof",4)==0) {
            token = PDBQ_TORSDOF;
        } else if (strncmp(c,"endt",4)==0) {
            token = PDBQ_ENDTORS;
        } else if (strncmp(c,"bran",4)==0) {
            token = PDBQ_BRANCH;
        } else if (strncmp(c,"endb",4)==0) {
            token = PDBQ_ENDBRANCH;
        } else if (strncmp(c,"cons",4)==0) {
            token = PDBQ_CONSTRAINT;
        } else if (strncmp(c,"begin_res",9)==0) {
            token = PDBQ_BEGIN_RES;
        } else if (strncmp(c,"end_res",7)==0) {
            token = PDBQ_END_RES;
        } else if (strncmp(c,"conect",6)==0) {
            token = PDBQ_CONECT;
        } else if (strncmp(c,"user",4)==0) {
            token = PDBQ_NULL;
        }
    } else {
		token = PDBQ_UNRECOGNIZED;
    }

    return( token );
}
/* EOF */
