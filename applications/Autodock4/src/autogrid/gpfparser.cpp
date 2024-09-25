/*

 $Id: gpfparser.cpp,v 1.12 2013/09/16 21:09:30 mp Exp $

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

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "gpftoken.h"
#include "autogrid.h"
#include "constants.h"

int gpfparser( char line[LINE_LEN] )

/******************************************************************************/
/*      Name: gpfparser                                                       */
/*  Function: Parse the AutoGrid parameter file line                          */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 02/01/95 (1-feb-1995)                                           */
/*----------------------------------------------------------------------------*/
/*    Inputs: line                                                            */
/*   Returns: integer token describing the keyword found.                     */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 02/01/95 GMM     Entered code.                                             */
/******************************************************************************/

{
    int l, i, token = -1 ;	       /* return -1 if nothing is recognized. */
    char c[LINE_LEN];

    l = (int)strindex(line, " ");
    if (l == -1) {
        l = (int)strindex(line, "\t");
        if (l == -1) {
            l = (int)strlen(line);
	}
    }
    for(i=0; i<l; i++) {
        c[i] = (char)tolower( (int)line[i] );
    }

    if ((c[0]=='\n')||(c[0]=='\0')) {
        token = GPF_NULL;

    } else if (c[0]=='#') {
        token = GPF_COMMENT;

    } else if (equal(c,"receptor_types",14)) {
        token = GPF_RECEPTOR_TYPES;

    } else if (equal(c,"receptor",8)) {
        token = GPF_RECEPTOR;

    } else if (equal(c,"gridfld",7)) {
        token = GPF_GRIDFLD;

    } else if (equal(c,"npts",4)) {
        token = GPF_NPTS;

    } else if (equal(c,"spacing",7)) {
        token = GPF_SPACING;

    } else if (equal(c,"gridcenter",10)) {
        token = GPF_GRIDCENTER;

    } else if (equal(c,"types",5)) {
        token = GPF_LIGAND_TYPES;

    } else if (equal(c,"ligand_types",12)) {
        token = GPF_LIGAND_TYPES;


    } else if (equal(c,"map",3)) {
        token = GPF_MAP;

    } else if (equal(c,"elecmap",7)) {
        token = GPF_ELECMAP;

    } else if (equal(c,"dsolvmap",8) || equal(c,"desolvmap",9)) {
        token = GPF_DSOLVMAP;

    } else if (equal(c,"covalentmap",11)) {
        token = GPF_COVALENTMAP;

    } else if (equal(c,"nbp_coeffs",10)) {
        token = GPF_NBP_COEFFS;

    } else if (equal(c,"nbp_r_eps",9)) {
        token = GPF_NBP_R_EPS;

    } else if (equal(c,"dielectric",10)) {
        token = GPF_DIEL;

    } else if (equal(c,"qasp",4)) {
        token = GPF_QASP;

    } else if (equal(c,"fmap",4)) {
        token = GPF_FMAP;

    } else if (equal(c,"disorder_h",10)) {
        token = GPF_DISORDER;

    } else if (equal(c,"smooth",6)) {
        token = GPF_SMOOTH;

    } else if (equal(c,"sol_par",7)) {
        token = GPF_SOL_PAR;

    } else if (equal(c,"constant",8)) {
        token = GPF_CONSTANT;

    } else if (equal(c,"parameter_file",14)) {
        token = GPF_PARAM_FILE;

    } else if (equal(c,"use_vina_potential",18)) {
        token = GPF_USE_VINA_POTENTIAL;

    }
    return(token);
}
/* EOF */
