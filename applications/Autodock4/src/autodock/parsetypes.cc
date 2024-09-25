/*

 $Id: parsetypes.cc,v 1.6 2010/08/27 00:05:08 mp Exp $

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

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "parsetypes.h"

int parsetypes(char  line[], char * words[], const int maxwords)
/*utility func for parsing types: caution: modifies "line" in place */
{
/******************************************************************************/
/*      Name: parsetypes                                                      */
/*  Function: Parse the AutoGrid types lines                                  */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 02/01/95 (1-feb-1995)                                           */
/*----------------------------------------------------------------------------*/
/*    Inputs: line, array of pointers, cut-off number of words                */
/*   Returns: integer, number of types found or -1 if too many.               */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 06/02/03 RH      Entered code.                                             */
/******************************************************************************/

    char *char_ptr = line;
    int num_types = 0;
    /*flag for first word which is always a keyword*/
    int found_keyword = 0;

    while(1) {
        /*skip spaces*/
        while(isspace(*char_ptr)){
            char_ptr++;
        }
        /*done parsing when get eol 'null' character*/
        /* could get null after a space*/
        if (*char_ptr == '\0'){
            /*return number of 'types' found*/
            return num_types;
        }
        /* the first word is the keyword not a type*/
        if(found_keyword==0){
            found_keyword++;
        } else {
            /*words is a list of indicies of beginning of 1 or 2 char types*/
            /*if there are too many types, return*/
            if(num_types >= maxwords) return -1;

            words[num_types++] = char_ptr;
        }
        /*once in a type, skip possible 2nd characters up to a space or null
         * character*/
        while(!isspace(*char_ptr) && *char_ptr!='\0'){
            char_ptr++;
        }
        /*done parsing when get eol 'null' character*/
        /* could get null after a character*/
        if(*char_ptr=='\0'){
            return num_types;
        }
        /*make each 'type' a null terminated string*/
        *char_ptr++ = '\0';
    }
}
/* EOF */
