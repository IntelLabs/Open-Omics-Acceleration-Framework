/*

 $Id: strindex.cc,v 1.6 2009/06/10 00:09:09 rhuey Exp $

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

#include "strindex.h"


int strindex( const char s[], const char t[] )

/******************************************************************************/
/*      Name: strindex                                                        */
/*  Function: return index of t in s, -1 if none.                             */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 11/30/92                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: s,t                                                             */
/*   Returns: i, index of t in s, -1 if not found.                            */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 11/30/92 GMM     Entered code.                                             */
/******************************************************************************/

{
    register int i,c1,c2;

    for (i=0; s[i] != '\0'; i++) {
        for (c1=i, c2=0; s[c1]==t[c2] && t[c2]!='\0'; c1++, c2++)
            ;
        if (c2>0 && t[c2]=='\0')
            return( i );
    }
    return( -1 );
}
/* EOF */
