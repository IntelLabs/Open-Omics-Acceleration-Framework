/*

 $Id: typedefs.h,v 1.5 2009/05/08 23:17:35 rhuey Exp $

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

#ifndef _TYPEDEFS_H
#define _TYPEDEFS_H

/******************************************************************************
 *      Name: typedefs.h                                                      *
 *  Function: Defines types used in Molecular Applications.                   *
 *Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
 *----------------------------------------------------------------------------*
 *    Author: Garrett Matthew Morris, The Scripps Research Institute          *
 *      Date: JAN/18/2003                                                     *
 *----------------------------------------------------------------------------*
 *    Inputs: none                                                            *
 *   Returns: nothing                                                         *
 *   Globals: none                                                            *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 01/18/03 GMM     This header added                                         *
 ******************************************************************************/


#ifdef USE_INT_AS_LONG
    typedef int  FourByteLong;
    typedef unsigned int UnsignedFourByteLong;
#else
    typedef long FourByteLong;
    typedef unsigned long UnsignedFourByteLong;
#endif

#ifdef USE_DOUBLE
    typedef double Real;
#else
    typedef float Real;
#endif

#ifdef USE_VELOCITY_ENGINE
typedef union
{
	vector float vec;
	float		 elements[4];
} Float4;
#endif

#endif
/* EOF */
