/*

 $Id: trjtokens.h,v 1.3 2009/05/08 23:02:18 rhuey Exp $

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

/******************************************************************************
 *      Name: trjtokens.h                                                     *
 *  Function: Defines tokens for parsing AUTODOCK trajectory files.           *
 *Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
 *----------------------------------------------------------------------------*
 *    Author: Garrett Matthew Morris, The Scripps Research Institute          *
 *      Date: 02/28/1995                                                      *
 *----------------------------------------------------------------------------*
 *    Inputs: none                                                            *
 *   Returns: nothing                                                         *
 *   Globals: none                                                            *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 02/28/95 GMM     This header added                                         *
 ******************************************************************************/


#define TRJ_NULL     0      /* Trajectory-file token for '\n' or '\0'. */
#define TRJ_NTOR     1      /* Trajectory-file token for "ntorsions" */
#define TRJ_RUN      2      /* Trajectory-file token for "run". */
#define TRJ_CYCLE    3      /* Trajectory-file token for "cycle". */
#define TRJ_STATE    4      /* Trajectory-file token for "state". */
#define TRJ_TEMP     5      /* Trajectory-file token for "temp". */
