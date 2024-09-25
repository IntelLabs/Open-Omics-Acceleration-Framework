/*

 $Id: autoglobal.h,v 1.20 2012/05/01 00:22:29 mp Exp $

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

#ifndef _AUTOGLOBAL
#define _AUTOGLOBAL

#include <sys/types.h>
#include <string.h>
#include <stdio.h>

#include "structs.h"

/******************************************************************************/
/*      Name: autoglobal.h                                                    */
/*  Function: Global variables for Autodock modules.                          */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett M. Morris                                               */
/*                                                                            */
/*            e-mail: garrett@scripps.edu				      */
/*                                                                            */
/*            The Scripps Research Institute                                  */
/*            Department of Molecular Biology, MB5                            */
/*            10666 North Torrey Pines Road                                   */
/*            La Jolla, CA 92037.                                             */
/*                                                                            */
/*      Date: 03/18/93                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: None.                                                           */
/*   Returns: Nothing.                                                        */
/*   Globals: programname, AutoDockHelp, command_mode,          */
/*            command_in_fp, command_out_fp, GPF, logFile.                    */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 03/18/93 GMM     Created.                                                  */
/******************************************************************************/


/*----------------------------------------------------------------------------*/
/* Globals,                                                                   */
/*----------------------------------------------------------------------------*/

char    *programname;

char    dock_param_fn[PATH_MAX];
char    grid_param_fn[PATH_MAX];

int     command_mode = FALSE;
int     debug = 0;
int	    ElecMap = 0;
int	    DesolvMap = 0;
int     ignore_errors = FALSE;
int     keepresnum = 1;
int     parse_tors_mode = FALSE;
int     write_stateFile = FALSE;

Real	idct = 1.0;
// For energy breakdown of non-bonded interactions

FILE    *command_in_fp;
FILE    *command_out_fp;
FILE    *parFile;
FILE    *GPF;
int     outlev;   // made global for debugging convenience
FILE    *logFile;   // made global for debugging convenience
FILE    *stateFile;

Unbound_Model ad4_unbound_model = Unbound_Default;

#endif /*_AUTOGLOBAL*/
/* EOF */
