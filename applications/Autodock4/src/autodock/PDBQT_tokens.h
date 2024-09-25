/*

 $Id: PDBQT_tokens.h,v 1.3 2009/05/08 23:02:10 rhuey Exp $

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

/* pdbqtokens.h */

/******************************************************************************
 *      Name: pdbqtokens.h                                                    *
 *  Function: Defines the tokens for PDBQ parsing.                            *
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


#define PDBQ_UNRECOGNIZED    (-1)      /* PDBQ-file token for unrecognized symbols. */
#define PDBQ_NULL    0      /* PDBQ-file token for '\n' or '\0'. */
#define PDBQ_ROOT    1      /* PDBQ-file token for "ROOT" */
#define PDBQ_ENDROOT 2      /* PDBQ-file token for "ENDROOT". */
#define PDBQ_ATOM    3      /* PDBQ-file token for "ATOM". */
#define PDBQ_HETATM  4      /* PDBQ-file token for "HETATM". */
#define PDBQ_TORS    5      /* PDBQ-file token for "TORS". */
#define PDBQ_BRANCH  6      /* PDBQ-file token for "BRANCH". */
#define PDBQ_ENDBRANCH    7 /* PDBQ-file token for "ENDBRANCH". */
#define PDBQ_ENDTORS 8      /* PDBQ-file token for "ENDTORS". */
#define PDBQ_REMARK  9      /* PDBQ-file token for "REMARK". */
#define PDBQ_CONSTRAINT 10  /* PDBQ-file token for "CONSTRAINT". */
#define PDBQ_BEGIN_RES 11   /* PDBQ-file token for "BEGIN_RES". */
#define PDBQ_END_RES 12     /* PDBQ-file token for "END_RES". */
#define PDBQ_TORSDOF 14     /* PDBQ-file token for "TORSDOF" - torsional degrees of freedom. */
#define PDBQ_CONECT  15     /* PDBQ-file token for "CONECT" - PDB connectivity record. */
