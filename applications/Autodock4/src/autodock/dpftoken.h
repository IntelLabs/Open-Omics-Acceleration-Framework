/*

 $Id: dpftoken.h,v 1.36 2012/10/30 20:15:10 mp Exp $

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
 *      Name: dpftoken.h                                                      *
 *  Function: Define tokens for parsing DPFs (docking parameter files)        *
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
 * 09/06/95 RSH     GA/SW tokens added                                        *
 * 02/28/95 GMM     This header added                                         *
 ******************************************************************************/

#ifndef DPF_TOKENS
#define DPF_TOKENS

enum DpfTokens {
  DPF_UNKNOWN                  = -1,
  DPF_NULL                     =  0,
  DPF_COMMENT                  ,
  DPF_BLANK_LINE               ,
  DPF_FLD                      ,
  DPF_MAP                      ,
  DPF_MOVE                     ,
  DPF_ABOUT                    ,
  DPF_TRAN0                    ,
  DPF_AXISANGLE0               ,
  DPF_NDIHE                    ,
  DPF_DIHE0                    ,
  DPF_TSTEP                    ,
  DPF_QSTEP                    ,
  DPF_DSTEP                    ,
  DPF_TRNRF                    ,
  DPF_QUARF                    ,
  DPF_DIHRF                    ,
  DPF_FLEX                     ,
  DPF_INTNBP_COEFFS            ,
  DPF_RT0	               ,
  DPF_RTRF                     ,
  DPF_RUNS                     ,
  DPF_CYCLES                   ,
  DPF_ACCS                     ,
  DPF_REJS                     ,
  DPF_SELECT                   ,
  DPF_OUTLEV                   ,
  DPF_RMSTOL                   ,
  DPF_TRJFRQ                   ,
  DPF_TRJBEG                   ,
  DPF_TRJEND                   ,
  DPF_TRJOUT                   ,
  DPF_TRJSEL                   ,
  DPF_EXTNRG                   ,
  DPF_CLUSTER                  ,
  DPF_CLUSALL                  ,
  DPF_RMSNOSYM                 ,
  DPF_SCHEDGEOMETRIC            ,
  DPF_SCHEDLIN                 ,
  DPF_RMSREF                   ,
  DPF_INTELEC                  ,
  DPF_SEED                     ,
  DPF_INTNBP_REQM_EPS          ,
  DPF_WATCH                    ,
  DPF_GAUSSTORCON              ,
  DPF_BARRIER                  ,
  DPF_SHOWTORPEN               ,
  DPF_HARDTORCON               ,
  DPF_E0MAX                    ,
  DPF_CHARMAP                  ,
  DPF_RAMP_VDW_REPULSION       ,
  DPF_SIMANNEAL	               ,
  DPF_GALS                     ,
  DPF_SET_GA                   ,
  DPF_SET_SW1                  ,
  DPF_SET_PSW1                 ,
  GA_pop_size                  ,
  GA_num_generations           ,
  GA_num_evals                 ,
  GA_window_size               ,
  GA_low                       ,
  GA_high	               ,
  GA_elitism                   ,
  GA_mutation_rate             ,
  GA_crossover_rate            ,
  GA_Cauchy_alpha              ,
  GA_Cauchy_beta               ,
  SW_max_its                   ,
  SW_max_succ                  ,
  SW_max_fail                  ,
  SW_rho                       ,
  SW_lb_rho                    ,
  LS_search_freq               ,
  DPF_LS                       ,
  DPF_GS                       ,
  DPF_ANALYSIS	               ,
  DPF_TORSDOF	               ,
  DPF_INVESTIGATE	       ,
  DPF_LIG_NOT_INHIB            ,
  DPF_HOLLOW_OUT               ,
  DPF_COLINY	               ,
  DPF_INCLUDE_1_4_INTERACTIONS ,
  DPF_PARAMETER_LIBRARY	       ,
  DPF_LIGAND_TYPES	       ,
  DPF_UNBOUND	               ,
  DPF_EPDB	               ,
  DPF_TERMINATION	       ,
  GA_CROSSOVER_MODE	       ,
  DPF_POPFILE                  ,
  DPF_SET_PATTERN              ,
  DPF_COMPUTE_UNBOUND_EXTENDED ,
  DPF_FLEXRES                  ,
  DPF_ELECMAP                  ,
  DPF_DESOLVMAP                ,
  DPF_UNBOUND_INTNBP_COEFFS    ,
  DPF_RMSATOMS                 ,
  DPF_CONFSAMPLER              ,
  DPF_REORIENT                 ,
  DPF_QUATERNION0              ,
  DPF_COPYRIGHT                ,
  DPF_WARRANTY                 ,
  DPF_QUAT0	               ,
  DPF_PARAMETER_VERSION        ,
  DPF_UNBOUND_MODEL            ,
  PSW_TRANS_SCALE              ,
  PSW_ROT_SCALE                ,
  PSW_TORS_SCALE               ,
  GA_PROPORTIONAL_SELECTION    ,
  GA_TOURNAMENT_SELECTION      ,
  GA_BOLTZMAN_SELECTION        ,
  PSO_W_START	                   , // the initializers are unnecessary so...
  PSO_W_END	                   ,
  PSO_TVMAX                    ,
  PSO_QVMAX                    ,
  PSO_RVMAX                    ,
  PSO_C1 		       ,
  PSO_C2 		       ,
  PSO_K                       ,
  PSO_NEIGHBORS_DYNAMIC	       ,
  PSO_NEIGHBORS_SYMMETRIC       ,
  PSO_RANDOM_BY_DIMENSION      ,
  PSO_ADAPTIVE_VELOCITY       ,
  DPF_PSO_CONSTRICTION	      , //this is the no longer supported PSO implemented in autodock so far (5/2011)
  PSO_REGENERATE_AT_LIMIT	,
  PSO_STAGE2CONSTRICTION	,
  PSO_INTERPOLATE_AS_SCALARS	,
  DPF_PSO_SSM		        ,
  DPF_PARSWARMOPT              ,  //this is the current PSO implemented in autodock (7/2011)
  GA_LINEAR_RANKING_SELECTION  ,
  DPF_RMSMODE                  ,
  DPF_SCALE_EINTERMOL          ,
  DPF_OUTPUT_POP_STATS         ,
  DPF_OUTPUT_RESNUM_AS         ,
  DPF_SMOOTH                   ,
};
#endif
/* EOF */
