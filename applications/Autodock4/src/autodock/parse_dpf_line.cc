/*

 $Id: parse_dpf_line.cc,v 1.43 2014/06/12 01:44:07 mp Exp $

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
#include "parse_dpf_line.h"



int parse_dpf_line( const char line[] )

/******************************************************************************/
/*      Name: parse_dpf_line                                                  */
/*  Function: Parse the docking parameter file line                           */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 19/05/94                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: line                                                            */
/*   Returns: integer token describing the keyword found, or -1 if unknown    */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 27 May 2010 M Pique     length of keyword must match length of token    */
/* 06/09/95 RSH     Changed to an array implementation                        */
/* 19/05/94 GMM     Entered code.                                             */
/******************************************************************************/

{
    int j, i, token;
    char c[LINE_LEN];

    const struct {
       const char *const lexeme;
       const int tokenvalue;
    } tokentable[] = {{"ligand", DPF_MOVE},  
                      {"fld", DPF_FLD}, 
                      {"map", DPF_MAP}, 
                      {"move", DPF_MOVE}, 
                      {"about", DPF_ABOUT}, 
                      {"tran0", DPF_TRAN0}, 
                      {"quat0", DPF_QUAT0}, 
                      {"ndihe", DPF_NDIHE}, 
                      {"dihe0", DPF_DIHE0}, 
                      {"torsdof", DPF_TORSDOF}, 
                      {"tstep", DPF_TSTEP}, 
                      {"qstep", DPF_QSTEP}, 
                      {"dstep", DPF_DSTEP}, 
                      {"trnrf", DPF_TRNRF}, 
                      {"quarf", DPF_QUARF}, 
                      {"dihrf", DPF_DIHRF}, 
                      {"flex", DPF_FLEX}, 
                      {"intnbp_coeffs", DPF_INTNBP_COEFFS}, 
                      {"rt0", DPF_RT0}, 
                      {"rtrf", DPF_RTRF}, 
                      {"runs", DPF_RUNS}, 
                      {"cycles", DPF_CYCLES}, 
                      {"accs", DPF_ACCS}, 
                      {"rejs", DPF_REJS}, 
                      {"select", DPF_SELECT}, 
                      {"outlev", DPF_OUTLEV}, 
                      {"rmstol", DPF_RMSTOL}, 
                      {"trjfrq", DPF_TRJFRQ}, 
                      {"trjbeg", DPF_TRJBEG}, 
                      {"trjend", DPF_TRJEND}, 
                      {"trjout", DPF_TRJOUT}, 
                      {"trjsel", DPF_TRJSEL}, 
                      {"extnrg", DPF_EXTNRG}, 
                      {"cluster", DPF_CLUSTER}, 
                      {"write_all", DPF_CLUSALL}, 
                      {"write_all_cluster_members", DPF_CLUSALL}, 
                      {"charmap", DPF_CHARMAP}, 
                      {"rmsnosym", DPF_RMSNOSYM}, 
                      {"rmsref", DPF_RMSREF}, 
                      {"rmsmode", DPF_RMSMODE}, 
                      {"watch", DPF_WATCH}, 
                      {"geometric_schedule", DPF_SCHEDGEOMETRIC}, 
                      {"linear_schedule", DPF_SCHEDLIN}, 
                      {"schedule_linear", DPF_SCHEDLIN}, 
                      {"linsched", DPF_SCHEDLIN}, 
                      {"schedlin", DPF_SCHEDLIN}, 
                      {"intelec", DPF_INTELEC}, 
                      {"smooth", DPF_SMOOTH}, 
                      {"seed", DPF_SEED}, 
                      {"e0max", DPF_E0MAX}, 
                      {"simanneal", DPF_SIMANNEAL}, 
                      {"hardtorcon", DPF_HARDTORCON}, 
                      {"intnbp_r_eps", DPF_INTNBP_REQM_EPS}, 
                      {"gausstorcon", DPF_GAUSSTORCON}, 
                      {"barrier", DPF_BARRIER}, 
                      {"showtorpen", DPF_SHOWTORPEN}, 
                      {"ga_run", DPF_GALS}, 
                      {"gals_run", DPF_GALS}, 
                      {"do_gals", DPF_GALS}, 
                      {"set_ga", DPF_SET_GA}, 
                      {"set_sw1", DPF_SET_SW1}, 
                      {"set_psw1", DPF_SET_PSW1}, 
                      {"analysis", DPF_ANALYSIS}, 
                      {"ga_pop_size", GA_pop_size}, 
                      {"ga_num_generations", GA_num_generations}, 
                      {"ga_num_evals", GA_num_evals}, 
                      {"ga_window_size", GA_window_size}, 
                      {"ga_low", GA_low}, 
                      {"ga_high", GA_high}, 
                      {"ga_elitism", GA_elitism}, 
                      {"ga_mutation_rate", GA_mutation_rate}, 
                      {"ga_crossover_rate", GA_crossover_rate}, 
                      {"ga_cauchy_alpha", GA_Cauchy_alpha}, 
                      {"ga_cauchy_beta", GA_Cauchy_beta}, 
                      {"sw_max_its", SW_max_its}, 
                      {"sw_max_succ", SW_max_succ}, 
                      {"sw_max_fail", SW_max_fail}, 
                      {"sw_rho", SW_rho}, 
                      {"sw_lb_rho", SW_lb_rho}, 
                      {"psw_trans_scale", PSW_TRANS_SCALE}, 
                      {"psw_rot_scale", PSW_ROT_SCALE}, 
                      {"psw_tors_scale", PSW_TORS_SCALE}, 
                      {"do_local_only", DPF_LS}, 
                      {"ls_run", DPF_LS}, 
                      {"do_global_only", DPF_GS}, 
                      {"ga_only_run", DPF_GS}, 
                      {"ls_search_freq", LS_search_freq}, 
                      {"bin_energies_by_rmsd", DPF_INVESTIGATE}, 
                      {"investigate", DPF_INVESTIGATE}, 
              {"ligand_is_not_inhibitor", DPF_LIG_NOT_INHIB}, 
              {"include_1_4_interactions", DPF_INCLUDE_1_4_INTERACTIONS}, 
	      {"scale_eintermol", DPF_SCALE_EINTERMOL},
              {"parameter_library", DPF_PARAMETER_LIBRARY}, 
              {"parameter_file", DPF_PARAMETER_LIBRARY} 
              , {"ligand_types", DPF_LIGAND_TYPES}      
              , {"unbound", DPF_UNBOUND}      
              , {"epdb", DPF_EPDB}      
              , {"ga_termination_criterion", DPF_TERMINATION}      
              , {"ga_termination", DPF_TERMINATION}      
              , {"ga_crossover_mode", GA_CROSSOVER_MODE}      
              , {"output_pop_file", DPF_POPFILE}      
              , {"output_population_statistics", DPF_OUTPUT_POP_STATS}      
              , {"output_resnum_as", DPF_OUTPUT_RESNUM_AS}      
              , {"set_pattern", DPF_SET_PATTERN}      
              , {"compute_unbound_extended", DPF_COMPUTE_UNBOUND_EXTENDED} 
              , {"set_unbound_energy", DPF_UNBOUND}      
              , {"flexible_residues", DPF_FLEXRES} 
              , {"flexres", DPF_FLEXRES} 
              , {"elecmap", DPF_ELECMAP} 
              , {"desolvmap", DPF_DESOLVMAP} 
              , {"dsolvmap", DPF_DESOLVMAP} 
              , {"unbound_intnbp_coeffs", DPF_UNBOUND_INTNBP_COEFFS} 
              , {"rmsatoms", DPF_RMSATOMS} 
              , {"confsampler", DPF_CONFSAMPLER} 
              , {"reorient", DPF_REORIENT} 
              , {"axisangle0", DPF_AXISANGLE0} 
              , {"quaternion0", DPF_QUATERNION0} 
              , {"copyright", DPF_COPYRIGHT} 
              , {"warranty", DPF_WARRANTY} 
              , {"autodock_parameter_version", DPF_PARAMETER_VERSION} 
              , {"unbound_model", DPF_UNBOUND_MODEL} 
              , {"unbound_energy", DPF_UNBOUND} 
              , {"ga_proportional_selection", GA_PROPORTIONAL_SELECTION}
              , {"ga_linear_ranking_selection", GA_LINEAR_RANKING_SELECTION}      
              , {"ga_tournament_selection", GA_TOURNAMENT_SELECTION}      
              , {"ga_boltzman_selection", GA_BOLTZMAN_SELECTION}
              , {"pso_c1", PSO_C1}
              , {"pso_c2", PSO_C2}
              , {"pso_neighbors", PSO_K}
              , {"pso_k", PSO_K} // syns
              , {"pso_neighbors_symmetric", PSO_NEIGHBORS_SYMMETRIC}
              , {"pso_neighbors_dynamic", PSO_NEIGHBORS_DYNAMIC}
              , {"pso_w_start", PSO_W_START}
              , {"pso_w_end", PSO_W_END}
              , {"pso_tvmax", PSO_TVMAX}
              , {"pso_qvmax", PSO_QVMAX}
              , {"pso_rvmax", PSO_RVMAX}
              , {"pso_random_by_dimension", PSO_RANDOM_BY_DIMENSION}
              , {"pso_adaptive_velocity", PSO_ADAPTIVE_VELOCITY}
              , {"pso_stage2constriction", PSO_STAGE2CONSTRICTION}
              , {"pso_interpolate_as_scalars", PSO_INTERPOLATE_AS_SCALARS}
              , {"pso_regenerate_at_limit", PSO_REGENERATE_AT_LIMIT}
              , {"do_pso", DPF_PARSWARMOPT} 
              , {"do_cpso", DPF_PARSWARMOPT}  // syns

#if defined(USING_COLINY)
              , {"coliny", DPF_COLINY}  
#endif
              };

    // build guaranteed-terminated copy of the input line's first token into 'c'
    for (i=0; line[i]!='\0' && isspace(line[i]); i++) ; // skip leading blanks
    for (j=0; j<(int)sizeof(c) && line[i]!='\0' && !isspace(line[i]); i++, j++) {
        c[j] = line[i];
    }
    c[j]='\0'; // assure termination
    token = DPF_UNKNOWN;               /* return -1 if nothing is recognized. */

    /*  Recognize one character tokens  */

    if ((c[0]=='\n') || (c[0]=='\0')) {
        token = DPF_BLANK_LINE;
    } else if (c[0]=='#') {
        token = DPF_COMMENT;
    } else for (i=0;  i<(int)(sizeof(tokentable)/sizeof(*tokentable)); i++) {
    /*  Recognize token strings  */
#ifdef DEBUG
extern FILE *logFile; // DEBUG only
        /*(void)fprintf(logFile,"i = %d, tokentable[i].lexeme = %s, tokentable[i].value = %d, c = %s\n",i,tokentable[i].lexeme,tokentable[i].tokenvalue,c);*/
#endif
        // short match version: if (strncasecmp(tokentable[i].lexeme, c, j) == 0) {
        if (strcasecmp(tokentable[i].lexeme, c) == 0) {
            token = tokentable[i].tokenvalue;
            break;
        }
   }
    return(token);
}
/* EOF */
