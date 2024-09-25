/*

 $Id: read_parameter_library.cc,v 1.30 2014/02/01 05:14:53 mp Exp $

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "autocomm.h"
#include "constants.h"
#include "parameters.h"
#include "openfile.h"
#include "parse_param_line.h"
#include "atom_parameter_manager.h"
#include "default_parameters.h"
#include "stop.h"


extern char *programname;
extern int debug;


static Boole string_begins_with(const char *const a, const char *const b);
static Boole string_ends_with(const char *const a, const char *const b);
#define streq(a,b) (0==strcasecmp(a,b))

static char parameter_library[MAX_CHARS];

void read_parameter_library(
	FILE *logFile,
	const int outlev,
        const char *const FN_parameter_library,
	Linear_FE_Model *AD4
        )
{
    static ParameterEntry thisParameter;
    FILE *parameter_library_file;
    char parameter_library_line[LINE_LEN];
    int nfields;
    int param_keyword = -1;
    int int_hbond_type = 0;
    char msg[1000]; // for error messages

    pr(logFile, "Using read_parameter_library() to try to open and read \"%s\".\n\n", FN_parameter_library);

    // Open and read the parameter library
    //
    if ((parameter_library_file = ad_fopen(FN_parameter_library, "r", logFile)) == NULL) {
         sprintf(msg,"Sorry, I can't find or open %s\n", FN_parameter_library);
         stop(msg); // exits
    }

    // remember this filename for report_parameter_library()
    snprintf(parameter_library, sizeof parameter_library, "from file: \"%s\"", FN_parameter_library);
    while (fgets(parameter_library_line, sizeof(parameter_library_line), parameter_library_file) != NULL) {
        param_keyword = parse_param_line(parameter_library_line);
        //@@what was this->??param_keyword = parse_param_line("huhu");
        if (debug > 0) {
        pr(logFile, "DEBUG: parameter_library_line = %sDEBUG: param_keyword          = %d\n", parameter_library_line, param_keyword);
        }

// define convenience macro for common processing
#define process(term, string)  \
  nfields = sscanf(parameter_library_line, "%*s %lf", &AD4->term); \
  if (nfields < 1) { \
                    sprintf( msg, "%s: ERROR:  Must supply a %s coefficient as a floating point number.\n\n", programname, string); \
                    stop(msg); \
                } \
		if( outlev >= LOGETABLES ) \
                pr( logFile, "Free energy coefficient for the %s term = \t%.4lf\n\n", string, AD4->term);

        switch (param_keyword) {
            case PAR_:
            case PAR_NULL:
            case PAR_COMMENT:
                break;

            case PAR_VDW:
		process( coeff_vdW, "van der Waals term")
                break;

            case PAR_HBOND:
		process( coeff_hbond, "H-bonding term")
                break;

            case PAR_ESTAT:
		process( coeff_estat, "electrostatic term")
                break;

            case PAR_DESOLV:
		process( coeff_desolv, "desolvation term")
                break;

            case PAR_TORS:
		process( coeff_tors, "torsional term")
                break;

            case PAR_UNBOUND:
                sprintf( msg, 
		"%s: ERROR: the unbound model cannot be specified in the parameter library file.\nUse the DPF parameter 'unbound_model' instead.\n", 
		programname);
		stop(msg);
                break;

            case PAR_ATOM_PAR:
                // Read in one line of atom parameters;
                // NB: scanf doesn't try to write missing fields
                nfields = sscanf(parameter_library_line, "%*s %s %lf %lf %lf %lf %lf %lf %d %d %d %d",
                                    thisParameter.autogrid_type,
                                    &thisParameter.Rij,
                                    &thisParameter.epsij,
                                    &thisParameter.vol,
                                    &thisParameter.solpar,
                                    &thisParameter.Rij_hb,
                                    &thisParameter.epsij_hb,
                                    &int_hbond_type,
                                    &thisParameter.rec_index,
                                    &thisParameter.map_index,
                                    &thisParameter.bond_index);
                if (nfields < 2) {
                    continue; // skip any parameter_library_line without enough info
                }

                if (int_hbond_type == 0) {
                    thisParameter.hbond = NON;
                } else if (int_hbond_type == 1) {
                    thisParameter.hbond = DS;
                } else if (int_hbond_type == 2) {
                    thisParameter.hbond = D1;
                } else if (int_hbond_type == 3) {
                    thisParameter.hbond = AS;
                } else if (int_hbond_type == 4) {
                    thisParameter.hbond = A1;
                } else if (int_hbond_type == 5) {
                    thisParameter.hbond = A2;
                } else {
                    thisParameter.hbond = NON;
                }

                thisParameter.epsij    *= AD4->coeff_vdW;
                thisParameter.epsij_hb *= AD4->coeff_hbond;

                apm_enter(thisParameter.autogrid_type, thisParameter);
		if(outlev >= LOGETABLES) {
                pr(logFile, "Parameters for the atom type \"%s\" were read in from \"%s\" as follows:\n\n", thisParameter.autogrid_type, FN_parameter_library);

                if (outlev > LOGETABLES) {
		    // high precision
                    pr(logFile, "\tR-eqm = %.6f Angstrom\n",
                            thisParameter.Rij);
                    pr(logFile, "\tweighted epsilon = %.8f\n",
                            thisParameter.epsij);
                    pr(logFile, "\tAtomic fragmental volume = %.6f\n",
                            thisParameter.vol);
                    pr(logFile, "\tAtomic solvation parameter = %.8f\n",
                            thisParameter.solpar);
                    pr(logFile, "\tH-bonding R-eqm = %.6f\n",
                            thisParameter.Rij_hb);
                    pr(logFile, "\tweighted H-bonding epsilon = %.8f\n",
                            thisParameter.epsij_hb);
                    pr(logFile, "\tH-bonding type = %d,  bond index = %d\n",
                            thisParameter.hbond, thisParameter.bond_index);
                } else {
                    pr(logFile, "\tR-eqm = %.2f Angstrom,  weighted epsilon = %.3f,\n\tAt.frag.vol. = %.3f,  At.solv.par. = %.3f,\n\tHb R-eqm = %.3f,  weighted Hb epsilon = %.3f,\n\tHb type = %d,  bond index = %d\n\n",
                            thisParameter.Rij, thisParameter.epsij, thisParameter.vol, thisParameter.solpar,
                            thisParameter.Rij_hb, thisParameter.epsij_hb, thisParameter.hbond, thisParameter.bond_index );
                }
		}
                break;

            default:
                break;
        } // switch
    } // while there is another line of parameters to read in
}

void setup_parameter_library( FILE *logFile, const int outlev, 
 const char *const model_text, const Unbound_Model unbound_model,
 Linear_FE_Model *AD4)
{
    static ParameterEntry thisParameter;
    char parameter_library_line[LINE_LEN];
    int nfields;
    int param_keyword = -1;
    int int_hbond_type = 0;
    register int counter = 0;
    char msg[1000]; // for error messages

    if ( outlev >= LOGETABLES )
    pr(logFile, "Setting up parameter library with AutoDock %s values.\n\n\n", 
                 model_text);

    // Default parameters
    //
    // These are set up in "default_parameters.h"
    // and stored in the param_string_VERSION_NUM[MAX_LINES] array
    // so far we have param_string_4_0 and param_string_4_1
    // remember this choice for report_parameter_library()

    const char *const *param_string = param_string_4_1; // default
    if (unbound_model==Extended) {
        param_string=param_string_4_0;
        strncpy(parameter_library, "'extended' [AutoDock 4.0 default]", sizeof parameter_library);
    }
    else if (unbound_model==Unbound_Same_As_Bound) {
        param_string=param_string_4_1;
        strncpy(parameter_library, "'same as bound' [AutoDock 4.2 default]", sizeof parameter_library);
    }
    else {
        sprintf(msg, "BUG: cannot determine %s parameter values \n",model_text);
        stop(msg);
    }


    while ( param_string[counter] != NULL) {
	const char* const s =  param_string[counter];
        param_keyword = parse_param_line(s);

        (void)strcpy(parameter_library_line, param_string[counter]);
        counter++;
        if (debug > 0) {
            pr(logFile, "DEBUG: parameter_library_line = %sDEBUG: param_keyword          = %d\n", parameter_library_line, param_keyword);
        }

        switch (param_keyword) {
            case PAR_:
            case PAR_NULL:
            case PAR_COMMENT:
                break;

            case PAR_VDW:
		process( coeff_vdW, "van der Waals term")
                break;

            case PAR_HBOND:
		process( coeff_hbond, "H-bonding term")
                break;

            case PAR_ESTAT:
		process( coeff_estat, "electrostatic term")
                break;

            case PAR_DESOLV:
		process( coeff_desolv, "desolvation term")
                break;

            case PAR_TORS:
		process( coeff_tors, "torsional term")
                break;

            case PAR_UNBOUND:
		sprintf(msg,
                "%s: ERROR: the unbound model cannot be specified in the parameter library file.\nUse the DPF parameter 'unbound_model' instead.\n",
		programname);
		stop(msg);
                break;

            case PAR_ATOM_PAR:
                // Read in one line of atom parameters;
                // NB: scanf doesn't try to write missing fields
                nfields = sscanf(parameter_library_line, "%*s %s %lf %lf %lf %lf %lf %lf %d %d %d %d",
                                    thisParameter.autogrid_type,
                                    &thisParameter.Rij,
                                    &thisParameter.epsij,
                                    &thisParameter.vol,
                                    &thisParameter.solpar,
                                    &thisParameter.Rij_hb,
                                    &thisParameter.epsij_hb,
                                    &int_hbond_type,
                                    &thisParameter.rec_index,
                                    &thisParameter.map_index,
                                    &thisParameter.bond_index);
                if (nfields < 2) {
                    continue; // skip any parameter_library_line without enough info
                }

                if (int_hbond_type == 0) {
                    thisParameter.hbond = NON;
                } else if (int_hbond_type == 1) {
                    thisParameter.hbond = DS;
                } else if (int_hbond_type == 2) {
                    thisParameter.hbond = D1;
                } else if (int_hbond_type == 3) {
                    thisParameter.hbond = AS;
                } else if (int_hbond_type == 4) {
                    thisParameter.hbond = A1;
                } else if (int_hbond_type == 5) {
                    thisParameter.hbond = A2;
                } else {
                    thisParameter.hbond = NON;
                }

                thisParameter.epsij    *= AD4->coeff_vdW;
                thisParameter.epsij_hb *= AD4->coeff_hbond;
#ifdef ADVINA // {
// experimental 
    struct xs_vdw_lookup { char name[4]; float value;};
    static struct xs_vdw_lookup xs_vdw_lookup[] = { // adapted from AutoDock Vina atom_constants.h
    {"C", 1.9},
    {"A", 1.9},
    {"N", 1.8},
    {"NA", 1.8},
    {"NS", 1.8},
    {"O", 1.7},
    {"OA", 1.7},
    {"OS", 1.7},
    {"S", 2.0},
    {"SA", 2.0},
    {"P", 2.1},
    {"F", 1.5},
    {"Mg", 1.2}, //?Metal Donor?
    //{"X", 1.2}, //SER-OG, THR-OG or TYR-OH
    //"Cl", 1.8,
    //"Br", 2.0,
    //"I", 2.2,
    //"Met_D", 1.2,
    {"", -1 }};

                {
                int xs_vdw_i = -1;
                for (int i = 0; xs_vdw_lookup[i].value>0;i++){
                    if (streq(thisParameter.autogrid_type, xs_vdw_lookup[i].name)){
                        xs_vdw_i = i;
                        break;
                    }
                }
                if (xs_vdw_i<0) {
                    pr( logFile, "%s: WARNING: atom type %s not found in the xs_vdw_lookup table.\n\n", programname, thisParameter.autogrid_type);
                    thisParameter.xs_radius = 1.8; // default if not found; should probably be error 
                    }
                else thisParameter.xs_radius = xs_vdw_lookup[xs_vdw_i].value;
                }
#endif // } ADVINA

                apm_enter(thisParameter.autogrid_type, thisParameter);
		if(outlev >= LOGETABLES) {
                pr(logFile, "Parameters for the atom type \"%s\" were initialised with the following default values:\n\n", thisParameter.autogrid_type);


                if (outlev > LOGETABLES) {
                    pr(logFile, "\tR-eqm = %5.2f Angstrom\n\tweighted epsilon = %5.3f\n\tAtomic fragmental volume = %5.3f\n\tAtomic solvation parameter = %5.3f\n\tH-bonding R-eqm = %5.3f\n\tweighted H-bonding epsilon = %5.3f\n\tH-bonding type = %d,  bond index = %d\n\n",
                            thisParameter.Rij, thisParameter.epsij, thisParameter.vol, thisParameter.solpar,
                            thisParameter.Rij_hb, thisParameter.epsij_hb, thisParameter.hbond, thisParameter.bond_index );
                } else {
                    pr(logFile, "\tR-eqm = %.2f Angstrom,  weighted epsilon = %.3f,\n\tAt.frag.vol. = %.3f,  At.solv.par. = %.3f,\n\tHb R-eqm = %.3f,  weighted Hb epsilon = %.3f,\n\tHb type = %d,  bond index = %d\n\n",
                            thisParameter.Rij, thisParameter.epsij, thisParameter.vol, thisParameter.solpar,
                            thisParameter.Rij_hb, thisParameter.epsij_hb, thisParameter.hbond, thisParameter.bond_index );
                }
		}
                break;

            default:
                break;
        } // switch
    } // while there is another line of parameters to read in
}

const char * report_parameter_library() {
    return parameter_library;
}

static inline Boole string_begins_with(const char *const a, const char *const b) {
    // does string a begin with b  (eg   a begins with "version 4" )
    int alen=strlen(a);
    int blen=strlen(b);
    return alen>=blen && 0==strncmp(b, a, blen) ;
}

static inline Boole string_ends_with(const char *const a, const char *const b) {
    // does string a end with b  (eg   a ends with .pdb)
    int alen=strlen(a);
    int blen=strlen(b);
    return alen>=blen && 0==strcmp(b, a+alen-blen) ;
}

/* EOF */
