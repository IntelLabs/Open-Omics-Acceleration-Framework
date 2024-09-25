/*

 $Id: writePDBQT.cc,v 1.44 2014/07/10 23:27:37 mp Exp $

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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "assert.h"
#include "writePDBQT.h"
#include "parse_PDBQT_line.h"
#include "calculateEnergies.h"

extern int keepresnum;
extern int write_stateFile;
extern FILE *stateFile;

void
writePDBQT(const int irun, const FourByteLong seed[2], 

		 const char *const smFileName, 
		 const char *const dpfFN, 
		 const Real sml_center[SPACE], 
		 /* not const */ State state, 
		 const int ntor, 
		 /* not const */ Real *const  Ptr_eintra, 
		 /* not const */ Real *const  Ptr_einter, 
		 const int natom, 
		 const char atomstuff[MAX_ATOMS][MAX_CHARS], 
		 /* not const */ Real crd[MAX_ATOMS][SPACE], 
		    EnergyComponent	peratomE[MAX_ATOMS],
		 const Real charge[MAX_ATOMS], 
		 const Real abs_charge[MAX_ATOMS], 
		 const Real qsp_abs_charge[MAX_ATOMS], 
		 const int ligand_is_inhibitor, 
		 const Real torsFreeEnergy, 
		 const Real vt[MAX_TORS][SPACE], 
		 const int tlist[MAX_TORS][MAX_ATOMS], 
		 const Real crdpdb[MAX_ATOMS][SPACE], 
		 const NonbondParam *const nonbondlist, 
         const EnergyTables *const ptr_ad_energy_tables, 
		 const int type[MAX_ATOMS],  // aka 'map_index' in 'ParameterEntry' structures
		 const int Nnb, 
		 int Nnb_array[3],
		 GroupEnergy *group_energy, 
		 const int true_ligand_atoms,
		 const Boole B_calcIntElec, 
         #include "map_declare.h"
		 const int ignore_inter[MAX_ATOMS], 
		 const Boole B_include_1_4_interactions, 
		 const Real scale_1_4, 
         const ParameterEntry parameterArray[MAX_ATOM_TYPES], 
		 const Real unbound_internal_FE, 

         const GridMapSetInfo *const info, 
         const int state_type,  // 0 means the state is unbound, 1 means the state is docked
         const char PDBQT_record[MAX_RECORDS][LINE_LEN], 
         const Boole B_use_non_bond_cutoff, 
         const Boole B_have_flexible_residues, 
         const Unbound_Model ad4_unbound_model,
	 const int outlev, 
	 FILE *logFile
         )

{
	int i = 0;
    EnergyBreakdown eb;

	EnergyComponent totalE;
	Real MaxValue = 99.99L;
	Real MinValue = -99.99L;

    char state_type_string[MAX_CHARS];
    char state_type_prefix_string[MAX_CHARS];
    char state_type_prefix_USER_string[MAX_CHARS];
    Real this_emap = 0.; // includes desolv
    Real this_elec = 0.;

    // Initialise various character strings
    if (state_type == 0) {
        strcpy(state_type_string, "UNBOUND");
        strcpy(state_type_prefix_string, "UNBOUND: ");
        strcpy(state_type_prefix_USER_string, "UNBOUND: USER    ");
    } else if (state_type == 1) {
        strcpy(state_type_string, "DOCKED");
        strcpy(state_type_prefix_string, "DOCKED: ");
        strcpy(state_type_prefix_USER_string, "DOCKED: USER    ");
    }

    // Write out the state variables

	// pass original center to printState - could be improved - MP 2010-05
	state.Center.x = sml_center[X];
	state.Center.y = sml_center[Y];
	state.Center.z = sml_center[Z];

	pr(logFile, "Detailed state: ");
	printState(logFile, state, 6); // detailed, include center, ntor, no newline
	pr(logFile, "\n");

        // "outlev" is the level of detail: >2 is high, 0 is low
	if (outlev >= LOGRUNVV) {
        pr(logFile, "QState:\t");
        printState(logFile, state, 5); // short format, as quaternion
        pr(logFile, "\n");
	}

	if (outlev >= LOGFORADT) {
	// MP I believe ADT write_models_from_states.py expects this line  2012-05
        pr(logFile, "State:\t"); // various possibly longer formats, as axis-angle
        printState(logFile, state, outlev);
        pr(logFile, "\n");
	}
        pr(logFile, "\n"); // separator before DOCKED: ...  lines

    // Convert state variables to x, y, z-coordinates
	cnv_state_to_coords( state, vt, tlist, ntor, crdpdb, crd, natom,
	  true_ligand_atoms, outlev, logFile);

    // Calculate the energy breakdown
    eb = calculateBindingEnergies( natom, ntor, unbound_internal_FE, torsFreeEnergy, B_have_flexible_residues, 
         crd, charge, abs_charge, type, map, info, 
         ignore_inter, peratomE, &totalE,
         nonbondlist, ptr_ad_energy_tables, 
	 Nnb, Nnb_array, group_energy, true_ligand_atoms, B_calcIntElec, 
         B_include_1_4_interactions, scale_1_4, qsp_abs_charge, B_use_non_bond_cutoff, ad4_unbound_model, outlev, logFile);

    // Set the total intramolecular energy (sum of intramolecular energies of ligand and of protein)
    if (ntor > 0) {
        // Add the intramolecular energy of the receptor, for the (moving, fixed) atom pairs // (2)
        if(Ptr_eintra!=NULL) *Ptr_eintra = 
	  group_energy->intra_moving_moving_lig.total + 
	   group_energy->intra_moving_moving_rec.total +
	     eb.e_intra_moving_fixed_rec;
    } else {
        if(Ptr_eintra!=NULL) *Ptr_eintra = 0.0;
    }

    // Set the total intermolecular energy
    if (state_type == 1) {
        // DOCKED
        // Set *Ptr_einter, the intermolecular energy, only for DOCKED states, not for UNBOUND states
        if(Ptr_einter!=NULL) *Ptr_einter = eb.e_inter;
    } else {
        // UNBOUND
        // "intermolecular" energy is meaningless for unbound state, so set this to zero
	static EnergyComponent zeroEC; // const
        if(Ptr_einter!=NULL) *Ptr_einter = 0.0;
        eb.e_inter = 0.0;
        totalE = zeroEC;
        eb.e_inter_moving_fixed = 0.0;
        eb.e_inter_moving_moving = 0.0;
    }

	if (outlev >= LOGFORADT ) {
		AxisAngle aa = QuatToAxisAngle( state.Q );
		// output of coordinates
        pr( logFile, "%s: MODEL     %4d\n", state_type_string, irun+1 );
        pr( logFile, "%s: USER    Run = %d\n", state_type_string, irun+1 );
        pr( logFile, "%s: USER    DPF = %s\n", state_type_string, dpfFN );
        pr( logFile, "%s: USER  \n", state_type_string );
        
	// see also main.cc and analysis.cc for similar code:
        printEnergies( &eb, state_type_prefix_USER_string, ligand_is_inhibitor, 
	totalE.vdW_Hb+totalE.desolv, totalE.elec, 
	 B_have_flexible_residues,  // next two terms are meaningful only if have flexible residues...
	 group_energy->inter_moving_moving.vdW_Hb + group_energy->inter_moving_moving.desolv,
	 group_energy->inter_moving_moving.elec,
	 ad4_unbound_model, outlev, logFile);

        // Write part of the "XML" state file
		if (write_stateFile) {
			pr(stateFile, "\n");
			pr(stateFile, "\t<run id=\"%4d\">\n", irun + 1);
			pr(stateFile, "\t\t<seed>" FBL_FMT " " FBL_FMT "</seed>\n", seed[0], seed[1]);
			pr(stateFile, "\t\t<dpf>%s</dpf>\n", dpfFN);
            printStateEnergies( &eb, state_type_prefix_USER_string, ligand_is_inhibitor, outlev, stateFile);
		} // End write state file

		(void) fprintf(logFile, "%s: USER    NEWDPF move %s\n", state_type_string, smFileName);
		(void) fprintf(logFile, "%s: USER    NEWDPF about %f %f %f\n", state_type_string, sml_center[X], sml_center[Y], sml_center[Z]);
		(void) fprintf(logFile, "%s: USER    NEWDPF tran0 %f %f %f\n", state_type_string, state.T.x, state.T.y, state.T.z);
		(void) fprintf(logFile, "%s: USER    NEWDPF quaternion0 %f %f %f %f\n", state_type_string, state.Q.x, state.Q.y, state.Q.z, state.Q.w);
		(void) fprintf(logFile, "%s: USER    NEWDPF axisangle0 %f %f %f %f\n", state_type_string, aa.nx, aa.ny, aa.nz, RadiansToDegrees(WrpRad(ModRad(aa.ang))));
		// note quat0 is deprecated, same as axis-angle
		(void) fprintf(logFile, "%s: USER    NEWDPF quat0 %f %f %f %f\n", state_type_string, aa.nx, aa.ny, aa.nz, RadiansToDegrees(WrpRad(ModRad(aa.ang))));
		if (ntor > 0) {
            // ndihe is deprecated; uses the number of torsions in the PDBQT's torsion tree
			// (void) fprintf(logFile, "%s: USER    NEWDPF ndihe %d\n", state_type_string, ntor);
			(void) fprintf(logFile, "%s: USER    NEWDPF dihe0 ", state_type_string);
			for (i = 0; i < ntor; i++) {
				(void) fprintf(logFile, "%.2f ", RadiansToDegrees(WrpRad(ModRad(state.tor[i]))));
			}
			(void) fprintf(logFile, "\n");

		}
        
        // Write remaining part of the "XML" state file
		if (write_stateFile) {
			pr(stateFile, "\t\t<move>%s</move>\n", smFileName);
			pr(stateFile, "\t\t<about>%f %f %f</about>\n", sml_center[X], sml_center[Y], sml_center[Z]);

			pr(stateFile, "\t\t<tran0>%f %f %f</tran0>\n", state.T.x, state.T.y, state.T.z);
			pr(stateFile, "\t\t<quaternion0>%f %f %f %f</quaternion0>\n", state.Q.x, state.Q.y, state.Q.z, state.Q.w);
			// quat0 is deprecated, same as axisangle0
			pr(stateFile, "\t\t<quat0>%f %f %f %f</quat0>\n", aa.nx, aa.ny, aa.nz, RadiansToDegrees(WrpRad(ModRad(aa.ang))));
			pr(stateFile, "\t\t<axisangle0>%f %f %f %f</axisangle0>\n", aa.nx, aa.ny, aa.nz, RadiansToDegrees(WrpRad(ModRad(aa.ang))));
			if (ntor > 0) {
				pr(stateFile, "\t\t<ndihe>%d</ndihe>\n", ntor);
				pr(stateFile, "\t\t<dihe0>");
				for (i = 0; i < ntor; i++) {
					(void) fprintf(stateFile, "%.2f ", RadiansToDegrees(WrpRad(ModRad(state.tor[i]))));
				}
				(void) fprintf(stateFile, "\n");
				pr(stateFile, "</dihe0>\n");
			}
			pr(stateFile, "\t</run>\n");
		} // End write state file


        (void) fprintf(logFile, "%s: USER  keepresnum = %d \n", state_type_string, keepresnum);
        (void) fprintf(logFile, "%s: USER  \n", state_type_string);

        // Count the number of non-NULL records in the PDBQT file
        int nrecord = 0;
        int r = 0;
        for (r = 0; PDBQT_record[r][0] != '\0'; r++) { }
        nrecord = r;

        int keyword_id = -1;
        int print_header = FALSE;

        // Zero the atom counter, 
        i = 0;
        for (r = 0; r < nrecord; r++) {
            // If this record is neither an ATOM nor a HETATM then print it, 
            // else print the new coordinates of this atom.
            keyword_id = parse_PDBQT_line(PDBQT_record[r]);
            if (keyword_id == PDBQ_ROOT) {
                // Print the header just before we print out the ROOT record
                print_header = TRUE;
            }
            if ((keyword_id == PDBQ_ATOM) || (keyword_id == PDBQ_HETATM)) {
                assert(i >= 0 && i < natom);
                // If the state_type is unbound, then ignore the per-atom intermolecular
                // emap and elec values; set these to 0.
                if (state_type == 1) {
                    // DOCKED
                    this_emap = (peratomE[i].vdW_Hb+peratomE[i].desolv >= 0.) ? 
		        min(peratomE[i].vdW_Hb+peratomE[i].desolv, MaxValue) 
		      : max(peratomE[i].vdW_Hb+peratomE[i].desolv, MinValue);
                    this_elec = (peratomE[i].elec >= 0.) ? 
		        min(peratomE[i].elec, MaxValue) : max(peratomE[i].elec, MinValue);
                } else {
                    // UNBOUND
                    this_emap = 0.;
                    this_elec = 0.;
                }
                if (keepresnum > 0) {
                    // Retain the original Residue Numbering (held in atomstuff)
                    print_PDBQT_atom_resstr(logFile, state_type_prefix_string, 
                                   i , atomstuff[i], crd,
                                   //i incremented in print_PDBQT_atom_resstr
                                   this_emap, this_elec, 
                                   charge[i], parameterArray[type[i]].autogrid_type, "\n" );
                } else {
                    // Change the residue number to the run number 
                    print_PDBQT_atom_resnum(logFile, state_type_prefix_string, 
                                   i , atomstuff[i], irun+1, crd,
                                   //i incremented in print_PDBQT_atom_resnum
                                   this_emap, this_elec, 
                                   charge[i], parameterArray[type[i]].autogrid_type, "\n" );

                }
                // Increment the atom counter
                i++;
            } else {
                if (print_header) {
                    (void) fprintf(logFile, "%s: USER                              x       y       z     vdW  Elec       q    Type\n", state_type_string);
                    (void) fprintf(logFile, "%s: USER                           _______ _______ _______ _____ _____    ______ ____\n", state_type_string);
                    // Make sure we only print the header once
                    print_header = FALSE;
                }
                (void) fprintf(logFile, "%s%s", state_type_prefix_string, PDBQT_record[r]);
            }
        } // r

        (void) fprintf(logFile, "%s: TER\n", state_type_string);
        (void) fprintf(logFile, "%s: ENDMDL\n", state_type_string);
        //(void) fprintf(logFile, UnderLine);
    } // outlev >= LOGFORADT
   (void) fflush(logFile);
} // writePDBQT()

void print_PDBQT( FILE *const logFile, 
                  const char *const prefix,
                  const int true_ligand_atoms,  // not necessarily 
                  const char atomstuff[MAX_ATOMS][MAX_CHARS], 
                  const Real crd[MAX_ATOMS][SPACE], 
                  const Real charge[MAX_ATOMS], 
                  const ParameterEntry parameterArray[MAX_ATOM_TYPES], 
                  const int type[MAX_ATOMS], 
                  const char *const suffix)
{ // Print out the coordinates
    for (int i=0; i<true_ligand_atoms; i++) {
        print_PDBQT_atom_resstr(logFile, prefix,
        i , atomstuff[i], crd,
        1., 0., 
       charge[i], parameterArray[type[i]].autogrid_type, suffix );
    }
    pr( logFile, "\n\n" );
} // end Print out the coordinates

void print_PDBQT_atom_resstr( FILE *const logFile, 
                  const char *const prefix, 
                  const int atom_num, // 0-origin 
                  const char *const atomstuff, 
                  const Real crd[MAX_ATOMS][SPACE], 
                  const Real vdW, 
                  const Real Elec, 
                  const Real charge, 
                  const char *const element, // 2-char AD type really eg HD, OA, Mg, Cl, Br ...
                  const char *const suffix //newline or empty
                  )
{
	char AtmNamResNamNumInsCode[20]; /* PDB record 0-origin indices 11-29 (from blank after serial_number to just before xcrd */
    sprintf(AtmNamResNamNumInsCode, "%-19.19s", &atomstuff[11]);
    // see constants.h for FORMAT_PDBQT_ATOM_RESSTR        
    pr(logFile, FORMAT_PDBQT_ATOM_RESSTR, prefix, atom_num+1, AtmNamResNamNumInsCode, 
       crd[atom_num][X], crd[atom_num][Y], crd[atom_num][Z], 
       vdW, Elec,  charge, element, suffix);
}

void print_PDBQT_atom_resnum( FILE *const logFile, 
                  const char *const prefix, 
                  const int atom_num, // 0-origin 
                  const char *const atomstuff, 
                  const int resnum, 
                  const Real crd[MAX_ATOMS][SPACE], 
                  const Real vdW, 
                  const Real Elec, 
                  const Real charge, 
                  const char *const element, 
                  const char *const suffix //newline or empty
                  )
{
	char AtmNamResNamNumInsCode[20]; /* PDB record 0-origin indices 11-29 (from blank after serial_number to just before xcrd */
    sprintf(AtmNamResNamNumInsCode, "%11.11s%4d%4.4s", &atomstuff[11], 
       resnum, &atomstuff[26]);
    pr(logFile, FORMAT_PDBQT_ATOM_RESSTR, prefix, atom_num+1, AtmNamResNamNumInsCode, 
       crd[atom_num][X], crd[atom_num][Y], crd[atom_num][Z], 
       vdW, Elec,  charge, element, suffix);
}




/* EOF */
