/*

 $Id: readPDBQT.cc,v 1.46 2014/07/10 19:42:02 mp Exp $

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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/param.h>
#include <time.h>
#include <ctype.h>		/* tolower */
#include <search.h>
#include "readPDBQT.h"
#include "PDBQT_tokens.h"
#include "structs.h"
#include "atom_parameter_manager.h"
#include "stack.h"
#include "timesys.h" // for struct tms

/*----------------------------------------------------------------------------*/

extern int      parse_tors_mode;
extern char    *programname;

/*----------------------------------------------------------------------------*/
Molecule readPDBQT(char input_line[LINE_LEN],
                    const int num_atom_maps,

                    int *const P_natom,
                    Real crdpdb[MAX_ATOMS][NTRN],
                    Real charge[MAX_ATOMS],
                    Boole *const P_B_haveCharges,
                    int map_index[MAX_ATOMS], //was:int type[MAX_ATOMS]
                    int bond_index[MAX_ATOMS],
                    char pdbaname[MAX_ATOMS][5],

                    char *const FN_ligand,
                    char *const FN_flexres,
                    const Boole B_have_flexible_residues,

                    char atomstuff[MAX_ATOMS][MAX_CHARS],
                    int *const P_n_heavy_atoms_in_ligand,
		    int *P_true_ligand_atoms, // set

                    Boole *const P_B_constrain,
                    int *const P_atomC1,
                    int *const P_atomC2,
                    Real *const P_sqlower,
                    Real *const P_squpper,

                    int *const P_ntor, // number of torsions in lig + flexres
                    int *const P_ntor_ligand,   // the number of torsions in the ligand (excluding the flexible residues in receptor)
                    int tlist[MAX_TORS][MAX_ATOMS],
                    Real vt[MAX_TORS][NTRN],

                    int *const P_Nnb,
		    int Nnb_array[3],
                    NonbondParam *const nonbondlist,

                    const Clock& jobStart,
                    const struct tms& tms_jobStart,
                    const char *const hostnm,
                    int *const P_ntorsdof,
                    int ignore_inter[MAX_ATOMS],
                    const int B_include_1_4_interactions,

                    /* not const */ char PDBQT_record[MAX_RECORDS][LINE_LEN],

                    /* not const */ int end_of_branch[MAX_TORS],
		    const int debug,
                    const int outlev,
		    FILE *logFile
                    )

{
	FILE           *FP_ligand;
	FILE           *FP_flexres;
	char     dummy[LINE_LEN];
	char     error_message[LINE_LEN];
	char     message[LINE_LEN];

	Real   aq = 0.;
	Real   lq = 0.;
	Real   total_charge_ligand = 0.;
	Real   uq = 0.;

    int             branch_last_piece[MAX_TORS+2];//0 unused, 1 means ROOT
	static int      atomnumber[MAX_RECORDS]; /* -1 for non-ATOM/HETATM */ // static to save stack space
	int pdbatomnumber[2][MAX_RECORDS]; /* two namespaces; each records's PDB atom number */
	int             iq = 0;
	int             natom = 0;
	static int      nbmatrix[MAX_ATOMS][MAX_ATOMS]; // static to save stack space
	int             nrecord = 0;
    int             nligand_record = 0;
	int             ntor = 0;
	static int      ntype[MAX_ATOMS]; // static to save stack space
	static int      rigid_piece[MAX_ATOMS]; // static to save stack space
    int             keyword_id = -1;
	int             nres = 0;
	int             nrigid_piece = 0;
	int             piece; 
	int		errorcode = 0;

	Boole           B_found_begin_res = FALSE; //tracks whether any flexres found
	Boole           B_found_first_flexres = FALSE; //tracks whether first flexres ROOT found
	Boole           B_has_conect_records = FALSE;
	Boole           B_is_in_branch = FALSE;
	Boole           B_is_in_flexres = FALSE;


	register int    i = 0;
	register int    j = 0;

	const static Real QTOL = 0.050;//rh increased from 0.005 4/2009

    // Definitions to help determine the end_of_branch array for "Branch Crossover Mode".
    // The "end_of_branch" array is like a dictionary, where the index corresponds to the key
    // and is the number of the torsion to be crossed over, while the value is the number
    // of the first torsion after the last torsion in the sub-tree (or "branch") being exchanged.
    stack s;  // Stack used to determine the values for the end_of_branch[] array
    int parent;  

	Molecule        mol;

	ParameterEntry  this_parameter_entry;

	for (i = 0; i < MAX_RECORDS; i++) {
		atomnumber[i] = -1;
		pdbatomnumber[0][i] = 9999999; // no PDB atom number 
		pdbatomnumber[1][i] = 9999999; // no PDB atom number 
	}

	for (i = 0; i < MAX_ATOMS; i++) {
		ntype[i] = 0;
		rigid_piece[i] = 0;
	}

	for (i = 0; i < MAX_TORS+2; i++) {
        branch_last_piece[i] = 0;//0 unused, 1 means ROOT
    }

    // Create the stack
    s = stack_create( MAX_TORS+2 );
    if (debug > 0) {
        pr(logFile, "DEBUG: 170:  stack_create(%d)\n", MAX_TORS+2 );
    }

    // Initialize the stack
    piece = -1; // "sentinel" value
    stack_push(s, piece);
    if (debug > 0) {
        pr(logFile, "DEBUG: 176:  stack_push(s, %d)\n", piece );
    }

    //  Attempt to open the ligand PDBQT file...
	sscanf(input_line, "%*s %s", FN_ligand);
	if (openFile(FN_ligand, "r", &FP_ligand, jobStart, tms_jobStart, TRUE, logFile)) {
		pr(logFile, "Ligand PDBQT file = \"%s\"\n\n", FN_ligand);
	}

    if (B_have_flexible_residues) {
        //  Attempt to open the flexible residue PDBQT file...
        if (openFile(FN_flexres, "r", &FP_flexres, jobStart, tms_jobStart, TRUE, logFile)) {
            if(outlev>=LOGMIN)
            pr(logFile, "Flexible Residues PDBQT file = \"%s\"\n\n", FN_flexres);
        }
    }

    //  Count the number of records in the Ligand PDBQT file first...
	nligand_record = 0;
	while (fgets(dummy, LINE_LEN, FP_ligand) != NULL) {
		++nligand_record;
	}
	(void) fclose(FP_ligand);

    nrecord = nligand_record;

    if (B_have_flexible_residues) {
        //  Count the number of records in the Flexible Residue PDBQT file next, 
        //  but don't reset the counter...
        while (fgets(dummy, LINE_LEN, FP_flexres) != NULL) {
            ++nrecord;
        }
        (void) fclose(FP_flexres);
    }

    // Set the nrecord-th entry of PDBQT_record to NULL, aka '\0'
    (void) strcpy(PDBQT_record[nrecord], "\0");

    // Read in the PDBQT file(s) if there are not too many records
	if (nrecord > MAX_RECORDS) {
		prStr(error_message, "ERROR: %d records read in, but only dimensioned for %d.\nChange \"MAX_RECORDS\" in \"constants.h\".", nrecord, MAX_RECORDS);
		stop(error_message);
	} else {
        // Read in the input Ligand PDBQT file...
		if (openFile(FN_ligand, "r", &FP_ligand, jobStart, tms_jobStart, TRUE, logFile)) {
			if(outlev>=LOGBASIC){
			  pr(logFile,   "INPUT LIGAND PDBQT FILE:");
			  pr(logFile, "\n________________________\n\n");
			 }
			for (i = 0; i < nligand_record; i++) {
				if (fgets(PDBQT_record[i], LINE_LEN, FP_ligand) != NULL) {
					if(outlev>=LOGBASIC)
					pr(logFile, "INPUT-LIGAND-PDBQT: %s", PDBQT_record[i]);
				}
			} // i
			if(outlev>=LOGBASIC) pr(logFile, UnderLine);
		} // if
		(void) fclose(FP_ligand);

        if (B_have_flexible_residues) {
            // Read in the input Flexible Residues PDBQT file...
            if (openFile(FN_flexres, "r", &FP_flexres, jobStart, tms_jobStart, TRUE, logFile)) {
		if(outlev>=LOGBASIC) {
                   pr(logFile,   "INPUT FLEXIBLE RESIDUES PDBQT FILE:");
                   pr(logFile, "\n___________________________________\n\n");
		   }
                for (i = nligand_record; i < nrecord; i++) {
                    if (fgets(PDBQT_record[i], LINE_LEN, FP_flexres) != NULL) {
                        if(outlev>=LOGBASIC) pr(logFile, "INPUT-FLEXRES-PDBQT: %s", PDBQT_record[i]);
                    }
                } // i
                if(outlev>=LOGBASIC) pr(logFile, UnderLine);
            } // if
            (void) fclose(FP_flexres);
        }

	} // if

	// Count the ATOMs and HETATMs; store the (x,y,z) coordinates...
	// Also, check for any BEGIN_RES records, for receptor flexibility...
	if(outlev>=LOGLIGREAD) {
	pr(logFile, "\nDetermining Atom Types and Parameters for the Moving Atoms\n");
	pr(logFile,   "__________________________________________________________\n\n");
	}
	natom = 0;
    //debug = 1;

    // Loop over all the lines in either the ligand or the "reconstructed" combined-ligand-flexible-residues file
	for (i = 0; i < nrecord; i++) {
		strncpy(input_line, PDBQT_record[i], (size_t) LINE_LEN);
		// Parse this line in the ligand file
		keyword_id = parse_PDBQT_line(input_line);

        switch ( keyword_id ) {
            case PDBQ_ATOM:
            case PDBQ_HETATM:
                if ( ! B_is_in_branch ) {
                    // Flag this as an error
                    // Incorrectly nested BRANCH/ENDBRANCH records
		    prStr(error_message,
                    "%s: ERROR:  All ATOM and HETATM records must be given before any nested BRANCHes; see line %d in PDBQT file \"%s\".\n\n", programname, i+1, FN_ligand);
                    stop(error_message);
                }
                
                // Check that the line is at least 78 characters long
                if (strlen(input_line) < 78) {
		    prStr(error_message,
                    "%s: FATAL ERROR: line %d is too short!\n", programname, i+1);
                    stop(error_message);
                }

                ParameterEntry * found_parm;
                int pdbanum;

                // Set up rigid_piece array by reading in the records of the PDBQT file;
                // each "rigid_piece" is a self-contained rigid entity.
                rigid_piece[natom] = piece;
                if (debug > 0) {
                    pr(logFile, "DEBUG: 289:  rigid_piece[%d] = %d (nrigid_piece)\n", natom, rigid_piece[natom] );
                    pr(logFile, "DEBUG: 290:  for natom %d piece=%d\n", natom, piece );
                }

                // Read the coordinates and store them in crdpdb[],
                // read the partial atomic charge and store it in charge[],
                // and read the parameters of this atom and store them in this_parameter_entry,
                // setting the "autogrid_type" in this_parameter_entry.
                readPDBQTLine(input_line, &pdbanum, crdpdb[natom], &charge[natom], &this_parameter_entry,outlev,logFile);

                /*
                // Verify the pdb number for this atom
                if ( pdbanum != (natom + 1) ) {
		    prStr(error_message,
                    "%s: ERROR:  ATOM and HETATM records must be numbered sequentially from 1.  See line %d in PDBQT file \"%s\".\n\n", programname, i+1, FN_ligand);
                    stop(error_message);
                }
                */

                // Set the serial 0-origin index atomnumber[i] for this record
                atomnumber[i] = natom;
		// Save the record's PDB  atom number, in one of the two
		// namespaces: 0 for true ligand, 1 for flex res
		pdbatomnumber[ B_found_begin_res ? 1 : 0 ][i] =  pdbanum;
#ifdef DEBUG
	fprintf(logFile, "DEBUG pdbatomnumber[%d][rec %2d]= %d\n",
	B_found_begin_res ? 1 : 0, i, pdbanum);
#endif

                for (j = 0; j < NTRN; j++) {
                    mol.crdpdb[natom][j] = crdpdb[natom][j];
                    mol.crd[natom][j] = crdpdb[natom][j];
                }

                if (!B_found_begin_res) {
                    // Only accumulate charges on the ligand...
                    total_charge_ligand += charge[natom];
                }
                *P_B_haveCharges = TRUE;

                strncpy(atomstuff[natom], input_line, (size_t) 30);
                atomstuff[natom][30] = '\0';
                strcpy(mol.atomstr[natom], atomstuff[natom]);

                //sscanf(&input_line[12], "%s", pdbaname[natom]);
		// safer:
		strncpy(pdbaname[natom], &input_line[12],  
		  sizeof(pdbaname[natom])-1);
		pdbaname[natom][sizeof(pdbaname[natom])-1]='\0';

                // "map_index" is used as an index into the AutoGrid "map" array to look up 
                // the correct energies in the current grid cell, thus:	map[][][][map_index[natom]]
                map_index[natom] = -1;

                // Find the AutoDock 4 atom type for the current atom.
                // "apm_find" is the new AutoDock 4 atom typing mechanism
                found_parm = apm_find(this_parameter_entry.autogrid_type);
                if (found_parm != NULL) {
                    // We found the atom type and its parameters
                    map_index[natom] = found_parm->map_index;
                    bond_index[natom] = found_parm->bond_index;
                    if (outlev >= LOGETABLES) {
                        (void) fprintf(logFile, "Found parameters for atom %d, atom type \"%s\", grid map index = %d\n",
                                   natom + 1, found_parm->autogrid_type, found_parm->map_index);
                    }
                } else {
                    // We could not find this parameter -- return an error
                    prStr(message, "\n%s: *** WARNING!  Unknown atom type \"%s\" found.  You should add parameters for it to the parameter library first! ***\n\n", programname, this_parameter_entry.autogrid_type);
                    pr_2x(stderr, logFile, message);
		    stop("unknown atom type");
                }

                if (map_index[natom] == -1) {
                    prStr(message, "%s: WARNING: the atom type (%s) of atom number %d could not be found;\n\tcheck that this atom type is listed after the \"ligand_types\" keyword in the DPF,\n\tand make sure to add a \"map\" keyword to the DPF for this atom type.\n\tNote that AutoDock will use the default atom type = 1, instead.\n\n", programname, this_parameter_entry.autogrid_type, natom+1);
                    pr_2x(stderr, logFile, message);
                    map_index[natom] = 0; // we are 0-based internally, 1-based in printed output
                }

                // Increment the number of atoms having this atomtype
                ++ntype[map_index[natom]];

                // Increment the number of atoms found in PDBQT file
                ++natom;

                // Increment the number of non-hydrogen atoms in the ligand
                if (!B_found_begin_res) {
                    //(void) fprintf(logFile, "DEBUG!!! Found parameters for atom %d, atom type \"%s\", grid map index = %d\n",
                               //natom, found_parm->autogrid_type, found_parm->map_index);
                    //(void) fprintf(logFile, "DEBUG!!!  (strcmp(t,\"H\")==%d) || (strcmp(t,\"HD\")==%d) || (strcmp(t,\"HS\")==%d))\n\n",
                                //strcmp(found_parm->autogrid_type,"H"), strcmp(found_parm->autogrid_type,"HD"), strcmp(found_parm->autogrid_type,"HS") );
                    if ( ! is_hydrogen_type( found_parm->autogrid_type ) ) {
                        //(void) fprintf(logFile, "DEBUG!!! is not a hydrogen atom.  Incrementing number of heavy atoms in ligand.\n");
                        ++(*P_n_heavy_atoms_in_ligand);
                    }
                }
                break;

            case PDBQ_ROOT:
                B_is_in_branch = TRUE;
                if (!B_is_in_flexres ){
                    ++nrigid_piece; //allocate a new rigid body for ligand root atoms 
                    piece = nrigid_piece; //CAUTION: ligand root must have value 1
                }
                else if  (!B_found_first_flexres) {
                    // allocate a new rigid body for FIRST flexres root atoms 
                    // Subsequent flexres root atoms will use this same piece
                    // number (so we do not run out of rigid body pieces (NTORS))
                    B_found_first_flexres = TRUE;
                    ++nrigid_piece; 
                    piece = nrigid_piece; //CAUTION: 
                }
                stack_push(s, piece);
                if (debug > 0) {
                    pr(logFile, "DEBUG: 382: ROOT piece=%d\n", piece);
                }
                break;

            case PDBQ_BRANCH:
                B_is_in_branch = TRUE;
                if (nrigid_piece>MAX_TORS){
		            prStr(error_message, "PDBQT ERROR: too many torsions, maximum number of torsions is %d", MAX_TORS);
		            stop(error_message);
                }
                stack_push(s, piece);//at this pt push the parent piece number
                if (debug > 0) {
                    pr(logFile, "DEBUG: 386:  stack_push(s, %d)\n", piece );
                }
                ++nrigid_piece; //allocate a new rigid body for this branch
                piece = nrigid_piece;
                if (debug > 0) {
                    pr(logFile, "DEBUG: 391:  piece = %d\n", piece );
                    pr(logFile, "DEBUG: 394:  nrigid_piece = %d\n", nrigid_piece );
                }
                break;

            case PDBQ_ENDROOT:
                B_is_in_branch = FALSE;
                (void) stack_pop(s);
                break;

            case PDBQ_ENDBRANCH:
                B_is_in_branch = FALSE;
                //end_of_branch[ piece ] = max( end_of_branch[ piece ], piece);
                branch_last_piece[ piece ] = max( branch_last_piece[ piece ], piece);
                if (debug > 0) {
                    pr(logFile, "DEBUG: 407:   END BRANCH: branch_last_piece[ piece %d ] = %d\n", piece, branch_last_piece[ piece ] );
                }
                parent = stack_pop(s);
                if (parent < 0) {
		            prStr(error_message, "PDBQT ERROR: Encountered end of branch without corresponding branch");
		            stop(error_message);
                }
                if (parent>1){ 
                    //if parent is 1, it is the 'ligand' root. In that case,we don't have an eob for it.
                    //end_of_branch[ parent ] = max( end_of_branch[ parent ], end_of_branch[ piece ] );
                    branch_last_piece[parent ] = max( branch_last_piece[parent ], branch_last_piece[piece ] );
                } 
                if (debug > 0) {
                    pr(logFile, "DEBUG: 425:   parent = %d, branch_last_piece[ parent ] = %d\n", parent, branch_last_piece[ parent ] );
                }
                if (debug > 0) {
                    pr(logFile, "DEBUG: 428:   branch_last_piece[ %d ] = %d\n", parent, branch_last_piece[ parent ] );
                }
                piece = parent;
                break;

            case PDBQ_NULL:
            case PDBQ_REMARK:
            case PDBQ_TORS:
            case PDBQ_ENDTORS:
            case PDBQ_TORSDOF:
            case PDBQ_CONSTRAINT:
                break;

            case PDBQ_CONECT:
                // At least some of the atoms in the "ligand" may have their connectivity specified
                // so we could set up their bonded entries. For future versions...
                B_has_conect_records = TRUE;
                break;

            case PDBQ_BEGIN_RES:
                // Then a flexible receptor sidechain was found in the PDBQ file.
                // Flag that we've found a BEGIN_RES record.
                B_found_begin_res = TRUE;
                if (debug > 0) {
                    pr(logFile, "DEBUG: 457: BEGIN_RES\n");
                }
                // Increment the number of residues
                nres++;
	            B_is_in_flexres = TRUE;
                break;

            case PDBQ_END_RES:
	            B_is_in_flexres = FALSE;
                break;

            case PDBQ_UNRECOGNIZED:
            default:
                pr(logFile, "%s: WARNING: Unrecognized PDBQT record type in line:\n", programname );
                pr(logFile, "%s: WARNING: %s\n", programname, input_line );
                break;

        } // end switch( keyword_id )

		if (!B_found_begin_res) {
			// No BEGIN_RES tag has been found yet.
            // Keep updating "true_ligand_atoms" until we find a "BEGIN_RES".
            // "true_ligand_atoms" is the number of atoms in the moving ligand, 
            // and excludes all atoms in the flexible sidechain residues of the receptor.
			*P_true_ligand_atoms = natom;
        }

	} // i, next record in PDBQT file 

    if(outlev>=LOGLIGREAD) {
    pr(logFile, "\nNumber of atoms in movable ligand = %d\n\n", *P_true_ligand_atoms);

    pr(logFile, "Number of non-hydrogen atoms in movable ligand = %d\n\n", *P_n_heavy_atoms_in_ligand);

	pr(logFile, "\nNumber of atoms found in flexible receptor sidechains (\"residues\") =\t%d atoms\n\n", natom - *P_true_ligand_atoms);

	pr(logFile, "Total number of atoms found =\t%d atoms\n\n", natom);

	pr(logFile, "Number of flexible residues found in the receptor =\t%d residues\n\n", nres);
   }


    //check for mismatched BRANCH/ENDBRANCH records
    // specifically:
    // stack_pop should return the sentinel ie -1
    if (B_is_in_branch || -1 != stack_pop(s)){
	    prStr(error_message, "ERROR: BRANCH statement without a corresponding ENDBRANCH\n\n");
        stop(error_message);
    }
    // create end_of_branch array
    // 0th element refers to first torsion
    // Each entry holds the index of the first torsion NOT in this branch
    // Since our piece numbers start at 1 for the ROOT, the atoms moved by
    // the first torsion are piece number 2. Thus we have to look 2 higher
    // in the branch_last_piece array. (Note the "-1" is here because 
    // end_of_branch is zero-based).
    for (j = 0; j < nrigid_piece-1; j++)  {
        end_of_branch[j] = branch_last_piece[j+2]-1; 
        if (debug > 0) {
            pr(logFile, "\n end_of_branch[%2d] := %2d (branch_last_piece[%2d]-1)",
                                          j,end_of_branch[j],j+2);
        };
        if (end_of_branch[j] <0) {
		    prStr(error_message, "ERROR: end_of_branch[%d] < 0 (%d)\n", j, end_of_branch[j]);
		    stop(error_message);
        }
    }

	if (natom > MAX_ATOMS) {
		prStr(error_message, "ERROR: Too many atoms found (i.e. %d); maximum allowed is %d.\nChange the \"#define MAX_ATOMS\" line in \"constants.h\"\n.", natom, MAX_ATOMS);
		stop(error_message);
	} else {
		*P_natom = natom;
		mol.natom = natom;
	}

	if(outlev>=LOGLIGREAD) {
	pr(logFile, "\nSummary of number of atoms of a given atom type:\n");
	pr(logFile, "------------------------------------------------\n\n");
	for (i = 0; i < num_atom_maps; i++) {
		pr(logFile, "Number of atoms with atom type %d = %2d\n", i + 1, ntype[i]);
	}

	pr(logFile, "\n\nSummary of total charge on ligand, residues and overall:\n");
	pr(logFile, "-------------------------------------------------------\n");
    }

    // Check total charge on ligand
	if(outlev>=LOGBASIC)
	pr(logFile, "Total charge on ligand                               =\t%+.3f e\n", total_charge_ligand);
	iq = (int) ((aq = fabs(total_charge_ligand)) + 0.5);
	lq = iq - QTOL;
	uq = iq + QTOL;
	if (!((aq >= lq) && (aq <= uq))) {
		prStr(message, "\n%s: *** Caution!  Non-integral total charge (%.3f e) on ligand may indicate a problem... ***\n\n", programname, total_charge_ligand);
		pr_2x(stderr, logFile, message);
	}

	/*
	 * Work out where the torsions are; and what they move...
	 *
	 * Also, detect which atoms we should ignore in the
	 * intermolecular energy calculation (ignore_inter[MAX_ATOMS]
	 * array)
	 */
	mkTorTree(atomnumber, pdbatomnumber, PDBQT_record, nrecord, tlist, &ntor, P_ntor_ligand, FN_ligand, pdbaname,
              P_B_constrain, P_atomC1, P_atomC2, P_sqlower, P_squpper, P_ntorsdof, ignore_inter,
	      *P_true_ligand_atoms, outlev, logFile);

	*P_ntor = ntor;

	mol.S.ntor = ntor;

	if (ntor > 0) {
        //  Create a list of internal non-bonds to be used in internal energy calculation...
	//int		nbonds[MAX_ATOMS];
	//int             bonded[MAX_ATOMS][MAX_NBONDS];
	int		*nbonds; // [natom];
	int             (*bonded) /* [natom] */ [MAX_NBONDS];
		if (debug > 0) {
			pr(logFile, "Finding bonds.\n\n");
		}
		// Initialise the bonded and nbonds arrays
	nbonds=(int *) calloc(natom, sizeof *nbonds);
	bonded=(int (*)[MAX_NBONDS])  malloc(natom * sizeof *bonded);
	if(nbonds==NULL || bonded==NULL) stop("memory alloc failure in readPDBQT");
        for (i = 0; i < natom; i++)  for (j = 0; j < MAX_NBONDS; j++) {
		bonded[i][j] = -1;
	}

        if (debug > 0) {
			printbonds(natom, nbonds, bonded, "\nDEBUG:  1. BEFORE getbonds, bonded[][] array is:\n\n", 1, outlev, logFile);
		}
        // find all the bonds in the ligand
		errorcode = getbonds(crdpdb, 0, *P_true_ligand_atoms, bond_index, nbonds, bonded, debug, outlev, logFile);
	   if(errorcode!=0)  {
			pr(logFile, " ERROR in ligand getbonds, code=%d\n", errorcode);
			stop(" ERROR in ligand getbonds");
			}

        if (B_have_flexible_residues) {
            // find all the bonds in the receptor
            errorcode = getbonds(crdpdb, *P_true_ligand_atoms, natom, bond_index, nbonds, bonded, debug, outlev, logFile);
	   if(errorcode!=0)  {
			pr(logFile, " ERROR in receptor getbonds, code=%d\n", errorcode);
			stop(" ERROR in receptor getbonds");
			}
        }

		if (debug > 0) {
			printbonds(natom, nbonds, bonded, "\nDEBUG:  2. AFTER getbonds, bonded[][] array is:\n\n", 0, outlev, logFile);
			pr(logFile, "Detecting all non-bonds.\n\n");
		}
		errorcode = nonbonds(crdpdb, nbmatrix, natom, bond_index, B_include_1_4_interactions, nbonds, bonded, debug, outlev, logFile);
		if(errorcode!=0)  {
			pr(logFile, " ERROR in nonbonds, code=%d\n", errorcode);
			stop(" ERROR in nonbonds");
			}

		if (debug > 0) {
			printbonds(natom, nbonds, bonded, "\nDEBUG:  4. AFTER nonbonds, bonded[][] array is:\n\n", 0, outlev, logFile);
			pr(logFile, "Weeding out non-bonds in rigid parts of the torsion tree.\n\n");
		}
		weedbonds(natom, pdbaname, rigid_piece, ntor, tlist, nbmatrix, P_Nnb, nonbondlist, Nnb_array, *P_true_ligand_atoms, map_index, debug, outlev, logFile);

		print_nonbonds(natom, pdbaname, rigid_piece, ntor, tlist, nbmatrix, *P_Nnb, nonbondlist, map_index, outlev, logFile);

        // Update the unit vectors for the torsion rotations
        update_torsion_vectors( crdpdb, ntor, tlist, vt, &mol, debug, 
	 outlev, logFile);

		flushLog;

	free(nbonds);
	free(bonded);
	} else {
		if(outlev>=LOGLIGREAD)
		fprintf(logFile, ">>> No torsions detected, so skipping \"nonbonds\", \"weedbonds\" and \"torNorVec\" <<<\n\n");
	}

    //  End program if just parsing torsions...
	if (parse_tors_mode) {
		prStr(message, "\n\n *** PARSE TORSIONS MODE - Stopping here ***\n\n");
		success(hostnm, jobStart, tms_jobStart, logFile);
		exit(EXIT_SUCCESS);
	}
	return mol;
}


/*----------------------------------------------------------------------------*/
/* readPDBQTLine.cc */

void
readPDBQTLine( char line[LINE_LEN],
               int  *const ptr_pdbatomnumber,
               Real crd[SPACE],
               Real *const ptr_q,
               ParameterEntry *const this_parameter_entry,
	       const int outlev,
	       FILE *logFile)
/*----------------------------------------------------------------------------*/
{
    char char8[9];
    char char6[7];
    char char5[6];
    char char2[3];
	char message[LINE_LEN];

    // Initialise char5
    (void) strcpy( char5, "    0" );
    char5[5] = '\0';

    // Initialise char8
    (void) strcpy( char8, "   0.000" );
    char8[8] = '\0';

    // Initialise char6
    (void) strcpy( char6, "  0.00" );
    char6[6] = '\0';

    // Initialise char2
    (void) strcpy( char2, "C " );
    char2[2] = '\0';

#define check_sscanf( str, fmt, val, fieldname )  if (1 != sscanf( str, fmt, val ))  {\
    prStr(message, "\n%s: WARNING! Could not read " fieldname " in PDBQT line \"%s\".\n", programname, line );\
    pr_2x(stderr, logFile, message); }

    // Read in the PDB record's number for this atom
    (void) strncpy( char5, &line[6], (size_t)5 );
    char5[5] = '\0';
    check_sscanf( char5, "%d", ptr_pdbatomnumber, "PDB atom number" );

	// Read in the X, Y, Z coordinates
    (void) strncpy( char8, &line[30], (size_t)8 );
    char8[8] = '\0';
    check_sscanf( char8, FDFMT, &crd[X], "x-coordinate" );

    (void) strncpy( char8, &line[38], (size_t)8 );
    char8[8] = '\0';
    check_sscanf( char8, FDFMT, &crd[Y], "y-coordinate" );

    (void) strncpy( char8, &line[46], (size_t)8 );
    char8[8] = '\0';
    check_sscanf( char8, FDFMT, &crd[Z], "z-coordinate" );

#ifdef DEBUG
	(void) fprintf(stderr, "readPDBQTLine: %s", line);
#endif				/* DEBUG */

    // partial charge, q
    (void) strncpy( char6, &line[70], (size_t)6 );
    char6[6] = '\0';
    check_sscanf( char6, FDFMT, ptr_q, "partial charge" );

    // atom type name
    (void) strncpy( char2, &line[77], (size_t)2 );
    char2[2] = '\0';
    check_sscanf( char2, "%s", this_parameter_entry->autogrid_type, "atom type" );

#undef check_sscanf

#ifdef DEBUG
	fprintf(stderr, "readPDBQTLine:  %d, %.3f, %.3f, %.3f, %.3f, %s\n", *ptr_pdbatomnumber, crd[X], crd[Y], crd[Z], *ptr_q,
    this_parameter_entry->autogrid_type);
#endif				/* DEBUG */
}

/* EOF */
