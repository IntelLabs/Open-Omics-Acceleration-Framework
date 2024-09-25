/*

 $Id: mkTorTree.cc,v 1.28 2014/06/12 01:44:07 mp Exp $

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
#include <ctype.h>
#include "mkTorTree.h"
#include "PDBQT_tokens.h"


extern char  *programname;
    
/* local functions: */
static
int check_atomnumber( const int number, 
 const int atomnumber[], const int nrecord, const int natoms, FILE *logFile);
//MP in progress #define check_atomnumber_ok(a) check_atomnumber((a), atomnumber, nrecord, natoms, logFile)

#define check_atomnumber_ok( a )  (((a) >= 0) && ((a) < natoms))

static int pdbatomnum_find_record( int pdbn, int pdbatomnumber[], int nrecord);

void mkTorTree( const int   atomnumber[ MAX_RECORDS ],
		const int   pdbatomnumber[2][MAX_RECORDS],
                const char  Rec_line[ MAX_RECORDS ][ LINE_LEN ],
                const int   nrecord,

/* not const */ int   tlist[ MAX_TORS+1][ MAX_ATOMS ],
/* not const */ int   *const P_ntor,
/* not const */ int   *const P_ntor_ligand,

          const char  *const smFileName,

          const char  pdbaname[ MAX_ATOMS ][ 5 ],
/* not const */ Boole *const P_B_constrain,
/* not const */ int   *const P_atomC1,
/* not const */ int   *const P_atomC2,
/* not const */ Real  *const P_sqlower,
/* not const */ Real  *const P_squpper,
/* not const */ int   *const P_ntorsdof,
/* not const */ int   ignore_inter[MAX_ATOMS],
	int true_ligand_atoms,
	int outlev,
	FILE *logFile)

{

    int   itor = 0;
    int   keyword_id = -1;
    int   nbranches = 0;
    int   ntor=0;
    int   tlistsort[ MAX_TORS+1][ MAX_ATOMS ];
    Boole   found_new_res = false;
    int   nres = 0;
    int   natoms_in_res = 0;
    int   natoms = 0;
    Boole   found_first_res=false;

    register int   i = 0;
    register int   j = 0;
    register int   k = 0;


    Real lower = 0.;
    Real temp  = 0.;
    Real upper = 0.01;

#ifdef DEBUG
    int   oo = 0;
    char  C = 'A';
#endif /* DEBUG */

    /* ntor = *P_ntor; */

    /* TorsionTree torstree; */

    for (i = 0; i  < MAX_TORS;  i++ ) {
        for (j = 0;  j < MAX_ATOMS;  j++ ) {
            tlistsort[ i ][ j ] = 0;
        }
    }

    i = 0;
    j = 0;

    /* ________________________________________________________________
      |  Work out the torsion angle tree.                              |
      |________________________________________________________________|
      |  tlist contains information on torsion angles:                 |
      |                                                                |
      |  tlist[ i ][ ATM1 ]           = | atom numbers defining torsion|
      |  tlist[ i ][ ATM2 ]             |                              |
      |  tlist[ i ][ NUM_ATM_MOVED ]  = number of atoms to be rotated  |
      |  tlist[ i ][ 3 ] and on       = atom IDs to be rotated.        |
      |                                                                |
      |  tlist[ ntor ][  ] holds atom identifiers of true ligand root  |
      |________________________________________________________________|
      | NOTE:  code does not explicitly check for keyword 'ROOT'.      |
      |________________________________________________________________|
    */

#ifdef DEBUG
    pr( logFile, "\n\n" );
    pr( logFile, "                                                  |Atoms|Total\n" );
    pr( logFile, "                                                  | 1| 2|Moved\n" );
    pr( logFile, "   Atom                                           |__|__|__|\n" );
    pr( logFile, "i  #  rec5 C atomlast nbranches j  ntor tlist[  ] [ 0| 1| 2| 3  4  5  6  7  8  9  10 11 12 13 14 15 ]\n" );
    pr( logFile, "__ __ ____ _ ________ _________ __ ____ _____________________________________________________________\n" );
#endif /* DEBUG */

    for (i = 0;  i < nrecord;  i++) {

        if ( (keyword_id = parse_PDBQT_line( Rec_line[ i ] )) == -1) {
            pr( logFile, "%s: Unrecognized keyword found while parsing PDBQT file, line:\n|%s|\n", programname, Rec_line[ i ] );
            continue;
        }

#ifdef DEBUG
        pr( logFile, "PDBQT-Line %d: %s", i+1, Rec_line[i] );
#endif /* DEBUG */

        switch( keyword_id ) {

    /*____________________________________________________________*/
            case PDBQ_REMARK:

		if(outlev>=LOGFORADT)
                pr( logFile, "%s", Rec_line[ i ] );
                break;

    /*____________________________________________________________*/
            case PDBQ_NULL:

                break;

    /*____________________________________________________________*/
            case PDBQ_ATOM: 
            case PDBQ_HETATM:

             /* This is an ATOM or HETATM. */

             if (found_new_res) {
                    /* We are in a residue. */
                    if (natoms_in_res < 2) {
                        /* flag the first two atoms in each new 
                         * residue to prevent them being
                         * included in the intermolecular
                         * energy calculation.  */
                         ignore_inter[natoms] = 1;
                    }
                    /* Keep counting the number of atoms in the residue. */
                    natoms_in_res++;
                } else {
                    /* We are not in a residue.
                     *
                     * "found_new_res" can only be reset to FALSE 
                     * if we encounter an "END_RES" record. 
                     * By default, found_new_res is set to FALSE. */
                     
                    /* reset the atom counter */
                    natoms_in_res = 0;
                }
                /* Increment atom counter for all atoms in PDBQT file */
                natoms++;

#ifdef DEBUG
                C = 'A';
                PrintDebugTors;
                PrintDebugTors2;
                pr( logFile, "]\n" );
#endif /* DEBUG */

                break;

    /*____________________________________________________________*/
            case PDBQ_BRANCH:
	    case PDBQ_TORS:  /* obsolete synonym for "branch" */
                if (ntor >= MAX_TORS) {
    		    char  error_message[ LINE_LEN ];
                    prStr( error_message, "ERROR: Too many torsions have been found (i.e. %d); maximum allowed is %d.\n Either: change the \"#define MAX_TORS\" line in constants.h\n Or:     edit \"%s\" to reduce the number of torsions defined.", (ntor+1), MAX_TORS, smFileName );
                    stop( error_message );
                }

		// read ATM1 and ATM2 for this BRANCH, convert from PDB to index
		for(int a=0;a<2;a++) { 
		  int pdbnum, recnum;
                    sscanf(Rec_line[ i ],(a==0)?"%*s %d %*d":"%*s %*d %d",
			&pdbnum );
		    recnum = pdbatomnum_find_record(pdbnum, 
		      (int *)pdbatomnumber[found_first_res?1:0], nrecord); 
#ifdef DEBUG
fprintf(logFile, "DEBUG find(pdbnum=%2d,a=%d nrecord=%3d) at recnum=%d (%s) anum=%d\n",
   pdbnum, a, nrecord, recnum, Rec_line[recnum], atomnumber[recnum]);
#endif /* DEBUG */
		    if (recnum<0) {
                        char  error_message[ LINE_LEN ];
			prStr(error_message, 
			"ERROR: line %d:\n%s\ntorsion %d, no atom numbered %d found in %s\n",
			 (i+1),Rec_line[i], (ntor+1), pdbnum, (a==0)?"ligand file":"flexres file");
			stop(error_message);
			}
                    tlist[ ntor ][ (a==0)?ATM1:ATM2]= atomnumber[recnum];
#ifdef DEBUG
fprintf(logFile, "    tlist[ntor=%d][a=%d %d] = %d \n", ntor,a,((a==0)?ATM1:ATM2),tlist[ntor][(a==0)?ATM1:ATM2]);
#endif /* DEBUG */
		}

		// check for degenerate torsion selection
                if ( tlist[ ntor ][ ATM2 ] == tlist[ ntor ][ ATM1 ]) {
                    char  error_message[ LINE_LEN ];
                    prStr( error_message, "ERROR: line %d:\n%sThe two atoms defining torsion %d are the same! (%d and %d)", (i+1), Rec_line[ i ], (ntor+1),
 tlist[ ntor ][ ATM1 ], tlist[ ntor ][ ATM2 ] );
                    stop( error_message );
                } /* endif */
		/* convert from 1-origin to internal 0-origin */
                //MP --tlist[ ntor ][ ATM1 ];
                //MP --tlist[ ntor ][ ATM2 ];
                nbranches = 0;

#ifdef DEBUG
                C = 'B';
                PrintDebugTors;
                PrintDebugTors2;
                pr( logFile, "]\n" );
#endif /* DEBUG */

		/* mark atoms within this branch and any internally nested branches.
		 * Note the atoms will be seen once for each branch they
		 * are within, from the "i" loop.
		 */
                for ( j = (i+1); j < nrecord; j++) {
		    int keyword;
#ifdef DEBUG
                    C = 'b';
                    PrintDebugTors;
                    PrintDebugTors2;
                    pr( logFile, "]\n" );
#endif /* DEBUG */
		    keyword= parse_PDBQT_line(Rec_line[j]);

                    if (keyword==PDBQ_ENDBRANCH && nbranches == 0)  break;
                    if (keyword==PDBQ_ENDBRANCH && nbranches != 0) --nbranches;
                    if (keyword==PDBQ_BRANCH) ++nbranches;
                    if (keyword==PDBQ_ATOM||keyword==PDBQ_HETATM) {
                        tlist[ ntor ][ tlist[ ntor ][ NUM_ATM_MOVED ] + 3 ] = atomnumber[ j ];
                        ++tlist[ ntor ][ NUM_ATM_MOVED ];
                    } /* endif */
                } /* j */
                ++ntor;

#ifdef DEBUG
                PrintDebugTors;
                PrintDebugTors2;
                pr( logFile, "]\n" );
#endif /* DEBUG */

                break;

    /*____________________________________________________________*/
            case PDBQ_CONSTRAINT:

                sscanf(Rec_line[ i ],"%*s %d %d " FDFMT2, P_atomC1, P_atomC2, &lower, &upper);

                *P_B_constrain = TRUE;

                upper = fabs( (double)upper );
                lower = fabs( (double)lower );

		if(outlev>=LOGLIGREAD)
                pr( logFile, "Constrain the distance between atom %d and atom %d to be within %.3f and %.3f Angstroms.\n\n", *P_atomC1, *P_atomC2, lower, upper);

                if (lower > upper) {
                    pr( logFile, "WARNING!  The lower distance constraint bound was larger than the upper bound. I will switch these around.\n");
                    temp = upper;
                    upper = lower;
                    lower = temp;

                } else if (lower == upper) {
                    pr( logFile, "WARNING!  The lower distance constraint bound is the same as the upper bound.\n");
                    upper += 0.01;
                }
                *P_sqlower = lower * lower;
                *P_squpper = upper * upper;
                break;

    /*____________________________________________________________*/
            case PDBQ_BEGIN_RES:
                found_new_res = true;
                // if this is the first BEGIN_RES tag, then set the number of torsions in the ligand
                if (!found_first_res) {
                    *P_ntor_ligand = ntor;
                }
                found_first_res=true;
                natoms_in_res = 0; /* reset number of atoms in this residue */
                break;

    /*____________________________________________________________*/
            case PDBQ_END_RES:
                found_new_res = false;
                nres++;
		if(outlev>=LOGLIGREAD)
                pr(logFile, "Residue number %d has %d moving atoms.\n\n", nres, natoms_in_res-2);
                break;

    /*____________________________________________________________*/
            case PDBQ_TORSDOF:
                sscanf(Rec_line[i], "%*s %d", P_ntorsdof);
		if(outlev>=LOGLIGREAD)
                pr( logFile, "\nTORSDOF record detected: number of torsional degress of freedom has been set to %d.\n", *P_ntorsdof );
                break;

    /*____________________________________________________________*/
            default:
                break;
    /*____________________________________________________________*/
        } /* switch -- finished parsing this line of PDBQT file*/
    } /* i --- do next record in PDBQT file... */

    /*
    \   Sort Torsion list on number of atoms moved,
     \______________________________________________________________
    */
    // Checked for above, as well as in readPDBQT, but this is extra insurance
    if (ntor > MAX_TORS) {
        char  error_message[ LINE_LEN ];
        prStr( error_message, "ERROR: Too many torsions have been found (i.e. %d); maximum allowed is %d.\n Either: change the \"#define MAX_TORS\" line in constants.h\n Or:     edit \"%s\" to reduce the number of torsions defined.", (ntor+1), MAX_TORS, smFileName );
        stop( error_message );
    } else {
        *P_ntor = ntor;
    }
    //if there are no flexible residues, still need to set P_ntor_ligand
    if (!found_first_res) *P_ntor_ligand = ntor;


    Boole B_atom_number_OK = TRUE;

    // M Pique - really need to check only ATM1 & ATM2 for PDB number existence
    //  since the range is determined by contiguous atom serial number range
    for (itor=0; itor<ntor; itor++ ) {
	char error_msg[LINE_LEN];
        B_atom_number_OK &= check_atomnumber_ok( tlist[ itor ][ ATM1 ] );
        B_atom_number_OK &= check_atomnumber_ok( tlist[ itor ][ ATM2 ] );
/** MPique in progress
        if (B_atom_number_OK) for (int i=0;  i < tlist[ itor ][ NUM_ATM_MOVED ]; i++ ) {
                B_atom_number_OK &= check_atomnumber_ok( tlist[ itor ][ 3+i ] );
            }
***/
        if (B_atom_number_OK) continue;
        prStr(error_msg, "%s: ERROR:  Torsion number %d between %s atom %d and atom %d has one or more atoms (out of %d atoms) that are out of range.\n\n",
 programname, itor+1, 
 ((itor<=*P_ntor_ligand)?"ligand":"flexres"),
  1+tlist[itor][ATM1],
  1+tlist[itor][ATM2],
  tlist[itor][NUM_ATM_MOVED] );
	stop(error_msg); // exits
    }

    itor = 0;
    for (i=0;  i<MAX_ATOMS;  i++) {
        for ( j=0; j<ntor; j++ ) {
            if (tlist[ j ][ NUM_ATM_MOVED ] == i) {
                for (k=0;  k<MAX_ATOMS;  k++) {
                    tlistsort[ itor ][ k ] = tlist[ j ][ k ];
                }
                ++itor;
            }
        }
    }
    for ( i=0; i<ntor; i++ ) {
        for (j=0;  j<MAX_ATOMS;  j++) {
            tlist[ i ][ j ] = tlistsort[ i ][ j ];
        }
    }

    /* fill in last entry in tlist - root atoms. 
     *  root atoms here mean only within the "true ligand".
     *  Entries for ATM1 and ATM2 are meaningless for the root.
     * M Pique spring 2012 
     */
    { // local block
    Boole B_is_root_atom[MAX_ATOMS];
    int nrootatoms; // assumed in root until proven otherwise
    // fill boolean array with TRUE up to true_ligand_atoms
    for(int i=0;i<true_ligand_atoms;i++)  B_is_root_atom[i]= TRUE;
    // set FALSE any atoms that participate in a torsion (including pivot atom)
    for(int tor=0;tor<ntor;tor++) {
        for(int i=0;i<tlist[tor][NUM_ATM_MOVED];i++) 
	  B_is_root_atom[ tlist[tor][i+NUM_ATM_MOVED+1] ] = FALSE;
	B_is_root_atom[ tlist[tor][ATM2] ] = FALSE; // pivot atom for tor
    }
    // count atoms still in root
    nrootatoms=0;
    for(int i=0;i<true_ligand_atoms;i++) if(B_is_root_atom[i]) {
    	tlist[ntor][NUM_ATM_MOVED+nrootatoms+1] = i;
	++nrootatoms;
	}
    tlist[ntor][ATM1] = tlist[ntor][ATM2] = 0; // meaningless
    tlist[ntor][NUM_ATM_MOVED] = nrootatoms;
    } // local block

    if (ntor > 0 ) {
        pr( logFile, "\n\nNumber of Rotatable Bonds in Small Molecule =\t%d torsions\n", ntor);
	if(outlev>=LOGLIGREAD)  {
        pr( logFile, "\n\nTORSION TREE\n____________\n\nSorted in order of increasing number of atoms moved:\n\n" );
     
        pr( logFile, "Torsion                    #\n" );
        pr( logFile, " #  Atom1--Atom2 Moved List of Atoms Moved\n" );
        pr( logFile, "___ ____________ _____ ________________________________________________________\n");
        for ( j=0; j<ntor; j++ ) {
           int   imax ;
            pr( logFile, "%2d  %5s--%-5s  %3d  ", j+1, pdbaname[ tlist[ j ][ ATM1 ] ], pdbaname[ tlist[ j ][ ATM2 ] ], tlist[ j ][ NUM_ATM_MOVED ] );
            imax = tlist[ j ][ NUM_ATM_MOVED ] + 2;
            for ( i = 3; i <= imax; i++ ) {
                pr( logFile, "%s%c", pdbaname[ tlist[ j ][ i ]], (i<imax)?',':'.' );
            }
            pr( logFile, "\n" );
        }
        pr( logFile, "\n" );
	} // end outlev
    } else { 
	if(outlev>=LOGBASIC)
        pr( logFile, "\n*** No Rotatable Bonds detected in Small Molecule. ***\n\n" );
    }
}

static
int check_atomnumber( const int number, const int *atomnumber, const int nrecord, 
const int natoms, FILE *logFile) {
	/* Check that the 0-origin serial number is within range and exists.
	 * Note that the number passed is zero-origin 
	 * as is the "atomnum" array content (which is -1 for non ATOM/HETATM)
 	 *
         *  Error return 0 for:
         *  negative
	 *  too large
	 *  entry not found in any record
         *  duplicate entry found 
	 */
	int i;
	int found_count;
	if(number<0) {
		pr(logFile, "ERROR: Atom number %d is < 0.\n", 1+number);
		return 0;
		}
	if(number>natoms) {
		pr(logFile, "ERROR: Atom number %d is > natoms (%d).\n", 1+number, natoms);
		return 0;
		}
	found_count=0;
	for(i=0;i<nrecord;i++) if(atomnumber[i]==number) found_count++;
	if(found_count==1) return 1; /* OK */
	if(found_count==0) {
		pr(logFile, "ERROR: Atom number %d not found in PDBQT file.\n", 1+number);
		}
	if(found_count>1) {
		pr(logFile, "ERROR: Atom number %d found multiple times (%d)in PDBQT file.\n", 1+number, found_count);
		}
	return 0;
	}

static int pdbatomnum_find_record( int pdbn, int pdbatomnumber[], int nrecord ) {
	/* return record number of atom (in passed name space) that has PDB
	 * atom number pdbn, or -1 if not found 
	 */
		
	for(int i=0;i<nrecord;i++) if(pdbn==pdbatomnumber[i]) return i;
	return -1;
	}

	
/* EOF */
