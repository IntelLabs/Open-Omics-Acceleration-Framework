/*

 $Id: weedbonds.cc,v 1.24 2014/06/12 01:44:08 mp Exp $

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
#include <stdlib.h>
#include "weedbonds.h"



void weedbonds( const int natom,
                const char pdbaname[MAX_ATOMS][5],
                const int rigid_piece[MAX_ATOMS],
                const int ntor,
                const int tlist[MAX_TORS][MAX_ATOMS],
      /* not const */ int nbmatrix[MAX_ATOMS][MAX_ATOMS],
      /* not const */ int *const Addr_Nnb,
      /* not const */ NonbondParam *nonbondlist,
      /* not const */ int Nnb_array[3],
		const int true_ligand_atoms,
                const int type[MAX_ATOMS],
		const int debug,
                const int outlev,
		FILE *logFile)

{
    int a11=0;
    int a12=0;
    int a21=0;
    int a22=0;
    int p11 = 0;
    int p12 = 0;
    int p21 = 0;
    int p22 = 0;
    int p = 0;
    int Nnb = 0;

    register int i = 0;
    register int j = 0;
    register int k = 0;

    char error_msg[LINE_LEN];


/*___________________________________________________________________________
|    ENDBRANCH---TORS---BRANCH---R O O T---BRANCH---ENDBRANCH                |
|                                  /              \                          |
|                                BRANCH            BRANCH--TORS---ENDBRANCH  |
|                                /                  \                        |
|                               ENDBRANCH            ENDBRANCH               |
|____________________________________________________________________________|
|  Eliminate all rigidly bonded atoms:                                       |
|                                     * those atoms which are at the ROOT;   |
|                                     * atoms between TORS and BRANCH;       |
|                                     * atoms between BRANCH and ENDBRANCH.  |
|  This is necessary for internal energy calculations.                       |
|____________________________________________________________________________|
| Weed out bonds in rigid pieces,                                            |
|____________________________________________________________________________|
*/        


    for (j = 0;  j < natom;  j++) {
        for (i = 0;  i < natom;  i++) {
            // Is atom "i" in the same rigid piece as atom "j"?
            if (rigid_piece[i] == rigid_piece[j]) {
                // Set the entry for atoms "i" and "j" in the
                //    nonbond matrix to 0 ;
                // 
                // Later on, we will not calculate the interaction energy
                //    between atoms "i" and "j"
                nbmatrix[j][i] = nbmatrix[i][j] = 0;
            }
        } // i
    } // j

    /* 
    \   Weed out bonds across torsions,
     \______________________________________________________________
    */
    // Loop over all "ntor" torsions, "i"
    for (i=0; i<ntor; i++) {
        nbmatrix[ tlist[i][ATM2] ][ tlist[i][ATM1] ] = 0;

        if (debug>0) {
            pr(logFile, "DEBUG: Inside weedbonds now; torsion i=%d, outlev is set to %d.\n\n", i, outlev);
        }
#ifdef DEBUG
    pr( logFile, "__2__ i=%-2d: nbmatrix[tlist[%-2d][ATM2]=%-2d][tlist[%-2d][ATM1]=%-2d]=%d\n",i,i,tlist[i][ATM2],i,tlist[i][ATM1],nbmatrix[tlist[i][ATM2]][tlist[i][ATM1]]);
#endif // DEBUG

    } // i

    /* 
    \  Weed out bonds from atoms directly connected to rigid pieces,
     \_ we think these are 1-3 interactions mp+rh, 10-2008______________________
    */
    for (i=0; i<ntor; i++) {

        a11 = tlist[i][ATM1];
        a21 = tlist[i][ATM2];
        p11 = rigid_piece[a11];
        p21 = rigid_piece[a21];

        for (j=0; j<ntor; j++) {

            a12 = tlist[j][ATM1];
            a22 = tlist[j][ATM2];
            p12 = rigid_piece[a12];
            p22 = rigid_piece[a22];

            if (p11 == p12)  { nbmatrix[ a22 ][ a21 ] = nbmatrix[ a21 ][ a22 ] = 0; }
            if (p11 == p22)  { nbmatrix[ a12 ][ a21 ] = nbmatrix[ a21 ][ a12 ] = 0; }
            if (p21 == p12)  { nbmatrix[ a22 ][ a11 ] = nbmatrix[ a11 ][ a22 ] = 0; }
            if (p21 == p22)  { nbmatrix[ a12 ][ a11 ] = nbmatrix[ a11 ][ a12 ] = 0; }

        } // j

        for (k = 0;  k < natom;  k++) {
            p = rigid_piece[k];
            if (p11 == p)  { nbmatrix[ k ][ a21 ] = nbmatrix[ a21 ][ k ] = 0; }
            if (p21 == p)  { nbmatrix[ k ][ a11 ] = nbmatrix[ a11 ][ k ] = 0; }
        } // k
    } // i

    /*
    for ( i = 0;  i < (natom-1);  i++ ) {
        for ( j = i+1;  j < natom;  j++ ) {
            if ( ((nbmatrix[i][j] == 1) && (nbmatrix[j][i] == 1)) || ((nbmatrix[i][j] == 4) && (nbmatrix[j][i] == 4)) ) {
                nonbondlist[Nnb][ATM1] = i;
                nonbondlist[Nnb][ATM2] = j;
                nonbondlist[Nnb][TYPE1] = type[i];
                nonbondlist[Nnb][TYPE2] = type[j];
                nonbondlist[Nnb][NBTYPE] = nbmatrix[i][j];
                ++Nnb;
            } else if ( (nbmatrix[i][j] != 0 && nbmatrix[j][i] == 0) || (nbmatrix[i][j] == 0 && nbmatrix[j][i] != 0) ) {
                prStr( error_msg, "BUG: ASSYMMETRY detected in Non-Bond Matrix at %d,%d\n", i,j);
		stop(error_msg);
            }
        } // j
    } // i
    */

    // intramolecular non-bonds for ligand
    for ( i = 0;  i < (true_ligand_atoms-1);  i++ ) {
        for ( j = i+1;  j < true_ligand_atoms;  j++ ) {
            if ( (nbmatrix[i][j] == 1 && nbmatrix[j][i] == 1) || (nbmatrix[i][j] == 4 && nbmatrix[j][i] == 4) ) {
                nonbondlist[Nnb].a1 = i;
                nonbondlist[Nnb].a2 = j;
                nonbondlist[Nnb].t1 = type[i];
                nonbondlist[Nnb].t2 = type[j];
                nonbondlist[Nnb].nonbond_type = nbmatrix[i][j];
                ++Nnb;
            } else if ( (nbmatrix[i][j] != 0 && nbmatrix[j][i] == 0) || (nbmatrix[i][j] == 0 && nbmatrix[j][i] != 0) ) {
                prStr( error_msg, "BUG: ASSYMMETRY detected in Non-Bond Matrix at %d,%d\n", i,j);
		stop(error_msg);
            }
        } // j
    } // i
    Nnb_array[INTRA_LIGAND] = Nnb;

    // intermolecular non-bonds for ligand,receptor
    for ( i = 0;  i < true_ligand_atoms;  i++ ) {
        for ( j = true_ligand_atoms;  j < natom;  j++ ) {
            if ( (nbmatrix[i][j] == 1 && nbmatrix[j][i] == 1) || (nbmatrix[i][j] == 4 && nbmatrix[j][i] == 4) ) {
                nonbondlist[Nnb].a1 = i;
                nonbondlist[Nnb].a2 = j;
                nonbondlist[Nnb].t1 = type[i];
                nonbondlist[Nnb].t2 = type[j];
                nonbondlist[Nnb].nonbond_type = nbmatrix[i][j];
                ++Nnb;
            } else if ( (nbmatrix[i][j] != 0 && nbmatrix[j][i] == 0) || (nbmatrix[i][j] == 0 && nbmatrix[j][i] != 0) ) {
                prStr( error_msg, "BUG: ASSYMMETRY detected in Non-Bond Matrix at %d,%d\n", i,j);
		stop(error_msg);
            }
        } // j
    } // i
    Nnb_array[INTER] = Nnb;

    // intramolecular non-bonds for receptor
    for ( i = true_ligand_atoms;  i < natom-1;  i++ ) {
        for ( j = i+1;  j < natom;  j++ ) {
            if ( (nbmatrix[i][j] == 1 && nbmatrix[j][i] == 1) || (nbmatrix[i][j] == 4 && nbmatrix[j][i] == 4) ) {
                nonbondlist[Nnb].a1 = i;
                nonbondlist[Nnb].a2 = j;
                nonbondlist[Nnb].t1 = type[i];
                nonbondlist[Nnb].t2 = type[j];
                nonbondlist[Nnb].nonbond_type = nbmatrix[i][j];
                ++Nnb;
            } else if ( (nbmatrix[i][j] != 0 && nbmatrix[j][i] == 0) || (nbmatrix[i][j] == 0 && nbmatrix[j][i] != 0) ) {
                prStr(error_msg, "BUG: ASSYMMETRY detected in Non-Bond Matrix at %d,%d\n", i,j);
		stop(error_msg);
            }
        } // j
    } // i
    Nnb_array[INTRA_RECEPTOR] = Nnb;

    if (Nnb > MAX_NONBONDS) {
        prStr( error_msg, "too many non-bonded interactions (%d) in small molecule\n\t(increase MAX_NONBONDS from %d).", Nnb, MAX_NONBONDS );
        stop( error_msg );
        exit( EXIT_FAILURE );
    } else {
        *Addr_Nnb = Nnb;
    }

    flushLog;

} // weedbonds


void print_nonbonds(
                const int natom,
                const char pdbaname[MAX_ATOMS][5],
                const int rigid_piece[MAX_ATOMS],
                const int ntor,
                const int tlist[MAX_TORS][MAX_ATOMS],
      /* not const */ int nbmatrix[MAX_ATOMS][MAX_ATOMS],
                const int Nnb,
                const NonbondParam *const nonbondlist,
                const int type[MAX_ATOMS],
                const int outlev,
		FILE *logFile)

{
    register int i = 0;
    register int j = 0;
    register int k = 0;
    int i_atmnum = 0;
    int j_atmnum = 0;
    int Nnbonds[MAX_ATOMS];
    const static int OUTNUMATM = 10;
    int repflag = FALSE;
    int n_a = 0;

    // Set the number of nonbonds of each atom "i" to 0
    for (i = 0;  i < natom;  i++) {
        Nnbonds[i] = 0;
    }

    if (outlev >= LOGLIGREAD)  {
        // Print out the matrix of non-bonded interactions
        if (ntor > 0) {
            pr( logFile, "\n\nMatrix of Non-Bonded Interactions:\n" );
            pr( logFile, "__________________________________\n\n" );
            pr( logFile, "Key:\nX = non-bonded interaction\n" );
            pr( logFile, "_ = 1,2 or 1,3 interaction\n\n" );
            pr( logFile, "\nAtom: ID: " ); for (j = 0;  j < natom;  j++)   fprintf( logFile, "%2d", (1+j) );
            pr( logFile, "\n_____ ___ " ); for (j = 0;  j < natom;  j++)   fprintf( logFile, "__" );
            for (j = 0;  j < natom;  j++) {
                pr( logFile, "\n%4s  %2d  ", pdbaname[j], 1+j );
                for (i = 0;  i < natom;  i++) {
                    pr( logFile, "|%c", (nbmatrix[j][i])?'X':'_' );
                } // i
            } // j
            pr( logFile, "\n\n" );
            flushLog;
        } //  endif 
    }

#ifdef DEBUG
        for (i = 0;  i < Nnb;  i++) {
            pr( logFile,"> nonbondlist[%2d][0,1] = %2d,%2d\n", i,nonbondlist[i].a1,nonbondlist[i].a2 );
        } //  i 
#endif // DEBUG

    for (i = 0;  i < Nnb;  i++) {

        i_atmnum = nonbondlist[i].a1;
        j_atmnum = nonbondlist[i].a2;

#ifdef DEBUG
            pr( logFile, "* Assigning nbmatrix[%d][%d] to %d *\n", Nnbonds[i_atmnum], i_atmnum, j_atmnum );
#endif // DEBUG

        nbmatrix[Nnbonds[i_atmnum]][i_atmnum] = j_atmnum;
        
        /* NOTE: re-utilizes the nbmatrix array;
        \        nbmatrix is `corrupted' after this...
         \ 
          \ Normally, nbmatrix contains 0s and 1s, but this
           \ assigns atom-numbers, which can be >1.
          */

#ifdef DEBUG
            pr( logFile, ">>--> i = %d, j_atmnum = %d, nbmatrix[ %d ][ %d ] = %d\n",
            i, j_atmnum, Nnbonds[i_atmnum],i_atmnum, nbmatrix[Nnbonds[i_atmnum]][i_atmnum] );
#endif // DEBUG

        ++Nnbonds[i_atmnum];
    } //  i 

#ifdef DEBUG
    if (ntor > 0) {
            pr( logFile, "\nCORRUPTED Matrix of Non-Bonded Interactions:\n" );
            pr( logFile, "____________________________________________\n\n" );
            pr( logFile, "\nAtom: ID: " ); for (j = 0;  j < natom;  j++)   fprintf( logFile, "%2d", j );
            pr( logFile, "\n_____ ___ " ); for (j = 0;  j < natom;  j++)   fprintf( logFile, "__" );
            for (i = 0;  i < natom;  i++) {
                pr( logFile, "\n%4s  %2d  ", pdbaname[i], i );
                for (j = 0;  j < natom;  j++) {
                        pr( logFile, "%2d", nbmatrix[j][i] );
                } // j
            } // i
            pr( logFile, "\n\n" );
            flushLog;
    } // endif
#endif // DEBUG

    if (outlev >= LOGLIGREAD) {
        // Print out a list of internal non-bonded interactions
        if (ntor > 0) {
            pr( logFile, "\n\nList of Internal Non-Bonded Interactions:\n" );
            pr( logFile, "_________________________________________\n\n" );
            pr( logFile, "First Second Non-Bonded\n" );
            pr( logFile, "Atom  Atom(s)\n" );
            pr( logFile, "_____ " );
        } // endif

        for (i = 0;  i < natom;  i++) {
            if (Nnbonds[i] != 0) {
                for ( j = 0;  j < OUTNUMATM;  j++ ) {
                    pr( logFile, "_____" );
                } // j
                break;
            } // endif
        } // i

        if (ntor > 0) {
            pr( logFile, "_\n" );
        }

        for (i = 0;  i < natom;  i++) {
            // loop over all atoms, "i"
            if (Nnbonds[i] != 0) {
                // this atom "i" has some nonbonds
                repflag = FALSE; // set the repeat(?) flag to FALSE.
                n_a = 0;  // number of atoms output on this line so far.
                pr( logFile, "%4s: ", pdbaname[i] );
                for ( j = 0; j < Nnbonds[i]; j++) {
                    // loop over all nonbonds "j" for this atom
                    ++n_a;
                    if (n_a >= OUTNUMATM) {
                        n_a = 0;
                        pr( logFile, "\n      " );
                    }
                    pr( logFile, "%s", pdbaname[ nbmatrix[j][i] ] );
                    for ( k = 0; k < j; k++ ) {
                        // loop over all "k" less than "j", the current nonbond
                        if (nbmatrix[k][i] == nbmatrix[j][i]) {
                            repflag = TRUE;
                            break;
                        }
                    } //  k 
                    if (repflag && (k != j)) {
                        pr( logFile,"(ERROR! %d & %d)", j, k );
                    }
                    pr( logFile, "%s", (j == (Nnbonds[i]-1))?".":", " );
                } //  j 
                pr( logFile,"\n\n");
            } //  endif 
        } //  i 
    } // outlev > -1

    if( outlev >= LOGLIGREAD ){
    pr( logFile, "\nInternal Non-bonded Interactions before,\t%d\n", (natom+1)*natom/2);
    pr( logFile, "                       and after weeding =\t%d\n\n", Nnb);

    flushLog;
    }
}
// EOF
