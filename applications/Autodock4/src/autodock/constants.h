/*

 $Id: constants.h,v 1.44 2014/06/12 01:44:07 mp Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.
 AutoGrid is a Trade Mark of The Scripps Research Institute.

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

#include "autocomm.h" /* PI, TWOPI. */
#include "dpftoken.h"
#include "partokens.h"


/******************************************************************************
 *      Name: constants.h                                                     *
 *  Function: Automated Docking of Small Molecule to Macromolecule            *
 *              Header file for Autodock modules.                             *
 *Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
 *----------------------------------------------------------------------------*
 *    Author: Garrett Matthew Morris*                                         *
 *                                                                            *
 *              The Scripps Research Institute                                *
 *              Department of Molecular Biology, MB5                          *
 *              10666 North Torrey Pines Road                                 *
 *              La Jolla, CA 92037.                                           *
 *              e-mail: garrett@scripps.edu                                   *
 *                                                                            *
 *      Date: 10/04/94                                                        *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 *          DSG     Quaternion rotations                                      *
 *          DSG     Generates torsions from annotated pdb list                *
 *          DSG     Generates internal energies                               *
 *          DSG     Performs a limited cluster analysis of conformations      *
 * 07/05/92 GMM     C translation                                             *
 * 14/05/92 GMM     Time-dependent seed in random-number generation           *
 * 29/10/92 GMM     Application Visualization System (AVS) readable grid      *
 *                    display file input.                                     *
 *                    [AVS is a trademark of Stardent Computer Inc.]          *
 * 06/11/92 GMM     Command line parsing, using Bruce S. Duncan's "setflags". *
 * 06/11/92 GMM     New command mode, allowing communication between Autodock *
 *                    and other, invoking programs.                           *
 ******************************************************************************/

// note: this file is shared by AutoDock and AutoGrid

#ifndef CONSTANTS
#define CONSTANTS

#define REJECT       0
#define ACCEPT       1

#define U            0        /* x-coordinate in trilinear interpolation. */
#define V            1        /* y-coordinate in trilinear interpolation. */
#define W            2        /* z-coordinate in trilinear interpolation. */
#define QX           3        /* quaternion x-coordinate */
#define QY           4        /* quaternion y-coordinate */
#define QZ           5        /* quaternion z-coordinate */
#define QW           6        /* quaternion w-rotation */
#define QUAT         7        /* No. of translation & rotations in quaternion */

#define NTRN         3        /* Number of translation axes */
#define NQTN         4        /* Number of quaternion components */

#define LOWER 0               /* Lower limit of torsion constraint. */
#define UPPER 1               /* Upper limit of torsion constraint. */

#define NLIG         4        /* Max. num. of ligands to be docked */
#define NLIGTOR      32       /* Max. num. of torsions in each ligand */

#define NRES         32       /* Max. num. of flexible residue */
#define NRESTOR      5        /* Max. num. of torsions in a residue */

#define NMOL         36       /* NMOL = NLIG + NRES, always!!! */

#define WORDLEN      9        /* Length of coordinates in characters. */

#ifdef USE_XCODE
// Xcode-gmm
#define MAX_RUNS     16      /* Maximum number of runs. */
#define MAX_ATOMS    128     /* Maximum number of atoms in Small Molecule. */
#define MAX_RECORDS  128     /* Maximum number of PDBQ records =~1.5*MAX_ATOMS */
#define MAX_NONBONDS 16384   /* 1048576   Maximum number of non-bonds in Small Molecule.*/
#else
#define MAX_RUNS     2000      /* Maximum number of runs. Used in Tests */
#define MAX_ATOMS    2048     /* Maximum number of atoms in Small Molecule. */
#define MAX_NBONDS    16     /* Maximum number of bonds per atom, in readPDBQT  */
#define MAX_RECORDS  2048     /* Maximum number of PDBQ records =~1.5*MAX_ATOMS */
#define MAX_NONBONDS 524288 /* 1048576   Maximum number of non-bonds in Small Molecule.*/
#endif

#define MAX_TREES    8        /* Maximum number of torsion trees */
#define MAX_TORS     32       /* Maximum number of torsions in Small Molecule. */
#define MAX_TORS_IN_ROTAMER 8 /* Maximum number of torsions in rotamer. */
#define MAX_TOR_CON  8        /* Maximum number of constraints per torsion. */

/* Note, torsion constraints will allocate (MAX_TORS *MAX_TOR_CON *2) Reals */

#define NTORDIVS     256      /* Integer # angular divisions in torsion profiles*/
#define NTORDIVS_2   128      /* Integer half of NTORDIVS */
#define ONE_RAD_IN_DIVS 40.74366543  /* 1 radian = NTORDIVS_2 / PI divisions */

/* #define ONE_RAD_IN_DIVS 81.48733086 512,256,1 radian = NTORDIVS / (2*PI) divisions */

/* Convert angle to index */
#define RadiansToDivs(radians) ((int)((radians)*ONE_RAD_IN_DIVS) + NTORDIVS_2)

/* Convert index to angle */
#define DivsToRadians(divs) (((Real)((divs) - NTORDIVS_2))/ONE_RAD_IN_DIVS)

#define TORBARMAX    65535    /* Torsion barrier-energy maximum value */
#define DEFHWDTH     10.      /* Default half-width (deg)for torsion-constraints*/

/*
 * Number of &tor0 entries in TOR_ARG_LIST must equal MAX_TORS...
 * It's a kludge, I know, but you only need to get this line right once.
 */
#define TOR_ARG_LIST        &sInit.tor[0], &sInit.tor[1], &sInit.tor[2], &sInit.tor[3], &sInit.tor[4], &sInit.tor[5], &sInit.tor[6], &sInit.tor[7], &sInit.tor[8], &sInit.tor[9], &sInit.tor[10], &sInit.tor[11], &sInit.tor[12], &sInit.tor[13], &sInit.tor[14], &sInit.tor[15], &sInit.tor[16], &sInit.tor[17], &sInit.tor[18], &sInit.tor[19], &sInit.tor[20], &sInit.tor[21], &sInit.tor[22], &sInit.tor[23], &sInit.tor[24], &sInit.tor[25], &sInit.tor[26], &sInit.tor[27], &sInit.tor[28], &sInit.tor[29], &sInit.tor[30], &sInit.tor[31]


#define ENERGY_CUTOFF 500.    /* Arbitrary intermolecular cutoff, above 
                                 which intramolecular energy is not calculated. 
				 (appears unused in AD code - perhaps great idea - MP 2012)
				 */
#define HI_NRG_JUMP_FACTOR 2. /* Scale up the range of random jumps by this when the 
                                 last energy was higher than ENERGY_CUTOFF. */

#ifdef USE_8A_NBCUTOFF

#define NBC         8.00      /* Non-bonded cutoff for internal energy calc./Ang*/ 
#define NEINT    2048         /* Number of values in internal energy table */
#define A_DIV     100.00      /* Used in distance look-up table. */

#else

#define NBC        64.00      /* Non-bonded cutoff for internal energy calc./Ang*/
#define NEINT  131072         /* Number of values in internal energy table */
#define A_DIV     100.00      /* Used in distance look-up table. */

#endif

#define NBC2     (NBC*NBC)    /* NBC^2, units: Angstrom^2 */
#define NEINT_1 (NEINT - 1)   /* index of last entry in internal energy table */
#define INV_A_DIV  (1./A_DIV)      /* Used in distance look-up table. i.e. every 1/100-th of an Angstrom */
#define INT_SQA_DIV   (NEINT/(int)NBC2)      /* Xcode-gmm */
#define SQA_DIV    (NEINT/NBC2)      /* Used in square-distance look-up table. */
#define INV_SQA_DIV (1./SQA_DIV) /* 0.03125   INV_SQA_DIV  =  1/SQA_DIV  =  NBC2 / NEINT   */


#define NDIEL 16384           /* Number of dielectric and desolvation values in lookup table.
                                 NDIEL is bigger than NEINT because electrostatic interactions are much
                                 longer-range than van der Waals interactions. */
#define NDIEL_1 (NDIEL - 1)   /* The last valid index in dielectric and desolv lookup tables, NDIEL minus 1 */

/*
 * Alternate Scheme:-              (Uses less memory; for smaller ligands, < 8.0 Ang.)
 *
 * #define NEINT       8192        Number of values in internal energy table
 * #define NEINT_1 (NEINT - 1)     index of last entry in internal energy table
 * #define A_DIV    128.00         Used in look-up table.
 * #define INV_A_DIV  0.0078125    INV_A_DIV = NBC2 / NEINT, i.e. every 7.8125/1000-th of an Angstrom
 * #define NBC        8.00         Non-bonded cutoff for internal energy calc./Ang
 * #define NBC2      64.00         NBC^2, units: Angstrom^2
 */

#define ATM1            0     /* Index for the first atom in a bond or non-bond */
#define ATM2            1     /* Index for the second atom in a bond or non-bond */

#define NUM_ATM_MOVED   2     /* Index for the number of atoms moved by a torsion-rotation */

#define TYPE1			2     /* Index for the atom type of first atom - used in nonbondlist  Xcode-gmm */
#define TYPE2			3     /* Index for the atom type of second atom - used in nonbondlist Xcode-gmm */
#define NBTYPE			4     /* Index for the type of nonbond */
#define MAX_NBDATA	    5     /* 5 elements of data for each nonbond: atm1, atm2, type1, type2 & whether 1,4 or other */

#define TINYDELTA 0.001       /* To nudge ligand into grid... */

/* Prevent internal electrostatic energy calculation from blowing up! */

#define RMIN_ELEC 0.5         /* if atoms closer than this in Angstroms, clamp distance in internal elec. calc.*/
#define RMIN_ELEC2  (RMIN_ELEC * RMIN_ELEC) /* if atoms closer than this, clamp square 
                                               of distance in internal elec. calc.*/

/* Set reasonable limits on Rij and epsilon-ij values*/

#define RIJ_MIN 0.9
#define RIJ_MAX 6.0
#define EPSIJ_MIN 0.0
#define EPSIJ_MAX 10.0


/* Statistics output frequency in GAs */
 
#define OUTLEV0_GENS 0   /* Never output based on number of generations */
#define OUTLEV1_GENS 100 /* output every 100 generations */
#define OUTLEV2_GENS 1 /* output every generation */

/* logFile outlev handling: low to high
 * This is AutoDock only, not AutoGrid so not certain it should be here 
 * in constants.h but that's where I'm putting it for now - MPique 2012
 */

#define LOGMIN -2
#define LOGMINCLUST -1
#define LOGBASIC 0
#define LOGFORADT 1
#define LOGRUNV 2
#define LOGLIGREAD 3
#define LOGRECREAD 4
#define LOGRUNVV 5
#define LOGRUNVVV 6
#define LOGETABLES 7
#define LOGNBINTE 8   // analysis.cc nonbond internal energy table
#define LOGNBINTEV 9   // analysis.cc;nbe.cc nonbond internal energy table verbose
const struct {
	int value;
	char *key;
	} outlev_lookup[] = { 
 {LOGMIN, "min"},
 {LOGMINCLUST, "mincluster"},
 {LOGBASIC, "basic"},
 {LOGFORADT, "adt"},
 {LOGFORADT, "default"},
 {LOGRUNV, "runv"},
 {LOGLIGREAD, "ligread"},
 {LOGRECREAD, "recread"},
 {LOGRUNVV, "runvv"},
 {LOGRUNVVV, "runvvv"},
 {LOGETABLES, "etables"},
 {LOGNBINTE, "nbinte"}, 
 {LOGNBINTEV, "nbintev"}, 
 {0, ""} /* last */ };
#endif  /* CONSTANTS */


/*----------------------------------------------------------------------------* 
 * Macros,                                                                    * 
 *----------------------------------------------------------------------------*/

#ifndef MACROS
#define MACROS

#define equal(a,b,n) ( strncmp(a,b,(size_t)(n)) == (int)0 )

#define max(x,y)     ( ((x) > (y)) ? (x) : (y) )
#define min(x,y)     ( ((x) < (y)) ? (x) : (y) )
#define clamp(x,lowerbound)        ( ((x) < (lowerbound)) ? (lowerbound) : (x) )

// Convert from degrees to radians
#define DegreesToRadians(deg)     ( (deg) * 0.01745329252 )
// Convert from radians to degrees
#define RadiansToDegrees(rad)     ( (rad) * 57.29577951 )

// Wrap the angle supplied to be in the range 0. to 360.
#define ModDeg(a)    fmod((double)(a),(double)360.)
// Wrap the angle supplied to be in the range 0. to 2*PI
#define ModRad(a)    fmod((double)(a),(double)TWOPI)

#define WrpDeg(a)    (((a)>180.)? ((a)-360.) :(((a)<-180.)? ((a)+360.) :(a)))
#define WrpRad(a)    (((a)> PI)? ((a)-TWOPI) :(((a)< -PI)?  ((a)+TWOPI) :(a)))
// #define WrpModRad(a) (((a)> PI)? (ModRad(a)-TWOPI) :(((a)< -PI)?  (ModRad(a)+TWOPI) :(a)))
#define WrpModRad(a) WrpRad( ModRad( (a) ) )

/*
 * #define        RedFac(s0,sN,N)                expf( logf((sN)/(s0)) / ((N)-1))
 * N.B. You must compile with ANSI (-Aa on HPPA) in order to use expf 
 * and logf, otherwise Reals are automatically promoted to doubles.
 */
#define        RedFac(s0,sN,N)    exp( log((sN)/(s0)) / ((N)-1))

#define arithmetic_mean( r1, r2 ) (0.5L * ((r1) + (r2)))
#define geometric_mean( e1, e2 )  sqrt((double)(e1) * (double)(e2))

#define sq(a)                     ( (a) * (a) )
#define SQ(a)                     ( (a) * (a) )

#define sqhypotenuse(x,y,z)       ( sq(x) + sq(y) + sq(z) )
#define hypotenuse(x,y,z)         (sqrt( sqhypotenuse((x), (y), (z)) ))
#define hypotenuse_F(x,y,z)       (sqrtf( sqhypotenuse((x), (y), (z)) ))

#define sqhypotenuse4(x,y,z,w)    ((double)(sq(x) + sq(y) + sq(z) + sq(w)))
#define hypotenuse4(x,y,z,w)      (sqrt( sqhypotenuse4((x), (y), (z), (w)) ))
#define hypotenuse4_F(x,y,z,w)    (sqrtf( sqhypotenuse4((x), (y), (z), (w)) ))


/* index_to_Ang converts from an array index to a distance */
#define index_to_Ang(i)           ( ( (Real) (i) ) * INV_A_DIV )

/* Ang_to_index converts from a distance to an array index */
#define Ang_to_index(r)           ( (int) ( (r) * A_DIV ) )

/* BoundedAng_to_index converts from a distance to an array index, but never returns an index out of bounds. */
#define BoundedAng_to_index(r)    ((((int)((r)*A_DIV)) > NEINT_1) ? NEINT_1 : ((int)((r)*A_DIV))

/* index_to_SqAng converts from an array index to the square of a distance */
#define index_to_SqAng(i)         ( ( (Real) (i) ) * INV_SQA_DIV )

/* SqAng_to_index converts from the square of a distance to an array index */
#define SqAng_to_index(r2)        ( (int) ( (r2) * SQA_DIV ) )
#define SqAng_to_index_NBC2(r2)   ( (int) ( ( min( r2, NBC2 ) ) * SQA_DIV ) )

/* BoundedSqAng_to_index converts from the square of a distance to an array index, but never returns an index out of bounds. */
#define BoundedSqAng_to_index(r2) ( (((int)((r2)*SQA_DIV)) > NEINT_1) ? NEINT_1 : ((int)((r2)*SQA_DIV)) )

/* BoundedNeint never returns an index greater than (NEINT - 1). */
#define BoundedNeint(i)           (((i) > NEINT_1) ? NEINT_1 : (i))

/* BoundedNdiel never returns an index greater than (NDIEL - 1). */
#define BoundedNdiel(i)           (((i) > NDIEL_1) ? NDIEL_1 : (i))

#define is_out_grid(x,y,z) (((x)<=(xlo)) || ((x)>=(xhi)) || ((y)<=(ylo)) || ((y)>=(yhi)) || ((z)<=(zlo)) || ((z)>=(zhi))) 

#define is_out_grid_info(x,y,z) (((x)<=(info->lo[X])) || ((x)>=(info->hi[X])) || ((y)<=(info->lo[Y])) || ((y)>=(info->hi[Y])) || ((z)<=(info->lo[Z])) || ((z)>=(info->hi[Z])))

#define is_zero(x) ((x) > -APPROX_ZERO) && ((x) < APPROX_ZERO) ? 1 : 0

/*----------------------------------------------------------------------------* 
 * Random numbers,                                                            * 
 *----------------------------------------------------------------------------*/

// This is platform-independent RNG-based.  All are wrappers around ignlgi() - I think! MPique 2014
#include "ranlib.h"
#define local_random()      genunf(0., 1.)

#ifdef sgi
    #include <ieeefp.h>
    #define        ISNAN(n)        isnand(n)
#else
    #define        ISNAN(n)        isnan(n)
#endif /* !sgi */

#define random_sign        ( (local_random() < 0.5) ?  (-1.0) : (+1.0) )
#define random_pm1()    ( 2.*local_random() - 1. ) /* ..."pm"="Plus or Minus" */
#define Randpm1                random_pm1()
#define RandpmPI        (PI * Randpm1)
#define random_pm(x)    ( (x) * Randpm1 )    /* ...random num. between +/- x; */
/* Return a random number in the range a to b...  */
#define random_range(a,b)    ( (a) + (((b)-(a))*local_random()) )

#endif  /* MACROS */


/*----------------------------------------------------------------------------* 
 * Debugging,                                                                 * 
.*----------------------------------------------------------------------------*/

#ifndef DEBUG_STUFF
#define DEBUG_STUFF

/*
 * Uncomment next line, and re-make executable, for general debugging...
 */
/* #define DEBUG        1         / * Switch on debugging output */

#ifdef DEBUG

/*  For debugging torsions... */

#define    PrintDebugTors    (void)fprintf(logFile, "%-2d %-2d %4s %c %-8d %-9d %-2d %-4d      [%-2d] [ ", i, atomnumber[i],rec5,C,atomlast,nbranches,j, ntor, ntor )
#define    PrintDebugTors2   for(oo=0; oo<3+tlist[ntor][2]; oo++)  (void)fprintf(logFile, "%-2d ", tlist[ ntor ][ oo ] )

#endif /* DEBUG */

#define DPrint(expr)        (void)printf(#expr " = %.4g\n", expr)

#endif /* DEBUG_STUFF */


/*----------------------------------------------------------------------------* 
 * For system timings,                                                        * 
 *----------------------------------------------------------------------------*/

#ifndef sgi
#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 60
#endif /* CLOCKS_PER_SEC */
#endif /* sgi */

#ifndef _CONST_INT
#define _CONST_INT

#define CONST_INT const int
#define CONST_FLOAT const Real

#endif /* _CONST_INT */

#ifndef _PDB_FORMATS
#define _PDB_FORMATS


/*----------------------------------------------------------------------------* 
 * Format for output                                                          * 
 *----------------------------------------------------------------------------*/

/* 
 * PDB
 * 
 * Standard PDB v2.1 format with segID, element and charge: 
 */
/* serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, segID, element, charge */
#define FORMAT_PDB2_ATOM       "ATOM  %5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s" 

/* 
 * PDBQ
 * 
 * Note: the initial %s is the prefix
 * ""          gives normal PDBQ.
 * "DOCKED: "  gives DOCKED prefix.
 * Note: no new line at end of some of these formats:
 */

/* 
 * PDBQT format
 */
//  The numbers in parentheses are the columns that this field spans, using 0-based counting.
//     e.g. the atom_name spans from column 12 to column 15 inclusive.
//  The digit immediately below the line of '|' characters is the length of the field.
//     e.g. the atom_name field is 4 characters long.
//
// record (0-5)
// |     serial_number (6-10)
// |     |     atom_name (12-15)
// |     |     |   alternate_locator (16)
// |     |     |   |residue_name (17-19)
// |     |     |   ||   chain_ID (21)
// |     |     |   ||   |residue_number (22-25)
// |     |     |   ||   ||   insertion_code (26)
// |     |     |   ||   ||   |   x (30-37)
// |     |     |   ||   ||   |   |       y (38-45)
// |     |     |   ||   ||   |   |       |       z (46-53)
// |     |     |   ||   ||   |   |       |       |       occupancy (54-59)
// |     |     |   ||   ||   |   |       |       |       |     temperature_factor (60-65)
// |     |     |   ||   ||   |   |       |       |       |     |         charge (70-75)
// |     |     |   ||   ||   |   |       |       |       |     |         |      atom_type (77-78)
// |     |     |   ||   ||   |   |       |       |       |     |         |      |                      
// 65432154321 43211321 143211   876543218765432187654321654321654321    654321 21
// ATOM     10  HO4 PGP     1      22.065  29.222  38.002  1.00  0.00     0.210 HD
#define FORMAT_PDBQT_ATOM_RESSTR         "%sATOM  %5d%-19.19s%8.3f%8.3f%8.3f%+6.2f%+6.2f    %+6.3f %-2s%s" // prefix + 79 chars + suffix

// #define FORMAT_PDBQT_ATOM_RESNUM         "%sATOM  %5d %.8s%5d     %8.3f%8.3f%8.3f%+6.2f%+6.2f    %+6.3f %-2s"
// #define FORMAT_PDBQT_ATOM_RANKRUN_STR      "ATOM  %5d %.13s     %8.3f%8.3f%8.3f%6d%6d    %+6.2f %8.3f %-2s\n"
// #define FORMAT_PDBQT_ATOM_RUN_NUM          "ATOM  %5d %.8s%5d     %8.3f%8.3f%8.3f%6d%+6.2f    %6.3f %-2s\n"
                                                                  
#endif /* _PDB_FORMATS */

/*
 * For nonbond calculations
 */
#define INTRA_LIGAND   0 // index of last of nonbonds in intramolecular energy calculation of ligand
#define INTER          1 // index of last nonbond in intermolecular energy calculation
#define INTRA_RECEPTOR 2 // index of last of nonbonds in intramolecular energy calculation of receptor

#ifndef __THERMODYNAMIC_CONSTANTS__
#define __THERMODYNAMIC_CONSTANTS__
/*
 * Thermodynamic constants
 */

const Real RJ = 8.31441;     // in J/K/mol, Gas Constant, Atkins Phys.Chem., 2/e
const Real Rcal = 1.9871917; // in cal/K/mol, Gas Constant, RJ/4.184
const Real T0K = 273.15;     // 0 degrees Celsius, in K
const Real TK = 298.15;      // Room temperature, in K
#endif

/*
 * Vectors
 */

#define Dot_product( a, b )  \
    ( (a)[X] * (b)[X]  +  (a)[Y] * (b)[Y]  +  (a)[Z] * (b)[Z] )

#define Cross_product( n, a, b )  \
    ( (n)[X] = (a)[Y] * (b)[Z]  -  (a)[Z] * (b)[Y],  \
      (n)[Y] = (a)[Z] * (b)[X]  -  (a)[X] * (b)[Z],  \
      (n)[Z] = (a)[X] * (b)[Y]  -  (a)[Y] * (b)[X] )

#define Length( a )  sqrt( Dot_product( (a), (a) ) )

#define Magnitude( a )  Length( (a) )

#define Angle_between( a, b ) \
    ( acos( Dot_product( (a), (b) ) / ( Length( (a) ) * Length( (b) ) ) ) )

#define Add_vectors( result, a, b )  \
    ( (result)[X] = (a)[X] + (b)[X], \
      (result)[Y] = (a)[Y] + (b)[Y], \
      (result)[Z] = (a)[Z] + (b)[Z] )

#define Subtract_vectors( result, a, b )  \
    ( (result)[X] = (a)[X] - (b)[X], \
      (result)[Y] = (a)[Y] - (b)[Y], \
      (result)[Z] = (a)[Z] - (b)[Z] )

#define Print_vector( fileptr, string, vector ) \
    (void) fprintf( (fileptr), "%s  [ %.3f, %.3f, %.3f ]\n", \
                    (string), (vector)[X], (vector)[Y], (vector)[Z] )


/*
 * For unbound and bound internal energy calculations,
 * used when calling intnbtable()
 */
#define UNBOUND_CALCULATION TRUE
#define BOUND_CALCULATION FALSE

#ifndef _UNBOUND_MODEL
#define _UNBOUND_MODEL
/*
 * Set the options for unbound model
 */
enum Unbound_Model { Unbound_Default=0, Unbound_Same_As_Bound=1, Extended=2, Compact=3, User=4 };
#endif

/*----------------------------------------------------------------------------* 
 * End of file                                                                * 
 *----------------------------------------------------------------------------*/
