/*

 $Id: eval.cc,v 1.39 2014/06/12 01:44:07 mp Exp $

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

/********************************************************************
     These are the functions associated with the evaluation object.

                                rsh 9/95
********************************************************************/


#include <math.h>
#include "eval.h"
#include "stateLibrary.h"
#include "assert.h"

#include <stdio.h>
#include <string.h>

#ifdef sgi
    #include <ieeefp.h>
#endif

#ifdef sun
    #include <ieeefp.h>
#endif

/*  The chromosome is assumed to have a layout like this -

       | x | y | z | qx | qy | qz | qw | tor1 | ... | tor N |

    where:
       x is the x translation
       y is the y translation
       z is the z translation
       qx, qy, qz, qw are the components of a 4D-normalized quaternion
       tor 1, ..., tor N are the ntor torsion angles
*/

void make_state_from_rep(const Representation *const *const rep, /* not const */ State *const stateNow, const int outlev, FILE *logFile) /* not a member function */
/*
    This routine modifies the various components of stateNow to correspond
    to the chromosome.  
*/
{
   register int i;

#ifdef DEBUG
   (void)fprintf(logFile, "eval.cc/make_state_from_rep(Representation **rep, State *stateNow)\n");
#endif /* DEBUG */

   //  Do the translations
   assert( !ISNAN( rep[0]->gene(0).real ) );
   stateNow->T.x = rep[0]->gene(0).real;
   assert( !ISNAN( rep[1]->gene(0).real ) );
   stateNow->T.y = rep[1]->gene(0).real;
   assert( !ISNAN( rep[2]->gene(0).real ) );
   stateNow->T.z = rep[2]->gene(0).real;

   //  Set up the quaternion
   assert( !ISNAN( rep[3]->gene(0).real ) );
   stateNow->Q.x = rep[3]->gene(0).real;
   assert( !ISNAN( rep[3]->gene(1).real ) );
   stateNow->Q.y = rep[3]->gene(1).real;
   assert( !ISNAN( rep[3]->gene(2).real ) );
   stateNow->Q.z = rep[3]->gene(2).real;
   assert( !ISNAN( rep[3]->gene(3).real ) );
   stateNow->Q.w = rep[3]->gene(3).real;

   //  Copy the angles
   for (i=0; i<stateNow->ntor; i++) {
      assert( !ISNAN( rep[4]->gene(i).real ) );
      stateNow->tor[i] = rep[4]->gene(i).real;
   }
}

double Eval::operator()(const Representation *const *const rep)
{
   make_state_from_rep(rep, &stateNow, outlev, logFile);
   return eval();
}

double Eval::operator()(const Representation *const *const rep, const int term)
{
   make_state_from_rep(rep, &stateNow, outlev, logFile);
   return eval(term);
}


double Eval::eval()
{
#ifdef DEBUG
   (void) fprintf(logFile,"eval.cc eval() calling eval(3)\n");
#endif /* DEBUG */
   return eval(3); // default is total energy
}


double Eval::eval(const int term)

// Use this method, eval(int term), to compute just one particular term of the total energy
//
// we define term=0 as total energy
//           term=1 as total non-bonded energy, i.e. vdW+Hb+desolv
//           term=2 as total electrostatic energy
//           term=3 as total energy if invoked by eval()

{
   register int i;
   int   B_outside = 0;
   int   I_tor = 0;
   int   indx = 0;
   double energy = 0.0L;
   double retval = 0.0L;

        EnergyComponent	totalE;        // total energy components


#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/double Eval::eval(int term=%d)\n", term);
#endif /* DEBUG */

#ifdef DEBUG
    if (is_out_grid_info(stateNow.T.x, stateNow.T.y, stateNow.T.z)) {
       (void)fprintf(logFile,"eval.cc/stateNow.T is outside grid!\n");
    }
#endif /* DEBUG */

#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/eval(int term)  Converting state to coordinates...\n");
    printState( logFile, stateNow, 2 );
#endif /* DEBUG */
 
   // Ligand could be inside or could still be outside, check all the atoms...
   // cnv_state_to_coords(stateNow, vt, tlist, stateNow.ntor, crdreo, crd, natom);
   cnv_state_to_coords(stateNow, vt, tlist, stateNow.ntor, crdpdb, crd, natom,
    true_ligand_atoms, outlev, logFile);

#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/Checking to see if all coordinates are inside grid...\n");
#endif /* DEBUG */

   //  Check to see if crd is valid
   for (i=0; (i<natom)&&(!B_outside); i++) {
      B_outside = is_out_grid_info(crd[i][0], crd[i][1], crd[i][2]);
   }

   // Use standard energy function

#ifdef DEBUG
    if(B_outside) (void)fprintf(logFile,"eval.cc/Some coordinates are outside grid...\n");
    else (void)fprintf(logFile,"eval.cc/All coordinates are inside grid...\n");
#endif /* DEBUG */

    if (B_compute_intermol_energy) {
        if(term==3) // do not need energy breakdown in this eval() case
        energy = scale_eintermol * trilinterp( 0, natom, crd, charge, abs_charge, type, map, 
	     info, ignore_inter, NULL, NULL,
	     NULL_ENERGY_BREAKDOWN);
        else
        energy = scale_eintermol * trilinterp( 0, natom, crd, charge, abs_charge, type, map, 
	     info, ignore_inter, NULL, &totalE,
	     NULL_ENERGY_BREAKDOWN);
    }
    
#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/double Eval::eval(int term=%d) after trilinterp, energy= %.5lf\n", term, energy);
#endif /* DEBUG */
    energy += eintcal( nonbondlist, ptr_ad_energy_tables, crd, Nnb,
		 Nnb_array, NULL, /* group_energy,  perhaps do not need energy breakdown MP TODO 2012 */
                 B_calcIntElec, B_include_1_4_interactions, 
                 scale_1_4, qsp_abs_charge,
                 B_use_non_bond_cutoff, B_have_flexible_residues,
		 outlev, logFile);
#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/double Eval::eval(int term=%d) after eintcal, energy= %.5lf\n", term, energy);
#endif /* DEBUG */
 
    if (B_isGaussTorCon) {
        for (I_tor = 0; I_tor <= stateNow.ntor; I_tor++) {
            if (B_isTorConstrained[I_tor] == 1) {
                indx = RadiansToDivs( WrpModRad(stateNow.tor[I_tor]) );
                if (B_ShowTorE) {
                    energy += (double)(US_TorE[I_tor] = US_torProfile[I_tor][indx]);
                } else {
                    energy += (double)US_torProfile[I_tor][indx];
                }
            }
        } // I_tor
    }/*if*/

   // increment evaluation counter only for "total energy" case
   if(term==3) num_evals++;

   if ((!finite(energy)) || ISNAN(energy)) {
      (void)fprintf( logFile, "eval.cc:  ERROR!  energy is %s!\n\n",
       (!finite(energy))?"infinite":"not a number");
#define DUMMYATOMSTUFFINF "ATOM  #####  C   INF X   1    "
#define DUMMYATOMSTUFFNAN "ATOM  #####  C   NAN X   1    "
      for (i=0; i<natom; i++) {
            print_PDBQT_atom_resstr( logFile, "", i,   DUMMYATOMSTUFFINF, crd, 
             0.0, 0.0, charge[i],"", "\n"); 
      } // i
   }
    switch (term) {
    default:
    case 0:
    case 3:
        // Return the total energy, scaled by e_intermol
        retval = energy;
        break;
    case 1:
        // Return the non-bonded energy, vdW+Hb+desolv, not scaled by e_intermol
        retval = (double)totalE.vdW_Hb+totalE.desolv;
        break;
    case 2:
        // Return the electrostatics energy, not scaled by e_intermol
        retval = (double)totalE.elec;
        break;
    }

#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/double Eval::eval(int term=%d) returns retval= %.5lf\n", term, retval);
#endif /*DEBUG*/
   return(retval);
}

int Eval::write(const Representation *const *const rep,
 const int true_ligand_atoms, const int outlev, FILE *logFile)
{
    int i=0, retval=0;
    //char rec14[14];

#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/int Eval::write(FILE *out_file, Representation **rep)\n");
#endif /*DEBUG*/

    make_state_from_rep(rep, &stateNow, outlev, logFile);
    cnv_state_to_coords(stateNow, vt, tlist, stateNow.ntor, crdpdb, crd, natom,
     true_ligand_atoms, outlev, logFile);
    for (i=0; i<natom; i++) {
            print_PDBQT_atom_resstr( logFile, "", i,   " C   RES     1 ", crd, 
             0.0, 0.0, charge[i],"", "\n"); 
    } // i
    return retval;
}

#ifdef USING_COLINY // {
double Eval::operator()(const double* const vec, const int len, const int outlev, FILE *logFile)
{
   make_state_from_rep(vec, len, &stateNow, outlev, logFile);
   return eval();
}


void make_state_from_rep(const double *const rep, const int n, /* not const */ State *const now, const int outlev, FILE *logFile)
{
#   ifdef DEBUG
    (void)fprintf(logFile, "eval.cc/make_state_from_rep(double *rep, int n, State *now)\n");
#   endif /* DEBUG */

    //  Do the translations
    now->T.x = rep[0];
    now->T.y = rep[1];
    now->T.z = rep[2];

    //  Set up the quaternion
    now->Q.x = rep[3];
    now->Q.y = rep[4];
    now->Q.z = rep[5];
    now->Q.w = rep[6];

    //  Copy the angles
    now->ntor = n - 7;
    for (int i=0, j=7; j<n; i++, j++) {
      now->tor[i] = rep[j];
    }

    //mkUnitQuat(&(now->Q));
}

/* next function is for Coliny only */
extern Eval evaluate;

double ADEvalFn(/* not const */ double *const x, const int n)
{
//
// Normalize the data
//
//
// Quaternion vector
/*
double sum=0.0;
if (x[3] < 0.0) x[3] = 1e-16;
if (x[4] < 0.0) x[4] = 1e-16;
if (x[5] < 0.0) x[5] = 1e-16;
*/
double sum = sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5]);
if (sum < 1e-8)
   x[3]=x[4]=x[5]=1.0L/sqrt(3.0L);
   else {
      x[3] /= sum;
      x[4] /= sum;
      x[5] /= sum;
      }

// torsion angles
for (int i=6; i<n; i++)
  x[i] = WrpModRad(x[i]);

return ::evaluate(x,n);
}
//
#endif // USING_COLINY // }
