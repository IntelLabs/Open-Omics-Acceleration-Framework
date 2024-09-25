/*

 $Id: prClusterHist.cc,v 1.13 2014/06/12 01:44:07 mp Exp $

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
#include "prClusterHist.h"
#include "constants.h"


void prClusterHist( const int ncluster,
                    const int irunmax,
                    ConstReal clus_rms_tol,
                    const int num_in_clu[MAX_RUNS],
                    const int cluster[MAX_RUNS][MAX_RUNS],
                    const Real econf[MAX_RUNS],
                    const Real clu_rms[MAX_RUNS][MAX_RUNS],
                    const Real ref_rms[MAX_RUNS],
		    const int outlev,
		    FILE *logFile)
{
    Real          etot = 0.,
                  eavg = 0.;

    register int  Rank=0,
                  j=0;

    int           num_multi_mem_clu=0,
                  ConfNum=0,
                  Rank1=1;

    Real          p_i = 0.0;
    Real          N = 1.0;
    Real          Information_entropy = 0.0;

    N = (Real) irunmax;

    (void)fprintf( logFile, "\n\n" );
    (void)fprintf( logFile, "________________________________________________________________________________\n\n" );

    (void)fprintf( logFile, "Number of distinct conformational clusters found = %d,  out of %d runs,\nUsing an rmsd-tolerance of %.1f A\n\n\n", ncluster, irunmax, (double)clus_rms_tol );
 
    (void)fprintf( logFile, "\tCLUSTERING HISTOGRAM\n" );
    (void)fprintf( logFile, "\t____________________\n\n");

    (void)fprintf( logFile, "________________________________________________________________________________\n");
    (void)fprintf( logFile, "     |           |     |           |     |                                    \n" );
    (void)fprintf( logFile, "Clus | Lowest    | Run | Mean      | Num | Histogram                          \n" );
    (void)fprintf( logFile, "-ter | Binding   |     | Binding   | in  |                                    \n" );
    (void)fprintf( logFile, "Rank | Energy    |     | Energy    | Clus|    5    10   15   20   25   30   35\n" );
    (void)fprintf( logFile, "_____|___________|_____|___________|_____|____:____|____:____|____:____|____:___");

    double U_internal_energy_summation = 0.0;
    double Q_partition_function = 0.0;
    /* */
    double U_internal_energy = 0.0;
    double A_free_energy = 0.0;
    double S_entropy = 0.0;
    /* */
    double this_energy = 0.0;
    double RT = Rcal * TK;

    // Loop over all the clusters
    for (Rank = 0;  Rank < ncluster;  Rank++) {
        Rank1 = Rank + 1;

        // Print the Cluster Rank
        (void)fprintf( logFile, "\n%4d |", Rank1);

        // Print the Lowest Binding Energy in this cluster
        if (econf[cluster[Rank][0]] > 999999.99) {
            (void)fprintf( logFile, "%+10.2e", econf[cluster[Rank][0]]);
        } else {
            (void)fprintf( logFile, "%+10.2f", econf[cluster[Rank][0]]);
        }

        // Print the number of the Run that found this docking
        (void)fprintf( logFile, " |%4d |", cluster[Rank][0]+1);

        // Calculate the proportion in this cluster, p_i
        p_i = (Real) num_in_clu[Rank] / N;

        // Add this cluster's contribution to the information entropy
        // using log to the base N, where N is the number of dockings, irunmax
        Information_entropy -= p_i * log( p_i ) / log( N );

        // Print the Mean Binding Energy
        if (num_in_clu[Rank] > 1) {
            // There is are more than one dockings in this cluster.

            // Calculate the average energy of this cluster
            etot = 0.;
            for (j = 0;  j < num_in_clu[Rank]; j++ ) {
                this_energy = econf[ cluster[Rank][j] ];

                U_internal_energy_summation += this_energy * exp( -this_energy / RT );
                Q_partition_function += exp( -this_energy / RT );

                etot += this_energy;
            }
            eavg = etot / (Real)num_in_clu[Rank];
            num_multi_mem_clu++;

            // Print the average energy of this cluster
            if (eavg > 999999.99) {
                (void)fprintf( logFile, "%+10.2e |", eavg );
            } else {
                (void)fprintf( logFile, "%+10.2f |", eavg );
            }
        } else {
            // There is only one docking in this cluster
            this_energy = econf[ cluster[Rank][0] ];

            U_internal_energy_summation += this_energy * exp( -this_energy / RT );
            Q_partition_function += exp( -this_energy / RT );

            // Print the average energy of this cluster
            if (this_energy > 999999.99) {
                (void)fprintf( logFile, "%+10.2e |", this_energy);
            } else {
                (void)fprintf( logFile, "%+10.2f |", this_energy);
            }
        }

        // Print the Number in Cluster
        (void)fprintf( logFile, "%4d |", num_in_clu[Rank] );

        // Print the Histogram
        // Print a '#' symbol for each member in this cluster
        for (j=0;  j<num_in_clu[Rank]; j++) {
            (void)fprintf( logFile, "#" );
        }/*j*/
    }/*Rank*/

    

    (void)fprintf( logFile, "\n" );
    (void)fprintf( logFile, "_____|___________|_____|___________|_____|______________________________________\n");
    (void)fprintf( logFile, "\n" );

    // Print the number of clusters which have more than one docking result
    if (num_multi_mem_clu > 0) {
        (void)fprintf( logFile, "\nNumber of multi-member conformational clusters found = %d, out of %d runs.\n\n", num_multi_mem_clu, irunmax );
    }


    // Print the table of RMSD values for each docking result
    (void)fprintf( logFile, "\tRMSD TABLE\n" );
    (void)fprintf( logFile, "\t__________\n\n");

    (void)fprintf( logFile, "_____________________________________________________________________\n");
    (void)fprintf( logFile, "     |      |      |           |         |                 |\n");
    (void)fprintf( logFile, "Rank | Sub- | Run  | Binding   | Cluster | Reference       | Grep\n");
    (void)fprintf( logFile, "     | Rank |      | Energy    | RMSD    | RMSD            | Pattern\n");
    (void)fprintf( logFile, "_____|______|______|___________|_________|_________________|___________\n" );

    for (Rank = 0;  Rank < ncluster;  Rank++) {
        Rank1 = Rank + 1;
        ConfNum = 0;
        for (j=0;  j<num_in_clu[Rank]; j++) {
            ++ConfNum;
            if (econf[cluster[Rank][j]] > 999999.99) {
                (void)fprintf( logFile, "%4d   %4d   %4d  %+10.2e  %8.2f  %8.2f           RANKING\n", Rank1, ConfNum,  cluster[Rank][j]+1, econf[cluster[Rank][j]], ((j==0)?(0.):(clu_rms[Rank][j])), ref_rms[cluster[Rank][j]] );
            } else {
                (void)fprintf( logFile, "%4d   %4d   %4d  %+10.2f  %8.2f  %8.2f           RANKING\n", Rank1, ConfNum,  cluster[Rank][j]+1, econf[cluster[Rank][j]], ((j==0)?(0.):(clu_rms[Rank][j])), ref_rms[cluster[Rank][j]] );
            }
        }   /*j*/
        //(void)fprintf( logFile, ".......................................................................\n" );
    }/*Rank*/
    (void)fprintf( logFile, "_______________________________________________________________________\n\n");


    // Print the Information Entropy value at this rmstol clustering tolerance
    (void)fprintf( logFile, "\n\n\tINFORMATION ENTROPY ANALYSIS FOR THIS CLUSTERING\n" );
    (void)fprintf( logFile,     "\t________________________________________________\n" );
    (void)fprintf( logFile, "\n\n" );
    (void)fprintf( logFile, "Information entropy for this clustering = %4.2f  (rmstol = %.2f Angstrom)\n", Information_entropy, clus_rms_tol );
    (void)fprintf( logFile, "\n" );
    (void)fprintf( logFile, "_______________________________________________________________________\n\n");


    // Finish the calculation of the internal energy, U, which depends on the partition function, Q:
    U_internal_energy = (1. / Q_partition_function) * U_internal_energy_summation;
    // Calculate the free energy, A
    A_free_energy = -RT * log( Q_partition_function );
    // Calculate the entropy, S
    S_entropy = -( A_free_energy - U_internal_energy ) / TK;

    // Print the Statistical Thermodyamics values 
    (void)fprintf( logFile, "\n\tSTATISTICAL MECHANICAL ANALYSIS\n" );
    (void)fprintf( logFile,   "\t_______________________________\n" );
    (void)fprintf( logFile, "\n\n" );
    (void)fprintf( logFile, "Partition function, Q = %8.2f            at Temperature, T = %.2f K\n", Q_partition_function, TK );
    (void)fprintf( logFile, "Free energy,        A ~ %8.2f kcal/mol   at Temperature, T = %.2f K\n", A_free_energy, TK );
    (void)fprintf( logFile, "Internal energy,    U = %8.2f kcal/mol   at Temperature, T = %.2f K\n", U_internal_energy, TK );
    (void)fprintf( logFile, "Entropy,            S = %8.2f kcal/mol/K at Temperature, T = %.2f K\n", S_entropy, TK );
    (void)fprintf( logFile, "\n" );
    (void)fprintf( logFile, "_______________________________________________________________________\n\n");


    fflush( logFile );
}
/* EOF */
