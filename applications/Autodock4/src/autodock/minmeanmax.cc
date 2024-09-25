/*

 $Id: minmeanmax.cc,v 1.9 2010/08/27 00:05:07 mp Exp $

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

/* minmeanmax.cc */

/* Calculates the min, mean and max of each gene for a given population */

/* N.B.: the macro is defined in "hybrids.h" */

#include <stdio.h>
#include "rep.h"
#include "support.h" /* Population defined here */

#define NUM_BINS 50

#define which_bin(x, min, max, nbins)  int((nbins)*((x)-(min))/((max)-(min)))
#define clamp_range(x, min, max) ((x) < (min)) ? (min) : ((x) > (max)) ? (max) : (x)

void minmeanmax( FILE *const fp, const Population &pop, const int num_generations, const GridMapSetInfo *const info )
{
#ifdef DEBUG
   fprintf( fp, "minmeanmax.cc/minmeanmax(); initialization\n" );
#endif
   register int i=0, g=0, b=0;
   int     num_indvs=1, num_genes=1;
   int last_bin = NUM_BINS - 1;
   // double  value=0., *minimum, *sum, *maximum, *best;
   double  value=0.;
   double minimum[ 7+MAX_TORS ], sum[ 7+MAX_TORS ], maximum[ 7+MAX_TORS ], best[ 7+MAX_TORS ];
   // double *bin_min, *bin_max;
   double bin_min[ 7+MAX_TORS ], bin_max[ 7+MAX_TORS ];
   double  best_energy=1e20;
   // double energy;
   double bmin=0., bmax=0.;
   Element temp;
   int best_individual = 0;
   double best_gene = 0.;
   char binchar = ' ';

#ifdef DEBUG
   fprintf( fp, "minmeanmax.cc/minmeanmax(); about to set num_indvs and num_genes\n" );
#endif
   num_indvs = pop.num_individuals();
   num_genes = pop[0].genotyp.num_genes();
#ifdef DEBUG
   fprintf( fp, "minmeanmax.cc/minmeanmax(); num_indvs= %d  num_genes=%d\n", num_indvs, num_genes);
#endif

   int num[ 7+MAX_TORS ][ NUM_BINS ];
#ifdef DEBUG
   fprintf( fp, "minmeanmax.cc/minmeanmax(); num[ %d ][ %d ] created.\n", 7+MAX_TORS, NUM_BINS);
#endif

   // Initialise arrays
   for (g=0; g < 7+MAX_TORS; g++) {
       minimum[ g ] = 0.;
       sum[ g ] = 0.;
       maximum[ g ] = 0.;
       best[ g ] = 0.;
   }

   // Set bin maxima and minima
   // Translation x,y,z components correspond to 0, 1, 2
   for (g=0; g < 3; g++) {
       bin_min[g] = info->lo[g];
       bin_max[g] = info->hi[g];
   }
   // Quaternion qx,qy,qz,qw components correspond to 3, 4, 5, 6
   for (g=3; g < 7; g++) {
       bin_min[g] = -1.;
       bin_max[g] = 1.;
   }
   // Torsion angles corresponds to 7 to (7+MAX_TORS - 1)
   for (g=7; g < 7+MAX_TORS; g++) {
       bin_min[g] = -PI;
       bin_max[g] = PI;
   }

#ifdef DEBUG
   fprintf( fp, "minmeanmax.cc/minmeanmax(); bin_min[] and bin_max[] have been set.\n");
#endif

   // energy  = best_energy;

   // Initialise the number in each bin.
   for (g=0; g < num_genes; g++) {
      for (i=0; i < NUM_BINS; i++) {
          num[g][i] = 0;
      }
   }
#ifdef DEBUG
   fprintf( fp, "minmeanmax.cc/minmeanmax(); num[][] set to zero.\n");
#endif

   // Set the initial values for the
   // minimum, sum, maximum and best arrays.
   for (g=0; g < num_genes; g++) {
#ifdef DEBUG
      fprintf( fp, "minmeanmax.cc/minmeanmax(); g= %d\n", g);
#endif
      temp       = pop[0].genotyp.gread(g);
      //value = temp.real;
#ifdef DEBUG
      fprintf( fp, "minmeanmax.cc/minmeanmax(); about to temp.real \n");
#endif
#ifdef DEBUG
      fprintf( fp, "minmeanmax.cc/minmeanmax(); temp.real= %lf\n", temp.real);
#endif
      minimum[g] = temp.real;
      sum[g]     = temp.real;
      maximum[g] = temp.real;
      best[g]    = temp.real;
      // build up the histogram, accumulating each gene's value in the appropriate bin:
      b = which_bin( temp.real, bin_min[g], bin_max[g], NUM_BINS );
#ifdef DEBUG
   fprintf( fp, "minmeanmax.cc/minmeanmax(); after which_bin, b= %d\n", b);
#endif
      b = clamp_range( b, 0, last_bin );
#ifdef DEBUG
   fprintf( fp, "minmeanmax.cc/minmeanmax(); after clamp_range, b= %d\n", b);
#endif
      ++num[ g ][ b ];
      if (pop[0].value(Normal_Eval) < best_energy) {
          best_energy = pop[0].value(Normal_Eval);
          best_individual = 0;
      }
   }

   // For all the individuals, i, in the population, 
   // except the 0-th individual:
   for (i=1; i < num_indvs; i++) {
      // For all the genes, g, in the individual, i 
      for (g=0; g < num_genes; g++) {
#ifdef DEBUG
         fprintf( fp, "minmeanmax.cc/minmeanmax(); g= %d\n", g);
#endif
         temp = pop[i].genotyp.gread(g);
         value = temp.real;
#ifdef DEBUG
      fprintf( fp, "minmeanmax.cc/minmeanmax(); value= %lf\n", value);
#endif
         if (minimum[g] > value) {
            minimum[g] = value;
         }
         if (maximum[g] < value) {
            maximum[g] = value;
      }
         sum[g] += value;
         // build up the histogram, accumulating each gene's value in the appropriate bin:
         b = which_bin( value, bin_min[g], bin_max[g], NUM_BINS );
#ifdef DEBUG
   fprintf( fp, "minmeanmax.cc/minmeanmax(); after which_bin, b= %d\n", b);
#endif
         b = clamp_range( b, 0, last_bin );
#ifdef DEBUG
   fprintf( fp, "minmeanmax.cc/minmeanmax(); after clamp_range, b= %d\n", b);
#endif
         ++num[ g ][ b ];
      } // g
      if (pop[i].value(Normal_Eval) < best_energy) {
          best_energy = pop[i].value(Normal_Eval);
          best_individual = i;
      }
   } // i

   fprintf( fp, "\nbest_energy= %5.2lf  best_individual= %d\n\n", best_energy, best_individual);

   // print out the statistics for each gene across the population
   for (g=0; g < num_genes; g++) {
#ifdef DEBUG
       fprintf( fp, "minmeanmax.cc/minmeanmax(); g= %d\n", g);
#endif
       // best individual's gene
       best_gene = pop[ best_individual ].genotyp.gread( g ).real;
       fprintf( fp, "\nnum_generations= %d  gene= %2d  best_gene= %5.2lf\n\n", num_generations, g, best_gene);
       // g-th gene in all individuals
       fprintf( fp, "num_generations= %d  gene= %2d    minimum= %8.3f  mean= %8.3f  maximum= %8.3f\n", num_generations, g, minimum[g], sum[g]/num_indvs, maximum[g] ); fflush( fp );
       for ( b=0; b < NUM_BINS; b++ ) {
          bmin = bin_min[g] + b * (bin_max[g] - bin_min[g]) / NUM_BINS;
          bmax = bin_min[g] + (b+1) * (bin_max[g] - bin_min[g]) / NUM_BINS;
          if ((bmin < best_gene) && (best_gene < bmax)) {
              binchar = '*';
          } else {
              binchar = ' ';
          }
          fprintf( fp, "gene= %2d    bin= %2d  binmin= %5.1lf  binmax= %5.1lf  num= %3d  %c|", g, b, bmin, bmax, num[ g ][ b ], binchar );
          for ( i=0; i<num[ g ][ b ]; i++ ) {
             fprintf( fp, "#" );
   }
          fprintf( fp, "\n" );
   }
   } // g

   return;
}
