/*

 $Id: support.cc,v 1.47 2014/07/18 05:43:18 mp Exp $

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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "eval.h"
#include "support.h"
#include "stateLibrary.h"
#include "structs.h"

//#define DEBUG 
extern FILE *logFile;  // DEBUG only
extern int outlev;

//extern class Eval evaluate;


//  These are the member functions for the support classes.


Population::Population(const Population &original)
: lhb(original.lhb), size(original.size), evaluate(original.evaluate)
{
   register int i;

//MP #ifdef DEBUG
   (void)fprintf(logFile, 
 "support.cc %d /Population::Population(Population &original) size=%d\n",
  __LINE__, size);
//MP #endif /* DEBUG */

   heap = new Individual[size];
   for (i=0; i<size; i++) {
      heap[i] = original.heap[i];
      heap[i].age = 0L; // gmm, 1998-07-14
   }
   // MPique 2014 TODO should this also copy the end_of_branch ?
}

/*  Heap Functions:  In this case, the heap condition means the maximal 
    element wrt fitness (i.e. the best) is at the top of the heap.  lhb 
    is the index of the last element to be inserted into the heap.  Some 
    the standard functions on the heap can be accomplished in the following 
    manner:
       Root = size - 1  (Note: that the root is fixed)
       LeftChild(I) = 2I - size
       RightChild(I) = 2I - size - 1
       Parent(I) = (I + size + 1)/2
    It is important to notice that the heap condition is maintained from
    lhb to size-1 *in reverse order*.
*/

void Population::swap(Individual &individual1, Individual &individual2) const
{
   Individual temp;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Population::swap(Individual &individual1, Individual &individual2)\n");
#endif /* DEBUG */


   temp = individual1;
   individual1 = individual2;
   individual2 = temp;
}

/*  This routine assumes that the heap condition is satisfied
    between lhb and size-1 and that the new individual is in
    position lhb-1.
*/
void Population::SiftUp(void) /* not const */ 
{
   int i, parent;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Population::SiftUp(void)\n");
#endif /* DEBUG */


   i = lhb-1;
   while (((parent=(i+size+1)/2)<size)&&(heap[parent].value(Normal_Eval)>heap[i].value(Normal_Eval))) {
      swap(heap[parent], heap[i]);
      i = parent;
   }
   lhb--;
}

/*  This routine assumes that the heap condition is satisfied
    between lhb & size-2 initially, and that the individual at size-1
    needs to be accomodated.*/
void Population::SiftDown(void) /* not const */
{
   int i, child;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Population::SiftDown(void)\n");
#endif /* DEBUG */


   i = size-1;
   while ((child=2*i-size)>=lhb) {
      if (child-1>=lhb) {
         if (heap[child-1].value(Normal_Eval)<heap[child].value(Normal_Eval)) {
            child--;
         }
      }
      /*  Now child holds the index of the best child of i  */
      if (heap[i].value(Normal_Eval)<heap[child].value(Normal_Eval)) {
         break;
      } else {
         swap(heap[child], heap[i]);
         i = child;
      }
   }
}

void Population::msort(const int m) /* not const */
{
   register int i;
   char error_message[200];

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Population::msort(int m=%d)\n", m);
#endif /* DEBUG */

   

   //  First make a heap of the whole array, i.e lhb = 0 & uhb = size
   lhb = size-1;
   while (lhb>0) {
      SiftUp();
   }

   assert(lhb==0);
 
   //  Now place the m best members at the beginning of the array
   if (m==size) return; //we're done
   if (m>size) {
        sprintf(error_message, "support.cc/Population::msort(int m=%d) -- ERROR!  m > size!\n\n", m);
        stop(error_message);
   }

   for (i=0; i<m && i<size-1; i++) {
#ifdef DEBUG
   (void)fprintf(stderr, "support.cc/placing %d of %d best of %d, lhb=%d )\n", i, m, size, lhb);
#endif /* DEBUG */
      swap(heap[i], heap[size-1]);
      lhb++;
#ifdef DEBUG
   (void)fprintf(stderr, "support.cc/calling SiftDown lhb=%d )\n", lhb);
#endif /* DEBUG */
      SiftDown();
   }
   
   //  Assert: heap[0..m-1] sorted
}

void Population::print(FILE *const output, const int num) const {
   register int i;

#ifdef DEBUG
   (void)fprintf(output, "support.cc/void Population::print(FILE *const output, int num=%d)\n", num);
#endif /* DEBUG */

   (void)fprintf( output, "The top %d individuals in the population:\n\n", num);
   for (i=0; i<num; i++) {
      (void)fprintf( output, "(%d):\t %8.2f\n", i+1, heap[i].value(Always_Eval));
   }
}

static int
doublecompare(const void *const p1, const void *const p2)
{
	const double i = *((double *)p1);
	const double j = *((double *)p2);
	if ( i > j ) return 1;
	if ( i < j ) return -1;
	return 0;
}

static double
compute_k_quantile(const int k, const int q, const double energy[], const int size)
{
 // return k-th (0-origin) "q"-tile of array (assumed sorted ascending)
 // example: median is k=1  q=2
 // example: first quartile is k=1 q=4
 //
 // If the computed index falls between two input sample values, the
 //   weighted average of the two is returned.
 //
 // M Pique December 2009 from wikipedia "Quantile" 
 //    "Estimating the quantiles of a population", "Weighted average"
 //    http://en.wikipedia.org/wiki/Quantile 
  const double p = (size-1) * k / q;
  double j_double; // integer part of p
  const double g = modf(p, &j_double); // fractional part of p
  const int ilow = (int) j_double;
  assert(ilow>=0 && ilow<size);
  if (g == 0.0) return energy[ilow];
  const int ihigh=ilow+1; // indices to be weighted
  assert(ihigh<size);
  return energy[ilow] + g*(energy[ihigh]-energy[ilow]);
 }


  
   

int Population::printPopulationStatistics(FILE *const output, const int level, const char suffix[]) /* not const changed in support.cc hack TODO */  {
// write best energy, etc, depending on level, followed by suffix string
// return 0 if OK, non-zero if error
// Code adapted from gs.cc  - M Pique  December 2009
int returnCode=0;
double sum, worst_e;
unsigned long oldest;
int oldestIndividual;
// best_e and best_i are set here, TODO move to better place M Pique 2010-03
   if(size<=0) return 2; // error, give up
   sum=best_e=worst_e=heap[best_i=0].value(Normal_Eval);
   oldest=heap[oldestIndividual=0].age;
   for (int i=1; i<size; i++) {
      double e = heap[i].value(Normal_Eval);
      sum+=e;
      if(e>worst_e) worst_e = e;
      if(e<best_e) {
          best_e = e;
	  best_i = i;
	  }
      if (heap[i].age >= oldest) {
          oldest = heap[i].age;
          oldestIndividual = i;
      }
   }

   // level == 1 - short output as gs.cc did previously: print oldest and best (in exactly same format)
   // level == 2 - add age info (supporting DEBUG3 option in gs.cc)
   // level == 3 - no age info, but add worst, mean, median, quartiles, 
   //    standard deviation to (1) in simpler format
    switch (level) {
    case 1:
	    // MP debug @@ next line added [best_i]
	    (void)fprintf(output, " Oldest's energy: %.3f    Lowest energy: %.3f [%d]", 
               heap[oldestIndividual].value(Normal_Eval), best_e, best_i);
	       break;
    case 2:
       (void)fprintf(output, " Oldest ind.: %u/%u, age: %lu, energy: %.3f    Lowest energy individual: %u/%u, age: %lu, energy: %.3f", 
               oldestIndividual+1, size, 
	       heap[oldestIndividual].age, 
               heap[oldestIndividual].value(Normal_Eval), 
	       best_i+1, size, 
               heap[best_i].age, best_e);
	       break;
    default:
	{
	double sum_squares, median, stddev;
	//double q14, q34; // 1st and 3rd quartiles
	//double q15, q45; // 1st and 4th quintiles
	double energy[size]; // array for sorting
	const double mean = sum/size;
	sum_squares=0;
	for (int i=0; i<size; i++) energy[i] = heap[i].value(Normal_Eval);
	for (int i=0; i<size; i++) sum_squares += (energy[i]-mean)*(energy[i]-mean);
	if(size<2) stddev=0; else stddev=sqrt(sum_squares/(size-1));
	// compute median and quartiles 
	//  Avoid using msort because reordering the actual population
	//  causes runs to produce different results depending on statistics output level
	qsort( (void *) energy, size, sizeof(*energy), doublecompare); // ascend

	if( 1 == (size%2) ) median = energy[size/2]; // odd size
	else median = (energy[size/2] + energy[size/2-1])/2.0; // even size

	(void) fprintf(output, "Lowest: %.3f Highest: %.3f Mean: %.3f Median: %.3f Std.Dev: %.3f", 
	  best_e, worst_e, mean, median, stddev);

	// quartiles and quintiles
	//q15 = compute_k_quantile(1, 5, energy, size);
	//q14 = compute_k_quantile(1, 4, energy, size);
	//q34 = compute_k_quantile(3, 4, energy, size);
	//q45 = compute_k_quantile(4, 5, energy, size);
#define quantile(k, q) \
	(void) fprintf(output, " Q%d/%d: %.3f", \
	  k, q, compute_k_quantile(k, q, energy, size))

	quantile(1,10);
	quantile(1, 5);
	quantile(1, 4);
	quantile(3, 4);
	quantile(4, 5);
	quantile(9,10);

#undef quantile
	  // debug print every energy:
	if(level>3) for(int i=0; i<size; i++) fprintf(output, " %.3f", energy[i]);
	}
	break;
    }
    if(suffix!=NULL) fprintf(output, "%s", suffix);
    return returnCode;
}

int Population::printPopulationStatisticsVerbose(FILE *const  output, 
 const unsigned int generations, const long int nevals, const int ntor, const char suffix[])  /* not const changed in support.cc hack TODO */ { /* print with generations & #evals */
int returnCode=0;
   // print "Population at Generation:" line with low/high/mean/median/stddev...
   (void) fprintf(output, "Population at Generation: %3u ", generations);

   // highest level, without newline at end, sets population best_i, best_e:
   returnCode= printPopulationStatistics(output, 3, "");

   // print "State0:" field with compact state for best individual.
   (void) fprintf(output," State0: ");
   heap[best_i].printIndividualsState( output, ntor, 4); // 4 means print compact state, no following newline

   (void) fprintf(output, " Num.evals: %ld%s", nevals, suffix );
    return returnCode; 
}
 	

void Population::printPopulationAsStates(FILE *const output, const int num, const int ntor) const {
   register int i;
#ifdef DEBUG2
   register int j;
   char resstr[LINE_LEN];
#endif /* DEBUG2 */
   double thisValue;

#ifdef DEBUG
   (void)fprintf(output, "support.cc/void Population::printPopulationAsStates(FILE *const output, int num=%d, int ntor=%d)\n", num, ntor);
#endif /* DEBUG */

   // Print an XML-like tag indicating this is a population, with attribute size
   // being the number of individuals in the population
   (void)fprintf( output, "<population size=\"%d\">\n", num);
   for (i=0; i<num; i++) {
      thisValue = heap[i].value(Always_Eval);
      (void)fprintf( output, "%4d\t%9.4lg\t", i+1, thisValue);
      (void)fprintf( output, "%4lu\t", heap[i].age );
      heap[i].printIndividualsState(output, ntor, 0);
      (void)fprintf( output, "\n");

#ifdef DEBUG2
      // to print only infinite or NaN structures // if (!finite(thisValue) || ISNAN(thisValue)) {//debug
      // Convert state to coords and print it out...//debug
      cnv_state_to_coords(heap[i].state(ntor), heap[i].mol->vt,  heap[i].mol->tlist,  ntor, heap[i].mol->crdpdb,  heap[i].mol->crd,  heap[i].mol->natom,
       true_ligand_atoms, outlev, output);//debug
      (void)fprintf(output, "MODEL     %4d\n", i+1);
      for (j=0; j<heap[i].mol->natom; j++) {//debug
            print_PDBQT_atom_resnum( output, "", j, DUMMYATOMSTUFFINF,   1,  // replace 1 with i+1 for incrementing residue numbers.
            heap[i].mol->crd, 0.0, 0.0, 0.0, "", "\n"); //debug
      }/*j*///debug
      (void)fprintf(output, "ENDMDL\n");
      // to print only infinite or NaN structures // }// thisValue is either infinite or not-a-number.//debug
#endif /* DEBUG2 */

   }// i
   (void)fprintf( output, "</population>\n");
}

void Population::printPopulationAsCoordsEnergies(FILE *const output, const int num, const int ntor) const {
   register int i;
   double thisValue;

#ifdef DEBUG
   (void)fprintf(output, "support.cc/void Population::printPopulationAsCoordsEnergies(FILE *const output, int num=%d, int ntor=%d)\n", num, ntor);
#endif // DEBUG

   if( output == NULL) stop("printPopulationAsCoordsEnergies received NULL file pointer");
   (void)fprintf( output, "The top %d individuals in the population:\n\n", num);
   for (i=0; i<num; i++) {

      // Print the number of this individual in the population (counting from 1, not 0)
      (void)fprintf( output, "%d\t", i+1);

      // Print the translation
      heap[i].printIndividualsState(output, ntor, 3);  // 3 means print just the translation
      //heap[i].printIndividualsState(output, ntor, 0);  // 0 means print the whole state
      (void)fprintf( output, "\n");

      // Print the energy
      thisValue = heap[i].value(Normal_Eval); // was Always_Eval before 13-Jan-2006
      (void)fprintf( output, "\t%9.2lf", thisValue);

      // Print the non-bonded energy, i.e. vdW+Hb+desolv
      thisValue = heap[i].value(Always_Eval_Nonbond);
      (void)fprintf( output, "\t%9.2lf", thisValue);
   
      // Print the electrostatic energy, i.e. elec
      thisValue = heap[i].value(Always_Eval_Elec);
      (void)fprintf( output, "\t%9.2lf", thisValue);
     
      // Write a newline at the end
      (void)fprintf( output, "\n");

      // We need the coordinates of this individual to compute the electrostatic and nonbond energies
      //cnv_state_to_coords( heap[i].state(ntor), heap[i].mol->vt,  heap[i].mol->tlist,  ntor, heap[i].mol->crdpdb,  heap[i].mol->crd,  heap[i].mol->natom,
       //true_ligand_atoms, outlev, output);

   }// i
   (void)fprintf( output, "\n");
}

void Population::set_eob(/* not const */ int init_end_of_branch[MAX_TORS])
// Set the end_of_branch[MAX_TORS] array
{
   for (register int i=0; i<MAX_TORS; i++) {
       end_of_branch[i] = init_end_of_branch[i];
   }
}

int Population::get_eob(const int init_tor) const
// Get the end_of_branch[] value for the supplied torsion number, init_tor
{
    if ((init_tor >= 0) && (init_tor < MAX_TORS)) {
        return end_of_branch[init_tor];
    } else {
	char error_message[200];
        (void)sprintf(error_message, "support.cc/Population::get_eob(int init_tor=%d) -- ERROR!  Attempt to access out-of-bounds torsion!\n", init_tor);
        stop(error_message);
	return(0); // dummy to keep lint happy NOTREACHED
    }
}

Genotype::Genotype(unsigned int init_number_of_vectors, Representation **
const init_rep_vector)
: number_of_vectors(init_number_of_vectors), rep_vector(init_rep_vector), 
  modified(0)
{
   register unsigned int i, j, k;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Genotype::Genotype(unsigned int init_number_of_vectors=%d, Representation **init_rep_vector)\n", init_number_of_vectors);
#endif /* DEBUG */


   number_of_genes = 0;
   for (i=0; i<number_of_vectors; i++) {
      number_of_genes += rep_vector[i]->number_of_points();
#ifdef DEBUG
      (void)fprintf(logFile, "support.cc/Genotype::Genotype(init_number_of_vectors=%d, **init_rep_vector) number_of_genes=%d   rep_vector[%d]->number_of_points()=%d\n", init_number_of_vectors, number_of_genes, i, rep_vector[i]->number_of_points());
#endif /* DEBUG */
   }
 
   i=0;
   lookup = new Lookup[number_of_genes];
   for (j=0; j<number_of_vectors; j++) {
      for (k=0; k<rep_vector[j]->number_of_points(); k++) {
         lookup[i].vector = j;
         lookup[i].index = k;
         i++;
      }
   }
}

Genotype::Genotype(const Genotype &original)
{
   register unsigned int i;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Genotype::Genotype(Genotype const &original)\n");
#endif /* DEBUG */

   number_of_genes = original.number_of_genes;
   number_of_vectors = original.number_of_vectors;
   modified = original.modified;
   if (original.rep_vector!=NULL) {
      rep_vector = new Representation*[number_of_vectors];
      lookup = new Lookup[number_of_genes];

      for (i=0; i<number_of_vectors; i++) {
         rep_vector[i] = original.rep_vector[i]->clone();
      }

      for (i=0; i<number_of_genes; i++) {
         lookup[i] = original.lookup[i];
      }
   } else {
      rep_vector = NULL;
      lookup = NULL;
   }
}

Genotype::~Genotype(void)
{
   register unsigned int i;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Genotype::~Genotype(void)\n");
#endif /* DEBUG */


   if (rep_vector!=NULL) {
      for (i=0; i<number_of_vectors; i++) {
         delete rep_vector[i];
      }
      delete [] rep_vector;
      delete [] lookup;
   }
}

Genotype &Genotype::operator=(const Genotype &original)
{
   unsigned int i;

#ifdef DEBUG
   (void)fprintf(logFile, "\nsupport.cc/Genotype &Genotype::operator=(const Genotype &original): this==original is %d\n\n", this==&original);
#endif /* DEBUG */

   if (this==&original){ //Prevent self assignment
      return *this;
   }

/*** MPique 2014 This next block at one point was dumping core
 * and probably should be looked at more carefully in context.
 * Its job is to remove old rep from "A" in "A=B" expression
***/

   if (rep_vector!=NULL) {
      for (i=0; i<number_of_vectors; i++) {
         delete rep_vector[i];
      }
      delete [] rep_vector;
      delete [] lookup;
   }


   number_of_vectors = original.number_of_vectors;
   number_of_genes = original.number_of_genes;
//   modified = original.modified;
   modified = 1;

   if (original.rep_vector!=NULL) {
      rep_vector = new Representation *[number_of_vectors];
      lookup = new Lookup[number_of_genes];

      for (i=0; i<number_of_vectors; i++) {
         rep_vector[i] = original.rep_vector[i]->clone();
      }

      for (i=0; i<number_of_genes; i++) {
         lookup[i] = original.lookup[i];
      }
   } else {
      rep_vector = NULL;
      lookup = NULL;
   }
   return(*this);
}


void Genotype::write(const Element& value, const int gene_number) /* not const */
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Genotype::write(Element value, int gene_number=%d)\n", gene_number);
#endif /* DEBUG */

   modified = 1;
   rep_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Genotype::write(const unsigned char& value, const int gene_number) /* not const */
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Genotype::write(unsigned char value=%c, int gene_number=%d)\n", value, gene_number);
#endif /* DEBUG */

   modified = 1;
   rep_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Genotype::write(const int& value, const int gene_number) /* not const */
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Genotype::write(int value=%ld, int gene_number=%d)\n", value, gene_number);
#endif /* DEBUG */

   modified = 1;
   rep_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Genotype::write(ConstDouble value, const int gene_number) /* not const */
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Genotype::write(double value=%lf, int gene_number=%d)\n", value, gene_number);
#endif /* DEBUG */

   modified = 1;
   rep_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Genotype::write(const Representation &value, const int gene_number)
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Genotype::write(const Representation &value, int gene_number=%d)\n", gene_number);
#endif /* DEBUG */

   modified = 1;
//   *rep_vector[lookup[gene_number].vector] = value;
   *(rep_vector[gene_number]) = value;
}

Quat Genotype::readQuat() const
{
    Quat q;
    q.x = gread(3).real;
    q.y = gread(4).real;
    q.z = gread(5).real;
    q.w = gread(6).real;
    return q;
}

void Genotype::writeQuat( const Quat& q ) /* not const */
{
    write( q.x, 3 );
    write( q.y, 4 );
    write( q.z, 5 );
    write( q.w, 6 );
}

Quat Phenotype::readQuat() const
{
    Quat q;
    q.x = gread(3).real;
    q.y = gread(4).real;
    q.z = gread(5).real;
    q.w = gread(6).real;
    return q;
}

void Phenotype::writeQuat( const Quat& q ) /* not const */
{
    write( q.x, 3 );
    write( q.y, 4 );
    write( q.z, 5 );
    write( q.w, 6 );
}

Phenotype::Phenotype(const unsigned int init_number_of_dimensions, Representation **const init_value_vector, Eval *init_pevaluate)
: number_of_dimensions(init_number_of_dimensions), value_vector(init_value_vector), 
  value(0.0), evalflag(0), pevaluate(init_pevaluate)
{
   register unsigned int i, j, k;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Phenotype::Phenotype(unsigned int init_number_of_dimensions=%d, Representation **init_value_vector)\n", init_number_of_dimensions);
#endif /* DEBUG */


   number_of_points = 0;
   for (i=0; i<number_of_dimensions; i++) {
      number_of_points += value_vector[i]->number_of_points();
   }

   i = 0;
   lookup = new Lookup[number_of_points];
   for (j=0; j<number_of_dimensions; j++) {
      for (k=0; k<value_vector[j]->number_of_points(); k++) {
         assert ( i < number_of_points ); // mp!
         lookup[i].vector = j;
         lookup[i].index = k;
         i++;
      }
   }
}

Phenotype::Phenotype(const Phenotype &original)
{
   register unsigned int i;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Phenotype::Phenotype(const Phenotype &original)\n");
#endif /* DEBUG */


   number_of_dimensions = original.number_of_dimensions;
   number_of_points = original.number_of_points;
   evalflag = original.evalflag;
   value = original.value;
   pevaluate= original.pevaluate;

   if (original.value_vector!=NULL) {
      value_vector = new Representation *[number_of_dimensions];
      lookup = new Lookup[number_of_points];

      for (i=0; i<number_of_dimensions; i++) {
         value_vector[i] = original.value_vector[i]->clone();
      }
      for (i=0; i<number_of_points; i++) {
         lookup[i] = original.lookup[i];
      }
   } else {
      value_vector = NULL;
      lookup = NULL;
   }
}

Phenotype &Phenotype::operator=(const Phenotype &original)
{
   register unsigned int i;

#ifdef DEBUG
   (void)fprintf(logFile, "\nsupport.cc/Phenotype &Phenotype::operator=(const Phenotype &original): this==original is %d\n\n", this==&original);
#endif /* DEBUG */

   if (this==&original){ //Prevent self assignment
      return *this;
   }

/** MPique 2014 - removed for now as is dumping core
*/
   //  Do the destructors get called on each element of value_vector?
   if (value_vector!=NULL) {
      for (i=0; i<number_of_dimensions; i++) {
         delete value_vector[i];
      }
      delete [] value_vector;
      delete [] lookup;
   }
/*
****/

   number_of_dimensions = original.number_of_dimensions;
   number_of_points = original.number_of_points;
   value = original.value;
   evalflag = original.evalflag;
   pevaluate = original.pevaluate;

   if (original.value_vector!=NULL) {
      value_vector = new Representation *[number_of_dimensions];
      lookup = new Lookup[number_of_points];
   } else {
      value_vector = NULL;
      lookup = NULL;
   }

   for (i=0; i<number_of_dimensions; i++) {
      value_vector[i] = original.value_vector[i]->clone();
   }

   for (i=0; i<number_of_points; i++) {
      lookup[i] = original.lookup[i];
   }

   return(*this);
}

Phenotype::~Phenotype(void)
{
   register unsigned int i;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Phenotype::~Phenotype(void)\n");
#endif /* DEBUG */


   if (value_vector!=NULL) {
      for (i=0; i<number_of_dimensions; i++) {  
         delete value_vector[i];
      }
      delete [] value_vector;
      delete [] lookup;
   }
}

void Phenotype::write(const Element& value, const int gene_number) /* not const */
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Phenotype::write(Element value, int gene_number=%d)\n", gene_number);
#endif /* DEBUG */

   evalflag = 0;
   value_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Phenotype::write(const unsigned char& value, const int gene_number) /* not const */
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Phenotype::write(unsigned char value=%c, int gene_number=%d)\n", value, gene_number);
#endif /* DEBUG */

   evalflag = 0;
   value_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Phenotype::write(const int& value, const int gene_number) /* not const */
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Phenotype::write(int value=%ld, int gene_number=%d)\n", value, gene_number);
#endif /* DEBUG */

   evalflag = 0;
   value_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Phenotype::write(ConstDouble value, const int gene_number) /* not const */
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Phenotype::write(double value=%lf, int gene_number=%d)\n", value, gene_number);
#endif /* DEBUG */

   evalflag = 0;
   value_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Phenotype::write(const Representation &value, const int gene_number) /* not const */
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Phenotype::write(const Representation &value, int gene_number=%d)\n", gene_number);
#endif /* DEBUG */

   evalflag = 0;
//   *(value_vector[lookup[gene_number].vector]) = value;
   *(value_vector[gene_number]) = value;
}

double Phenotype::evaluate(const EvalMode& mode) /* not const */
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/double Phenotype::evaluate(EvalMode mode)\n");
#endif /* DEBUG */

// MPique - the pevaluate is an "Eval *" member of the Phenotype object
//  and this code calls its () operator, see eval.cc,
//  which calls make_state_from_ind and then its eval() function.
   switch(mode)
   {
      case Always_Eval:
         //value = ::evaluate(value_vector);
         value = (*pevaluate)(value_vector);
         evalflag = 1;
         break;
      case Always_Eval_Nonbond:
         //value = ::evaluate(value_vector, 1); // term=1 as total non-bonded energy, i.e. vdW+Hb+desolv
         value = (*pevaluate)(value_vector, 1); // term=1 as total non-bonded energy, i.e. vdW+Hb+desolv
         evalflag = 0;
         break;
      case Always_Eval_Elec:
         //value = ::pevaluate(value_vector, 2); // term=2 as total electrostatic energy
         value = (*pevaluate)(value_vector, 2); // term=2 as total electrostatic energy
         evalflag = 0;
         break;
      case Normal_Eval:
         if (!evalflag) {
            //value = ::evaluate(value_vector);
            value = (*pevaluate)(value_vector);
            evalflag = 1;
         }
         break;
      case Reset:
         evalflag = 0;
         break;
      default:
         stop("BUGCHECK: Phenotype::evaluate(const EvalMode& mode) Unknown mode!\n");
         break;
   }

   return(value);
}

State Phenotype::make_state(const int ntor) const
{
   State retval;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/State Phenotype::make_state(int ntor=%d)\n", ntor);
#endif /* DEBUG */


   retval.ntor = ntor;
   make_state_from_rep(value_vector, &retval, outlev, logFile);
   //make_state_from_rep(value_vector, &retval);
   return(retval);
}

Individual &Population::operator[](const int ind_num) const
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Individual &Population::operator[](int ind_num=%d)\n", ind_num);
#endif /* DEBUG */

   if ((ind_num<0)||(ind_num>=size)) {
	char error_message[200];
      prStr( error_message,  "ERROR: BUGCHECK support.cc/Trying to access %d, an out of bounds individual! (0<i<%d)\n", ind_num, size);
      //stop(error_message);
      return(heap[0]); // notreached
   } else {
      return(heap[ind_num]);
   }
}

State Individual::state(const int ntor) const
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/State Individual::state(int ntor=%d)\n", ntor);
#endif /* DEBUG */

   return(phenotyp.make_state(ntor));
}

void Individual::getMol(Molecule * /* not const */ returnedMol) const		// Steffen : bug candidate here
{
// Converts phenotype to mol's state and returns this individual's mol data.

    State molState;
    Molecule molcopy;

    molState = phenotyp.make_state(mol->S.ntor);
    molcopy = copyStateToMolecule(&molState, mol);
    returnedMol = &molcopy;							// Steffen : or here, molcopy is not returned
}

void Individual::printIndividualsState(FILE *const filePtr, const int ntor, const int detail) const
{
#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Individual::printIndividualsState(FILE *filePtr, int ntor=%d, int detail=%d)\n", ntor, detail);
#endif /* DEBUG */

    printState( filePtr, state(ntor), detail ); 
}

void Individual::incrementAge(void) /* not const */
{
    ++age;
}

Population &Population::operator=(const Population &original)
{
   register int i;

#ifdef DEBUG
   (void)fprintf(logFile, "\nsupport.cc/Population &Population::operator=(const Population &original):this=original is %d\n\n", this==&original);
#endif /* DEBUG */

   if (this==&original){ //Prevent self assignement
      return *this;
   }

   if (heap!=NULL) {
      delete [] heap;
   }

   size = original.size;
   heap = new Individual[size];
   lhb = original.lhb;
   for (i=0; i<size; i++) {
      heap[i] = original.heap[i];
   }
   // MPique 2014 TODO should this also copy the end_of_branch info?

   return(*this);
}
