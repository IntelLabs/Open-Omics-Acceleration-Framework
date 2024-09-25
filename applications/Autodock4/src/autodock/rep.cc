/*

 $Id: rep.cc,v 1.21 2014/06/12 01:44:08 mp Exp $

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
     The methods associated with the Representation class hierarchy.

				rsh 9/95
********************************************************************/
 
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include "rep.h"
#include "ranlib.h"
#include "structs.h"
#include "stop.h"

extern FILE *logFile; // DEBUG and bug check messages only
extern int debug;

//______________________________________________________________________________
//
//  Initializations
int IntVector::low = -INT_MAX/4;
int IntVector::high = INT_MAX/4;
/* A nonstatic data member cannot be defined outside its class:
 * Real RealVector::low = REALV_LOW;
 * Real RealVector::high = REALV_HIGH;
 */
//  For now assume that normalize handles this constraint
Real ConstrainedRealVector::low = REALV_LOW;
Real ConstrainedRealVector::high = REALV_HIGH;
double ConstrainedRealVector::sum = 1.0;
Real BitVector::one_prob = 0.5;

//______________________________________________________________________________
//
//  The member functions for the canonical base classes
//______________________________________________________________________________


//______________________________________________________________________________
//
//  This constructor is used to generate the initial (random) instances
//  of an integer vector.
IntVector::IntVector(const int number_of_els)
: Representation(number_of_els)
{
   register int i;

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/IntVector::IntVector(int number_of_els=%d) \n",number_of_els);
#endif /* DEBUG */

   mytype = T_IntV;
   vector = new int[number_of_els];
   for (i=0; i<number_of_els; i++) {
      vector[i] = ignuin(low, high);
   }
}

//______________________________________________________________________________
//
IntVector::IntVector(const int num_els, const int init_low, const int init_high)
: Representation(num_els)
{
   register int i;

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/IntVector::IntVector(int num_els=%d, int init_low=%d, int init_high=%d) \n",num_els,init_low,init_high);
#endif /* DEBUG */


   mytype = T_IntV;
   vector = new int[num_els];
   for (i=0; i<num_els; i++) {
      vector[i] = ignuin(init_low, init_high);
   }
}

//______________________________________________________________________________
//
//  This constructor does an actual copy of the vector.
//  We could make gains by doing reference counting, but
//  that's for the future.
IntVector::IntVector(const IntVector &original)
: Representation(original.number_of_pts)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/IntVector::IntVector(const IntVector &original) \n");
#endif /* DEBUG */

   mytype = T_IntV;
   if (original.vector!=NULL) {
      vector = new int[number_of_pts];
   } else {
      vector = NULL;
   }

   for (register unsigned int i=0; i<number_of_pts; i++) {
      vector[i] = original.vector[i];
   }
}

//______________________________________________________________________________
//
void IntVector::write(const unsigned char& value, const int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void IntVector::write(unsigned char value=%c, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing a Bit to an Int!\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= \"%c\", gene= %d\n", value, gene); // used to be "stderr"
}

//______________________________________________________________________________
//
void IntVector::write(const int& value, const int gene) /* not const */ 
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void IntVector::write(int value=%d, int gene=%d)` \n",value,gene);
#endif /* DEBUG */

   if (value<low) {
      vector[gene] = low;
   } else if (value>high) {
      vector[gene] = high;
   } else {
      vector[gene] = value;
   }
}

//______________________________________________________________________________
//
void IntVector::write(ConstDouble value, const int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void IntVector::write(double value=%lf, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing a Real to an Int!\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= %lf, gene= %d\n", value, gene); // used to be "stderr"
}

//______________________________________________________________________________
//
void IntVector::write(const Element& value, const int gene) /* not const */
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void IntVector::write(const Element value, int gene=%d) \n",gene);
#endif /* DEBUG */

   if (value.integer<low) {
      (void)fprintf(logFile,"Writing out-of-bounds Int!\n"); // used to be "stderr"
      vector[gene] = low;
   } else if (value.integer>high) {
      (void)fprintf(logFile,"Writing out-of-bounds Int!\n"); // used to be "stderr"
      vector[gene] = high;
   } else {
      vector[gene] = value.integer;
   }
}


//______________________________________________________________________________
//
const Element IntVector::gene(const unsigned int gene_number) const
{
   Element retval;

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const Element IntVector::gene(unsigned int gene_number=%d) const \n",gene_number);
#endif /* DEBUG */


   if (gene_number>=number_of_pts) {
      char error_message[200];
      (void)sprintf(error_message, "ERROR: BUGCHECK: Trying to access an out-of-bounds IntVector gene! (gene_number=%d >= number_of_pts=%d)\n", gene_number, number_of_pts); // used to be "stderr"
      stop(error_message);
      return(retval); // NOTREACHED
   } else {
      retval.integer = vector[gene_number];
      return(retval);  // typecast int as Element
   }
}

//______________________________________________________________________________
//
const void *IntVector::internals(void) const
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const void *IntVector::internals(void) const \n");
#endif /* DEBUG */

   return((void *)(&vector[0]));
}

//______________________________________________________________________________
//
Representation &IntVector::operator=(const Representation &original)
{
   register unsigned int i;

#ifdef DEBUG
    (void)fprintf(logFile, "\nrep.cc/Representation &IntVector::operator=(const Representation &original) \n");
#endif /* DEBUG */


   const int *const array = (int *)original.internals();
   if (original.type()==T_IntV) {
      number_of_pts = original.number_of_points();
      if (vector!=NULL) {
         delete [] vector;
      }

      if (array!=NULL) {
         vector = new int[number_of_pts];
      } else {
         vector = NULL;
      }

      for (i=0; i<number_of_pts; i++) {
         vector[i] = array[i];
      }
   } else {
      stop("Unable to invoke Representation &IntVector operator= because Representations don't match!\n");
   }

   return(*this);
}

//______________________________________________________________________________
//
//  This constructor is used to initialize the starting population, 
//  with random values between REALV_LOW and REALV_HIGH.
//
RealVector::RealVector(/* not const */ int num_els)
: Representation(num_els)
{
#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/RealVector::RealVector(int num_els=%d) \n",num_els);
#endif /* DEBUG */
   mytype = T_RealV;
   low = REALV_LOW;
   high = REALV_HIGH;
   vector = new double[num_els];
   for (; --num_els>=0;) {
      vector[num_els] = double(genunf(low, high));
#ifdef DEBUG
      (void)fprintf(logFile, "rep.cc/RealVector::RealVector(num_els)   vector[num_els] = %.3f\n", vector[num_els] );
#endif /* DEBUG */
   }
} // RealVector::RealVector(int num_els)

//______________________________________________________________________________
//
//  This constructor is used to initialize the starting population, 
//  with elements of the vector set to random values between 
//  the user-specified values, init_low and init_high.
//
RealVector::RealVector(/* not const */ int num_els, ConstDouble init_low, ConstDouble init_high)
: Representation(num_els)
{
#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/RealVector::RealVector(int num_els=%d, double init_low=%lf, double init_high=%lf) \n",num_els,init_low,init_high);
#endif /* DEBUG */
   mytype = T_RealV;
   low = init_low;
   high = init_high;
   vector = new double[num_els];
   for (; --num_els>=0;) {
      vector[num_els] = double(genunf(init_low, init_high));
#ifdef DEBUG
      (void)fprintf(logFile, "rep.cc/RealVector::RealVector(num_els, init_low, init_high)   vector[num_els] = %.3f\n", vector[num_els] );
#endif /* DEBUG */
   }
} // RealVector::RealVector(int num_els, double init_low, double init_high)

//______________________________________________________________________________
//
//  This constructor is used to initialize the starting population, 
//  setting the value of the second and all remaining elements in the 
//  vector (if any) to a uniformly-distributed random number between 
//  the user-specified bounds, init_low and init_high;
//  but the first value in this vector is set to the value supplied as the last argument.
//  This is useful for specifying an initial axis-angle rotation angle.
//
RealVector::RealVector(const int num_els, ConstDouble init_low, ConstDouble init_high, ConstDouble init_first_value)
: Representation(num_els)
{
#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/RealVector::RealVector(int num_els=%d, double init_low=%lf, double init_high=%lf, double init_first_value=%lf) \n",num_els,init_low,init_high,init_first_value);
#endif /* DEBUG */
   register int i=0;
   mytype = T_RealV;
   low = init_low;
   high = init_high;
   vector = new double[num_els];
   // Set the first element in the vector to "init_first_value":
   vector[i] = init_first_value;
#ifdef DEBUG
      (void)fprintf(logFile, "rep.cc/RealVector::RealVector(i, init_low, init_high, init_first_value)   vector[i] = %.3f\n", vector[i] );
#endif /* DEBUG */
   // Set the second and remaining elements in the vector, if any:
   for (i=1; i<num_els; i++) {
      vector[i] = double(genunf(init_low, init_high));
#ifdef DEBUG
      (void)fprintf(logFile, "rep.cc/RealVector::RealVector(i, init_low, init_high, init_first_value)   vector[i] = %.3f\n", vector[i] );
#endif /* DEBUG */
   }
} // RealVector::RealVector(int num_els, double init_low, double init_high, double init_first_value)

//______________________________________________________________________________
//
//  This constructor is used to initialize the starting population, 
//  with elements of a vector of length 3 being set to the user specified values
//  nx, ny, and nz.
//  This is useful for specifying an initial orientation's axis components
//
RealVector::RealVector( const int num_els, ConstDouble init_low, ConstDouble init_high, ConstDouble nx, ConstDouble ny, ConstDouble nz )
: Representation(num_els)
{
#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/RealVector::RealVector(double nx=%lf, double ny=%lf, double nz=%lf, int num_els=%d) \n", nx, ny, nz, num_els);
#endif /* DEBUG */

   mytype = T_RealV;
   low = init_low;
   high = init_high;
   vector = new double[3];
   // Set the unit vector.
   vector[0] = nx;
   vector[1] = ny;
   vector[2] = nz;
#ifdef DEBUG
   (void)fprintf(logFile, "rep.cc/RealVector::RealVector(nx,ny,nz,num_els)   vector[0] = %.3f\n", vector[0] );
   (void)fprintf(logFile, "rep.cc/RealVector::RealVector(nx,ny,nz,num_els)   vector[1] = %.3f\n", vector[1] );
   (void)fprintf(logFile, "rep.cc/RealVector::RealVector(nx,ny,nz,num_els)   vector[2] = %.3f\n", vector[2] );
#endif /* DEBUG */
} // RealVector::RealVector(double nx, double ny, double nz, int num_els)

//______________________________________________________________________________
//
//  This constructor is used to initialize the starting population, 
//  with elements of a vector of length 4 being set to the user specified values
//  w, x, y, z
//  This is useful for specifying an initial quaternion rotation' unit vector.
//
RealVector::RealVector( const int num_els, ConstDouble init_low, ConstDouble init_high, ConstDouble x, ConstDouble y, ConstDouble z, ConstDouble w)
: Representation(num_els)
{
#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/RealVector::RealVector(int num_els=%d, double x=%lf, double y=%lf, double z=%lf, double w=%lf) \n", num_els, x, y, z, w );
#endif /* DEBUG */

   mytype = T_RealV;
   low = init_low;
   high = init_high;
   vector = new double[4];
   // Set the x,y,z,w components of the quaternion.
   vector[0] = x;
   vector[1] = y;
   vector[2] = z;
   vector[3] = w;
#ifdef DEBUG
   (void)fprintf(logFile, "rep.cc/RealVector::RealVector(num_els,x,y,z,w)   vector[0] = %.3f\n", vector[0] );
   (void)fprintf(logFile, "rep.cc/RealVector::RealVector(num_els,x,y,z,w)   vector[1] = %.3f\n", vector[1] );
   (void)fprintf(logFile, "rep.cc/RealVector::RealVector(num_els,x,y,z,w)   vector[2] = %.3f\n", vector[2] );
   (void)fprintf(logFile, "rep.cc/RealVector::RealVector(num_els,x,y,z,w)   vector[3] = %.3f\n", vector[3] );
#endif /* DEBUG */
}
//______________________________________________________________________________
//
//  Do a deep copy of the original
//
RealVector::RealVector(const RealVector &original)
: Representation(original.number_of_pts)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/RealVector::RealVector(const RealVector &original) \n");
#endif /* DEBUG */

   mytype = T_RealV;
   low =  original.low;
   high = original.high;
   if (original.vector!=NULL) {
      vector = new double[original.number_of_pts];
   } else {
      vector = NULL;
   }

   for (register unsigned int i=0; i<original.number_of_pts; i++) {
      vector[i] = original.vector[i];
#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/i=%d, original.number_of_pts=%d, vector[%d]= %.3f\n",i, original.number_of_pts, i, vector[i]);
#endif /* DEBUG */
   }
}

//______________________________________________________________________________
//
void RealVector::write(const unsigned char& value, const int gene) /* not const ... in sibling classes */
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void RealVector::write(unsigned char value=%c, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing a Bit to a Real!\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= \"%c\", gene= %d\n", value, gene); // used to be "stderr"
}

//______________________________________________________________________________
//
void RealVector::write(const int& value, const int gene) /* not const ... in sibling classes */
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void RealVector::write(int value=%d, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing an Int to a Real!\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= %ld, gene= %d\n", (long)value, gene); // used to be "stderr"
}

//______________________________________________________________________________
//
void RealVector::write(ConstDouble value, const int gene) /* not const */
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void RealVector::write(double value=%lf, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   if (value<low) {
      // if (debug > 0) {
          // (void)fprintf(logFile,"WARNING:  Writing out of bounds Real!  value (%lf) too low (%lf)\n",value,low); // used to be "stderr"
      // }
      vector[gene] = low;
   } else if (value>high) {
      // if (debug > 0) {
          // (void)fprintf(logFile,"WARNING:  Writing out of bounds Real!  value (%lf) too high (%lf)\n",value,high); // used to be "stderr"
      // }
      vector[gene] = high;
   } else {
      vector[gene] = value;
   }
}


//______________________________________________________________________________
//
void RealVector::write(const Element& value, const int gene) /* not const */
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void RealVector::write(const Element value, int gene=%d) \n",gene);
#endif /* DEBUG */

   if (value.real<low) {
      vector[gene] = low;
   } else if (value.real>high) {
      vector[gene] = high;
   } else {
      vector[gene] = value.real;
   }
}

//______________________________________________________________________________
//
const Element RealVector::gene(const unsigned int gene_number) const
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const Element RealVector::gene(unsigned int gene_number=%d) const \n",gene_number);
#endif /* DEBUG */

   Element retval;

   if (gene_number>=number_of_pts) {
      (void)fprintf(logFile,"ERROR: Trying to access an out-of-bounds RealVector gene (gene_number=%d >= number_of_pts=%d)\n", gene_number, number_of_pts); // used to be "stderr"
      retval.real = 0.0;
      return(retval);
   } else {
#ifdef DEBUG
      pr( logFile, "rep.cc /  retval.real = vector[gene_number=%d] = %lf\n", gene_number, vector[gene_number] );
#endif
      retval.real = vector[gene_number];
#ifdef DEBUG
      pr( logFile, "rep.cc / retval.real = %lf\n", retval.real );
#endif
      return(retval);
   }
}

//______________________________________________________________________________
//
const void *RealVector::internals(void) const
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const void *RealVector::internals(void) const \n");
#endif /* DEBUG */

   return((void *)(&vector[0]));
}

//______________________________________________________________________________
//
Representation &RealVector::operator=(const Representation &original)
{

#ifdef DEBUG
    (void)fprintf(logFile, "\nrep.cc/Representation &RealVector::operator=(const Representation &original) \n");
#endif /* DEBUG */

   register unsigned int i;
   double *array;

   if (original.type()==T_RealV) {
      low = ((const RealVector &)original).low;
      high = ((const RealVector &)original).high;
      array = (double *)original.internals();
      number_of_pts = original.number_of_points();
      if (vector!=NULL) {
         delete [] vector;
      }

      if (array!=NULL) {
         vector = new double[number_of_pts];
      } else {
         vector = NULL;
      }

      for (i=0; i<number_of_pts; i++) {
         vector[i] = array[i];
      }
   } else {
      (void)fprintf(logFile,"Unable to invoke operator= because Representations don't match!\n"); // used to be "stderr"
   }

   return(*this);
}

//______________________________________________________________________________
//
ConstrainedRealVector::ConstrainedRealVector(/* not const */ int num_els)
:  Representation(num_els)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/ConstrainedRealVector::ConstrainedRealVector(int num_els) \n");
#endif /* DEBUG */

   mytype = T_CRealV;
   normalized = 0;
   vector = new double[num_els];
   for (; --num_els>=0;) {
      vector[num_els] = double(genunf(low, high));
   }

   normalize();
}

//______________________________________________________________________________
//
ConstrainedRealVector::ConstrainedRealVector(/* not const */ int num_els, ConstDouble init_low, ConstDouble init_high)
:  Representation(num_els)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/ConstrainedRealVector::ConstrainedRealVector(int num_els=%d, double init_low=%lf, double init_high=%lf) \n",num_els,init_low,init_high);
#endif /* DEBUG */

   mytype = T_CRealV;
   normalized = 0;
   low = init_low;
   high = init_high;
   vector = new double[num_els];
   for (; --num_els>=0;) {
      vector[num_els] = double(genunf(init_low, init_high));
   }

   normalize();
}

//______________________________________________________________________________
//
ConstrainedRealVector::ConstrainedRealVector(const ConstrainedRealVector &original)
:  Representation(original.number_of_pts)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/ConstrainedRealVector::ConstrainedRealVector(const ConstrainedRealVector &original) \n");
#endif /* DEBUG */

   mytype = T_CRealV;
   normalized = original.normalized;
   if (original.vector != NULL) {
      vector = new double[original.number_of_pts];
   } else {
      vector = NULL;
   }

   for (register unsigned int i=0; i<original.number_of_pts; i++) {
      vector[i] = original.vector[i];
   }
}

//______________________________________________________________________________
//
void ConstrainedRealVector::write(const unsigned char& value, const int gene) /* not const, inherited */
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void ConstrainedRealVector::write(unsigned char value=%c, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing a Bit to a Constrained Real\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= \"%c\",  gene= %d\n", value, gene); // used to be "stderr"
}

//______________________________________________________________________________
//
void ConstrainedRealVector::write(const int& value, const int gene) 
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void ConstrainedRealVector::write(int value=%ld, int gene=%d) \n",(long)value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing an Integer to a Constrained Real\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= %ld, gene= %d\n",(long)value,gene); // used to be "stderr"
}

//______________________________________________________________________________
//
void ConstrainedRealVector::write(ConstDouble value, const int gene) /* not const */
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void ConstrainedRealVector::write(double value=%lf, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   if (value<low) {
      (void)fprintf(logFile,"Writing out-of-bounds Constrained Real\n"); // used to be "stderr"
      vector[gene] = low;
   } else if (value>high) {
      (void)fprintf(logFile,"Writing out-of-bounds Constrained Real\n"); // used to be "stderr"
      vector[gene] = high;
   } else {
      vector[gene] = value;
   }

   normalized = 0;
}


//______________________________________________________________________________
//
void ConstrainedRealVector::write(ConstDouble a, ConstDouble b, ConstDouble c, ConstDouble d)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void ConstrainedRealVector::write( double a=%lf, double b=%lf, double c=%lf, double d=%lf ) \n", a, b, c, d );
#endif /* DEBUG */

#define clamp_and_set_vector( value, gene ) \
   if (value<low) { \
      (void)fprintf(logFile,"Writing out-of-bounds Constrained Real\n"); \
      vector[gene] = low; \
   } else if (value>high) { \
      (void)fprintf(logFile,"Writing out-of-bounds Constrained Real\n"); \
      vector[gene] = high; \
   } else { \
      vector[gene] = value; \
   }

   assert( number_of_pts >= 4 );

   clamp_and_set_vector( a, 0 );
   clamp_and_set_vector( b, 1 );
   clamp_and_set_vector( c, 2 );
   clamp_and_set_vector( d, 3 );

   normalize();
}
//______________________________________________________________________________
//
void ConstrainedRealVector::write(const Element& value, const int gene) /* not const */
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void ConstrainedRealVector::write(const Element value, int gene=%d) \n",gene);
#endif /* DEBUG */

   if (value.real<low) {
      vector[gene] = low;
   } else if (value.real>high) {
      vector[gene] = high;
   } else {
      vector[gene] = value.real;
   }

   normalized = 0;
}

//______________________________________________________________________________
//
void ConstrainedRealVector::normalize(void) /* not const */
{
#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void ConstrainedRealVector::normalize(void) const \n");
#endif /* DEBUG */

   if (!normalized) {
      register unsigned int i;
      register double tempsum = 0.0, hypotenuse;

      for (i=0; i<number_of_pts; i++) {
         tempsum += vector[i] * vector[i];
      }

      if ((tempsum - sum  >  ACCURACY) || (sum - tempsum  >  ACCURACY)) {
         hypotenuse = sqrt(tempsum);
         // normalize the vector[]
         for (i=0; i<number_of_pts; i++) {
            vector[i] /= hypotenuse;
         }
      }

      //normalized = 1;
      set_normalized_true();
   }
}


//______________________________________________________________________________
//
const Element ConstrainedRealVector::gene(const unsigned int gene_number) const
{
   Element retval;

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const Element ConstrainedRealVector::gene(unsigned int gene_number=%d) const \n",gene_number);
#endif /* DEBUG */


   if (gene_number>=number_of_pts) {
      (void)fprintf(logFile,"ERROR: Trying to access an out-of-bounds ConstrainedRealVector gene (gene_number=%d >= number_of_pts=%d)\n", gene_number, number_of_pts); // used to be "stderr"
      retval.real = 0.0;
      return(retval);
   } else {
      // normalize();  // cannot normalize because gene(int) is const
      retval.real = vector[gene_number];
      return(retval);
   }
}

//______________________________________________________________________________
//
const void *ConstrainedRealVector::internals(void) const
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const void *ConstrainedRealVector::internals(void) const \n");
#endif /* DEBUG */

   return((void *)(&vector[0]));
}

//______________________________________________________________________________
//
Representation &ConstrainedRealVector::operator=(const Representation &original)
{
   register unsigned int i;
   double *array;

#ifdef DEBUG
    (void)fprintf(logFile, "\nrep.cc/Representation &ConstrainedRealVector::operator=(const Representation &original) \n");
#endif /* DEBUG */


   array = (double *)original.internals();
   if (original.type()==T_CRealV) {
      number_of_pts = original.number_of_points();
      normalized = original.is_normalized();
      if (vector!=NULL) {
         delete [] vector;
      }
      
      if (array!=NULL) {
         vector = new double[number_of_pts];
      } else {
         vector = NULL;
      }
      
      for (i=0; i<number_of_pts; i++) {
         vector[i] = array[i];
      }
   } else {
      (void)fprintf(logFile,"Unable to invoke operator= because Representations don't match!\n"); // used to be "stderr"
   }

   return(*this);
}

//______________________________________________________________________________
//
//  This constructor is used to initialize the first
//  generation of any particular bitvector.  Right
//  now bits are assumed to be unsigned chars.
BitVector::BitVector(/* not const */ int num_els)
: Representation(num_els)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/BitVector::BitVector(int num_els=%d) \n",num_els);
#endif /* DEBUG */

   mytype = T_BitV;
   vector = new unsigned char[num_els];
   for (; --num_els>=0;) {
      vector[num_els] = ((ranf()<one_prob)? 1 : 0);
   }
}

//______________________________________________________________________________
//
BitVector::BitVector(/* not const */ int num_els, ConstReal prob)
: Representation(num_els)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/BitVector::BitVector(int num_els=%d, Real prob=%f) \n",num_els,prob);
#endif /* DEBUG */

   mytype = T_BitV;
   vector = new unsigned char[num_els];
   for (; --num_els>=0;) {
      vector[num_els] = ((ranf()<prob)? 1 : 0);
   }
}

//______________________________________________________________________________
//
//  There are probably better ways of doing this, e.g.
//  using memcpy()
BitVector::BitVector(const BitVector &original)
: Representation(original.number_of_pts)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/BitVector::BitVector(const BitVector &original) \n");
#endif /* DEBUG */

   mytype = T_BitV;
   if (original.vector!=NULL) {
      vector = new unsigned char[number_of_pts];
   } else {
      vector = NULL;
   }

   for (register unsigned int i=0; i<number_of_pts; i++) {
      vector[i] = original.vector[i];
   }
}

//______________________________________________________________________________
//
void BitVector::write(const unsigned char& value, const int gene) /* not const */
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void BitVector::write(unsigned char value=%c, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   vector[gene] = value;
}

//______________________________________________________________________________
//
void BitVector::write(const int& value, const int gene) 
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void BitVector::write(int value=%d, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing Int to Bit!\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= %ld, gene= %d\n",(long)value,gene); // used to be "stderr"
}

//______________________________________________________________________________
//
void BitVector::write(ConstDouble value, const int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void BitVector::write(double value=%lf, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing Real to Bit!\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= %lf, gene= %d\n",value,gene); // used to be "stderr"
}


//______________________________________________________________________________
//
void BitVector::write(const Element& value, const int gene) /* not const */
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void BitVector::write(const Element value, int gene=%d) \n",gene);
#endif /* DEBUG */

   vector[gene] = value.bit;
}


//______________________________________________________________________________
//
const Element BitVector::gene(const unsigned int gene_number) const
{
   Element retval;

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const Element BitVector::gene(unsigned int gene_number=%d) const \n",gene_number);
#endif /* DEBUG */


   if (gene_number>=number_of_pts) {
      (void)fprintf(logFile,"ERROR: Trying to access an out-of-bounds BitVector gene (gene_number=%d >= number_of_pts=%d)\n", gene_number, number_of_pts); // used to be "stderr"
      retval.bit = 0;
      return(retval);
   } else {
      retval.bit = vector[gene_number];
      return(retval);
   }
}

//______________________________________________________________________________
//
const void *BitVector::internals(void) const
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const void *BitVector::internals(void) const \n");
#endif /* DEBUG */

   return((void *)(&vector[0]));
}

//______________________________________________________________________________
//
Representation &BitVector::operator=(const Representation &original)
{
   register unsigned int i;
   unsigned char *array;

#ifdef DEBUG
    (void)fprintf(logFile, "\nrep.cc/Representation &BitVector::operator=(const Representation &original) \n");
#endif /* DEBUG */


   array = (unsigned char *)original.internals();
   if (original.type()==T_BitV) {
      if (vector!=NULL) {
         delete [] vector;
      }

      number_of_pts = original.number_of_points();
      if (array!=NULL) {
         vector = new unsigned char[number_of_pts];
      } else {
         vector = NULL;
      }

      for (i=0; i<number_of_pts; i++) {
         vector[i] = array[i];
      }
   } else {
      (void)fprintf(logFile,"Unable to invoke operator= because Representations don't match!\n"); // used to be "stderr"
   }

   return(*this);
}

//______________________________________________________________________________
//
