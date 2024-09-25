/*

 $Id: rep.h,v 1.20 2014/06/12 01:44:08 mp Exp $

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

//  These are the classes  associated with the Representation class
//  hierarchy.  The Representation class is meant to be a generic
//  place holder for any type of Representation that a user might
//  need to build a problem out of.  By the way, these derived 
//  (Representation) classes look like perfect candidates for templates 
//  We need to make sure that the const pointers
//  are never used to indirectly change values!
//  rsh 07/08/95
//  
//  The type Element was added to take care of all of the void *
//  rsh 02/23/96

#ifndef _REP_H
#define _REP_H

#include "rep_constants.h"

#include <stdio.h>
#include "structs.h"
#include "typedefs.h"

#if 0
#define ACCURACY      0.001
#define REALV_LOW  -999.999 //gmm, 2003-11-11
#define REALV_HIGH  999.999 //gmm, 2003-11-11
#else
const Real ACCURACY   =    0.001;
const Real REALV_LOW  = -999.999; //gmm, 2003-11-11
const Real REALV_HIGH  = 999.999; //gmm, 2003-11-11
#endif

// const Real REALV_LOW  = -100.0 //gmm, 2003-10-13
// const Real REALV_HIGH =  100.0 //gmm, 2003-10-13
// const Real REALV_LOW  = -3.14159265358979323846 //gmm, 1998-07-08
// const Real REALV_HIGH =  3.14159265358979323846 //gmm, 1998-07-08

enum RepType { T_BASE, T_IntV, T_RealV, T_CRealV, T_BitV, T_Orientation };

typedef union 
{
   double real;
   int integer;
   unsigned char bit;
} Element;

class Representation
{
   protected:
      unsigned int number_of_pts;
      RepType mytype;
      unsigned char normalized; // =1 means the vector's normalized

   public:
      Representation(void);
      Representation(const unsigned int);
      virtual ~Representation(void);
      virtual Representation &operator=(const Representation &) = 0;
      unsigned int number_of_points(void) const;
      int is_normalized(void) const;
      void set_normalized_true(void);
      void set_normalized_false(void);
      virtual RepType type(void) const; 
      virtual void write(const unsigned char& value, const int) = 0; /* not const, e.g. in bitvector */
      virtual void write(const int& value,  const int) = 0;
      virtual void write(ConstDouble value,          const int) = 0; /* not const, e.g. in RealVector */
      virtual void write(const Element& value,       const int) = 0;
      virtual const Element gene(const unsigned int) const = 0;
      virtual Representation *clone(void) const = 0;
      virtual const void *internals(void) const = 0;
};

class IntVector : public Representation
{
//   friend void debug(IntVector &);
   protected:
      static int low, high;
      int *vector;

      const void *internals(void) const;
      Representation *clone(void) const;

   public:
      IntVector(void);
      IntVector(const int);
      IntVector(const int, int *const);
      IntVector(const int, const int, const int);
      IntVector(const IntVector &);
      ~IntVector(void);
      void write(const unsigned char& value, const int); /* not const - inherited */
      void write(const int&  value, const int); /* not const */
      void write(ConstDouble          value, const int);
      void write(const Element&       value, const int); /* not const */
      Representation &operator=(const Representation &);
      const Element gene(const unsigned int) const;
};

class RealVector : public Representation
{
//   friend void debug(RealVector &);
   protected:
      Real high, low;
      double *vector;

      const void *internals(void) const;
      Representation *clone(void) const;

   public:
      RealVector(void)
      : Representation(0), 
        high(REALV_HIGH), low(REALV_LOW), 
        vector((double *)NULL)
        { 
            mytype = T_RealV;
        }

      RealVector(/* not const */ int);
      RealVector(const int, double *const);
      RealVector(/* not const */ int, ConstDouble, ConstDouble);
      // Use this to set the first value in the vector--useful for random quaternions
      RealVector(const int, ConstDouble, ConstDouble, ConstDouble);
      RealVector(const int, ConstDouble, ConstDouble, ConstDouble, ConstDouble, ConstDouble, ConstDouble);  // sets a quaternion's x,y,z,w values
      // Use this to create a vector of length 3 with these values--useful for random quaternions
      RealVector(const int, ConstDouble, ConstDouble, ConstDouble, ConstDouble, ConstDouble);
      RealVector(const RealVector &);
      ~RealVector(void);
      void write(const unsigned char&, const int) /* not const - inherited */;
      void write(const int&, const int); /* not const, e.g. in IntVector */
//    void write(const void *, int);
      void write(const Element&, const int) /* not const */;
      Representation &operator=(const Representation &);
//    const void *gene(unsigned int) const;
#ifdef DEBUG
      /*
      inline const Element gene(unsigned int gene_number) const
      {
          Element retval;
          retval.real = vector[gene_number];
          return retval;
      }
      inline void write(double value, int gene)
      {
          value = value < low  ? low  : 
                  value > high ? high :
                  value;
          vector[gene] = value;
      }
      */
      // non-inlined versions, possibly with range checking
      void write(ConstDouble value, const int gene) /* not const */;
      const Element gene(const unsigned int) const;
#else
      // non-inlined versions, possibly with range checking
      void write(ConstDouble value, const int gene) /* not const */;
      const Element gene(unsigned int) const;
#endif
};

//  Maybe this should be derived from RealVector
class ConstrainedRealVector : public Representation
{
//   friend debug(ConstrainedRealVector &);
   protected:
      static Real high, low;
      static double sum;
      double *vector;

      const void *internals(void) const;
      Representation *clone(void) const;
      void normalize(void); /* not const */

   public:
      ConstrainedRealVector(void);
      ConstrainedRealVector(/* not const */ int);
      ConstrainedRealVector(int, double *const);
      ConstrainedRealVector(/* not const */ int, ConstDouble, ConstDouble);
      ConstrainedRealVector(const ConstrainedRealVector &);
      ~ConstrainedRealVector(void);
      void write(const unsigned char& value, const int gene); /* not const - inherited */
      void write(const int& value,  const int gene); /* not const - e.g. in IntVector */
      void write(ConstDouble value,          const int gene); /* not const */
      void write(const Element&,             const int gene); /* not const */
      void write(ConstDouble a, ConstDouble b, ConstDouble c, ConstDouble d);
      Representation &operator=(const Representation &);
      const Element gene(const unsigned int) const;
};

class BitVector : public Representation
{
//   friend void debug(BitVector &);
   protected:
      static Real one_prob;
      unsigned char *vector;

      const void *internals(void) const;
      Representation *clone(void) const;

   public:
      BitVector(void);
      BitVector(/* not const */ int);
      BitVector(int, unsigned char *const);
      BitVector(/* not const */ int, ConstReal);
      BitVector(const BitVector &);
      ~BitVector(void);
      void write(const unsigned char&, const int); /* not const */
      void write(const int&, const int); /* not const - inherited */
      void write(ConstDouble, const int); /* not const - inherited */
//      void write(const void *, int);
      void write(const Element&, const int); /* not const */
      Representation &operator=(const Representation &);
//      const void *gene(unsigned int) const;
      const Element gene(const unsigned int) const;
};

/**************************************************************************
      Inline Functions
**************************************************************************/

inline Representation::Representation(void)
    : number_of_pts(0) 
{
    mytype = T_BASE;
}

inline Representation::Representation(const unsigned int pts)
    : number_of_pts(pts) 
{
}

inline Representation::~Representation(void)
{
}

inline unsigned int Representation::number_of_points(void) const
{
   return(number_of_pts);
}

inline int Representation::is_normalized(void) const
{
   return(normalized);
}

inline void Representation::set_normalized_true(void)
{
    normalized = 1;
}

inline void Representation::set_normalized_false(void)
{
    normalized = 0;
}

inline RepType Representation::type(void) const
{
   return(mytype);
}

inline Representation *Representation::clone(void) const
{
   return(NULL);
}

inline IntVector::IntVector(void)
: Representation(0)
{
   vector = (int *)NULL;
   mytype = T_IntV;
}

//  This constructor does a shallow copy of the array
inline IntVector::IntVector(const int num_els, int *const array)
: Representation(num_els), vector(array)
{
    mytype = T_IntV;
}

inline IntVector::~IntVector(void)
{
   if(vector!=(int *)NULL)
   {
      delete [] vector;
   }
}

inline Representation *IntVector::clone(void) const
{
   return(new IntVector(*this));
}

//  This performs a shallow copy of the array
inline RealVector::RealVector(const int num_els, double *const array)
: Representation(num_els), high(REALV_HIGH), low(REALV_LOW), 
        vector(array)
{
    mytype = T_RealV;
}

inline RealVector::~RealVector(void)
{
   if(vector!=(double *)NULL)
   {
      delete [] vector;
   }
}

inline Representation *RealVector::clone(void) const
{
   return(new RealVector(*this));
}

inline ConstrainedRealVector::ConstrainedRealVector(void)
: Representation(0)
{
   normalized = 1;
   vector = (double *)NULL;
   mytype = T_CRealV;
}

inline ConstrainedRealVector::ConstrainedRealVector(const int num_els, double *const array)
:  Representation(num_els)
{
   normalized = 0;
   vector = array;
   mytype = T_CRealV;
}

inline ConstrainedRealVector::~ConstrainedRealVector(void)
{
   if(vector!=(double *)NULL)
   { 
      delete [] vector;
   }
}

inline Representation *ConstrainedRealVector::clone(void) const
{
   return(new ConstrainedRealVector(*this));
}

inline BitVector::BitVector(void)
: Representation(0)
{
   vector = (unsigned char *)NULL;
   mytype = T_BitV;
}

//  Do a shallow copy
inline BitVector::BitVector(int num_els, unsigned char *const array)
: Representation(num_els)
{
   vector = array;
   mytype = T_BitV;
}

inline BitVector::~BitVector(void)
{
   if(vector!=(unsigned char *)NULL)
   {
      delete [] vector;
   }
}

inline Representation *BitVector::clone(void) const
{
   return(new BitVector(*this));
}

#endif
