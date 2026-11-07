/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_SideIndex
#define included_pdat_SideIndex

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

/**
 * Class SideIndex implements a simple n-dimensional integer
 * vector for side centered variables.  Side indices contain an integer
 * index location in AMR index space along with the designated side axis
 * (X=0, Y=1, or Z=2).  See the side box geometry class for more information
 * about the mapping between the AMR index space and the side indices.
 *
 * @see hier::Index
 * @see SideData
 * @see SideGeometry
 * @see SideIterator
 */

class SideIndex:public hier::Index
{
public:
   /**
    * The default constructor for a side index creates an uninitialized index.
    */
   explicit SideIndex(
      const tbox::Dimension& dim);

   /**
    * Construct a side index from a regular index, axis, and side.  The axis
    * can be one of SideIndex::X (0), SideIndex::Y (1), or
    * SideIndex::Z (2). The side argument can be one of the constants
    * SideIndex::Lower (0) or SideIndex::Upper (1).
    */
   SideIndex(
      const hier::Index& rhs,
      const int axis,
      const int side);

   /**
    * The copy constructor creates a side index equal to the argument.
    */
   SideIndex(
      const SideIndex& rhs);

   /**
    * The assignment operator sets the side index equal to the argument.
    *
    * @pre getDim() == rhs.getDim()
    */
   SideIndex&
   operator = (
      const SideIndex& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator = (rhs);
      d_axis = rhs.d_axis;
      return *this;
   }

   /**
    * The side index destructor does nothing interesting.
    */
   ~SideIndex();

   /**
    * Get the axis for which this side index is defined (X=0, Y=1, Z=2).
    */
   int
   getAxis() const
   {
      return d_axis;
   }

   /**
    * Set the side axis (X=0, Y=1, Z=2).
    */
   void
   setAxis(
      const int axis)
   {
      d_axis = axis;
   }

   /**
    * Convert the side index into the index on the left hand side
    * (argument side == 0) or the right hand side (argument side == 1).
    */
   hier::Index
   toCell(
      const int side) const;

   /**
    * Plus-equals operator for a side index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   SideIndex&
   operator += (
      const hier::IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator += (rhs);
      return *this;
   }

   /**
    * Plus-equals operator for a side index and a side index.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getAxis() == rhs.getAxis()
    */
   SideIndex&
   operator += (
      const SideIndex& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_axis == rhs.d_axis);
      hier::Index::operator += (rhs);
      return *this;
   }

   /**
    * Plus operator for a side index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   SideIndex
   operator + (
      const hier::IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      SideIndex tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * Plus operator for a side index and a side index.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getAxis() == rhs.getAxis()
    */
   SideIndex
   operator + (
      const SideIndex& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_axis == rhs.d_axis);
      SideIndex tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * Plus-equals operator for a side index and an integer.
    */
   SideIndex&
   operator += (
      const int rhs)
   {
      hier::Index::operator += (rhs);
      return *this;
   }

   /**
    * Plus operator for a side index and an integer.
    */
   SideIndex
   operator + (
      const int rhs) const
   {
      SideIndex tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * Minus-equals operator for a side index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   SideIndex&
   operator -= (
      const hier::IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator -= (rhs);
      return *this;
   }

   /**
    * Minus operator for a side index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   SideIndex
   operator - (
      const hier::IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      SideIndex tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * Minus-equals operator for a side index and a side index.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getAxis() == rhs.getAxis()
    */
   SideIndex&
   operator -= (
      const SideIndex& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_axis == rhs.d_axis);
      hier::Index::operator -= (rhs);
      return *this;
   }

   /**
    * Minus operator for a side index and a side index.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getAxis() == rhs.getAxis()
    */
   SideIndex
   operator - (
      const SideIndex& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_axis == rhs.d_axis);
      SideIndex tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * Minus-equals operator for a side index and an integer.
    */
   SideIndex&
   operator -= (
      const int rhs)
   {
      hier::Index::operator -= (rhs);
      return *this;
   }

   /**
    * Minus operator for a side index and an integer.
    */
   SideIndex
   operator - (
      const int rhs) const
   {
      SideIndex tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * Times-equals operator for a side index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   SideIndex&
   operator *= (
      const hier::IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator *= (rhs);
      return *this;
   }

   /**
    * Times operator for a side index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   SideIndex
   operator * (
      const hier::IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      SideIndex tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * Times-equals operator for a side index and an integer.
    */
   SideIndex&
   operator *= (
      const int rhs)
   {
      hier::Index::operator *= (rhs);
      return *this;
   }

   /**
    * Times operator for a side index and an integer.
    */
   SideIndex
   operator * (
      const int rhs) const
   {
      SideIndex tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * Returns true if two side index objects are equal.  All components
    * and the corresponding side axes must be the same for equality.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator == (
      const SideIndex& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      return ((hier::Index *)this)->operator == (rhs) && (d_axis == rhs.d_axis);
   }

   /**
    * Returns true if two side index objects are not equal.  Any of
    * the components or axes may be different for inequality.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator != (
      const SideIndex& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      return ((hier::Index *)this)->operator != (rhs) || (d_axis != rhs.d_axis);
   }

   enum {
      X = 0,
      Y = 1,
      Z = 2,
      Lower = 0,
      Upper = 1
   };

private:
   int d_axis;
};

}
}

#endif
