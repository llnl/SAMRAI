/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_FaceIndex
#define included_pdat_FaceIndex

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

/**
 * Class FaceIndex implements a simple n-dimensional integer
 * vector for face centered variables.  Face indices contain an integer
 * index location in AMR index space along with the designated face axis
 * (X=0, Y=1, or Z=2).  See the face box geometry class for more information
 * about the mapping between the AMR index space and the face indices.
 *
 * @see hier::Index
 * @see FaceData
 * @see FaceGeometry
 * @see FaceIterator
 */

class FaceIndex:public hier::Index
{
public:
   /**
    * The default constructor for a face index creates an uninitialized index.
    */
   explicit FaceIndex(
      const tbox::Dimension& dim);

   /**
    * Construct a face index from a regular index, axis, and face.  The axis
    * can be one of FaceIndex::X (0), FaceIndex::Y (1), or
    * FaceIndex::Z (2). The face argument can be one of the constants
    * FaceIndex::Lower (0) or FaceIndex::Upper (1).
    */
   FaceIndex(
      const hier::Index& rhs,
      const int axis,
      const int face);

   /**
    * The copy constructor creates a face index equal to the argument.
    */
   FaceIndex(
      const FaceIndex& rhs);

   /**
    * The assignment operator sets the face index equal to the argument.
    *
    * @pre getDim() == rhs.getDim()
    */
   FaceIndex&
   operator = (
      const FaceIndex& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator = (rhs);
      d_axis = rhs.d_axis;
      return *this;
   }

   /**
    * The face index destructor does nothing interesting.
    */
   ~FaceIndex();

   /**
    * Get the axis for which this face index is defined (X=0, Y=1, Z=2).
    */
   int
   getAxis() const
   {
      return d_axis;
   }

   /**
    * Set the face axis (X=0, Y=1, Z=2).
    */
   void
   setAxis(
      const int axis)
   {
      d_axis = axis;
   }

   /**
    * Convert the face index into the index on the left hand face
    * (argument face == 0) or the right hand face (argument face == 1).
    */
   hier::Index
   toCell(
      const int face) const;

   /**
    * Plus-equals operator for a face index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   FaceIndex&
   operator += (
      const hier::IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator += (rhs);
      return *this;
   }

   /**
    * Plus-equals operator for a face index and a face index.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getAxis() == rhs.getAxis()
    */
   FaceIndex&
   operator += (
      const FaceIndex& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_axis == rhs.d_axis);
      hier::Index::operator += (rhs);
      return *this;
   }

   /**
    * Plus operator for a face index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   FaceIndex
   operator + (
      const hier::IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      FaceIndex tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * Plus operator for a face index and a face index.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getAxis() == rhs.getAxis()
    */
   FaceIndex
   operator + (
      const FaceIndex& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_axis == rhs.d_axis);
      FaceIndex tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * Plus-equals operator for a face index and an integer.
    */
   FaceIndex&
   operator += (
      const int rhs)
   {
      hier::Index::operator += (rhs);
      return *this;
   }

   /**
    * Plus operator for a face index and an integer.
    */
   FaceIndex
   operator + (
      const int rhs) const
   {
      FaceIndex tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * Minus-equals operator for a face index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   FaceIndex&
   operator -= (
      const hier::IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator -= (rhs);
      return *this;
   }

   /**
    * Minus operator for a face index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   FaceIndex
   operator - (
      const hier::IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      FaceIndex tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * Minus-equals operator for a face index and a face index.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getAxis() == rhs.getAxis()
    */
   FaceIndex&
   operator -= (
      const FaceIndex& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_axis == rhs.d_axis);
      hier::Index::operator -= (rhs);
      return *this;
   }

   /**
    * Minus operator for a face index and a face index.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getAxis() == rhs.getAxis()
    */
   FaceIndex
   operator - (
      const FaceIndex& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_axis == rhs.d_axis);
      FaceIndex tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * Minus-equals operator for a face index and an integer.
    */
   FaceIndex&
   operator -= (
      const int rhs)
   {
      hier::Index::operator -= (rhs);
      return *this;
   }

   /**
    * Minus operator for a face index and an integer.
    */
   FaceIndex
   operator - (
      const int rhs) const
   {
      FaceIndex tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * Times-equals operator for a face index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   FaceIndex&
   operator *= (
      const hier::IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator *= (rhs);
      return *this;
   }

   /**
    * Times operator for a face index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   FaceIndex
   operator * (
      const hier::IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      FaceIndex tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * Times-equals operator for a face index and an integer.
    */
   FaceIndex&
   operator *= (
      const int rhs)
   {
      hier::Index::operator *= (rhs);
      return *this;
   }

   /**
    * Times operator for a face index and an integer.
    */
   FaceIndex
   operator * (
      const int rhs) const
   {
      FaceIndex tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * Returns true if two face index objects are equal.  All components
    * and the corresponding face axes must be the same for equality.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator == (
      const FaceIndex& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      return ((hier::Index *)this)->operator == (rhs) && (d_axis == rhs.d_axis);
   }

   /**
    * Returns true if two face index objects are not equal.  Any of
    * the components or axes may be different for inequality.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator != (
      const FaceIndex& rhs) const
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
