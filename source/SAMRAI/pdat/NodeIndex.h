/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_NodeIndex
#define included_pdat_NodeIndex

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Utilities.h"

#include <vector>

namespace SAMRAI {
namespace pdat {

/**
 * Class NodeIndex implements a simple n-dimensional integer
 * vector for node centered variables.  Given a hier::Box in the AMR abstract
 * index space, the index space for a node-centered variable runs from the
 * lower corner of the box to the upper corner of the box plus one in each
 * direction.  See the node box geometry class for more information about
 * the mapping between the AMR index space and the node indices.
 *
 * @see hier::Index
 * @see NodeData
 * @see NodeGeometry
 * @see NodeIterator
 */

class NodeIndex:public hier::Index
{
public:
   /**
    * The Corner enumerated type is used when converting from a cell centered
    * index to a node centered index.  In 1d, use Left and Right.  In 2d, use
    * LowerLeft, LowerRight, UpperLeft, and UpperRight.  In 3d, the naming is
    * less intuitive, and use names LLL through UUU, where L means lower and
    * U means upper.  Therefore, to get the box upper in X, lower in Y, and
    * lower in Z, use corner name ULL.
    */
   enum Corner {
      Left = 0, Right = 1,
      LowerLeft = 0, LowerRight = 1, UpperLeft = 2, UpperRight = 3,
      LLL = 0, ULL = 1, LUL = 2, UUL = 3, LLU = 4, ULU = 5, LUU = 6, UUU = 7
   };

   /**
    * The default constructor for a node index creates an uninitialized index.
    */
   explicit NodeIndex(
      const tbox::Dimension& dim);

   /**
    * Construct a node index from a regular index and a corner.
    *
    * The Corner enumerated type is only defined for 3D or lower, so use
    * the next constructor with an hier::IntVector argument when using higher
    * dimensions.
    *
    * @pre getDim() == rhs.getDim()
    */
   NodeIndex(
      const hier::Index& rhs,
      const Corner corner);

   /**
    * Construct a node index from a regular index and an hier::IntVector.  For
    * each direction, if the hier::IntVector contains a 0, the node index will
    * represent a lower bound in that direction, and if 1 will represent an
    * upper bound in that direction.
    *
    * @pre getDim() == rhs.getDim()
    * @pre for each element of corner, e, e == 0 || e == 1
    */
   NodeIndex(
      const hier::Index& rhs,
      const hier::IntVector& corner);

   /**
    * The copy constructor creates a node index equal to the argument.
    *
    * @pre getDim() == rhs.getDim()
    */
   NodeIndex(
      const NodeIndex& rhs);

   /**
    * The assignment operator sets the node index equal to the argument.
    *
    * @pre getDim() == rhs.getDim()
    */
   NodeIndex&
   operator = (
      const NodeIndex& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator = (rhs);
      return *this;
   }

   /**
    * The node index destructor does nothing interesting.
    */
   ~NodeIndex();

   /**
    * Plus-equals operator for a node index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   NodeIndex&
   operator += (
      const hier::IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator += (rhs);
      return *this;
   }

   /**
    * Plus-equals operator for a node index and a node index.
    *
    * @pre getDim() == rhs.getDim()
    */
   NodeIndex&
   operator += (
      const NodeIndex& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator += (rhs);
      return *this;
   }

   /**
    * Plus operator for a node index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   NodeIndex
   operator + (
      const hier::IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      NodeIndex tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * Plus operator for a node index and a node index.
    *
    * @pre getDim() == rhs.getDim()
    */
   NodeIndex
   operator + (
      const NodeIndex& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      NodeIndex tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * Plus-equals operator for a node index and an integer.
    */
   NodeIndex&
   operator += (
      const int rhs)
   {
      hier::Index::operator += (rhs);
      return *this;
   }

   /**
    * Plus operator for a node index and an integer.
    */
   NodeIndex
   operator + (
      const int rhs) const
   {
      NodeIndex tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * Minus-equals operator for a node index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   NodeIndex&
   operator -= (
      const hier::IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator -= (rhs);
      return *this;
   }

   /**
    * Minus operator for a node index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   NodeIndex
   operator - (
      const hier::IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      NodeIndex tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * Minus-equals operator for a node index and a node index.
    *
    * @pre getDim() == rhs.getDim()
    */
   NodeIndex&
   operator -= (
      const NodeIndex& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator -= (rhs);
      return *this;
   }

   /**
    * Minus operator for a node index and a node index.
    *
    * @pre getDim() == rhs.getDim()
    */
   NodeIndex
   operator - (
      const NodeIndex& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      NodeIndex tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * Minus-equals operator for a node index and an integer.
    */
   NodeIndex&
   operator -= (
      const int rhs)
   {
      hier::Index::operator -= (rhs);
      return *this;
   }

   /**
    * Minus operator for a node index and an integer.
    */
   NodeIndex
   operator - (
      const int rhs) const
   {
      NodeIndex tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * Times-equals operator for a node index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   NodeIndex&
   operator *= (
      const hier::IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator *= (rhs);
      return *this;
   }

   /**
    * Times operator for a node index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   NodeIndex
   operator * (
      const hier::IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      NodeIndex tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * Times-equals operator for a node index and an integer.
    */
   NodeIndex&
   operator *= (
      const int rhs)
   {
      hier::Index::operator *= (rhs);
      return *this;
   }

   /**
    * Times operator for a node index and an integer.
    */
   NodeIndex
   operator * (
      const int rhs) const
   {
      NodeIndex tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * Returns true if two node index objects are equal.
    * All components must be the same for equality.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator == (
      const NodeIndex& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      return ((hier::Index *)this)->operator == (rhs);
   }

   /**
    * Returns true if two node index objects are not equal.
    * Any of the components may be different for inequality.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator != (
      const NodeIndex& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      return ((hier::Index *)this)->operator != (rhs);
   }

private:
   /*
    * Initializes the offsets if it has not yet been done
    */
   void
   setOffsets();

   static std::vector<hier::IntVector> s_offsets[SAMRAI::MAX_DIM_VAL];
   static bool s_offsets_are_set[SAMRAI::MAX_DIM_VAL];
};

}
}

#endif
