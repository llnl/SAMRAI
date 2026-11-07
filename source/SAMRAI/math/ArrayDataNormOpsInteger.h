/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated norm operations for real data arrays.
 *
 ************************************************************************/

#ifndef included_math_ArrayDataNormOpsInteger
#define included_math_ArrayDataNormOpsInteger

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/ArrayData.h"

namespace SAMRAI {
namespace math {

/**
 * Class ArrayDataNormOpsInteger provides a set of common norm
 * operations that may be applied to arrays of integer data values
 * maintained as pdat::ArrayData<TYPE> objects.  The intent of this class
 * is to provide a single implementation of these operations as they are needed
 * by objects that perform these operations on the standard array-based patch
 * data types (i.e., cell-centered, face-centered, node-centered).
 * Note that each operation is performed on the intersection of the box in
 * the function argument list and the boxes associated with all
 * pdat::ArrayData<TYPE> objects.  Currently, the only norm operation implemented
 * in this class is the absolute value operation.
 *
 * @see pdat::ArrayData
 */

class ArrayDataNormOpsInteger
{
public:
   /**
    * Empty constructor and destructor.
    */
   ArrayDataNormOpsInteger();

   ~ArrayDataNormOpsInteger();

   /**
    * Set destination component to absolute value of source component.
    * That is, each destination entry is set to \f$d_i = \| s_i \|\f$.
    *
    * @pre (dst.getDim() == src.getDim()) && (dst.getDim() == box.getDim())
    * @pre dst.getDepth() == src.getDepth()
    */
   void
   abs(
      pdat::ArrayData<int>& dst,
      const pdat::ArrayData<int>& src,
      const hier::Box& box) const;

private:
   // The following are not implemented:
   ArrayDataNormOpsInteger(
      const ArrayDataNormOpsInteger&);
   ArrayDataNormOpsInteger&
   operator = (
      const ArrayDataNormOpsInteger&);
};

}
}
#endif
