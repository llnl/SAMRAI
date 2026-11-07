/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer face-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchFaceDataOpsInteger
#define included_math_PatchFaceDataOpsInteger

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/math/PatchFaceDataBasicOps.h"
#include "SAMRAI/math/ArrayDataNormOpsInteger.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iostream>
#include <memory>

namespace SAMRAI {
namespace math {

/**
 * Class PatchFaceDataOpsInteger provides a collection of operations
 * that may be used to manipulate integer face-centered patch data.  The
 * operations include basic arithmetic, min, max, etc.  With the assertion
 * of a few basic routines, this class inherits its interface (and
 * thus its functionality) from the base class PatchFaceDataBasicOps
 * from which it is derived.
 *
 * A more extensive set of operations is implemented for real (double and
 * float) and complex patch data in the classes PatchFaceDataOpsReal
 * and PatchFaceDataOpsComplex, repsectively.
 *
 * @see PatchFaceDataBasicOps
 */

class PatchFaceDataOpsInteger:
   public PatchFaceDataBasicOps<int>
{
public:
   /**
    * Empty constructor and destructor.
    */
   PatchFaceDataOpsInteger();

   virtual ~PatchFaceDataOpsInteger();

   /**
    * Return the number of data values for the face-centered data object
    * in the given box.  Note that it is assumed that the box refers to
    * the cell-centered index space corresponding to the patch hierarchy.
    *
    * @pre data
    * @pre data->getDim() == box.getDim()
    */
   int
   numberOfEntries(
      const std::shared_ptr<pdat::FaceData<int> >& data,
      const hier::Box& box) const;

   /**
    * Copy dst data to src data over given box.
    *
    * @pre dst && src
    * @pre (dst->getDim() == src->getDim()) && (dst->getDim() == box.getDim())
    */
   void
   copyData(
      const std::shared_ptr<pdat::FaceData<int> >& dst,
      const std::shared_ptr<pdat::FaceData<int> >& src,
      const hier::Box& box) const;

   /**
    * Swap pointers for patch data objects.  Objects are checked for
    * consistency of depth, box, and ghost box.
    *
    * @pre patch
    * @pre patch->getPatchData(data1_id) is actually a std::shared_ptr<pdat::FaceData<int> >
    * @pre patch->getPatchData(data2_id) is actually a std::shared_ptr<pdat::FaceData<int> >
    * @pre patch->getPatchData(data1_id)->getDepth() == patch->getPatchData(data2_id)->getDepth()
    * @pre patch->getPatchData(data1_id)->getBox().isSpatiallyEqual(patch->getPatchData(data2_id)->getBox())
    * @pre patch->getPatchData(data1_id)->getGhostBox().isSpatiallyEqual(patch->getPatchData(data2_id)->getGhostBox())
    */
   void
   swapData(
      const std::shared_ptr<hier::Patch>& patch,
      const int data1_id,
      const int data2_id) const;

   /**
    * Print data entries over given box to given output stream.
    *
    * @pre data
    * @pre data->getDim() == box.getDim()
    */
   void
   printData(
      const std::shared_ptr<pdat::FaceData<int> >& data,
      const hier::Box& box,
      std::ostream& s = tbox::plog) const;

   /**
    * Initialize data to given scalar over given box.
    *
    * @pre dst
    * @pre dst->getDim() == box.getDim()
    */
   void
   setToScalar(
      const std::shared_ptr<pdat::FaceData<int> >& dst,
      const int& alpha,
      const hier::Box& box) const
   {
      TBOX_ASSERT(dst);
      TBOX_ASSERT_OBJDIM_EQUALITY2(*dst, box);
      dst->fillAll(alpha, box);
   }

   /**
    * Set destination component to absolute value of source component.
    * That is, each destination entry is set to \f$d_i = \| s_i \|\f$.
    *
    * @pre dst && src
    * @pre (dst->getDim() == src->getDim()) && (dst->getDim() == box.getDim())
    */
   void
   abs(
      const std::shared_ptr<pdat::FaceData<int> >& dst,
      const std::shared_ptr<pdat::FaceData<int> >& src,
      const hier::Box& box) const;

private:
   // The following are not implemented:
   PatchFaceDataOpsInteger(
      const PatchFaceDataOpsInteger&);
   PatchFaceDataOpsInteger&
   operator = (
      const PatchFaceDataOpsInteger&);

   ArrayDataNormOpsInteger d_array_ops;

};

}
}

#endif
