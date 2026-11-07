/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex side-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchSideDataOpsComplex
#define included_math_PatchSideDataOpsComplex

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/math/PatchSideDataBasicOps.h"
#include "SAMRAI/math/PatchSideDataNormOpsComplex.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iostream>
#include <memory>

namespace SAMRAI {
namespace math {

/**
 * Class PatchSideDataOpsComplex provides a collection of operations
 * that may be used to manipulate complex side-centered patch data.  The
 * operations include basic arithmetic and norms.  With the
 * assertion of a few basic routines, this class inherits its interface (and
 * thus its functionality) from the base classes PatchSideDataBasicOps,
 * PatchSideDataNormOpsComplex from which it is derived.  The
 * name of each of these base classes is indicative of the set of
 * side-centered patch data operations that it provides.
 *
 * A similar set of operations is implemented for real (double and float) and
 * integer patch data in the classes PatchSideDataOpsReal and
 * PatchSideDataOpsInteger, repsectively.
 *
 * @see PatchSideDataBasicOps
 * @see PatchSideDataNormOpsComplex
 */

class PatchSideDataOpsComplex:
   public PatchSideDataBasicOps<dcomplex>,
   public PatchSideDataNormOpsComplex
{
public:
   /**
    * Empty constructor and destructor.
    */
   PatchSideDataOpsComplex();

   virtual ~PatchSideDataOpsComplex();

   /**
    * Copy dst data to src data over given box.
    *
    * @pre dst && src
    * @pre dst->getDirectionVector() == src->getDirectionVector()
    * @pre (dst->getDim() == src->getDim()) && (dst->getDim() == box.getDim())
    */
   void
   copyData(
      const std::shared_ptr<pdat::SideData<dcomplex> >& dst,
      const std::shared_ptr<pdat::SideData<dcomplex> >& src,
      const hier::Box& box) const;

   /**
    * Swap pointers for patch data objects.  Objects are checked for
    * consistency of depth, box, and ghost box.
    *
    * @pre patch
    * @pre patch->getPatchData(data1_id) is actually a std::shared_ptr<pdat::SideData<dcomplex> >
    * @pre patch->getPatchData(data2_id) is actually a std::shared_ptr<pdat::SideData<dcomplex> >
    * @pre patch->getPatchData(data1_id)->getDepth() == patch->getPatchData(data2_id)->getDepth()
    * @pre patch->getPatchData(data1_id)->getBox().isSpatiallyEqual(patch->getPatchData(data2_id)->getBox())
    * @pre patch->getPatchData(data1_id)->getDirectionVector() == patch->getPatchData(data2_id)->getDirectionVector()
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
      const std::shared_ptr<pdat::SideData<dcomplex> >& data,
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
      const std::shared_ptr<pdat::SideData<dcomplex> >& dst,
      const dcomplex& alpha,
      const hier::Box& box) const
   {
      TBOX_ASSERT(dst);
      TBOX_ASSERT_OBJDIM_EQUALITY2(*dst, box);
      dst->fillAll(alpha, box);
   }

private:
   // The following are not implemented:
   PatchSideDataOpsComplex(
      const PatchSideDataOpsComplex&);
   PatchSideDataOpsComplex&
   operator = (
      const PatchSideDataOpsComplex&);

};

}
}

#endif
