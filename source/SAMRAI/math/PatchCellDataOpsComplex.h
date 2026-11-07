/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex cell-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchCellDataOpsComplex
#define included_math_PatchCellDataOpsComplex

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/PatchCellDataBasicOps.h"
#include "SAMRAI/math/PatchCellDataNormOpsComplex.h"
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
 * Class PatchCellDataOpsComplex provides a collection of operations
 * that may be used to manipulate complex cell-centered patch data.  The
 * operations include basic arithmetic and norms.  With the
 * assertion of a few basic routines, this class inherits its interface (and
 * thus its functionality) from the base classes PatchCellDataBasicOps,
 * PatchCellDataNormOpsComplex from which it is derived.  The
 * name of each of these base classes is indicative of the set of
 * cell-centered patch data operations that it provides.
 *
 * A similar set of operations is implemented for real (double and float) and
 * integer patch data in the classes PatchCellDataOpsReal and
 * PatchCellDataOpsInteger, respectively.
 *
 * @see PatchCellDataBasicOps
 * @see PatchCellDataNormOpsComplex
 */

class PatchCellDataOpsComplex:
   public PatchCellDataBasicOps<dcomplex>,
   public PatchCellDataNormOpsComplex
{
public:
   /**
    * Empty constructor and destructor.
    */
   PatchCellDataOpsComplex();

   virtual ~PatchCellDataOpsComplex();

   /**
    * Copy dst data to src data over given box.
    *
    * @pre dst && src
    * @pre (dst->getDim() == src->getDim()) && (dst->getDim() == box.getDim())
    */
   void
   copyData(
      const std::shared_ptr<pdat::CellData<dcomplex> >& dst,
      const std::shared_ptr<pdat::CellData<dcomplex> >& src,
      const hier::Box& box) const
   {
      TBOX_ASSERT(dst && src);
      TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);
      dst->getArrayData().copy(src->getArrayData(), box);
   }

   /**
    * Swap pointers for patch data objects.  Objects are checked for
    * consistency of depth, box, and ghost box.
    *
    * @pre patch->getPatchData(data1_id) is actually a std::shared_ptr<pdat::CellData<dcomplex> >
    * @pre patch->getPatchData(data2_id) is actually a std::shared_ptr<pdat::CellData<dcomplex> >
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
      const std::shared_ptr<pdat::CellData<dcomplex> >& data,
      const hier::Box& box,
      std::ostream& s = tbox::plog) const
   {
      TBOX_ASSERT(data);
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);
      s << "Data box = " << box << std::endl;
      data->print(box, s);
      s << "\n";
   }

   /**
    * Initialize data to given scalar over given box.
    *
    * @pre dst
    * @pre dst->getDim() == box.getDim()
    */
   void
   setToScalar(
      const std::shared_ptr<pdat::CellData<dcomplex> >& dst,
      const dcomplex& alpha,
      const hier::Box& box) const
   {
      TBOX_ASSERT(dst);
      TBOX_ASSERT_OBJDIM_EQUALITY2(*dst, box);
      dst->fillAll(alpha, box);
   }

private:
   // The following are not implemented:
   PatchCellDataOpsComplex(
      const PatchCellDataOpsComplex&);
   PatchCellDataOpsComplex&
   operator = (
      const PatchCellDataOpsComplex&);

};

}
}

#endif
