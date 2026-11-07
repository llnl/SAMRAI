/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer cell-centered patch data.
 *
 ************************************************************************/
#include "SAMRAI/math/PatchCellDataOpsInteger.h"

namespace SAMRAI {
namespace math {

PatchCellDataOpsInteger::PatchCellDataOpsInteger()
{
}

PatchCellDataOpsInteger::~PatchCellDataOpsInteger()
{
}

/*
 *************************************************************************
 *
 * General operations for integer cell-centered patch data.
 *
 *************************************************************************
 */

void
PatchCellDataOpsInteger::swapData(
   const std::shared_ptr<hier::Patch>& patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   std::shared_ptr<pdat::CellData<int> > d1(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
         patch->getPatchData(data1_id)));
   std::shared_ptr<pdat::CellData<int> > d2(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
         patch->getPatchData(data2_id)));

   TBOX_ASSERT(d1 && d2);
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));

   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

}
}
