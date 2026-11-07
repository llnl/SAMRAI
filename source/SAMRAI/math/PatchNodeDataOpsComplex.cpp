/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex node-centered patch data.
 *
 ************************************************************************/
#include "SAMRAI/math/PatchNodeDataOpsComplex.h"
#include "SAMRAI/pdat/NodeGeometry.h"

namespace SAMRAI {
namespace math {

PatchNodeDataOpsComplex::PatchNodeDataOpsComplex()
{
}

PatchNodeDataOpsComplex::~PatchNodeDataOpsComplex()
{
}

/*
 *************************************************************************
 *
 * General operations for complex node-centered patch data.
 *
 *************************************************************************
 */

void
PatchNodeDataOpsComplex::swapData(
   const std::shared_ptr<hier::Patch>& patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   std::shared_ptr<pdat::NodeData<dcomplex> > d1(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<dcomplex>, hier::PatchData>(
         patch->getPatchData(data1_id)));
   std::shared_ptr<pdat::NodeData<dcomplex> > d2(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<dcomplex>, hier::PatchData>(
         patch->getPatchData(data2_id)));

   TBOX_ASSERT(d1 && d2);
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));

   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

void
PatchNodeDataOpsComplex::printData(
   const std::shared_ptr<pdat::NodeData<dcomplex> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

}
}
