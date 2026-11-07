/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex edge-centered patch data.
 *
 ************************************************************************/
#include "SAMRAI/math/PatchEdgeDataOpsComplex.h"
#include "SAMRAI/pdat/EdgeGeometry.h"

namespace SAMRAI {
namespace math {

PatchEdgeDataOpsComplex::PatchEdgeDataOpsComplex()
{
}

PatchEdgeDataOpsComplex::~PatchEdgeDataOpsComplex()
{
}

/*
 *************************************************************************
 *
 * General operations for complex edge-centered patch data.
 *
 *************************************************************************
 */

void
PatchEdgeDataOpsComplex::swapData(
   const std::shared_ptr<hier::Patch>& patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   std::shared_ptr<pdat::EdgeData<dcomplex> > d1(
      SAMRAI_SHARED_PTR_CAST<pdat::EdgeData<dcomplex>, hier::PatchData>(
         patch->getPatchData(data1_id)));
   std::shared_ptr<pdat::EdgeData<dcomplex> > d2(
      SAMRAI_SHARED_PTR_CAST<pdat::EdgeData<dcomplex>, hier::PatchData>(
         patch->getPatchData(data2_id)));

   TBOX_ASSERT(d1 && d2);
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));

   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

void
PatchEdgeDataOpsComplex::printData(
   const std::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

void
PatchEdgeDataOpsComplex::copyData(
   const std::shared_ptr<pdat::EdgeData<dcomplex> >& dst,
   const std::shared_ptr<pdat::EdgeData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   int dimVal = box.getDim().getValue();
   for (int d = 0; d < dimVal; ++d) {
      dst->getArrayData(d).copy(src->getArrayData(d),
         pdat::EdgeGeometry::toEdgeBox(box, d));
   }
}

}
}
