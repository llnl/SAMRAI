/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer side-centered patch data.
 *
 ************************************************************************/
#include "SAMRAI/math/PatchSideDataOpsInteger.h"
#include "SAMRAI/pdat/SideGeometry.h"

namespace SAMRAI {
namespace math {

PatchSideDataOpsInteger::PatchSideDataOpsInteger()
{
}

PatchSideDataOpsInteger::~PatchSideDataOpsInteger()
{
}

/*
 *************************************************************************
 *
 * Compute the number of data entries on a patch in the given box.
 *
 *************************************************************************
 */

int
PatchSideDataOpsInteger::numberOfEntries(
   const std::shared_ptr<pdat::SideData<int> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   int dimVal = box.getDim().getValue();
   int retval = 0;
   const hier::Box ibox = box * data->getGhostBox();
   const hier::IntVector& directions = data->getDirectionVector();
   const int data_depth = data->getDepth();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         retval +=
            static_cast<int>((pdat::SideGeometry::toSideBox(ibox, d).size()) * data_depth);
      }
   }
   return retval;
}

/*
 *************************************************************************
 *
 * General operations for integer side-centered patch data.
 *
 *************************************************************************
 */

void PatchSideDataOpsInteger::swapData(
   const std::shared_ptr<hier::Patch>& patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   std::shared_ptr<pdat::SideData<int> > d1(
      SAMRAI_SHARED_PTR_CAST<pdat::SideData<int>, hier::PatchData>(
         patch->getPatchData(data1_id)));
   std::shared_ptr<pdat::SideData<int> > d2(
      SAMRAI_SHARED_PTR_CAST<pdat::SideData<int>, hier::PatchData>(
         patch->getPatchData(data2_id)));

   TBOX_ASSERT(d1 && d2);
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getDirectionVector() == d2->getDirectionVector());
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));

   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

void
PatchSideDataOpsInteger::printData(
   const std::shared_ptr<pdat::SideData<int> >& data,
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
PatchSideDataOpsInteger::copyData(
   const std::shared_ptr<pdat::SideData<int> >& dst,
   const std::shared_ptr<pdat::SideData<int> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();
   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         dst->getArrayData(d).copy(src->getArrayData(d),
            pdat::SideGeometry::toSideBox(box, d));
      }
   }
}
void
PatchSideDataOpsInteger::abs(
   const std::shared_ptr<pdat::SideData<int> >& dst,
   const std::shared_ptr<pdat::SideData<int> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   int dimVal = box.getDim().getValue();
   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         d_array_ops.abs(dst->getArrayData(d),
            src->getArrayData(d),
            pdat::SideGeometry::toSideBox(box, d));
      }
   }
}

}
}
