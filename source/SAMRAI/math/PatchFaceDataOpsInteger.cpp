/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer face-centered patch data.
 *
 ************************************************************************/
#include "SAMRAI/math/PatchFaceDataOpsInteger.h"
#include "SAMRAI/pdat/FaceGeometry.h"

namespace SAMRAI {
namespace math {

PatchFaceDataOpsInteger::PatchFaceDataOpsInteger()
{
}

PatchFaceDataOpsInteger::~PatchFaceDataOpsInteger()
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
PatchFaceDataOpsInteger::numberOfEntries(
   const std::shared_ptr<pdat::FaceData<int> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   int dimVal = box.getDim().getValue();
   int retval = 0;
   const hier::Box ibox = box * data->getGhostBox();
   const int data_depth = data->getDepth();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      retval +=
         static_cast<tbox::Dimension::dir_t>((pdat::FaceGeometry::toFaceBox(ibox,
                                                 d).size()) * data_depth);
   }
   return retval;
}

/*
 *************************************************************************
 *
 * General operations for integer face-centered patch data.
 *
 *************************************************************************
 */

void
PatchFaceDataOpsInteger::swapData(
   const std::shared_ptr<hier::Patch>& patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   std::shared_ptr<pdat::FaceData<int> > d1(
      SAMRAI_SHARED_PTR_CAST<pdat::FaceData<int>, hier::PatchData>(
         patch->getPatchData(data1_id)));
   std::shared_ptr<pdat::FaceData<int> > d2(
      SAMRAI_SHARED_PTR_CAST<pdat::FaceData<int>, hier::PatchData>(
         patch->getPatchData(data2_id)));

   TBOX_ASSERT(d1 && d2);
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));

   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

void
PatchFaceDataOpsInteger::printData(
   const std::shared_ptr<pdat::FaceData<int> >& data,
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
PatchFaceDataOpsInteger::copyData(
   const std::shared_ptr<pdat::FaceData<int> >& dst,
   const std::shared_ptr<pdat::FaceData<int> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   int dimVal = box.getDim().getValue();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      dst->getArrayData(d).copy(src->getArrayData(d),
         pdat::FaceGeometry::toFaceBox(box, d));
   }
}

void
PatchFaceDataOpsInteger::abs(
   const std::shared_ptr<pdat::FaceData<int> >& dst,
   const std::shared_ptr<pdat::FaceData<int> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   int dimVal = box.getDim().getValue();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      d_array_ops.abs(dst->getArrayData(d),
         src->getArrayData(d),
         pdat::FaceGeometry::toFaceBox(box, d));
   }
}

}
}
