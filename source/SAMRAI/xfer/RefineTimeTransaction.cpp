/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Communication transaction for time interpolation during data refining
 *
 ************************************************************************/
#include "SAMRAI/xfer/RefineTimeTransaction.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <typeinfo>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Initialization, set/unset functions for static array of refine items
 * and interpolation time.
 *
 *************************************************************************
 */

double RefineTimeTransaction::s_time = 0.0;

/*
 *************************************************************************
 *
 * Constructor sets state of transaction.
 *
 *************************************************************************
 */
RefineTimeTransaction::RefineTimeTransaction(
   const std::shared_ptr<hier::PatchLevel>& dst_level,
   const std::shared_ptr<hier::PatchLevel>& src_level,
   const std::shared_ptr<hier::BoxOverlap>& overlap,
   const hier::Box& dst_box,
   const hier::Box& src_box,
   const hier::Box& box,
   const RefineClasses::Data** refine_data,
   int item_id):
   d_dst_patch(),
   d_dst_patch_rank(dst_box.getOwnerRank()),
   d_src_patch(),
   d_src_patch_rank(src_box.getOwnerRank()),
   d_overlap(overlap),
   d_box(box),
   d_refine_data(refine_data),
   d_item_id(item_id)
{
   TBOX_ASSERT(dst_level);
   TBOX_ASSERT(src_level);
   TBOX_ASSERT(overlap);
   TBOX_ASSERT(dst_box.getLocalId() >= 0);
   TBOX_ASSERT(src_box.getLocalId() >= 0);
   TBOX_ASSERT(item_id >= 0);
   TBOX_ASSERT(refine_data[item_id] != 0);
   TBOX_ASSERT_OBJDIM_EQUALITY5(*dst_level,
      *src_level,
      dst_box,
      src_box,
      box);

   if (d_dst_patch_rank == dst_level->getBoxLevel()->getMPI().getRank()) {
      d_dst_patch = dst_level->getPatch(dst_box.getGlobalId());
   }
   if (d_src_patch_rank == dst_level->getBoxLevel()->getMPI().getRank()) {
      d_src_patch = src_level->getPatch(src_box.getGlobalId());
   }
}

RefineTimeTransaction::~RefineTimeTransaction()
{
}

/*
 *************************************************************************
 *
 * Functions overridden in tbox::Transaction base class.
 *
 *************************************************************************
 */

bool
RefineTimeTransaction::canEstimateIncomingMessageSize()
{
   bool can_estimate = false;
   if (d_src_patch) {
      can_estimate =
         d_src_patch->getPatchData(d_refine_data[d_item_id]->d_src_told)
         ->canEstimateStreamSizeFromBox();
   } else {
      can_estimate =
         d_dst_patch->getPatchData(d_refine_data[d_item_id]->d_scratch)
         ->canEstimateStreamSizeFromBox();
   }
   return can_estimate;
}

size_t
RefineTimeTransaction::computeIncomingMessageSize()
{
   d_incoming_bytes =
      d_dst_patch->getPatchData(d_refine_data[d_item_id]->d_scratch)
      ->getDataStreamSize(*d_overlap);
   return d_incoming_bytes;
}

size_t
RefineTimeTransaction::computeOutgoingMessageSize()
{
   d_outgoing_bytes =
      d_src_patch->getPatchData(d_refine_data[d_item_id]->d_src_told)
      ->getDataStreamSize(*d_overlap);
   return d_outgoing_bytes;
}

int
RefineTimeTransaction::getSourceProcessor() {
   return d_src_patch_rank;
}

int
RefineTimeTransaction::getDestinationProcessor() {
   return d_dst_patch_rank;
}

void
RefineTimeTransaction::packStream(
   tbox::MessageStream& stream)
{
   const hier::PatchData& src_told_data =
      *(d_src_patch->getPatchData(d_refine_data[d_item_id]->d_src_told));

   const double& told = src_told_data.getTime();
   if (tbox::MathUtilities<double>::equalEps(s_time, told)) {
      src_told_data.packStream(stream, *d_overlap);
   } else {

      hier::Box temporary_box(d_box.getDim());
      temporary_box.initialize(d_box,
                               d_src_patch->getBox().getLocalId(),
                               tbox::SAMRAI_MPI::getInvalidRank());

      hier::Patch temporary_patch(
         temporary_box,
         d_src_patch->getPatchDescriptor());

      std::shared_ptr<hier::PatchData> temporary_patch_data(
         d_src_patch->getPatchDescriptor()
         ->getPatchDataFactory(d_refine_data[d_item_id]->d_src_told)
         ->allocate(temporary_patch));
      temporary_patch_data->setTime(s_time);

      timeInterpolate(
         *temporary_patch_data,
         *d_overlap,
         src_told_data,
         d_src_patch->getPatchData(d_refine_data[d_item_id]->d_src_tnew));

      temporary_patch_data->packStream(stream, *d_overlap);
   }
}

void
RefineTimeTransaction::unpackStream(
   tbox::MessageStream& stream)
{
   d_dst_patch->getPatchData(d_refine_data[d_item_id]->d_scratch)
   ->unpackStream(stream, *d_overlap);
}

void
RefineTimeTransaction::copyLocalData()
{
   hier::PatchData& scratch_data =
      *(d_dst_patch->getPatchData(d_refine_data[d_item_id]->d_scratch));
   const hier::PatchData& src_told_data =
      *(d_src_patch->getPatchData(d_refine_data[d_item_id]->d_src_told));
   const std::shared_ptr<hier::PatchData>& src_tnew_data =
      d_src_patch->getPatchData(d_refine_data[d_item_id]->d_src_tnew);

   const double& told = src_told_data.getTime(); 

   if (tbox::MathUtilities<double>::equalEps(s_time, told)) {

      /*
       * If the destination time is same as told, do a regular copy.
       */

      scratch_data.copy(src_told_data, *d_overlap);

   } else if (d_overlap->getSourceOffset() ==
       hier::IntVector::getZero(d_box.getDim()) &&
       d_overlap->getTransformation().getRotation() ==
       hier::Transformation::NO_ROTATE) {


      /*
       * If there is no offset between the source and destination, then
       * time interpolate directly to the destination patchdata.  Otherwise,
       * time interpolate into a temporary patchdata and copy the result
       * to the destination patchdata.
       */

      timeInterpolate(scratch_data, *d_overlap, src_told_data, src_tnew_data);

   } else {

      hier::Box temporary_box(d_box.getDim());
      temporary_box.initialize(d_box,
                               d_src_patch->getBox().getLocalId(),
                               tbox::SAMRAI_MPI::getInvalidRank());

      hier::Patch temporary_patch(
         temporary_box,
         d_src_patch->getPatchDescriptor());

      std::shared_ptr<hier::PatchData> temp(
         d_src_patch->getPatchDescriptor()
         ->getPatchDataFactory(d_refine_data[d_item_id]->d_src_told)
         ->allocate(temporary_patch));

      temp->setTime(s_time);

      timeInterpolate(*temp, *d_overlap, src_told_data, src_tnew_data);

      scratch_data.copy(*temp, *d_overlap);

   }

}

void
RefineTimeTransaction::timeInterpolate(
   hier::PatchData& pd_dst,
   const hier::BoxOverlap& overlap,
   const hier::PatchData& pd_old,
   const std::shared_ptr<hier::PatchData>& pd_new)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(pd_dst, pd_old);
   TBOX_ASSERT(tbox::MathUtilities<double>::equalEps(pd_dst.getTime(), s_time));

   if (tbox::MathUtilities<double>::equalEps(pd_old.getTime(), s_time)) {
      d_refine_data[d_item_id]->
      d_optime->timeInterpolate(pd_dst, d_box, overlap, pd_old, pd_old);
   } else {

      TBOX_ASSERT(pd_new);
      TBOX_ASSERT_OBJDIM_EQUALITY2(pd_dst, *pd_new);
      TBOX_ASSERT(pd_old.getTime() < s_time);
      TBOX_ASSERT(pd_new->getTime() >= s_time);

      d_refine_data[d_item_id]->
      d_optime->timeInterpolate(pd_dst, d_box, overlap, pd_old, *pd_new);
   }
}

/*
 *************************************************************************
 *
 * Function to print state of transaction.
 *
 *************************************************************************
 */

void
RefineTimeTransaction::printClassData(
   std::ostream& stream) const
{
   stream << "Refine Time Transaction" << std::endl;
   stream << "   transaction time:        " << s_time << std::endl;
   stream << "   refine item array:        "
          << (RefineClasses::Data *)d_refine_data[d_item_id] << std::endl;
   stream << "   destination patch rank:        " << d_dst_patch_rank
          << std::endl;
   stream << "   source patch rank:             " << d_src_patch_rank
          << std::endl;
   stream << "   time interpolation box:  " << d_box << std::endl;
   stream << "   refine item id :  " << d_item_id << std::endl;
   if (d_refine_data) {
      auto& optime = *d_refine_data[d_item_id]->d_optime;
      stream << "   destination patch data id:  "
             << d_refine_data[d_item_id]->d_scratch << std::endl;
      stream << "   source (old) patch data id: "
             << d_refine_data[d_item_id]->d_src_told << std::endl;
      stream << "   source (new) patch data id: "
             << d_refine_data[d_item_id]->d_src_tnew << std::endl;
      stream << "   time interpolation name id: "
             << typeid(optime).name() << std::endl;
   }
   stream << "   incoming bytes:          " << d_incoming_bytes << std::endl;
   stream << "   outgoing bytes:          " << d_outgoing_bytes << std::endl;
   stream << "   destination patch:           "
          << d_dst_patch.get() << std::endl;
   stream << "   source level:           "
          << d_src_patch.get() << std::endl;
   stream << "   overlap:                 " << std::endl;
   d_overlap->print(stream);
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
