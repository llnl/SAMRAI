/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Communication transaction for data copies during data refining
 *
 ************************************************************************/
#include "SAMRAI/xfer/RefineCopyTransaction.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

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
 * Constructor sets state of transaction.
 *
 *************************************************************************
 */

RefineCopyTransaction::RefineCopyTransaction(
   const std::shared_ptr<hier::PatchLevel>& dst_level,
   const std::shared_ptr<hier::PatchLevel>& src_level,
   const std::shared_ptr<hier::BoxOverlap>& overlap,
   const hier::Box& dst_box,
   const hier::Box& src_box,
   const RefineClasses::Data** refine_data,
   int item_id):
   d_dst_patch_rank(dst_box.getOwnerRank()),
   d_src_patch_rank(src_box.getOwnerRank()),
   d_overlap(overlap),
   d_refine_data(refine_data),
   d_item_id(item_id),
   d_incoming_bytes(0),
   d_outgoing_bytes(0)
{
   TBOX_ASSERT(dst_level);
   TBOX_ASSERT(src_level);
   TBOX_ASSERT(overlap);
   TBOX_ASSERT(dst_box.getLocalId() >= 0);
   TBOX_ASSERT(src_box.getLocalId() >= 0);
   TBOX_ASSERT(item_id >= 0);
   TBOX_ASSERT(refine_data[item_id] != 0);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst_level,
      *src_level,
      dst_box,
      src_box);

   // Note: s_num_coarsen_items cannot be used at this point!

   if (d_dst_patch_rank == dst_level->getBoxLevel()->getMPI().getRank()) {
      d_dst_patch = dst_level->getPatch(dst_box.getGlobalId());
   }
   if (d_src_patch_rank == dst_level->getBoxLevel()->getMPI().getRank()) {
      d_src_patch = src_level->getPatch(src_box.getGlobalId());
   }
}

RefineCopyTransaction::~RefineCopyTransaction()
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
RefineCopyTransaction::canEstimateIncomingMessageSize()
{
   bool can_estimate = false;
   if (d_src_patch) {
      can_estimate =
         d_src_patch->getPatchData(d_refine_data[d_item_id]->d_src)
         ->canEstimateStreamSizeFromBox();
   } else {
      can_estimate =
         d_dst_patch->getPatchData(d_refine_data[d_item_id]->d_scratch)
         ->canEstimateStreamSizeFromBox();
   }
   return can_estimate;
}

size_t
RefineCopyTransaction::computeIncomingMessageSize()
{
   d_incoming_bytes =
      d_dst_patch->getPatchData(d_refine_data[d_item_id]->d_scratch)
      ->getDataStreamSize(*d_overlap);
   return d_incoming_bytes;
}

size_t
RefineCopyTransaction::computeOutgoingMessageSize()
{
   d_outgoing_bytes =
      d_src_patch->getPatchData(d_refine_data[d_item_id]->d_src)
      ->getDataStreamSize(*d_overlap);
   return d_outgoing_bytes;
}

int
RefineCopyTransaction::getSourceProcessor() {
   return d_src_patch_rank;
}

int
RefineCopyTransaction::getDestinationProcessor() {
   return d_dst_patch_rank;
}

void
RefineCopyTransaction::packStream(
   tbox::MessageStream& stream)
{
   d_src_patch->getPatchData(d_refine_data[d_item_id]->d_src)
   ->packStreamFuseable(stream, *d_overlap);
}

void
RefineCopyTransaction::unpackStream(
   tbox::MessageStream& stream)
{
   d_dst_patch->getPatchData(d_refine_data[d_item_id]->d_scratch)
   ->unpackStreamFuseable(stream, *d_overlap);
}

void
RefineCopyTransaction::copyLocalData()
{
   hier::PatchData& dst_data =
      *d_dst_patch->getPatchData(d_refine_data[d_item_id]->d_scratch);

   const hier::PatchData& src_data =
      *d_src_patch->getPatchData(d_refine_data[d_item_id]->d_src);

   dst_data.copyFuseable(src_data, *d_overlap);
}

/*
 *************************************************************************
 *
 * Function to print state of transaction.
 *
 *************************************************************************
 */

void
RefineCopyTransaction::printClassData(
   std::ostream& stream) const
{
   stream << "Refine Copy Transaction" << std::endl;
   stream << "   refine item:        "
          << (RefineClasses::Data *)d_refine_data[d_item_id] << std::endl;
   stream << "   destination patch rank:       " << d_dst_patch_rank
          << std::endl;
   stream << "   source patch rank:            " << d_src_patch_rank
          << std::endl;
   if (d_refine_data) {
      stream << "   destination patch data id: "
             << d_refine_data[d_item_id]->d_scratch << std::endl;
      stream << "   source patch data id:      "
             << d_refine_data[d_item_id]->d_src << std::endl;
   }
   stream << "   incoming bytes:         " << d_incoming_bytes << std::endl;
   stream << "   outgoing bytes:         " << d_outgoing_bytes << std::endl;
   stream << "   destination patch:           "
          << d_dst_patch.get() << std::endl;
   stream << "   source patch:           "
          << d_src_patch.get() << std::endl;
   stream << "   overlap:                " << std::endl;
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
