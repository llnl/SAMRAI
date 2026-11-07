/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Concrete factory to create standard copy and time transactions
 *                for refine schedules.
 *
 ************************************************************************/
#include "SAMRAI/xfer/StandardRefineTransactionFactory.h"

#include "SAMRAI/xfer/RefineCopyTransaction.h"
#include "SAMRAI/xfer/RefineTimeTransaction.h"


namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Default constructor and destructor.
 *
 *************************************************************************
 */

StandardRefineTransactionFactory::StandardRefineTransactionFactory()
{
}

StandardRefineTransactionFactory::~StandardRefineTransactionFactory()
{
}

void
StandardRefineTransactionFactory::setTransactionTime(
   double fill_time)
{
   RefineTimeTransaction::setTransactionTime(fill_time);
}

/*
 *************************************************************************
 *
 * Allocate appropriate transaction object.
 *
 *************************************************************************
 */

std::shared_ptr<tbox::Transaction>
StandardRefineTransactionFactory::allocate(
   const std::shared_ptr<hier::PatchLevel>& dst_level,
   const std::shared_ptr<hier::PatchLevel>& src_level,
   const std::shared_ptr<hier::BoxOverlap>& overlap,
   const hier::Box& dst_box,
   const hier::Box& src_box,
   const RefineClasses::Data** refine_data,
   int item_id,
   const hier::Box& box,
   bool use_time_interpolation) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY5(*dst_level,
      *src_level,
      dst_box,
      src_box,
      box);

   if (use_time_interpolation) {

      return std::make_shared<RefineTimeTransaction>(
                dst_level,
                src_level,
                overlap,
                dst_box,
                src_box,
                box,
                refine_data,
                item_id);

   } else {

      return std::make_shared<RefineCopyTransaction>(
                dst_level,
                src_level,
                overlap,
                dst_box,
                src_box,
                refine_data,
                item_id);

   }
}

void
StandardRefineTransactionFactory::preprocessScratchSpace(
   const std::shared_ptr<hier::PatchLevel>& level,
   double fill_time,
   const hier::ComponentSelector& preprocess_vector) const
{
   NULL_USE(level);
   NULL_USE(fill_time);
   NULL_USE(preprocess_vector);
}

}
}
