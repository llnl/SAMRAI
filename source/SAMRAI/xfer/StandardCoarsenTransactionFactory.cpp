/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Concrete factory for create standard copy transactions
 *                for coarsen schedules.
 *
 ************************************************************************/
#include "SAMRAI/xfer/StandardCoarsenTransactionFactory.h"

#include "SAMRAI/xfer/CoarsenCopyTransaction.h"


namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Default constructor and destructor.
 *
 *************************************************************************
 */

StandardCoarsenTransactionFactory::StandardCoarsenTransactionFactory()
{
}

StandardCoarsenTransactionFactory::~StandardCoarsenTransactionFactory()
{
}

/*
 *************************************************************************
 *
 * Allocate appropriate transaction object.
 *
 *************************************************************************
 */

std::shared_ptr<tbox::Transaction>
StandardCoarsenTransactionFactory::allocate(
   const std::shared_ptr<hier::PatchLevel>& dst_level,
   const std::shared_ptr<hier::PatchLevel>& src_level,
   const std::shared_ptr<hier::BoxOverlap>& overlap,
   const hier::Box& dst_box,
   const hier::Box& src_box,
   const CoarsenClasses::Data** coarsen_data,
   int item_id) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst_level,
      *src_level,
      dst_box,
      src_box);

   return std::make_shared<CoarsenCopyTransaction>(
             dst_level,
             src_level,
             overlap,
             dst_box,
             src_box,
             coarsen_data,
             item_id);
}

}
}
