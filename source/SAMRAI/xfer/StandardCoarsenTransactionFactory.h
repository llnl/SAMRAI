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

#ifndef included_xfer_StandardCoarsenTransactionFactory
#define included_xfer_StandardCoarsenTransactionFactory

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Transaction.h"
#include "SAMRAI/xfer/CoarsenClasses.h"
#include "SAMRAI/xfer/CoarsenTransactionFactory.h"

#include <memory>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Concrete subclass of CoarsenTransactionFactory base class that
 * allocates CoarsenCopyTransaction objects for a CoarsenSchedule object.
 *
 * @see CoarsenCopyTransaction
 */

class StandardCoarsenTransactionFactory:public CoarsenTransactionFactory
{
public:
   /*!
    * @brief Default constructor.
    */
   StandardCoarsenTransactionFactory();

   /*!
    * @brief Virtual destructor.
    */
   virtual ~StandardCoarsenTransactionFactory();

   /*!
    * @brief Allocate a CoarsenCopyTransaction object.
    *
    * @param dst_level      std::shared_ptr to destination patch level.
    * @param src_level      std::shared_ptr to source patch level.
    * @param overlap        std::shared_ptr to overlap region between
    *                       patches.
    * @param dst_box        Destination Box in destination patch level.
    * @param src_box        Source Box in source patch level.
    * @param coarsen_data   Pointer to array of coarsen data items
    * @param item_id        Integer index of CoarsenClass::Data item associated
    *                       with transaction.
    *
    * @pre (dst_level->getDim() == src_level->getDim()) &&
    *      (dst_level->getDim() == dst_box.getDim()) &&
    *      (dst_level->getDim() == src_box.getDim())
    */
   virtual std::shared_ptr<tbox::Transaction>
   allocate(
      const std::shared_ptr<hier::PatchLevel>& dst_level,
      const std::shared_ptr<hier::PatchLevel>& src_level,
      const std::shared_ptr<hier::BoxOverlap>& overlap,
      const hier::Box& dst_box,
      const hier::Box& src_box,
      const CoarsenClasses::Data ** coarsen_data,
      int item_id) const;

private:
   // The following two functions are not implemented
   StandardCoarsenTransactionFactory(
      const StandardCoarsenTransactionFactory&);
   StandardCoarsenTransactionFactory&
   operator = (
      const StandardCoarsenTransactionFactory&);

};

}
}
#endif
