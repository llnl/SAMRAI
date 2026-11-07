/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface for factory objects that create transactions for
 *                oarsen schedules.
 *
 ************************************************************************/

#ifndef included_xfer_CoarsenTransactionFactory
#define included_xfer_CoarsenTransactionFactory

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Transaction.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/CoarsenClasses.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Abstract base class defining the interface for all concrete
 * transaction factory objects that generate data transaction objects used with
 * a CoarsenSchedule object.  A concrete subclass will allocate new transaction
 * objects.  This class is an example of the ``Abstract Factory'' method
 * described in the Design Patterns book by Gamma, et al.
 *
 * To add a new type of Transaction object MyCoarsenTransaction:
 *
 * -# Implement a concrete CoarsenTransactionFactory object as a subclass
 *       that is derived from this CoarsenTransactionFactory base class.
 *       Implement the abstract virtual functions as appropriate for the
 *       concrete subclass; in particular, the allocate() function must return
 *       a new instance of the desired transaction object.
 * -# The type of the transaction allocated by the concrete factory is a
 *       Transaction.  Thus, the new transaction object must be derived
 *       from the Transaction base class and implement the abstract
 *       virtual functions declared by the base class.
 *
 * @see tbox::Transaction
 */

class CoarsenTransactionFactory
{
public:
   /*!
    * @brief Default constructor.
    */
   CoarsenTransactionFactory();

   /*!
    * @brief Virtual destructor.
    */
   virtual ~CoarsenTransactionFactory();

   /*!
    * @brief Pure virtual function to allocate a concrete coarsen transaction
    * object.  This routine is called by the coarsen schedule during
    * construction of the schedule.
    *
    * @param dst_level      std::shared_ptr to destination patch level.
    * @param src_level      std::shared_ptr to source patch level.
    * @param overlap        std::shared_ptr to overlap region between patches.
    * @param dst_box        Destination Box in destination patch level.
    * @param src_box        Source Box in source patch level.
    * @param coarsen_data   Pointer to array of coarsen data items
    * @param item_id        Integer index of CoarsenClass::Data item associated
    *                       with transaction.
    */
   virtual std::shared_ptr<tbox::Transaction>
   allocate(
      const std::shared_ptr<hier::PatchLevel>& dst_level,
      const std::shared_ptr<hier::PatchLevel>& src_level,
      const std::shared_ptr<hier::BoxOverlap>& overlap,
      const hier::Box& dst_box,
      const hier::Box& src_box,
      const CoarsenClasses::Data ** coarsen_data,
      int item_id) const = 0;

private:
   CoarsenTransactionFactory(
      const CoarsenTransactionFactory&);                  // not implemented
   CoarsenTransactionFactory&
   operator = (
      const CoarsenTransactionFactory&);                  // not implemented

};

}
}
#endif
