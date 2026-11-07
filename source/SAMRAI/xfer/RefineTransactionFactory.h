/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface for factory objects that create transactions for
 *                refine schedules.
 *
 ************************************************************************/

#ifndef included_xfer_RefineTransactionFactory
#define included_xfer_RefineTransactionFactory

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Transaction.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineClasses.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Abstract base class defining the interface for all concrete
 * transaction factory objects that generate data transaction objects used with
 * a RefineSchedule object.  A concrete subclass will allocate new transaction
 * objects.  This class is an example of the ``Abstract Factory'' method
 * described in the Design Patterns book by Gamma, et al.
 *
 * To add a new type of Transaction object MyRefineTransaction:
 *
 * -# Implement a concrete RefineTransactionFactory object as a subclass
 *       that is derived from this RefineTransactionFactory base class.
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

class RefineTransactionFactory
{
public:
   /*!
    * @brief Default constructor.
    */
   RefineTransactionFactory();

   /*!
    * @brief Virtual destructor.
    */
   virtual ~RefineTransactionFactory();

   /*!
    * @brief Pure virtual function to allocate a concrete refine transaction
    * object.  This routine is called by the refine schedule during
    * construction of the schedule.
    *
    * @param dst_level      std::shared_ptr to destination patch level.
    * @param src_level      std::shared_ptr to source patch level.
    * @param overlap        std::shared_ptr to overlap region between
    *                       patches.
    * @param dst_box        Destination Box in destination patch level.
    * @param src_box        Source Box in source patch level.
    * @param refine_data    Pointer to array of refine data items
    * @param item_id        Integer index of RefineClass::Data item associated
    *                       with transaction.
    * @param box            Optional const reference to box defining region of
    *                       refine transaction.  Default is an empty box.
    * @param use_time_interpolation  Optional boolean flag indicating whether
    *                       the refine transaction involves time interpolation.
    *                       Default is false.
    */
   virtual std::shared_ptr<tbox::Transaction>
   allocate(
      const std::shared_ptr<hier::PatchLevel>& dst_level,
      const std::shared_ptr<hier::PatchLevel>& src_level,
      const std::shared_ptr<hier::BoxOverlap>& overlap,
      const hier::Box& dst_box,
      const hier::Box& src_box,
      const RefineClasses::Data** refine_data,
      int item_id,
      const hier::Box& box,
      bool use_time_interpolation = false) const = 0;

   std::shared_ptr<tbox::Transaction>
   allocate(
      const std::shared_ptr<hier::PatchLevel>& dst_level,
      const std::shared_ptr<hier::PatchLevel>& src_level,
      const std::shared_ptr<hier::BoxOverlap>& overlap,
      const hier::Box& dst_box,
      const hier::Box& src_box,
      const RefineClasses::Data** refine_data,
      int item_id) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY4(*dst_level,
         *src_level,
         dst_box,
         src_box);
      return allocate(
         dst_level,
         src_level,
         overlap,
         dst_box,
         src_box,
         refine_data,
         item_id,
         hier::Box::getEmptyBox(src_level->getDim()),
         false);
   }

   /*!
    * @brief Virtual function to set simulation time for transaction objects.
    * This operation is called by the refine schedule during the execution of
    * the RefineSchedule::fillData() routine before data communication
    * operations begin.  This function is optional for the concrete transaction
    * factory object.  The default implementation is a no-op.
    */
   virtual void
   setTransactionTime(
      double fill_time);

   /*!
    * @brief Virtual function allowing transaction factory to preprocess
    * scratch space data before transactactions use it if they need to.  This
    * function is optional for the concrete transaction factory object.
    * The default implementation is a no-op.
    *
    * @param level        std::shared_ptr to patch level holding scratch data.
    * @param fill_time    Double value of simulation time corresponding to
    *                     RefineSchedule operations.
    * @param preprocess_vector Const reference to ComponentSelector that
    *                     indicates patch data array indices of scratch patch
    *                     data objects to preprocess.
    */
   virtual void
   preprocessScratchSpace(
      const std::shared_ptr<hier::PatchLevel>& level,
      double fill_time,
      const hier::ComponentSelector& preprocess_vector) const = 0;

private:
   // The following two functions are not implemented
   RefineTransactionFactory(
      const RefineTransactionFactory&);
   RefineTransactionFactory&
   operator = (
      const RefineTransactionFactory&);

};

}
}

#endif
