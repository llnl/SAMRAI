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

#ifndef included_xfer_StandardRefineTransactionFactory
#define included_xfer_StandardRefineTransactionFactory

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Transaction.h"
#include "SAMRAI/xfer/RefineClasses.h"
#include "SAMRAI/xfer/RefineTransactionFactory.h"

#include <memory>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Concrete subclass of RefineTransactionFactory base class that
 * allocates RefineCopyTransaction and RefineTimeTransaction objects for a
 * RefineSchedule object.
 *
 * @see RefineCopyTransaction
 * @see RefineTimeTransaction
 * @see RefineTransactionFactory
 */

class StandardRefineTransactionFactory:public RefineTransactionFactory
{
public:
   /*!
    * @brief Default constructor.
    */
   StandardRefineTransactionFactory();

   /*!
    * @brief Virtual destructor.
    */
   virtual ~StandardRefineTransactionFactory();

   /*!
    * @brief Set simulation time used by the refine time transaction objects.
    */
   virtual void
   setTransactionTime(
      double fill_time);

   /*!
    * @brief Allocate an appropriate refine copy or time transaction object.
    * When time interpolation flag is passed as true a RefineTimeTransaction
    * object will be created.  Otherwise, a RefineCopyTransaction aill be
    * created.
    *
    * @param dst_level      std::shared_ptr to destination patch level.
    * @param src_level      std::shared_ptr to source patch level.
    * @param overlap        std::shared_ptr to overlap region between
    *                       patches.
    * @param dst_box        Destination Box in destination patch level.
    * @param src_box        Source Box in source patch level.
    * @param refine_data    Pointer to array of refine data items.
    * @param item_id        Integer index of RefineClass::Data item associated
    *                       with transaction.
    * @param box            Optional const reference to box defining region of
    *                       refine transaction.  Default is an empty box.
    * @param use_time_interpolation  Optional boolean flag indicating whether
    *                       the refine transaction involves time interpolation.
    *                       Default is false.
    *
    * @pre (dst_level->getDim() == src_level->getDim()) &&
    *      (dst_level->getDim() == dst_box.getDim()) &&
    *      (dst_level->getDim() == src_box.getDim()) &&
    *      (dst_level->getDim() == box.getDim())
    */
   virtual std::shared_ptr<tbox::Transaction>
   allocate(
      const std::shared_ptr<hier::PatchLevel>& dst_level,
      const std::shared_ptr<hier::PatchLevel>& src_level,
      const std::shared_ptr<hier::BoxOverlap>& overlap,
      const hier::Box& dst_box,
      const hier::Box& src_box,
      const RefineClasses::Data ** refine_data,
      int item_id,
      const hier::Box& box,       // Default in v 2.x  = hier::Box()
      bool use_time_interpolation = false) const;

   /*!
    * @brief Virtual function allowing transaction factory to preprocess
    * scratch space data before transactactions use it if they need to.  This
    * function is optional for the concrete transaction factory object.
    * The default implementation is a no-op.
    *
    * @param level        std::shared_ptr to patch level holding scratch
    *                     data.
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
      const hier::ComponentSelector& preprocess_vector) const;

private:
   // The following two functions are not implemented
   StandardRefineTransactionFactory(
      const StandardRefineTransactionFactory&);
   StandardRefineTransactionFactory&
   operator = (
      const StandardRefineTransactionFactory&);

};

}
}
#endif
