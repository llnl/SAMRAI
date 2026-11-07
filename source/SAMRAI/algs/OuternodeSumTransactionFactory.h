/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Factory for creating outernode sum transaction objects
 *
 ************************************************************************/

#ifndef included_algs_OuternodeSumTransactionFactory
#define included_algs_OuternodeSumTransactionFactory

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineClasses.h"
#include "SAMRAI/xfer/RefineTransactionFactory.h"

#include <memory>

namespace SAMRAI {
namespace algs {

/*!
 * @brief Concrete subclass of the xfer::RefineTransactionFactory base class
 * that allocates transaction outernode sum objects for a xfer::RefineSchedule
 * object.
 *
 * @see xfer::RefineTransactionFactory
 * @see xfer::OuternodeSumTransaction
 */

class OuternodeSumTransactionFactory:public xfer::RefineTransactionFactory
{
public:
   /*!
    * @brief Default constructor.
    */
   OuternodeSumTransactionFactory();

   /*!
    * @brief Virtual destructor for base class.
    */
   virtual ~OuternodeSumTransactionFactory();

   /*!
    * @brief Set the array of xfer::RefineClasses::Data items used by the
    * transactions.
    */
   void
   setRefineItems(
      const xfer::RefineClasses::Data * const * refine_items);

   /*!
    * @brief Clear the array of xfer::RefineClasses::Data items used by the
    * transactions.
    */
   void
   unsetRefineItems();

   /*!
    * @brief Allocate an OuternodeSumTransaction object.
    *
    * @param dst_level      std::shared_ptr to destination patch level.
    * @param src_level      std::shared_ptr to source patch level.
    * @param overlap        std::shared_ptr to overlap region between
    *                       patches.
    * @param dst_node       Destination Box in destination patch level.
    * @param src_node       Source Box in source patch level.
    * @param refine_data    Pointer to array of refine data items
    * @param item_id        Integer index of xfer::RefineClasses::Data item
    *                       associated with transaction.
    * @param box            Const reference to box defining region of
    *                       refine transaction.  Use following allocate method
    *                       if not needed.
    * @param use_time_interpolation  Optional boolean flag indicating whether
    *                       the refine transaction involves time interpolation.
    *                       Default is false.
    *
    * @pre dst_level
    * @pre src_level
    * @pre overlap
    * @pre dst_node.getLocalId() >= 0
    * @pre src_node.getLocalId() >= 0
    * @pre item_id >= 0
    * @pre (dst_level->getDim() == src_level->getDim()) &&
    *      (dst_level->getDim() == dst_node.getDim()) &&
    *      (dst_level->getDim() == src_node.getDim())
    */
   std::shared_ptr<tbox::Transaction>
   allocate(
      const std::shared_ptr<hier::PatchLevel>& dst_level,
      const std::shared_ptr<hier::PatchLevel>& src_level,
      const std::shared_ptr<hier::BoxOverlap>& overlap,
      const hier::Box& dst_node,
      const hier::Box& src_node,
      const xfer::RefineClasses::Data ** refine_data,
      int item_id,
      const hier::Box& box,
      bool use_time_interpolation = false) const;

   /*!
    * @brief Allocate an OuternodeSumTransaction object.
    *
    * Same as previous allocate routine but with default empty box and no
    * timer interpolation.
    *
    * @pre dst_level
    * @pre src_level
    * @pre overlap
    * @pre dst_node.getLocalId() >= 0
    * @pre src_node.getLocalId() >= 0
    * @pre item_id >= 0
    * @pre (dst_level->getDim() == src_level->getDim()) &&
    *      (dst_level->getDim() == dst_node.getDim()) &&
    *      (dst_level->getDim() == src_node.getDim())
    */
   std::shared_ptr<tbox::Transaction>
   allocate(
      const std::shared_ptr<hier::PatchLevel>& dst_level,
      const std::shared_ptr<hier::PatchLevel>& src_level,
      const std::shared_ptr<hier::BoxOverlap>& overlap,
      const hier::Box& dst_node,
      const hier::Box& src_node,
      const xfer::RefineClasses::Data ** refine_data,
      int item_id) const;

   /*!
    * @brief Function to initialize scratch space data for the sum transactions
    * (patch data components indicated by the component selector) to zero.
    *
    * @param level        std::shared_ptr to patch level holding scratch
    *                     data.
    * @param fill_time    Double value of simulation time at which preprocess
    *                     operation is called.
    * @param preprocess_vector Const reference to hier::ComponentSelector
    *                     indicating patch data array indices of scratch patch
    *                     data objects to preprocess.
    *
    * @pre level
    */
   void
   preprocessScratchSpace(
      const std::shared_ptr<hier::PatchLevel>& level,
      double fill_time,
      const hier::ComponentSelector& preprocess_vector) const;

private:
   // The following two functions are not implemented
   OuternodeSumTransactionFactory(
      const OuternodeSumTransactionFactory&);
   OuternodeSumTransactionFactory&
   operator = (
      const OuternodeSumTransactionFactory&);

};

}
}
#endif
