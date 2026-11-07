/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Container of loads for TreeLoadBalancer.
 *
 ************************************************************************/

#ifndef included_mesh_TransitLoad
#define included_mesh_TransitLoad

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/MappingConnector.h"
#include "SAMRAI/tbox/MessageStream.h"

#include <iostream>
#include <string>

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Base class for container of work in transit during
 * partitioning.
 *
 * TransitLoad base class follows the prototype design pattern.
 * Subclasses must implement clone(), initialize(), and a copy
 * constructor according to this design pattern.
 *
 * TransitLoad has dual responsibilities.
 *
 * - First, it's a container of loads that move around.  The
 * implementation may represent the work however it wants and should
 * be able to shift load from one container to another.  See
 * adjustLoad().
 *
 * - Second, it generates the mappings between the pre- and
 * post-balance BoxLevels.  It must have this responsibility because
 * the implementation alone knows how the work is represented.  See
 * assignToLocalAndPopulateMaps().
 *
 * For usage examples, see TreeLoadBalancer.
 */

class TransitLoad
{

public:
   TransitLoad();

   TransitLoad(
      const TransitLoad& other);

   virtual ~TransitLoad() {
   }

   //@{
   //! @name Methods required by prototype design pattern

   //! @brief Clone object according to design pattern.
   virtual TransitLoad *
   clone() const = 0;

   //! @brief Initialize object according to design pattern.
   virtual void
   initialize() = 0;

   //@}

   //@{
   //! @brief Container-like characteristics

   //! @brief Return the total load contained.
   virtual double
   getSumLoad() const = 0;

   //! @brief Insert all boxes from the given BoxContainer.
   virtual void
   insertAll(
      const hier::BoxContainer& other) = 0;

   /*! @brief Insert all boxes while applying a minimum load value
    *
    * This virtual method allows child classes to implement a minimum load
    * which forces boxes to have this value as their minimum load even if
    * their computed load is smaller.  As this is an optional feature that
    * can be implemented in child classes, the default implementation
    * ignores the minimum load value and is equivalent to insertAll().
    *
    * @param other  box container with boxes to inserted into this object.
    * @param minimum_load  artificial minimum load value
    *                      (ignored in default implementation)
    */
   virtual void
   insertAllWithArtificialMinimum(
      const hier::BoxContainer& other,
      double minimum_load)
   {
      // Child classes can implement usage of the minimum_load.  This
      // default implementation will ignore it.
      NULL_USE(minimum_load);
      insertAll(other);
   }

   /*!
    * @brief Insert all boxes from the given TransitLoad.
    *
    * Changes to other as an implementation side-effect is allowed.
    * This and other are guaranteed to be the same concrete type.
    *
    * @param other [in] Other TransitLoad container whose
    * implementation matches this one.
    */
   virtual void
   insertAll(
      TransitLoad& other) = 0;

   /*!
    * @brief Insert Boxes from the given container but preserve
    * load values currently held by the TransitLoad.
    */
   virtual void
   insertAllWithExistingLoads(
      const hier::BoxContainer& box_container) = 0;

   /*!
    * @brief Set Workload values in the the TransitLoad object
    *
    * @param patch_level  Level holding workload data
    * @param work_data_id Patch data id for workload data
    */
   virtual void
   setWorkload(
      const hier::PatchLevel& patch_level,
      const int work_data_id) = 0;

   //! @brief Return number of items in this container.
   virtual size_t
   getNumberOfItems() const = 0;

   //! @brief Return number of processes contributing to the contents.
   virtual size_t
   getNumberOfOriginatingProcesses() const = 0;

   //! @brief Whether container is empty.
   virtual bool empty() const {
      return getNumberOfItems() == 0;
   }

   //! @brief Empty the container.
   virtual void
   clear() = 0;

   //@}

   //@{
   //! @name Packing/unpacking for communication.

   //! @brief Put content into MessageStream.
   virtual void
   putToMessageStream(
      tbox::MessageStream& msg) const = 0;

   //! @brief Add to content from MessageStream.
   virtual void
   getFromMessageStream(
      tbox::MessageStream& msg) = 0;

   friend tbox::MessageStream& operator << (
      tbox::MessageStream& msg, const TransitLoad& transit_load) {
      transit_load.putToMessageStream(msg);
      return msg;
   }

   friend tbox::MessageStream& operator >> (
      tbox::MessageStream& msg, TransitLoad& transit_load) {
      transit_load.getFromMessageStream(msg);
      return msg;
   }
   //@}

   /*!
    * @brief Adjust the load in this TransitSet by moving work
    * between it and another TransitSet.
    *
    * @param[in,out] hold_bin Holding bin for reserve load.  hold_bin
    * implementation is guaranteed to match this one.
    *
    * @param[in] ideal_load The load that this bin should have.
    *
    * @param[in] low_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @param[in] high_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @return Net load added to this TransitSet.  If negative, load
    * decreased.
    */
   virtual double
   adjustLoad(
      TransitLoad& hold_bin,
      double ideal_load,
      double low_load,
      double high_load) = 0;

   /*!
    * @brief Assign contents to local process.
    *
    * This method may use communication.
    *
    * @param [in,out] balanced_box_level Empty BoxLevel to populate with
    * the contents of this TransitLoad.
    *
    * @param[in] unbalanced_box_level
    *
    * @param [in] flexible_load_tol
    *
    * @param [in] alt_mpi Alternate SAMRAI_MPI to use for communication if
    * needed, overriding that in balanced_box_level.
    */
   virtual void
   assignToLocal(
      hier::BoxLevel& balanced_box_level,
      const hier::BoxLevel& unbalanced_box_level,
      double flexible_load_tol = 0.0,
      const tbox::SAMRAI_MPI& alt_mpi = tbox::SAMRAI_MPI(MPI_COMM_NULL)) = 0;

   /*!
    * @brief Assign contents to local process and populate the
    * balanced<==>unbalanced maps.
    *
    * This method may use communication.
    *
    * @param [in,out] balanced_box_level Empty BoxLevel to populate with
    * the contents of this TransitLoad.
    *
    * @param [in,out] balanced_to_unbalanced Empty Connector to populate
    * with the balanced--->unbalanced edges.
    *
    * @param [in,out] unbalanced_to_balanced Empty Connector to populate
    * with the unbalanced--->balanced edges.
    *
    * @param [in] flexible_load_tol
    *
    * @param [in] alt_mpi Alternate SAMRAI_MPI to use for communication if
    * needed, overriding that in balanced_box_level.
    */
   virtual void
   assignToLocalAndPopulateMaps(
      hier::BoxLevel& balanced_box_level,
      hier::MappingConnector& balanced_to_unbalanced,
      hier::MappingConnector& unbalanced_to_balanced,
      double flexible_load_tol = 0.0,
      const tbox::SAMRAI_MPI& alt_mpi = tbox::SAMRAI_MPI(MPI_COMM_NULL)) = 0;

   //@{
   //! @name Parameters in box breaking

   //! @brief Whether to allow box breaking.
   void setAllowBoxBreaking(bool allow_box_breaking) {
      d_allow_box_breaking = allow_box_breaking;
   }

   //! @brief Whether object may break up boxes.
   bool getAllowBoxBreaking() const {
      return d_allow_box_breaking;
   }

   /*!
    * @brief Set threshold for resisting small and thin boxes.
    *
    * This is not a hard limit, but implementations should use it to
    * determine when a box is getting too small in any coordinate
    * direction.  Default initial threshold is very small (practically
    * unlimited).
    */
   void setThresholdWidth(double threshold_width) {
      d_threshold_width = threshold_width;
   }

   //! @brief Return the threshold width for resisting small and thin boxes.
   double getThresholdWidth() const {
      return d_threshold_width;
   }

   //@}

   //@{
   //! @name Diagnostic support

   //! @brief Set prefix of internal diagnostic timer names.
   virtual void
   setTimerPrefix(
      const std::string&) {
   }

   //! @brief Object printing method to aid in debugging.
   virtual void
   recursivePrint(
      std::ostream& co = tbox::plog,
      const std::string& border = std::string(),
      int detail_depth = 1) const;

   //@}

private:
   bool d_allow_box_breaking;
   double d_threshold_width;

};

}
}

#endif
