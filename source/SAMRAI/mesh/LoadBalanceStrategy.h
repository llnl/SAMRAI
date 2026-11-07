/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface for box load balancing routines.
 *
 ************************************************************************/

#ifndef included_mesh_LoadBalanceStrategy
#define included_mesh_LoadBalanceStrategy

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/RankGroup.h"

#include <memory>

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Class LoadBalanceStrategy is an abstract base class that
 * defines a Strategy pattern interface for operations that load
 * balance patches on a single AMR patch hierarchy level.  Typically,
 * such operations are invoked after the domain of a new hierarchy
 * level is determined (e.g., via some error estimation procedure) and
 * is applied to the collection of boxes that describe the domain.  The
 * load balancing process produces a set of boxes from which patches on
 * the new level are created and a processor mapping describing how the
 * new patches are mapped to processors.
 *
 * @see hier::PatchLevel
 */

class LoadBalanceStrategy
{
public:
   /*!
    * This virtual destructor does nothing interesting.
    */
   virtual ~LoadBalanceStrategy();

   /*!
    * Indicate whether load balancing procedure for given level depends on
    * patch data on mesh.  This can be used to determine whether a level
    * needs to be rebalanced although its box configuration is unchanged.
    *
    * @return Boolean value of true if load balance routines for level
    *        depend on patch data; false otherwise.
    *
    * @param level_number Integer level number.
    */
   virtual bool
   getLoadBalanceDependsOnPatchData(
      int level_number) const = 0;

   /*!
    * @brief Given a BoxLevel, representing the domain of a specified
    * level in the AMR hierarchy, generate a new BoxLevel from which the
    * patches for the level may be formed and update two Connectors
    * incident on the changed BoxLevel.
    *
    * The union of the boxes in the balance_box_level is the same
    * before and after the the method call.
    *
    * @param[in,out] balance_box_level Input BoxLevel.  On input, this is the pre-balance
    *  BoxLevel.  On output, it is the balanced BoxLevel.
    *
    * @param[in,out] balance_to_anchor Connector between the balance_box_level and
    *  some given "anchor box_level".
    *  This must be accurate on input.  On putput, connects the newly
    *  balanced balance_box_level to the anchor box_level.
    *  If set to NULL, this parameter is ignored and Connector update
    *  is skipped.
    *
    * @param[in] hierarchy The hierarchy where the work distribution
    * data lives.
    *
    * @param[in] level_number The number of the level where the work
    * distribution data lives.
    *
    * @param[in] min_size hier::IntVector representing mimimum box size.
    *
    * @param[in] max_size hier::IntVector representing maximum box size.
    *
    * @param[in] domain_box_level Description of the domain.
    *
    * @param[in] bad_interval
    *  hier::IntVector indicating the length of an interval
    *  of cells along each side of the box where chopping
    *  the box may produce boxes with certain "bad" properties.
    *  For example, this is primiarily used to avoid generating
    *  ghost regions for patches that intersect the domain
    *  boundary in ways that may it difficult for a use to
    *  provide boundary values.  Thus, it is typically related
    *  to the maximum ghost cell width in the problem.  See
    *  hier::BoxUtilities header file for more information.
    *
    * @param[in] cut_factor
    *  hier::IntVector indicating factor for chopping
    *  each side of a box; i.e., after chopping a box,
    *  the number of cells along each direction of each
    *  piece must be an integer multiple of the corresponding
    *  entry in the cut factor vector.  For example, the
    *  cut factor may be related to the coarsen ratio between
    *  levels in the hierarchy in which case it may be used
    *  to produce boxes that can be coarsened by a certain
    *  factor if needed.  See hier::BoxUtilities header file
    *  for more information.
    *
    * @param[in] rank_group
    *  Input tbox::RankGroup indicating a set of ranks on which all boxes
    *  in the output balance_box_level will be restricted.  Some
    *  child classes may not make use of this argument.
    */
   virtual void
   loadBalanceBoxLevel(
      hier::BoxLevel& balance_box_level,
      hier::Connector * balance_to_anchor,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const hier::IntVector& min_size,
      const hier::IntVector& max_size,
      const hier::BoxLevel& domain_box_level,
      const hier::IntVector& bad_interval,
      const hier::IntVector& cut_factor, // Default v 2.x.x = 1
      const tbox::RankGroup& rank_group = tbox::RankGroup()) const = 0;

   virtual void
   setWorkloadPatchDataIndex(
      int data_id,
      int level_number = -1) = 0;

protected:
   /*!
    * Construct load balance strategy object.
    */
   LoadBalanceStrategy();

   /*!
    * @brief Write load data to log for postprocessing later.
    *
    * For development only.  Not for general use.
    *
    * @param[in] rank
    *
    * @param[in] load
    *
    * @param[in] nbox
    *
    * TODO: This method does not belong in a strategy base class and
    * should be moved elsewhere.
    */
   static void
   markLoadForPostprocessing(
      int rank,
      double load,
      int nbox);

private:
   // The following are not implemented:
   LoadBalanceStrategy(
      const LoadBalanceStrategy&);

   LoadBalanceStrategy&
   operator = (
      const LoadBalanceStrategy&);

   static int s_sequence_number;

};

}
}
#endif
