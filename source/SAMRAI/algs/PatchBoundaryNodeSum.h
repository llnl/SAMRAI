/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Routines for summing node data at patch boundaries
 *
 ************************************************************************/

#ifndef included_algs_PatchBoundaryNodeSum
#define included_algs_PatchBoundaryNodeSum

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/OuternodeVariable.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/RefineTransactionFactory.h"
#include "SAMRAI/tbox/Utilities.h"

#include <string>
#include <memory>

namespace SAMRAI {
namespace algs {

/*!
 *  @brief Class PatchBoundaryNodeSum provides operations for summing node data
 *  values at nodes that are shared by multiple patches on a single level or
 *  across multiple hierarchy levels.
 *
 *  Usage of a patch boundry node sum involves the following sequence of steps:
 *
 *  -# Construct a patch boundry node sum object.  For example,
 *     \verbatim
 *         PatchBoundaryNodeSum my_node_sum("My Node Sum");
 *     \endverbatim
 *  -# Register node data quantities to sum.  For example,
 *     \verbatim
 *         my_node_sum.registerSum(node_data_id1);
 *         my_node_sum.registerSum(node_data_id2);
 *         etc...
 *     \endverbatim
 *  -# Setup the sum operations for either single level or a range of levels
 *     in a patch hierarchy.  For example,
 *     \verbatim
 *         my_node_sum.setupSum(level);    // single level
 *         -- or --
 *         my_node_sum.setupSum(hierarchy, coarsest_ln, finest_ln);  // multiple levels
 *     \endverbatim
 *  -# Execute the sum operation.  For example,
 *     \verbatim
 *         my_node_sum.computeSum()
 *     \endverbatim
 *
 *  The result of these operations is that each node patch data value
 *  associated with the registered ids at patch boundaries, on either the
 *  single level or range of hierarchy levels, is replaced by the sum of all
 *  data values at the node.
 *
 *  Note that only one of the setupSum() functions may be called once a
 *  PatchBoundaryNodeSum object is created.
 */

class PatchBoundaryNodeSum
{
public:
   /*!
    *  @brief Static function used to predetermine number of patch data
    *         slots ahared among all PatchBoundaryNodeSum
    *         objects (i.e., static members).  To get a correct count,
    *         this routine should only be called once.
    *
    *  @return integer number of internal patch data slots required
    *          to perform sum.
    *  @param max_variables_to_register integer value indicating
    *          maximum number of patch data ids that will be registered
    *          with node sum objects.TO
    */
   static int
   getNumSharedPatchDataSlots(
      int max_variables_to_register)
   {
      // node boundary sum requires two internal outernode variables
      // (source and destination) for each registered variable.
      return 2 * max_variables_to_register;
   }

   /*!
    *  @brief Static function used to predetermine number of patch data
    *         slots unique to each PatchBoundaryNodeSum
    *         object (i.e., non-static members).  To get a correct count,
    *         this routine should be called exactly once for each object
    *         that will be constructed.
    *
    *  @return integer number of internal patch data slots required
    *          to perform sum.
    *  @param max_variables_to_register integer value indicating
    *          maximum number of patch data ids that will be registered
    *          with node sum objects.
    */
   static int
   getNumUniquePatchDataSlots(
      int max_variables_to_register)
   {
      NULL_USE(max_variables_to_register);
      // all patch data slots used by node boundary sum are static
      // and shared among all objects.
      return 0;
   }

   /*!
    *  @brief Constructor initializes object to default (mostly undefined)
    *  state.
    *
    *  @param object_name const std::string reference for name of object used
    *  in error reporting.
    *
    *  @pre !object_name.empty()
    */
   explicit PatchBoundaryNodeSum(
      const std::string& object_name);

   /*!
    *  @brief Destructor for the schedule releases all internal storage.
    */
   ~PatchBoundaryNodeSum();

   /*!
    *  @brief Register node data with given patch data identifier for summing.
    *
    *  @param node_data_id  integer patch data index for node data to sum
    *
    *  @pre !d_setup_called
    *  @pre node_data_id >= 0
    *  @pre hier::VariableDatabase::getDatabase()->getPatchDescriptor()->getPatchDataFactory(node_data_id) is actually a std::shared_ptr<pdat::NodeDataFactory<double> >
    */
   void
   registerSum(
      int node_data_id);

   /*!
    *  @brief Set up summation operations for node data across shared nodes
    *         on a single level.
    *
    *  If the other setupSum() function for a range of hierarchy levels has
    *  been called previously for this object, an error will result.
    *
    *  @param level         pointer to level on which to perform node sum
    *
    *  @pre level
    *  @pre !d_hierarchy_setup_called
    */
   void
   setupSum(
      const std::shared_ptr<hier::PatchLevel>& level);

   /*!
    *  @brief Set up for summation operations for node data at shared nodes
    *         across a range of hierarchy levels.
    *
    *  If the other setupSum() function for a single level has been called
    *  previously for this object, an error will result.
    *
    *  @param hierarchy      pointer to hierarchy on which to perform node sum
    *  @param coarsest_level coarsest level number for node sum
    *  @param finest_level   finest level number for node sum
    *
    *  @pre hierarchy
    *  @pre (coarsest_level >= 0) && (finest_level >= coarsest_level) &&
    *       (finest_level <= hierarchy->getFinestLevelNumber())
    *  @pre !d_hierarchy_setup_called
    */
   void
   setupSum(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level);

   /*!
    *  @brief Compute sum of node values at each shared node and replace
    *         each such node value with the corresponding sum.
    *
    *  At the end of this method, all values at shared node locations on
    *  patch boundaries (on levels indicated by the call to one of the
    *  setupSum() routines) will have the same value.
    *
    *  When the setupSum() method taking a range of patch levels in a
    *  hierarchy is called, this method will compute the sum of nodal
    *  quantities at all the specified patch boundaries.  For nodes at a
    *  coarse-fine boundary, nodal sums will only be performed where the
    *  coarse and fine nodes overlap.  A node on a fine level that is not
    *  also a node on the next coarser level (a so-called "hanging node")
    *  will not be summed.
    *
    *  The boolean "fill_hanging_nodes" argument specifies whether the
    *  the hanging nodes should be filled using linearly interpolated values
    *  from neighboring non-hanging nodes (i.e. those that overlap nodes on
    *  a coarse level). The correct steps required to deal with hanging
    *  nodes is algorithm dependent so, if left unspecified, values at the
    *  hanging nodes will not be adjusted.  However, because many algorithms
    *  average hanging nodes we provide the capability to do it here.  Note
    *  that the hanging node interpolation provided does not take into
    *  consideration the spatial location of the nodes.  So the interpolation
    *  may not be correct for coordinate systems other than standard Cartesian
    *  grid geometry.
    *
    *  @param fill_hanging_nodes Optional boolean value specifying whether
    *         hanging node values should be set to values interpolated from
    *         neighboring non-hanging node values.  The default is false.
    */
   void
   computeSum(
      const bool fill_hanging_nodes = false) const;

   /*!
    * @brief Returns the object name.
    *
    * @return The object name.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

private:
   /*!
    * @brief Perform node sum across single level.
    *
    * Called from computeSum().
    *
    * @pre level
    */
   void
   doLevelSum(
      const std::shared_ptr<hier::PatchLevel>& level) const;

   /*!
    * @ Sum node data on a coarse-fine boundary
    *
    * A fine level and a coarse level are given as arguments, with the
    * coarse level being a coarsened representation of the fine level.
    * This method modifies node data on the coarse-fine boundary of the
    * fine level by summing the existing node data values with outernode
    * data values from the coarse level.
    *
    * The data to modify are specified by the node_data_id and onode_data_id
    * arrays.  Each entry in node_data_id identifies data that will be
    * summed with data identified by the corresponding entry in onode_data_id.
    *
    * If the boolean fill_hanging_nodes is false, only data on the nodes
    * coincident between the fine and coarse levels will be modified.  If true,
    * linear interpolation will be used to fill the remaining fine nodes
    * on the coarse-fine boundary.  See documentation of method computeSum()
    * for more information on this argument.
    *
    * @param fine_level            Level where data will be modified
    * @param coarsened_fine_level  Coarsened version of fine_level
    * @param node_data_id   Vector of data ids specifying data to modify
    * @param onode_data_id  Vector of data ids specifying data to use in sums
    * @param fill_hanging_nodes    Tells whether to fill fine data on
    *                              intermediate fine nodes.
    *
    * @pre fine_level
    * @pre coarsened_fine_level
    * @pre fine_level->getDim() == coarsened_fine_level->getDim()
    * @pre node_data_id.size() == onode_data_id.size()
    * @pre for each member, i, of node_data_id fine_level->checkAllocated(i)
    * @pre for each member, i, of onode_data_id coarsened_fine_level->checkAllocated(i)
    */
   void
   doLocalCoarseFineBoundarySum(
      const std::shared_ptr<hier::PatchLevel>& fine_level,
      const std::shared_ptr<hier::PatchLevel>& coarsened_fine_level,
      const std::vector<int>& node_data_id,
      const std::vector<int>& onode_data_id,
      bool fill_hanging_nodes) const;

   /*
    * @brief Copy node data to outernode data on all patches of level
    *
    * Data specified in the node_data_id array will be copied to data
    * specified by the onode_data_id array on all patches.
    *
    * @param level
    * @param node_data_id   Vector of data ids for NodeData source
    * @param onode_data_id  Vector of data ids for OuternodeData destination
    */
   void
   copyNodeToOuternodeOnLevel(
      const std::shared_ptr<hier::PatchLevel>& level,
      const std::vector<int>& node_data_id,
      const std::vector<int>& onode_data_id) const;

   /*
    * @brief Copy outernode data to node data on all patches of level
    *
    * Data specified in the onode_data_id array will be copied to data
    * specified by the node_data_id array on all patches.
    *
    * @param level
    * @param onode_data_id Vector of data ids for OuternodeData source
    * @param node_data_id  Vector of data ids for NodeData destination
    */
   void
   copyOuternodeToNodeOnLevel(
      const std::shared_ptr<hier::PatchLevel>& level,
      const std::vector<int>& onode_data_id,
      const std::vector<int>& node_data_id) const;

   /*
    * Static members for managing shared temporary data among multiple
    * PatchBoundaryNodeSum objects.
    */
   static int s_instance_counter;
   // These arrays are indexed [data depth][number of variables with depth]
   static std::vector<std::vector<int> > s_onode_src_id_array;
   static std::vector<std::vector<int> > s_onode_dst_id_array;

   enum PATCH_BDRY_NODE_SUM_DATA_ID { ID_UNDEFINED = -1 };

   std::string d_object_name;
   bool d_setup_called;

   int d_num_reg_sum;

   // These arrays are indexed [variable registration sequence number]
   std::vector<int> d_user_node_data_id;
   std::vector<int> d_user_node_depth;

   // These arrays are indexed [data depth]
   std::vector<int> d_num_registered_data_by_depth;

   /*
    * Node-centered variables and patch data indices used as internal work
    * quantities.
    */
   // These arrays are indexed [variable registration sequence number]
   std::vector<std::shared_ptr<hier::Variable> > d_tmp_onode_src_variable;
   std::vector<std::shared_ptr<hier::Variable> > d_tmp_onode_dst_variable;

   // These arrays are indexed [variable registration sequence number]
   std::vector<int> d_onode_src_id;
   std::vector<int> d_onode_dst_id;

   /*
    * Sets of indices for temporary variables to expedite allocation and
    * deallocation.
    */
   hier::ComponentSelector d_onode_src_data_set;
   hier::ComponentSelector d_onode_dst_data_set;

   std::shared_ptr<hier::PatchLevel> d_level;

   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;
   int d_coarsest_level;
   int d_finest_level;

   bool d_level_setup_called;
   bool d_hierarchy_setup_called;

   std::shared_ptr<xfer::RefineTransactionFactory> d_sum_transaction_factory;

   std::vector<std::shared_ptr<xfer::RefineSchedule> >
   d_single_level_sum_schedule;
   std::vector<std::shared_ptr<xfer::RefineSchedule> >
   d_cfbdry_copy_schedule;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule> >
   d_sync_coarsen_schedule;

   // A coarsened version of each fine level.
   std::vector<std::shared_ptr<hier::PatchLevel> > d_cfbdry_tmp_level;

   std::vector<std::shared_ptr<hier::CoarseFineBoundary> >
   d_coarse_fine_boundary;

};

}
}

#endif
