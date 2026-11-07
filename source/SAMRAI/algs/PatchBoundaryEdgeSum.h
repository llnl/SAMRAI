/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Routines for summing edge data at patch boundaries
 *
 ************************************************************************/

#ifndef included_algs_PatchBoundaryEdgeSum
#define included_algs_PatchBoundaryEdgeSum

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/pdat/OuteredgeVariable.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/RefineTransactionFactory.h"
#include "SAMRAI/tbox/Utilities.h"

#include <string>
#include <vector>
#include <memory>

namespace SAMRAI {
namespace algs {

/*!
 *  @brief Class PatchBoundaryEdgeSum provides operations summing edge data
 *  values at edges that are shared by multiple patches on a single level.
 *  Note that this utility only works on a SINGLE patch level, not on a
 *  multiple levels in an AMR patch hierarchy like the PatchBoundaryNodeSum
 *  class.   Unlike node data, edge data at coarse-fine boundaries are not
 *  co-located, so the sum operation is not clearly defined.
 *
 *  Usage of a patch boundry edge sum involves the following sequence of steps:
 *
 *  -# Construct a patch boundry edge sum object.  For example,
 *     \verbatim
 *         PatchBoundaryEdgeSum my_edge_sum("My Edge Sum");
 *     \endverbatim
 *  -# Register edge data quantities to sum.  For example,
 *     \verbatim
 *         my_edge_sum.registerSum(edge_data_id1);
 *         my_edge_sum.registerSum(edge_data_id2);
 *         etc...
 *     \endverbatim
 *  -# Setup the sum operations for a single level.  For example,
 *     \verbatim
 *         my_edge_sum.setupSum(level);
 *     \endverbatim
 *  -# Execute the sum operation.  For example,
 *     \verbatim
 *         my_edge_sum.computeSum()
 *     \endverbatim
 *
 *  The result of these operations is that each edge patch data value
 *  associated with the registered ids at patch boundaries on the level is
 *  replaced by the sum of all data values at the edge.
 */

class PatchBoundaryEdgeSum
{
public:
   /*!
    *  @brief Static function used to predetermine number of patch data
    *         slots ahared among all PatchBoundaryEdgeSum
    *         objects (i.e., static members).  To get a correct count,
    *         this routine should only be called once.
    *
    *  @return integer number of internal patch data slots required
    *          to perform sum.
    *  @param max_variables_to_register integer value indicating
    *          maximum number of patch data ids that will be registered
    *          with edge sum objects.
    */
   static int
   getNumSharedPatchDataSlots(
      int max_variables_to_register)
   {
      // edge boundary sum requires two internal outeredge variables
      // (source and destination) for each registered variable.
      return 2 * max_variables_to_register;
   }

   /*!
    *  @brief Static function used to predetermine number of patch data
    *         slots unique to each PatchBoundaryEdgeSum
    *         object (i.e., non-static members).  To get a correct count,
    *         this routine should be called exactly once for each object
    *         that will be constructed.
    *
    *  @return integer number of internal patch data slots required
    *          to perform sum.
    *  @param max_variables_to_register integer value indicating
    *          maximum number of patch data ids that will be registered
    *          with edge sum objects.
    */
   static int
   getNumUniquePatchDataSlots(
      int max_variables_to_register)
   {
      NULL_USE(max_variables_to_register);
      // all patch data slots used by edge boundary sum are static
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
   explicit PatchBoundaryEdgeSum(
      const std::string& object_name);

   /*!
    *  @brief Destructor for the schedule releases all internal storage.
    */
   ~PatchBoundaryEdgeSum();

   /*!
    *  @brief Register edge data with given patch data identifier for summing.
    *
    *  @param edge_data_id  integer patch data index for edge data to sum
    *
    *  @pre !d_setup_called
    *  @pre edge_data_id >= 0
    *  @pre hier::VariableDatabase::getDatabase()->getPatchDescriptor()->getPatchDataFactory(edge_data_id) is actually a std::shared_ptr<pdat::EdgeDataFactory<double> >
    */
   void
   registerSum(
      int edge_data_id);

   /*!
    *  @brief Set up summation operations for edge data across shared edges
    *         on a single level.
    *
    *  @param level         pointer to level on which to perform edge sum
    *
    *  @pre level
    */
   void
   setupSum(
      const std::shared_ptr<hier::PatchLevel>& level);

   /*!
    *  @brief Compute sum of edge values at each shared edge and replace
    *         each such edge value with the corresponding sum.
    *
    *  At the end of this method, all values at shared edge locations on
    *  patch boundaries will have the same value.
    */
   void
   computeSum() const;

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
   /*
    * Private member function to perform edge sum across single level --
    * called from computeSum()
    */
   void
   doLevelSum(
      const std::shared_ptr<hier::PatchLevel>& level) const;

   /*
    * Static members for managing shared temporary data among multiple
    * PatchBoundaryEdgeSum objects.
    */
   static int s_instance_counter;
   // These arrays are indexed [data depth][number of variables with depth]
   static std::vector<std::vector<int> > s_oedge_src_id_array;
   static std::vector<std::vector<int> > s_oedge_dst_id_array;

   enum PATCH_BDRY_EDGE_SUM_DATA_ID { ID_UNDEFINED = -1 };

   std::string d_object_name;
   bool d_setup_called;

   int d_num_reg_sum;

   // These arrays are indexed [variable registration sequence number]
   std::vector<int> d_user_edge_data_id;
   std::vector<int> d_user_edge_depth;

   // These arrays are indexed [data depth]
   std::vector<int> d_num_registered_data_by_depth;

   /*
    * Edge-centered variables and patch data indices used as internal work
    * quantities.
    */
   // These arrays are indexed [variable registration sequence number]
   std::vector<std::shared_ptr<hier::Variable> > d_tmp_oedge_src_variable;
   std::vector<std::shared_ptr<hier::Variable> > d_tmp_oedge_dst_variable;

   // These arrays are indexed [variable registration sequence number]
   std::vector<int> d_oedge_src_id;
   std::vector<int> d_oedge_dst_id;

   /*
    * Sets of indices for temporary variables to expedite allocation and
    * deallocation.
    */
   hier::ComponentSelector d_oedge_src_data_set;
   hier::ComponentSelector d_oedge_dst_data_set;

   std::shared_ptr<hier::PatchLevel> d_level;

   std::shared_ptr<xfer::RefineTransactionFactory> d_sum_transaction_factory;

   std::shared_ptr<xfer::RefineSchedule> d_single_level_sum_schedule;

};

}
}

#endif
