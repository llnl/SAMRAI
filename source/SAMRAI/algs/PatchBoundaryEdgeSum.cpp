/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Routines for summing edge data at patch boundaries
 *
 ************************************************************************/
#include "SAMRAI/algs/PatchBoundaryEdgeSum.h"

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/EdgeDataFactory.h"
#include "SAMRAI/pdat/OuteredgeData.h"
#include "SAMRAI/algs/OuteredgeSumTransactionFactory.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/tbox/Collectives.h"

namespace SAMRAI {
namespace algs {

/*
 *************************************************************************
 *
 * Initialize the static data members.
 *
 *************************************************************************
 */

int PatchBoundaryEdgeSum::s_instance_counter = 0;

std::vector<std::vector<int> > PatchBoundaryEdgeSum::s_oedge_src_id_array =
   std::vector<std::vector<int> >(0);
std::vector<std::vector<int> > PatchBoundaryEdgeSum::s_oedge_dst_id_array =
   std::vector<std::vector<int> >(0);

/*
 *************************************************************************
 *
 * Constructor patch boundary edge sum objects initializes data members
 * to default (undefined) states.
 *
 *************************************************************************
 */

PatchBoundaryEdgeSum::PatchBoundaryEdgeSum(
   const std::string& object_name):
   d_setup_called(false),
   d_num_reg_sum(0),
   d_sum_transaction_factory(std::make_shared<OuteredgeSumTransactionFactory>())
{
   TBOX_ASSERT(!object_name.empty());

   d_object_name = object_name;

   ++s_instance_counter;
}

/*
 *************************************************************************
 *
 * Destructor removes temporary outeredge patch data ids from
 * variable database, if defined.
 *
 *************************************************************************
 */

PatchBoundaryEdgeSum::~PatchBoundaryEdgeSum()
{

   --s_instance_counter;
   if (s_instance_counter == 0) {
      const int arr_length_depth =
         static_cast<int>(s_oedge_src_id_array.size());
      for (int id = 0; id < arr_length_depth; ++id) {
         const int arr_length_nvar =
            static_cast<int>(s_oedge_src_id_array[id].size());

         for (int iv = 0; iv < arr_length_nvar; ++iv) {

            if (s_oedge_src_id_array[id][iv] >= 0) {
               hier::VariableDatabase::getDatabase()->
               removeInternalSAMRAIVariablePatchDataIndex(
                  s_oedge_src_id_array[id][iv]);
            }
            if (s_oedge_dst_id_array[id][iv] >= 0) {
               hier::VariableDatabase::getDatabase()->
               removeInternalSAMRAIVariablePatchDataIndex(
                  s_oedge_dst_id_array[id][iv]);
            }

            s_oedge_src_id_array[id].resize(0);
            s_oedge_dst_id_array[id].resize(0);

         }
      }

      s_oedge_src_id_array.resize(0);
      s_oedge_dst_id_array.resize(0);
   }

}

/*
 *************************************************************************
 *
 * Register edge patch data index for summation.
 *
 *************************************************************************
 */

void
PatchBoundaryEdgeSum::registerSum(
   int edge_data_id)
{
   if (d_setup_called) {
      TBOX_ERROR("PatchBoundaryEdgeSum::register error..."
         << "\nobject named " << d_object_name
         << "\nCannot call registerSum with this PatchBoundaryEdgeSum"
         << "\nobject since it has already been used to create communication"
         << "\nschedules; i.e., setupSum() has been called."
         << std::endl);
   }

   if (edge_data_id < 0) {
      TBOX_ERROR("PatchBoundaryEdgeSum register error..."
         << "\nobject named " << d_object_name
         << "\n edge_data_id = " << edge_data_id
         << " is an invalid patch data identifier." << std::endl);
   }

   hier::VariableDatabase* var_db = hier::VariableDatabase::getDatabase();

   std::shared_ptr<pdat::EdgeDataFactory<double> > edge_factory(
      SAMRAI_SHARED_PTR_CAST<pdat::EdgeDataFactory<double>, hier::PatchDataFactory>(
         var_db->getPatchDescriptor()->getPatchDataFactory(edge_data_id)));

   TBOX_ASSERT(edge_factory);

   const tbox::Dimension& dim(edge_factory->getDim());

   static std::string tmp_oedge_src_variable_name(
      "PatchBoundaryEdgeSum__internal-oedge-src");
   static std::string tmp_oedge_dst_variable_name(
      "PatchBoundaryEdgeSum__internal-oedge-dst");

   const int reg_sum_id = d_num_reg_sum;

   ++d_num_reg_sum;

   d_user_edge_data_id.resize(d_num_reg_sum);
   d_user_edge_data_id[reg_sum_id] = ID_UNDEFINED;
   d_user_edge_depth.resize(d_num_reg_sum);
   d_user_edge_depth[reg_sum_id] = ID_UNDEFINED;
   d_tmp_oedge_src_variable.resize(d_num_reg_sum);
   d_tmp_oedge_dst_variable.resize(d_num_reg_sum);
   d_oedge_src_id.resize(d_num_reg_sum);
   d_oedge_src_id[reg_sum_id] = ID_UNDEFINED;
   d_oedge_dst_id.resize(d_num_reg_sum);
   d_oedge_dst_id[reg_sum_id] = ID_UNDEFINED;

   const int data_depth = edge_factory->getDepth();
   const int array_by_depth_size = data_depth + 1;

   if (static_cast<int>(d_num_registered_data_by_depth.size()) <
       array_by_depth_size) {
      const int old_size =
         static_cast<int>(d_num_registered_data_by_depth.size());
      const int new_size = array_by_depth_size;

      d_num_registered_data_by_depth.resize(new_size);
      for (int i = old_size; i < new_size; ++i) {
         d_num_registered_data_by_depth[i] = 0;
      }
   }

   const int data_depth_id = d_num_registered_data_by_depth[data_depth];
   const int num_data_at_depth = data_depth_id + 1;

   if (static_cast<int>(s_oedge_src_id_array.size()) < array_by_depth_size) {
      s_oedge_src_id_array.resize(array_by_depth_size);
      s_oedge_dst_id_array.resize(array_by_depth_size);
   }

   if (static_cast<int>(s_oedge_src_id_array[data_depth].size()) <
       num_data_at_depth) {
      const int old_size =
         static_cast<int>(s_oedge_src_id_array[data_depth].size());
      const int new_size = num_data_at_depth;

      s_oedge_src_id_array[data_depth].resize(new_size);
      s_oedge_dst_id_array[data_depth].resize(new_size);
      for (int i = old_size; i < new_size; ++i) {
         s_oedge_src_id_array[data_depth][i] = ID_UNDEFINED;
         s_oedge_dst_id_array[data_depth][i] = ID_UNDEFINED;
      }
   }

   std::string var_suffix = tbox::Utilities::intToString(data_depth_id, 4)
      + "__depth=" + tbox::Utilities::intToString(data_depth);

   std::string toedge_src_var_name = tmp_oedge_src_variable_name + var_suffix;
   d_tmp_oedge_src_variable[reg_sum_id] = var_db->getVariable(
         toedge_src_var_name);
   if (!d_tmp_oedge_src_variable[reg_sum_id]) {
      d_tmp_oedge_src_variable[reg_sum_id].reset(
         new pdat::OuteredgeVariable<double>(dim,
            toedge_src_var_name,
            data_depth));
   }

   std::string toedge_dst_var_name = tmp_oedge_dst_variable_name + var_suffix;
   d_tmp_oedge_dst_variable[reg_sum_id] = var_db->getVariable(
         toedge_dst_var_name);
   if (!d_tmp_oedge_dst_variable[reg_sum_id]) {
      d_tmp_oedge_dst_variable[reg_sum_id].reset(
         new pdat::OuteredgeVariable<double>(dim,
            toedge_dst_var_name,
            data_depth));
   }

   if (s_oedge_src_id_array[data_depth][data_depth_id] < 0) {
      s_oedge_src_id_array[data_depth][data_depth_id] =
         var_db->registerInternalSAMRAIVariable(
            d_tmp_oedge_src_variable[reg_sum_id],
            hier::IntVector::getZero(dim));
   }
   if (s_oedge_dst_id_array[data_depth][data_depth_id] < 0) {
      s_oedge_dst_id_array[data_depth][data_depth_id] =
         var_db->registerInternalSAMRAIVariable(
            d_tmp_oedge_dst_variable[reg_sum_id],
            hier::IntVector::getZero(dim));
   }

   d_user_edge_data_id[reg_sum_id] = edge_data_id;
   d_user_edge_depth[reg_sum_id] = data_depth;

   d_num_registered_data_by_depth[data_depth] = num_data_at_depth;

   d_oedge_src_id[reg_sum_id] =
      s_oedge_src_id_array[data_depth][data_depth_id];
   d_oedge_dst_id[reg_sum_id] =
      s_oedge_dst_id_array[data_depth][data_depth_id];

   d_oedge_src_data_set.setFlag(d_oedge_src_id[reg_sum_id]);
   d_oedge_dst_data_set.setFlag(d_oedge_dst_id[reg_sum_id]);

}

/*
 *************************************************************************
 *
 * Set up schedule to sum edge data around patch boundaries
 * on a single level.
 *
 *************************************************************************
 */

void
PatchBoundaryEdgeSum::setupSum(
   const std::shared_ptr<hier::PatchLevel>& level)
{
   TBOX_ASSERT(level);

   d_setup_called = true;

   d_level = level;

   // Communication algorithm for summing outeredge values on a level
   xfer::RefineAlgorithm single_level_sum_algorithm;

   for (int i = 0; i < d_num_reg_sum; ++i) {
      single_level_sum_algorithm.registerRefine(d_oedge_dst_id[i],  // dst data
         d_oedge_src_id[i],                                         // src data
         d_oedge_dst_id[i],                                         // scratch data
         std::shared_ptr<hier::RefineOperator>());
   }

   d_single_level_sum_schedule =
      single_level_sum_algorithm.createSchedule(
         d_level,
         0,
         d_sum_transaction_factory);

}

/*
 *************************************************************************
 *
 * Perform patch boundary edge sum across single level or multiple
 * hierarchy levels depending on how object was initialized.
 *
 *************************************************************************
 */

void
PatchBoundaryEdgeSum::computeSum() const
{
   d_level->allocatePatchData(d_oedge_src_data_set);
   d_level->allocatePatchData(d_oedge_dst_data_set);

   doLevelSum(d_level);

   d_level->deallocatePatchData(d_oedge_src_data_set);
   d_level->deallocatePatchData(d_oedge_dst_data_set);

}

/*
 *************************************************************************
 *
 * Private member function that performs edge sum across single level.
 *
 * 1. Copy edge data to local outeredge data.
 * 2. Transfer and sum outeredge data.
 * 3. Copy outeredge data back to edge data.
 *
 *************************************************************************
 */

void
PatchBoundaryEdgeSum::doLevelSum(
   const std::shared_ptr<hier::PatchLevel>& level) const
{
   TBOX_ASSERT(level);

   for (hier::PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
      const std::shared_ptr<hier::Patch>& patch = *ip;

      int array_size = static_cast<int>(d_user_edge_data_id.size());
      for (int i = 0; i < array_size; ++i) {
         std::shared_ptr<pdat::EdgeData<double> > edge_data(
            SAMRAI_SHARED_PTR_CAST<pdat::EdgeData<double>, hier::PatchData>(
               patch->getPatchData(d_user_edge_data_id[i])));
         std::shared_ptr<pdat::OuteredgeData<double> > oedge_data(
            SAMRAI_SHARED_PTR_CAST<pdat::OuteredgeData<double>, hier::PatchData>(
               patch->getPatchData(d_oedge_src_id[i])));

         TBOX_ASSERT(edge_data);
         TBOX_ASSERT(oedge_data);

         oedge_data->copy(*edge_data);
      }
   }
#if defined(HAVE_RAJA)
   tbox::parallel_synchronize();
#endif

   d_single_level_sum_schedule->fillData(0.0, false);

   for (hier::PatchLevel::iterator ip2(level->begin());
        ip2 != level->end(); ++ip2) {
      const std::shared_ptr<hier::Patch>& patch = *ip2;

      int array_size = static_cast<int>(d_user_edge_data_id.size());
      for (int i = 0; i < array_size; ++i) {
         std::shared_ptr<pdat::EdgeData<double> > edge_data(
            SAMRAI_SHARED_PTR_CAST<pdat::EdgeData<double>, hier::PatchData>(
               patch->getPatchData(d_user_edge_data_id[i])));
         std::shared_ptr<pdat::OuteredgeData<double> > oedge_data(
            SAMRAI_SHARED_PTR_CAST<pdat::OuteredgeData<double>, hier::PatchData>(
               patch->getPatchData(d_oedge_dst_id[i])));

         TBOX_ASSERT(edge_data);
         TBOX_ASSERT(oedge_data);

         oedge_data->copy2(*edge_data);

      }
   }
#if defined(HAVE_RAJA)
   tbox::parallel_synchronize();
#endif
}

}
}
