/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Refine algorithm for data transfer between AMR levels
 *
 ************************************************************************/
#include "SAMRAI/xfer/RefineAlgorithm.h"

#include "SAMRAI/xfer/BoxGeometryVariableFillPattern.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"
#include "SAMRAI/xfer/StandardRefineTransactionFactory.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Default constructor creates a new RefineClasses object.
 *
 *************************************************************************
 */

RefineAlgorithm::RefineAlgorithm():
   d_refine_classes(std::make_shared<RefineClasses>()),
   d_schedule_created(false)
{
}

/*
 *************************************************************************
 *
 * The destructor implicitly deletes the list storage associated with
 * the refine algorithm.
 *
 *************************************************************************
 */

RefineAlgorithm::~RefineAlgorithm()
{
}

/*
 *************************************************************************
 *
 * Register a refine operation that will not require time interpolation.
 *
 *************************************************************************
 */

void
RefineAlgorithm::registerRefine(
   const int dst,
   const int src,
   const int scratch,
   const std::shared_ptr<hier::RefineOperator>& oprefine,
   const std::shared_ptr<VariableFillPattern>& var_fill_pattern,
   const std::vector<int>& work_ids)
{
   if (d_schedule_created) {
      TBOX_ERROR("RefineAlgorithm::registerRefine error..."
         << "\nCannot call registerRefine with a RefineAlgorithm"
         << "\nobject that has already been used to create a schedule."
         << std::endl);
   }

   RefineClasses::Data data;

   data.d_dst = dst;
   data.d_src = src;
   data.d_src_told = -1;
   data.d_src_tnew = -1;
   data.d_scratch = scratch;
   data.d_fine_bdry_reps_var = hier::VariableDatabase::getDatabase()->
      getPatchDescriptor()->getPatchDataFactory(dst)->
      fineBoundaryRepresentsVariable();
   data.d_time_interpolate = false;
   data.d_oprefine = oprefine;
   data.d_optime.reset();
   data.d_tag = -1;
   if (var_fill_pattern) {
      data.d_var_fill_pattern = var_fill_pattern;
   } else {
      data.d_var_fill_pattern.reset(new BoxGeometryVariableFillPattern());
   }
   data.d_work = work_ids;

   d_refine_classes->insertEquivalenceClassItem(data);
}

/*
 *************************************************************************
 *
 * Register a refine operation that will require time interpolation.
 *
 *************************************************************************
 */

void
RefineAlgorithm::registerRefine(
   const int dst,
   const int src,
   const int src_told,
   const int src_tnew,
   const int scratch,
   const std::shared_ptr<hier::RefineOperator>& oprefine,
   const std::shared_ptr<hier::TimeInterpolateOperator>& optime,
   const std::shared_ptr<VariableFillPattern>& var_fill_pattern,
   const std::vector<int>& work_ids)
{
   TBOX_ASSERT(optime);

   if (d_schedule_created) {
      TBOX_ERROR("RefineAlgorithm::registerRefine error..."
         << "\nCannot call registerRefine with a RefineAlgorithm object"
         << "\nthat has already been used to create a schedule."
         << std::endl);
   }

   RefineClasses::Data data;

   data.d_dst = dst;
   data.d_src = src;
   data.d_src_told = src_told;
   data.d_src_tnew = src_tnew;
   data.d_scratch = scratch;
   data.d_fine_bdry_reps_var = hier::VariableDatabase::getDatabase()->
      getPatchDescriptor()->getPatchDataFactory(dst)->
      fineBoundaryRepresentsVariable();
   data.d_time_interpolate = true;
   data.d_oprefine = oprefine;
   data.d_optime = optime;
   data.d_tag = -1;
   if (var_fill_pattern) {
      data.d_var_fill_pattern = var_fill_pattern;
   } else {
      data.d_var_fill_pattern.reset(new BoxGeometryVariableFillPattern());
   }
   data.d_work = work_ids;

   d_refine_classes->insertEquivalenceClassItem(data);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that will copy data from the
 * interiors of the specified level into the ghost cells and
 * interiors of the same level.
 *
 *************************************************************************
 */

std::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const std::shared_ptr<hier::PatchLevel>& level,
   RefinePatchStrategy* patch_strategy,
   const std::shared_ptr<RefineTransactionFactory>& transaction_factory)
{
   TBOX_ASSERT(level);

   d_schedule_created = true;

   std::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   std::shared_ptr<PatchLevelFullFillPattern> fill_pattern(
      std::make_shared<PatchLevelFullFillPattern>());

   return std::make_shared<RefineSchedule>(
             fill_pattern,
             level,
             level,
             d_refine_classes,
             trans_factory,
             patch_strategy);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that will copy data from the
 * interiors of the specified level into the ghost cells and
 * interiors of the same level.
 *
 *************************************************************************
 */

std::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const std::shared_ptr<PatchLevelFillPattern>& fill_pattern,
   const std::shared_ptr<hier::PatchLevel>& level,
   RefinePatchStrategy* patch_strategy,
   const std::shared_ptr<RefineTransactionFactory>& transaction_factory)
{
   TBOX_ASSERT(level);

   d_schedule_created = true;

   std::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   return std::make_shared<RefineSchedule>(
             fill_pattern,
             level,
             level,
             d_refine_classes,
             trans_factory,
             patch_strategy);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that will copy data from the
 * interiors of the source level into the ghost cell and interiors
 * of the destination level.
 *
 *************************************************************************
 */

std::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const std::shared_ptr<hier::PatchLevel>& dst_level,
   const std::shared_ptr<hier::PatchLevel>& src_level,
   RefinePatchStrategy* patch_strategy,
   bool use_time_refinement,
   const std::shared_ptr<RefineTransactionFactory>& transaction_factory)
{
   // TBOX_ERROR("Untried method!  I think this method should work, but it's never been excercised.  When code crashes here, remove this line and rerun.  If problem continues, it could well be due to excercising this code.  --BTNG");

   TBOX_ASSERT(dst_level);
   TBOX_ASSERT(src_level);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *src_level);

   d_schedule_created = true;

   std::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   std::shared_ptr<PatchLevelFullFillPattern> fill_pattern(
      std::make_shared<PatchLevelFullFillPattern>());

   return std::make_shared<RefineSchedule>(
             fill_pattern,
             dst_level,
             src_level,
             d_refine_classes,
             trans_factory,
             patch_strategy,
             use_time_refinement);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that will copy data from the
 * interiors of the source level into the ghost cell and interiors
 * of the destination level.
 *
 *************************************************************************
 */

std::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const std::shared_ptr<PatchLevelFillPattern>& fill_pattern,
   const std::shared_ptr<hier::PatchLevel>& dst_level,
   const std::shared_ptr<hier::PatchLevel>& src_level,
   RefinePatchStrategy* patch_strategy,
   bool use_time_refinement,
   const std::shared_ptr<RefineTransactionFactory>& transaction_factory)
{
   TBOX_ASSERT(dst_level);
   TBOX_ASSERT(src_level);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *src_level);

   d_schedule_created = true;

   std::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   return std::make_shared<RefineSchedule>(
             fill_pattern,
             dst_level,
             src_level,
             d_refine_classes,
             trans_factory,
             patch_strategy,
             use_time_refinement);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that copies data from the interiors
 * of the same level and coarser levels into the interior and boundary
 * cells of the given level.
 *
 *************************************************************************
 */

std::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const std::shared_ptr<hier::PatchLevel>& level,
   const int next_coarser_level,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   RefinePatchStrategy* patch_strategy,
   bool use_time_refinement,
   const std::shared_ptr<RefineTransactionFactory>& transaction_factory)
{

   // Do we all agree on the destination box_level?
   TBOX_ASSERT(level);
   TBOX_ASSERT((next_coarser_level == -1) || hierarchy);
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   if (hierarchy) {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*level, *hierarchy);
   }
#endif

   d_schedule_created = true;

   std::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   std::shared_ptr<PatchLevelFullFillPattern> fill_pattern(
      std::make_shared<PatchLevelFullFillPattern>());

   return std::make_shared<RefineSchedule>(
             fill_pattern,
             level,
             level,
             next_coarser_level,
             hierarchy,
             d_refine_classes,
             trans_factory,
             patch_strategy,
             use_time_refinement);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that copies data from the interiors
 * of the same level and coarser levels into the interior and boundary
 * cells of the given level.
 *
 *************************************************************************
 */

std::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const std::shared_ptr<PatchLevelFillPattern>& fill_pattern,
   const std::shared_ptr<hier::PatchLevel>& level,
   const int next_coarser_level,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   RefinePatchStrategy* patch_strategy,
   bool use_time_refinement,
   const std::shared_ptr<RefineTransactionFactory>& transaction_factory)
{

   // Do we all agree on the destination box_level?
   TBOX_ASSERT(level);
   TBOX_ASSERT((next_coarser_level == -1) || hierarchy);
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   if (hierarchy) {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*level, *hierarchy);
   }
#endif

   d_schedule_created = true;

   std::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   return std::make_shared<RefineSchedule>(
             fill_pattern,
             level,
             level,
             next_coarser_level,
             hierarchy,
             d_refine_classes,
             trans_factory,
             patch_strategy,
             use_time_refinement);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that copies data from the interiors
 * of the old level and coarser levels into the ghost cells and interior
 * cells of the given new level.
 *
 *************************************************************************
 */

std::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const std::shared_ptr<hier::PatchLevel>& dst_level,
   const std::shared_ptr<hier::PatchLevel>& src_level,
   const int next_coarser_level,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   RefinePatchStrategy* patch_strategy,
   bool use_time_refinement,
   const std::shared_ptr<RefineTransactionFactory>& transaction_factory)
{
   NULL_USE(use_time_refinement);

   TBOX_ASSERT(dst_level);
   TBOX_ASSERT((next_coarser_level == -1) || hierarchy);
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   if (src_level) {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *src_level);
   }
   if (hierarchy) {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *hierarchy);
   }
#endif

   // Do we all agree on the destination box_level?
   if (src_level) {
      if (next_coarser_level >= 0) {
      }
   }

   d_schedule_created = true;

   std::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   std::shared_ptr<PatchLevelFullFillPattern> fill_pattern(
      std::make_shared<PatchLevelFullFillPattern>());

   return std::make_shared<RefineSchedule>(
             fill_pattern,
             dst_level,
             src_level,
             next_coarser_level,
             hierarchy,
             d_refine_classes,
             trans_factory,
             patch_strategy,
             false);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that copies data from the interiors
 * of the old level and coarser levels into the ghost cells and interior
 * cells of the given new level.
 *
 *************************************************************************
 */

std::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const std::shared_ptr<PatchLevelFillPattern>& fill_pattern,
   const std::shared_ptr<hier::PatchLevel>& dst_level,
   const std::shared_ptr<hier::PatchLevel>& src_level,
   const int next_coarser_level,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   RefinePatchStrategy* patch_strategy,
   bool use_time_refinement,
   const std::shared_ptr<RefineTransactionFactory>& transaction_factory)
{
   NULL_USE(use_time_refinement);

   TBOX_ASSERT(dst_level);
   TBOX_ASSERT((next_coarser_level == -1) || hierarchy);
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   if (src_level) {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *src_level);
   }
   if (hierarchy) {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *hierarchy);
   }
#endif

   // Do we all agree on the destination box_level?
   if (src_level) {
      if (next_coarser_level >= 0) {
      }
   }

   d_schedule_created = true;

   std::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   return std::make_shared<RefineSchedule>(
             fill_pattern,
             dst_level,
             src_level,
             next_coarser_level,
             hierarchy,
             d_refine_classes,
             trans_factory,
             patch_strategy,
             false);
}

/*
 **************************************************************************
 *
 * Reconfigure refine schedule to perform operations in this algorithm.
 *
 **************************************************************************
 */

bool
RefineAlgorithm::checkConsistency(
   const std::shared_ptr<RefineSchedule>& schedule) const
{
   TBOX_ASSERT(schedule);
   return d_refine_classes->classesMatch(schedule->getEquivalenceClasses());
}

void RefineAlgorithm::resetSchedule(
   const std::shared_ptr<RefineSchedule>& schedule) const
{
   TBOX_ASSERT(schedule);
   if (d_refine_classes->classesMatch(schedule->getEquivalenceClasses())) {
      schedule->reset(d_refine_classes);
   } else {
      TBOX_ERROR("RefineAlgorithm::resetSchedule error..."
         << "\n Items in RefineClasses object passed to reset"
         << "\n routine does not match that in existing schedule."
         << std::endl);
   }
}

/*
 *************************************************************************
 *
 * Print refine algorithm data to the specified output stream.
 *
 *************************************************************************
 */

void
RefineAlgorithm::printClassData(
   std::ostream& stream) const
{
   stream << "RefineAlgorithm::printClassData()" << std::endl;
   stream << "----------------------------------------" << std::endl;
   d_refine_classes->printClassData(stream);
}

}
}
