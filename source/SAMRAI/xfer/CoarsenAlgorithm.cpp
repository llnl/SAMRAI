/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Coarsening algorithm for data transfer between AMR levels
 *
 ************************************************************************/
#include "SAMRAI/xfer/CoarsenAlgorithm.h"

#include "SAMRAI/xfer/BoxGeometryVariableFillPattern.h"
#include "SAMRAI/xfer/StandardCoarsenTransactionFactory.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/VariableDatabase.h"


namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * The constructor creates a new CoarsenClasses object
 * and caches a boolean indiating whether to copy data to the
 * destination space on the coarse level before coarsening.
 *
 *************************************************************************
 */

CoarsenAlgorithm::CoarsenAlgorithm(
   const tbox::Dimension& dim,
   bool fill_coarse_data):
   d_dim(dim),
   d_coarsen_classes(std::make_shared<CoarsenClasses>()),
   d_fill_coarse_data(fill_coarse_data),
   d_schedule_created(false)
{
}

/*
 *************************************************************************
 *
 * The destructor implicitly deallocates the list data.
 *
 *************************************************************************
 */

CoarsenAlgorithm::~CoarsenAlgorithm()
{
}

/*
 *************************************************************************
 *
 * Register a coarsening operation with the coarsening algorithm.
 *
 *************************************************************************
 */

void
CoarsenAlgorithm::registerCoarsen(
   const int dst,
   const int src,
   const std::shared_ptr<hier::CoarsenOperator>& opcoarsen,
   const hier::IntVector& gcw_to_coarsen,
   const std::shared_ptr<VariableFillPattern>& var_fill_pattern)
{
   if (d_schedule_created) {
      TBOX_ERROR(
         "CoarsenAlgorithm::registerCoarsen error..."
         << "\nCannot call registerCoarsen with this coarsen algorithm"
         << "\nobject since it has already been used to create a coarsen schedule."
         << std::endl);
   }

   CoarsenClasses::Data data(d_dim);

   data.d_dst = dst;
   data.d_src = src;
   data.d_fine_bdry_reps_var = hier::VariableDatabase::getDatabase()->
      getPatchDescriptor()->getPatchDataFactory(dst)->
      fineBoundaryRepresentsVariable();
   data.d_gcw_to_coarsen = gcw_to_coarsen;
   data.d_opcoarsen = opcoarsen;
   data.d_tag = -1;
   if (var_fill_pattern) {
      data.d_var_fill_pattern = var_fill_pattern;
   } else {
      data.d_var_fill_pattern.reset(new BoxGeometryVariableFillPattern());
   }

   d_coarsen_classes->insertEquivalenceClassItem(data);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that will coarsen data from fine
 * patch level to the coarse patch level.
 *
 *************************************************************************
 */

std::shared_ptr<CoarsenSchedule>
CoarsenAlgorithm::createSchedule(
   const std::shared_ptr<hier::PatchLevel>& crse_level,
   const std::shared_ptr<hier::PatchLevel>& fine_level,
   CoarsenPatchStrategy* patch_strategy,
   const std::shared_ptr<CoarsenTransactionFactory>& transaction_factory)
{
   TBOX_ASSERT(crse_level);
   TBOX_ASSERT(fine_level);
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY2(d_dim, *crse_level, *fine_level);

   d_schedule_created = true;

   std::shared_ptr<CoarsenTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardCoarsenTransactionFactory());
   }

   return std::make_shared<CoarsenSchedule>(
             crse_level,
             fine_level,
             d_coarsen_classes,
             trans_factory,
             patch_strategy,
             d_fill_coarse_data);
}

void
CoarsenAlgorithm::resetSchedule(
   const std::shared_ptr<CoarsenSchedule>& schedule) const
{

   TBOX_ASSERT(schedule);

   if (d_coarsen_classes->classesMatch(schedule->getEquivalenceClasses())) {
      schedule->reset(d_coarsen_classes);
   } else {
      TBOX_ERROR("CoarsenAlgorithm::resetSchedule error..."
         << "\n CoarsenClasses object passed to reset routine"
         << "\n does not match that owned by existing schedule."
         << std::endl);
   }
}

/*
 *************************************************************************
 *
 * Print coarsen algorithm data to the specified output stream.
 *
 *************************************************************************
 */

void
CoarsenAlgorithm::printClassData(
   std::ostream& stream) const
{
   stream << "CoarsenAlgorithm::printClassData()" << std::endl;
   stream << "----------------------------------------" << std::endl;
   stream << "d_fill_coarse_data = " << d_fill_coarse_data << std::endl;

   d_coarsen_classes->printClassData(stream);
}

}
}
