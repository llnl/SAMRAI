/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Manager class for patch data communication tests.
 *
 ************************************************************************/

#include "CommTester.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/xfer/CompositeBoundaryAlgorithm.h"

namespace SAMRAI {


/*
 *************************************************************************
 *
 * The constructor initializes object state.  The destructor is empty.
 *
 *************************************************************************
 */

CommTester::CommTester(
   const std::string& object_name,
   const tbox::Dimension& dim,
   std::shared_ptr<tbox::Database> main_input_db,
   PatchDataTestStrategy* data_test,
   bool do_refine,
   bool do_coarsen,
   const std::string& refine_option):
   RefinePatchStrategy(),
   CoarsenPatchStrategy(),
   d_dim(dim),
   d_fill_source_algorithm(),
   d_refine_algorithm(),
   d_coarsen_algorithm(dim),
   d_reset_refine_algorithm(),
   d_reset_coarsen_algorithm(dim)
{
   NULL_USE(main_input_db);

   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(main_input_db);
   TBOX_ASSERT(data_test != 0);

   d_object_name = object_name;
   d_data_test_strategy = data_test;

   d_fake_time = 0.0;

   d_fake_cycle = 0;

   d_is_reset = false;

   d_do_refine = do_refine;
   d_do_coarsen = false;
   if (!do_refine) {
      d_do_coarsen = do_coarsen;
   }

   d_refine_option = refine_option;
   if (!((d_refine_option == "INTERIOR_FROM_SAME_LEVEL")
         || (d_refine_option == "INTERIOR_FROM_COARSER_LEVEL"))) {
      TBOX_ERROR(object_name << " input error: illegal refine_option = "
                             << d_refine_option << std::endl);
   }

   d_patch_data_components.clrAllFlags();
   d_fill_source_schedule.resize(0);
   d_refine_schedule.resize(0);
   d_coarsen_schedule.resize(0);

   d_source =
      hier::VariableDatabase::getDatabase()->getContext("SOURCE");
   d_destination =
      hier::VariableDatabase::getDatabase()->getContext("DESTINATION");
   d_refine_scratch =
      hier::VariableDatabase::getDatabase()->getContext("REFINE_SCRATCH");

   d_reset_source =
      hier::VariableDatabase::getDatabase()->getContext("SOURCE");
   d_reset_destination =
      hier::VariableDatabase::getDatabase()->getContext("DESTINATION");
   d_reset_refine_scratch =
      hier::VariableDatabase::getDatabase()->getContext("REFINE_SCRATCH");

   d_data_test_strategy->registerVariables(this);

}

CommTester::~CommTester()
{

}

/*
 *************************************************************************
 *
 * Add variable with associated attributes to set of test variables.
 *
 *************************************************************************
 */

void CommTester::registerVariable(
   const std::shared_ptr<hier::Variable> src_variable,
   const std::shared_ptr<hier::Variable> dst_variable,
   const hier::IntVector& src_ghosts,
   const hier::IntVector& dst_ghosts,
   const std::shared_ptr<hier::BaseGridGeometry> xfer_geom,
   const std::string& operator_name)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(src_ghosts, dst_ghosts);

   TBOX_ASSERT(src_variable);
   TBOX_ASSERT(dst_variable);
   TBOX_ASSERT(xfer_geom);
   TBOX_ASSERT(!operator_name.empty());

   const tbox::Dimension dim(src_ghosts.getDim());

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   int src_id = variable_db->registerVariableAndContext(src_variable,
         d_source,
         src_ghosts);

   int dst_id = variable_db->registerVariableAndContext(dst_variable,
         d_destination,
         dst_ghosts);

   TBOX_ASSERT(src_id != -1);
   TBOX_ASSERT(dst_id != -1);

   d_patch_data_components.setFlag(src_id);
   d_patch_data_components.setFlag(dst_id);

   if (d_do_refine) {
      std::shared_ptr<hier::RefineOperator> refine_operator(
         xfer_geom->lookupRefineOperator(src_variable, operator_name));

      hier::IntVector scratch_ghosts(hier::IntVector::max(src_ghosts,
                                        dst_ghosts));
      scratch_ghosts.max(hier::IntVector(scratch_ghosts.getDim(), 1));
      if (refine_operator) {
         scratch_ghosts.max(refine_operator->getStencilWidth(dim));
      }
      int scratch_id =
         variable_db->registerVariableAndContext(src_variable,
            d_refine_scratch,
            scratch_ghosts);
      TBOX_ASSERT(scratch_id != -1);

      d_patch_data_components.setFlag(scratch_id);

      d_refine_algorithm.registerRefine(dst_id,
         src_id,
         scratch_id,
         refine_operator);

      if (src_ghosts >= scratch_ghosts) {
         d_fill_source_algorithm.registerRefine(src_id,
            src_id,
            src_id,
            refine_operator);
      }
   } else if (d_do_coarsen) {
      std::shared_ptr<hier::CoarsenOperator> coarsen_operator(
         xfer_geom->lookupCoarsenOperator(src_variable, operator_name));
      d_coarsen_algorithm.registerCoarsen(dst_id,
         src_id,
         coarsen_operator);
   }

   registerVariableForReset(src_variable, dst_variable,
      src_ghosts, dst_ghosts, xfer_geom,
      operator_name);
}

void CommTester::registerVariableForReset(
   const std::shared_ptr<hier::Variable> src_variable,
   const std::shared_ptr<hier::Variable> dst_variable,
   const hier::IntVector& src_ghosts,
   const hier::IntVector& dst_ghosts,
   const std::shared_ptr<hier::BaseGridGeometry> xfer_geom,
   const std::string& operator_name)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(src_ghosts, dst_ghosts);

   TBOX_ASSERT(src_variable);
   TBOX_ASSERT(dst_variable);
   TBOX_ASSERT(xfer_geom);
   TBOX_ASSERT(!operator_name.empty());

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   int src_id = variable_db->registerVariableAndContext(src_variable,
         d_reset_source,
         src_ghosts);

   int dst_id = variable_db->registerVariableAndContext(dst_variable,
         d_reset_destination,
         dst_ghosts);

   d_patch_data_components.setFlag(src_id);
   d_patch_data_components.setFlag(dst_id);

   if (d_do_refine) {
      std::shared_ptr<hier::RefineOperator> refine_operator(
         xfer_geom->lookupRefineOperator(src_variable, operator_name));

      hier::IntVector scratch_ghosts(hier::IntVector::max(src_ghosts,
                                        dst_ghosts));
      scratch_ghosts.max(hier::IntVector(scratch_ghosts.getDim(), 1));
      if (refine_operator) {
         scratch_ghosts.max(refine_operator->getStencilWidth(scratch_ghosts.getDim()));
      }
      int scratch_id =
         variable_db->registerVariableAndContext(src_variable,
            d_reset_refine_scratch,
            scratch_ghosts);

      d_patch_data_components.setFlag(scratch_id);

      d_reset_refine_algorithm.registerRefine(dst_id,
         src_id,
         scratch_id,
         refine_operator);

   } else if (d_do_coarsen) {
      std::shared_ptr<hier::CoarsenOperator> coarsen_operator(
         xfer_geom->lookupCoarsenOperator(src_variable, operator_name));
      d_reset_coarsen_algorithm.registerCoarsen(dst_id,
         src_id,
         coarsen_operator);
   }

}

/*
 *************************************************************************
 *
 * Create refine and coarsen communication schedules for hierarchy.
 *
 *************************************************************************
 */

void CommTester::createRefineSchedule(
   const int level_number)
{
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= d_patch_hierarchy->getFinestLevelNumber()));

   std::shared_ptr<hier::PatchLevel> level(
      d_patch_hierarchy->getPatchLevel(level_number));

   if (d_do_refine) {

      d_fill_source_schedule.resize(d_patch_hierarchy->getNumberOfLevels());
      d_fill_source_schedule[level_number].reset();
      d_refine_schedule.resize(d_patch_hierarchy->getNumberOfLevels());
      d_refine_schedule[level_number].reset();

      const hier::Connector& peer_cnect =
         d_patch_hierarchy->getPatchLevel(level_number)->findConnector(
            *d_patch_hierarchy->getPatchLevel(level_number),
            d_patch_hierarchy->getRequiredConnectorWidth(level_number, level_number, true),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            false);
      const hier::Connector* cnect_to_coarser = level_number > 0 ?
         &d_patch_hierarchy->getPatchLevel(level_number)->findConnectorWithTranspose(
            *d_patch_hierarchy->getPatchLevel(level_number - 1),
            d_patch_hierarchy->getRequiredConnectorWidth(level_number, level_number - 1, true),
            d_patch_hierarchy->getRequiredConnectorWidth(level_number - 1, level_number),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            false) : 0;

      if (0) {
         // These are expensive checks.
         peer_cnect.assertOverlapCorrectness();
         if (cnect_to_coarser) {
            cnect_to_coarser->assertOverlapCorrectness();
            if (cnect_to_coarser->hasTranspose()) {
               cnect_to_coarser->getTranspose().assertOverlapCorrectness();
            }
         }
         d_patch_hierarchy->recursivePrint(tbox::plog, "", 3);
      }

      d_fill_source_schedule[level_number] =
         d_fill_source_algorithm.createSchedule(level, this);

      if ((level_number == 0) ||
          (d_refine_option == "INTERIOR_FROM_SAME_LEVEL")) {
         d_refine_schedule[level_number] =
            d_refine_algorithm.createSchedule(level,
               level_number - 1,
               d_patch_hierarchy,
               this);
      } else if (d_refine_option == "INTERIOR_FROM_COARSER_LEVEL") {
         d_refine_schedule[level_number] =
            d_refine_algorithm.createSchedule(level,
               std::shared_ptr<hier::PatchLevel>(),
               level_number - 1,
               d_patch_hierarchy,
               this);
      }

   }

}

void CommTester::resetRefineSchedule(
   const int level_number)
{
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= d_patch_hierarchy->getFinestLevelNumber()));

   if (d_do_refine) {

      d_reset_refine_algorithm.resetSchedule(d_refine_schedule[level_number]);

   }

   d_is_reset = true;
}

void CommTester::createCoarsenSchedule(
   const int level_number)
{
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= d_patch_hierarchy->getFinestLevelNumber()));

   if (d_do_coarsen && (level_number > 0)) {

      d_coarsen_schedule.resize(d_patch_hierarchy->getNumberOfLevels());
      d_coarsen_schedule[level_number].reset();

      std::shared_ptr<hier::PatchLevel> level(
         d_patch_hierarchy->getPatchLevel(level_number));
      std::shared_ptr<hier::PatchLevel> coarser_level(
         d_patch_hierarchy->getPatchLevel(level_number - 1));

      d_coarsen_schedule[level_number] =
         d_coarsen_algorithm.createSchedule(coarser_level,
            level,
            this);

   }

}

void CommTester::resetCoarsenSchedule(
   const int level_number)
{
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= d_patch_hierarchy->getFinestLevelNumber()));

   if (d_do_coarsen && (level_number > 0)) {

      d_reset_coarsen_algorithm.resetSchedule(d_coarsen_schedule[level_number]);

   }

   d_is_reset = true;
}

/*
 *************************************************************************
 *
 * Perform data refine and coarsen operations.
 *
 *************************************************************************
 */

void CommTester::performRefineOperations(
   const int level_number)
{
   if (d_do_refine) {
      if (d_fill_source_schedule[level_number] &&
          level_number < static_cast<int>(d_fill_source_schedule.size()) - 1) {
         d_data_test_strategy->setDataContext(d_source);
         d_fill_source_schedule[level_number]->fillData(d_fake_time);
         // synchronize is covered by RefineSchedule::recursiveFill at a finer grain
      }
      if (d_is_reset) {
         d_data_test_strategy->setDataContext(d_reset_refine_scratch);
      } else {
         d_data_test_strategy->setDataContext(d_refine_scratch);
      }
      if (d_refine_schedule[level_number]) {
         d_refine_schedule[level_number]->fillData(d_fake_time);
         // synchronize is covered by RefineSchedule::recursiveFill at a finer grain
      }
      d_data_test_strategy->clearDataContext();
   }
}

void CommTester::performCoarsenOperations(
   const int level_number)
{
   if (d_do_coarsen) {
      if (d_is_reset) {
         d_data_test_strategy->setDataContext(d_reset_source);
      } else {
         d_data_test_strategy->setDataContext(d_source);
      }
      if (d_coarsen_schedule[level_number]) {
         d_coarsen_schedule[level_number]->coarsenData();
         // synchronize is provided at a finer grain in coarsenData after communicate
      }
      d_data_test_strategy->clearDataContext();
   }
}

bool CommTester::performCompositeBoundaryComm(
      const int level_number)
{
   bool test_failed = false;

   d_data_test_strategy->setDataContext(d_destination);
   std::list<int> dst_ids;
   d_data_test_strategy->setDataIds(dst_ids);
   if (!dst_ids.empty()) {
      SAMRAI::xfer::CompositeBoundaryAlgorithm cba(d_patch_hierarchy, 3);
      for (std::list<int>::const_iterator itr = dst_ids.begin();
           itr != dst_ids.end(); ++itr) {
         cba.addDataId(*itr);
      }
      std::shared_ptr<SAMRAI::xfer::CompositeBoundarySchedule> cbsched =
         cba.createSchedule(level_number);
      cbsched->fillData(d_fake_time);
      // synchronize is covered at a finer grain

      std::shared_ptr<hier::PatchLevel> level(
         d_patch_hierarchy->getPatchLevel(level_number));

      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;

         for (std::list<int>::const_iterator itr = dst_ids.begin();
              itr != dst_ids.end(); ++itr) {
            const std::vector<std::shared_ptr<hier::PatchData> >& bdry_data =
               cbsched->getBoundaryPatchData(*patch, *itr);

            bool result = d_data_test_strategy->verifyCompositeBoundaryData(
               *patch, d_patch_hierarchy, *itr, level_number, bdry_data);

            if (!result) {
               test_failed = true;
            }
         }
      }
   }

   d_data_test_strategy->clearDataContext();

   return !test_failed;
}


/*
 *************************************************************************
 *
 * Verify results of communication operations.
 *
 *************************************************************************
 */

bool CommTester::verifyCommunicationResults() const
{
   bool tests_pass = true;
   if (d_is_reset) {
      d_data_test_strategy->setDataContext(d_reset_destination);
   } else {
      d_data_test_strategy->setDataContext(d_destination);
   }
   for (int ln = 0;
        ln <= d_patch_hierarchy->getFinestLevelNumber(); ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_patch_hierarchy->getPatchLevel(ln));

      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;

         tests_pass = tests_pass &&
            d_data_test_strategy->verifyResults(*patch, d_patch_hierarchy, ln);
      }
   }
   d_data_test_strategy->clearDataContext();

   return tests_pass;
}

/*
 *************************************************************************
 *
 * Cell tagging and patch level data initialization routines declared
 * in the GradientDetectorStrategy interface.  They are used to
 * construct the hierarchy initially.
 *
 *************************************************************************
 */

void CommTester::initializeLevelData(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const double time,
   const bool can_be_refined,
   const bool initial_time,
   const std::shared_ptr<hier::PatchLevel>& old_level,
   const bool allocate_data)
{
   NULL_USE(can_be_refined);
   NULL_USE(initial_time);
   NULL_USE(old_level);
   NULL_USE(allocate_data);

   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
   TBOX_ASSERT(level_number >= 0);

   hier::PatchLevel& level =
      (hier::PatchLevel &) * hierarchy->getPatchLevel(level_number);

   level.allocatePatchData(d_patch_data_components, time);

   for (hier::PatchLevel::iterator p(level.begin());
        p != level.end(); ++p) {
      hier::Patch& patch = **p;

      d_data_test_strategy->setDataContext(d_source);
      d_data_test_strategy->initializeDataOnPatch(patch,
         hierarchy,
         level.getLevelNumber(),
         's');
      d_data_test_strategy->clearDataContext();

      d_data_test_strategy->setDataContext(d_reset_source);
      d_data_test_strategy->initializeDataOnPatch(patch,
         hierarchy,
         level.getLevelNumber(),
         's');
      d_data_test_strategy->clearDataContext();

      if (d_do_coarsen) {

         /*
          * TODO: Why are we initializing destination data?  Shouldn't
          * it be uninitialized so we can check whether coarsening
          * actually set the data correctly?  Maybe we initialize it
          * so the un-refined coarse cells (which are unchanged by
          * coarsening) won't trigger a false positive during
          * verification.  If that is the problem, it should be fixed
          * by limiting the verification to refined coarse cells.
          * --BTNG
          */

         d_data_test_strategy->setDataContext(d_destination);
         d_data_test_strategy->initializeDataOnPatch(patch,
            hierarchy,
            level.getLevelNumber(),
            'd');
         d_data_test_strategy->clearDataContext();

         d_data_test_strategy->setDataContext(d_reset_destination);
         d_data_test_strategy->initializeDataOnPatch(patch,
            hierarchy,
            level.getLevelNumber(),
            'd');
         d_data_test_strategy->clearDataContext();

      }

   }

}

void CommTester::resetHierarchyConfiguration(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level)
{
   NULL_USE(hierarchy);
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);
}

/*
 *************************************************************************
 *
 * Physical boundary condition and user-defined coarsen and refine
 * operations declared in RefinePatchStrategy and CoarsenPatchStrategy.
 * They are passed off to patch data test object.
 *
 *************************************************************************
 */

void CommTester::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw)
{
   NULL_USE(time);
   d_data_test_strategy->setPhysicalBoundaryConditions(patch,
      d_fake_time,
      gcw);
}

hier::IntVector CommTester::getRefineOpStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getOne(dim);
}

void CommTester::preprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{
   d_data_test_strategy->preprocessRefine(fine, coarse, fine_box, ratio);
}

void CommTester::postprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{
   d_data_test_strategy->postprocessRefine(fine, coarse, fine_box, ratio);
}

hier::IntVector CommTester::getCoarsenOpStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getZero(dim);
}

void CommTester::preprocessCoarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio)
{
   d_data_test_strategy->preprocessCoarsen(coarse, fine, coarse_box, ratio);
}

void CommTester::postprocessCoarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio)
{
   d_data_test_strategy->postprocessCoarsen(coarse, fine, coarse_box, ratio);
}

/*
 *************************************************************************
 *
 * Create and configure gridding objects used to build the hierarchy.
 * Then, create hierarchy and initialize data.  Note this routine
 * must be called after variables are registered with this tester object.
 *
 *************************************************************************
 */

void CommTester::setupHierarchy(
   std::shared_ptr<tbox::Database> main_input_db,
   std::shared_ptr<mesh::StandardTagAndInitialize> cell_tagger)
{
   TBOX_ASSERT(main_input_db);

   d_patch_hierarchy.reset(
      new hier::PatchHierarchy("PatchHierarchy",
         d_data_test_strategy->getGridGeometry(),
         main_input_db->getDatabase("PatchHierarchy")));

   std::shared_ptr<mesh::BergerRigoutsos> box_generator(
      new mesh::BergerRigoutsos(d_dim,
         main_input_db->getDatabase("BergerRigoutsos")));

   std::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
      new mesh::TreeLoadBalancer(
         d_dim,
         "TreeLoadBalancer",
         main_input_db->getDatabase("TreeLoadBalancer")));
   load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

   std::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm(
      new mesh::GriddingAlgorithm(
         d_patch_hierarchy,
         "GriddingAlgorithm",
         main_input_db->getDatabase("GriddingAlgorithm"),
         cell_tagger,
         box_generator,
         load_balancer));

   int fake_tag_buffer = 0;

   gridding_algorithm->makeCoarsestLevel(d_fake_time);

   bool initial_cycle = true;
   for (int ln = 0; d_patch_hierarchy->levelCanBeRefined(ln); ++ln) {
      gridding_algorithm->makeFinerLevel(
         fake_tag_buffer,
         initial_cycle,
         d_fake_cycle,
         d_fake_time);
   }

   /*
    * Clear the PersistentOverlapConnectors so we can test whether
    * the communication schedules can work properly in the mode where
    * there are no precomputed Connectors.
    */
   for (int ln = 0; ln < d_patch_hierarchy->getNumberOfLevels(); ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_patch_hierarchy->getPatchLevel(ln));
      level->getBoxLevel()->clearPersistentOverlapConnectors();
   }

   if (0) {
      tbox::plog << "h:  generated hierarchy:\n";
      d_patch_hierarchy->recursivePrint(tbox::plog, "h:  ", 3);
      tbox::plog << "h:  box_level hierarchy:\n";
      d_patch_hierarchy->recursivePrint(tbox::plog,
         "",
         3);
   }

}

void
CommTester::putCoordinatesToDatabase(
   std::shared_ptr<tbox::Database>& coords_db,
   const hier::Patch& patch,
   const hier::Box& box)
{
   NULL_USE(box);

   std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
/*
   if (pgeom) {
      pgeom->putBlueprintCoords(coords_db, patch.getBox());
   }
*/
   const tbox::Dimension& dim(patch.getDim());

   pdat::NodeData<double> coords(patch.getBox(), dim.getValue(),
                                 hier::IntVector::getZero(dim));

   const hier::Index& box_lo = patch.getBox().lower();
   const double* x_lo = pgeom->getXLower();
   const double* dx = pgeom->getDx();
   
   pdat::NodeIterator nend = pdat::NodeGeometry::end(patch.getBox());
   for (pdat::NodeIterator itr(pdat::NodeGeometry::begin(patch.getBox())); itr != nend; ++itr) {
      const pdat::NodeIndex& ni = *itr;
      for (int d = 0; d < dim.getValue(); ++d) {
         coords(ni, d) = x_lo[d] + (ni(d)-box_lo(d))*dx[d];
      }
   }

   coords_db->putString("type", "explicit");

   std::shared_ptr<tbox::Database> values_db = coords_db->putDatabase("values");

   int data_size = coords.getArrayData().getBox().size();

   values_db->putDoubleArray("x", coords.getPointer(0), data_size);
   if (dim.getValue() > 1) {
      values_db->putDoubleArray("y", coords.getPointer(1), data_size);
   }
   if (dim.getValue() > 2) {
      values_db->putDoubleArray("z", coords.getPointer(2), data_size);
   }
}


}
