/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   AMR hierarchy generation and regridding routines.
 *
 ************************************************************************/
#include "SAMRAI/mesh/GriddingAlgorithm.h"

#include "SAMRAI/tbox/IEEE.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/math/PatchCellDataBasicOps.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/pdat/CellIntegerConstantRefine.h"
#include "SAMRAI/pdat/CellConstantRefine.h"
#include "SAMRAI/xfer/PatchInteriorVariableFillPattern.h"
#include "SAMRAI/xfer/PatchLevelInteriorFillPattern.h"
#include "SAMRAI/tbox/AllocatorDatabase.h"
#include "SAMRAI/tbox/Collectives.h"
#include "SAMRAI/tbox/NVTXUtilities.h"
#include "SAMRAI/tbox/Utilities.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <iomanip>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

const int GriddingAlgorithm::ALGS_GRIDDING_ALGORITHM_VERSION = 3;

const int GriddingAlgorithm::EXISTING_FINE_PATCHES_TAG_VAL = 1;
const int GriddingAlgorithm::NEW_FINE_PATCHES_TAG_VAL = -4;
const int GriddingAlgorithm::BUFFER_TAG_VAL = -9;

/*
 *************************************************************************
 *
 * Constructor and destructor for GriddingAlgorithm.
 *
 *************************************************************************
 */
GriddingAlgorithm::GriddingAlgorithm(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const std::string& object_name,
   const std::shared_ptr<tbox::Database>& input_db,
   const std::shared_ptr<TagAndInitializeStrategy>& tag_init_strategy,
   const std::shared_ptr<BoxGeneratorStrategy>& generator,
   const std::shared_ptr<LoadBalanceStrategy>& balancer,
   const std::shared_ptr<LoadBalanceStrategy>& balancer0):
   GriddingAlgorithmStrategy(),
   d_hierarchy(hierarchy),
   d_connector_width_requestor(),
   d_buf_tag_ghosts(hierarchy->getDim(), 1),
   d_object_name(object_name),
   d_tag_init_strategy(tag_init_strategy),
   d_box_generator(generator),
   d_load_balancer(balancer),
   d_load_balancer0(balancer0 ? balancer0 : balancer),
   d_mb_tagger_strategy(0),
   d_true_tag(1),
   d_false_tag(0),
   d_from_fine_pretag(EXISTING_FINE_PATCHES_TAG_VAL),
   d_from_fine_tag(NEW_FINE_PATCHES_TAG_VAL),
   d_buffer_tag(BUFFER_TAG_VAL),
   d_base_ln(-1),
   d_tag_to_cluster_width(),
   d_check_nonrefined_tags('w'),
   d_check_overlapping_patches('w'),
   d_check_nonnesting_user_boxes('e'),
   d_check_boundary_proximity_violation('e'),
   d_sequentialize_patch_indices(true),
   d_log_metadata_statistics(false),
   d_oca(),
   d_mca(),
   d_blcu(),
   d_oca0(),
   d_mca0(),
   d_blcu0(),
   d_enforce_proper_nesting(true),
   d_extend_to_domain_boundary(true),
   d_load_balance(true),
   d_save_tag_data(false),
   d_barrier_and_time(false),
   d_check_overflow_nesting(false),
   d_check_proper_nesting(false),
   d_check_connectors(false),
   d_print_steps(false)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(tag_init_strategy);
   TBOX_ASSERT(generator);
   TBOX_ASSERT(balancer);

   tbox::RestartManager::getManager()->
   registerRestartItem(d_object_name, this);

   d_hierarchy->registerConnectorWidthRequestor(
      d_connector_width_requestor);

   const tbox::Dimension& dim = hierarchy->getDim();

   /*
    * Construct integer tag variables and add to variable database.
    * Note that variables and patch data indices are shared among
    * all GriddingAlgorithm instances.  The VariableDatabase
    * holds the variables, once contructed and registered via the
    * VariableDatabase::makeInternalSAMRAIWorkVariablePatchDataIndex()
    * function call.  Note that variables are registered and patch data
    * indices are made only for the first time through the constructor.
    */
   hier::VariableDatabase* var_db = hier::VariableDatabase::getDatabase();

   std::string tag_interior_variable_name("GriddingAlgorithm__tag-interior");
   std::string tag_saved_variable_name("GriddingAlgorithm__tag-saved");
   std::string tag_algorithm_variable_name("GriddingAlgorithm__tag-algorithm");
   std::string tag_buffer_variable_name("GriddingAlgorithm__tag-buffer");

   std::ostringstream dim_extension;
   dim_extension << "_" << dim.getValue();

   tag_interior_variable_name += dim_extension.str();
   tag_saved_variable_name += dim_extension.str();
   tag_algorithm_variable_name += dim_extension.str();
   tag_buffer_variable_name += dim_extension.str();

   d_user_tag = std::dynamic_pointer_cast<pdat::CellVariable<int>, hier::Variable>(
         var_db->getVariable(tag_interior_variable_name));
   if (!d_user_tag) {
      d_user_tag.reset(
         new pdat::CellVariable<int>(
           dim,
           tag_interior_variable_name,
           tbox::AllocatorDatabase::getDatabase()->getTagAllocator(),
           1));
   }

   d_saved_tag = std::dynamic_pointer_cast<pdat::CellVariable<int>, hier::Variable>(
         var_db->getVariable(tag_saved_variable_name));
   if (!d_saved_tag) {
      d_saved_tag.reset(
         new pdat::CellVariable<int>(dim,
           tag_saved_variable_name,
           tbox::AllocatorDatabase::getDatabase()->getTagAllocator(),
           1));
   }

   d_boolean_tag = std::dynamic_pointer_cast<pdat::CellVariable<int>, hier::Variable>(
         var_db->getVariable(tag_algorithm_variable_name));
   if (!d_boolean_tag) {
#if defined(HAVE_UMPIRE)
      umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
      tbox::ResourceAllocator host_alloc = rm.getAllocator(umpire::resource::Host);
#else
      tbox::ResourceAllocator host_alloc;
#endif

      d_boolean_tag.reset(
         new pdat::CellVariable<int>(dim,
           tag_algorithm_variable_name,
           host_alloc,
           1));
   }

   d_buf_tag = std::dynamic_pointer_cast<pdat::CellVariable<int>, hier::Variable>(
         var_db->getVariable(tag_buffer_variable_name));
   if (!d_buf_tag) {
      d_buf_tag.reset(new pdat::CellVariable<int>(dim,
            tag_buffer_variable_name,
            tbox::AllocatorDatabase::getDatabase()->getTagAllocator(),
            1));
   }

   d_user_tag_indx = var_db->registerInternalSAMRAIVariable(d_user_tag,
            hier::IntVector::getZero(dim));
   d_saved_tag_indx = var_db->registerInternalSAMRAIVariable(d_saved_tag,
            hier::IntVector::getZero(dim));
   d_boolean_tag_indx = var_db->registerInternalSAMRAIVariable(d_boolean_tag,
            hier::IntVector::getZero(dim));
   d_buf_tag_indx = var_db->registerInternalSAMRAIVariable(d_buf_tag,
            hier::IntVector::getOne(dim));
   d_buf_tag_ghosts = hier::IntVector::getOne(dim);

   if (d_hierarchy->getGridGeometry()->getNumberBlocks() > 1) {
      d_mb_tagger_strategy = new MultiblockGriddingTagger();
      d_mb_tagger_strategy->setScratchTagPatchDataIndex(d_buf_tag_indx);
   }

   /*
    * Initialize communication algorithm for exchanging buffered tag data.
    */
   d_bdry_fill_tags.reset(new xfer::RefineAlgorithm());

   d_bdry_fill_tags->registerRefine(d_buf_tag_indx,
      d_buf_tag_indx,
      d_buf_tag_indx,
      std::shared_ptr<hier::RefineOperator>());

   d_fill_saved_tags.reset(new xfer::RefineAlgorithm());

   d_fill_saved_tags->registerRefine(d_saved_tag_indx,
      d_user_tag_indx,
      d_saved_tag_indx,
      std::shared_ptr<hier::RefineOperator>(new pdat::CellConstantRefine<int>()),
      std::shared_ptr<xfer::VariableFillPattern>(new xfer::PatchInteriorVariableFillPattern(dim)));

   d_proper_nesting_complement.resize(d_hierarchy->getMaxNumberOfLevels());
   d_to_nesting_complement.resize(d_hierarchy->getMaxNumberOfLevels());

   /*
    * Initialize object with data read from input and restart databases.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }

   getFromInput(input_db, is_from_restart);

   if (d_hierarchy->getMaxNumberOfLevels() > 1) {
      std::vector<hier::IntVector> ratio_to_coarser(
         d_hierarchy->getMaxNumberOfLevels(),
         hier::IntVector(dim,1));
      for (int ln = 0; ln < d_hierarchy->getMaxNumberOfLevels(); ++ln) {
         ratio_to_coarser[ln] = d_hierarchy->getRatioToCoarserLevel(ln);
      }
      d_tag_init_strategy->checkCoarsenRatios(ratio_to_coarser);
      // Note: checkCoarsenRatios() must precede getErrorCoarsenRatio().
   }
   if (d_tag_init_strategy->getErrorCoarsenRatio() > 1) {
      std::shared_ptr<StandardTagAndInitialize> std_tag_init(
         std::dynamic_pointer_cast<StandardTagAndInitialize, TagAndInitializeStrategy>(
            d_tag_init_strategy));
      if (std_tag_init) {
         d_hierarchy->registerConnectorWidthRequestor(
            std_tag_init->getConnectorWidthRequestor());
      }
   }

   d_bdry_sched_tags.resize(d_hierarchy->getMaxNumberOfLevels());

   d_oca.setSAMRAI_MPI(d_hierarchy->getDomainBoxLevel().getMPI(), true);
   d_mca.setSAMRAI_MPI(d_hierarchy->getDomainBoxLevel().getMPI(), true);
   d_oca0.setSAMRAI_MPI(d_hierarchy->getDomainBoxLevel().getMPI(), true);
   d_mca0.setSAMRAI_MPI(d_hierarchy->getDomainBoxLevel().getMPI(), true);

   warnIfDomainTooSmallInPeriodicDir();

   allocateTimers();

#ifdef GA_RECORD_STATS

   d_boxes_stat.resize(d_hierarchy->getMaxNumberOfLevels());
   d_cells_stat.resize(d_hierarchy->getMaxNumberOfLevels());
   d_timestamp_stat.resize(d_hierarchy->getMaxNumberOfLevels());

   for (int ln = 0; ln < d_hierarchy->getMaxNumberOfLevels(); ++ln) {
      std::string ln_text = tbox::Utilities::intToString(ln, 2);
      const std::string num_boxes_str = std::string("GA_BoxesL") + ln_text;
      const std::string num_cells_str = std::string("GA_CellsL") + ln_text;
      const std::string timestamp_str = std::string("GA_TimeL") + ln_text;
      d_boxes_stat[ln] = tbox::Statistician::getStatistician()->
         getStatistic(num_boxes_str, "PROC_STAT");
      d_cells_stat[ln] = tbox::Statistician::getStatistician()->
         getStatistic(num_cells_str, "PROC_STAT");
      d_timestamp_stat[ln] = tbox::Statistician::getStatistician()->
         getStatistic(timestamp_str, "PROC_STAT");
   }

   d_oca.setTimerPrefix("mesh::GriddingAlgorithm");
   d_mca.setTimerPrefix("mesh::GriddingAlgorithm");
   d_blcu.setTimerPrefix("mesh::GriddingAlgorithm");
   d_oca0.setTimerPrefix("mesh::GriddingAlgorithm0");
   d_mca0.setTimerPrefix("mesh::GriddingAlgorithm0");
   d_blcu0.setTimerPrefix("mesh::GriddingAlgorithm0");

#endif

}

/*
 *************************************************************************
 *
 * Destructor tells the tbox::RestartManager to remove this object from
 * the list of restart items.
 *
 *************************************************************************
 */
GriddingAlgorithm::~GriddingAlgorithm()
{
   tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
   delete d_mb_tagger_strategy;
   d_mb_tagger_strategy = 0;
}

std::shared_ptr<TagAndInitializeStrategy>
GriddingAlgorithm::getTagAndInitializeStrategy() const
{
   return d_tag_init_strategy;
}

std::shared_ptr<LoadBalanceStrategy>
GriddingAlgorithm::getLoadBalanceStrategy() const
{
   return d_load_balancer;
}

std::shared_ptr<LoadBalanceStrategy>
GriddingAlgorithm::getLoadBalanceStrategyZero() const
{
   return d_load_balancer0;
}

std::shared_ptr<hier::PatchHierarchy>
GriddingAlgorithm::getPatchHierarchy() const
{
   return d_hierarchy;
}

/*
 *************************************************************************
 *
 * Construct coarsest level in the patch hierarchy (i.e., level 0).
 * This routine operates in two modes:
 *
 * (1) If level 0 does not exist in the hierarchy, then a new level
 *  will be created based on the domain description provided by
 *  the grid geometry associated with the hierarchy.
 * (2) If level 0 exists in the hierarchy, then a new level is made
 *  by re-applying the load balancing routines to the existing level.
 *  The pre-existing level will be removed from the hierarchy.
 *
 * In either case, the new level is placed in the hierarchy as level 0.
 * If level 0 can be refined (i.e., it will be subject to tagging cells
 * for refinement), the domain boxes are checked against the constraints
 * of regridding process.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::makeCoarsestLevel(
   const double level_time)
{
   if (d_print_steps) {
      tbox::plog << "GriddingAlgorithm::makeCoarsestLevel: entered with getNumberOfLevels = "
                 << d_hierarchy->getNumberOfLevels() << "\n";
   }

   const tbox::Dimension& dim = d_hierarchy->getDim();

   if (d_barrier_and_time) {
      t_make_coarsest->barrierAndStart();
   }

   if (d_tag_to_cluster_width.empty()) {
      computeTagToClusterWidths();
   }

   const hier::IntVector& zero_vec(hier::IntVector::getZero(dim));
   const size_t num_blocks = d_hierarchy->getGridGeometry()->getNumberBlocks();

   TBOX_ASSERT(d_hierarchy->getMaxNumberOfLevels() > 0);

   t_make_domain->start();

   const int ln = 0;

   bool level_zero_exists = false;
   if ((d_hierarchy->getNumberOfLevels() > 0)) {
      if (d_hierarchy->getPatchLevel(ln)) {
         level_zero_exists = true;
      }
   }

   hier::IntVector smallest_patch(dim);
   hier::IntVector largest_patch(dim);
   hier::IntVector extend_ghosts(dim);
   {
      hier::IntVector smallest_box_to_refine(dim, 0, num_blocks);
      // "false" argument: for_building_finer level = false
      getGriddingParameters(
         smallest_patch,
         smallest_box_to_refine,
         largest_patch,
         extend_ghosts,
         ln,
         false);
   }

   /*
    * If making coarsest level for the first time, check domain boxes
    * for violations of user constraints.
    */
   if (!level_zero_exists) {
      for (unsigned int b = 0; b < num_blocks; ++b) {
         hier::BoxContainer domain_boxes(
            d_hierarchy->getGridGeometry()->getPhysicalDomain(),
            hier::BlockId(b));
         checkDomainBoxes(domain_boxes);
      }
   }

   const hier::BoxLevel& domain_box_level(d_hierarchy->getDomainBoxLevel());

   TBOX_ASSERT(!domain_box_level.getMPI().hasReceivableMessage());    // Check errant messages.

   t_make_domain->stop();

   /*
    * domain_to_domain is a temporary connector for bridging the
    * balanced BoxLevel to itself.  Since domain is small, owned
    * by just proc 0 and already globalized, computing
    * domain_to_domain is fast and does not require communication.
    */

   std::shared_ptr<hier::Connector> domain_to_domain;

   if (d_print_steps) {
      tbox::plog << "GriddingAlgorithm::makeCoarsestLevel: finding domain<==>domain.\n";
   }

   d_oca0.findOverlaps(domain_to_domain,
      domain_box_level,
      domain_box_level,
      hier::IntVector::max(
         d_hierarchy->getRequiredConnectorWidth(0, 0, true),
         hier::IntVector::getOne(dim)));

   if (d_barrier_and_time) {
      t_load_balance0->barrierAndStart();
   }

   hier::IntVector patch_cut_factor(dim, d_tag_init_strategy->getErrorCoarsenRatio());

   /*
    * TODO: The code for generating the coarsest level's boxes is not
    * scalable because we use the domain BoxLevel in bridges and
    * Connector modifications, forcing proc 0 (which owns all domain
    * nodes) to do all of the searches.  This problem will show up in
    * performance data if we re-partition the coarsest level
    * regularly.
    */

   std::shared_ptr<hier::BoxLevel> new_box_level(
      std::make_shared<hier::BoxLevel>(domain_box_level));
   hier::Connector domain_to_new(*domain_to_domain);
   domain_to_new.setHead(*new_box_level, true);
   hier::Connector new_to_domain(*domain_to_domain);
   new_to_domain.setBase(*new_box_level, true);
   new_to_domain.setTranspose(&domain_to_new, false);

   if (d_print_steps) {
      tbox::plog << "GriddingAlgorithm::makeCoarsestLevel: partitioning domain.\n";
   }

   TBOX_ASSERT(!new_box_level->getMPI().hasReceivableMessage());    // Check errant messages.

   std::shared_ptr<hier::Connector> new_to_zero;
   std::shared_ptr<hier::Connector> zero_to_new;
   if (level_zero_exists) {
      /*
       * If level zero exists, that means we are rebalancing level zero.
       *
       * Most calls to loadBalanceBoxLevel pass as the second argument
       * a pointer to a "balance_to_reference" Connector, which is a
       * Connector between the BoxLevel to be load balanced (here,
       * new_box_levle) and a reference level -- that is a BoxLevel that
       * will not be changed during the loadBalanceBoxLevel call.  When load
       * balancing a finer level in the AMR hierarchy, the reference level
       * is typically the current next coarser level, level N-1 when balancing
       * level N.
       *
       * As there is no level N-1 when balancing level zero, here we use
       * the current level zero as the reference level and construct the
       * needed Connector.
       *
       * If level zero doesn't exist, this section is skipped and a null
       * pointer is passed to loadBalanceBoxLevel.  A null pointer is a
       * valid argument when creating the first instance of level zero.
       */

      hier::IntVector connector_width(
         hier::IntVector::max(
            d_hierarchy->getRequiredConnectorWidth(0, 0),
            hier::IntVector::getOne(dim)));

      new_to_zero = std::make_shared<hier::Connector>(
         *new_box_level,
         *(d_hierarchy->getPatchLevel(0)->getBoxLevel()),
         connector_width);
      d_oca0.findOverlaps_assumedPartition(*new_to_zero);

      zero_to_new = std::make_shared<hier::Connector>(dim);
      zero_to_new->computeTransposeOf(*new_to_zero,
                                      new_box_level->getMPI());

      new_to_zero->setTranspose(zero_to_new.get(), false);
   }

   d_load_balancer0->loadBalanceBoxLevel(
      *new_box_level,
      new_to_zero.get(),
      d_hierarchy,
      ln,
      smallest_patch,
      largest_patch,
      domain_box_level,
      extend_ghosts,
      patch_cut_factor);

   if (d_barrier_and_time) {
      t_load_balance0->stop();
   }

   if (d_sequentialize_patch_indices) {
      if (d_print_steps) {
         tbox::plog << "GriddingAlgorithm::makeCoarsestLevel: renumbering boxes.\n";
      }
      renumberBoxes(*new_box_level, 0, false, true);
   }

   if (d_print_steps) {
      tbox::plog << "GriddingAlgorithm::makeCoarsestLevel: adding periodic images.\n";
   }
   d_blcu0.addPeriodicImages(
      *new_box_level,
      d_hierarchy->getGridGeometry()->getDomainSearchTree(),
      hier::IntVector::max(
         d_hierarchy->getRequiredConnectorWidth(0, 0, true),
         hier::IntVector::getOne(dim)));

   TBOX_ASSERT(!new_box_level->getMPI().hasReceivableMessage());    // Check errant messages.

   std::shared_ptr<hier::Connector> new_to_new;
   if (domain_box_level.getLocalNumberOfBoxes(0) ==
       domain_box_level.getGlobalNumberOfBoxes()) {
      /*
       * If proc 0 owns all new boxes, it is faster find new<==>new by
       * globalizing the new boxes.
       *
       * The standard approach of bridging basically does the same,
       * but forces proc 0 to compute all the overlaps and send that
       * info to each processor one at a time.
       */

      if (d_barrier_and_time) {
         t_find_new_to_new->barrierAndStart();
      }

      if (d_print_steps) {
         tbox::plog << "GriddingAlgorithm::makeCoarsestLevel: finding new<==>new." << std::endl;
      }
      new_to_new.reset(new hier::Connector(*new_box_level, *new_box_level,
            d_hierarchy->getRequiredConnectorWidth(0, 0, true)));
      d_oca0.findOverlaps_assumedPartition(*new_to_new);

      if (d_barrier_and_time) {
         t_find_new_to_new->stop();
      }

   } else {

      if (d_barrier_and_time) {
         t_bridge_new_to_new->barrierAndStart();
      }

      if (d_print_steps) {
         tbox::plog << "GriddingAlgorithm::makeCoarsestLevel: bridging for domain<==>domain.\n";
      }
      d_oca0.bridgeWithNesting(
         new_to_new,
         new_to_domain,
         domain_to_new,
         zero_vec,
         zero_vec,
         d_hierarchy->getRequiredConnectorWidth(0, 0, true),
         false);

      if (d_barrier_and_time) {
         t_bridge_new_to_new->stop();
      }

      TBOX_ASSERT(new_to_new->getConnectorWidth() ==
         d_hierarchy->getRequiredConnectorWidth(0, 0));
      TBOX_ASSERT(&new_to_new->getBase() == new_box_level.get());
      TBOX_ASSERT(&new_to_new->getHead() == new_box_level.get());

   }

   if (d_check_overlapping_patches != 'i') {
      checkOverlappingPatches(*new_to_new);
   }

   if (d_print_steps) {
      tbox::plog << "GriddingAlgorithm::makeCoarsestLevel: making level 0.\n";
   }
   t_make_new->start();
   if (!level_zero_exists) {

      d_hierarchy->makeNewPatchLevel(ln, new_box_level);
      /*
       * Add computed Connectors to new level's collection of
       * persistent overlap Connectors.
       */
      if (!new_box_level->hasConnector(new_to_new->getHead(),
             new_to_new->getConnectorWidth())) {
         new_box_level->cacheConnector(new_to_new);
      }

      d_hierarchy->getGridGeometry()->adjustMultiblockPatchLevelBoundaries(
         *d_hierarchy->getPatchLevel(ln));

      if (d_print_steps) {
         tbox::plog << "GriddingAlgorithm::makeCoarsestLevel: initializing data on level 0 (a).\n";
      }
      // "true" argument: const bool initial_time = true;
      d_tag_init_strategy->initializeLevelData(d_hierarchy,
         ln,
         level_time,
         d_hierarchy->levelCanBeRefined(ln),
         true);

   } else {

      /*
       * Save old data before they are overwritten by the new BoxLevel.
       */
      std::shared_ptr<hier::PatchLevel> old_level(
         d_hierarchy->getPatchLevel(ln));

      /*
       * Compute old<==>new.  Doing it this way is not scalable, but
       * we only do this for the coarsest level.  The old approach of
       * bridging across the domain BoxLevel is probably not
       * scalable anyway, because the domain is usually owned by just
       * one processor.
       */
      old_level->getBoxLevel()->createConnectorWithTranspose(*new_box_level,
         d_hierarchy->getRequiredConnectorWidth(0, 0, true),
         d_hierarchy->getRequiredConnectorWidth(0, 0, true));

      d_tag_init_strategy->processHierarchyBeforeAddingNewLevel(d_hierarchy,
         ln,
         new_box_level);

      d_hierarchy->removePatchLevel(ln);

      d_hierarchy->makeNewPatchLevel(ln, new_box_level);

      d_hierarchy->getGridGeometry()->adjustMultiblockPatchLevelBoundaries(
         *d_hierarchy->getPatchLevel(ln));

      if (d_print_steps) {
         tbox::plog << "GriddingAlgorithm::makeCoarsestLevel: initializing data on level 0 (b).\n";
      }
      // "false" argument: const bool initial_time = false;
      d_tag_init_strategy->initializeLevelData(d_hierarchy,
         ln,
         level_time,
         d_hierarchy->levelCanBeRefined(ln),
         false,
         old_level);

      if (d_hierarchy->getNumberOfLevels() > 1) {
         const hier::Connector& new_to_old =
            d_hierarchy->getPatchLevel(ln)->findConnector(
               *old_level,
               d_hierarchy->getRequiredConnectorWidth(0, 0),
               hier::CONNECTOR_IMPLICIT_CREATION_RULE,
               true);

         const hier::Connector& old_to_fine =
            old_level->findConnector(
               *(d_hierarchy->getPatchLevel(ln+1)),
               d_hierarchy->getRequiredConnectorWidth(0, 1),
               hier::CONNECTOR_IMPLICIT_CREATION_RULE,
               true);

         std::shared_ptr<hier::Connector> new_to_fine;
         d_oca.bridgeWithNesting(
            new_to_fine,
            new_to_old,
            old_to_fine,
            zero_vec,
            zero_vec,
            d_hierarchy->getRequiredConnectorWidth(0, 1),
            true);

         std::shared_ptr<hier::Connector> new_to_new;
         d_oca.bridgeWithNesting(
            new_to_new,
            new_to_old,
            new_to_old.getTranspose(),
            zero_vec,
            zero_vec,
            d_hierarchy->getRequiredConnectorWidth(0, 0),
            false);

         d_hierarchy->getPatchLevel(0)->cacheConnector(new_to_fine);
         d_hierarchy->getPatchLevel(0)->cacheConnector(new_to_new);

      }

      old_level.reset();
   }

   std::shared_ptr<hier::PatchLevel> new_level(
      d_hierarchy->getPatchLevel(ln));

   new_level->allocatePatchData(d_saved_tag_indx);
   fillTags(d_false_tag, new_level, d_saved_tag_indx);

   if (d_print_steps) {
      tbox::plog << "GriddingAlgorithm::makeCoarsestLevel: made level 0.\n";
   }
   t_make_new->stop();

   if (d_barrier_and_time) {
      t_reset_hier->barrierAndStart();
   }
   d_tag_init_strategy->resetHierarchyConfiguration(d_hierarchy, ln, ln);
   if (d_barrier_and_time) {
      t_reset_hier->stop();
   }

   if (d_log_metadata_statistics) {
      d_hierarchy->logMetadataStatistics("makeCoarsestLevel", 0, 0, level_time, false, false);
   }

#ifdef GA_RECORD_STATS
   recordStatistics(level_time);
#endif

   if (d_barrier_and_time) {
      t_make_coarsest->stop();
   }

   if (d_print_steps) {
      tbox::plog << "GriddingAlgorithm::makeCoarsestLevel: returning.\n";
   }
}

/*
 *************************************************************************
 *
 * Perform error estimation process on the finest hierarchy level to
 * determine if and where a new finest level is needed.  If it is
 * determined  that a new finest level should exist, it is created and
 * its patch data is allocated and initialized.  The algorithm is
 * summarized as follows:
 *
 * (1) Compute boxes whose union constitutes the region within the level
 *  in which the next finer level may reside (i.e., proper nesting).
 *
 * (2) Set tags on the level to indicate which cells should be refined.
 *
 * (3) Buffer the tags.  This prevents disturbances from moving off
 *  refined regions before the next remesh occurs.
 *
 * (4) Determine boxes for new patches that will cover the tagged cells.
 *  This step includes load balancing of the patches.
 *
 * (5) If there exist some regions to refine, construct the next finer
 *  level and insert it in the hierarchy.  Then, initialize the data
 *  on that level.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::makeFinerLevel(
   const int tag_buffer,
   const bool initial_cycle,
   const int cycle,
   const double level_time,
   const double regrid_start_time)
{
   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::makeFinerLevel: entered with finest ln = "
      << d_hierarchy->getFinestLevelNumber() << "\n";
   }

   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT(d_hierarchy->getPatchLevel(
         d_hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(tag_buffer >= 0);

   if (d_barrier_and_time) {
      t_make_finer->barrierAndStart();
   }

   if (d_tag_to_cluster_width.empty()) {
      computeTagToClusterWidths();
   }

   const tbox::Dimension& dim = d_hierarchy->getDim();

   const hier::IntVector& zero_vector(hier::IntVector::getZero(dim));

   const int tag_ln = d_hierarchy->getFinestLevelNumber();
   const int new_ln = tag_ln + 1;

   if (d_hierarchy->levelCanBeRefined(tag_ln)) {

      /*
       * d_base_ln is used by private methods during regrid.
       */
      d_base_ln = tag_ln;

      /*
       * Compute nesting data at d_base_ln for use in constructing
       * level d_base_ln+1;
       */
      computeProperNestingData(d_base_ln, d_oca);

      const std::shared_ptr<hier::PatchLevel> tag_level(
         d_hierarchy->getPatchLevel(tag_ln));

      std::shared_ptr<hier::BoxLevel> new_box_level;
      std::shared_ptr<hier::Connector> tag_to_new;
      hier::Connector* new_to_tag = 0;

      /*
       * The boolean "do_tagging" specifies whether or not tagging will
       * be performed. This will be true except in two circumstances:
       *   1) only user supplied refine boxes are used
       *   2) the boxes are read from a previously dumped file.
       *
       * If either of these circumstances is true, tagging operations
       * are NOT necessary so do_tagging will be set to false.
       */
      bool do_tagging = true;
      if (d_tag_init_strategy->refineUserBoxInputOnly(cycle, level_time)) {
         do_tagging = false;
      }

      /*
       * Tag cells, determine refine boxes from tagged regions, and
       * load balance in preparation for constructing new refined level.
       */
      if (do_tagging) {

         /*
          * Initialize integer tag arrays on level to false.
          */
         tag_level->allocatePatchData(d_user_tag_indx, level_time);
         fillTags(d_false_tag, tag_level, d_user_tag_indx);

         /*
          * Perform pre-processing of error estimation data.
          */
         d_tag_init_strategy->preprocessErrorEstimation(d_hierarchy,
            tag_ln,
            cycle,
            level_time,
            regrid_start_time,
            initial_cycle);

         /*
          * Determine cells needing refinement on level and set tags to true.
          * Because we are constructing a new level, not regridding the level,
          * the coarsest_sync_level argument is always false.
          */
         bool coarsest_sync_level = false;
         d_tag_init_strategy->tagCellsForRefinement(d_hierarchy,
            tag_ln,
            cycle,
            level_time,
            d_user_tag_indx,
            initial_cycle,
            coarsest_sync_level,
            d_hierarchy->levelCanBeRefined(tag_ln),
            regrid_start_time);

         /*
          * Set algorithmic tags to false/true values that will be understood
          * in the buffering and clustering steps.
          */
         tag_level->allocatePatchData(d_boolean_tag_indx, level_time);
         setBooleanTagData(tag_level, false);

         /*
          * Check for user-tagged cells that violate proper nesting.
          * except if user specified that the violating tags be ignored.
          */
         if (d_check_nonrefined_tags != 'i') {
            checkNonrefinedTags(*tag_level, tag_ln, d_oca);
         }

         /*
          * Buffer true tagged cells by specified amount which should be
          * sufficient to keep disturbance on refined region until next regrid
          * of the level occurs.
          */
         hier::IntVector max_descriptor_ghosts(
            d_hierarchy->getPatchDescriptor()->getMaxGhostWidth(dim));

         /*
          * If the tag buffer value passed into this method is greater than the
          * current ghost width of the data that handles tag buffering, then
          * the call to resetTagBufferingData resets that data to have a ghost
          * width equal to the tag buffer.
          */

         if (tag_buffer > d_buf_tag_ghosts.max()) {
            resetTagBufferingData(tag_buffer);
         }

         /*
          * Create communication schedule for buffer tags on this level.
          */
         t_bdry_fill_tags_create->start();
         d_bdry_sched_tags[tag_ln] =
            d_bdry_fill_tags->createSchedule(tag_level, d_mb_tagger_strategy);
         t_bdry_fill_tags_create->stop();

         tag_level->allocatePatchData(d_buf_tag_indx, level_time);
         bufferTagsOnLevel(d_true_tag, tag_level, tag_buffer);

         tag_level->deallocatePatchData(d_buf_tag_indx);

         /*
          * We cannot leave this method with the tag buffering data having
          * ghosts greater than any other data managed by the patch descriptor,
          * so if that is the case, we reset it to the default value of 1.
          */

         if (needResetTagBuffer(max_descriptor_ghosts)) {
            resetTagBufferingData(1);
         }

         /*
          * Determine Boxes for new fine level.
          */
         findRefinementBoxes(new_box_level,
            tag_to_new,
            tag_ln);

         /*
          * Callback to user code to check user tags on the tag level.
          */
         d_tag_init_strategy->checkUserTagData(d_hierarchy,
            tag_ln,
            cycle,
            level_time,
            d_user_tag_indx);

         if (new_box_level && new_box_level->isInitialized()) {

            new_to_tag = &tag_to_new->getTranspose();

            if (d_check_proper_nesting) {
               /*
                * Check that the new BoxLevel nests inside the tag level.
                *
                * SAMRAI convention (my best understanding of it) supported
                * (or should have been supported) by grid generation:
                * - L0 must be equivalent to the domain.
                * - L1 must nest in L0 by the max ghost width in L1 index space,
                *   except where L1 touches physical boundary.
                * - L(n) must nest in L(n-1) by getProperNestingBuffer(n-1),
                *   except where L(n) touches physical boundary.
                */

               hier::IntVector required_nesting(dim);
               if (tag_ln > 0) {
                  required_nesting = hier::IntVector(dim,
                        d_hierarchy->getProperNestingBuffer(tag_ln));
               } else {
                  required_nesting = max_descriptor_ghosts;
               }

               bool locally_nests = false;
               const bool new_nests_in_tag =
                  d_blcu0.baseNestsInHead(
                     &locally_nests,
                     *new_box_level,
                     *d_hierarchy->getBoxLevel(tag_ln),
                     required_nesting,
                     zero_vector,
                     zero_vector,
                     &d_hierarchy->getGridGeometry()->getPeriodicDomainSearchTree());

               if (!new_nests_in_tag) {
                  std::shared_ptr<hier::BoxLevel> violator;
                  std::shared_ptr<hier::MappingConnector> new_to_violator;
                  t_compute_external_parts->start();
                  d_blcu0.computeExternalParts(
                     violator,
                     new_to_violator,
                     *new_to_tag,
                     hier::IntVector(dim, -d_hierarchy->getProperNestingBuffer(tag_ln)),
                     d_hierarchy->getGridGeometry()->getDomainSearchTree());
                  t_compute_external_parts->stop();

                  TBOX_ERROR(
                     "Internal library error: Failed to produce proper nesting."
                     << "GriddingAlgorithm::makeFinerLevel:\n"
                     << "tag_ln=" << tag_ln << ":\n"
                     << "new box_level does not properly nest\n"
                     << "in tag box_level by the required nesting buffer of "
                     << d_hierarchy->getProperNestingBuffer(tag_ln)
                     << ".\nLocal nestingness: " << locally_nests
                     << ".\nProper nesting violation with new_box_level of\n"
                     << new_box_level->format("N->", 2)
                     << "Proper nesting violation with tag box_level of\n"
                     << d_hierarchy->getBoxLevel(tag_ln)->format("F->", 2)
                     << "tag_to_new:\n" << tag_to_new->format("N->", 2)
                     << "new_to_tag:\n" << new_to_tag->format("N->", 2)
                     << "violator:\n" << violator->format("N->", 2)
                     << "new_to_violator:\n" << new_to_violator->format("N->", 2)
                     << std::endl);
               }
            }

         }

         /*
          * Deallocate algorithm tag arrays and schedule--no longer needed.
          */
         tag_level->deallocatePatchData(d_boolean_tag_indx);
         d_bdry_sched_tags[tag_ln].reset();

      } else { /* do_tagging == false */

         /*
          * If tagging is not necessary (do_tagging = false) we simply
          * need to access the level boxes, either from a dumped file or
          * from user-supplied refine boxes, and load balance them before
          * constructing the finer level.
          */
         bool remove_old_fine_level = false;
         readLevelBoxes(new_box_level,
            tag_to_new,
            tag_ln,
            level_time,
            cycle,
            remove_old_fine_level);

         new_to_tag = &tag_to_new->getTranspose();

         /*
          * Check for user-specified boxes that violate nesting requirements.
          */

         if (d_check_nonnesting_user_boxes != 'i') {
            checkNonnestingUserBoxes(
               *new_to_tag,
               new_to_tag->getRatio() * d_hierarchy->getProperNestingBuffer(tag_ln));
         }

         if (d_check_boundary_proximity_violation != 'i') {
            checkBoundaryProximityViolation(tag_ln, *new_box_level);
         }

      } /* do_tagging == false */

      if (new_box_level && new_box_level->isInitialized()) {

         /*
          * If we made a new_box_level, proceed to make a
          * PatchLevel from it.
          */

         // Bridge for new<==>new.
         if (d_barrier_and_time) {
            t_bridge_new_to_new->barrierAndStart();
         }
         std::shared_ptr<hier::Connector> new_to_new;
         d_oca.bridgeWithNesting(
            new_to_new,
            *new_to_tag,
            *tag_to_new,
            zero_vector,
            zero_vector,
            d_hierarchy->getRequiredConnectorWidth(new_ln, new_ln, true),
            false);
         if (d_barrier_and_time) {
            t_bridge_new_to_new->barrierAndStop();
         }

         new_box_level->cacheConnector(new_to_new);
         tag_level->cacheConnector(tag_to_new);

         if (d_check_overlapping_patches != 'i') {
            checkOverlappingPatches(*new_to_new);
         }

         d_tag_init_strategy->processHierarchyBeforeAddingNewLevel(d_hierarchy,
            new_ln,
            new_box_level);

         d_hierarchy->makeNewPatchLevel(new_ln, new_box_level);

         d_hierarchy->getGridGeometry()->adjustMultiblockPatchLevelBoundaries(
            *d_hierarchy->getPatchLevel(new_ln));

         d_tag_init_strategy->initializeLevelData(d_hierarchy,
            new_ln,
            level_time,
            d_hierarchy->levelCanBeRefined(new_ln),
            initial_cycle);

         t_reset_hier->barrierAndStart();
         d_tag_init_strategy->resetHierarchyConfiguration(d_hierarchy,
            new_ln,
            new_ln);
         t_reset_hier->stop();

         /*
          * Optionally refine and store user tag values on the new level.
          */
         if (new_ln > 0 && do_tagging && d_save_tag_data) {

            std::shared_ptr<hier::PatchLevel> new_level(
               d_hierarchy->getPatchLevel(new_ln));

            new_level->allocatePatchData(d_saved_tag_indx);

            d_saved_tags_sched =
               d_fill_saved_tags->createSchedule(
                  std::shared_ptr<xfer::PatchLevelFillPattern>(
                     new xfer::PatchLevelInteriorFillPattern()),
                  new_level,
                  std::shared_ptr<hier::PatchLevel>(),
                  new_ln-1,
                  d_hierarchy);

            d_saved_tags_sched->fillData(level_time, false);

            d_tag_init_strategy->checkNewLevelTagData(d_hierarchy,
               new_ln,
               d_saved_tag_indx);

         }

         if (d_log_metadata_statistics) {
            d_hierarchy->logMetadataStatistics("makeFinerLevel",
               d_hierarchy->getFinestLevelNumber(), cycle, level_time, false, true);
         }
      }

      if (do_tagging) {
         tag_level->deallocatePatchData(d_user_tag_indx);
      }

      d_base_ln = -1;

   }  // if level cannot be refined, the routine drops through...

   if (d_barrier_and_time) {
      t_make_finer->stop();
   }

}

/*
 *************************************************************************
 *
 * Regrid each level in the hierarchy which is finer than the
 * specified level.  If the regridding procedure employs time
 * integration, we perform any pre-processing necessary to regrid the
 * levels.  Then, each level finer than the specified level is
 * regridded from fine to coarse.  The recursive regridding procedure
 * is performed by the function regridFinerLevel().  Finally, after
 * the new hierarchy configuration is set, the application-specific
 * operations for resetting hierarchy-dependent infomation is called.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::regridAllFinerLevels(
   const int level_number,
   const std::vector<int>& tag_buffer,
   const int cycle,
   const double level_time,
   const std::vector<double>& regrid_start_time,
   const bool level_is_coarsest_sync_level)
{
   SAMRAI_CALI_CXX_MARK_FUNCTION;
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= d_hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(d_hierarchy->getPatchLevel(level_number));
   TBOX_ASSERT(static_cast<int>(tag_buffer.size()) >= level_number + 1);
#ifdef DEBUG_CHECK_ASSERTIONS
   int array_size = static_cast<int>(tag_buffer.size());
   for (int i = 0; i < array_size; ++i) {
      TBOX_ASSERT(tag_buffer[i] >= 0);
   }
#endif

   if (d_barrier_and_time) {
      t_regrid_all_finer->barrierAndStart();
   }

   if (d_tag_to_cluster_width.empty()) {
      computeTagToClusterWidths();
   }

   if (d_hierarchy->levelCanBeRefined(level_number)) {

      if (d_print_steps) {
         tbox::plog
         << "GriddingAlgorithm::regridAllFinerLevels: regridding finer than "
         << level_number << std::endl;
      }

      /*
       * d_base_ln is used by private methods during regrid.
       */
      d_base_ln = level_number;

      t_process_error->start();
      /*
       * Perform pre-processing of error estimation data, if
       * appropriate.
       */
      if (d_tag_init_strategy->usesTimeIntegration(cycle, level_time)) {
         for (int ln = level_number;
              ln <= d_hierarchy->getFinestLevelNumber(); ++ln) {
            if (d_hierarchy->levelCanBeRefined(ln)) {
               bool initial_time = false;
               double level_regrid_start_time = 0.;
               if (static_cast<int>(regrid_start_time.size()) < ln + 1) {
                  TBOX_ERROR("GriddingAlgorithm::regridAllFinerLevels()...\n"
                     << "no regrid_start_time specified for level " << ln
                     << "." << std::endl);
               } else {
                  level_regrid_start_time = regrid_start_time[ln];
               }

               d_tag_init_strategy->preprocessErrorEstimation(d_hierarchy,
                  ln,
                  cycle,
                  level_time,
                  level_regrid_start_time,
                  initial_time);
            }
         }
      }
      t_process_error->stop();

      /*
       * Recursively regrid each finer level.
       */
      const int finest_level_not_regridded = level_number;
      regridFinerLevel(
         level_number,
         level_time,
         cycle,
         finest_level_not_regridded,
         level_is_coarsest_sync_level,
         tag_buffer,
         regrid_start_time);

      /*
       * Invoke application-specific routines to reset information for those
       * levels which have been modified.
       */

      SAMRAI_CALI_MARK_BEGIN("regrid");
      if (d_hierarchy->getFinestLevelNumber() >= (level_number + 1)) {
         if (d_barrier_and_time) {
            t_reset_hier->barrierAndStart();
         }
         d_tag_init_strategy->
         resetHierarchyConfiguration(d_hierarchy,
            level_number + 1,
            d_hierarchy->getFinestLevelNumber());
         if (d_barrier_and_time) {
            t_reset_hier->stop();
         }
      }
      SAMRAI_CALI_MARK_END("regrid");

      d_base_ln = -1;

      if (d_print_steps) {
         tbox::plog
         << "GriddingAlgorithm::regridAllFinerLevels: regridded finer than "
         << level_number << std::endl;
      }

   } //  if level cannot be refined, the routine drops through...
   else {
      if (d_print_steps) {
         tbox::plog
         << "GriddingAlgorithm::regridAllFinerLevels: level "
         << level_number << " cannot be refined." << std::endl;
      }
   }

#ifdef GA_RECORD_STATS
   // Verified that this does not use much time.
   recordStatistics(level_time);
#endif

   if (d_barrier_and_time) {
      t_regrid_all_finer->stop();
   }

}

/*
 *************************************************************************
 *
 * Recursively, regrid each AMR hierarchy level finer than the
 * specified level (indexed by tag_ln).  The process is as follows:
 *
 * (1) Initialize tags to false on the level.
 *
 * (2) If a finer level exists, set tag to true on level for each cell
 * that is refined.
 *
 * (3) Tag cells for refinement on level by applying application-
 * specific error estimation routines.
 *
 * (4) If a finer level exists, invoke process recursively (i.e.,
 * invoke step 1 on next finer level).
 *
 * (5) (Note we have popped out of recursion at this point).  Buffer
 * true tags on current level to keep disturbances on fine grids until
 * regridding occurs next.
 *
 * (6) If level tag_ln+2 exists, tag under its boxes to ensure level
 * tag_ln+1 properly nests it.
 *
 * (7) Determine box configuration for new finer level, by calling
 * findRefinementBoxes() function.
 *
 * (8) If a finer level should exist in the hierarchy, create its
 * patches from the box description and initialize its data.  If
 * necessary, discard old level.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::regridFinerLevel(
   const int tag_ln,
   const double regrid_time,
   const int regrid_cycle,
   const int finest_level_not_regridded,
   const bool level_is_coarsest_sync_level,
   const std::vector<int>& tag_buffer,
   const std::vector<double>& regrid_start_time)
{
   TBOX_ASSERT((tag_ln >= 0)
      && (tag_ln <= d_hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(d_hierarchy->getPatchLevel(tag_ln));
   TBOX_ASSERT(finest_level_not_regridded >= 0
      && finest_level_not_regridded <= tag_ln);
   TBOX_ASSERT(static_cast<int>(tag_buffer.size()) >= tag_ln + 1);
#ifdef DEBUG_CHECK_ASSERTIONS
   int array_size = static_cast<int>(tag_buffer.size());
   for (int i = 0; i < array_size; ++i) {
      TBOX_ASSERT(tag_buffer[i] >= 0);
   }
#endif

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::regridFinerLevel: entered with tag_ln = "
      << tag_ln << "\n";
   }

   if (d_hierarchy->levelCanBeRefined(tag_ln)) {

      int new_ln = tag_ln + 1;

      std::shared_ptr<hier::PatchLevel> tag_level(
         d_hierarchy->getPatchLevel(tag_ln));

      /*
       * Compute nesting data at tag_ln for use in constructing
       * level tag_ln+1;
       */
      computeProperNestingData(tag_ln, d_oca);

      /*
       * The boolean "do_tagging" specifies whether or not tagging will
       * be performed. This will be true except in two circumstances:
       *   1) only user supplied refine boxes are used
       *   2) the boxes are read from a previously dumped file.
       *
       * If either of these circumstances is true, tagging operations
       * are NOT necessary so do_tagging will be set to false.
       *
       * The old level is generally removed when regridding, but
       * some circumstances may warrant keeping the old level.  For
       * example, if the refine region has not changed, there is no
       * need to regenerate the finer level.  The boolean
       * "remove_old_fine_level" specifies if the old level should
       * be removed.
       */
      bool do_tagging = true;
      if (d_tag_init_strategy->refineUserBoxInputOnly(regrid_cycle,
             regrid_time)) {
         do_tagging = false;
      }
      bool remove_old_fine_level = true;

      std::shared_ptr<hier::BoxLevel> new_box_level;
      std::shared_ptr<hier::Connector> tag_to_new;

      /*
       * tag_to_finer is [tag_ln]->[tag_ln+2].
       *
       * These are declared in this scope, computed and cached if
       * [tag_ln+2] exists.  They are used to tag the footprint of
       * [tag_ln+2] on level tag_ln.  Later on, they are used to
       * bridge for [new_ln] <-> [tag_ln+2].
       */
      std::shared_ptr<hier::Connector> tag_to_finer;

      /*
       * Tag cells, determine refine boxes from tagged regions, and
       * load balance in preparation for constructing new refined level.
       */
      if (do_tagging) {
         /*
          * Tagging stuff have been factored out to shorten this method.
          */
        RANGE_PUSH("doTaggingBefore", 3);
         regridFinerLevel_doTaggingBeforeRecursiveRegrid(
            tag_ln,
            level_is_coarsest_sync_level,
            regrid_start_time,
            regrid_time,
            regrid_cycle);
         RANGE_POP;

         /*
          * Perform regridding recursively on finer levels, if appropriate.
          */
         if (d_hierarchy->finerLevelExists(tag_ln)
             && d_hierarchy->levelCanBeRefined(new_ln)) {

            if (d_print_steps) {
               tbox::plog
               << "GriddingAlgorithm::regridFinerLevel: recursing to tag_ln = "
               << new_ln << "\n";
            }

            regridFinerLevel(
               new_ln,
               regrid_time,
               regrid_cycle,
               finest_level_not_regridded,
               false,
               tag_buffer,
               regrid_start_time);

            if (d_print_steps) {
               tbox::plog
               << "GriddingAlgorithm::regridFinerLevel: recursion returned to tag_ln = "
               << tag_ln << "\n";
            }

         }

         /*
          * Tagging stuff have been factored out to shorten this method.
          */
         RANGE_PUSH("doTaggingAfter", 3);
         regridFinerLevel_doTaggingAfterRecursiveRegrid(
            tag_to_finer,
            tag_ln,
            tag_buffer,
            regrid_time);
         RANGE_POP;

         /*
          * Determine boxes containing cells on level with a true tag
          * value.
          */
         RANGE_PUSH("findRefinemetBoxes", 3);
         findRefinementBoxes(
            new_box_level,
            tag_to_new,
            tag_ln);
         RANGE_POP;

         d_tag_init_strategy->checkUserTagData(d_hierarchy,
            tag_ln,
            regrid_cycle,
            regrid_time,
            d_user_tag_indx);

         if (d_print_steps) {
            if (new_box_level && new_box_level->isInitialized()) {
               tbox::plog
               << "GriddingAlgorithm::regridFinerLevel: got inititalized new_box_level\n";
            } else {
               tbox::plog
               << "GriddingAlgorithm::regridFinerLevel: got un-inititalized new_box_level\n";
            }
         }

         /*
          * Deallocate tag arrays and schedule; no longer needed on current
          * level.
          */

         tag_level->deallocatePatchData(d_boolean_tag_indx);
         d_bdry_sched_tags[tag_ln].reset();

      } else { /* do_tagging == false */

         if (d_hierarchy->finerLevelExists(tag_ln)
             && d_hierarchy->levelCanBeRefined(new_ln)) {
            regridFinerLevel(
               new_ln,
               regrid_time,
               regrid_cycle,
               finest_level_not_regridded,
               false,
               tag_buffer,
               regrid_start_time);
         }

         /*
          * If tagging is not necessary (do_tagging = false) we simply
          * need to access the user-supplied refine boxes, and load
          * balance them before constructing the finer level.
          */
         readLevelBoxes(new_box_level,
            tag_to_new,
            tag_ln,
            regrid_time,
            regrid_cycle,
            remove_old_fine_level);

      } /* end do_tagging == false */

      /*
       * Make new finer level (new_ln) if necessary, or remove
       * next finer level if it is no longer needed.
       */

      RANGE_PUSH("new_box_level", 3);
      if (new_box_level && new_box_level->isInitialized()) {

         /*
          * Create the new PatchLevel from the new_box_level.
          */
         regridFinerLevel_createAndInstallNewLevel(
            tag_ln,
            regrid_time,
            tag_to_new,
            new_box_level);

         if (d_log_metadata_statistics) {
            // Don't log the coarse Connector, if the coarse level will be updated.
            d_hierarchy->logMetadataStatistics("regridFinerLevel",
               new_ln,
               regrid_cycle,
               regrid_time,
               new_ln < d_hierarchy->getFinestLevelNumber(),
               true);
         }

      } else {

         /*
          * The new level has no boxes, so we don't generate it.
          * Remove the pre-existing fine level if it exists.
          */

         if (d_hierarchy->finerLevelExists(tag_ln)
             && remove_old_fine_level) {
            d_tag_init_strategy->processLevelBeforeRemoval(
               d_hierarchy,
               new_ln,
               d_hierarchy->getPatchLevel(new_ln));
            d_hierarchy->removePatchLevel(new_ln);
         }

      } // if we are not re-regenerating level new_ln.
      RANGE_POP;

      if (do_tagging) {
         tag_level->deallocatePatchData(d_user_tag_indx);
      }


   } //  if level cannot be refined, the routine drops through...

}

/*
 *************************************************************************
 * Various tagging stuff done before recursively regridding a finer level.
 *************************************************************************
 */
void
GriddingAlgorithm::regridFinerLevel_doTaggingBeforeRecursiveRegrid(
   const int tag_ln,
   const bool level_is_coarsest_sync_level,
   const std::vector<double>& regrid_start_time,
   const double regrid_time,
   const int regrid_cycle)
{
   if (d_print_steps) {
      tbox::plog
      <<
      "GriddingAlgorithm::regridFinerLevel_doTaggingBeforeRecursiveRegrid: entered with tag_ln = "
      << tag_ln << "\n";
   }

   t_regrid_finer_do_tagging_before->start();

   const hier::IntVector& zero_vec(hier::IntVector::getZero(d_hierarchy->getDim()));

   const std::shared_ptr<hier::PatchLevel>& tag_level(
      d_hierarchy->getPatchLevel(tag_ln));

   /*
    * Create communication schedule for buffer tags and set tags to
    * false.
    */

   tag_level->allocatePatchData(d_user_tag_indx, regrid_time);
   fillTags(d_false_tag, tag_level, d_user_tag_indx);

   /*
    * Set tags to true for cells that currently cover next finer level.
    * Note that this is not needed for all regridding strategies.  But
    * knowledge of currently refined cells is generally useful to avoid
    * repeated refining and coarsening of cells near boundaries of
    * refined regions where error estimates may hover around the error
    * tolerance. For example, the regridding scheme may require that
    * the error in a cell that is currently refined fall below a
    * certain tolerance  (generally different than the tolerance to
    * refine the cell in the first place) before the cell will be
    * de-refined.
    */

   if (d_hierarchy->finerLevelExists(tag_ln)) {
      fillTagsFromBoxLevel(
         d_from_fine_pretag,
         tag_level,
         d_user_tag_indx,
         d_hierarchy->getPatchLevel(tag_ln)->findConnector(
            *d_hierarchy->getPatchLevel(tag_ln + 1),
            d_hierarchy->getRequiredConnectorWidth(tag_ln, tag_ln + 1, true),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            false),
         true,
         zero_vec,
         false);
   }

   /*
    * Determine cells needing refinement according to a specific
    * error estimation procedure and set to true.
    *
    * The "level_is_coarsest_sync_level" is provided as an argument
    * to this method.  Provide the additional check of whether the
    * level is not the coarsest level and that it is not a new level
    * in the hierarchy.  If all three conditions are true, the
    * "coarsest_sync_level" argument passed into the tagCells method
    * will be true.  Otherwise, it will be false.
    */

   bool coarsest_sync_level =
      level_is_coarsest_sync_level &&
      tag_ln > 0 &&
      tag_ln <= d_base_ln;

   bool initial_time = false;
   double level_regrid_start_time = 0.;
   if (d_tag_init_strategy->usesTimeIntegration(regrid_cycle, regrid_time)) {
      if (static_cast<int>(regrid_start_time.size()) < tag_ln + 1) {
         TBOX_ERROR("GriddingAlgorithm::regridFinerLevel_doTaggingBeforeRecursiveRegrid:\n"
            << "no regrid_start_time specified for level " << tag_ln
            << "." << std::endl);
      } else {
         level_regrid_start_time = regrid_start_time[tag_ln];
      }
   }

   if (d_barrier_and_time) {
      t_tag_cells_for_refinement->barrierAndStart();
   }
   d_tag_init_strategy->tagCellsForRefinement(d_hierarchy,
      tag_ln,
      regrid_cycle,
      regrid_time,
      d_user_tag_indx,
      initial_time,
      coarsest_sync_level,
      d_hierarchy->levelCanBeRefined(tag_ln),
      level_regrid_start_time);
   if (d_barrier_and_time) {
      t_tag_cells_for_refinement->stop();
   }

   /*
    * Check for user-tagged cells that violate proper nesting,
    * except if user specified that the violating tags be ignored.
    */
   if (d_check_nonrefined_tags != 'i') {
      checkNonrefinedTags(*tag_level, tag_ln, d_oca);
   }

   tag_level->allocatePatchData(d_boolean_tag_indx, regrid_time);
   setBooleanTagData(tag_level, false);

   t_regrid_finer_do_tagging_before->stop();
}

/*
 *************************************************************************
 * Various tagging stuff done after recursively regridding a finer level.
 *************************************************************************
 */
void
GriddingAlgorithm::regridFinerLevel_doTaggingAfterRecursiveRegrid(
   std::shared_ptr<hier::Connector>& tag_to_finer,
   const int tag_ln,
   const std::vector<int>& tag_buffer,
   double regrid_time)
{
   if (d_print_steps) {
      tbox::plog
      <<
      "GriddingAlgorithm::regridFinerLevel_doTaggingAfterRecursiveRegrid: entered with tag_ln = "
      << tag_ln << "\n";
   }

   t_regrid_finer_do_tagging_after->start();

   const int new_ln = tag_ln + 1;
   const std::shared_ptr<hier::PatchLevel>& tag_level(
      d_hierarchy->getPatchLevel(tag_ln));

   /*
    * Buffer true tagged cells by specified amount which should be
    * sufficient to keep disturbance on refined region until next
    * regrid of the level occurs.
    */

   hier::IntVector max_descriptor_ghosts(
      d_hierarchy->getPatchDescriptor()->getMaxGhostWidth(d_hierarchy->getDim()));

   /*
    * If the tag buffer value passed into this method is greater than the
    * current ghost width of the data that handles tag buffering, then the
    * call to resetTagBufferingData resets that data to have a ghost width
    * equal to the tag buffer.
    */
   if (tag_buffer[tag_ln] > d_buf_tag_ghosts.max()) {
      resetTagBufferingData(tag_buffer[tag_ln]);
   }

   t_bdry_fill_tags_create->start();
   d_bdry_sched_tags[tag_ln] =
      d_bdry_fill_tags->createSchedule(tag_level, d_mb_tagger_strategy);
   t_bdry_fill_tags_create->stop();

   tag_level->allocatePatchData(d_buf_tag_indx, regrid_time);
   bufferTagsOnLevel(d_true_tag, tag_level, tag_buffer[tag_ln]);

   /*
    * We cannot leave this method with the tag buffering data having ghosts
    * greater than any other data managed by the patch descriptor, so if that
    * is the case, we reset it to the default value of 1.
    */
   if (needResetTagBuffer(max_descriptor_ghosts)) {
      resetTagBufferingData(1);
   }

   if (d_hierarchy->finerLevelExists(new_ln)) {

      /*
       * Add tags to ensure new_ln properly nests level
       * new_ln+1.  This means tagging under level new_ln+1
       * plus an appropriate nesting buffer equal to how
       * much new_ln+1 should nest in new_ln.
       *
       * To determine where to tag, we compute Connector
       * [tag_ln] ---> [new_ln+1], aka tag_to_finer.  Then use
       * its neighbor data to determine where to tag.
       */

      if (d_barrier_and_time) {
         t_second_finer_tagging->start();
      }

      const hier::Connector& tag_to_old =
         d_hierarchy->getPatchLevel(tag_ln)->findConnectorWithTranspose(
            *d_hierarchy->getPatchLevel(new_ln),
            d_hierarchy->getRequiredConnectorWidth(tag_ln, new_ln, true),
            d_hierarchy->getRequiredConnectorWidth(new_ln, tag_ln),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            false);
      const hier::Connector& old_to_finer =
         d_hierarchy->getPatchLevel(new_ln)->findConnectorWithTranspose(
            *d_hierarchy->getPatchLevel(new_ln + 1),
            d_hierarchy->getRequiredConnectorWidth(new_ln, new_ln + 1, true),
            d_hierarchy->getRequiredConnectorWidth(new_ln + 1, new_ln),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            false);
      d_oca.bridge(
         tag_to_finer,
         tag_to_old,
         old_to_finer,
         true);

      // Nesting buffer in resolution of level new_ln+1.
      const hier::IntVector nesting_buffer =
         d_hierarchy->getRatioToCoarserLevel(new_ln + 1)
         * d_hierarchy->getProperNestingBuffer(tag_ln + 1);

      TBOX_ASSERT(
         tag_to_finer->getConnectorWidth()
         * d_hierarchy->getRatioToCoarserLevel(tag_ln + 1)
         * d_hierarchy->getRatioToCoarserLevel(new_ln + 1) >= nesting_buffer);

      // Add periodic relationships in tag_to_finer.
      hier::BoxLevel dummy_finer_box_level =
         old_to_finer.getTranspose().getBase();
      d_blcu.addPeriodicImagesAndRelationships(
         dummy_finer_box_level,
         tag_to_finer->getTranspose(),
         d_hierarchy->getGridGeometry()->getDomainSearchTree(),
         d_hierarchy->getPatchLevel(tag_ln)->findConnector(
            *d_hierarchy->getPatchLevel(tag_ln),
            d_hierarchy->getRequiredConnectorWidth(tag_ln, tag_ln, true),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            false));

      fillTagsFromBoxLevel(
         d_from_fine_tag,
         tag_level,
         d_user_tag_indx,
         *tag_to_finer,
         true,
         nesting_buffer,
         true);
      /*
       * New tags may have been added to user tags via the above step.
       * Here we set algorithmic tags at those new locations without
       * overwriting any algorithmic tags that already exist.
       */
      setBooleanTagData(tag_level, true);

      if (d_barrier_and_time) {
         t_second_finer_tagging->stop();
      }

   } // End tagging under level new_ln+1.

   tag_level->deallocatePatchData(d_buf_tag_indx);

   t_regrid_finer_do_tagging_after->stop();
}

/*
 *************************************************************************
 * After the metadata describing the new level is computed,
 * this method creates and installs new PatchLevel in the hierarchy.
 *************************************************************************
 */
void
GriddingAlgorithm::regridFinerLevel_createAndInstallNewLevel(
   const int tag_ln,
   const double regrid_time,
   std::shared_ptr<hier::Connector>& tag_to_new,
   std::shared_ptr<hier::BoxLevel> new_box_level)
{
   TBOX_ASSERT(tag_to_new && tag_to_new->hasTranspose());
   TBOX_ASSERT(new_box_level);

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::regridFinerLevel_createAndInstallNewLevel: starting\n";
   }

   t_regrid_finer_create_and_install->start();

   hier::Connector& new_to_tag = tag_to_new->getTranspose();

   /*
    * Compute self-Connector for the new level.
    */

   const tbox::Dimension& dim = d_hierarchy->getDim();

   const int new_ln = tag_ln + 1;
   const std::shared_ptr<hier::PatchLevel>& tag_level(
      d_hierarchy->getPatchLevel(tag_ln));

   const hier::IntVector& zero_vector(hier::IntVector::getZero(dim));

   std::shared_ptr<hier::Connector> new_to_new;

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::regridFinerLevel_createAndInstallNewLevel: bridging for new<==>new\n";
   }

   if (d_barrier_and_time) {
      t_bridge_new_to_new->barrierAndStart();
   }
   d_oca.bridgeWithNesting(
      new_to_new,
      new_to_tag,
      *tag_to_new,
      zero_vector,
      zero_vector,
      d_hierarchy->getRequiredConnectorWidth(new_ln, new_ln, true),
      false);
   if (d_barrier_and_time) {
      t_bridge_new_to_new->barrierAndStop();
   }

   TBOX_ASSERT(new_to_new->getConnectorWidth() ==
      d_hierarchy->getRequiredConnectorWidth(new_ln, new_ln));

   if (d_check_overlapping_patches != 'i') {
      checkOverlappingPatches(*new_to_new);
   }

   if (d_barrier_and_time) {
      new_box_level->getMPI().Barrier();
   }
   t_regrid_finer_create->start();

   /*
    * Either remove pre-existing fine level from hierarchy and make
    * a new level, or just make a new fine level for hierarchy.
    */

   /*
    * Save references to old objects before hierarchy removes them.
    * We need this while installing new objects.
    */

   std::shared_ptr<hier::PatchLevel> old_fine_level;
   std::shared_ptr<const hier::BoxLevel> old_box_level;
   const hier::Connector* old_to_tag = 0;

   hier::IntVector ratio(tag_level->getRatioToLevelZero()
                         * d_hierarchy->getRatioToCoarserLevel(new_ln));

   if (d_hierarchy->finerLevelExists(tag_ln)) {

      old_box_level = d_hierarchy->getBoxLevel(new_ln);
      old_to_tag =
         &d_hierarchy->getPatchLevel(new_ln)->findConnectorWithTranspose(
            *d_hierarchy->getPatchLevel(tag_ln),
            d_hierarchy->getRequiredConnectorWidth(new_ln, tag_ln, true),
            d_hierarchy->getRequiredConnectorWidth(tag_ln, new_ln),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            false);

      old_fine_level = d_hierarchy->getPatchLevel(new_ln);
      TBOX_ASSERT(ratio == old_fine_level->getRatioToLevelZero());

   }


   if (d_check_proper_nesting && d_hierarchy->levelExists(new_ln + 1)) {
      /*
       * Check that the new_box_level nests the next finer
       * level (new_ln+1).
       */

      hier::IntVector required_nesting(dim, d_hierarchy->getProperNestingBuffer(new_ln));
      required_nesting *= d_hierarchy->getRatioToCoarserLevel(new_ln + 1);

      bool locally_nests = false;
      const bool finer_nests_in_new =
         d_blcu.baseNestsInHead(
            &locally_nests,
            *d_hierarchy->getBoxLevel(new_ln + 1),
            *new_box_level,
            required_nesting,
            zero_vector,
            zero_vector,
            &d_hierarchy->getGridGeometry()->getPeriodicDomainSearchTree());

      if (!finer_nests_in_new) {

         tbox::perr
            << "GriddingAlgorithm::regridFinerLevel_createAndInstallNewLevel: new box_level\n"
            << "at ln=" << new_ln
            << " does not properly nest\n"
            << "existing finer box_level at ln="
            << new_ln + 1
            << " by the required nesting buffer of "
            << required_nesting << " in fine resolution.\n"
            << "Local nestingness: " << locally_nests
            << ".\nWriting BoxLevels out to log file."
            << std::endl;
         tbox::plog
            << "Proper nesting violation with new_box_level of\n"
            << new_box_level->format("N->", 2)
            << "Proper nesting violation with finer box_level of\n"
            << d_hierarchy->getBoxLevel(new_ln + 1)->format("F->", 2);

         std::shared_ptr<hier::BoxLevel> external;
         std::shared_ptr<hier::MappingConnector> finer_to_external;
         std::shared_ptr<hier::Connector> finer_to_new;
         d_oca.findOverlaps(finer_to_new,
                            *d_hierarchy->getBoxLevel(new_ln + 1),
                            *new_box_level,
                            required_nesting);
         tbox::plog << "Finer to new:\n" << finer_to_new->format("FN->", 3);
         d_blcu.computeExternalParts(
            external,
            finer_to_external,
            *finer_to_new,
            -required_nesting,
            d_hierarchy->getGridGeometry()->getDomainSearchTree());
         tbox::plog << "External parts:\n" << finer_to_external->format("FE->", 3);

         TBOX_ERROR(
            "Internal library error: Failed to produce proper nesting."
            << std::endl);

      } /* !finer_nests_in_new */

   } /* d_check_proper_nesting */

   /*
    * Cache Connectors for new level.
    */
   new_box_level->cacheConnector(new_to_new);
   tag_level->cacheConnector(tag_to_new);

   std::shared_ptr<hier::Connector> old_to_new;
   if (old_box_level) {

      /*
       * Connect old to new by bridging.
       *
       * Cache these Connectors for use when creating schedules to
       * transfer data from old to new.
       */

      if (d_print_steps) {
         tbox::plog
         <<
         "GriddingAlgorithm::regridFinerLevel_createAndInstallNewLevel: bridging for new<==>old\n";
      }

      if (d_barrier_and_time) {
         t_bridge_new_to_old->barrierAndStart();
      }
      d_oca.bridgeWithNesting(
         old_to_new,
         *old_to_tag,
         d_hierarchy->getPatchLevel(tag_ln)->getBoxLevel()->findConnectorWithTranspose(
            *new_box_level,
            d_hierarchy->getRequiredConnectorWidth(tag_ln, tag_ln + 1, true),
            d_hierarchy->getRequiredConnectorWidth(tag_ln + 1, tag_ln),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            false),
         zero_vector,
         zero_vector,
         d_hierarchy->getRequiredConnectorWidth(new_ln, new_ln, true),
         true);
      if (d_barrier_and_time) {
         t_bridge_new_to_old->barrierAndStop();
      }

      old_fine_level->cacheConnector(old_to_new);

   }

   if (d_hierarchy->levelExists(new_ln + 1)) {
      /*
       * There is a level finer than new_ln.  Connect the new level to
       * the finer level.
       */

      if (d_print_steps) {
         tbox::plog
         <<
         "GriddingAlgorithm::regridFinerLevel_createAndInstallNewLevel: bridging for new<==>finer\n";
      }

      const hier::Connector& old_to_finer =
         d_hierarchy->getPatchLevel(new_ln)->getBoxLevel()->findConnectorWithTranspose(
            *d_hierarchy->getPatchLevel(new_ln + 1)->getBoxLevel(),
            d_hierarchy->getRequiredConnectorWidth(new_ln, new_ln + 1),
            d_hierarchy->getRequiredConnectorWidth(new_ln + 1, new_ln),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE);
      std::shared_ptr<hier::Connector> new_to_finer;

      if (d_barrier_and_time) {
         t_bridge_new_to_finer->barrierAndStart();
      }
      d_oca.bridgeWithNesting(
         new_to_finer,
         old_to_new->getTranspose(),
         old_to_finer,
         -hier::IntVector::getOne(dim),
         zero_vector,
         d_hierarchy->getRequiredConnectorWidth(new_ln, new_ln + 1, true),
         true);
      if (d_barrier_and_time) {
         t_bridge_new_to_finer->barrierAndStop();
      }

#ifdef DEBUG_CHECK_ASSERTIONS
      hier::Connector& finer_to_new = new_to_finer->getTranspose();
#endif

      TBOX_ASSERT(
         new_to_finer->getConnectorWidth() ==
         d_hierarchy->getRequiredConnectorWidth(new_ln, new_ln + 1));
      TBOX_ASSERT(
         finer_to_new.getConnectorWidth() ==
         d_hierarchy->getRequiredConnectorWidth(new_ln + 1, new_ln));

#ifdef DEBUG_CHECK_ASSERTIONS
      new_to_finer->assertOverlapCorrectness(false, true, true);
      finer_to_new.assertOverlapCorrectness(false, true, true);
#endif

      std::shared_ptr<hier::PatchLevel> finer_level(
         d_hierarchy->getPatchLevel(new_ln + 1));

      new_box_level->cacheConnector(new_to_finer);
   }

   d_tag_init_strategy->processHierarchyBeforeAddingNewLevel(d_hierarchy,
      new_ln,
      new_box_level);

   if (old_box_level) {
      d_hierarchy->removePatchLevel(new_ln);
   }
   d_hierarchy->makeNewPatchLevel(new_ln, new_box_level);

   d_hierarchy->getGridGeometry()->adjustMultiblockPatchLevelBoundaries(
      *d_hierarchy->getPatchLevel(new_ln));

   t_regrid_finer_create->stop();

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::regridFinerLevel_createAndInstallNewLevel: initializing level data\n";
   }

   if (d_barrier_and_time) {
      t_initialize_level_data->barrierAndStart();
   }
   // "false" argument": const bool initial_time = false;
   d_tag_init_strategy->initializeLevelData(d_hierarchy,
      new_ln,
      regrid_time,
      d_hierarchy->levelCanBeRefined(new_ln),
      false,
      old_fine_level);
   if (d_barrier_and_time) {
      t_initialize_level_data->barrierAndStop();
   }

   if (new_ln > 0 && d_save_tag_data) {
      std::shared_ptr<hier::PatchLevel> new_level(
         d_hierarchy->getPatchLevel(new_ln));

      new_level->allocatePatchData(d_saved_tag_indx);

      d_saved_tags_sched =
         d_fill_saved_tags->createSchedule(
            std::shared_ptr<xfer::PatchLevelFillPattern>(
               new xfer::PatchLevelInteriorFillPattern()),
            d_hierarchy->getPatchLevel(new_ln),
            std::shared_ptr<hier::PatchLevel>(),
            new_ln-1,
            d_hierarchy);

      d_saved_tags_sched->fillData(regrid_time, false);

      d_tag_init_strategy->checkNewLevelTagData(d_hierarchy,
         new_ln,
         d_saved_tag_indx);

   }

   /*
    * Destroy old patch level, if such a level existed prior to regrid.
    */
   old_fine_level.reset();

   t_regrid_finer_create_and_install->stop();

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::regridFinerLevel_createAndInstallNewLevel: returning\n";
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void
GriddingAlgorithm::computeTagToClusterWidths()
{
   TBOX_ASSERT(d_tag_to_cluster_width.empty());    // Never recompute.

   const tbox::Dimension& dim = d_hierarchy->getDim();

   d_tag_to_cluster_width.resize(d_hierarchy->getMaxNumberOfLevels() - 1,
      hier::IntVector::getZero(dim));

   const size_t nblocks = d_hierarchy->getGridGeometry()->getNumberBlocks();

   for (int ln = d_hierarchy->getMaxNumberOfLevels() - 2; ln >= 0; --ln) {
      /*
       * Construct list of boxes covering the true tags on the level.
       * Note that box list will be contained in the bounding box
       * but will not be contained in the list of proper nesting boxes,
       * in general.  So we intersect the box list against the list of
       * nesting boxes.  Note that this may produce boxes which are too
       * small.  Thus, boxes are regrown later.
       */

      hier::IntVector smallest_patch(dim);
      hier::IntVector smallest_box_to_refine(dim, 0, nblocks);
      hier::IntVector largest_patch(dim);
      hier::IntVector extend_ghosts(dim);
      // "true" argument: for_building_finer level = true
      getGriddingParameters(smallest_patch,
         smallest_box_to_refine,
         largest_patch,
         extend_ghosts,
         ln + 1,
         true);
      const hier::IntVector extend_ghosts_in_tag_space =
         hier::IntVector::ceilingDivide(extend_ghosts,
            d_hierarchy->getRatioToCoarserLevel(ln + 1));

      /*
       * Compute the width for tag<==>cluster.  This width be wide enough to
       * guarantee completeness for tag<==>new when we massage the
       * cluster into the new level.  In the massage step, we grow the new boxes.
       * The width of tag<==>cluster must be big enough to see new overlaps
       * generated by the growths.
       */
      d_tag_to_cluster_width[ln] =
         d_hierarchy->getRequiredConnectorWidth(ln, ln + 1);

      // For width of d_tag_to_cluster_width[ln+1] in bridge new<==>tag<==>new
      if (ln + 1 < static_cast<int>(d_tag_to_cluster_width.size())) {
         d_tag_to_cluster_width[ln].max(
            hier::IntVector::ceilingDivide(d_tag_to_cluster_width[ln + 1],
               d_hierarchy->getRatioToCoarserLevel(ln + 1)));
      }

      if (d_extend_to_domain_boundary) {
         // For extending boxes to domain boundary by amount of extend_ghosts_in_tag_space.
         d_tag_to_cluster_width[ln] += extend_ghosts_in_tag_space;
      }
      // For growing within domain by amount of smallest_box_to_refine.
      d_tag_to_cluster_width[ln] += smallest_box_to_refine;
   }

   d_connector_width_requestor.setTagToClusterWidth(d_tag_to_cluster_width);

   // Commit to computing the widths required by the hierarchy.
   d_hierarchy->getRequiredConnectorWidth(0, 0, true);
}

/*
 *************************************************************************
 * Check for boundary proximity violations.  Boxes may not be within
 * extend_ghosts of a physical boundary without touching the physical
 * boundary.
 *************************************************************************
 */
size_t
GriddingAlgorithm::checkBoundaryProximityViolation(
   const hier::BoxLevel& box_level,
   const hier::IntVector& extend_ghosts) const
{
   /*
    * 1. Compute the boundary regions of the boxes.
    *    a. Grow temporary boxes by the max ghost width.
    *    b. Remove orig boxes from the grown boxes to get ghost regions.
    *    c. Remove domain from the ghost region.
    * 2. Physical boundary regions which should not
    *    be less wide than the max ghost width.  If they are,
    *    it means they are partially inside the domain.
    */

   const hier::BaseGridGeometry& grid_geometry(
      *d_hierarchy->getGridGeometry());

   const hier::BoxContainer& periodic_domain_search_tree(
      grid_geometry.getPeriodicDomainSearchTree());

   hier::BoxContainer refined_periodic_domain_search_tree(
      periodic_domain_search_tree);
   refined_periodic_domain_search_tree.refine(box_level.getRefinementRatio());
   refined_periodic_domain_search_tree.makeTree(&grid_geometry);

   size_t nerr(0);

   const tbox::Dimension& dim = d_hierarchy->getDim();
   for (hier::RealBoxConstIterator bi(box_level.getBoxes().realBegin());
        bi != box_level.getBoxes().realEnd(); ++bi) {

      hier::BoxContainer external_parts(*bi);
      external_parts.grow(extend_ghosts);
      external_parts.removeIntersections(
         box_level.getRefinementRatio(),
         refined_periodic_domain_search_tree);

      for (hier::BoxContainer::iterator bli = external_parts.begin();
           bli != external_parts.end();
           ++bli) {
         hier::IntVector extend(
            extend_ghosts.getBlockVector(bli->getBlockId()));
         hier::IntVector leftover_size(bli->numberCells());
         for (int d = 0; d < dim.getValue(); ++d) {
            if (leftover_size(d) != 0 && leftover_size(d) < extend(d)) {
               ++nerr;
               TBOX_WARNING("GriddingAlgorithm::checkBoundaryProximityViolation:\n"
                  << "User-specified box (refined) " << *bi
                  << " violates boundary proximity.\n"
                  << "In direction " << d << ", it is "
                  << extend(d) - leftover_size(d)
                  << " cells from a physical domain boundary.\n"
                  << "All boxes must be at least " << extend
                  << " from physical boundaries or touching the physical boundary."
                  << std::endl);
            }
         }
      }

   }

   return nerr;
}

/*
 *******************************************************************
 * Check domain boxes for violations of user constraints.
 *******************************************************************
 */
void
GriddingAlgorithm::checkDomainBoxes(const hier::BoxContainer& domain_boxes) const {
   const tbox::Dimension& dim = d_hierarchy->getDim();

   const size_t nblocks = d_hierarchy->getGridGeometry()->getNumberBlocks();

   hier::IntVector smallest_patch(dim);
   hier::IntVector largest_patch(dim);
   hier::IntVector extend_ghosts(dim);
   {
      hier::IntVector smallest_box_to_refine(dim, 0, nblocks);
      // "false" argument: for_building_finer level = false
      getGriddingParameters(
         smallest_patch,
         smallest_box_to_refine,
         largest_patch,
         extend_ghosts,
         0,
         false);
   }

   /*
    * Check minimum size violations.
    */
   int i = 0;
   for (hier::BoxContainer::const_iterator itr = domain_boxes.begin();
        itr != domain_boxes.end(); ++itr, ++i) {

      hier::Box test_box = *itr;
      for (tbox::Dimension::dir_t dir = 0; dir < dim.getValue(); ++dir) {

         if (test_box.numberCells(dir) < smallest_patch(dir)) {

            int error_coarsen_ratio =
               d_tag_init_strategy->getErrorCoarsenRatio();
            if (error_coarsen_ratio > 1) {
               TBOX_ERROR(
                  d_object_name << ": " << "\ndomain Box " << i << ", " << test_box
                                << ", violates the minimum patch size constraints."
                                << "\nVerify that boxes are larger than"
                                << "the maximum ghost width and/or"
                                << "\nthe specified minimum patch size."
                                << "\nNOTE: to assure the constraints are"
                                << "properly enforced during coarsening for"
                                << "\nerror computation, the minimum patch"
                                << "size is the smallest patch size multiplied"
                                << "\nby the error coarsen ratio, which is "
                                << error_coarsen_ratio
                                << " in this case."
                                << std::endl);
            } else {
               TBOX_ERROR(
                  d_object_name << ": "
                                << "\ndomain Box " << i << ", " << test_box
                                << ", violates the minimum patch size constraints."
                                << "\nVerify that boxes are larger than"
                                << " the maximum ghost width and/or"
                                << "\nthe specified minimum patch size, "
                                << smallest_patch << "."
                                << std::endl);
            }
         }
      }
   }

   /*
    * Check for overlapping boxes.
    */
   if (domain_boxes.boxesIntersect()) {
      TBOX_ERROR(d_object_name << ":  "
                               << "Boxes specified for coarsest level "
                               << "contain intersections with each other!"
                               << std::endl);
   }

   /*
    * Check for violations of implementation of TagAndInitStrategy.
    */
   if ((d_hierarchy->getMaxNumberOfLevels() > 1)
       && (!d_tag_init_strategy->coarsestLevelBoxesOK(domain_boxes))) {
      TBOX_ERROR(d_object_name << ":  "
                               << "level gridding strategy encountered"
                               << " a problem with the domain boxes!"
                               << std::endl);
   }
}

/*
 *******************************************************************
 * Check for non-nesting user-specified boxes.
 *******************************************************************
 */
void
GriddingAlgorithm::checkNonnestingUserBoxes(
   const hier::Connector& new_to_tag,
   const hier::IntVector& nesting_buffer) const
{

   const hier::BoxLevel& new_box_level(new_to_tag.getBase());

   std::shared_ptr<hier::BoxLevel> violating_parts;
   std::shared_ptr<hier::MappingConnector> new_to_violating_parts;

   d_blcu.computeExternalParts(
      violating_parts,
      new_to_violating_parts,
      new_to_tag,
      -nesting_buffer,
      d_hierarchy->getGridGeometry()->getDomainSearchTree());

   if (violating_parts->getGlobalNumberOfBoxes() > 0) {

      tbox::perr
      << "GriddingAlgorihtm::checkNonnestingUserBoxes: user-specified refinement boxes\n"
      << "violates nesting requirement.  Diagnostics will be\n"
      << "writen to log files." << std::endl;
      const std::string left_margin("ERR: ");
      tbox::plog
      << left_margin << "Tag BoxLevel:\n" << new_to_tag.getHead().format(left_margin, 2)
      << left_margin << "User-specified boxes:\n" << new_box_level.format(left_margin, 2)
      << left_margin << "Violating parts:\n" << violating_parts->format(left_margin, 2)
      << left_margin << "User-specified boxes and their violating parts:\n"
      << new_to_violating_parts->format(left_margin, 2);

      if (d_check_nonnesting_user_boxes == 'e') {
         TBOX_ERROR("Exiting due to above error" << std::endl);
      }
      if (d_check_nonnesting_user_boxes == 'w') {
         TBOX_WARNING("Proceeding with nesting violation as requested.\n"
            << "SAMRAI is not guaranteed to work with nesting"
            << "violations!" << std::endl);
      }

   }
}

/*
 *******************************************************************
 * Check for non-nesting user-specified boxes.
 *******************************************************************
 */
void
GriddingAlgorithm::checkBoundaryProximityViolation(
   const int tag_ln,
   const hier::BoxLevel& new_box_level) const
{
   const size_t nblocks = d_hierarchy->getGridGeometry()->getNumberBlocks();

   const tbox::Dimension& dim = d_hierarchy->getDim();
   hier::IntVector extend_ghosts(dim);
   hier::IntVector smallest_patch(dim);
   hier::IntVector smallest_box_to_refine(dim, 0, nblocks);
   hier::IntVector largest_patch(dim);
   getGriddingParameters(smallest_patch,
      smallest_box_to_refine,
      largest_patch,
      extend_ghosts,
      tag_ln,
      false);

   const size_t nerr(
      checkBoundaryProximityViolation(
         new_box_level,
         extend_ghosts));
   if (nerr > 0 && d_check_boundary_proximity_violation == 'e') {
      TBOX_ERROR("GriddingAlgorithm::checkBoundaryProximityViolation: User error:\n"
         << "Making level " << tag_ln + 1 << ".\n"
         << "New boxes violate boundary proximity.\n"
         << "All boxes must be at least " << extend_ghosts
         << " from physical boundaries or touching the physical boundary."
         << std::endl);
   }
}

/*
 *************************************************************************
 *                                                                       *
 *                                                                       *
 *************************************************************************
 */
void
GriddingAlgorithm::recordStatistics(
   double current_time)
{
#ifdef GA_RECORD_STATS
// GA_RECORD_STATS is defined in GriddingAlgorithm.h
/*
 * For statistics, record number of cells and patches on new level.
 */
   for (int ln = 0; ln < d_hierarchy->getMaxNumberOfLevels(); ++ln) {
      int level_gridcells = 0;
      int level_local_patches = 0;
      if (ln < d_hierarchy->getNumberOfLevels()) {
         const std::shared_ptr<hier::PatchLevel>& patch_level =
            d_hierarchy->getPatchLevel(ln);
         level_gridcells = patch_level->getLocalNumberOfCells();
         level_local_patches = patch_level->getLocalNumberOfPatches();
      }
      d_boxes_stat[ln]->recordProcStat(double(level_local_patches));
      d_cells_stat[ln]->recordProcStat(double(level_gridcells));
      d_timestamp_stat[ln]->recordProcStat(double(current_time));
   }
#endif
}

/*
 *************************************************************************
 *                                                                       *
 *                                                                       *
 *************************************************************************
 */
void
GriddingAlgorithm::printStatistics(
   std::ostream& s) const
{
#ifdef GA_RECORD_STATS
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());
   /*
    * Output statistics.
    */
   // Collect statistic on mesh size.
   tbox::Statistician* statn = tbox::Statistician::getStatistician();

   statn->finalize(false);
   // statn->printLocalStatData(s);
   if (mpi.getRank() == 0) {
      // statn->printAllGlobalStatData(s);
      for (int ln = 0; ln < d_hierarchy->getMaxNumberOfLevels(); ++ln) {
         tbox::Statistic& cstat = *d_cells_stat[ln];
         tbox::Statistic& bstat = *d_boxes_stat[ln];
         tbox::Statistic& tstat = *d_timestamp_stat[ln];
         s << "statistic " << cstat.getName() << ":" << std::endl;
         if (0) {
            s << "Global: \n";
            statn->printGlobalProcStatDataFormatted(cstat.getInstanceId(), s);
         }
         s
         <<
         "Seq#   SimTime           C-Sum   C-Avg   C-Min ->      C-Max  C-Max/Avg     B-Sum    B-Avg B-Min -> B-Max B-Max/Avg  C/B-Avg\n";
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
         for (int sn = 0; sn < cstat.getStatSequenceLength(); ++sn) {
            const double csum = statn->getGlobalProcStatSum(cstat.getInstanceId(), sn);
            const double cmax = statn->getGlobalProcStatMax(cstat.getInstanceId(), sn);
            const double cmin = statn->getGlobalProcStatMin(cstat.getInstanceId(), sn);
            const double cavg = csum / mpi.getSize();
            const double cmaxnorm = cavg != 0 ? cmax / cavg : 1;
            const double bsum = statn->getGlobalProcStatSum(bstat.getInstanceId(), sn);
            const double bmax = statn->getGlobalProcStatMax(bstat.getInstanceId(), sn);
            const double bmin = statn->getGlobalProcStatMin(bstat.getInstanceId(), sn);
            const double bavg = bsum / mpi.getSize();
            const double bmaxnorm = bavg != 0 ? bmax / bavg : 1;
            const double stime = statn->getGlobalProcStatMin(
                  tstat.getInstanceId(), sn);
            s << std::setw(4) << sn
              << " " << std::scientific << std::setprecision(6) << std::setw(12) << stime
              << " " << std::fixed << std::setprecision(0) << std::setw(12) << csum
              << " " << std::setw(7) << cavg
              << " " << std::setw(7) << cmin
              << " -> " << std::setw(10) << cmax
              << "  " << std::setw(9) << std::setprecision(2) << cmaxnorm
              << " " << std::fixed << std::setprecision(0)
              << std::setw(9) << bsum
              << "  " << std::fixed << std::setprecision(2)
              << std::setw(7) << bavg
              << "   " << std::fixed << std::setprecision(0)
              << std::setw(3) << bmin
              << " -> " << std::setw(5) << bmax
              << "  " << std::setw(8) << std::setprecision(2) << bmaxnorm
              << "   " << std::setw(6) << std::setprecision(0) << (bsum != 0 ? csum / bsum : 0)
              << std::endl;
         }
      }
   }
#else
   s << "GriddingAlgorithm statistics is disabled.  See GA_RECORD_STATS in GriddingAlgorithm.h\n";
#endif
}

/*
 *************************************************************************
 * All tags reside in the tag level.  But due to nesting restrictions,
 * not all cells in the level may be refined.
 *
 * We look for the portions of the level that would violate nesting if
 * refined.  Any tags there are nonnesting tags.
 *
 * TODO: remove level from this interface.  tag_ln is sufficient.
 *************************************************************************
 */
void
GriddingAlgorithm::checkNonrefinedTags(
   const hier::PatchLevel& level,
   int tag_ln,
   const hier::OverlapConnectorAlgorithm& oca) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   const tbox::Dimension& dim = d_hierarchy->getDim();
#endif

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dim, level);

   const hier::BoxLevel& tag_box_level = *d_hierarchy->getBoxLevel(tag_ln);
   std::shared_ptr<hier::BoxLevel> violator;
   std::shared_ptr<hier::MappingConnector> tag_to_violator;
   const hier::Connector& tag_box_level_to_self =
      d_hierarchy->getPatchLevel(tag_ln)->findConnector(
         *d_hierarchy->getPatchLevel(tag_ln),
         d_hierarchy->getRequiredConnectorWidth(tag_ln, tag_ln, true),
         hier::CONNECTOR_IMPLICIT_CREATION_RULE,
         false);
   computeNestingViolator(
      violator,
      tag_to_violator,
      tag_box_level,
      tag_box_level_to_self,
      tag_ln,
      oca);

   /*
    * Check for user-tagged cells in the violating parts of the tag level.
    */
   math::PatchCellDataBasicOps<int> dataop;
   int maxval = 0;
   for (hier::Connector::ConstNeighborhoodIterator ei = tag_to_violator->begin();
        ei != tag_to_violator->end(); ++ei) {
      const hier::BoxId& box_id = *ei;
      std::shared_ptr<pdat::CellData<int> > tag_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            level.getPatch(box_id)->getPatchData(d_user_tag_indx)));
      TBOX_ASSERT(tag_data);
      for (hier::Connector::ConstNeighborIterator na = tag_to_violator->begin(ei);
           na != tag_to_violator->end(ei); ++na) {
         const hier::Box& vio_box = *na;
         maxval = dataop.max(tag_data, vio_box);
         if (maxval > 0) {
            break;
         }
      }
      if (maxval > 0) {
         break;
      }
   }
   const tbox::SAMRAI_MPI mpi(tag_box_level.getMPI());
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&maxval, 1, MPI_MAX);
   }

   if (maxval > 0) {
      if (d_check_nonrefined_tags == 'w') {
         TBOX_WARNING("User code has tagged cells in\n"
            << "violation of nesting requirements.\n"
            << "Violating tags will be discarded.\n"
            << "See GriddingAlgorithm::checkNonrefinedTags()\n");
      } else if (d_check_nonrefined_tags == 'e') {
         TBOX_ERROR("User code has tagged cells in\n"
            << "violation of nesting requirements.\n"
            << "See GriddingAlgorithm::checkNonrefinedTags()\n");
      }
   }
}

/*
 *******************************************************************
 * Reset tag buffering data to be able to handle given buffer size.
 *******************************************************************
 */

void GriddingAlgorithm::resetTagBufferingData(const int tag_buffer)
{
   const tbox::Dimension& dim = d_hierarchy->getDim();

   d_buf_tag_ghosts = hier::IntVector(dim, tag_buffer);

   d_bdry_fill_tags.reset();

   hier::VariableDatabase* var_db = hier::VariableDatabase::getDatabase();

   /*
    * Remove d_buf_tag from the VariableDatabase and re-register it with
    * the new ghost width.
    */
   var_db->removeInternalSAMRAIVariablePatchDataIndex(d_buf_tag_indx);

   d_buf_tag_indx = var_db->registerInternalSAMRAIVariable(d_buf_tag,
         d_buf_tag_ghosts);

   if (d_hierarchy->getGridGeometry()->getNumberBlocks() > 1) {
      TBOX_ASSERT(d_mb_tagger_strategy);
      d_mb_tagger_strategy->setScratchTagPatchDataIndex(d_buf_tag_indx);
   }

   d_bdry_fill_tags.reset(new xfer::RefineAlgorithm());

   d_bdry_fill_tags->
   registerRefine(d_buf_tag_indx,
      d_buf_tag_indx,
      d_buf_tag_indx,
      std::shared_ptr<hier::RefineOperator>());
}

/*
 * *************************************************************************
 * *************************************************************************
 */
bool
GriddingAlgorithm::needResetTagBuffer(
   const hier::IntVector& max_ghosts) const
{
   bool ret_val = false;

   int buf_max = d_buf_tag_ghosts.max();
   int ghost_max = max_ghosts.max();
   if (buf_max > ghost_max) {
      ret_val = true;
   } else {
      int ghost_min = max_ghosts.min();

      if (ghost_min != ghost_max) {
         const tbox::Dimension& dim = d_hierarchy->getDim();
         int nblocks = max_ghosts.getNumBlocks();

         for (int b = 0; b < nblocks; ++b) {
            for (int d = 0; d < dim.getValue(); ++d) {
               if (d_buf_tag_ghosts[d] > max_ghosts(b,d)) {
                  ret_val = true;
                  break;
               }
            }
         }
      }
   }
   return ret_val;
}

/*
 *************************************************************************
 *************************************************************************
 */
void
GriddingAlgorithm::checkOverlappingPatches(
   const hier::BoxLevel& box_level) const
{
   hier::Connector box_level_to_self(
      box_level,
      box_level,
      hier::IntVector::getZero(box_level.getDim()));
   hier::OverlapConnectorAlgorithm oca;
   oca.findOverlaps_assumedPartition(box_level_to_self);
   checkOverlappingPatches(box_level_to_self);
}

/*
 *************************************************************************
 *************************************************************************
 */
void
GriddingAlgorithm::checkOverlappingPatches(
   const hier::Connector& box_level_to_self) const
{
   bool has_overlap = false;
   const hier::BoxLevel& box_level = box_level_to_self.getBase();
   const hier::BaseGridGeometry& grid_geom = *box_level.getGridGeometry();
   const hier::IntVector& ratio = box_level.getRefinementRatio();

   for (hier::Connector::ConstNeighborhoodIterator ei = box_level_to_self.begin();
        ei != box_level_to_self.end() && !has_overlap; ++ei) {

      const hier::Box& box = *box_level.getBoxStrict(*ei);

      for (hier::Connector::ConstNeighborIterator na = box_level_to_self.begin(ei);
           na != box_level_to_self.end(ei) && !has_overlap;
           ++na) {
         const hier::Box& nabr = *na;

         if (!nabr.isIdEqual(box)) {
            if (nabr.getBlockId() == box.getBlockId()) {
               has_overlap = nabr.intersects(box);
            } else {
               hier::Box nabr_box(nabr);
               grid_geom.transformBox(nabr_box,
                  ratio,
                  box.getBlockId(),
                  nabr.getBlockId());
               has_overlap = nabr_box.intersects(box);
            }
         }
      }
   }

   if (has_overlap) {
      if (d_check_overlapping_patches == 'w') {
         TBOX_WARNING(
            "PatchLevel has patches which overlap in index space\n"
            << "See GriddingAlgorithm::checkOverlappingPatches().\n"
            <<
            "Note that setting allow_patches_smaller_than_minimum_size_to_prevent_overlaps = FALSE\n"
            <<
            "in the PatchHierarchy can allow some patches to violate min size in order to prevent overlap.\n");
      } else if (d_check_overlapping_patches == 'e') {
         TBOX_ERROR(
            "PatchLevel has patches which overlap in index space\n"
            << "See GriddingAlgorithm::checkOverlappingPatches().\n"
            <<
            "Note that setting allow_patches_smaller_than_minimum_size_to_prevent_overlaps = FALSE\n"
            <<
            "in the PatchHierarchy can allow some patches to violate min size in order to prevent overlap.\n");
      }
   }
}

/*
 *************************************************************************
 *
 * For cases where tagging is not performed read the new level boxes
 * either from user input or from stored level boxes.
 *
 *************************************************************************
 */
void
GriddingAlgorithm::readLevelBoxes(
   std::shared_ptr<hier::BoxLevel>& new_box_level,
   std::shared_ptr<hier::Connector>& coarser_to_new,
   const int tag_ln,
   const double regrid_time,
   const int regrid_cycle,
   bool& remove_old_fine_level)
{
   TBOX_ASSERT((tag_ln >= 0) &&
      (tag_ln <= d_hierarchy->getFinestLevelNumber()));

   const tbox::Dimension& dim = d_hierarchy->getDim();

   const hier::BoxLevel& coarser_box_level = *d_hierarchy->getBoxLevel(tag_ln);

   int fine_level_number = tag_ln + 1;
   hier::BoxContainer boxes_to_refine;

   /*
    * Access the user supplied refine boxes.  The
    * "new_level_has_new_boxes" boolean specifies whether the
    * level boxes have changed from the last time
    * getUserSuppliedRefineBoxes() was called.  If they have changed,
    * it returns true.  If they are unchanged, it returns false.
    */
   bool new_level_has_new_boxes = true;
   if (d_tag_init_strategy->refineUserBoxInputOnly(regrid_cycle, regrid_time)) {

      new_level_has_new_boxes = d_tag_init_strategy->
         getUserSuppliedRefineBoxes(boxes_to_refine,
            tag_ln,
            regrid_cycle,
            regrid_time);

   }

   /*
    * If "new_level_has_new_boxes" is false we wish to keep the
    * existing fine level intact.  Avoid further work by setting
    * the parameter "compute_load_balanced_level_boxes" to false
    * and indicate that we want to avoid removing the old fine level
    * by setting "remove_old_fine_level" to false.
    */
   bool compute_load_balanced_level_boxes = true;
   if (!new_level_has_new_boxes) {
      compute_load_balanced_level_boxes = false;
      remove_old_fine_level = false;
   }

   /*
    * If we are using the nonuniform load balance option, we
    * still need to redo the load balance and construct a new level,
    * even if the level boxes have not changed.
    */

   if (d_load_balancer->getLoadBalanceDependsOnPatchData(fine_level_number)
       && !boxes_to_refine.empty()) {
      compute_load_balanced_level_boxes = true;
      remove_old_fine_level = true;
   }

   /*
    * If the boxes_to_refine are empty, this implies that no
    * refinement is desired so a new finer level will NOT be
    * constructed.  In this case, avoid load balance steps and
    * specify that we want to remove the old fine level.
    */
   if (boxes_to_refine.empty()) {
      compute_load_balanced_level_boxes = false;
      remove_old_fine_level = true;
   }

   if (compute_load_balanced_level_boxes) {

      hier::BoxLevel unbalanced_box_level(
         coarser_box_level.getRefinementRatio(),
         coarser_box_level.getGridGeometry(),
         d_hierarchy->getMPI(),
         hier::BoxLevel::GLOBALIZED);
      hier::LocalId i(0);
      for (hier::BoxContainer::iterator itr = boxes_to_refine.begin();
           itr != boxes_to_refine.end(); ++itr, ++i) {
         hier::Box unbalanced_box(*itr, i, 0);
         unbalanced_box_level.addBox(unbalanced_box);
      }

      const hier::IntVector& ratio =
         d_hierarchy->getRatioToCoarserLevel(fine_level_number);

      new_box_level.reset(new hier::BoxLevel(unbalanced_box_level));
      d_oca0.findOverlapsWithTranspose(coarser_to_new,
         coarser_box_level,
         *new_box_level,
         d_hierarchy->getRequiredConnectorWidth(tag_ln, tag_ln + 1, true),
         hier::IntVector::ceilingDivide(
            d_hierarchy->getRequiredConnectorWidth(tag_ln + 1, tag_ln, true), ratio));

      hier::Connector& new_to_coarser = coarser_to_new->getTranspose();

      size_t num_blocks = new_box_level->getGridGeometry()->getNumberBlocks();
      hier::IntVector smallest_patch(dim);
      hier::IntVector largest_patch(dim);
      hier::IntVector extend_ghosts(dim);
      {
         hier::IntVector smallest_box_to_refine(dim, 0, num_blocks);
         // "false" argument: for_building_finer level = false
         getGriddingParameters(smallest_patch,
            smallest_box_to_refine,
            largest_patch,
            extend_ghosts,
            fine_level_number,
            false);
      }

      refineNewBoxLevel(*new_box_level,
         *coarser_to_new,
         ratio);

      hier::IntVector patch_cut_factor(
         dim,
         d_tag_init_strategy->getErrorCoarsenRatio(),
         num_blocks);
      patch_cut_factor.max(ratio);

      t_load_balance0->start();
      d_load_balancer0->loadBalanceBoxLevel(
         *new_box_level,
         &new_to_coarser,
         d_hierarchy,
         tag_ln,
         smallest_patch,
         largest_patch,
         d_hierarchy->getDomainBoxLevel(),
         extend_ghosts,
         patch_cut_factor);
      t_load_balance0->stop();

      if (d_sequentialize_patch_indices) {
         renumberBoxes(*new_box_level, coarser_to_new.get(), false, true);
      }

      d_blcu0.addPeriodicImages(
         *new_box_level,
         d_hierarchy->getGridGeometry()->getDomainSearchTree(),
         new_to_coarser.getConnectorWidth());
      d_oca0.findOverlaps(*coarser_to_new);
      d_oca0.findOverlaps(new_to_coarser);

   }
}

/*
 *************************************************************************
 * Set all tags on a level to tag_value.
 *************************************************************************
 */

void
GriddingAlgorithm::fillTags(
   const int tag_value,
   const std::shared_ptr<hier::PatchLevel>& tag_level,
   const int tag_index) const
{
   TBOX_ASSERT((tag_value == d_true_tag) || (tag_value == d_false_tag));
   TBOX_ASSERT(tag_level);
   TBOX_ASSERT(tag_index == d_user_tag_indx || tag_index == d_buf_tag_indx
               || tag_index == d_boolean_tag_indx
               || tag_index == d_saved_tag_indx);

   t_fill_tags->start();

   for (hier::PatchLevel::iterator ip(tag_level->begin());
        ip != tag_level->end(); ++ip) {

      const std::shared_ptr<hier::Patch>& patch = *ip;
      std::shared_ptr<pdat::CellData<int> > tag_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch->getPatchData(tag_index)));
      TBOX_ASSERT(tag_data);

      tag_data->fill(tag_value);

   }
#if defined(HAVE_RAJA)
   tbox::parallel_synchronize();
#endif
   t_fill_tags->stop();
}

/*
 *************************************************************************
 *
 * Set each integer value in specified tag array to tag_value where
 * patch level intersects given box array.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::fillTagsFromBoxLevel(
   const int tag_value,
   const std::shared_ptr<hier::PatchLevel>& tag_level,
   const int tag_index,
   const hier::Connector& tag_level_to_fill_box_level,
   const bool interior_only,
   const hier::IntVector& fill_box_growth,
   const bool preserve_existing_tags) const
{
   TBOX_ASSERT((tag_value == d_from_fine_tag) ||
               (tag_value == d_from_fine_pretag));
   TBOX_ASSERT(tag_level);
   TBOX_ASSERT(tag_index == d_user_tag_indx || tag_index == d_buf_tag_indx ||
               tag_index == d_boolean_tag_indx
               || tag_index == d_saved_tag_indx);

   /*
    * This method assumes fill is finer than tag, but that is easy to
    * change, if needed.
    */
   TBOX_ASSERT(tag_level_to_fill_box_level.getHeadCoarserFlag() == false);

   t_fill_tags->start();

   const std::shared_ptr<const hier::BaseGridGeometry>& grid_geom(
      d_hierarchy->getGridGeometry());

   const hier::IntVector& ratio = tag_level_to_fill_box_level.getRatio();

   hier::IntVector growth_in_tag_resolution(fill_box_growth,
                                            grid_geom->getNumberBlocks());
   growth_in_tag_resolution.ceilingDivide(
      tag_level_to_fill_box_level.getRatio());

   for (hier::PatchLevel::iterator ip(tag_level->begin());
        ip != tag_level->end(); ++ip) {
      const std::shared_ptr<hier::Patch>& patch = *ip;

      std::shared_ptr<pdat::CellData<int> > tag_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch->getPatchData(tag_index)));

      TBOX_ASSERT(tag_data);

      const hier::BoxId& box_id(patch->getBox().getBoxId());

      NeighborSet neighbors;

      d_oca.extractNeighbors(
         neighbors,
         tag_level_to_fill_box_level,
         box_id,
         growth_in_tag_resolution);

      for (NeighborSet::const_iterator
           ni = neighbors.begin(); ni != neighbors.end(); ++ni) {
         const hier::Box& neighbor(*ni);
         hier::Box box = neighbor;
         box.grow(fill_box_growth);
         box.coarsen(ratio);
         if (neighbor.getBlockId() != patch->getBox().getBlockId()) {
            grid_geom->transformBox(box,
               tag_level->getLevelNumber(),
               patch->getBox().getBlockId(),
               neighbor.getBlockId());
         }
         if (interior_only) {
            box *= tag_data->getBox();
         }
         if (!preserve_existing_tags) {
            tag_data->fill(tag_value, box);
         } else {
            pdat::CellIterator icend(pdat::CellGeometry::end(box));
            for (pdat::CellIterator ic(pdat::CellGeometry::begin(box));
                 ic != icend; ++ic) {
               if ((*tag_data)(*ic) == d_false_tag) {
                  (*tag_data)(*ic) = tag_value;
               }
            }
         }
      }
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
   }
   t_fill_tags->stop();
}


/*
 *************************************************************************
 * Interprets user tags and sets up algorithmic tags, which will have only
 * values of d_false_tag or d_true_tag, as required by box generator
 * algorithms.
 *************************************************************************
 */
void GriddingAlgorithm::setBooleanTagData(
   const std::shared_ptr<hier::PatchLevel>& tag_level,
   bool preserve_existing_tags) const
{
   for (hier::PatchLevel::iterator ip(tag_level->begin());
        ip != tag_level->end(); ++ip) {
      const std::shared_ptr<hier::Patch>& patch = *ip;

      std::shared_ptr<pdat::CellData<int> > user_tag_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch->getPatchData(d_user_tag_indx)));
      std::shared_ptr<pdat::CellData<int> > boolean_tag_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch->getPatchData(d_boolean_tag_indx)));

      if (!preserve_existing_tags) {
         boolean_tag_data->fillAll(d_false_tag);
      }
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      size_t data_length = boolean_tag_data->getGhostBox().size();
      TBOX_ASSERT(data_length == user_tag_data->getGhostBox().size());

      const int* user_tag_ptr = user_tag_data->getPointer();
      int* boolean_tag_ptr = boolean_tag_data->getPointer();

     for (unsigned int i = 0; i < data_length; ++i) {
        if (user_tag_ptr[i] != d_false_tag) {
           boolean_tag_ptr[i] = d_true_tag;
        }
     }
   }
}

/*
 *************************************************************************
 *
 * Buffer each integer tag with given value on the patch level by the
 * specified buffer size.  Note that the patch data indexed by
 * d_buf_tag_indx is used temporarily to buffer the tag data. The
 * communication of ghost cell (i.e., buffer) information forces all
 * tags on all patch interiors to represent a consistent buffering of
 * the original configuration of tagged cells.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::bufferTagsOnLevel(
   const int tag_value,
   const std::shared_ptr<hier::PatchLevel>& level,
   const int buffer_size) const
{
   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::bufferTagsOnLevel: entered with tag_ln = "
      << level->getLevelNumber() << "\n";
   }

   const tbox::Dimension& dim = d_hierarchy->getDim();

   TBOX_ASSERT((tag_value == d_true_tag) || (tag_value == d_false_tag));
   TBOX_ASSERT(level);
   TBOX_ASSERT(buffer_size >= 0);
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dim, *level);
   /*
    * Start timer for this method.
    */
   t_buffer_tags->start();

   SAMRAI_CALI_MARK_BEGIN("setBuffer");
   /*
    * Set temporary buffered tags based on buffer width and
    * distance from actual tags.
    */
   const int not_tag = ((tag_value == d_true_tag) ? d_false_tag : d_true_tag);
   for (hier::PatchLevel::iterator ip1(level->begin());
        ip1 != level->end(); ++ip1) {
      const std::shared_ptr<hier::Patch>& patch = *ip1;

      std::shared_ptr<pdat::CellData<int> > buf_tag_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch->getPatchData(d_buf_tag_indx)));
      std::shared_ptr<pdat::CellData<int> > boolean_tag_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch->getPatchData(d_boolean_tag_indx)));

      TBOX_ASSERT(buf_tag_data);
      TBOX_ASSERT(boolean_tag_data);

      buf_tag_data->fillAll(not_tag);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      const hier::Box& interior(patch->getBox());

     pdat::CellIterator icend(pdat::CellGeometry::end(interior));
     for (pdat::CellIterator ic(pdat::CellGeometry::begin(interior));
          ic != icend; ++ic) {
        if ((*boolean_tag_data)(*ic) == tag_value) {
           (*buf_tag_data)(*ic) = d_true_tag;
        }
     }
   }
   SAMRAI_CALI_MARK_END("setBuffer");

   SAMRAI_CALI_MARK_BEGIN("bufferBCdata");
   /*
    * Communicate boundary data for buffered tag array so that tags
    * near patch boundaries will become buffered properly.
    */
   const double dummy_time = 0.0;

   t_bdry_fill_tags_comm->start();
   d_bdry_sched_tags[level->getLevelNumber()]->fillData(dummy_time, false);
   t_bdry_fill_tags_comm->stop();

   /*
    * Buffer tags on patch interior according to buffered tag data.
    */
   for (hier::PatchLevel::iterator ip2(level->begin());
        ip2 != level->end(); ++ip2) {
      const std::shared_ptr<hier::Patch>& patch = *ip2;

      std::shared_ptr<pdat::CellData<int> > buf_tag_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch->getPatchData(d_buf_tag_indx)));
      std::shared_ptr<pdat::CellData<int> > boolean_tag_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch->getPatchData(d_boolean_tag_indx)));

      TBOX_ASSERT(buf_tag_data);
      TBOX_ASSERT(boolean_tag_data);

      const hier::Box& tag_box(boolean_tag_data->getBox());
      const hier::BlockId& tag_box_block_id = tag_box.getBlockId();
      hier::Box buf_tag_box(tag_box);
      buf_tag_box.grow(hier::IntVector(dim, buffer_size));

      boolean_tag_data->fillAll(not_tag);

      pdat::CellIterator icend(pdat::CellGeometry::end(buf_tag_box));
      for (pdat::CellIterator ic(pdat::CellGeometry::begin(buf_tag_box));
           ic != icend; ++ic) {
         if ((*buf_tag_data)(*ic) == d_true_tag) {
            hier::Box buf_box(*ic - buffer_size,
                              *ic + buffer_size,
                              tag_box_block_id);
            boolean_tag_data->getArrayData().fillSequential(tag_value, buf_box);
         }
      }
   }
   SAMRAI_CALI_MARK_END("bufferBCdata");

   SAMRAI_CALI_MARK_BEGIN("bufferAdjust");
   /*
    * If a cell has a true boolean tag and a false user tag, the tag is
    * a result of buffering and is set to the d_buffer_tag value in the
    * user tags.
    */
   for (hier::PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
      const std::shared_ptr<hier::Patch>& patch = *ip;

      std::shared_ptr<pdat::CellData<int> > user_tag_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch->getPatchData(d_user_tag_indx)));
      std::shared_ptr<pdat::CellData<int> > boolean_tag_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch->getPatchData(d_boolean_tag_indx)));

      size_t data_length = boolean_tag_data->getGhostBox().size();
      TBOX_ASSERT(data_length == user_tag_data->getGhostBox().size());

      int* user_tag_ptr = user_tag_data->getPointer();
      const int* boolean_tag_ptr = boolean_tag_data->getPointer();

     for (unsigned int i = 0; i < data_length; ++i) {
        if (boolean_tag_ptr[i] == d_true_tag &&
            user_tag_ptr[i] == d_false_tag) {
           user_tag_ptr[i] = d_buffer_tag;
        }
     }
   }
   SAMRAI_CALI_MARK_END("bufferAdjust");
   t_buffer_tags->stop();
}

/*
 *************************************************************************
 *
 * Given a patch level, determine an appropriate array of boxes from
 * which a new finer level may be constructed.  That is, find an array
 * of boxes that covers all tags having the specified tag value.  Note
 * that it is assumed that the integer tag arrays have been set
 * properly; i.e., cells have been tagged through error estimation and
 * the tags have been buffered to ensure disturbances remain on fine
 * level until next regrid occurs.  Note that load balancing is
 * performed once an appropriate list of boxes containing the tags is
 * found.  This procedure massages the list of boxes further and then
 * assigns each to a single processor (i.e., the mapping).
 *
 *************************************************************************
 */

void
GriddingAlgorithm::findRefinementBoxes(
   std::shared_ptr<hier::BoxLevel>& new_box_level,
   std::shared_ptr<hier::Connector>& tag_to_new,
   const int tag_ln) const
{
   TBOX_ASSERT((tag_ln >= 0) &&
      (tag_ln <= d_hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(d_hierarchy->getPatchLevel(tag_ln));

   const tbox::Dimension& dim = d_hierarchy->getDim();

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::findRefinementBoxes: entered with tag_ln = "
      << tag_ln << "\n";
   }

   /*
    * Start timer for this method.
    */
   if (d_barrier_and_time) {
      t_find_refinement->barrierAndStart();
   }

   const hier::BoxLevel& tag_box_level = *d_hierarchy->getBoxLevel(tag_ln);

   const int new_ln = tag_ln + 1;

   /*
    * Construct list of boxes covering the true tags on the level.
    * Note that box list will be contained in the bounding box
    * but will not be contained in the list of proper nesting boxes,
    * in general.  So we intersect the box list against the list of
    * nesting boxes.  Note that this may produce boxes which are too
    * small.  Thus, boxes are regrown later.
    */
   const size_t nblocks = d_hierarchy->getGridGeometry()->getNumberBlocks();

   hier::IntVector smallest_patch(dim);
   hier::IntVector smallest_box_to_refine(dim, 0, nblocks);
   hier::IntVector largest_patch(dim);
   hier::IntVector extend_ghosts(dim);
   // "true" argument: for_building_finer level = true
   getGriddingParameters(smallest_patch,
      smallest_box_to_refine,
      largest_patch,
      extend_ghosts,
      new_ln,
      true);

   const hier::IntVector extend_ghosts_in_tag_space =
      hier::IntVector::ceilingDivide(
         extend_ghosts,
         d_hierarchy->getRatioToCoarserLevel(new_ln));

   std::shared_ptr<hier::PatchLevel> level(
      d_hierarchy->getPatchLevel(tag_ln));

   if (d_print_steps) {
      tbox::plog << "GriddingAlgorithm::findRefinementBoxes: clustering\n";
   }

   t_find_boxes_containing_tags->barrierAndStart();
   hier::IntVector ratio = d_hierarchy->getRatioToCoarserLevel(new_ln);

   hier::BoxContainer bounding_container;
   for (hier::BlockId::block_t bn = 0; bn < nblocks; ++bn) {
      hier::BlockId block_id(bn);
      hier::Box bounding_box(dim);
      bounding_box =
         d_hierarchy->getBoxLevel(tag_ln)->getGlobalBoundingBox(block_id);
      if (!bounding_box.empty()) {
         bounding_container.pushBack(bounding_box);
      }
   }

   if (!bounding_container.empty()) {
      d_box_generator->setMinimumCellRequest(
         d_hierarchy->getMinimumCellRequest(tag_ln+1));
      d_box_generator->setRatioToNewLevel(ratio);
      d_box_generator->findBoxesContainingTags(
         new_box_level,
         tag_to_new,
         level, d_boolean_tag_indx, d_true_tag, bounding_container,
         smallest_box_to_refine,
         d_tag_to_cluster_width[tag_ln]);
   }
   t_find_boxes_containing_tags->stop();

   const hier::PeriodicShiftCatalog& shift_catalog =
      d_hierarchy->getGridGeometry()->getPeriodicShiftCatalog();

   if (new_box_level && new_box_level->getGlobalNumberOfBoxes() > 0) {

      if (d_check_connectors) {
         /*
          * At this stage, there are no edges to periodic images yet, so
          * don't check for them.
          */
         if (d_print_steps) {
            tbox::plog
            <<
            "GriddingAlgorithm::findRefinementBoxes: checking new-->tag from findBoxesContainingTags\n";
         }
         TBOX_ASSERT(tag_to_new->getTranspose().checkOverlapCorrectness(false, true, true) == 0);
         if (d_print_steps) {
            tbox::plog
            <<
            "GriddingAlgorithm::findRefinementBoxes: checking tag-->new from findBoxesContainingTags\n";
         }
         TBOX_ASSERT(tag_to_new->checkOverlapCorrectness(false, true, true) == 0);
      }

      enforceOverflowNesting(*new_box_level, *tag_to_new);

      /*
       * If clustering implementation didn't provide the requested width,
       * recompute the tag<==>new with the right width now.
       *
       * The bridge generates some periodic edges that we don't need just
       * yet, so remove them.
       */
      if (tag_to_new->getConnectorWidth() != d_tag_to_cluster_width[tag_ln]) {
         t_fix_zero_width_clustering->barrierAndStart();
         const hier::Connector& tag_to_tag =
            tag_box_level.findConnectorWithTranspose(
               tag_box_level,
               d_tag_to_cluster_width[tag_ln],
               d_tag_to_cluster_width[tag_ln],
               hier::CONNECTOR_IMPLICIT_CREATION_RULE);
         d_oca.bridgeWithNesting(
            tag_to_new,
            tag_to_tag,
            hier::Connector(*tag_to_new),
            hier::IntVector::getZero(dim),
            hier::IntVector::getZero(dim),
            d_tag_to_cluster_width[tag_ln],
            true);
         if (shift_catalog.isPeriodic()) {
            tag_to_new->removePeriodicRelationships();
            tag_to_new->getTranspose().removePeriodicRelationships();
         }
         t_fix_zero_width_clustering->barrierAndStop();
      }

      if (d_enforce_proper_nesting) {
         enforceProperNesting(*new_box_level, *tag_to_new, tag_ln);
      }

      if (d_extend_to_domain_boundary) {
         extendBoxesToDomainBoundary(
            *new_box_level,
            *tag_to_new,
            level->getPhysicalDomainArray(),
            extend_ghosts_in_tag_space);
      }

      bool allow_patches_smaller_than_minimum_size_to_prevent_overlaps =
         d_hierarchy->allowPatchesSmallerThanMinimumSize();

      // BTNG: these if-else blocks can be significantly simplified by factoring.
      if (!allow_patches_smaller_than_minimum_size_to_prevent_overlaps) {
         if (d_print_steps) {
            tbox::plog
            << "GriddingAlgorithm::findRefinementBoxes: growing boxes\n";
         }
         t_extend_within_domain->start();
         growBoxesWithinNestingDomain(
            *new_box_level,
            *tag_to_new,
            smallest_box_to_refine,
            tag_ln);
         t_extend_within_domain->stop();
      } else {
         const hier::IntVector periodic_dirs(
            d_hierarchy->getGridGeometry()->getPeriodicShift(
               hier::IntVector::getOne(dim)));

         bool need_to_grow = false;
         hier::IntVector min_size(hier::IntVector::getOne(dim));
         for (int i = 0; i < dim.getValue(); ++i) {
            if (periodic_dirs(i)) {
               need_to_grow = true;
               min_size(i) =
                  smallest_box_to_refine(0,i);
            }
         }

         if (need_to_grow) {
            t_extend_within_domain->start();
            growBoxesWithinNestingDomain(
               *new_box_level,
               *tag_to_new,
               min_size,
               tag_ln);
            t_extend_within_domain->stop();
         } else {
            /*
             * Had we need to grow, the growth would have shrunken the widths.
             * We must manually shrink the widths as if we used the growing map
             * (matching the result of an empty map).
             */
            tag_to_new->shrinkWidth(
               tag_to_new->getConnectorWidth() - smallest_box_to_refine);
            tag_to_new->getTranspose().shrinkWidth(
               tag_to_new->getTranspose().getConnectorWidth() -
                  smallest_box_to_refine);
         }
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      std::set<int> new_local_ids;
      const hier::BoxContainer& new_boxes = new_box_level->getBoxes();
      for (hier::BoxContainer::const_iterator new_itr = new_boxes.begin();
           new_itr != new_boxes.end(); ++new_itr) {
         new_local_ids.insert(new_itr->getBoxId().getLocalId().getValue());
      }
      TBOX_ASSERT(static_cast<int>(new_local_ids.size()) == new_boxes.size());
#endif

      /*
       * We have been working with new_box_level in the
       * tag_box_level's index space.  Now, refine it so we can
       * build the new level.
       */
      refineNewBoxLevel(*new_box_level,
         *tag_to_new,
         ratio);

      if (d_check_connectors) {
         TBOX_ASSERT(tag_to_new->getTranspose().checkOverlapCorrectness(false, true, true) == 0);
         TBOX_ASSERT(tag_to_new->checkOverlapCorrectness(false, true, true) == 0);
      }

      if (d_load_balance) {
         if (d_print_steps) {
            tbox::plog
            << "GriddingAlgorithm::findRefinementBoxes: load balancing\n";
         }

         t_load_balance->barrierAndStart();

         const hier::IntVector& patch_cut_factor = ratio;

         d_load_balancer->loadBalanceBoxLevel(
            *new_box_level,
            &tag_to_new->getTranspose(),
            d_hierarchy,
            new_ln,
            smallest_patch,
            largest_patch,
            d_hierarchy->getDomainBoxLevel(),
            extend_ghosts,
            patch_cut_factor);

         t_load_balance->stop();

         if (d_check_connectors) {
            tbox::plog << "GriddingAlgorithm::findRefinementBoxes: checking new-tag" << std::endl;
            int errs = 0;
            if (tag_to_new->getTranspose().checkOverlapCorrectness(false, true, true)) {
               ++errs;
               tbox::perr << "Overlap error found in new_to_tag!\n";
            }
            if (tag_to_new->checkOverlapCorrectness(false, true, true)) {
               ++errs;
               tbox::perr << "Overlap error found in tag_to_new!\n";
            }
            if (tag_to_new->getTranspose().checkTransposeCorrectness(*tag_to_new, true)) {
               ++errs;
               tbox::perr << "Transpose error found in new-tag transpose!\n";
            }
            if (errs != 0) {
               TBOX_ERROR(
                  "Errors found after using load balance map."
                  << "new_box_level:\n" << new_box_level->format("", 2)
                  << "tag_box_level:\n" << tag_box_level.format("", 2)
                  << "new_to_tag:\n" << tag_to_new->getTranspose().format("", 2)
                  << "tag_to_new:\n" << tag_to_new->format("", 2)
                  << std::endl);
            }
         }
      }

      if (d_sequentialize_patch_indices) {
         if (d_print_steps) {
            tbox::plog << "GriddingAlgorithm::findRefinementBoxes: begin sorting boxes."
                       << std::endl;
         }
         renumberBoxes(*new_box_level, tag_to_new.get(), false, true);
         if (d_print_steps) {
            tbox::plog << "GriddingAlgorithm::findRefinementBoxes: end sorting boxes." << std::endl;
         }
      }

      /*
       * Add periodic image Boxes to new_box_level and add edges
       * incident on those nodes.
       */
      const hier::Connector& tag_to_tag =
         d_hierarchy->getPatchLevel(tag_ln)->findConnector(
            *d_hierarchy->getPatchLevel(tag_ln),
            d_hierarchy->getRequiredConnectorWidth(tag_ln, tag_ln, true),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            false);
      if (d_print_steps) {
         tbox::plog << "GriddingAlgorithm::findRefinementBoxes: begin adding periodic images."
                    << std::endl;
      }
      d_blcu.addPeriodicImagesAndRelationships(
         *new_box_level,
         tag_to_new->getTranspose(),
         d_hierarchy->getGridGeometry()->getDomainSearchTree(),
         tag_to_tag);
      if (d_print_steps) {
         tbox::plog << "GriddingAlgorithm::findRefinementBoxes: finished adding periodic images."
                    << std::endl;
      }

      if (d_check_connectors) {
         TBOX_ASSERT(tag_to_new->getTranspose().checkOverlapCorrectness() == 0);
         TBOX_ASSERT(tag_to_new->checkOverlapCorrectness() == 0);
      }

   } else if (new_box_level && new_box_level->getGlobalNumberOfBoxes() == 0) {

      /*
       * On return, new_box_level should be initialized if we
       * generated boxes, deallocated if we didnt.
       */
      new_box_level->clear();

   }

   d_hierarchy->getMPI().Barrier();

   if (d_barrier_and_time) {
      t_find_refinement->stop();
   }

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::findRefinementBoxes: leaving with tag_ln = "
      << tag_ln << "\n";
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
GriddingAlgorithm::renumberBoxes(
   hier::BoxLevel& new_box_level,
   hier::Connector* ref_to_new,
   bool sort_by_corners,
   bool sequentialize_global_indices) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   const tbox::Dimension& dim = d_hierarchy->getDim();
#endif
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dim, new_box_level);

   t_renumber_boxes->barrierAndStart();

   const hier::OverlapConnectorAlgorithm* oca = &d_oca;
   const hier::MappingConnectorAlgorithm* mca = &d_mca;
   const hier::BoxLevelConnectorUtils* blcu = &d_blcu;
   if (ref_to_new && &ref_to_new->getBase() == &d_hierarchy->getDomainBoxLevel()) {
      oca = &d_oca0;
      mca = &d_mca0;
      blcu = &d_blcu0;
   }

   std::shared_ptr<hier::MappingConnector> sorting_map;
   std::shared_ptr<hier::BoxLevel> seq_box_level;
   blcu->makeSortingMap(
      seq_box_level,
      sorting_map,
      new_box_level,
      sort_by_corners,
      sequentialize_global_indices);

   /*
    * Modify works well in most cases, but if ref is actually the
    * domain its cost is O(N^2).  In such a case, recompute instead
    * of modify.
    */
   if (ref_to_new == 0) {
      hier::BoxLevel::swap(new_box_level, *seq_box_level);
   } else if (&ref_to_new->getBase() != &d_hierarchy->getDomainBoxLevel()) {
      mca->modify(*ref_to_new, *sorting_map, &new_box_level);
   } else {
      hier::BoxLevel::swap(new_box_level, *seq_box_level);
      ref_to_new->clearNeighborhoods();
      ref_to_new->getTranspose().clearNeighborhoods();
      ref_to_new->setHead(new_box_level, true);
      ref_to_new->getTranspose().setBase(new_box_level, true);
      oca->findOverlaps_assumedPartition(*ref_to_new);
      ref_to_new->removePeriodicRelationships();
      hier::Connector* new_to_ref = new hier::Connector(ref_to_new->getHead(),
            ref_to_new->getBase(),
            hier::Connector::convertHeadWidthToBase(
               ref_to_new->getHead().getRefinementRatio(),
               ref_to_new->getBase().getRefinementRatio(),
               ref_to_new->getConnectorWidth()));
      oca->findOverlaps(*new_to_ref);
      ref_to_new->setTranspose(new_to_ref, true);
   }

   t_renumber_boxes->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
GriddingAlgorithm::refineNewBoxLevel(
   hier::BoxLevel& new_box_level,
   hier::Connector& tag_to_new,
   const hier::IntVector& ratio) const
{
   TBOX_ASSERT(tag_to_new.hasTranspose());

   hier::Connector& new_to_tag = tag_to_new.getTranspose();

   new_box_level.refineBoxes(new_box_level,
      ratio,
      new_box_level.getRefinementRatio() * ratio);
   new_box_level.finalize();

   hier::IntVector new_to_tag_width =
      ratio * new_to_tag.getConnectorWidth();
   new_to_tag.setBase(new_box_level);
   new_to_tag.setWidth(new_to_tag_width, true);

   tag_to_new.setHead(new_box_level, true);
   tag_to_new.refineLocalNeighbors(ratio);
}

/*
 *************************************************************************
 *************************************************************************
 */

void
GriddingAlgorithm::enforceProperNesting(
   hier::BoxLevel& new_box_level,
   hier::Connector& tag_to_new,
   int tag_ln) const
{
   if (d_barrier_and_time) {
      t_enforce_proper_nesting->barrierAndStart();
   }

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::enforceProperNesting: entered.\n";
   }

   const int new_ln = tag_ln + 1;

   std::shared_ptr<hier::MappingConnector> unnested_to_nested;
   std::shared_ptr<hier::BoxLevel> nested_box_level;

   makeProperNestingMap(
      nested_box_level,
      unnested_to_nested,
      new_box_level,
      tag_to_new.getTranspose(),
      new_ln,
      d_oca);

   if (d_print_steps) {
      tbox::plog << "GriddingAlgorithm::enforceProperNesting: applying proper nesting map.\n";
   }

   t_use_nesting_map->start();
   d_mca.modify(tag_to_new,
      *unnested_to_nested,
      &new_box_level);
   t_use_nesting_map->stop();

   if (tag_ln == d_base_ln && d_check_proper_nesting) {
      /*
       * Tag level will be regridded when we exit the current
       * recursion if tag_ln is not d_base_ln, so do not check
       * proper nesting in that case.
       *
       * Check that the new box_level nest in the tag
       * level (tag_ln).
       */
      hier::IntVector required_nesting(d_hierarchy->getDim());
      if (tag_ln > 0) {
         required_nesting =
            hier::IntVector(d_hierarchy->getDim(), d_hierarchy->getProperNestingBuffer(tag_ln));
      } else {
         required_nesting =
            d_hierarchy->getPatchDescriptor()->getMaxGhostWidth(d_hierarchy->getDim());
      }
      bool locally_nests = false;
      const bool new_nests_in_tag =
         d_blcu.baseNestsInHead(
            &locally_nests,
            new_box_level,
            tag_to_new.getBase(),
            required_nesting,
            hier::IntVector::getZero(d_hierarchy->getDim()),
            hier::IntVector::getZero(d_hierarchy->getDim()),
            &d_hierarchy->getGridGeometry()->getPeriodicDomainSearchTree());
      if (!new_nests_in_tag) {
         tbox::perr << "GriddingAlgorithm" "enforceProperNesting: new BoxLevel\n"
                    << "at ln=" << new_ln
                    << " does not properly nest in\n"
                    << "tag level at tag_ln=" << tag_ln
                    << " by the required nesting buffer of "
                    << required_nesting
                    << ".\nLocal nestingness: " << locally_nests
                    << ".\nWriting BoxLevels out to log file."
                    << std::endl;
         tbox::plog
         << "Proper nesting violation with new BoxLevel of\n"
         << new_box_level.format("N->", 2)
         << "Proper nesting violation with tag BoxLevel of\n"
         << tag_to_new.getBase().format("T->", 2);
         std::shared_ptr<hier::BoxLevel> external;
         std::shared_ptr<hier::Connector> tmp_new_to_tag;
         d_oca.findOverlaps(tmp_new_to_tag,
            new_box_level,
            tag_to_new.getBase(),
            required_nesting);
         tbox::plog << "tmp_new_to_tag:\n" << tmp_new_to_tag->format("NT->", 3);
         std::shared_ptr<hier::MappingConnector> new_to_external;
         d_blcu.computeExternalParts(
            external,
            new_to_external,
            *tmp_new_to_tag,
            -required_nesting,
            d_hierarchy->getGridGeometry()->getDomainSearchTree());
         tbox::plog << "External parts:\n" << new_to_external->format("NE->", 3);
         TBOX_ERROR(
            "Internal library error: Failed to produce proper nesting."
            << std::endl);
      }
   }

   if (d_barrier_and_time) {
      t_enforce_proper_nesting->barrierAndStop();
   }
}

/*
 *************************************************************************
 * Extend Boxes to domain boundary if they are too close.
 *************************************************************************
 */

void
GriddingAlgorithm::extendBoxesToDomainBoundary(
   hier::BoxLevel& new_box_level,
   hier::Connector& tag_to_new,
   const std::vector<hier::BoxContainer>& physical_domain_array,
   const hier::IntVector& extend_ghosts) const
{
   TBOX_ASSERT(tag_to_new.hasTranspose());

   t_extend_to_domain_boundary->barrierAndStart();

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::extendBoxesToDomainBoundary: extending boxes to boundary\n";
   }

   tbox::SAMRAI_MPI mpi(new_box_level.getMPI());

   /*
    * Extend boxes to domain boundary if they are too close to one.  I
    * think there is no need to modify connectivities when a box is
    * extended to the domain boundary.  There is no increased overlap
    * to a finer level, because the finer level already nests in the
    * nodes being extended.  There is increased overlap with the
    * coarser level, but it should not result in additional relationships,
    * because there would not be enough room for an unseen coarse Box
    * to live in the small gap across which the Box is being extended.
    */
   const hier::BoxContainer& before_nodes =
      new_box_level.getBoxes();

   hier::BoxLevel after_box_level(
      new_box_level.getRefinementRatio(),
      new_box_level.getGridGeometry(),
      new_box_level.getMPI());

   hier::MappingConnector before_to_after(
      new_box_level,
      after_box_level,
      extend_ghosts);

   for (hier::BoxContainer::const_iterator nn = before_nodes.begin();
        nn != before_nodes.end(); ++nn) {
      const hier::Box& before_box = *nn;
      hier::Box after_box = before_box;
      hier::BoxUtilities::extendBoxToDomainBoundary(
         after_box,
         physical_domain_array[before_box.getBlockId().getBlockValue()],
         extend_ghosts);
      after_box_level.addBox(after_box);
      before_to_after.insertLocalNeighbor(
         after_box,
         before_box.getBoxId());
   }

   d_mca.modify(tag_to_new,
      before_to_after,
      &new_box_level);

   t_extend_to_domain_boundary->barrierAndStop();
}

/*
 *************************************************************************
 * Enforce overflow nesting by removing new cells that are outside the
 * tag level.
 *************************************************************************
 */

void
GriddingAlgorithm::enforceOverflowNesting(
   hier::BoxLevel& new_box_level,
   hier::Connector& tag_to_new) const
{
   if (d_barrier_and_time) {
      t_enforce_overflow_nesting->barrierAndStart();
   }

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::enforceOverflowNesting: enforcing overflow nesting\n";
   }

   /*
    * Do not allow the new box_level to overflow the tag box_level.
    * If we want to allow the overflow, we have to add the overflow
    * ammount to width of tag->new.  Such additions may make
    * clustering slower.
    */
   std::shared_ptr<hier::BoxLevel> nested_box_level;
   std::shared_ptr<hier::MappingConnector> unnested_to_nested;
   makeOverflowNestingMap(
      nested_box_level,
      unnested_to_nested,
      new_box_level,
      tag_to_new.getTranspose());
   t_use_overflow_map->start();
   d_mca.modify(tag_to_new,
      *unnested_to_nested,
      &new_box_level);
   t_use_overflow_map->stop();

   if (d_check_overflow_nesting) {
      if (d_print_steps) {
         tbox::plog << "GriddingAlgorithm::enforceOverflowNesting: checking overflow."
                    << std::endl;
      }
      bool locally_nested = false;
      bool nested = d_blcu.baseNestsInHead(
            &locally_nested,
            new_box_level,
            tag_to_new.getBase(),
            hier::IntVector::getZero(d_hierarchy->getDim()),
            hier::IntVector::getZero(d_hierarchy->getDim()),
            hier::IntVector::getZero(d_hierarchy->getDim()),
            &d_hierarchy->getGridGeometry()->getDomainSearchTree());
      if (!nested) {
         TBOX_ERROR(
            "Failed overflow nesting: new box_level does not nest in tagged box_level.\n"
            << "Local nestedness = " << locally_nested << std::endl
            << "tag_box_level:\n" << tag_to_new.getBase().format("", 2)
            << "new_box_level:\n" << new_box_level.format("", 2)
            << "tag_to_new:\n" << tag_to_new.format("", 2)
            << "new_to_tag:\n" << tag_to_new.getTranspose().format("", 2)
            << std::endl);
      }
   }

   if (d_barrier_and_time) {
      t_enforce_overflow_nesting->barrierAndStop();
   }
}

/*
 *************************************************************************
 * Make a map that can be used to enforce overflow nesting.
 *************************************************************************
 */

void
GriddingAlgorithm::makeOverflowNestingMap(
   std::shared_ptr<hier::BoxLevel>& nested_box_level,
   std::shared_ptr<hier::MappingConnector>& unnested_to_nested,
   const hier::BoxLevel& unnested_box_level,
   const hier::Connector& unnested_to_reference) const
{
#ifndef DEBUG_CHECK_ASSERTIONS
   NULL_USE(unnested_box_level);
#endif

   const tbox::Dimension& dim = d_hierarchy->getDim();

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dim, unnested_box_level);

   t_make_overflow_map->start();

   std::shared_ptr<hier::BoxLevel> violator_box_level;
   std::shared_ptr<hier::MappingConnector> unnested_to_violator;
   if (d_print_steps) {
      tbox::plog << " GriddingAlgorithm::makeOverflowNestingMap computing external parts."
                 << std::endl;
   }
   t_compute_external_parts->start();
   d_blcu.computeExternalParts(
      violator_box_level,
      unnested_to_violator,
      unnested_to_reference,
      hier::IntVector::getZero(dim),
      d_hierarchy->getGridGeometry()->getDomainSearchTree());
   t_compute_external_parts->stop();

   TBOX_ASSERT(unnested_to_violator->isLocal());

   if (d_print_steps) {
      tbox::plog << " GriddingAlgorithm::makeOverflowNestingMap making remainer map." << std::endl;
   }
   d_blcu.makeRemainderMap(
      nested_box_level,
      unnested_to_nested,
      *unnested_to_violator);
   t_make_overflow_map->stop();

   if (d_print_steps) {
      tbox::plog << " GriddingAlgorithm::makeOverflowNestingMap finished." << std::endl;
   }
}

/*
 *************************************************************************
 * Make a map that can be used to enforce proper nesting.
 * @param unnested_ln Level number for refinenement ratio of unnested_box_level.
 *************************************************************************
 */

void
GriddingAlgorithm::makeProperNestingMap(
   std::shared_ptr<hier::BoxLevel>& nested_box_level,
   std::shared_ptr<hier::MappingConnector>& unnested_to_nested,
   const hier::BoxLevel& unnested_box_level,
   const hier::Connector& unnested_to_hierarchy,
   const int unnested_ln,
   const hier::OverlapConnectorAlgorithm& oca) const
{
   TBOX_ASSERT(unnested_to_hierarchy.hasTranspose());
#ifdef DEBUG_CHECK_ASSERTIONS
   const tbox::Dimension& dim = d_hierarchy->getDim();
#endif

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dim, unnested_box_level);

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::makeProperNesingMap: entered.\n";
   }

   t_make_nesting_map->start();

   std::shared_ptr<hier::MappingConnector> unnested_to_violator;
   std::shared_ptr<hier::BoxLevel> violator;
   computeNestingViolator(
      violator,
      unnested_to_violator,
      unnested_box_level,
      unnested_to_hierarchy,
      unnested_ln - 1,
      oca);

   /*
    * unnested_to_violator is the Connector from the nodes
    * that violate nesting to their violating parts.
    * Convert it to the mapping from unnested to nested.
    */
   d_blcu.makeRemainderMap(
      nested_box_level,
      unnested_to_nested,
      *unnested_to_violator);

   t_make_nesting_map->stop();

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::makeProperNesingMap: exiting.\n";
   }
}

/*
 *************************************************************************
 * Make a map from a BoxLevel to parts of that BoxLevel
 * that violate proper nesting.
 *
 * The violating Boxes are found by comparing candidate
 * Boxes to d_to_nesting_complement's head BoxLevel.
 * Boxes inside the nesting complement violate nesting.
 *************************************************************************
 */
void
GriddingAlgorithm::computeNestingViolator(
   std::shared_ptr<hier::BoxLevel>& violator,
   std::shared_ptr<hier::MappingConnector>& candidate_to_violator,
   const hier::BoxLevel& candidate,
   const hier::Connector& candidate_to_hierarchy,
   const int tag_ln,
   const hier::OverlapConnectorAlgorithm& oca) const
{
   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::computeNestingViolator: entered.\n";
   }

   const tbox::Dimension& dim = d_hierarchy->getDim();

   TBOX_ASSERT(candidate_to_hierarchy.hasTranspose());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dim, candidate);

   hier::IntVector zero_vector(hier::IntVector::getZero(dim));

   // Check requirements on arguments.
   TBOX_ASSERT(candidate_to_hierarchy.getRatio() == 1);
   TBOX_ASSERT(candidate_to_hierarchy.getTranspose().getRatio() == 1);

   const hier::BaseGridGeometry& grid_geometry(
      *d_hierarchy->getGridGeometry());

   t_compute_nesting_violator->start();

   const hier::BoxContainer& candidate_boxes = candidate.getBoxes();
   /*
    * Bridge candidate to d_nesting_complement.  Any part of the
    * candidate BoxLevel that overlaps d_nesting_complement
    * violates nesting.
    */
   std::shared_ptr<hier::Connector> candidate_to_complement;

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::computeNestingViolator: bridging for candidate_to_complement.\n";
   }
   oca.bridge(candidate_to_complement,
      candidate_to_hierarchy,
      *d_to_nesting_complement[tag_ln],
      false);
   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::computeNestingViolator: bridged for candidate_to_complement.\n";
   }

   d_blcu.computeInternalParts(
      violator,
      candidate_to_violator,
      *candidate_to_complement,
      zero_vector,
      grid_geometry.getDomainSearchTree());
   /*
    * Above step ignored the domain complement components of nesting
    * definition (by necessity).  Where the candidate falls outside
    * the domain, it violates nesting.
    */

   hier::BoxContainer refined_domain_search_tree(
      d_hierarchy->getGridGeometry()->getDomainSearchTree());
   refined_domain_search_tree.refine(candidate.getRefinementRatio());
   refined_domain_search_tree.makeTree(&grid_geometry);

   for (hier::BoxContainer::const_iterator ni = candidate_boxes.begin();
        ni != candidate_boxes.end(); ++ni) {
      const hier::Box& cmb = *ni;
      hier::BoxContainer addl_violators(cmb);
      addl_violators.removeIntersections(
         candidate.getRefinementRatio(),
         refined_domain_search_tree);
      if (!addl_violators.empty()) {
         /*
          * Non-periodic BoxId needed for NeighborhoodSet::find()
          */
         hier::BoxId cmb_non_per_id(cmb.getGlobalId(),
                                    hier::PeriodicId::zero());
         if (candidate_to_violator->hasNeighborSet(cmb_non_per_id)) {
            /*
             * Remove parts that we already know, through
             * candidate_to_violator, are non-nesting.  Leftovers are
             * non-nesting parts not found using
             * candidate_to_complement.
             */
            hier::Connector::NeighborhoodIterator base_box_itr =
               candidate_to_violator->makeEmptyLocalNeighborhood(cmb_non_per_id);
            hier::Connector::ConstNeighborhoodIterator current_violators =
               candidate_to_violator->find(cmb_non_per_id);
            for (hier::Connector::ConstNeighborIterator na =
                    candidate_to_violator->begin(current_violators);
                 na != candidate_to_violator->end(current_violators) && !addl_violators.empty();
                 ++na) {
               addl_violators.removeIntersections(*na);
            }
            if (!addl_violators.empty()) {
               for (hier::BoxContainer::iterator bi = addl_violators.begin();
                    bi != addl_violators.end(); ++bi) {
                  hier::BoxContainer::const_iterator new_violator = violator->addBox(
                        *bi, cmb.getBlockId());
                  candidate_to_violator->insertLocalNeighbor(*new_violator,
                     base_box_itr);
               }
            }
         }
      }
   }

   t_compute_nesting_violator->stop();

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::computeNestingViolator: exiting.\n";
   }
}

/*
 *************************************************************************
 * Precompute data used to define proper nesting.  Data is associated
 * with level number ln, to be used for constructing level number ln+1.
 *
 * Data computed: d_proper_nesting_complement[ln],
 * d_to_proper_nesting_complement[ln] and its transpose.
 *
 * If ln > d_base_ln, assume data at ln-1 is already set.
 *************************************************************************
 */

void
GriddingAlgorithm::computeProperNestingData(
   const int ln,
   const hier::OverlapConnectorAlgorithm& oca)
{
   t_compute_proper_nesting_data->start();

   TBOX_ASSERT(d_base_ln >= 0 && ln >= d_base_ln);

   const tbox::Dimension& dim = d_hierarchy->getDim();

   if (ln == d_base_ln) {
      /*
       * At the base level, nesting domain is level d_base_ln,
       * shrunken by d_proper_nesting_buffer[d_base_ln].
       */
      const hier::Connector& self_connector =
         d_hierarchy->getPatchLevel(ln)->findConnector(
            *d_hierarchy->getPatchLevel(ln),
            d_hierarchy->getRequiredConnectorWidth(ln, ln, true),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            false);

      // This assert shoud pass due to GriddingAlgorithmConnectorWidthRequestor.
      TBOX_ASSERT(self_connector.getConnectorWidth() >=
         hier::IntVector(dim, -d_hierarchy->getProperNestingBuffer(ln)));

      std::shared_ptr<hier::MappingConnector> to_nesting_complement;
      d_blcu.computeExternalParts(
         d_proper_nesting_complement[ln],
         to_nesting_complement,
         self_connector,
         hier::IntVector(dim, -d_hierarchy->getProperNestingBuffer(ln)),
         d_hierarchy->getGridGeometry()->getDomainSearchTree());

      /*
       * Change to_nesting_complement from a mapping to an overlap Connector
       * by adding trivial edges normally omitted from mapping Connectors.
       */
      const hier::BoxContainer& tag_level_boxes = self_connector.getBase().getBoxes();
      for (hier::BoxContainer::const_iterator bi = tag_level_boxes.begin();
           bi != tag_level_boxes.end();
           ++bi) {
         if (!to_nesting_complement->hasNeighborSet(bi->getBoxId())) {
            to_nesting_complement->insertLocalNeighbor(*bi, bi->getBoxId());
         }
      }

      d_to_nesting_complement[ln] =
         std::static_pointer_cast<hier::Connector, hier::MappingConnector>(
            to_nesting_complement);
      d_to_nesting_complement[ln]->setTranspose(
         d_to_nesting_complement[ln]->createLocalTranspose(),
         true);

   } else {

      TBOX_ASSERT(d_to_nesting_complement[ln - 1]->isFinalized());

      /*
       * How to build d_proper_nesting_complement[ln] and connect it to level ln:
       *
       * In the left column are the BoxLevels in the hierarchy.
       * In the right are their proper nesting complements for that level number.
       *
       *                   (new
       *                 Connector)
       *                     |
       *             Mapped  |   Proper
       *               box   |   nesting
       *             levels  | complements
       *             ======  | ===========
       *                     v
       *               ln <-----> ln
       *                ^        ^
       *                |       /
       *                |      /
       *                |     /
       * (existing      |    / <--(temporary
       *  Connector)--> |   /      Connector)
       *                |  /
       *                | /
       *                |/
       *                v
       *             ln-1 <-----> ln-1
       *                     ^
       *                     |
       *                 (existing
       *                  Connector)
       *
       * We have existing Connectors between ln and ln-1 and also between
       * level ln-1 and the nesting complement at ln-1.
       *
       * 1. Build the complement at ln (d_proper_nesting_complement[ln])
       *    from the complement at ln-1 (d_proper_nesting_complement[ln-1])
       *    by refining and growing the complement boxes.
       *
       * 2. Build the temporary Connector from level ln-1 to
       *    complements at ln by using the fact that the complement at ln
       *    is similar to the one at ln-1.
       *
       * 3. Bridge for the new Connector from level ln to the
       *    complement at ln, using the temporary Connector.
       */

      /*
       * 1. Build the complement at ln (d_proper_nesting_complement[ln])
       *    from the complement at ln-1 (d_proper_nesting_complement[ln-1]).
       */
      d_proper_nesting_complement[ln].reset(new hier::BoxLevel(
            d_hierarchy->getBoxLevel(ln)->getRefinementRatio(),
            d_hierarchy->getGridGeometry(),
            d_to_nesting_complement[ln - 1]->getMPI()));
      const hier::BoxContainer& lnm1_complement_boxes =
         d_proper_nesting_complement[ln - 1]->getBoxes();
      for (hier::BoxContainer::const_iterator ni =
              lnm1_complement_boxes.begin();
           ni != lnm1_complement_boxes.end(); ++ni) {
         hier::Box tmp_box = *ni;
         TBOX_ASSERT(!tmp_box.isPeriodicImage());
         tmp_box.refine(d_hierarchy->getRatioToCoarserLevel(ln));
         tmp_box.grow(
            hier::IntVector(dim, d_hierarchy->getProperNestingBuffer(ln)));
         d_proper_nesting_complement[ln]->addBox(tmp_box);
      }

      /*
       * 2. Temporarily connect level ln-1 and d_proper_nesting_complement[ln].
       */
      hier::Connector lnm1_to_ln_complement(*d_hierarchy->getBoxLevel(ln - 1),
                                            *d_proper_nesting_complement[ln],
                                            d_to_nesting_complement[ln - 1]->getConnectorWidth());
      for (hier::Connector::ConstNeighborhoodIterator ei =
              d_to_nesting_complement[ln - 1]->begin();
           ei != d_to_nesting_complement[ln - 1]->end(); ++ei) {
         for (hier::Connector::ConstNeighborIterator na =
                 d_to_nesting_complement[ln - 1]->begin(ei);
              na != d_to_nesting_complement[ln - 1]->end(ei); ++na) {
            hier::Box tmp_box = *na;
            tmp_box.refine(d_hierarchy->getRatioToCoarserLevel(ln));
            tmp_box.grow(
               hier::IntVector(dim, d_hierarchy->getProperNestingBuffer(ln)));
            lnm1_to_ln_complement.insertLocalNeighbor(tmp_box, *ei);
         }
      }
      hier::Connector& from_nesting_complement =
         d_to_nesting_complement[ln - 1]->getTranspose();
      hier::Connector ln_complement_to_lnm1(from_nesting_complement);
      ln_complement_to_lnm1.setBase(*d_proper_nesting_complement[ln]);
      ln_complement_to_lnm1.setHead(*d_hierarchy->getBoxLevel(ln - 1));
      ln_complement_to_lnm1.setWidth(
         from_nesting_complement.getConnectorWidth()
         * d_hierarchy->getRatioToCoarserLevel(ln),
         true);
      lnm1_to_ln_complement.setTranspose(&ln_complement_to_lnm1, false);

      /*
       * 3. Bridge for Connector between level ln and d_proper_nesting_complement[ln].
       */
      oca.bridge(
         d_to_nesting_complement[ln],
         d_hierarchy->getPatchLevel(ln)->findConnectorWithTranspose(
            *d_hierarchy->getPatchLevel(ln - 1),
            d_hierarchy->getRequiredConnectorWidth(ln, ln - 1, true),
            d_hierarchy->getRequiredConnectorWidth(ln - 1, ln),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            false),
         lnm1_to_ln_complement,
         d_hierarchy->getRequiredConnectorWidth(ln - 1, ln, true),
         true);
   }

   t_compute_proper_nesting_data->stop();
}

/*
 *************************************************************************
 * Make a mapping Connector that can be used to grow boxes within
 * nesting domain by the minimum amount needed to make all boxes in a
 * BoxLevel satisfy the min_size requirement.
 *
 * Apply the map.
 *************************************************************************
 */

void
GriddingAlgorithm::growBoxesWithinNestingDomain(
   hier::BoxLevel& new_box_level,
   hier::Connector& tag_to_new,
   const hier::IntVector& min_size,
   const int tag_ln) const
{
   TBOX_ASSERT(tag_to_new.hasTranspose());

   const tbox::Dimension& dim = d_hierarchy->getDim();

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY2(dim,
      new_box_level,
      min_size);

   hier::Connector& new_to_tag = tag_to_new.getTranspose();

   const hier::BaseGridGeometry& grid_geometry(
      *d_hierarchy->getGridGeometry());
   const size_t nblocks = grid_geometry.getNumberBlocks();
   hier::IntVector current_min_size(dim, tbox::MathUtilities<int>::getMax());
   for (hier::BlockId::block_t bn = 0; bn < nblocks; ++bn) {
      current_min_size.min(new_box_level.getGlobalMinBoxSize(hier::BlockId(bn)));
   }

   if (min_size <= current_min_size) {
      /*
       * No box growing is needed.  Just shrink the Connector widths
       * to mimic expected the side-effect of applying a map with a
       * width of min_size.  The code should give the same result
       * without this special bypass.
       */
      tag_to_new.shrinkWidth(tag_to_new.getConnectorWidth() - min_size);
      new_to_tag.shrinkWidth(new_to_tag.getConnectorWidth() - min_size);
      return;
   }

   const hier::BoxContainer& new_boxes = new_box_level.getBoxes();

   const hier::Connector& tag_to_nesting_complement =
      *d_to_nesting_complement[tag_ln];

   /*
    * Connect new_box_level to the nesting complement so we
    * know where it cannot exist.  Use a Connector width of zero
    * because we don't need it any bigger and we don't want
    * inter-block neighbors.  (The new Boxes are already confined to
    * their own blocks before entering this method, so a zero
    * Connector width eliminates inter-block neighbors.)
    */

   std::shared_ptr<hier::Connector> new_to_nesting_complement;
   d_oca.bridge(
      new_to_nesting_complement,
      new_to_tag,
      tag_to_nesting_complement,
      min_size - 1,
      false);

   /*
    * Set up the empty grown_box_level to be populated as we
    * determine whether each box needs to be grown.
    */
   hier::BoxLevel grown_box_level(
      new_box_level.getRefinementRatio(),
      new_box_level.getGridGeometry(),
      new_box_level.getMPI());

   // Create the mapping Connector from new to grown.
   hier::MappingConnector new_to_grown(
      new_box_level,
      grown_box_level,
      min_size);

   hier::BoxContainer refined_domain_search_tree(
      grid_geometry.getDomainSearchTree());
   refined_domain_search_tree.refine(new_box_level.getRefinementRatio());
   refined_domain_search_tree.makeTree(&grid_geometry);

   std::vector<hier::Box> tmp_box_vector;
   tmp_box_vector.reserve(10);

   /*
    * Loop through the new Boxes and grow if needed.  For each
    * Box, determine a sufficient view of the domain where it
    * can grow into.  This view includes domain boxes overlapping the
    * Box minus parts removed to satisfy nesting requirements.
    */

   for (hier::BoxContainer::const_iterator ni = new_boxes.begin();
        ni != new_boxes.end(); ++ni) {
      const hier::Box& omb = *ni;
      TBOX_ASSERT(!omb.isPeriodicImage());

      if (omb.numberCells() >= min_size.getBlockVector(omb.getBlockId())) {
         // This box does not need growing.
         grown_box_level.addBox(omb);
         continue;
      }

      if (new_to_nesting_complement->hasNeighborSet(omb.getBoxId())) {
         // Box omb is near the nesting boundary and may touch it.
         hier::BoxContainer nearby_nesting_boundary;
         new_to_nesting_complement->getNeighborBoxes(omb.getBoxId(), nearby_nesting_boundary);
         hier::Box grown_box = omb;
         hier::BoxUtilities::growBoxWithinDomain(
            grown_box,
            nearby_nesting_boundary,
            min_size.getBlockVector(omb.getBlockId()));
         grown_box_level.addBox(grown_box);
         new_to_grown.insertLocalNeighbor(grown_box, omb.getBoxId());
      } else {
         // Box omb is not near the nesting boundary and need not be grown.
         grown_box_level.addBox(omb);
      }

   }

   /*
    * Use the mapping Connector.
    */

   d_mca.modify(tag_to_new,
      new_to_grown,
      &new_box_level);
}

void
GriddingAlgorithm::getGriddingParameters(
   hier::IntVector& smallest_patch,
   hier::IntVector& smallest_box_to_refine,
   hier::IntVector& largest_patch,
   hier::IntVector& extend_ghosts,
   const int level_number,
   const bool for_building_finer) const
{
   const tbox::Dimension& dim = d_hierarchy->getDim();

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY4(dim,
      smallest_patch,
      smallest_box_to_refine,
      largest_patch,
      extend_ghosts);

   TBOX_ASSERT((level_number >= 0) && (level_number < d_hierarchy->getMaxNumberOfLevels()));

   /*
    * Determine maximum ghost cell width needed over all variables
    * currently known to the patch descriptor, and set the smallest
    * patch size.  The maximum number of ghosts is multiplied by the
    * error coarsen ratio (which should always be 1 unless regridding
    * uses error estimation).  This assures that when levels are
    * coarsened during error estimation, the coarser level patches
    * will meet the ghost cell constraint.
    */
   bool allow_patches_smaller_than_ghost_width =
      d_hierarchy->allowPatchesSmallerThanGhostWidth();

   hier::IntVector max_ghosts(
      d_hierarchy->getPatchDescriptor()->getMaxGhostWidth(dim));
   max_ghosts = max_ghosts * d_tag_init_strategy->getErrorCoarsenRatio();
   smallest_patch = d_hierarchy->getSmallestPatchSize(level_number);
   if (!allow_patches_smaller_than_ghost_width) {
      smallest_patch.max(max_ghosts);
   } else {
      const hier::IntVector periodic_dirs(
         d_hierarchy->getGridGeometry()->getPeriodicShift(
            hier::IntVector::getOne(dim)));

      for (int i = 0; i < dim.getValue(); ++i) {
         if (periodic_dirs(i)) {
            smallest_patch(i) =
               tbox::MathUtilities<int>::Max(smallest_patch(i), max_ghosts(i));
         }
      }
   }

   /*
    * Set largest patch size.
    */
   largest_patch = d_hierarchy->getLargestPatchSize(level_number);

   /*
    * Following if-check prevents changing a negative largest_patch
    * because a negative value dissables the upper limit on patch size.
    */
   if (largest_patch > hier::IntVector::getZero(dim)) {
      largest_patch.max(smallest_patch);
   }

   /*
    * Set the smallest box to refine based on the number of cells that
    * coarsened patches must accomodate to meet ghost cell needs of variables.
    * On the finest level, the smallest box to refine is the smallest patch.
    * On coarser levels, it is a function of the error coarsen ratio and
    * the ratio to the next finer level.
    *
    * If we are accessing gridding parameters for a level that is being
    * reconstructed, the smallest box to refine is not applicable so we
    * set it to -1 to indicate an invalid entry in case it is used.
    */
   if (for_building_finer) {

      smallest_box_to_refine.setAll(smallest_patch);

      /*
       * Shouldn't this division be rounded up because it
       * represents a coarsening?  BTNG.
       */
      smallest_box_to_refine /=
         d_hierarchy->getRatioToCoarserLevel(level_number);
      /*
       * den = ratio from level_number to Richardson-coarsened
       * version of level_number+1.
       */
      const hier::IntVector den(
         d_hierarchy->getRatioToCoarserLevel(level_number)
         / d_tag_init_strategy->getErrorCoarsenRatio());
      /*
       * sz = max ghosts on Richardson-coarsened level_number+1, as
       * seen on level_number.
       */
      hier::IntVector multi_max_ghosts(max_ghosts,
                                       d_hierarchy->getNumberBlocks());

      const hier::IntVector sz(hier::IntVector::ceilingDivide(
                                  multi_max_ghosts, den));
      smallest_box_to_refine.max(sz);

   } else {

      smallest_box_to_refine.setAll(hier::IntVector(dim, -1));

   }

   /*
    * Determine number of cells box may be extended to physical
    * domain boundary to accomodate ghost cells.
    */
   if (extend_ghosts.getNumBlocks() == 1 &&
       d_hierarchy->getNumberBlocks() > 1) {
      extend_ghosts = hier::IntVector(max_ghosts,
                                      d_hierarchy->getNumberBlocks());
   } else {
      extend_ghosts.setAll(max_ghosts);
   }
}

/*
 *************************************************************************
 *************************************************************************
 */

void
GriddingAlgorithm::warnIfDomainTooSmallInPeriodicDir() const
{
   const tbox::Dimension& dim = d_hierarchy->getDim();

   const hier::PeriodicShiftCatalog& shift_catalog =
      d_hierarchy->getGridGeometry()->getPeriodicShiftCatalog();

   if (shift_catalog.isPeriodic()) {

      hier::IntVector periodic_shift(
         d_hierarchy->getGridGeometry()->getPeriodicShift(
            hier::IntVector::getOne(dim)));

      hier::IntVector domain_bounding_box_size(
         d_hierarchy->getDomainBoxLevel().
         getGlobalBoundingBox(hier::BlockId(0)).numberCells());

      for (int ln = 0; ln < d_hierarchy->getNumberOfLevels(); ++ln) {

         if (ln > 0) {
            periodic_shift *= d_hierarchy->getRatioToCoarserLevel(ln);
            domain_bounding_box_size *= d_hierarchy->getRatioToCoarserLevel(ln);
         }

         hier::IntVector smallest_patch_size(dim);
         hier::IntVector largest_patch_size(dim);
         hier::IntVector extend_ghosts(dim);
         hier::IntVector smallest_box_to_refine(dim, 0);
         // "false" argument: for_building_finer level = false
         getGriddingParameters(
            smallest_patch_size,
            smallest_box_to_refine,
            largest_patch_size,
            extend_ghosts,
            ln,
            false);

         for (int d = 0; d < dim.getValue(); ++d) {
            if (periodic_shift(d) > 0 &&
                domain_bounding_box_size(d) < smallest_patch_size(d)) {
               TBOX_WARNING(
                  "GriddingAlgorithm::warnIfDomainTooSmallInPeriodicDir: domain bounding box size\n"
                  << domain_bounding_box_size << " is smaller\n"
                  << "than the smallest patch size "
                  << smallest_patch_size << " on level "
                  << ln << " in direction " << d << "\n");
               break;
            }
         }

      }

   }
}

/*
 *************************************************************************
 *
 * Print out all attributes of class instance for debugging.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::printClassData(
   std::ostream& os) const
{
   os << "\nGriddingAlgorithm::printClassData..." << std::endl;
   os << "GriddingAlgorithm: this = "
      << (GriddingAlgorithm *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_tag_init_strategy = "
      << d_tag_init_strategy.get() << std::endl;
   os << "d_box_generator = "
      << d_box_generator.get() << std::endl;
   os << "d_load_balancer = "
      << d_load_balancer.get() << std::endl;
   os << "d_load_balancer0 = "
      << d_load_balancer0.get() << std::endl;
   os << "d_user_tag = " << d_user_tag.get() << std::endl;
   os << "d_user_tag_indx = " << d_user_tag_indx << std::endl;
   os << "d_saved_tag = " << d_user_tag.get() << std::endl;
   os << "d_saved_tag_indx = " << d_user_tag_indx << std::endl;
   os << "d_boolean_tag = " << d_boolean_tag.get() << std::endl;
   os << "d_boolean_tag_indx = " << d_boolean_tag_indx << std::endl;
   os << "d_buf_tag_indx = " << d_buf_tag_indx << std::endl;
   os << "d_true_tag = " << d_true_tag << std::endl;
   os << "d_false_tag = " << d_false_tag << std::endl;
}

/*
 *************************************************************************
 *
 * Write out class version number and data members to restart database.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("ALGS_GRIDDING_ALGORITHM_VERSION",
      ALGS_GRIDDING_ALGORITHM_VERSION);

   restart_db->putBool("check_overflow_nesting", d_check_overflow_nesting);
   restart_db->putBool("check_proper_nesting", d_check_proper_nesting);
   restart_db->putBool("DEV_check_connectors", d_check_connectors);
   restart_db->putBool("DEV_print_steps", d_print_steps);
   restart_db->putBool("DEV_log_metadata_statistics", d_log_metadata_statistics);

   restart_db->putChar("check_nonrefined_tags", d_check_nonrefined_tags);
   restart_db->putChar("check_overlapping_patches",
      d_check_overlapping_patches);
   restart_db->putChar("check_nonnesting_user_boxes",
      d_check_nonnesting_user_boxes);
   restart_db->putChar("DEV_check_boundary_proximity_violation",
      d_check_boundary_proximity_violation);

   restart_db->putBool("sequentialize_patch_indices",
      d_sequentialize_patch_indices);

   restart_db->putBool("enforce_proper_nesting", d_enforce_proper_nesting);
   restart_db->putBool("DEV_extend_to_domain_boundary",
      d_extend_to_domain_boundary);
   restart_db->putBool("DEV_load_balance", d_load_balance);

   restart_db->putBool("DEV_barrier_and_time", d_barrier_and_time);

   restart_db->putBool("save_tag_data", d_save_tag_data);
}

/*
 *************************************************************************
 *
 * If simulation is not from restart, read data from input database.
 * Otherwise, override data members initialized from restart with
 * values in the input database.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::getFromInput(
   const std::shared_ptr<tbox::Database>& input_db,
   bool is_from_restart)
{
   if (input_db) {
      if (!is_from_restart) {

         d_check_overflow_nesting =
            input_db->getBoolWithDefault("check_overflow_nesting", false);
         d_check_proper_nesting =
            input_db->getBoolWithDefault("check_proper_nesting", false);
         d_check_connectors =
            input_db->getBoolWithDefault("DEV_check_connectors", false);
         d_print_steps =
            input_db->getBoolWithDefault("DEV_print_steps", false);
         d_log_metadata_statistics =
            input_db->getBoolWithDefault("DEV_log_metadata_statistics", false);

         std::string tmp_str;

         tmp_str =
            input_db->getStringWithDefault("check_nonrefined_tags", std::string("WARN"));
         if (!(tmp_str == "IGNORE" || tmp_str == "WARN" || tmp_str == "ERROR")) {
            INPUT_VALUE_ERROR("check_nonrefined_tags");
         }
         d_check_nonrefined_tags = char(tolower(*tmp_str.c_str()));

         tmp_str =
            input_db->getStringWithDefault("check_overlapping_patches", std::string("IGNORE"));
         if (!(tmp_str == "IGNORE" || tmp_str == "WARN" || tmp_str == "ERROR")) {
            INPUT_VALUE_ERROR("check_overlapping_patches");
         }
         d_check_overlapping_patches = char(tolower(*tmp_str.c_str()));

         tmp_str =
            input_db->getStringWithDefault("check_nonnesting_user_boxes", std::string("ERROR"));
         if (!(tmp_str == "IGNORE" || tmp_str == "WARN" || tmp_str == "ERROR")) {
            INPUT_VALUE_ERROR("check_nonnesting_user_boxes");
         }
         d_check_nonnesting_user_boxes = char(tolower(*tmp_str.c_str()));

         tmp_str =
            input_db->getStringWithDefault("DEV_check_boundary_proximity_violation",
               std::string("ERROR"));
         if (!(tmp_str == "IGNORE" || tmp_str == "WARN" || tmp_str == "ERROR")) {
            INPUT_VALUE_ERROR("DEV_check_boundary_proximity_violation");
         }
         d_check_boundary_proximity_violation = char(tolower(*tmp_str.c_str()));

         d_sequentialize_patch_indices =
            input_db->getBoolWithDefault("sequentialize_patch_indices", true);

         d_enforce_proper_nesting =
            input_db->getBoolWithDefault("enforce_proper_nesting", true);
         d_extend_to_domain_boundary =
            input_db->getBoolWithDefault("DEV_extend_to_domain_boundary", true);
         d_load_balance =
            input_db->getBoolWithDefault("DEV_load_balance", true);

         d_barrier_and_time =
            input_db->getBoolWithDefault("DEV_barrier_and_time", false);

         d_save_tag_data =
            input_db->getBoolWithDefault("save_tag_data", false);

      } else {
         bool read_on_restart =
            input_db->getBoolWithDefault("read_on_restart", false);
         if (!read_on_restart) {
            return;
         }

         d_check_overflow_nesting =
            input_db->getBoolWithDefault("check_overflow_nesting",
               d_check_overflow_nesting);
         d_check_proper_nesting =
            input_db->getBoolWithDefault("check_proper_nesting",
               d_check_proper_nesting);
         d_check_connectors =
            input_db->getBoolWithDefault("DEV_check_connectors",
               d_check_connectors);
         d_print_steps =
            input_db->getBoolWithDefault("DEV_print_steps", d_print_steps);
         d_log_metadata_statistics =
            input_db->getBoolWithDefault("DEV_log_metadata_statistics",
               d_log_metadata_statistics);

         std::string tmp_str;

         if (input_db->keyExists("check_nonrefined_tags")) {
            tmp_str = input_db->getString("check_nonrefined_tags");
            d_check_nonrefined_tags = char(tolower(*tmp_str.c_str()));
            if (d_check_nonrefined_tags != 'i' &&
                d_check_nonrefined_tags != 'w' &&
                d_check_nonrefined_tags != 'e') {
               TBOX_ERROR(
                  "GriddingAlgorithm::getFromInput: input parameter check_nonrefined_tags\n"
                  << "can only be \"IGNORE\", \"WARN\" or \"ERROR\""
                  << std::endl);
            }
         }

         if (input_db->keyExists("check_overlapping_patches")) {
            tmp_str = input_db->getString("check_overlapping_patches");
            d_check_overlapping_patches = char(tolower(*tmp_str.c_str()));
            if (d_check_overlapping_patches != 'i' &&
                d_check_overlapping_patches != 'w' &&
                d_check_overlapping_patches != 'e') {
               TBOX_ERROR(
                  "GriddingAlgorithm::getFromInput: input parameter check_overlapping_patches\n"
                  << "can only be \"IGNORE\", \"WARN\" or \"ERROR\""
                  << std::endl);
            }
         }

         if (input_db->keyExists("check_nonnesting_user_boxes")) {
            tmp_str = input_db->getString("check_nonnesting_user_boxes");
            d_check_nonnesting_user_boxes = char(tolower(*tmp_str.c_str()));
            if (d_check_nonnesting_user_boxes != 'i' &&
                d_check_nonnesting_user_boxes != 'w' &&
                d_check_nonnesting_user_boxes != 'e') {
               TBOX_ERROR(
                  "GriddingAlgorithm::getFromInput: input parameter check_nonnesting_user_boxes\n"
                  << "can only be \"IGNORE\", \"WARN\" or \"ERROR\""
                  << std::endl);
            }
         }

         if (input_db->keyExists("DEV_check_boundary_proximity_violation")) {
            tmp_str =
               input_db->getString("DEV_check_boundary_proximity_violation");
            d_check_boundary_proximity_violation =
               char(tolower(*tmp_str.c_str()));
            if (d_check_boundary_proximity_violation != 'i' &&
                d_check_boundary_proximity_violation != 'w' &&
                d_check_boundary_proximity_violation != 'e') {
               TBOX_ERROR(
                  "GriddingAlgorithm::getFromInput: input parameter check_boundary_proximity_violation\n"
                  << "can only be \"IGNORE\", \"WARN\" or \"ERROR\""
                  << std::endl);
            }
         }

         d_sequentialize_patch_indices =
            input_db->getBoolWithDefault("sequentialize_patch_indices",
               d_sequentialize_patch_indices);

         d_enforce_proper_nesting =
            input_db->getBoolWithDefault("enforce_proper_nesting",
               d_enforce_proper_nesting);
         d_extend_to_domain_boundary =
            input_db->getBoolWithDefault("DEV_extend_to_domain_boundary",
               d_extend_to_domain_boundary);
         d_load_balance =
            input_db->getBoolWithDefault("DEV_load_balance", d_load_balance);

         d_barrier_and_time =
            input_db->getBoolWithDefault("DEV_barrier_and_time",
               d_barrier_and_time);

         d_save_tag_data =
            input_db->getBoolWithDefault("save_tag_data", false);
      }
   }
}

/*
 *************************************************************************
 *
 * Gets the database in the root database that corresponds to the object
 * name.  This method then checks to make sure that the version number
 * of the class is that same as the version number in the restart file.
 * If these values are equal, the data members are read in from the
 * restart database.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::getFromRestart()
{
   std::shared_ptr<tbox::Database> root_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in restart file." << std::endl);
   }
   std::shared_ptr<tbox::Database> db(root_db->getDatabase(d_object_name));

   int ver = db->getInteger("ALGS_GRIDDING_ALGORITHM_VERSION");
   if (ver != ALGS_GRIDDING_ALGORITHM_VERSION) {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "Restart file version different than class version."
                       << std::endl);
   }

   d_check_overflow_nesting = db->getBool("check_overflow_nesting");
   d_check_proper_nesting = db->getBool("check_proper_nesting");
   d_check_connectors = db->getBool("DEV_check_connectors");
   d_print_steps = db->getBool("DEV_print_steps");
   d_log_metadata_statistics = db->getBool("DEV_log_metadata_statistics");

   d_check_nonrefined_tags = db->getChar("check_nonrefined_tags");
   d_check_overlapping_patches = db->getChar("check_overlapping_patches");
   d_check_nonnesting_user_boxes = db->getChar("check_nonnesting_user_boxes");
   d_check_boundary_proximity_violation =
      db->getChar("DEV_check_boundary_proximity_violation");

   d_sequentialize_patch_indices = db->getBool("sequentialize_patch_indices");

   d_enforce_proper_nesting = db->getBool("enforce_proper_nesting");
   d_extend_to_domain_boundary = db->getBool("DEV_extend_to_domain_boundary");
   d_load_balance = db->getBool("DEV_load_balance");

   d_barrier_and_time = db->getBool("DEV_barrier_and_time");

   d_save_tag_data = db->getBool("save_tag_data");
}

/*
 *************************************************************************
 *************************************************************************
 */
void
GriddingAlgorithm::allocateTimers()
{
   /*
    * Timers:  for gathering performance information about box
    * calculus and other regridding operations.
    */
   t_load_balance = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::load_balance");
   t_load_balance0 = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::load_balance0");
   t_bdry_fill_tags_create = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::bdry_fill_tags_create");
   t_make_coarsest = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeCoarsestLevel()");
   t_make_finer = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeFinerLevel()");
   t_regrid_all_finer = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::regridAllFinerLevels()");
   t_regrid_finer_do_tagging_before = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::regridFinerLevel_doTaggingBeforeRecursiveRegrid()");
   t_regrid_finer_do_tagging_after = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::regridFinerLevel_doTaggingAfterRecursiveRegrid()");
   t_regrid_finer_create_and_install = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::regridFinerLevel_createAndInstallNewLevel()");
   t_regrid_finer_create = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::regridFinerLevel()_create");
   t_fill_tags = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::fillTagsFromBoxLevel()");
   t_tag_cells_for_refinement = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::tag_cells_for_refinement");
   t_initialize_level_data = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::initialize_level_data");
   t_buffer_tags = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::bufferTagsOnLevel()");
   t_second_finer_tagging = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::second_finer_tagging");
   t_bdry_fill_tags_comm = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::bdry_fill_tags_comm");
   t_find_refinement = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::findRefinementBoxes()");
   t_find_boxes_containing_tags = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::find_boxes_containing_tags");
   t_fix_zero_width_clustering = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::fix_zero_width_clustering");
   t_compute_proper_nesting_data = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::computeProperNestingData()");
   t_enforce_proper_nesting = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::enforceProperNesting()");
   t_make_nesting_map = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeProperNestingMap()");
   t_use_nesting_map = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::use_nesting_map");
   t_make_overflow_map = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeOverflowNestingMap()");
   t_use_overflow_map = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::use_overflow_map");
   t_compute_external_parts = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::compute_external_parts");
   t_compute_nesting_violator = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::computeNestingViolator()");
   t_extend_to_domain_boundary = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::extendBoxesToDomainBoundary()");
   t_extend_within_domain = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::extend_within_domain");
   t_grow_boxes_within_domain = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::grow_boxes_within_domain");
   t_renumber_boxes = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::renumberBoxes()");
   t_find_new_to_new = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::find_new_to_new");
   t_bridge_new_to_new = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::bridge_new_to_new");
   t_bridge_new_to_finer = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::bridge_new_to_finer");
   t_bridge_new_to_old = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::bridge_new_to_old");
   t_make_domain = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeCoarsestLevel()_make_domain");
   t_make_new = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeCoarsestLevel()_make_new");
   t_process_error = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::process_error");
   t_reset_hier = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::reset_hierarchy_config");
   t_enforce_overflow_nesting = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::enforceOverflowNesting()");
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
