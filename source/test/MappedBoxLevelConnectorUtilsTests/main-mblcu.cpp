/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test program for performance of tree search algorithm.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/BoxContainerSingleBlockIterator.h"
#include "SAMRAI/geom/GridGeometry.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/HDFDatabase.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

#include <algorithm>
#include <vector>
#include <iomanip>

using namespace SAMRAI;
using namespace tbox;

/*
 * Partition boxes in the given BoxLevel.  This method is meant
 * to partition and also to create a bunch of small boxes from the
 * user-input boxes in order to set up a non-trivial mesh
 * configuration.
 */
void
partitionBoxes(
   hier::BoxLevel& box_level,
   hier::BoxLevel& domain_box_level,
   const hier::IntVector& max_box_size,
   const hier::IntVector& min_box_size);

/*
 * Shrink a BoxLevel by the given amount.
 */
void
shrinkBoxLevel(
   std::shared_ptr<hier::BoxLevel>& small_box_level,
   const hier::BoxLevel& big_box_level,
   const hier::IntVector& shrinkage,
   const std::vector<hier::BlockId::block_t>& unshrunken_blocks);

/*
 * Refine a BoxLevel by the given ratio.
 */
void
refineBoxLevel(
   hier::BoxLevel& box_level,
   const hier::IntVector& ratio);

/*
 ************************************************************************
 *
 * This is an correctness test for the
 * computeExternalParts and
 * computeInternalParts method in
 * BoxLevelConnectorUtils.
 *
 * 1. Set up GridGeometry.
 *
 * 2. Build a big BoxLevel.
 *
 * 3. Build a small BoxLevel by shrinking the big one slightly.
 *
 * 4. Partition the two BoxLevels.
 *
 * 5. Check the internal and external parts the big and small
 *    BoxLevels with respect to each other.
 *
 * Inputs:
 *
 *   GridGeometry { ... }
 *
 *   Main {
 *      // Domain dimension
 *      dim = 2
 *
 *      // Base name of output files generated.
 *      base_name = "..."
 *
 *      // Logging option
 *      log_all_nodes  = TRUE
 *
 *      // For breaking up user-specified domain boxes.
 *      max_box_size = 7, 7
 *      min_box_size = 2, 2
 *
 *      // Index space to exclude from the big BoxLevel.
 *      exclude2 = [(0,0,4), (19,19,7)]
 *
 *      // Blocks not to be shrunken when generating small BoxLevel
 *      // This allows some small blocks to touch the domain boundary.
 *      unshrunken_boxes = 1, 3
 *
 *      // Refinement ratio of big BoxLevel.
 *      big_refinement_ratio = 2, 2
 *
 *      // Refinement ratio of small BoxLevel.
 *      small_refinement_ratio = 6, 6
 *   }
 */

int main(
   int argc,
   char* argv[])
{
   /*
    * Initialize MPI, SAMRAI.
    */

   SAMRAI_MPI::init(&argc, &argv);
   SAMRAIManager::initialize();
   SAMRAIManager::startup();
   tbox::SAMRAI_MPI mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   int fail_count = 0;

   {

      /*
       * Process command line arguments.  For each run, the input
       * filename must be specified.  Usage is:
       *
       * executable <input file name>
       */
      std::string input_filename;

      if (argc != 2) {
         TBOX_ERROR("USAGE:  " << argv[0] << " <input file> \n"
                               << "  options:\n"
                               << "  none at this time" << std::endl);
      } else {
         input_filename = argv[1];
      }

      /*
       * Create input database and parse all data in input file.
       */

      std::shared_ptr<InputDatabase> input_db(new InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Retrieve "Main" section from input database.
       * The main database is used only in main().
       * The base_name variable is a base name for
       * all name strings in this program.
       */

      std::shared_ptr<Database> main_db(input_db->getDatabase("Main"));

      const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      std::string base_name = "unnamed";
      base_name = main_db->getStringWithDefault("base_name", base_name);

      /*
       * Start logging.
       */
      const std::string log_filename = base_name + ".log";
      bool log_all_nodes = false;
      log_all_nodes = main_db->getBoolWithDefault("log_all_nodes",
            log_all_nodes);
      if (log_all_nodes) {
         PIO::logAllNodes(log_filename);
      } else {
         PIO::logOnlyNodeZero(log_filename);
      }

      /*
       * Generate the grid geometry.
       */
      if (!input_db->keyExists("GridGeometry")) {
         TBOX_ERROR("BoxLevelConnectorUtils test: could not find entry GridGeometry"
            << "\nin input.");
      }
      std::shared_ptr<hier::BaseGridGeometry> grid_geometry(
         new geom::GridGeometry(
            dim,
            "GridGeometry",
            input_db->getDatabase("GridGeometry")));
      grid_geometry->printClassData(tbox::plog);

      /*
       * Empty blocks are blocks not to have any Boxes.
       */
      std::vector<int> empty_blocks;
      if (main_db->isInteger("empty_blocks")) {
         empty_blocks = main_db->getIntegerVector("empty_blocks");
      }

      /*
       * Unshrunken blocks are blocks in which the small BoxLevel
       * has the same index space as the big BoxLevel.
       */
      std::vector<hier::BlockId::block_t> unshrunken_blocks;
      if (main_db->isInteger("unshrunken_blocks")) {
         std::vector<int> input_unshrunken =
            main_db->getIntegerVector("unshrunken_blocks");
         for (std::vector<int>::const_iterator itr = input_unshrunken.begin();
              itr != input_unshrunken.end(); ++itr) {
            unshrunken_blocks.push_back(
               static_cast<hier::BlockId::block_t>(*itr));
         }
      }

      plog << "Input database after initialization..." << std::endl;
      input_db->printClassData(plog);

      hier::IntVector one_vector(hier::IntVector::getOne(dim));
      hier::IntVector zero_vector(hier::IntVector::getZero(dim));

      /*
       * How much to shrink the big BoxLevel to get the small one.
       */
      const hier::IntVector shrinkage(hier::IntVector::getOne(dim));

      hier::BoxLevelConnectorUtils mblcu;

      /*
       * Set up the domain.
       */

      hier::BoxContainer domain_boxes;
      grid_geometry->computePhysicalDomain(domain_boxes, one_vector);
      tbox::plog << "domainboxes:\n"
                 << domain_boxes.format()
                 << std::endl;

      domain_boxes.makeTree(grid_geometry.get());

      /*
       * Construct the big_box_level.  It is a refinement of
       * the domain but without exclude* boxes.
       */

      hier::BoxLevel big_box_level(
         one_vector,
         grid_geometry,
         tbox::SAMRAI_MPI::getSAMRAIWorld());
      const std::string exclude("exclude");
      for (hier::BlockId::block_t bn = 0; bn < grid_geometry->getNumberBlocks(); ++bn) {

         const hier::BlockId block_id(bn);

         const std::string exclude_boxes_name = exclude + tbox::Utilities::intToString(static_cast<int>(bn));
         if (main_db->keyExists(exclude_boxes_name)) {

            /*
             * Block bn boxes in big_box_level are the
             * block_domain \ exclude_boxes.
             */

            hier::BoxContainer block_domain;
            grid_geometry->computePhysicalDomain(block_domain,
               one_vector,
               block_id);

            std::vector<tbox::DatabaseBox> db_box_vector =
               main_db->getDatabaseBoxVector(exclude_boxes_name);
            hier::BoxContainer exclude_boxes(db_box_vector);
            for (hier::BoxContainer::iterator itr = exclude_boxes.begin();
                 itr != exclude_boxes.end(); ++itr) {
               itr->setBlockId(block_id);
            }
            block_domain.unorder();
            block_domain.removeIntersections(exclude_boxes);

            hier::LocalId last_local_id(-1);
            for (hier::BoxContainer::iterator bi = block_domain.begin();
                 bi != block_domain.end(); ++bi) {
               big_box_level.addBoxWithoutUpdate(
                  hier::Box(*bi,
                     ++last_local_id,
                     0));
            }

         } else {

            /*
             * Block bn boxes in big_box_level are the same as
             * block bn domain.
             */
            for (hier::BoxContainerSingleBlockIterator bi(domain_boxes.begin(block_id));
                 bi != domain_boxes.end(block_id); ++bi) {
               big_box_level.addBoxWithoutUpdate(*bi);
            }

         }

      }
      big_box_level.finalize();

      const hier::BoxContainer& big_boxes(big_box_level.getBoxes());

      hier::BoxLevel big_domain_level = big_box_level;

      /*
       * Generate the "small" BoxLevel by shrinking the big one
       * back at its boundary.
       */
      std::shared_ptr<hier::BoxLevel> small_box_level;
      shrinkBoxLevel(small_box_level,
         big_box_level,
         shrinkage,
         unshrunken_blocks);

      hier::BoxLevel small_domain_level = *small_box_level;

      std::vector<hier::IntVector> refinement_ratios(
         1, hier::IntVector(dim, 1, grid_geometry->getNumberBlocks()));

      /*
       * Refine Boxlevels as user specified.
       */
      if (main_db->isInteger("big_refinement_ratio")) {
         hier::IntVector big_refinement_ratio(dim);
         main_db->getIntegerArray("big_refinement_ratio",
            &big_refinement_ratio[0],
            dim.getValue());
         refineBoxLevel(big_box_level, big_refinement_ratio);
      }
      refinement_ratios.push_back(big_box_level.getRefinementRatio());

      if (main_db->isInteger("small_refinement_ratio")) {
         hier::IntVector small_refinement_ratio(dim);
         main_db->getIntegerArray("small_refinement_ratio",
            &small_refinement_ratio[0],
            dim.getValue());
         refineBoxLevel(*small_box_level, small_refinement_ratio);
      }
      TBOX_ASSERT(small_box_level->getRefinementRatio() >
                  big_box_level.getRefinementRatio());
      refinement_ratios.push_back(small_box_level->getRefinementRatio() /
                                  big_box_level.getRefinementRatio());

      /*
       * These steps are usually handled by PatchHierarchy, but this
       * test does not use PatchHierarchy.
       */
      grid_geometry->setUpRatios(refinement_ratios);

      /*
       * Partition the big and small BoxLevels.
       *
       * Limit box sizes to make the configuration more complex than
       * the domain description.  Default is not to limit box size.
       */
      hier::IntVector max_box_size(dim, tbox::MathUtilities<int>::getMax());
      if (main_db->isInteger("max_box_size")) {
         main_db->getIntegerArray("max_box_size",
            &max_box_size[0],
            dim.getValue());
      }
      hier::IntVector min_box_size(dim, 2);
      if (main_db->isInteger("min_box_size")) {
         main_db->getIntegerArray("min_box_size",
            &min_box_size[0],
            dim.getValue());
      }
      partitionBoxes(*small_box_level, small_domain_level,
         max_box_size, min_box_size);
      partitionBoxes(big_box_level, big_domain_level,
         max_box_size, min_box_size);

      big_box_level.cacheGlobalReducedData();
      small_box_level->cacheGlobalReducedData();

      tbox::plog << "\nbig_box_level:\n"
                 << big_box_level.format("", 2)
                 << '\n'
                 << "small_box_level:\n"
                 << small_box_level->format("", 2)
                 << '\n'
      ;

      const hier::BoxContainer& small_boxes(small_box_level->getBoxes());

      const hier::BoxContainer small_box_tree(
         small_box_level->getGlobalizedVersion().getGlobalBoxes());
      small_box_tree.makeTree(grid_geometry.get());

      /*
       * Connectors between big and small BoxLevels.
       */

      const hier::Connector& small_to_big(
         small_box_level->createConnectorWithTranspose(big_box_level,
            refinement_ratios.back() * shrinkage,
            shrinkage));
      small_to_big.cacheGlobalReducedData();

      const hier::Connector& big_to_small = small_to_big.getTranspose();
      big_to_small.cacheGlobalReducedData();

      tbox::plog << "\nsmall_to_big:\n"
                 << small_to_big.format("", 2)
                 << '\n'
                 << "big_to_small:\n"
                 << big_to_small.format("", 2)
                 << '\n';

      /*
       * Setup is complete.  Begin testing.
       */

      {
         /*
          * small_box_level nests inside big_box_level
          * by the shrinkage amount.  Verify that
          * computeInternalParts finds all
          * small_box_level to be internal to
          * big_box_level and that
          * computeExternalParts finds none of
          * small_box_level to be external to
          * big_box_level.
          *
          * small_box_level's internal parts should include
          * everything.  Thus small_to_everything should only map
          * small_box_level to its own index space (no more and
          * no less).
          *
          * small_box_level's external parts should include
          * nothing.  Thus small_to_nothing should map
          * small_box_level to nothing.
          */
         std::shared_ptr<hier::BoxLevel> everything;
         std::shared_ptr<hier::BoxLevel> nothing;
         std::shared_ptr<hier::MappingConnector> small_to_everything,
                                                   small_to_nothing;
         mblcu.computeExternalParts(
            nothing,
            small_to_nothing,
            small_to_big,
            -shrinkage,
            domain_boxes);
         mblcu.computeInternalParts(
            everything,
            small_to_everything,
            small_to_big,
            -shrinkage,
            domain_boxes);
         tbox::plog << "\nsmall_to_nothing:\n"
                    << small_to_nothing->format("", 2) << '\n'
                    << "\nnothing:\n"
                    << nothing->format("", 2) << '\n'
                    << "small_to_everything:\n"
                    << small_to_everything->format("", 2) << '\n'
                    << "\neverything:\n"
                    << everything->format("", 2) << '\n'
         ;

         for (hier::BoxContainer::const_iterator bi = small_boxes.begin();
              bi != small_boxes.end(); ++bi) {
            const hier::Box& small_box = *bi;

            if (small_to_everything->hasNeighborSet(small_box.getBoxId())) {
               hier::Connector::ConstNeighborhoodIterator neighbors =
                  small_to_everything->find(small_box.getBoxId());

               hier::BoxContainer neighbor_box_list;
               for (hier::Connector::ConstNeighborIterator na =
                       small_to_everything->begin(neighbors);
                    na != small_to_everything->end(neighbors); ++na) {
                  if (!(*na).empty()) {
                     neighbor_box_list.pushBack(*na);

                     if (!small_box.contains(*na)) {
                        tbox::perr << "Mapping small_to_everyting erroneously mapped "
                                   << small_box << " to:\n" << *na
                                   << " which is outside itself.\n";
                        ++fail_count;
                     }
                  }
               }

               hier::BoxContainer tmp_box_list(small_box);
               tmp_box_list.removeIntersections(neighbor_box_list);
               if (tmp_box_list.size() != 0) {
                  tbox::perr << "Mapping small_to_everything erroneously mapped "
                             << small_box << " to something less than itself:\n";
                  small_to_everything->writeNeighborhoodToStream(
                     tbox::perr,
                     small_box.getBoxId());
               }

            }

            if (small_to_nothing->hasNeighborSet(small_box.getBoxId())) {
               if (!small_to_nothing->isEmptyNeighborhood(
                      small_box.getBoxId())) {
                  tbox::perr << "Mapping small_to_nothing erroneously mapped " << small_box
                             << " to:\n";
                  small_to_nothing->writeNeighborhoodToStream(
                     tbox::perr,
                     small_box.getBoxId());
                  tbox::perr << "\nIt should be mapped to nothing\n";
                  ++fail_count;
               }
            } else {
               tbox::perr << "Mapping small_to_nothing is missing a map from "
                          << small_box << " to nothing.\n";
               ++fail_count;
            }

         }
      }

      {
         /*
          * Compute the parts of big_box_level that are
          * internal to small_box_level and check for correctness.
          *
          * To verify that the internal parts are correctly computed:
          *
          * - check that the small_box_level and the internal
          * parts of big_box_level have the same index space.
          */

         std::shared_ptr<hier::BoxLevel> internal_box_level;
         std::shared_ptr<hier::MappingConnector> big_to_internal;
         mblcu.computeInternalParts(
            internal_box_level,
            big_to_internal,
            big_to_small,
            zero_vector);
         const hier::BoxContainer& internal_boxes(internal_box_level->getBoxes());
         tbox::plog << "internal_box_level:\n"
                    << internal_box_level->format("", 2)
                    << '\n'
                    << "big_to_internal:\n"
                    << big_to_internal->format("", 2);

         hier::BoxContainer internal_box_tree(
            internal_box_level->getGlobalizedVersion().getGlobalBoxes());
         internal_box_tree.makeTree(grid_geometry.get());

         for (hier::BoxContainer::const_iterator ni = small_boxes.begin();
              ni != small_boxes.end(); ++ni) {
            hier::BoxContainer tmp_box_list(*ni);
            small_to_big.getHeadCoarserFlag() ?
            tmp_box_list.coarsen(small_to_big.getRatio()) :
            tmp_box_list.refine(small_to_big.getRatio());
            tmp_box_list.removeIntersections(
               big_box_level.getRefinementRatio(),
               internal_box_tree);
            if (tmp_box_list.size() > 0) {
               tbox::perr << "Small box " << *ni << " should fall within "
                          << "the internal index space, but it doesn't." << std::endl;
               ++fail_count;
            }
         }

         for (hier::BoxContainer::const_iterator ni = internal_boxes.begin();
              ni != internal_boxes.end(); ++ni) {
            hier::BoxContainer tmp_box_list(*ni);
            big_to_small.getHeadCoarserFlag() ?
            tmp_box_list.coarsen(big_to_small.getRatio()) :
            tmp_box_list.refine(big_to_small.getRatio());
            tmp_box_list.removeIntersections(
               small_box_level->getRefinementRatio(),
               small_box_tree);
            if (tmp_box_list.size() > 0) {
               tbox::perr << "Internal box " << *ni << " should fall within "
                          << "the small index space, but it doesn't." << std::endl;
               ++fail_count;
            }
         }

      }

      {

         /*
          * Compute parts of big_box_level that are external to
          * small_box_level.
          *
          * Verify that the external parts are correctly computed:
          *
          * - check that external parts do not overlap small_box_level.
          *
          * - check that small_box_level does not overlap external parts.
          *
          * - check that big_box_level \ { small_box_level, external parts }
          *   is empty.
          */
         std::shared_ptr<hier::BoxLevel> external_box_level;
         std::shared_ptr<hier::MappingConnector> big_to_external;
         mblcu.computeExternalParts(
            external_box_level,
            big_to_external,
            big_to_small,
            zero_vector,
            hier::BoxContainer());
         const hier::BoxContainer& external_boxes(external_box_level->getBoxes());
         tbox::plog << "\nexternal_box_level:\n"
                    << external_box_level->format("", 2)
                    << '\n'
                    << "big_to_external:\n"
                    << big_to_external->format("", 2);

         hier::BoxContainer external_box_tree(
            external_box_level->getGlobalizedVersion().getGlobalBoxes());
         external_box_tree.makeTree(grid_geometry.get());

         for (hier::BoxContainer::const_iterator ni = external_boxes.begin();
              ni != external_boxes.end(); ++ni) {
            hier::BoxContainer tmp_box_list(*ni);
            big_to_small.getHeadCoarserFlag() ?
            tmp_box_list.coarsen(big_to_small.getRatio()) :
            tmp_box_list.refine(big_to_small.getRatio());
            tmp_box_list.intersectBoxes(
               small_box_level->getRefinementRatio(),
               small_box_tree);
            if (tmp_box_list.size() != 0) {
               tbox::perr << "External box " << *ni << " should not\n"
                          << "intersect small_box_level but does.\n"
                          << "Intersections:\n";
               tmp_box_list.print(tbox::perr);
               tbox::perr << std::endl;
               ++fail_count;
            }
         }

         for (hier::BoxContainer::const_iterator ni = small_boxes.begin();
              ni != small_boxes.end(); ++ni) {
            hier::BoxContainer tmp_box_list(*ni);
            small_to_big.getHeadCoarserFlag() ?
            tmp_box_list.coarsen(small_to_big.getRatio()) :
            tmp_box_list.refine(small_to_big.getRatio());
            tmp_box_list.intersectBoxes(
               big_box_level.getRefinementRatio(),
               external_box_tree);
            if (tmp_box_list.size() != 0) {
               tbox::perr << "Small box " << *ni << " should not intersect "
                          << "the external parts but is does.\n"
                          << "Intersections:\n";
               tmp_box_list.print(tbox::perr);
               tbox::perr << std::endl;
               ++fail_count;
            }
         }

         for (hier::BoxContainer::const_iterator ni = big_boxes.begin();
              ni != big_boxes.end(); ++ni) {
            hier::BoxContainer tmp_box_list(*ni);
            big_to_small.getHeadCoarserFlag() ?
            tmp_box_list.coarsen(big_to_small.getRatio()) :
            tmp_box_list.refine(big_to_small.getRatio());
            tmp_box_list.removeIntersections(
               small_box_level->getRefinementRatio(),
               small_box_tree);
            small_to_big.getHeadCoarserFlag() ?
            tmp_box_list.coarsen(small_to_big.getRatio()) :
            tmp_box_list.refine(small_to_big.getRatio());
            tmp_box_list.removeIntersections(
               big_box_level.getRefinementRatio(),
               external_box_tree);
            if (tmp_box_list.size() > 0) {
               tbox::perr << "Big box " << *ni << " should not be inside "
                          << "the small BoxLevel and the external parts but is not.\n"
                          << "Outside parts:\n";
               tmp_box_list.print(tbox::perr);
               tbox::perr << std::endl;
               ++fail_count;
            }
         }

      }

      if (fail_count == 0) {
         tbox::pout << "\nPASSED:  BoxLevelConnector test" << std::endl;
      }

      input_db.reset();
      main_db.reset();

      /*
       * Exit properly by shutting down services in correct order.
       */
      tbox::plog << "\nShutting down..." << std::endl;

   }

   /*
    * Shut down.
    */
   SAMRAIManager::shutdown();
   SAMRAIManager::finalize();
   SAMRAI_MPI::finalize();

   return fail_count;
}

/*
 * Partition boxes in the given BoxLevel.  This method is meant
 * to partition and also to create a bunch of small boxes from the
 * user-input boxes in order to set up a non-trivial mesh
 * configuration.
 */
void partitionBoxes(
   hier::BoxLevel& box_level,
   hier::BoxLevel& domain_box_level,
   const hier::IntVector& max_box_size,
   const hier::IntVector& min_box_size) {

   const tbox::Dimension& dim(box_level.getDim());

   domain_box_level.setParallelState(hier::BoxLevel::GLOBALIZED);

   mesh::TreeLoadBalancer load_balancer(box_level.getDim(),
                                        "TreeLoaadBalaancer");

   hier::Connector* dummy_connector = 0;

   const hier::IntVector bad_interval(dim, 1);
   const hier::IntVector cut_factor(hier::IntVector::getOne(dim));

   load_balancer.loadBalanceBoxLevel(
      box_level,
      dummy_connector,
      std::shared_ptr<hier::PatchHierarchy>(),
      0,
      min_box_size,
      max_box_size,
      domain_box_level,
      bad_interval,
      cut_factor);
}

void shrinkBoxLevel(
   std::shared_ptr<hier::BoxLevel>& small_box_level,
   const hier::BoxLevel& big_box_level,
   const hier::IntVector& shrinkage,
   const std::vector<hier::BlockId::block_t>& unshrunken_blocks)
{
   const std::shared_ptr<const hier::BaseGridGeometry>& grid_geometry(
      big_box_level.getGridGeometry());

   const int local_rank = big_box_level.getMPI().getRank();

   const hier::BoxContainer& big_boxes(big_box_level.getBoxes());

   const hier::Connector& big_to_big(
      big_box_level.createConnector(big_box_level, shrinkage));

   hier::BoxContainer visible_boxes(big_boxes);
   for (hier::Connector::ConstNeighborhoodIterator mi = big_to_big.begin();
        mi != big_to_big.end(); ++mi) {
      for (hier::Connector::ConstNeighborIterator ni = big_to_big.begin(mi);
           ni != big_to_big.end(mi); ++ni) {
         visible_boxes.insert(*ni);
      }
   }

   hier::BoxContainer boundary_boxes = visible_boxes;

   hier::BoxLevelConnectorUtils mblcu;

   mblcu.computeBoxesAroundBoundary(
      boundary_boxes,
      big_box_level.getRefinementRatio(),
      big_box_level.getGridGeometry());

   tbox::plog << "shrinkBoxLevel: Boundary plain boxes:\n"
              << boundary_boxes.format("\n", 2);

   /*
    * Construct the complement of the small_box_level by
    * growing the boundary boxes.
    */

   hier::BoxContainer complement_boxes;

   hier::LocalId last_local_id(-1);
   for (hier::BoxContainer::const_iterator bi = boundary_boxes.begin();
        bi != boundary_boxes.end(); ++bi) {
      hier::Box box(*bi);
      box.grow(shrinkage);
      hier::Box complement_box(
         box, ++last_local_id, local_rank);
      complement_boxes.insert(complement_box);
   }

   complement_boxes.makeTree(grid_geometry.get());

   /*
    * Construct the small_box_level.
    */

   small_box_level.reset(new hier::BoxLevel(
         big_box_level.getRefinementRatio(),
         grid_geometry,
         big_box_level.getMPI()));
   last_local_id = -1;
   for (hier::BoxContainer::const_iterator bi = big_boxes.begin();
        bi != big_boxes.end(); ++bi) {

      const hier::Box& box = *bi;

      int ix;
      for (ix = 0; ix < static_cast<int>(unshrunken_blocks.size()); ++ix) {
         if (box.getBlockId() == unshrunken_blocks[ix]) {
            break;
         }
      }

      if (ix < static_cast<int>(unshrunken_blocks.size())) {
         /*
          * This block should be excluded from shrinking.
          */
         small_box_level->addBoxWithoutUpdate(box);
      } else {

         hier::BoxContainer shrunken_boxes(box);

         shrunken_boxes.removeIntersections(
            big_box_level.getRefinementRatio(),
            complement_boxes);
         shrunken_boxes.simplify();

         for (hier::BoxContainer::iterator li = shrunken_boxes.begin();
              li != shrunken_boxes.end(); ++li) {
            const hier::Box shrunken_box(
               *li,
               ++last_local_id,
               box.getOwnerRank());
            TBOX_ASSERT(shrunken_box.getBlockId() ==
               box.getBlockId());

            small_box_level->addBoxWithoutUpdate(shrunken_box);
         }
      }

   }
   small_box_level->finalize();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void refineBoxLevel(hier::BoxLevel& box_level,
                    const hier::IntVector& ratio)
{
   box_level.refineBoxes(
      box_level,
      ratio,
      box_level.getRefinementRatio() * ratio);
   box_level.finalize();
}
