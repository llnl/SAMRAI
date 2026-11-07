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
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
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
 * Break up boxes in the given BoxLevel.  This method is meant
 * to create a bunch of small boxes from the user-input boxes in order
 * to set up a non-trivial mesh configuration.
 *
 * The database may have the following items: min_box_size,
 * max_box_size, refinement_ratio.
 *
 * The box_level will be refined by refinement_ratio and
 * partitioned to generate a non-trivial configuration for testing.
 */
void
breakUpBoxes(
   hier::BoxLevel& box_level,
   const hier::BoxLevel& domain_box_level,
   const std::shared_ptr<tbox::Database>& database);

void
alterAndGenerateMapping(
   std::shared_ptr<hier::BoxLevel>& box_level_c,
   std::shared_ptr<hier::MappingConnector>& b_to_c,
   const hier::BoxLevel& box_level_b,
   const std::shared_ptr<tbox::Database>& database);

/*
 ************************************************************************
 *
 * This is an accuracy test for the MappingConnectorAlgorithm class:
 *
 * 1. Read in user-specified GridGeometry.
 *
 * 2. Build a domain BoxLevel from GridGeometry.
 *
 * 3. Build BoxLevel A by refining and partitioning domain.
 *
 * 4. Build BoxLevel B by refining and partitioning domain.
 *    Compute overlap Connectors A<==>B.
 *
 * 5. Build BoxLevel C by changing B based on some simple formula.
 *    Generate mapping Connectors B<==>C.
 *
 * 6. Apply mapping B<==>C to update A<==>B.
 *
 * 7. Check correctness of updated A<==>B.
 *
 *************************************************************************
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

   size_t fail_count = 0;

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

      plog << "Input database after initialization..." << std::endl;
      input_db->printClassData(plog);

      /*
       * Generate the grid geometry.
       */
      if (!main_db->keyExists("GridGeometry")) {
         TBOX_ERROR("Multiblock tree search test: could not find entry GridGeometry"
            << "\nin input.");
      }
      std::shared_ptr<hier::BaseGridGeometry> grid_geometry(
         new geom::GridGeometry(
            dim,
            "GridGeometry",
            main_db->getDatabase("GridGeometry")));

      /*
       * Print input database again to fully show usage.
       */
      plog << "Input database after running..." << std::endl;
      input_db->printClassData(plog);

      const hier::IntVector& one_vector(hier::IntVector::getOne(dim));
      const hier::IntVector& zero_vector(hier::IntVector::getZero(dim));

      hier::BoxLevel domain_box_level(
         one_vector,
         grid_geometry,
         tbox::SAMRAI_MPI::getSAMRAIWorld(),
         hier::BoxLevel::GLOBALIZED);
      grid_geometry->computePhysicalDomain(
         domain_box_level,
         one_vector);
      domain_box_level.finalize();

      std::vector<hier::IntVector> refinement_ratios(
         1, domain_box_level.getRefinementRatio());

      /*
       * Generate BoxLevel A from the multiblock domain description
       * using input database BoxLevelA.
       */
      hier::BoxLevel box_level_a(domain_box_level);
      std::shared_ptr<Database> a_db(main_db->getDatabase("BoxLevelA"));
      breakUpBoxes(box_level_a, domain_box_level, a_db);
      box_level_a.cacheGlobalReducedData();
      refinement_ratios.push_back(box_level_a.getRefinementRatio());

      /*
       * Generate BoxLevel B from the multiblock domain description
       * using input database BoxLevelB.
       */
      hier::BoxLevel box_level_b(domain_box_level);
      std::shared_ptr<Database> b_db(main_db->getDatabase("BoxLevelB"));
      breakUpBoxes(box_level_b, domain_box_level, b_db);
      box_level_b.cacheGlobalReducedData();
      refinement_ratios.push_back(box_level_b.getRefinementRatio() /
                                  box_level_a.getRefinementRatio());

      /*
       * These steps are usually handled by PatchHierarchy, but this
       * test does not use PatchHierarchy.
       */
      grid_geometry->setUpRatios(refinement_ratios);

      /*
       * Generate Connector A<==>B, to be modified by the mapping
       * operation.
       */

      hier::IntVector base_width_a(zero_vector);
      if (main_db->isInteger("base_width_a")) {
         main_db->getIntegerArray("base_width_a", &base_width_a[0], dim.getValue());
      }
      hier::IntVector width_a(base_width_a);
      hier::IntVector width_b(
         hier::Connector::convertHeadWidthToBase(
            box_level_b.getRefinementRatio(),
            box_level_a.getRefinementRatio(),
            width_a));

      std::shared_ptr<hier::Connector> a_to_b;

      hier::OverlapConnectorAlgorithm oca;
      oca.findOverlapsWithTranspose(a_to_b,
         box_level_a,
         box_level_b,
         width_a,
         width_b);
      hier::Connector& b_to_a = a_to_b->getTranspose();
      // tbox::pout << "a_to_b:\n" << a_to_b->format("AB: ",2) << std::endl;
      // tbox::pout << "b_to_a:\n" << b_to_a.format("BA: ",2) << std::endl;

      a_to_b->checkConsistencyWithBase();
      a_to_b->checkConsistencyWithHead();
      b_to_a.checkConsistencyWithBase();
      b_to_a.checkConsistencyWithHead();

      a_to_b->checkOverlapCorrectness();
      b_to_a.checkOverlapCorrectness();

      /*
       * Generate BoxLevel C by altering B based on a simple formula.
       * Generate the mapping Connectors B<==>C.
       */

      std::shared_ptr<hier::BoxLevel> box_level_c;
      std::shared_ptr<hier::MappingConnector> b_to_c;
      std::shared_ptr<Database> alteration_db(
         main_db->getDatabase("Alteration"));

      alterAndGenerateMapping(
         box_level_c,
         b_to_c,
         box_level_b,
         alteration_db);
      box_level_c->cacheGlobalReducedData();
      // tbox::pout << "box level c:\n" << box_level_c->format("C: ",2) << std::endl;
      // tbox::pout << "b_to_c:\n" << b_to_c->format("BC: ",2) << std::endl;
      // tbox::pout << "c_to_b:\n" << b_to_c->getTranspose().format("CB: ",2) << std::endl;

      hier::MappingConnectorAlgorithm mca;
      mca.modify(*a_to_b,
         *b_to_c,
         &box_level_b,
         box_level_c.get());
      // tbox::pout << "box level b after modify:\n" << box_level_b.format("B: ",2) << std::endl;

      // tbox::pout << "checking a--->b consistency with base:" << std::endl;
      a_to_b->checkConsistencyWithBase();
      // tbox::pout << "checking a--->b consistency with head:" << std::endl;
      a_to_b->checkConsistencyWithHead();
      // tbox::pout << "checking b--->a consistency with base:" << std::endl;
      b_to_a.checkConsistencyWithBase();
      // tbox::pout << "checking b--->a consistency with head:" << std::endl;
      b_to_a.checkConsistencyWithHead();

      tbox::pout << "Checking for a--->b overlap correctness:" << std::endl;
      const int a_to_b_errors = a_to_b->checkOverlapCorrectness();
      if (a_to_b_errors) {
         tbox::pout << "... " << a_to_b_errors << " errors." << std::endl;
      } else {
         tbox::pout << "... OK." << std::endl;
      }

      tbox::pout << "Checking for b--->a overlap correctness:" << std::endl;
      const int b_to_a_errors = b_to_a.checkOverlapCorrectness();
      if (b_to_a_errors) {
         tbox::pout << "... " << b_to_a_errors << " errors." << std::endl;
      } else {
         tbox::pout << "... OK." << std::endl;
      }

      fail_count += a_to_b_errors + b_to_a_errors;

      if (fail_count == 0) {
         tbox::pout << "\nPASSED:  Connector modify" << std::endl;
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

   return int(fail_count);
}

/*
 * Break up boxes in the given BoxLevel.  This method is meant
 * to create a bunch of small boxes from the user-input boxes in order
 * to set up a non-trivial mesh configuration.
 *
 * 1. Refine the boxes according to refinement_ratio in the database.
 * 2. Partition according to min and max box sizes in the database.
 */
void breakUpBoxes(
   hier::BoxLevel& box_level,
   const hier::BoxLevel& domain_box_level,
   const std::shared_ptr<tbox::Database>& database) {

   const tbox::Dimension& dim(box_level.getDim());

   hier::IntVector refinement_ratio(hier::IntVector::getOne(dim));
   if (database->isInteger("refinement_ratio")) {
      database->getIntegerArray("refinement_ratio", &refinement_ratio[0], dim.getValue());
   }

   if (refinement_ratio != hier::IntVector::getOne(dim)) {
      box_level.refineBoxes(box_level,
         refinement_ratio,
         box_level.getRefinementRatio() * refinement_ratio);
      box_level.finalize();
   }

   hier::IntVector max_box_size(dim, tbox::MathUtilities<int>::getMax());
   if (database->isInteger("max_box_size")) {
      database->getIntegerArray("max_box_size", &max_box_size[0], dim.getValue());
   }

   hier::IntVector min_box_size(hier::IntVector::getOne(dim));
   if (database->isInteger("min_box_size")) {
      database->getIntegerArray("min_box_size", &min_box_size[0], dim.getValue());
   }

   mesh::TreeLoadBalancer load_balancer(box_level.getDim(),
                                        "TreeLoadBalancer");

   const int level_number(0);

   hier::Connector* dummy_connector = 0;

   const hier::IntVector bad_interval(hier::IntVector::getOne(dim));
   const hier::IntVector cut_factor(hier::IntVector::getOne(dim));

   load_balancer.loadBalanceBoxLevel(
      box_level,
      dummy_connector,
      std::shared_ptr<hier::PatchHierarchy>(),
      level_number,
      min_box_size,
      max_box_size,
      domain_box_level,
      bad_interval,
      cut_factor);
}

/*
 * Generate BoxLevel C by altering B based on a simple formula.
 * Generate the mapping Connectors B<==>C.
 */
void alterAndGenerateMapping(
   std::shared_ptr<hier::BoxLevel>& box_level_c,
   std::shared_ptr<hier::MappingConnector>& b_to_c,
   const hier::BoxLevel& box_level_b,
   const std::shared_ptr<tbox::Database>& database)
{
   const tbox::Dimension dim(box_level_b.getDim());

   /*
    * Increment for changing the LocalIds.
    * Set to zero to disable.
    */
   const int local_id_increment =
      database->getIntegerWithDefault("local_id_increment", 0);

   const hier::BoxContainer boxes_b(box_level_b.getBoxes());

   box_level_c.reset(new hier::BoxLevel(box_level_b.getRefinementRatio(),
         box_level_b.getGridGeometry(),
         box_level_b.getMPI()));

   b_to_c.reset(new hier::MappingConnector(box_level_b,
         *box_level_c,
         hier::IntVector::getZero(dim)));
   hier::MappingConnector* c_to_b =
      new hier::MappingConnector(*box_level_c,
         box_level_b,
         hier::IntVector::getZero(dim));
   b_to_c->setTranspose(c_to_b, true);
   for (hier::BoxContainer::const_iterator bi = boxes_b.begin();
        bi != boxes_b.end(); ++bi) {
      const hier::Box& box_b(*bi);
      hier::Box box_c(box_b,
                      box_b.getLocalId() + local_id_increment,
                      box_b.getOwnerRank(),
                      box_b.getPeriodicId());
      box_level_c->addBoxWithoutUpdate(box_c);
      b_to_c->insertLocalNeighbor(box_c, box_b.getBoxId());
      c_to_b->insertLocalNeighbor(box_b, box_c.getBoxId());
   }

   box_level_c->finalize();

   b_to_c->checkConsistencyWithBase();
   b_to_c->checkConsistencyWithHead();
   c_to_b->checkConsistencyWithBase();
   c_to_b->checkConsistencyWithHead();

   b_to_c->assertMappingValidity();
   c_to_b->assertMappingValidity();
}
