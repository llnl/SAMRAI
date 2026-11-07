/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Example program to demonstrate boundary utilities.
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include <string>
#include <memory>

// Headers for basic SAMRAI objects used in this code.
#include "SAMRAI/tbox/SAMRAIManager.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Utilities.h"

// Headers for classes specific to this example
#include "BoundaryDataTester.h"

using namespace SAMRAI;

int main(
   int argc,
   char* argv[])
{
   int fail_count = -1;

   /*
    * Initialize tbox::MPI and SAMRAI, enable logging, and process command line.
    * Note this example is set up to run in serial only.
    */

   tbox::SAMRAI_MPI::init(&argc, &argv);

   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   {

      if (argc != 2) {
         TBOX_ERROR(
            "USAGE:  " << argv[0] << " <input filename> "
                       << "<restart dir> <restore number> [options]\n"
                       << "  options:\n"
                       << "  none at this time"
                       << std::endl);
         return -1;
      }

      std::string input_filename = argv[1];

      /*
       * Create input database and parse all data in input file.
       */

      std::shared_ptr<tbox::InputDatabase> input_db(
         new tbox::InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Retrieve "GlobalInputs" section of the input database and set
       * values accordingly.
       */

      if (input_db->keyExists("GlobalInputs")) {
         std::shared_ptr<tbox::Database> global_db(
            input_db->getDatabase("GlobalInputs"));
         if (global_db->keyExists("call_abort_in_serial_instead_of_exit")) {
            bool flag = global_db->
               getBool("call_abort_in_serial_instead_of_exit");
            tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(flag);
         }
      }

      /*
       * Read "Main" input data.
       */

      std::shared_ptr<tbox::Database> main_db(input_db->getDatabase("Main"));

      const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      std::string log_file_name = "boundary.log";
      if (main_db->keyExists("log_file_name")) {
         log_file_name = main_db->getString("log_file_name");
      }
      tbox::PIO::logOnlyNodeZero(log_file_name);

      hier::IntVector num_boxes(dim, 1);
      if (main_db->keyExists("num_domain_boxes")) {
         int* tmp_arr = &num_boxes[0];
         main_db->getIntegerArray("num_domain_boxes", tmp_arr, dim.getValue());
      }

      /*
       * Create objects used in boundary data test.  Then, print out
       * state of BoundaryDataTester to log file for checking.
       */

      std::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
         new geom::CartesianGridGeometry(
            dim,
            "CartesianGridGeometry",
            input_db->getDatabase("CartesianGridGeometry")));

      std::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
         new hier::PatchHierarchy(
            "PatchHierarchy",
            grid_geometry));

      BoundaryDataTester* btester =
         new BoundaryDataTester(
            "BoundaryDataTester",
            dim,
            input_db->getDatabase("BoundaryDataTester"),
            grid_geometry);

      tbox::plog
      << "\nPRINTING BoundaryDataTester object state after initialization..."
      << std::endl;
      btester->printClassData(tbox::plog);

      /*
       * For simplicity, we manually create a hierachy with a single patch level.
       */

      tbox::plog << "\nBuilding patch hierarchy..." << std::endl;

      const hier::BoxContainer& domain = grid_geometry->getPhysicalDomain();
      hier::BoxContainer boxes(domain);
      boxes.unorder();
      if ((domain.size() == 1) &&
          (num_boxes != hier::IntVector(dim, 1))) {
         const hier::Box& dbox = domain.front();
         hier::IntVector max_size(dbox.numberCells());
         hier::IntVector min_size(dbox.numberCells() / num_boxes);
         hier::IntVector cut_factor(dim, 1);
         hier::IntVector bad_interval(dim, 1);
         hier::BoxUtilities::chopBoxes(boxes,
            max_size,
            min_size,
            cut_factor,
            bad_interval,
            domain);
      }

      hier::BoxLevelConnectorUtils edge_utils;
      std::shared_ptr<hier::BoxLevel> layer0(
         std::make_shared<hier::BoxLevel>(
            hier::IntVector(dim, 1), grid_geometry));
      hier::BoxContainer::const_iterator domain_boxes = domain.begin();
      int rank = mpi.getRank();
      int size = mpi.getSize();
      for (hier::LocalId ib(0); ib < boxes.size(); ++ib, ++domain_boxes) {
         if (ib % size == rank) {
            layer0->addBox(hier::Box(*domain_boxes, ib, rank));
         }
      }
      edge_utils.addPeriodicImages(
         *layer0,
         patch_hierarchy->getGridGeometry()->getDomainSearchTree(),
         hier::IntVector(dim, 2));

      patch_hierarchy->makeNewPatchLevel(0, layer0);

      // Add Connector required for schedule construction.
      std::shared_ptr<hier::PatchLevel> level0(
         patch_hierarchy->getPatchLevel(0));
      level0->createConnector(*level0, hier::IntVector(dim, 2));

      /*
       * Allocate data on hierarchy and set variable data on patch interiors
       * to input values.
       */

      tbox::plog << "\nAllocate and initialize data on patch hierarchy..."
                 << std::endl;

      btester->initializeDataOnPatchInteriors(patch_hierarchy, 0);

      tbox::plog << "Performing tests..." << std::endl;

      fail_count = btester->runBoundaryTest(patch_hierarchy, 0);

      tbox::plog << "\n\n\nDone." << std::endl;

      /*
       * At conclusion of test, deallocate objects.
       */
      patch_hierarchy.reset();
      grid_geometry.reset();

      if (btester) delete btester;

      if (fail_count == 0) {
         tbox::pout << "\nPASSED:  boundary test" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return fail_count;
}
