/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test program to test OverlapConnectorAlgorithm class
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/AssumedPartition.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/geom/GridGeometry.h"



using namespace SAMRAI;
using namespace hier;

/*!
 * @brief Primitive BoxGenerator (independent of mesh package)
 * creating boxes using an AssumedPartition followed by an index
 * filter to keep a subset of boxes.
 */
struct PrimitiveBoxGen {
   std::shared_ptr<BaseGridGeometry> d_geom;
   AssumedPartition d_ap;
   // Index filtering parameters.
   enum IndexFilter { ALL = 0 /* Keep all boxes */,
                      INTERVAL = 1      /* Keep d_num_keep, discard d_num_discard */
   };
   int d_index_filter;
   int d_num_keep;
   int d_num_discard;
   PrimitiveBoxGen(
      tbox::Database& database,
      const std::shared_ptr<BaseGridGeometry>& geom):
      d_index_filter(ALL),
      d_num_keep(1),
      d_num_discard(0)
   {
      d_geom = geom;
      getFromInput(database);
   }
   PrimitiveBoxGen(
      const PrimitiveBoxGen& other):
      d_geom(other.d_geom),
      d_ap(other.d_ap),
      d_index_filter(other.d_index_filter),
      d_num_keep(other.d_num_keep),
      d_num_discard(other.d_num_discard) {
   }
   void
   getFromInput(
      tbox::Database& input_db);
   void
   getBoxes(
      BoxContainer& boxes,
      int rank);
   void
   populateBoxLevel(
      BoxLevel& box_level);
};

struct CommonTestParams {
   PrimitiveBoxGen d_gen1;
   PrimitiveBoxGen d_gen2;
   CommonTestParams(
      const tbox::Dimension& dim);
   CommonTestParams(
      const CommonTestParams& other);
};

CommonTestParams
getTestParametersFromDatabase(
   tbox::Database& database);

int main(
   int argc,
   char* argv[])
{
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   const tbox::SAMRAI_MPI& mpi = tbox::SAMRAI_MPI::getSAMRAIWorld();

   int fail_count = 0;

   /*
    * Process command line arguments.  For each run, the input
    * filename must be specified.  Usage is:
    *
    * executable <input file name>
    */
   std::string input_filename;

   if (argc < 2) {
      TBOX_ERROR("USAGE:  " << argv[0] << " <input file> [case name]\n"
                            << "  options:\n"
                            << "  none at this time" << std::endl);
   } else {
      input_filename = argv[1];
   }

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {
      /*
       * Create input database and parse all data in input file.
       */
      std::shared_ptr<tbox::MemoryDatabase> input_db(new tbox::MemoryDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      std::shared_ptr<tbox::Database> main_db = input_db->getDatabase("Main");

      std::string base_name = "unnamed";
      base_name = main_db->getStringWithDefault("base_name", base_name);

      /*
       * Modify basename for this particular run.
       * Add the number of processes and the case name.
       */
      std::string base_name_ext = base_name;
      base_name_ext = base_name_ext + '-'
         + tbox::Utilities::nodeToString(mpi.getSize());

      /*
       * Start logging.
       */
      const std::string log_file_name = base_name_ext + ".log";
      bool log_all_nodes = false;
      log_all_nodes = main_db->getBoolWithDefault("log_all_nodes", log_all_nodes);
      if (log_all_nodes) {
         tbox::PIO::logAllNodes(log_file_name);
      } else {
         tbox::PIO::logOnlyNodeZero(log_file_name);
      }

      const int rank = mpi.getRank();

      {

         const tbox::Dimension dim(static_cast<tbox::Dimension::dir_t>(main_db->getInteger("dim")));

         if (!input_db->isDatabase("BlockGeometry")) {
            TBOX_ERROR(
               "getTestParametersFromDatabase: You must specify \"BlockGeometry\" in input database.");
         }
         // Note: Using GridGeometry only because BaseGridGeometry can't be instanstiated.
         std::shared_ptr<BaseGridGeometry> grid_geom =
            std::make_shared<geom::GridGeometry>(
               dim,
               "BlockGeometry",
               input_db->getDatabase("BlockGeometry"));

         int test_number = 0;
         while (true) {

            std::string test_name("Test");
            test_name += tbox::Utilities::intToString(test_number++, 2);

            std::shared_ptr<tbox::Database> test_db =
               input_db->getDatabaseWithDefault(test_name, std::shared_ptr<tbox::Database>());

            if (!test_db) {
               break;
            }

            const std::string nickname =
               test_db->getStringWithDefault("nickname", test_name);

            tbox::plog << "\n\n";
            tbox::pout << "Running " << test_name << " (" << nickname << ")\n";

            PrimitiveBoxGen pb1(*test_db->getDatabase("PrimitiveBoxGen1"), grid_geom);
            BoxContainer boxes1;
            pb1.getBoxes(boxes1, rank);
            BoxLevel l1(boxes1, IntVector::getOne(pb1.d_geom->getDim()), pb1.d_geom);
            l1.cacheGlobalReducedData();

            PrimitiveBoxGen pb2(*test_db->getDatabase("PrimitiveBoxGen2"), grid_geom);
            BoxContainer boxes2;
            pb2.getBoxes(boxes2, rank);
            BoxLevel l2(boxes2, IntVector::getOne(pb2.d_geom->getDim()), pb2.d_geom);
            l2.cacheGlobalReducedData();

            /*
             * Set up edges in l1_to_l2 by the contrivance specified
             * in the test database.  Then check transpose
             * correctness.
             */

            Connector l1_to_l2(l1, l2, IntVector::getZero(dim));
            Connector l2_to_l1(l2, l1, IntVector::getZero(dim));
            OverlapConnectorAlgorithm oca;
            oca.findOverlaps_assumedPartition(l1_to_l2);
            oca.findOverlaps_assumedPartition(l2_to_l1);

            tbox::plog << "Testing with:"
                       << "\nl1:\n" << l1.format("\t")
                       << "\nl2:\n" << l2.format("\t")
                       << "\nl1_to_l2:\n" << l1_to_l2.format("\t")
                       << "\nl2_to_l1:\n" << l2_to_l1.format("\t")
                       << std::endl;

            size_t fail_count_1 = l1_to_l2.checkOverlapCorrectness();
            size_t fail_count_2 = l2_to_l1.checkOverlapCorrectness();

            if (fail_count_1) {
               tbox::plog << "Error l1_to_l2 of " << test_name << " (" << nickname << ')'
                          << std::endl;
            }
            if (fail_count_2) {
               tbox::plog << "Error l2_to_l1 of " << test_name << " (" << nickname << ')'
                          << std::endl;
            }

            if (fail_count_1 || fail_count_2) {
               tbox::pout << "FAILED: " << test_name << " (" << nickname << ')' << std::endl;
               tbox::plog << "FAILED: " << test_name << " (" << nickname << ')' << std::endl;
            } else {
               tbox::plog << "PASSED: " << test_name << " (" << nickname << ')' << std::endl;
            }

            fail_count += static_cast<int>(fail_count_1 + fail_count_2);
         }

      }

      input_db->printClassData(tbox::plog);

   }

   if (fail_count == 0) {
      tbox::pout << "\nPASSED:  Connector" << std::endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();
   return fail_count;
}

/*
 *************************************************************************
 *************************************************************************
 */
void PrimitiveBoxGen::getFromInput(tbox::Database& database)
{
   int rank_begin = 0;
   int rank_end = tbox::SAMRAI_MPI::getSAMRAIWorld().getSize();
   int index_begin = database.getIntegerWithDefault("index_begin", 0);
   double parts_per_rank = database.getDoubleWithDefault("parts_per_rank", 1.0);
   d_ap.partition(d_geom->getPhysicalDomain(), rank_begin, rank_end, index_begin, parts_per_rank);
   tbox::plog << "PrimitiveBoxGen::getFromInput() generated AssumedPartition:\n";
   d_ap.recursivePrint(tbox::plog, "\t", 3);
   if (d_ap.selfCheck()) {
      TBOX_ERROR("Error in setting up AssumedPartition d_ap (selfCheck failed).\n");
   }

   std::string index_filter = database.getStringWithDefault("index_filter", "ALL");
   if (index_filter == "ALL") {
      d_index_filter = PrimitiveBoxGen::ALL;
   } else if (index_filter == "INTERVAL") {
      d_index_filter = PrimitiveBoxGen::INTERVAL;
   }
   d_num_keep = database.getIntegerWithDefault("num_keep", d_num_keep);
   d_num_discard = database.getIntegerWithDefault("num_discard", d_num_discard);
}

/*
 *************************************************************************
 *************************************************************************
 */
void PrimitiveBoxGen::getBoxes(BoxContainer& boxes, int rank)
{
   if (d_index_filter == ALL) {
      int idbegin = d_ap.beginOfRank(rank);
      int idend = d_ap.endOfRank(rank);
      for (int id = idbegin; id < idend; ++id) {
         boxes.push_back(d_ap.getBox(id));
      }
   } else if (d_index_filter == INTERVAL) {
      int idbegin = d_ap.beginOfRank(rank);
      int idend = d_ap.endOfRank(rank);
      int interval = d_num_keep + d_num_discard;
      for (int id = idbegin; id < idend; ++id) {
         int interval_id = id % interval;
         if (interval_id < d_num_keep) {
            boxes.push_back(d_ap.getBox(id));
         }
      }
   } else {
      TBOX_ERROR("Invalid value of index_filter: " << d_index_filter);
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void PrimitiveBoxGen::populateBoxLevel(BoxLevel& box_level)
{
   BoxContainer boxes;
   getBoxes(boxes, box_level.getMPI().getRank());
   for (BoxContainer::const_iterator bi = boxes.begin(); bi != boxes.end(); ++bi) {
      box_level.addBoxWithoutUpdate(*bi);
   }
   box_level.finalize();
}
