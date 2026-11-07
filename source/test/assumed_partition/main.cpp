/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test program to test the AssumedPartition classes
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/AssumedPartitionBox.h"
#include "SAMRAI/hier/AssumedPartition.h"
#include "SAMRAI/geom/GridGeometry.h"


using namespace SAMRAI;

struct CommonTestParams {
   hier::Box box;
   std::shared_ptr<hier::BaseGridGeometry> geometry;
   int rank_begin;
   int rank_end;
   int index_begin;
   double avg_parts_per_rank;
   CommonTestParams(
      const tbox::Dimension& dim):
      box(dim),
      rank_begin(0),
      rank_end(tbox::SAMRAI_MPI::getSAMRAIWorld().getSize()),
      index_begin(0),
      avg_parts_per_rank(1.0) {
   }
   CommonTestParams(
      const CommonTestParams& other):
      box(other.box),
      geometry(other.geometry),
      rank_begin(other.rank_begin),
      rank_end(other.rank_end),
      index_begin(other.index_begin),
      avg_parts_per_rank(other.avg_parts_per_rank) {
   }
};

CommonTestParams
getTestParametersFromDatabase(
   tbox::Database& test_db);

int
main(
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
      log_all_nodes = main_db->getBoolWithDefault("log_all_nodes",
            log_all_nodes);
      if (log_all_nodes) {
         tbox::PIO::logAllNodes(log_file_name);
      } else {
         tbox::PIO::logOnlyNodeZero(log_file_name);
      }

      {

         int test_number = 0;

         while (true) {

            std::string test_name("Test");
            test_name += tbox::Utilities::intToString(test_number, 2);

            std::shared_ptr<tbox::Database> test_db =
               input_db->getDatabaseWithDefault(test_name, std::shared_ptr<tbox::Database>());

            if (!test_db) {
               break;
            }

            const std::string nickname =
               test_db->getStringWithDefault("nickname", test_name);

            tbox::plog << "\n\n";
            tbox::pout << "Running " << test_name << " (" << nickname << ")\n";

            CommonTestParams ctp = getTestParametersFromDatabase(*test_db);

            size_t test_fail_count = 0;
            if (ctp.geometry) {
               // Test multi-box assumed partitions.
               hier::AssumedPartition ap(ctp.geometry->getPhysicalDomain(),
                                         ctp.rank_begin,
                                         ctp.rank_end,
                                         ctp.index_begin,
                                         ctp.avg_parts_per_rank);
               tbox::plog << "AssumedPartition:\n";
               ap.recursivePrint(tbox::plog, "\t");
               test_fail_count = ap.selfCheck();
            } else {
               // Test single-box assumed partitions.
               hier::AssumedPartitionBox apb(ctp.box,
                                             ctp.rank_begin,
                                             ctp.rank_end,
                                             ctp.index_begin,
                                             ctp.avg_parts_per_rank);
               tbox::plog << "AssumedPartitionBox:\n";
               apb.recursivePrint(tbox::plog, "\t");
               test_fail_count = apb.selfCheck();
            }

            fail_count += static_cast<int>(test_fail_count);
            if (test_fail_count) {
               tbox::pout << "FAILED: selfCheck found " << fail_count << " problems in test "
                          << test_name << " (" << nickname << ')' << std::endl;
            }

            ++test_number;

         }
      }

   }

   if (fail_count == 0) {
      tbox::pout << "\nPASSED:  assumed_partition" << std::endl;
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
CommonTestParams
getTestParametersFromDatabase(
   tbox::Database& test_db)
{
   const tbox::Dimension dim(static_cast<tbox::Dimension::dir_t>(test_db.getInteger("dim")));
   CommonTestParams ctp(dim);
   if (test_db.isDatabaseBox("box")) {
      ctp.box = test_db.getDatabaseBox("box");
      ctp.box.setBlockId(hier::BlockId(0));
   } else if (test_db.isDatabase("BlockGeometry")) {
      // Note: Using GridGeometry only because BaseGridGeometry can't be instanstiated.
      ctp.geometry.reset(
         new geom::GridGeometry(
            dim,
            "BlockGeometry",
            test_db.getDatabase("BlockGeometry")));
   } else {
      TBOX_ERROR("getTestParametersFromDatabase: You must specify either \"box\"\n"
         << "or \"BlockGeometry\" in each test.");
   }
   ctp.rank_begin = test_db.getIntegerWithDefault("rank_begin", ctp.rank_begin);
   ctp.rank_end = test_db.getIntegerWithDefault("rank_end", ctp.rank_end);
   ctp.index_begin = test_db.getIntegerWithDefault("index_begin", ctp.index_begin);
   ctp.avg_parts_per_rank = test_db.getDoubleWithDefault("avg_parts_per_rank",
         ctp.avg_parts_per_rank);
   return ctp;
}
