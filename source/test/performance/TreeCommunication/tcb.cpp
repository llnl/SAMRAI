/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Code and input for benchmarking and experimentation with tree-based communication.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include <iomanip>

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/tbox/CenteredRankTree.h"
#include "SAMRAI/tbox/Clock.h"
#include "SAMRAI/tbox/CommGraphWriter.h"
#include "SAMRAI/tbox/BreadthFirstRankTree.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/RankTreeStrategy.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/TimerManager.h"

using namespace SAMRAI;
using namespace tbox;

/*
 ************************************************************************
 * Tree communication benchmark code.
 *************************************************************************
 */

struct CommonTestSwitches {
   std::string tree_name;
   std::string message_pattern;
   int msg_length;
   int first_data_length;
   int processing_cost[2];
   bool randomize_processing_cost;
   bool verify_data;
   int repetition;
   bool barrier_after_each_repetition;
   int mpi_tags[2];
   CommonTestSwitches():
      tree_name(),
      message_pattern(),
      msg_length(1024),
      first_data_length(1),
      randomize_processing_cost(false),
      verify_data(false),
      repetition(1),
      barrier_after_each_repetition(false)
   {
      processing_cost[0] = 0;
      processing_cost[1] = 0;
      mpi_tags[0] = 0;
      mpi_tags[1] = 1;
   }
};

std::shared_ptr<RankTreeStrategy>
getTreeForTesting(
   const std::string& tree_name,
   Database& test_db,
   const SAMRAI_MPI& mpi);

void
setupAsyncComms(
   AsyncCommStage& child_stage,
   AsyncCommPeer<int> *& child_comms,
   AsyncCommStage& parent_stage,
   AsyncCommPeer<int> *& parent_comm,
   const SAMRAI_MPI& mpi,
   const RankTreeStrategy& rank_tree,
   const CommonTestSwitches& cts);

void
destroyAsyncComms(
   AsyncCommPeer<int> *& child_comms,
   AsyncCommPeer<int> *& parent_comm);

int
testUp(
   const RankTreeStrategy& rank_tree,
   const SAMRAI_MPI& mpi,
   const CommonTestSwitches& cts,
   Database& test_db);

int
testDown(
   const RankTreeStrategy& rank_tree,
   const SAMRAI_MPI& mpi,
   const CommonTestSwitches& cts,
   Database& test_db);

int
testUpThenDown(
   const RankTreeStrategy& rank_tree,
   const SAMRAI_MPI& mpi,
   const CommonTestSwitches& cts,
   Database& test_db);

int
testTreeLB(
   const RankTreeStrategy& rank_tree,
   const SAMRAI_MPI& mpi,
   const CommonTestSwitches& cts,
   Database& test_db);

void
setMessageData(
   std::vector<int>& msg,
   int rank);

int
verifyReceivedData(
   const AsyncCommPeer<int>& peer_comm);

void
getCommonTestSwitchesFromDatabase(
   CommonTestSwitches& cts,
   Database& test_db);

void
simulateDataProcessing(
   const CommonTestSwitches& cts);

SAMRAI_MPI::Comm
getRotatedMPI(
   const SAMRAI_MPI::Comm& old_comm);
SAMRAI_MPI::Comm
getSmallerMPI(
   const SAMRAI_MPI::Comm& old_comm);

// Replacement for usleep, which is broken on some machines.
void
my_usleep(
   size_t num_usec);
void
calibrate_my_usleep();
unsigned long dummy() {
   static unsigned long i = 0;
   i = (i + 1) % 1000;
   return i;
}
double my_usleep_calls_per_usec = 0;

/*
 ********************************************************************************
 *
 * Performance testing for tree communication.
 *
 * 1.
 *
 * 2.
 *
 * Input File:
 *
 * Test## { // ## is a 2-digit integer, sequentially from 0
 *
 *   nickname = "Foobar" // Optional name for this test.
 *
 *   tree_name = "BalancedDepthFirstTree" // BalancedDepthFirstTree || CenteredRankTree || ...
 *
 *   BalancedDepthFirstTree { // Parameters for BalancedDepthFirstTree in getTreeForTesting()
 *     do_left_leaf_switch = TRUE
 *   }
 *
 *   // Pattern of message travel:
 *   // "UP", "DOWN": Up or down the tree.
 *   // "UP_THEN_DOWN", "DOWN_THEN_UP": Self-explanatory
 *   // "TreeLB": Simulate the communication of the TreeLoadBalancer.
 *   message_pattern = "DOWN"
 *
 *   msg_length = 1024 // Message length (units of integer)
 *   first_data_length = 1 // One int.  see AsyncCommPeer::limitFirstDataLength().
 *
 *   repetition = 1 // Number of times to run.  0 disables test.
 *   barrier_after_each_repetition = FALSE // Whether to barrier after each rep.
 *
 *   processing_cost = 5, 2 // Simulate processing with this 5 microseconds per message and 2 per msg_length unit.
 *   randomize_processing_cost = FALSE // Multiply a random factor in [0,1] to the processing cost.
 *   verify_data = TRUE // Verify correctness of received data.
 *
 *   mpi_tags = 1, 2 // Array of 2 ints, see AsyncCommPeer::setMPITag().
 *
 *   // Specify the dependency for the down-message, as a funcion of MPI rank:
 *   // 1: down message depends only on parent
 *   // 2: down message depends on grandparent
 *   // 0: there is no down message
 *   // First value is for rank 0, second is for rank 1, and so on.
 *   // The array is repeated for ranks higher than specifed.
 *   // (Rank r has dependency according to index r%L, where L is
 *   // the length of down_message_dependency.)
 *   down_message_dependency = 1, 2, 0, 1, 2, 0
 * }
 *
 ********************************************************************************
 */

int main(
   int argc,
   char* argv[])
{
#ifndef HAVE_MPI
   // This test doesn't make sense without MPI because it cannot avoid MPI interfaces.
   NULL_USE(argc);
   std::cout << "PASSED: " << argv[0] << std::endl;
   return 0;

#else

   /*
    * The switch --alter_mpi=<string> is an option to alter SAMRAI's
    * MPI group.  It must be processed before starting SAMRAI because
    * it affects SAMRAI's initial MPI.  Valid values for <string> and
    * what they mean are seen in the if-else blocks below.
    */
   std::string arg1(argv[1]);
   std::string arg1value;
   if (arg1.find("--alter_mpi=", 0) < arg1.size()) {
      arg1value = arg1.c_str() + 12;
      --argc;
      for (int i = 1; i < argc; ++i) {
         argv[i] = argv[i + 1];
      }
   }

   // Start MPI first to we can generate a special communicator for SAMRAI.
   MPI_Init(&argc, &argv);

   MPI_Comm communicator = MPI_COMM_WORLD;
   {
      int rank;
      MPI_Comm_rank(communicator, &rank);
      if (arg1value == "drop1") {
         communicator = getSmallerMPI(communicator);
         if (rank == 0) {
            std::cout << "Dropped rank 1 from communicator." << std::endl;
         }
      } else if (arg1value == "rotate") {
         communicator = getRotatedMPI(communicator);
         if (rank == 0) {
            std::cout << "Rotated ranks in communicator." << std::endl;
         }
      } else if (arg1value.empty()) {
         if (rank == 0) {
            std::cout << "No change in communicator." << std::endl;
         }
      } else {
         if (rank == 0) {
            std::cout << "alter_mpi of " << arg1value << " unrecognized." << std::endl;
            MPI_Finalize();
            return 1;
         }
      }
   }

   if (communicator == MPI_COMM_NULL) {
      // This process has been excluded from the test.
      MPI_Finalize();
      return 0;
   }

   /*
    * Initialize SAMRAI.
    */
   SAMRAI_MPI::init(communicator);
   SAMRAIManager::initialize();
   SAMRAIManager::startup();
   tbox::SAMRAI_MPI samrai_mpi(SAMRAI_MPI::getSAMRAIWorld());
   tbox::SAMRAI_MPI world_mpi(MPI_COMM_WORLD);

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

   std::string case_name;
   if (argc > 2) {
      case_name = argv[2];
   }

   /*
    * Randomize the random number generator to avoid funny looking
    * first random number.
    */
   srand48(samrai_mpi.getRank() + 10 + samrai_mpi.getSize());

   int total_err_count = 0;

   {
      /*
       * Scope to force destruction of objects that would otherwise
       * leave allocated memory reported by the memory test.
       */

      /*
       * Create input database and parse all data in input file.
       */

      std::shared_ptr<InputDatabase> input_db(new InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Set up the timer manager.
       */
      if (input_db->isDatabase("TimerManager")) {
         TimerManager::createManager(input_db->getDatabase("TimerManager"));
      }
      /*
       * Create timers in the following orders so that they match
       * across processors.  That means creating some timers that are
       * not needed in the current scope.
       */
      TimerManager::getManager()->getTimer("apps::main::child_send");
      TimerManager::getManager()->getTimer("apps::main::child_recv");
      std::shared_ptr<tbox::Timer> t_child_wait = TimerManager::getManager()->getTimer(
            "apps::main::child_wait");
      TimerManager::getManager()->getTimer("apps::main::parent_send");
      TimerManager::getManager()->getTimer("apps::main::parent_recv");
      std::shared_ptr<tbox::Timer> t_parent_wait = TimerManager::getManager()->getTimer(
            "apps::main::parent_wait");
      TimerManager::getManager()->getTimer("apps::main::equiv_MPI");
      std::shared_ptr<tbox::Timer> t_processing = TimerManager::getManager()->getTimer(
            "apps::main::processing");

      /*
       * Retrieve "Main" section from input database.
       * The main database is used only in main().
       * The base_name variable is a base name for
       * all name strings in this program.
       */

      std::shared_ptr<Database> main_db = input_db->getDatabase("Main");

      std::string base_name = "unnamed";
      base_name = main_db->getStringWithDefault("base_name", base_name);

      /*
       * Modify basename for this particular run.
       * Add the number of processes and the case name.
       */
      if (!case_name.empty()) {
         base_name = base_name + '-' + case_name;
      }
      base_name = base_name + '-' + tbox::Utilities::intToString(
            samrai_mpi.getSize(),
            5);
      tbox::plog << "Added case name (" << case_name << ") and nprocs ("
                 << samrai_mpi.getSize() << ") to base name -> '"
                 << base_name << "'\n";

      if (!case_name.empty()) {
         tbox::plog << "Added case name (" << case_name << ") and nprocs ("
                    << samrai_mpi.getSize() << ") to base name -> '"
                    << base_name << "'\n";
      }

      /*
       * Start logging.
       */
      const std::string log_file_name = base_name + ".log";
      bool log_all_nodes = false;
      log_all_nodes = main_db->getBoolWithDefault("log_all_nodes",
            log_all_nodes);
      if (log_all_nodes) {
         PIO::logAllNodes(log_file_name);
      } else {
         PIO::logOnlyNodeZero(log_file_name);
      }
      plog << "This is SAMRAI process " << samrai_mpi.getRank() << " of " << samrai_mpi.getSize()
           << " and world process " << world_mpi.getRank() << " of " << world_mpi.getSize()
           << std::endl;

      calibrate_my_usleep();

      if (samrai_mpi.getCommunicator() != MPI_COMM_NULL) {
         /*
          * Run each test as it is pulled out of input_db.
          */

         int test_number = 0;

         CommGraphWriter comm_graph_writer;

         while (true) {

            std::string test_name("Test");
            test_name += Utilities::intToString(test_number, 2);

            std::shared_ptr<Database> test_db =
               input_db->getDatabaseWithDefault(test_name, std::shared_ptr<Database>());

            if (!test_db) {
               break;
            }

            const std::string nickname =
               test_db->getStringWithDefault("nickname", test_name);

            plog << "\n\n\nStarting test " << test_name << " (" << nickname << ")\n";

            CommonTestSwitches cts;
            getCommonTestSwitchesFromDatabase(cts, *test_db);

            const std::shared_ptr<RankTreeStrategy> rank_tree =
               getTreeForTesting(cts.tree_name, *test_db, samrai_mpi);

            // Write local part of tree to log.
            plog << "Tree " << cts.tree_name << ":\n"
                 << "  Root rank: " << rank_tree->getRootRank() << '\n'
                 << "  Child number: " << rank_tree->getChildNumber() << '\n'
                 << "  Generation number: " << rank_tree->getGenerationNumber() << '\n'
                 << "  Number of children: " << rank_tree->getNumberOfChildren() << '\n'
                 << "  " << rank_tree->getParentRank() << " <- " << rank_tree->getRank() << " ->";
            for (unsigned int i = 0; i < rank_tree->getNumberOfChildren(); ++i) {
               plog << ' ' << rank_tree->getChildRank(i);
            }
            plog << std::endl;

            int test_err_count = 0;

            if (cts.message_pattern == "UP") {
               test_err_count = testUp(*rank_tree, samrai_mpi, cts, *test_db);
            } else if (cts.message_pattern == "DOWN") {
               test_err_count = testDown(*rank_tree, samrai_mpi, cts, *test_db);
            } else if (cts.message_pattern == "UP_THEN_DOWN") {
               test_err_count = testUpThenDown(*rank_tree, samrai_mpi, cts, *test_db);
            } else if (cts.message_pattern == "TreeLB") {
               test_err_count = testTreeLB(*rank_tree, samrai_mpi, cts, *test_db);
            } else {
               TBOX_ERROR("Test message_pattern '" << cts.message_pattern << "' is not supported.");
            }

            total_err_count += test_err_count;

            plog << "\n";
            pout << "Completed " << test_name << " (" << nickname << ") with "
                 << test_err_count << " errs, total of " << total_err_count << "\n";

            if (test_err_count != 0) {
               perr << "Test FAILED.\n";
            }

            comm_graph_writer.addRecord(samrai_mpi, size_t(1 + rank_tree->getDegree()), size_t(1));

            for (unsigned int cn = 0; cn < rank_tree->getDegree(); ++cn) {
               comm_graph_writer.setEdgeInCurrentRecord(
                  size_t(cn),
                  "child_wait",
                  t_child_wait->getTotalWallclockTime() / cts.repetition,
                  CommGraphWriter::FROM,
                  rank_tree->getChildRank(cn));
            }

            comm_graph_writer.setEdgeInCurrentRecord(
               size_t(0 + rank_tree->getDegree()),
               "parent_wait",
               t_parent_wait->getTotalWallclockTime() / cts.repetition,
               CommGraphWriter::FROM,
               rank_tree->getParentRank());

            comm_graph_writer.setNodeValueInCurrentRecord(
               size_t(0),
               "processing_time",
               t_processing->getTotalWallclockTime() / cts.repetition);

            comm_graph_writer.writeGraphToTextStream(
               comm_graph_writer.getNumberOfRecords() - 1, plog);

            // Output timer results then reset for next test.
            tbox::TimerManager::getManager()->print(tbox::plog);
            tbox::TimerManager::getManager()->resetAllTimers();

            ++test_number;

         }

         if (samrai_mpi != SAMRAI_MPI::getSAMRAIWorld()) {
            samrai_mpi.freeCommunicator();
         }

         tbox::plog << "TreeCommunicationBenchmark completed " << test_number << " tests."
                    << std::endl;

      }

   }

   /*
    * Print input database again to fully show usage.
    */
   plog << "Input database after running..." << std::endl;
   tbox::InputManager::getManager()->getInputDatabase()->printClassData(plog);

   if (total_err_count == 0) {
      tbox::pout << "\nPASSED:  TreeCommunicationBenchmark" << std::endl;
   }

   /*
    * Exit properly by shutting down services in correct order.
    */
   tbox::plog << "\nShutting down..." << std::endl;

   /*
    * Shut down.
    */
   SAMRAIManager::shutdown();
   SAMRAIManager::finalize();
   SAMRAI_MPI::finalize();
   MPI_Finalize();

   return fail_count;

#endif
}

/*
 *************************************************************************
 * Get common test switches from the test database.
 *************************************************************************
 */
void getCommonTestSwitchesFromDatabase(
   CommonTestSwitches& cts,
   Database& test_db)
{
   cts.tree_name = test_db.getString("tree_name");
   cts.message_pattern = test_db.getString("message_pattern");
   cts.msg_length = test_db.getIntegerWithDefault("msg_length", 1);
   cts.first_data_length = test_db.getIntegerWithDefault("first_data_length", 1);
   cts.verify_data = test_db.getBoolWithDefault("verify_data", false);
   cts.repetition = test_db.getIntegerWithDefault("repetition", 1);
   cts.barrier_after_each_repetition = test_db.getBoolWithDefault("barrier_after_each_repetition",
         false);
   if (test_db.isInteger("mpi_tags")) {
      test_db.getIntegerArray("mpi_tags", cts.mpi_tags, 2);
   }
   if (test_db.isInteger("processing_cost")) {
      test_db.getIntegerArray("processing_cost", cts.processing_cost, 2);
   }
   cts.randomize_processing_cost = test_db.getBoolWithDefault("randomize_processing_cost", false);
}

/*
 ****************************************************************************
 * Test sending up the tree.
 ****************************************************************************
 */
int testUp(
   const RankTreeStrategy& rank_tree,
   const SAMRAI_MPI& mpi,
   const CommonTestSwitches& cts,
   Database& test_db)
{
   plog << "Test database:\n";
   test_db.printClassData(plog);

   int err_count = 0;

   AsyncCommStage child_stage;
   AsyncCommPeer<int>* child_comms;
   AsyncCommStage parent_stage;
   AsyncCommPeer<int>* parent_comm;
   setupAsyncComms(child_stage, child_comms, parent_stage, parent_comm,
      mpi, rank_tree, cts);

   std::vector<int> msg(cts.msg_length, 1);
   setMessageData(msg, mpi.getRank());

   std::shared_ptr<Timer> repetitions_timer = TimerManager::getManager()->getTimer(
         "apps::main::repetitions");
   std::shared_ptr<Timer> processing_timer = TimerManager::getManager()->getTimer(
         "apps::main::processing");
   std::shared_ptr<Timer> verify_timer = TimerManager::getManager()->getTimer(
         "apps::main::verify");
   std::shared_ptr<Timer> equiv_mpi_timer = TimerManager::getManager()->getTimer(
         "apps::main::equiv_MPI");

   mpi.Barrier();
   repetitions_timer->start();
   for (int r = 0; r < cts.repetition; ++r) {

      for (unsigned int ic = 0; ic < rank_tree.getNumberOfChildren(); ++ic) {
         child_comms[ic].beginRecv();
      }

      child_stage.advanceAll();
      child_stage.clearCompletionQueue();

      if (cts.verify_data) {
         verify_timer->start();
         for (unsigned int ic = 0; ic < rank_tree.getNumberOfChildren(); ++ic) {
            err_count += verifyReceivedData(child_comms[ic]);
         }
         verify_timer->stop();
      }

      processing_timer->start();
      simulateDataProcessing(cts);
      processing_timer->stop();

      if (rank_tree.getParentRank() != RankTreeStrategy::getInvalidRank()) {
         if (!parent_comm->isDone()) {
            parent_comm->completeCurrentOperation(); // From previous repetition.
         }
         parent_comm->beginSend(&msg[0], cts.msg_length);
      }

      if (cts.barrier_after_each_repetition) {
         mpi.Barrier();
      }
   }
   repetitions_timer->stop();

   parent_stage.advanceAll(); // Make sure all is completed before destroying.
   destroyAsyncComms(child_comms, parent_comm);

   /*
    * Run equivalent MPI operations for comparison.
    */

   std::vector<int> msgr(cts.msg_length, 1);
   equiv_mpi_timer->start();
   for (int r = 0; r < cts.repetition; ++r) {
      mpi.Reduce(&msg[0], &msgr[0], cts.msg_length, MPI_INT, MPI_SUM,
         rank_tree.getRootRank());
      if (cts.barrier_after_each_repetition) {
         mpi.Barrier();
      }
   }
   equiv_mpi_timer->stop();

   return err_count;
}

/*
 ****************************************************************************
 * Test sending down the tree.
 ****************************************************************************
 */
int testDown(
   const RankTreeStrategy& rank_tree,
   const SAMRAI_MPI& mpi,
   const CommonTestSwitches& cts,
   Database& test_db)
{
   plog << "Test database:\n";
   test_db.printClassData(plog);

   int err_count = 0;

   AsyncCommStage child_stage;
   AsyncCommPeer<int>* child_comms;
   AsyncCommStage parent_stage;
   AsyncCommPeer<int>* parent_comm;
   setupAsyncComms(child_stage, child_comms, parent_stage, parent_comm,
      mpi, rank_tree, cts);

   std::vector<int> msg(cts.msg_length, 1);
   setMessageData(msg, mpi.getRank());

   std::shared_ptr<Timer> repetitions_timer = TimerManager::getManager()->getTimer(
         "apps::main::repetitions");
   std::shared_ptr<Timer> processing_timer = TimerManager::getManager()->getTimer(
         "apps::main::processing");
   std::shared_ptr<Timer> verify_timer = TimerManager::getManager()->getTimer(
         "apps::main::verify");
   std::shared_ptr<Timer> equiv_mpi_timer = TimerManager::getManager()->getTimer(
         "apps::main::equiv_MPI");

   mpi.Barrier();
   repetitions_timer->start();
   for (int r = 0; r < cts.repetition; ++r) {

      if (rank_tree.getParentRank() != RankTreeStrategy::getInvalidRank()) {
         parent_comm->beginRecv();
         parent_comm->completeCurrentOperation();

         if (cts.verify_data) {
            verify_timer->start();
            err_count += verifyReceivedData(*parent_comm);
            verify_timer->stop();
         }
      }

      processing_timer->start();
      simulateDataProcessing(cts);
      processing_timer->stop();

      for (unsigned int ic = 0; ic < rank_tree.getNumberOfChildren(); ++ic) {
         if (!child_comms[ic].isDone()) {
            child_comms[ic].completeCurrentOperation(); // From previous repetition.
         }
         child_comms[ic].beginSend(&msg[0], cts.msg_length);
      }

      child_stage.advanceAll();
      child_stage.clearCompletionQueue();

      if (cts.barrier_after_each_repetition) {
         mpi.Barrier();
      }
   }
   repetitions_timer->stop();

   child_stage.advanceAll(); // Make sure all is completed before destroying.
   destroyAsyncComms(child_comms, parent_comm);

   /*
    * Run equivalent MPI operations for comparison.
    */

   equiv_mpi_timer->start();
   for (int r = 0; r < cts.repetition; ++r) {
      mpi.Bcast(&msg[0], cts.msg_length, MPI_INT, rank_tree.getRootRank());
      if (cts.barrier_after_each_repetition) {
         mpi.Barrier();
      }
   }
   equiv_mpi_timer->stop();

   return err_count;
}

/*
 ****************************************************************************
 * Test sending up then down the tree.
 ****************************************************************************
 */
int testUpThenDown(
   const RankTreeStrategy& rank_tree,
   const SAMRAI_MPI& mpi,
   const CommonTestSwitches& cts,
   Database& test_db)
{
   plog << "Test database:\n";
   test_db.printClassData(plog);

   int err_count = 0;

   AsyncCommStage child_stage;
   AsyncCommPeer<int>* child_comms;
   AsyncCommStage parent_stage;
   AsyncCommPeer<int>* parent_comm;
   setupAsyncComms(child_stage, child_comms, parent_stage, parent_comm,
      mpi, rank_tree, cts);

   std::vector<int> msg(cts.msg_length, 1);
   setMessageData(msg, mpi.getRank());

   std::shared_ptr<Timer> repetitions_timer = TimerManager::getManager()->getTimer(
         "apps::main::repetitions");
   std::shared_ptr<Timer> processing_timer = TimerManager::getManager()->getTimer(
         "apps::main::processing");
   std::shared_ptr<Timer> verify_timer = TimerManager::getManager()->getTimer(
         "apps::main::verify");
   std::shared_ptr<Timer> equiv_mpi_timer = TimerManager::getManager()->getTimer(
         "apps::main::equiv_MPI");

   mpi.Barrier();
   repetitions_timer->start();
   for (int r = 0; r < cts.repetition; ++r) {

      for (unsigned int ic = 0; ic < rank_tree.getNumberOfChildren(); ++ic) {
         child_comms[ic].beginRecv();
      }

      child_stage.advanceAll();
      child_stage.clearCompletionQueue();

      if (cts.verify_data) {
         verify_timer->start();
         for (unsigned int ic = 0; ic < rank_tree.getNumberOfChildren(); ++ic) {
            err_count += verifyReceivedData(child_comms[ic]);
         }
         verify_timer->stop();
      }

      processing_timer->start();
      simulateDataProcessing(cts);
      processing_timer->stop();

      if (rank_tree.getParentRank() != RankTreeStrategy::getInvalidRank()) {
         parent_comm->beginSend(&msg[0], cts.msg_length);
         parent_comm->completeCurrentOperation();

         parent_comm->beginRecv();
         parent_comm->completeCurrentOperation();

         if (cts.verify_data) {
            verify_timer->start();
            err_count += verifyReceivedData(*parent_comm);
            verify_timer->stop();
         }

         processing_timer->start();
         simulateDataProcessing(cts);
         processing_timer->stop();
      }

      for (unsigned int ic = 0; ic < rank_tree.getNumberOfChildren(); ++ic) {
         if (!child_comms[ic].isDone()) {
            TBOX_ERROR("Test code error: Should never be here.");
         }
         child_comms[ic].beginSend(&msg[0], cts.msg_length);
      }

      child_stage.advanceAll();
      child_stage.clearCompletionQueue();

      if (cts.barrier_after_each_repetition) {
         mpi.Barrier();
      }
   }
   repetitions_timer->stop();

   child_stage.advanceAll(); // Make sure all is completed before destroying.
   destroyAsyncComms(child_comms, parent_comm);

   /*
    * Run equivalent MPI operations for comparison.
    */

   std::vector<int> msgr(cts.msg_length, 1);
   mpi.Barrier();
   equiv_mpi_timer->start();
   for (int r = 0; r < cts.repetition; ++r) {
      mpi.Allreduce(&msg[0], &msgr[0], cts.msg_length, MPI_INT, MPI_SUM);
      if (cts.barrier_after_each_repetition) {
         mpi.Barrier();
      }
   }
   equiv_mpi_timer->stop();

   return err_count;
}

/*
 ****************************************************************************
 * Simulation the communication of the TreeLoadBalancer.
 * Test sending up then down the tree.  Some ranks expect a down
 * message dependent on the parent, some expect one dependent
 * on the grandparent and some expect no down message at all.
 * See the extra input parameter down_message_dependency.
 ****************************************************************************
 */
int testTreeLB(
   const RankTreeStrategy& rank_tree,
   const SAMRAI_MPI& mpi,
   const CommonTestSwitches& cts,
   Database& test_db)
{

   std::vector<int> down_message_dependency =
      test_db.getIntegerVector("down_message_dependency");
   const int dl = static_cast<int>(down_message_dependency.size());

#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < dl; ++i) {
      TBOX_ASSERT(down_message_dependency[i] >= 0 &&
         down_message_dependency[i] < 3);
   }
#endif

   plog << "Test database:\n";
   test_db.printClassData(plog);

   int err_count = 0;

   AsyncCommStage child_stage;
   AsyncCommPeer<int>* child_comms;
   AsyncCommStage parent_stage;
   AsyncCommPeer<int>* parent_comm;
   setupAsyncComms(child_stage, child_comms, parent_stage, parent_comm,
      mpi, rank_tree, cts);

   std::vector<int> msg(cts.msg_length, 1);
   setMessageData(msg, mpi.getRank());

   std::shared_ptr<Timer> repetitions_timer = TimerManager::getManager()->getTimer(
         "apps::main::repetitions");
   std::shared_ptr<Timer> processing_timer = TimerManager::getManager()->getTimer(
         "apps::main::processing");
   std::shared_ptr<Timer> verify_timer = TimerManager::getManager()->getTimer(
         "apps::main::verify");

   mpi.Barrier();
   repetitions_timer->start();
   for (int r = 0; r < cts.repetition; ++r) {

      for (unsigned int ic = 0; ic < rank_tree.getNumberOfChildren(); ++ic) {
         child_comms[ic].beginRecv();
      }

      child_stage.advanceAll();
      child_stage.clearCompletionQueue();

      if (cts.verify_data) {
         verify_timer->start();
         for (unsigned int ic = 0; ic < rank_tree.getNumberOfChildren(); ++ic) {
            err_count += verifyReceivedData(child_comms[ic]);
         }
         verify_timer->stop();
      }

      processing_timer->start();
      simulateDataProcessing(cts);
      processing_timer->stop();

      if (rank_tree.getParentRank() != RankTreeStrategy::getInvalidRank()) {
         parent_comm->beginSend(&msg[0], cts.msg_length);
      }

      // Send down-messages to children that are independent of their grandparents.
      for (unsigned int ic = 0; ic < rank_tree.getNumberOfChildren(); ++ic) {
         const int child_rank = rank_tree.getChildRank(ic);
         if (down_message_dependency[child_rank % dl] == 1) {
            if (!child_comms[ic].isDone()) {
               TBOX_ERROR("Test code error: Should never be here.");
            }
            child_comms[ic].beginSend(&msg[0], cts.msg_length);
         }
      }

      if (rank_tree.getParentRank() != RankTreeStrategy::getInvalidRank()) {
         parent_comm->completeCurrentOperation();

         if (down_message_dependency[rank_tree.getRank() % dl] != 0) {
            // This process expects a message from its parent.
            parent_comm->beginRecv();
            parent_comm->completeCurrentOperation();

            if (cts.verify_data) {
               verify_timer->start();
               err_count += verifyReceivedData(*parent_comm);
               verify_timer->stop();
            }

            processing_timer->start();
            simulateDataProcessing(cts);
            processing_timer->stop();
         }
      }

      // Send down-messages to children that depend on their grandparents.
      for (unsigned int ic = 0; ic < rank_tree.getNumberOfChildren(); ++ic) {
         const int child_rank = rank_tree.getChildRank(ic);
         if (down_message_dependency[child_rank % dl] == 2) {
            if (!child_comms[ic].isDone()) {
               TBOX_ERROR("Test code error: Should never be here.");
            }
            child_comms[ic].beginSend(&msg[0], cts.msg_length);
         }
      }

      child_stage.advanceAll();
      child_stage.clearCompletionQueue();

      if (cts.barrier_after_each_repetition) {
         mpi.Barrier();
      }
   }
   repetitions_timer->stop();

   child_stage.advanceAll(); // Make sure all is completed before destroying.
   destroyAsyncComms(child_comms, parent_comm);

   return err_count;
}

/*
 ****************************************************************************
 * Simulate the data_processing.
 ****************************************************************************
 */
void simulateDataProcessing(const CommonTestSwitches& cts)
{
   double randf = cts.randomize_processing_cost ? tbox::MathUtilities<double>::Rand(0.0, 1.0) : 1.0;
   unsigned long microsecs =
      static_cast<unsigned long>(randf
                                 * (cts.processing_cost[0] + cts.msg_length
                                    * cts.processing_cost[1]));
   my_usleep(microsecs);
}

/*
 ****************************************************************************
 * Sleep for some count of microseconds, a substitute for usleep.
 ****************************************************************************
 */
void my_usleep(size_t num_usec)
{
   size_t num_calls = static_cast<size_t>(double(num_usec) * my_usleep_calls_per_usec);
   for (size_t i = 0; i < num_calls; ++i) {
      dummy();
   }
}

/*
 ****************************************************************************
 * Calibration for my_usleep, a substitute for usleep.
 ****************************************************************************
 */
void calibrate_my_usleep()
{
   tbox::plog << "Calibrating my_usleep.\n";

   const size_t num_samples = 1000000;

   clock_t user_start, system_start;
   double wallclock_start;
   tbox::Clock::timestamp(user_start, system_start, wallclock_start);

   for (size_t i = 0; i < num_samples; ++i) {
      dummy();
   }

   clock_t user_stop, system_stop;
   double wallclock_stop;
   tbox::Clock::timestamp(user_stop, system_stop, wallclock_stop);

   double wall_time = wallclock_stop - wallclock_start;
   my_usleep_calls_per_usec = 1e-6 * num_samples / wall_time;

   tbox::plog << "Calibration completed.  my_usleep_calls_per_usec = "
              << my_usleep_calls_per_usec;

}

/*
 ****************************************************************************
 * Initialize data.
 * Fill array starting with value of rank, then increasing.
 ****************************************************************************
 */
void setMessageData(
   std::vector<int>& msg,
   int rank)
{
   for (size_t i = 0; i < msg.size(); ++i) {
      msg[i] = rank + static_cast<int>(i);
   }
}

/*
 ****************************************************************************
 * Verify that data matches what is set by setMessageData.
 ****************************************************************************
 */
int verifyReceivedData(
   const AsyncCommPeer<int>& peer_comm)
{
   TBOX_ASSERT(peer_comm.isDone());
   int err_count = 0;
   const int* msg = peer_comm.getRecvData();
   const int msg_length = peer_comm.getRecvSize();
   const int rank = peer_comm.getPeerRank();
   for (int i = 0; i < msg_length; ++i) {
      if (msg[i] != rank + i) {
         ++err_count;
      }
   }
   return err_count;
}

/*
 ****************************************************************************
 * Get the RankTreeStrategy implementation for the test Database.
 ****************************************************************************
 */
std::shared_ptr<RankTreeStrategy> getTreeForTesting(
   const std::string& tree_name,
   Database& test_db,
   const SAMRAI_MPI& mpi)
{
   std::shared_ptr<tbox::RankTreeStrategy> rank_tree;

   if (tree_name == "BalancedDepthFirstTree") {

      BalancedDepthFirstTree * bdfs(new BalancedDepthFirstTree());

      if (test_db.isDatabase("BalancedDepthFirstTree")) {
         std::shared_ptr<tbox::Database> tmp_db = test_db.getDatabase("BalancedDepthFirstTree");
         bool do_left_leaf_switch = tmp_db->getBoolWithDefault("do_left_leaf_switch", true);
         bdfs->setLeftLeafSwitching(do_left_leaf_switch);
      }

      bdfs->setupTree(RankGroup(mpi), mpi.getRank());
      rank_tree.reset(bdfs);

   } else if (tree_name == "CenteredRankTree") {

      CenteredRankTree * crt(new tbox::CenteredRankTree());

      if (test_db.isDatabase("CenteredRankTree")) {
         std::shared_ptr<tbox::Database> tmp_db = test_db.getDatabase("CenteredRankTree");
         bool make_first_rank_the_root = tmp_db->getBoolWithDefault("make_first_rank_the_root",
               true);
         crt->makeFirstRankTheRoot(make_first_rank_the_root);
      }

      crt->setupTree(RankGroup(mpi), mpi.getRank());
      rank_tree.reset(crt);

   } else if (tree_name == "BreadthFirstRankTree") {

      BreadthFirstRankTree * dft(new tbox::BreadthFirstRankTree());

      if (test_db.isDatabase("BreadthFirstRankTree")) {
         std::shared_ptr<tbox::Database> tmp_db = test_db.getDatabase("BreadthFirstRankTree");
         const int tree_degree = tmp_db->getIntegerWithDefault("tree_degree", true);
         dft->setTreeDegree(static_cast<unsigned short>(tree_degree));
      }

      dft->setupTree(RankGroup(mpi), mpi.getRank());
      rank_tree.reset(dft);

   } else {
      TBOX_ERROR("Unrecognized RankTreeStrategy " << tree_name);
   }

   return rank_tree;
}

/*
 ****************************************************************************
 *
 ****************************************************************************
 */
void setupAsyncComms(
   AsyncCommStage& child_stage,
   AsyncCommPeer<int> *& child_comms,
   AsyncCommStage& parent_stage,
   AsyncCommPeer<int> *& parent_comm,
   const SAMRAI_MPI& mpi,
   const RankTreeStrategy& rank_tree,
   const CommonTestSwitches& cts)
{
   child_comms = parent_comm = 0;

   const int num_children = rank_tree.getNumberOfChildren();

   if (num_children > 0) {

      child_comms = new tbox::AsyncCommPeer<int>[num_children];

      for (int child_num = 0; child_num < num_children; ++child_num) {

         const int child_rank = rank_tree.getChildRank(child_num);

         child_comms[child_num].initialize(&child_stage);
         child_comms[child_num].setPeerRank(child_rank);
         child_comms[child_num].setMPI(mpi);
         child_comms[child_num].setMPITag(cts.mpi_tags[0], cts.mpi_tags[1]);
         child_comms[child_num].limitFirstDataLength(cts.first_data_length);
         child_comms[child_num].setSendTimer(
            TimerManager::getManager()->getTimer("apps::main::child_send"));
         child_comms[child_num].setRecvTimer(
            TimerManager::getManager()->getTimer("apps::main::child_recv"));
         child_comms[child_num].setWaitTimer(
            TimerManager::getManager()->getTimer("apps::main::child_wait"));
      }
   }

   if (rank_tree.getParentRank() != tbox::RankTreeStrategy::getInvalidRank()) {

      const int parent_rank = rank_tree.getParentRank();

      parent_comm = new tbox::AsyncCommPeer<int>;
      parent_comm->initialize(&parent_stage);
      parent_comm->setPeerRank(parent_rank);
      parent_comm->setMPI(mpi);
      parent_comm->setMPITag(cts.mpi_tags[0], cts.mpi_tags[1]);
      parent_comm->limitFirstDataLength(cts.first_data_length);
      parent_comm->setSendTimer(
         TimerManager::getManager()->getTimer("apps::main::parent_send"));
      parent_comm->setRecvTimer(
         TimerManager::getManager()->getTimer("apps::main::parent_recv"));
      parent_comm->setWaitTimer(
         TimerManager::getManager()->getTimer("apps::main::parent_wait"));

   }

   parent_stage.setCommunicationWaitTimer(TimerManager::getManager()->getTimer(
         "apps::main::parent_wait"));
   child_stage.setCommunicationWaitTimer(TimerManager::getManager()->getTimer(
         "apps::main::child_wait"));
}

/*
 *************************************************************************
 *************************************************************************
 */
void destroyAsyncComms(
   AsyncCommPeer<int> *& child_comms,
   AsyncCommPeer<int> *& parent_comm)
{
   if (child_comms != 0) {
      delete[] child_comms;
   }
   if (parent_comm != 0) {
      delete parent_comm;
   }
   child_comms = parent_comm = 0;
}

/*
 *************************************************************************
 * Create a new communicator with ranks rotated downward by 1.
 *************************************************************************
 */
SAMRAI_MPI::Comm getRotatedMPI(const SAMRAI_MPI::Comm& old_comm)
{
#ifdef HAVE_MPI
   int old_rank, old_size;
   MPI_Comm_rank(old_comm, &old_rank);
   MPI_Comm_size(old_comm, &old_size);

   if (old_size == 1) {
      return old_comm;
   }

   MPI_Group old_group;
   MPI_Comm_group(old_comm, &old_group);

   MPI_Group new_group;
   int ranges[2][3];
   ranges[0][0] = 1;
   ranges[0][1] = old_size - 1;
   ranges[0][2] = 1;
   ranges[1][0] = 0;
   ranges[1][1] = 0;
   ranges[1][2] = 1;
   MPI_Group_range_incl(old_group, 2, ranges, &new_group);

   SAMRAI_MPI::Comm new_comm;
   MPI_Comm_create(old_comm, new_group, &new_comm);
   return new_comm;

#else
   NULL_USE(old_comm);
   return SAMRAI_MPI::Comm();

#endif
}

/*
 *************************************************************************
 * Create a new communicator with process 1 removed.
 *************************************************************************
 */
SAMRAI_MPI::Comm getSmallerMPI(const SAMRAI_MPI::Comm& old_comm)
{
#ifdef HAVE_MPI
   int old_rank, old_size;
   MPI_Comm_rank(old_comm, &old_rank);
   MPI_Comm_size(old_comm, &old_size);

   if (old_size == 1) {
      return old_comm;
   }

   MPI_Group old_group;
   MPI_Comm_group(old_comm, &old_group);

   MPI_Group new_group;
   int ranges[2][3];
   ranges[0][0] = 0;
   ranges[0][1] = 0;
   ranges[0][2] = 1;
   if (old_size == 2) {
      MPI_Group_range_incl(old_group, 1, ranges, &new_group);
   } else {
      ranges[1][0] = 2;
      ranges[1][1] = old_size - 1;
      ranges[1][2] = 1;
      MPI_Group_range_incl(old_group, 2, ranges, &new_group);
   }

   SAMRAI_MPI::Comm new_comm;
   MPI_Comm_create(old_comm, new_group, &new_comm);
   return old_rank == 1 ? MPI_COMM_NULL : new_comm;

#else
   NULL_USE(old_comm);
   return SAMRAI_MPI::Comm();

#endif
}
