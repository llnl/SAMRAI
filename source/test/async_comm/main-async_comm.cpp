/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test program for asynchronous communication classes
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/AsyncCommGroup.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <iomanip>

using namespace SAMRAI;
using namespace tbox;

/*
 ************************************************************************
 *
 * This program tests the asynchronous communication classes:
 * AsyncCommGroup
 * AsyncCommStage
 *
 * 1. Group the processors.  See code for heuristic rule for
 * defining groups.
 *
 * 2. Perform asynchronous communication within each group.
 *
 * 3. Check results.
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

   const int rank = mpi.getRank();
   int fail_count = 0;

   {

      /*
       * Process command line arguments.  For each run, the input
       * filename must be specified.  Usage is:
       *
       * executable <input file name>
       *
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
       * Make sure all processes are alive and well before running test.
       */
      SAMRAI_MPI::getSAMRAIWorld().Barrier();
      plog << "Process " << std::setw(5) << rank << " is ready." << std::endl;

      /*
       * Created a separate communicator for testing,
       * to avoid possible interference with other communications
       * by SAMRAI library.
       */
      tbox::SAMRAI_MPI isolated_mpi(MPI_COMM_NULL);
      isolated_mpi.dupCommunicator(SAMRAI_MPI::getSAMRAIWorld());
      plog << "Process " << std::setw(5) << rank
           << " duplicated Communicator." << std::endl;

      /*
       * Create input database and parse all data in input file.
       */

      std::shared_ptr<InputDatabase> input_db(new InputDatabase("input_db"));
      InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Set up the timer manager.
       */
      if (input_db->isDatabase("TimerManager")) {
         TimerManager::createManager(input_db->getDatabase("TimerManager"));
      }

      /*
       * Retrieve "Main" section from input database.
       * The main database is used only in main().
       * The base_name variable is a base name for
       * all name strings in this program.
       */

      std::shared_ptr<Database> main_db(input_db->getDatabase("Main"));
      std::string base_name = "unnamed";
      base_name = main_db->getStringWithDefault("base_name", base_name);

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

      plog << "********************* Note! *********************\n"
           << "* The asychronous communication classes are meant for\n"
           << "* large processor counts.\n"
           << "*\n"
           << "* For this test to be significant, you should run it on\n"
           << "* lots of processors.  I recommend running on a\n"
           << "* 'massively parallel' processor count (however you\n"
           << "* would like to define massively parallel).\n"
           << "*\n"
           << "* This program should not be used for performance\n"
           << "* testing because performance may be context-sensitive\n"
           << "* and these groups are rather contrived.";

      plog << "\n\n\n";

      const int sync_bcast_cycles =
         main_db->getIntegerWithDefault("sync_bcast_cycles", 1);
      const int sync_sumreduce_cycles =
         main_db->getIntegerWithDefault("sync_sumreduce_cycles", 1);
      const int asyncany_bcast_cycles =
         main_db->getIntegerWithDefault("asyncany_bcast_cycles", 1);
      const int asyncany_sumreduce_cycles =
         main_db->getIntegerWithDefault("asyncany_sumreduce_cycles", 1);
      const int asyncsome_bcast_cycles =
         main_db->getIntegerWithDefault("asyncsome_bcast_cycles", 1);
      const int asyncsome_sumreduce_cycles =
         main_db->getIntegerWithDefault("asyncsome_sumreduce_cycles", 1);

      int sync_bcast_count = 0;
      int sync_sumreduce_count = 0;
      int asyncany_bcast_count = 0;
      int asyncany_sumreduce_count = 0;
      int asyncsome_bcast_count = 0;
      int asyncsome_sumreduce_count = 0;

      const int def_num_groups = (mpi.getSize() + 1) / 2;
      plog << "Default num groups: " << def_num_groups << std::endl;
      const int num_groups =
         main_db->getIntegerWithDefault("num_groups", def_num_groups);
      plog << "Num groups: " << num_groups << std::endl;
      const int num_children =
         main_db->getIntegerWithDefault("num_children", 2);

      if (mpi.getRank() == 0) {
         plog << "Running num_groups = " << num_groups << std::endl;
         plog << "Running num_children = " << num_children << std::endl;
      }

      int pass_count = 0;

      int gi; // Group index.
      int ai; // Active group index.

      int count = 0;
      while ((sync_bcast_count < sync_bcast_cycles) ||
             (sync_sumreduce_count < sync_sumreduce_cycles) ||
             (asyncany_bcast_count < asyncany_bcast_cycles) ||
             (asyncany_sumreduce_count < asyncany_sumreduce_cycles) ||
             (asyncsome_bcast_count < asyncsome_bcast_cycles) ||
             (asyncsome_sumreduce_count < asyncsome_sumreduce_cycles)) {

         if (mpi.getRank() == 0) {
            plog << " Starting cycle number " << count << std::endl;
         }

         plog << "\n\n\n***************** Beginning Cycle Number "
              << count << " *******************\n\n";

         std::vector<std::vector<int> > group_ids(num_groups);
         std::vector<int> owners(num_groups);
         std::vector<int> active_flags(num_groups);

         std::vector<int> active_groups(num_groups);
         int num_active_groups = 0;

         /*
          * Define groups.
          * Group n includes all process whose rank divisible by n+1 -- with
          * a slight variation.  With each testing cycle, the "rank" is
          * increased by one.
          * Set owner of each group to roughly the processor in the middle
          * of the group.
          */
         for (int n = 0; n < num_groups; ++n) {

            int gsize = (mpi.getSize() + n) / (n + 1);
            group_ids[n].resize(gsize);
            active_flags[n] = false;
            for (int i = 0; i < gsize; ++i) {
               group_ids[n][i] = i * (n + 1);
               group_ids[n][i] =
                  (group_ids[n][i] + count) % mpi.getSize();
               if (group_ids[n][i] == rank) {
                  active_groups[num_active_groups++] = n;
                  active_flags[n] = true;
               }
            }

            owners[n] = group_ids[n][gsize / 2];

         }
         active_groups.resize(num_active_groups);

         /*
          * Write out group data.
          */
         plog << "Group definitions (" << num_groups << " groups):\n";
         plog << "(Groups with '*' contains the local process, "
              << rank << ".)\n";
         plog << "(Groups with '**' is owned by the local process\n\n";
         plog << " ID  size owner members...\n";
         for (int n = 0; n < num_groups; ++n) {
            plog << std::setw(3) << n
                 << std::setw(5) << group_ids[n].size()
                 << std::setw(4) << owners[n]
                 << (active_flags[n] ? '*' : ' ')
                 << (owners[n] == rank ? '*' : ' ') << ':';
            for (int i = 0; i < static_cast<int>(group_ids[n].size()); ++i) {
               plog << "  " << group_ids[n][i];
            }
            plog << '\n';
         }
         plog << '\n';

         plog << "Active groups (" << num_active_groups << " groups):";
         for (ai = 0; ai < num_active_groups; ++ai) {
            plog << "  " << active_groups[ai];
         }
         plog << "\n\n";

         /*
          * Initialize data for sum-reduce tests.
          * Compute the correct sum for comparison.
          */

         std::vector<int> sum(num_active_groups);
         std::vector<int> correct_sum(num_active_groups);

         for (ai = 0; ai < num_active_groups; ++ai) {
            sum[ai] = 1 + rank;
            correct_sum[ai] = 0;
            std::vector<int>& g = group_ids[active_groups[ai]];
            for (int j = 0; j < static_cast<int>(g.size()); ++j) {
               correct_sum[ai] += 1 + g[j];
            }
         }

         /*
          * Initialize data for broadcast test.
          * Broadcast data is 1001 + the group index.
          */
         std::vector<int> bcdata(num_active_groups);
         std::vector<int> correct_bcdata(num_active_groups);
         for (ai = 0; ai < num_active_groups; ++ai) {
            gi = active_groups[ai];
            bcdata[ai] = rank == owners[gi] ? 1001 + gi : -1;
            correct_bcdata[ai] = 1001 + gi;
         }

         /*
          * Create the communication stage and groups.
          * Each group uses its group index as the MPI tag.
          */
         AsyncCommStage comm_stage;
         AsyncCommGroup* comm_groups = new AsyncCommGroup[num_active_groups];
         for (ai = 0; ai < num_active_groups; ++ai) {
            gi = active_groups[ai];
            plog << "Initializing group " << gi << "\n";
            comm_groups[ai].initialize(num_children,
               &comm_stage);
            comm_groups[ai].setGroupAndRootRank(isolated_mpi,
               &group_ids[gi][0],
               static_cast<int>(group_ids[gi].size()),
               owners[gi]);
            comm_groups[ai].setMPITag(1000000 * count + gi);
            comm_groups[ai].setUseBlockingSendToParent(false);
            comm_groups[ai].setUseBlockingSendToChildren(false);
         }

         if (sync_bcast_count < sync_bcast_cycles) {
            /*
             * For the synchronous (groupwise) broadcast test,
             * each group broadcasts the its group id.
             */
            TBOX_ASSERT(!comm_stage.hasCompletedMembers());
            plog << "\n\n\n*********** Synchronous Broadcast "
                 << sync_bcast_count << " ************\n";
            for (ai = 0; ai < num_active_groups; ++ai)
               if (rank != owners[active_groups[ai]]) bcdata[ai] = -1;
            plog << "Job Group Result Correct  Note\n";
            for (ai = 0; ai < num_active_groups; ++ai) {
               AsyncCommGroup& comm_group = comm_groups[ai];
               comm_group.beginBcast(&bcdata[ai], 1);
               comm_group.completeCurrentOperation();
               gi = active_groups[ai];
               plog << std::setw(3) << ai
                    << std::setw(5) << gi
                    << std::setw(8) << bcdata[ai]
                    << std::setw(8) << correct_bcdata[ai]
               ;
               plog << "  Bcast difference = "
                    << bcdata[ai] - correct_bcdata[ai];
               if (bcdata[ai] != correct_bcdata[ai]) {
                  plog << "  Error!";
                  tbox::pout << "Error in bcast result for group "
                             << gi << std::endl;
                  ++fail_count;
               } else ++pass_count;
               plog << std::endl;
            }
            for (ai = 0; ai < num_active_groups; ++ai) {
               TBOX_ASSERT(comm_groups[ai].isDone());
            }
            TBOX_ASSERT(!comm_stage.hasPendingRequests());
            ++sync_bcast_count;
         }

         if (sync_sumreduce_count < sync_sumreduce_cycles) {
            TBOX_ASSERT(!comm_stage.hasCompletedMembers());
            /*
             * For the sum advanceSome reduce test,
             * each group sums up the ranks of its members, plus 1.
             */
            plog << "\n\n\n*********** Synchronous Sum Reduce "
                 << sync_sumreduce_count << " ************\n";
            for (ai = 0; ai < num_active_groups; ++ai) sum[ai] = 1 + rank;
            plog << "Job Group Result Correct  Note\n";
            for (ai = 0; ai < num_active_groups; ++ai) {
               AsyncCommGroup& comm_group = comm_groups[ai];
               comm_group.beginSumReduce(&sum[ai], 1);
               comm_group.completeCurrentOperation();
               TBOX_ASSERT(comm_group.isDone());
               gi = active_groups[ai];
               plog << std::setw(3) << ai
                    << std::setw(5) << gi
                    << std::setw(8) << sum[ai]
                    << std::setw(8) << correct_sum[ai]
               ;
               if (rank == owners[gi]) {
                  plog << "  Sum reduce difference = "
                       << sum[ai] - correct_sum[ai];
                  if (sum[ai] != correct_sum[ai]) {
                     plog << "  Error!";
                     tbox::pout << "Error in sum reduce result for group "
                                << gi << std::endl;
                     ++fail_count;
                  } else ++pass_count;
               } else {
                  plog << "  Not owner (not checking)";
               }
               plog << std::endl;
            }
            for (ai = 0; ai < num_active_groups; ++ai) {
               TBOX_ASSERT(comm_groups[ai].isDone());
            }
            TBOX_ASSERT(!comm_stage.hasPendingRequests());
            ++sync_sumreduce_count;
         }

         if (asyncany_bcast_count < asyncany_bcast_cycles) {
            TBOX_ASSERT(!comm_stage.hasCompletedMembers());
            /*
             * For the advanceSome broadcast test,
             * each group broadcasts the its group id.
             */
            plog << "\n\n\n*********** advanceAny Broadcast "
                 << asyncany_bcast_count << " ************\n";
            for (ai = 0; ai < num_active_groups; ++ai)
               if (rank != owners[active_groups[ai]]) bcdata[ai] = -1;
            plog << "Job Group Result Correct  Note\n";
            ai = 0;
            int counter = 0;
            while (counter < num_active_groups ||
                   comm_stage.hasPendingRequests()) {
               if (counter < num_active_groups) {
                  if (comm_groups[counter].beginBcast(&bcdata[counter], 1)) {
                     TBOX_ASSERT(comm_groups[counter].isDone());
                  }
                  ++counter;
               } else {
                  comm_stage.advanceAny();
               }
               if (comm_stage.hasCompletedMembers()) {
                  AsyncCommGroup* completed_group =
                     CPP_CAST<AsyncCommGroup *>(comm_stage.popCompletionQueue());
                  TBOX_ASSERT(completed_group);
                  ai = static_cast<int>(completed_group - comm_groups);
                  gi = active_groups[ai];
                  plog << std::setw(3) << ai
                       << std::setw(5) << gi
                       << std::setw(8) << bcdata[ai]
                       << std::setw(8) << correct_bcdata[ai]
                  ;
                  plog << "  Bcast difference = "
                       << bcdata[ai] - correct_bcdata[ai];
                  if (bcdata[ai] != correct_bcdata[ai]) {
                     plog << "  Error!";
                     tbox::pout << "Error in bcast result for group "
                                << gi << std::endl;
                     ++fail_count;
                  } else ++pass_count;
                  plog << std::endl;
                  TBOX_ASSERT(comm_groups[ai].isDone());
               }
               ++ai;
            }
            for (ai = 0; ai < num_active_groups; ++ai) {
               TBOX_ASSERT(comm_groups[ai].isDone());
            }
            TBOX_ASSERT(!comm_stage.hasPendingRequests());
            ++asyncany_bcast_count;
         }

         if (asyncany_sumreduce_count < asyncany_sumreduce_cycles) {
            TBOX_ASSERT(!comm_stage.hasCompletedMembers());
            /*
             * For the advanceSome broadcast test,
             * each group broadcasts the its group id.
             */
            plog << "\n\n\n*********** advanceAny Sum Reduce "
                 << asyncany_sumreduce_count << " ************\n";
            for (ai = 0; ai < num_active_groups; ++ai) sum[ai] = 1 + rank;
            plog << "Job Group Result Correct  Note\n";
            ai = 0;
            int counter = 0;
            while (counter < num_active_groups ||
                   comm_stage.hasPendingRequests()) {
               if (counter < num_active_groups) {
                  if (comm_groups[counter].beginSumReduce(&sum[counter], 1)) {
                     TBOX_ASSERT(comm_groups[counter].isDone());
                     comm_groups[counter].pushToCompletionQueue();
                  }
                  ++counter;
               }
               if (!comm_stage.hasCompletedMembers()) {
                  comm_stage.advanceAny();
                  TBOX_ASSERT(comm_stage.numberOfCompletedMembers() < 2);
               }
               TBOX_ASSERT(comm_stage.numberOfCompletedMembers() < 2);
               if (comm_stage.hasCompletedMembers()) {
                  AsyncCommGroup* completed_group =
                     CPP_CAST<AsyncCommGroup *>(comm_stage.popCompletionQueue());
                  TBOX_ASSERT(completed_group);
                  ai = static_cast<int>(completed_group - comm_groups);
                  gi = active_groups[ai];
                  plog << std::setw(3) << ai
                       << std::setw(5) << gi
                       << std::setw(8) << sum[ai]
                       << std::setw(8) << correct_sum[ai]
                  ;
                  if (rank == owners[gi]) {
                     plog << "  Sum reduce difference = "
                          << sum[ai] - correct_sum[ai];
                     if (sum[ai] != correct_sum[ai]) {
                        plog << "  Error!";
                        tbox::pout << "Error in sum reduce result for group "
                                   << gi << std::endl;
                        ++fail_count;
                     } else ++pass_count;
                  } else {
                     plog << "  Not owner (not checking)";
                  }
                  plog << std::endl;
                  TBOX_ASSERT(comm_groups[ai].isDone());
               }
               ++ai;
            }
            for (ai = 0; ai < num_active_groups; ++ai) {
               TBOX_ASSERT(comm_groups[ai].isDone());
            }
            TBOX_ASSERT(!comm_stage.hasPendingRequests());
            ++asyncany_sumreduce_count;
         }

         if (asyncsome_bcast_count < asyncsome_bcast_cycles) {
            TBOX_ASSERT(!comm_stage.hasCompletedMembers());
            /*
             * For the advanceSome broadcast test,
             * each group broadcasts the its group id.
             */
            plog << "\n\n\n*********** advanceSome Broadcast "
                 << asyncsome_bcast_count << " ************\n";
            for (ai = 0; ai < num_active_groups; ++ai)
               if (rank != owners[active_groups[ai]]) bcdata[ai] = -1;
            for (ai = 0; ai < num_active_groups; ++ai) {
               AsyncCommGroup& comm_group = comm_groups[ai];
               comm_group.beginBcast(&bcdata[ai], 1);
               if (comm_group.isDone()) {
                  comm_group.pushToCompletionQueue();
               }
            }
            plog << "Job Group Result Correct  Note\n";
            while (comm_stage.hasCompletedMembers() ||
                   comm_stage.advanceSome()) {
               AsyncCommGroup* completed_group =
                  CPP_CAST<AsyncCommGroup *>(comm_stage.popCompletionQueue());
               TBOX_ASSERT(completed_group != 0);
               ai = static_cast<int>(completed_group - comm_groups);
               gi = active_groups[ai];
               plog << std::setw(3) << ai
                    << std::setw(5) << gi
                    << std::setw(8) << bcdata[ai]
                    << std::setw(8) << correct_bcdata[ai]
               ;
               plog << "  Bcast difference = "
                    << bcdata[ai] - correct_bcdata[ai];
               if (bcdata[ai] != correct_bcdata[ai]) {
                  plog << "  Error!";
                  tbox::pout << "Error in bcast result for group "
                             << gi << std::endl;
                  ++fail_count;
               } else ++pass_count;
               plog << std::endl;
               TBOX_ASSERT(comm_groups[ai].isDone());
            }
            for (ai = 0; ai < num_active_groups; ++ai) {
               TBOX_ASSERT(comm_groups[ai].isDone());
            }
            TBOX_ASSERT(!comm_stage.hasPendingRequests());
            ++asyncsome_bcast_count;
         }

         if (asyncsome_sumreduce_count < asyncsome_sumreduce_cycles) {
            TBOX_ASSERT(!comm_stage.hasCompletedMembers());
            /*
             * For the sum advanceSome reduce test,
             * each group sums up the ranks of its members, plus 1.
             */
            plog << "\n\n\n*********** advanceSome Sum Reduce "
                 << asyncsome_sumreduce_count << " ************\n";
            for (ai = 0; ai < num_active_groups; ++ai) sum[ai] = 1 + rank;
            for (ai = 0; ai < num_active_groups; ++ai) {
               AsyncCommGroup& comm_group = comm_groups[ai];
               comm_group.beginSumReduce(&sum[ai], 1);
               if (comm_group.isDone()) {
                  comm_group.pushToCompletionQueue();
               }
            }
            plog << "Job Group Result Correct  Note\n";
            while (comm_stage.hasCompletedMembers() ||
                   comm_stage.advanceSome()) {
               AsyncCommGroup* completed_group =
                  CPP_CAST<AsyncCommGroup *>(comm_stage.popCompletionQueue());
               TBOX_ASSERT(completed_group != 0);
               ai = static_cast<int>(completed_group - comm_groups);
               gi = active_groups[ai];
               plog << std::setw(3) << ai
                    << std::setw(5) << gi
                    << std::setw(8) << sum[ai]
                    << std::setw(8) << correct_sum[ai]
               ;
               if (rank == owners[gi]) {
                  plog << "  Sum reduce difference = "
                       << sum[ai] - correct_sum[ai];
                  if (sum[ai] != correct_sum[ai]) {
                     plog << "  Error!";
                     tbox::pout << "Error in sum reduce result for group "
                                << gi << std::endl;
                     ++fail_count;
                  } else ++pass_count;
               } else {
                  plog << "  Not owner (not checking)";
               }
               plog << std::endl;
               TBOX_ASSERT(comm_groups[ai].isDone());
            }
            for (ai = 0; ai < num_active_groups; ++ai) {
               TBOX_ASSERT(comm_groups[ai].isDone());
            }
            TBOX_ASSERT(!comm_stage.hasPendingRequests());
            ++asyncsome_sumreduce_count;
         }

         ++count;
         delete[] comm_groups;
      }

      plog << '\n';
      plog << "pass_count = " << pass_count << std::endl;
      plog << "fail_count = " << fail_count << std::endl;
      plog << "\n************** Test completed **************\n" << std::endl;
      input_db->printClassData(tbox::plog);

      /*
       * Clean up and exit.
       */

      TimerManager::getManager()->print(plog);

   }

   if (fail_count == 0) {
      tbox::pout << "\nPASSED:  group_comm" << std::endl;
   }

   SAMRAIManager::shutdown();
   SAMRAIManager::finalize();
   SAMRAI_MPI::finalize();

   return fail_count;
}
