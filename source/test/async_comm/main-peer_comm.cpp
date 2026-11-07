/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test program for asynchronous peer communication classes
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <iomanip>
#include <vector>

using namespace SAMRAI;
using namespace tbox;

/*
 ************************************************************************
 * Define data to communicate, based on sender, recipient and count.
 ************************************************************************
 */
template<class TYPE>
class TypeIndependentTester
{

public:
   void
   runTest(
      int& pass_count,
      int& fail_count,
      tbox::SAMRAI_MPI::Comm isolated_communicator,
      int max_first_data_length,
      int use_advance_some,
      int num_cycles,
      int group_rel_first,
      int group_rel_last);
   void setSendData(
      int fr,
      int to,
      int count,
      std::vector<TYPE>& send_data) {
      NULL_USE(to);
      send_data.resize(count);
      for (int i = 0; i < count; ++i) send_data[i] = (TYPE)(fr + count + i);
   }
   bool checkRecvData(
      int fr,
      int to,
      int count,
      int recv_size,
      const TYPE* recv_data,
      std::string& size_correct,
      std::string& data_correct) {
      NULL_USE(to);

      bool rval = false;
      if (recv_size == count) size_correct = "SIZE OK";
      else {
         size_correct = "WRONG SIZE";
         rval = true;
      }
      data_correct = "DATA OK";
      for (int i = 0; i < count; ++i) {
         if (!tbox::MathUtilities<double>::equalEps(recv_data[i], fr + count
                + i)) {
            data_correct = " WRONG DATA";
            rval = true;
            break;
         }
      }
      return rval;
   }
};

template<class TYPE>
void TypeIndependentTester<TYPE>::runTest(
   int& pass_count,
   int& fail_count,
   tbox::SAMRAI_MPI::Comm isolated_communicator,
   int max_first_data_length,
   int use_advance_some,
   int num_cycles,
   int group_rel_first,
   int group_rel_last)
{
   tbox::SAMRAI_MPI mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   const int nproc = mpi.getSize();
   const int iproc = mpi.getRank();

   /*
    * Compute group size from the range given by group_rel_first
    * and group_rel_last, but to avoid duplicating members, don't
    * let group get bigger than nproc.
    */
   if (group_rel_last < group_rel_first) {
      TBOX_ERROR("Invalid input: group_rel_last < group_rel_first");
   }
   const int group_size =
      tbox::MathUtilities<int>::Min(group_rel_last - group_rel_first + 1, nproc);

   plog << num_cycles << " cycles." << std::endl;
   plog << "num procs = " << nproc << std::endl;
   plog << "group size = " << group_size << std::endl;

   pass_count = 0;
   fail_count = 0;

   if (group_size < 2) {
      plog << "Bypassing TypeIndependentTester::runTest due to trivial group size.\n";
      return;
   }

   /*
    * Allocate 2*group_size communication objects.  The group_size objects
    * are senders, the others are receivers.
    */
   AsyncCommPeer<TYPE>* peer_comms = new AsyncCommPeer<TYPE>[2 * group_size];
   AsyncCommStage stage;

   // Counter for number of completions of each communication group.
   std::vector<int> completion_counter(2 * group_size);

   for (int i = 0; i < 2 * group_size; ++i) {
      peer_comms[i].initialize(&stage);
      peer_comms[i].setMPI(SAMRAI_MPI(isolated_communicator));
      peer_comms[i].setMPITag(0, 1);
      peer_comms[i].limitFirstDataLength(max_first_data_length);
      completion_counter[i] = 0;
   }

   /*
    * First half of peer_comms is used for sending.
    * Second half is used for receiving.
    * Set peer ranks for each object in peer_comms.
    */

   for (int i = 0; i < group_size; ++i) {
      int peer_rank = iproc + group_rel_first + i;
      while (peer_rank < 0) peer_rank += nproc;
      peer_rank %= nproc;
      peer_comms[i].setPeerRank(peer_rank);
      plog << "Proc " << std::setw(3) << iproc
           << " peer_comms[" << std::setw(3) << i
           << "] send to " << std::setw(3) << peer_rank << std::endl;
   }

   for (int i = 0; i < group_size; ++i) {
      int peer_rank = iproc - group_rel_last + i;
      while (peer_rank < 0) peer_rank += nproc;
      peer_rank %= nproc;
      peer_comms[i + group_size].setPeerRank(peer_rank);
      plog << "Proc " << std::setw(3) << iproc
           << " peer_comms[" << std::setw(3) << i + group_size
           << "] recv fr " << std::setw(3) << peer_rank << std::endl;
   }

   /*
    * Test loop.  Each process will send and receive from every
    * member in its group num_cycles times.
    *
    * completion_counter tracks how many
    * of these operations have completed.
    */
   int count = 0;
   while (count < 2 * group_size * num_cycles) {

      /*
       * Find a list of completed members using either advanceSome
       * or advanceAny, as controlled by input file.
       */
      if (use_advance_some) {
         stage.advanceSome();
      } else {
         stage.advanceAny();
      }

      /*
       * Check completed members for correctness.
       */
      while (stage.hasCompletedMembers()) {

         AsyncCommStage::Member* completed_member = stage.popCompletionQueue();

         /*
          * If there has been a completed prop, process it.
          * Else process everything that has not reached
          * the required number of test cycles.
          */

         AsyncCommPeer<TYPE>* completed_comm_ =
            CPP_CAST<AsyncCommPeer<TYPE> *>(completed_member);
         TBOX_ASSERT(completed_comm_);
         AsyncCommPeer<TYPE>& completed_comm = *completed_comm_;

         int completed_comm_index = static_cast<int>(completed_comm_ - peer_comms);

         /*
          * Whether completed_comm is a sender or receiver is based on its index in peer_comms.
          */
         if (completed_comm_index < group_size) {
            // completed_comm is a sender.  No accuracy checks needed.
            plog << "comm_peer[" << std::setw(3) << completed_comm_index
                 << "] finished send # " << std::setw(3)
                 << completion_counter[completed_comm_index]
                 << " to " << std::setw(3) << completed_comm.getPeerRank() << "."
                 << std::endl;
         } else {
            // completed_comm is a receiver.  Do accuracy check on received data.
            std::string size_correct, data_correct;
            bool fail = checkRecvData(completed_comm.getPeerRank(),
                  iproc,
                  completion_counter[completed_comm_index],
                  completed_comm.getRecvSize(),
                  completed_comm.getRecvData(),
                  size_correct,
                  data_correct);
            plog << "comm_peer[" << std::setw(3) << completed_comm_index
                 << "] finished recv # " << std::setw(3)
                 << completion_counter[completed_comm_index]
                 << " fr " << std::setw(3) << completed_comm.getPeerRank() << ": "
                 << std::setw(5) << size_correct << ' '
                 << std::setw(5) << data_correct << ' '
                 << std::endl;
            if (fail) ++fail_count;
            else ++pass_count;
         }

         // Count number of completions for the current AsyncCommPeer.
         ++completion_counter[completed_comm_index];

         // Count number of completions for the whole test.
         ++count;

      }

      /*
       * Launch another cycle for members that are done but have not
       * completed all required cycles.
       */
      for (int i = 0; i < 2 * group_size; ++i) {

         AsyncCommPeer<TYPE>& peer_comm = peer_comms[i];

         if (completion_counter[i] < num_cycles &&
             peer_comm.isDone()) {

            if (i < group_size) {
               // This is a sender.
               peer_comm.setMPITag(2 * completion_counter[i],
                  2 * completion_counter[i] + 1);
               std::vector<TYPE> send_data;
               setSendData(iproc,
                  peer_comm.getPeerRank(),
                  completion_counter[i],
                  send_data);
               peer_comm.beginSend(send_data.size() > 0 ? &send_data[0] : 0,
                  static_cast<int>(send_data.size()));
               /*
                * Check if the new communication is done (because if it is,
                * the stage won't detect it--stage only detects non-NULL request.
                */
               if (peer_comm.isDone()) {
                  plog << "comm_peer[" << std::setw(3) << i
                       << "] finished send # " << std::setw(3)
                       << completion_counter[i]
                       << " items of size " << sizeof(TYPE) << " to "
                       << std::setw(3) << peer_comm.getPeerRank() << "."
                       << std::endl;
                  ++completion_counter[i];
                  ++count;
               }
            } else {
               // This is a receiver.
               peer_comm.setMPITag(2 * completion_counter[i],
                  2 * completion_counter[i] + 1);
               peer_comm.beginRecv();
               /*
                * Check if the new communication is done (because if it is,
                * the stage won't detect it--stage only detects non-NULL request.
                */
               if (peer_comm.isDone()) {
                  std::string size_correct, data_correct;
                  bool fail = checkRecvData(peer_comm.getPeerRank(),
                        iproc,
                        completion_counter[i],
                        peer_comm.getRecvSize(),
                        peer_comm.getRecvData(),
                        size_correct,
                        data_correct);
                  plog << "comm_peer[" << std::setw(3) << i
                       << "] finished recv # " << std::setw(3)
                       << completion_counter[i]
                       << " items of size " << sizeof(TYPE) << " fr "
                       << std::setw(3) << peer_comm.getPeerRank() << ": "
                       << std::setw(5) << size_correct << ' '
                       << std::setw(5) << data_correct << ' '
                       << std::endl;
                  if (fail) ++fail_count;
                  else ++pass_count;
                  ++completion_counter[i];
                  ++count;
               }
            }

         }

      }

   }

   delete[] peer_comms;
}

template class TypeIndependentTester<int>;
template class TypeIndependentTester<float>;
template class TypeIndependentTester<double>;
template class TypeIndependentTester<char>;

/*
 ************************************************************************
 *
 * This program tests the asynchronous communication classes:
 * AsyncCommPeer
 * AsyncCommStage
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

   const int iproc = mpi.getRank();
   int total_fail_count = 0;

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
      plog << "Process " << std::setw(5) << iproc << " is ready." << std::endl;

      /*
       * Created a separate communicator for testing,
       * to avoid possible interference with other communications
       * by SAMRAI library.
       */
      tbox::SAMRAI_MPI::Comm isolated_communicator(MPI_COMM_NULL);
      if (tbox::SAMRAI_MPI::usingMPI()) {
         tbox::SAMRAI_MPI::getSAMRAIWorld().Comm_dup(&isolated_communicator);
      }
      tbox::SAMRAI_MPI isolated_mpi(isolated_communicator);
      plog << "Process " << std::setw(5) << isolated_mpi.getRank()
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

      std::shared_ptr<Database> main_db = input_db->getDatabase("Main");
      std::string base_name = "unnamed";
      base_name = main_db->getStringWithDefault("base_name", base_name);

      int max_first_data_length =
         main_db->getIntegerWithDefault("max_first_data_length", 1);

      bool use_advance_some = false;
      use_advance_some = main_db->getBoolWithDefault("use_advance_some",
            use_advance_some);

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

      const int num_cycles =
         main_db->getIntegerWithDefault("num_cycles", 1);

      const int group_rel_first =
         main_db->getIntegerWithDefault("group_rel_first", 0);
      const int group_rel_last =
         main_db->getIntegerWithDefault("group_rel_last", 1);

      int pass_count, fail_count;

      plog << "\nTesting passing integer data:\n";
      plog << "sizeof(int) = " << sizeof(int) << "\n";
      TypeIndependentTester<int> tester_int;
      tester_int.runTest(pass_count,
         fail_count,
         isolated_communicator,
         max_first_data_length,
         use_advance_some,
         num_cycles,
         group_rel_first,
         group_rel_last);
      plog << "pass_count = " << pass_count << std::endl;
      plog << "bad_count  = " << fail_count << std::endl;
      total_fail_count += fail_count;

      plog << "\nTesting passing float data:\n";
      plog << "sizeof(float) = " << sizeof(float) << "\n";
      TypeIndependentTester<float> tester_float;
      tester_float.runTest(pass_count,
         fail_count,
         isolated_communicator,
         max_first_data_length,
         use_advance_some,
         num_cycles,
         group_rel_first,
         group_rel_last);
      plog << "pass_count = " << pass_count << std::endl;
      plog << "bad_count  = " << fail_count << std::endl;
      total_fail_count += fail_count;

      plog << "\nTesting passing double data:\n";
      plog << "sizeof(double) = " << sizeof(double) << "\n";
      TypeIndependentTester<double> tester_double;
      tester_double.runTest(pass_count,
         fail_count,
         isolated_communicator,
         max_first_data_length,
         use_advance_some,
         num_cycles,
         group_rel_first,
         group_rel_last);
      plog << "pass_count = " << pass_count << std::endl;
      plog << "bad_count  = " << fail_count << std::endl;
      total_fail_count += fail_count;

      plog << "\nTesting passing char data:\n";
      plog << "sizeof(char) = " << sizeof(char) << "\n";
      TypeIndependentTester<char> tester_char;
      tester_char.runTest(pass_count,
         fail_count,
         isolated_communicator,
         max_first_data_length,
         use_advance_some,
         num_cycles,
         group_rel_first,
         group_rel_last);
      plog << "pass_count = " << pass_count << std::endl;
      plog << "bad_count  = " << fail_count << std::endl;
      total_fail_count += fail_count;

      plog << "\n************** Test completed **************\n" << std::endl;
      input_db->printClassData(tbox::plog);

      /*
       * Clean up and exit.
       */

      TimerManager::getManager()->print(plog);

#if defined(HAVE_MPI)
      MPI_Comm_free(&isolated_communicator);
#endif

      plog << "Process " << std::setw(5) << iproc << " got " << fail_count
           << " failures and " << pass_count << " successes." << std::endl;

   }

   if (total_fail_count == 0) {
      tbox::pout << "\nPASSED:  peer_comm" << std::endl;
   }

   SAMRAIManager::shutdown();
   SAMRAIManager::finalize();
   SAMRAI_MPI::finalize();

   return total_fail_count;
}
