/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test program to test Connector class
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/AssumedPartition.h"
#include "SAMRAI/geom/GridGeometry.h"



using namespace SAMRAI;
using namespace hier;

/*!
 * @brief Primitive BoxGenerator (independent of mesh package)
 * creating boxes using an AssumedPartition followed by an index
 * filter to keep a subset of boxes.
 */
struct PrimitiveBoxGen {
   std::shared_ptr<hier::BaseGridGeometry> d_geom;
   hier::AssumedPartition d_ap;
   // Index filtering parameters.
   enum IndexFilter { ALL = 0 /* Keep all boxes */,
                      INTERVAL = 1 /* Keep d_num_keep, discard d_num_discard */,
                      LOWER = 2 /* Keep indices below d_frac */,
                      UPPER = 3 /* Keep indices above d_frac */
   };
   std::string d_nickname;
   double d_avg_parts_per_rank;
   int d_index_filter;
   int d_num_keep;
   int d_num_discard;
   double d_frac;
   PrimitiveBoxGen(
      tbox::Database& database,
      const std::shared_ptr<hier::BaseGridGeometry>& geom):
      d_avg_parts_per_rank(1.0)
   {
      d_geom = geom;
      getFromInput(database);
   }
   PrimitiveBoxGen(
      const PrimitiveBoxGen& other):
      d_geom(other.d_geom),
      d_ap(other.d_ap),
      d_avg_parts_per_rank(other.d_avg_parts_per_rank),
      d_index_filter(other.d_index_filter),
      d_num_keep(other.d_num_keep),
      d_num_discard(other.d_num_discard),
      d_frac(other.d_frac) {
   }
   void
   getFromInput(
      tbox::Database& input_db);
   void
   getBoxes(
      hier::BoxContainer& boxes);
};

// Parameters in a specific test.
struct CommonTestParams {
   std::string d_nickname;
   int d_base_num;
   int d_head_num;
   std::string d_method;
   // For method mod:
   int d_denom;
   // For method bracket:
   int d_begin_shift;
   int d_end_shift;
   int d_inc;
   // For method overlap: (no parameter)
   CommonTestParams(
      tbox::Database& test_db);
   void
   contriveConnector(
      Connector& conn,
      const std::vector<PrimitiveBoxGen>& boxgens,
      const std::vector<hier::BoxLevel>& levels,
      const hier::IntVector& connector_width);
};

CommonTestParams
getTestParametersFromDatabase(
   tbox::Database& test_db);

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

         const hier::IntVector& refinement_ratio = hier::IntVector::getOne(dim);

         hier::IntVector connector_width(dim);
         main_db->getIntegerArray("connector_width",
            &connector_width[0],
            connector_width.getDim().getValue());

         if (!input_db->isDatabase("BlockGeometry")) {
            TBOX_ERROR(
               "getTestParametersFromDatabase: You must specify \"BlockGeometry\" in input database.");
         }
         // Note: Using GridGeometry only because BaseGridGeometry can't be instanstiated.
         std::shared_ptr<hier::BaseGridGeometry> grid_geom =
            std::make_shared<geom::GridGeometry>(
               dim,
               "BlockGeometry",
               input_db->getDatabase("BlockGeometry"));

         hier::BoxLevelConnectorUtils blcu;

         /*
          * Generate BoxLevels and associated PrimitiveBoxGen objects.
          */

         std::vector<hier::BoxLevel> levels;
         std::vector<PrimitiveBoxGen> boxgens;

         while (true) {

            std::string level_name("PrimitiveBoxGen");
            level_name += tbox::Utilities::intToString(static_cast<int>(levels.size()), 1);

            std::shared_ptr<tbox::Database> level_db =
               input_db->getDatabaseWithDefault(level_name, std::shared_ptr<tbox::Database>());

            if (!level_db) {
               break;
            }

            PrimitiveBoxGen boxgen(*level_db, grid_geom);
            BoxContainer boxes;
            boxgen.getBoxes(boxes);
            BoxLevel box_level(boxes, refinement_ratio, grid_geom);
            blcu.addPeriodicImages(box_level,
               grid_geom->getPeriodicDomainSearchTree(), connector_width);
            box_level.cacheGlobalReducedData();

            boxgens.push_back(boxgen);
            levels.push_back(box_level);

            const BoxLevel& globalized = levels.back().getGlobalizedVersion();
            if (rank == 0) {
               tbox::plog << "Globalized version of BoxLevel #" << levels.size() - 1 << ":\n"
                          << globalized.format("\t");
            }
         }

         /*
          * Read in and run tests.
          */

         int test_number = 0;
         while (true) {

            std::string test_name("Test");
            test_name += tbox::Utilities::intToString(test_number++, 2);

            std::shared_ptr<tbox::Database> test_db =
               input_db->getDatabaseWithDefault(test_name, std::shared_ptr<tbox::Database>());

            if (!test_db) {
               break;
            }

            CommonTestParams testparams(*test_db);

            tbox::plog << "\n\n";
            tbox::pout << "Running " << test_name << " (" << testparams.d_nickname << ")\n";

            /*
             * Set up edges in forward by the contrivance specified
             * in the test database.  Then check transpose
             * correctness.
             */

            hier::Connector forward(dim);
            testparams.contriveConnector(forward, boxgens, levels, connector_width);

            tbox::plog << "Testing with:"
                       << "\nbase:\n" << forward.getBase().format("\t")
                       << "\nhead:\n" << forward.getHead().format("\t")
                       << "\nforward:\n" << forward.format("\t")
                       << std::endl;

            hier::Connector reverse(forward.getHead(), forward.getBase(), connector_width);
            reverse.computeTransposeOf(forward);
            tbox::plog << "Computed:\nreverse:\n" << reverse.format("\t")
                       << std::endl;

            size_t test_fail_count = forward.checkTransposeCorrectness(reverse);
            fail_count += static_cast<int>(test_fail_count);
            if (test_fail_count) {
               tbox::pout << "FAILED: " << test_name << " (" << testparams.d_nickname << ')'
                          << std::endl;
               tbox::plog << "FAILED: " << test_name << " (" << testparams.d_nickname << ')'
                          << std::endl;
            } else {
               tbox::plog << "PASSED: " << test_name << " (" << testparams.d_nickname << ')'
                          << std::endl;
            }
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
void PrimitiveBoxGen::getFromInput(tbox::Database& test_db)
{
   int rank_begin = 0;
   int rank_end = tbox::SAMRAI_MPI::getSAMRAIWorld().getSize();
   int index_begin = test_db.getIntegerWithDefault("index_begin", 0);

   d_avg_parts_per_rank = test_db.getDoubleWithDefault("avg_parts_per_rank", d_avg_parts_per_rank);

   d_nickname = test_db.getString("nickname");

   std::string index_filter = test_db.getStringWithDefault("index_filter", "ALL");
   if (index_filter == "ALL") {
      d_index_filter = PrimitiveBoxGen::ALL;
   } else if (index_filter == "INTERVAL") {
      d_index_filter = PrimitiveBoxGen::INTERVAL;
      d_num_keep = test_db.getInteger("num_keep");
      d_num_discard = test_db.getInteger("num_discard");
   } else if (index_filter == "LOWER") {
      d_index_filter = PrimitiveBoxGen::LOWER;
      d_frac = test_db.getDouble("frac");
   } else if (index_filter == "UPPER") {
      d_index_filter = PrimitiveBoxGen::UPPER;
      d_frac = test_db.getDouble("frac");
   }

   d_ap.partition(d_geom->getPhysicalDomain(),
      rank_begin, rank_end, index_begin,
      d_avg_parts_per_rank);

   tbox::plog << "PrimitiveBoxGen::getFromInput():\n"
              << "d_nickname = " << d_nickname << "\n"
              << "d_avg_parts_per_rank = " << d_avg_parts_per_rank << "\n"
              << "rank begin, end = " << rank_begin << "  " << rank_end << '\n'
              << "index_begin = " << index_begin << "\n"
              << "index_filter = " << index_filter << "\n"
              << "d_num_keep = " << d_num_keep << "\n"
              << "d_num_discard = " << d_num_discard << "\n"
              << "d_frac = " << d_frac << "\n"
              << "generated AssumedPartition:\n";
   d_ap.recursivePrint(tbox::plog, "\t", 3);

   if (d_ap.selfCheck()) {
      TBOX_ERROR("Error in setting up AssumedPartition d_ap (selfCheck failed).\n");
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void PrimitiveBoxGen::getBoxes(hier::BoxContainer& boxes)
{
   if (d_index_filter == ALL) {
      int idbegin = d_ap.begin();
      int idend = d_ap.end();
      for (int id = idbegin; id < idend; ++id) {
         boxes.push_back(d_ap.getBox(id));
      }
   } else if (d_index_filter == INTERVAL) {
      int idbegin = d_ap.begin();
      int idend = d_ap.end();
      int interval = d_num_keep + d_num_discard;
      for (int id = idbegin; id < idend; ++id) {
         int interval_id = id % interval;
         if (interval_id < d_num_keep) {
            boxes.push_back(d_ap.getBox(id));
         }
      }
   } else if (d_index_filter == LOWER) {
      int threshold = d_ap.begin() + static_cast<int>(d_frac * (d_ap.end() - d_ap.begin()));
      int idbegin = d_ap.begin();
      int idend = tbox::MathUtilities<int>::Min(threshold, d_ap.end());
      for (int id = idbegin; id < idend; ++id) {
         boxes.push_back(d_ap.getBox(id));
      }
   } else if (d_index_filter == UPPER) {
      int threshold = d_ap.begin() + static_cast<int>(d_frac * (d_ap.end() - d_ap.begin()));
      int idbegin = tbox::MathUtilities<int>::Max(threshold, d_ap.begin());
      int idend = d_ap.end();
      for (int id = idbegin; id < idend; ++id) {
         boxes.push_back(d_ap.getBox(id));
      }
   } else {
      TBOX_ERROR("Invalid value of index_filter: " << d_index_filter);
   }
}

CommonTestParams::CommonTestParams(
   tbox::Database& test_db):
   d_denom(0),
   d_begin_shift(0),
   d_end_shift(0),
   d_inc(0)
{
   d_nickname = test_db.getString("nickname");
   d_method = test_db.getString("method");

   int level_num[2];
   test_db.getIntegerArray("levels", level_num, 2);
   d_base_num = level_num[0];
   d_head_num = level_num[1];

   if (d_method == "mod") {
      d_denom = test_db.getInteger("denom");
   } else if (d_method == "bracket") {
      d_begin_shift = test_db.getInteger("begin_shift");
      d_end_shift = test_db.getInteger("end_shift");
      d_inc = test_db.getInteger("inc");
   } else if (d_method == "overlap") {
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void CommonTestParams::contriveConnector(
   Connector& conn,
   const std::vector<PrimitiveBoxGen>& boxgens,
   const std::vector<hier::BoxLevel>& levels,
   const hier::IntVector& connector_width)
{
   const PrimitiveBoxGen& boxgen_base = boxgens[d_base_num];
   const PrimitiveBoxGen& boxgen_head = boxgens[d_head_num];
   const hier::BoxLevel& base = levels[d_base_num];
   const hier::BoxLevel& head = levels[d_head_num];

   conn = Connector(base, head, connector_width);

   const int rank = conn.getBase().getMPI().getRank();

   if (d_method == "mod") {
      for (int i = boxgen_base.d_ap.beginOfRank(rank); i < boxgen_base.d_ap.endOfRank(rank); ++i) {
         hier::Box l1box = boxgen_base.d_ap.getBox(i);
         TBOX_ASSERT(l1box.getOwnerRank() == rank);
         for (int j = boxgen_head.d_ap.begin(); j < boxgen_head.d_ap.end(); ++j) {
            hier::Box l2box = boxgen_head.d_ap.getBox(j);
            TBOX_ASSERT(l2box.getOwnerRank() >= 0 &&
               l2box.getOwnerRank() < conn.getBase().getMPI().getSize());
            if ((i + j) % d_denom == 0) {
               conn.insertLocalNeighbor(l2box, l1box.getBoxId());
            }
         }
      }
      tbox::plog << "Contrived connector using 'mod':"
                 << "  denom=" << d_denom
                 << std::endl;
   } else if (d_method == "bracket") {
      for (int i = boxgen_base.d_ap.beginOfRank(rank); i < boxgen_base.d_ap.endOfRank(rank); ++i) {
         hier::Box l1box = boxgen_base.d_ap.getBox(i);
         TBOX_ASSERT(l1box.getOwnerRank() == rank);

         int begin = l1box.getLocalId().getValue() + d_begin_shift;
         int end = l1box.getLocalId().getValue() + d_end_shift;
         begin = tbox::MathUtilities<int>::Max(begin, boxgen_head.d_ap.begin());
         end = tbox::MathUtilities<int>::Min(end, boxgen_head.d_ap.end());

         for (int j = begin; j < end; j += d_inc) {
            hier::Box l2box = boxgen_head.d_ap.getBox(j);
            TBOX_ASSERT(l2box.getOwnerRank() >= 0 &&
               l2box.getOwnerRank() < conn.getBase().getMPI().getSize());
            conn.insertLocalNeighbor(l2box, l1box.getBoxId());
         }
      }
      tbox::plog << "Contrived connector using 'bracket':"
                 << "  begin_shift=" << d_begin_shift
                 << "  end_shift=" << d_end_shift
                 << "  inc=" << d_inc
                 << std::endl;
   } else if (d_method == "overlap") {
      hier::OverlapConnectorAlgorithm oca;
      oca.findOverlaps(conn);
      tbox::plog << "Contrived connector using 'overlap':"
                 << std::endl;
   } else {
      TBOX_ERROR("Contrivance method must be one of these: mod, bracket.");
   }
}
