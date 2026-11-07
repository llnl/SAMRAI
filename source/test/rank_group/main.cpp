/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test program for RankGroup with TreeLoadBalancer.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include <iomanip>

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/BoxLevelStatistics.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/ConnectorStatistics.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/appu/VisItDataWriter.h"

#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include <vector>
#include <algorithm>

#include <cmath>

using namespace SAMRAI;
using namespace tbox;

/*
 ************************************************************************
 *
 *
 *************************************************************************
 */

void
generatePrebalanceByUserBoxes(
   std::shared_ptr<tbox::Database> database,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const hier::IntVector& min_size,
   const hier::IntVector& max_gcw,
   std::shared_ptr<hier::BoxLevel>& balance_box_level,
   const hier::BoxLevel& anchor_box_level,
   std::shared_ptr<hier::Connector>& anchor_to_balance);

void
generatePrebalanceByUserShells(
   const std::shared_ptr<tbox::Database>& database,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const hier::IntVector& min_size,
   const hier::IntVector& max_gcw,
   std::shared_ptr<hier::BoxLevel>& balance_box_level,
   const std::shared_ptr<hier::BoxLevel>& anchor_box_level,
   std::shared_ptr<hier::Connector>& anchor_to_balance,
   int tag_level_number);

void
sortNodes(
   hier::BoxLevel& new_box_level,
   hier::Connector& tag_to_new,
   bool sort_by_corners,
   bool sequentialize_global_indices);

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
   tbox::SAMRAI_MPI mpi(SAMRAI_MPI::getSAMRAIWorld());

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

   {
      /*
       * Scope to force destruction of objects that would otherwise
       * leave allocated memory reported by the memory test.
       */

      /*
       * Create input database and parse all data in input file.
       */

      std::shared_ptr<InputDatabase> input_db(
         new InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

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

      const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));
      const int dimval = dim.getValue();

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
            mpi.getSize(),
            5);
      tbox::plog << "Added case name (" << case_name << ") and nprocs ("
                 << mpi.getSize() << ") to base name -> '"
                 << base_name << "'\n";

      if (!case_name.empty()) {
         tbox::plog << "Added case name (" << case_name << ") and nprocs ("
                    << mpi.getSize() << ") to base name -> '"
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

      plog << "Input database after initialization..." << std::endl;
      input_db->printClassData(plog);

      // tbox::TimerManager *tm = tbox::TimerManager::getManager();
      // std::shared_ptr<tbox::Timer> t_search_tree_for_set =
      // tm->getTimer("apps::main::search_tree_for_set");

      /*
       * Parameters.  Some of these should be specified by input deck.
       */
      hier::IntVector ghost_cell_width(dim, 2);
      if (main_db->isInteger("ghost_cell_width")) {
         main_db->getIntegerArray("ghost_cell_width", &ghost_cell_width[0], dimval);
      }

      hier::IntVector min_size(dim, 8);
      if (main_db->isInteger("min_size")) {
         main_db->getIntegerArray("min_size", &min_size[0], dimval);
      }
      hier::IntVector max_size(dim, tbox::MathUtilities<int>::getMax());
      if (main_db->isInteger("max_size")) {
         main_db->getIntegerArray("max_size", &max_size[0], dimval);
      }
      hier::IntVector bad_interval(dim, 2);
      hier::IntVector cut_factor(dim, 1);

      hier::OverlapConnectorAlgorithm oca;

      /*
       * Set up the domain from input.
       */
      std::vector<tbox::DatabaseBox> db_box_vector =
         main_db->getDatabaseBoxVector("domain_boxes");
      hier::BoxContainer domain_boxes(db_box_vector);
      for (hier::BoxContainer::iterator itr = domain_boxes.begin();
           itr != domain_boxes.end(); ++itr) {
         itr->setBlockId(hier::BlockId(0));
      }

      /*
       * Create hierarchy we can satisfy the load balancing
       * interface and dump visit output.
       *
       * anchor_box_level is used for level 0.
       * balance_box_level is used for level 1.
       */
      std::vector<double> xlo(dimval);
      std::vector<double> xhi(dimval);
      for (int i = 0; i < dimval; ++i) {
         xlo[i] = 0.0;
         xhi[i] = 1.0;
      }
      std::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
         new geom::CartesianGridGeometry(
            "GridGeometry",
            &xlo[0],
            &xhi[0],
            domain_boxes));

      std::shared_ptr<hier::PatchHierarchy> hierarchy(
         new hier::PatchHierarchy(
            "Hierarchy",
            grid_geometry));

      hierarchy->setMaxNumberOfLevels(2);

      hier::BoxLevel domain_box_level(
         hier::IntVector(dim, 1),
         grid_geometry,
         mpi,
         hier::BoxLevel::GLOBALIZED);
      for (hier::BoxContainer::iterator domain_boxes_itr = domain_boxes.begin();
           domain_boxes_itr != domain_boxes.end(); ++domain_boxes_itr) {
         domain_box_level.addBox(*domain_boxes_itr);
      }

      /*
       * Set up the load balancers.
       */

      mesh::ChopAndPackLoadBalancer cut_and_pack_lb(
         dim,
         "ChopAndPackLoadBalancer",
         input_db->getDatabaseWithDefault("ChopAndPackLoadBalancer",
            std::shared_ptr<tbox::Database>()));

      mesh::TreeLoadBalancer tree_lb(
         dim,
         "TreeLoadBalancer",
         input_db->getDatabaseWithDefault("TreeLoadBalancer",
            std::shared_ptr<tbox::Database>()));
      tree_lb.setSAMRAI_MPI(mpi);

      mesh::LoadBalanceStrategy* lb = 0;

      std::string load_balancer_type;
      if (main_db->isString("load_balancer_type")) {
         load_balancer_type = main_db->getString("load_balancer_type");
         if (load_balancer_type == "TreeLoadBalancer") {
            lb = &tree_lb;
         } else if (load_balancer_type == "ChopAndPackLoadBalancer") {
            lb = &cut_and_pack_lb;
         }
      }
      if (lb == 0) {
         TBOX_ERROR(
            "Missing or bad load_balancer specification in Main database.\n"
            << "Specify load_balancer_type = STRING, where STRING can be\n"
            << "\"ChopAndPackLoadBalancer\" or \"TreeLoadBalancer\".");
      }

      /*
       * Set up data used by TreeLoadBalancer.
       */
      std::shared_ptr<hier::BoxLevel> anchor_box_level(
         std::make_shared<hier::BoxLevel>(
            hier::IntVector(dim, 1), grid_geometry));
      std::shared_ptr<hier::BoxLevel> balance_box_level;
      std::shared_ptr<hier::Connector> anchor_to_balance;
      std::shared_ptr<hier::Connector> balance_to_balance;

      {
         std::vector<tbox::DatabaseBox> db_box_vector =
            main_db->getDatabaseBoxVector("anchor_boxes");
         hier::BoxContainer anchor_boxes(db_box_vector);
         const int boxes_per_proc =
            (anchor_boxes.size() + anchor_box_level->getMPI().getSize() - 1)
            / anchor_box_level->getMPI().getSize();
         const int my_boxes_start = anchor_box_level->getMPI().getRank()
            * boxes_per_proc;
         const int my_boxes_stop =
            tbox::MathUtilities<int>::Min(my_boxes_start + boxes_per_proc,
               anchor_boxes.size());
         hier::BoxContainer::iterator anchor_boxes_itr = anchor_boxes.begin();
         for (int i = 0; i < my_boxes_start; ++i) {
            if (anchor_boxes_itr != anchor_boxes.end()) {
               ++anchor_boxes_itr;
            }
         }
         for (int i = my_boxes_start; i < my_boxes_stop; ++i, ++anchor_boxes_itr) {
            anchor_box_level->addBox(*anchor_boxes_itr, hier::BlockId::zero());
         }
      }

      {
         /*
          * Load balance the anchor box_level, using the domain as its anchor.
          *
          * This is not a part of the performance test because does not
          * reflect the load balancer use in real apps.  We just neeed a
          * distributed anchor for the real loac balancing performance test.
          */
         std::shared_ptr<hier::Connector> domain_to_anchor;
         oca.findOverlapsWithTranspose(domain_to_anchor,
            domain_box_level,
            *anchor_box_level,
            hier::IntVector(dim, 2),
            hier::IntVector(dim, 2));
         hier::Connector* anchor_to_domain = &domain_to_anchor->getTranspose();

         tbox::plog << "\n\n\ninitial anchor loads:\n";
         mesh::BalanceUtilities::reduceAndReportLoadBalance(
            std::vector<double>(1, static_cast<double>(anchor_box_level->getLocalNumberOfCells())),
            anchor_box_level->getMPI());

         const int nnodes = mpi.getSize();
         std::vector<int> active_ranks;
         if (nnodes == 1) {
            active_ranks.resize(1);
            active_ranks[0] = 0;
         } else {
            active_ranks.resize(nnodes / 2);
            for (int i = 0; i < nnodes / 2; ++i) {
               active_ranks[i] = (i + 1) % (nnodes / 2);
            }
            std::sort(&active_ranks[0],
               &active_ranks[0] + static_cast<int>(active_ranks.size()));
         }
         tbox::RankGroup rank_group_0(active_ranks, mpi);

         lb->loadBalanceBoxLevel(
            *anchor_box_level,
            anchor_to_domain,
            hierarchy,
            0,
            min_size,
            max_size,
            domain_box_level,
            bad_interval,
            cut_factor,
            rank_group_0);

         sortNodes(*anchor_box_level,
            *domain_to_anchor,
            false,
            true);

         anchor_to_domain->assertOverlapCorrectness();
         domain_to_anchor->assertOverlapCorrectness();

         anchor_box_level->cacheGlobalReducedData();

         tbox::plog << "\n\n\nfinal anchor loads:\n";
         mesh::BalanceUtilities::reduceAndReportLoadBalance(
            std::vector<double>(1, static_cast<double>(anchor_box_level->getLocalNumberOfCells())),
            anchor_box_level->getMPI());
      }

      {
         std::string box_gen_method("PrebalanceByUserBoxes");
         box_gen_method = main_db->getStringWithDefault("box_gen_method",
               box_gen_method);
         if (box_gen_method == "PrebalanceByUserBoxes") {
            generatePrebalanceByUserBoxes(
               main_db->getDatabase("PrebalanceByUserBoxes"),
               hierarchy,
               min_size,
               ghost_cell_width,
               balance_box_level,
               *anchor_box_level,
               anchor_to_balance);
         } else if (box_gen_method == "PrebalanceByUserShells") {
            generatePrebalanceByUserShells(
               main_db->getDatabase("PrebalanceByUserShells"),
               hierarchy,
               min_size,
               ghost_cell_width,
               balance_box_level,
               anchor_box_level,
               anchor_to_balance,
               0);
         } else {
            TBOX_ERROR("Bad box_gen_method: '" << box_gen_method << "'");
         }
      }
      hier::Connector* balance_to_anchor = &anchor_to_balance->getTranspose();

      {
         /*
          * Output "before" data.
          */
         balance_box_level->cacheGlobalReducedData();
         tbox::plog << "\n\n\nBefore:\n";
         mesh::BalanceUtilities::reduceAndReportLoadBalance(
            std::vector<double>(1, static_cast<double>(balance_box_level->getLocalNumberOfCells())),
            balance_box_level->getMPI());

         hier::BoxLevelStatistics anchor_stats(*anchor_box_level);
         tbox::plog << "Anchor box_level node stats:\n";
         anchor_stats.printBoxStats(tbox::plog, "AL-> ");
         tbox::plog << "Anchor box_level:\n";
         anchor_box_level->recursivePrint(tbox::plog, "AL-> ", 2);

         hier::BoxLevelStatistics balance_stats(*balance_box_level);
         tbox::plog << "Balance box_level node stats:\n";
         balance_stats.printBoxStats(tbox::plog, "BL-> ");
         tbox::plog << "Balance box_level:\n";
         balance_box_level->recursivePrint(tbox::plog, "BL-> ", 2);

         hier::ConnectorStatistics balance_anchor_stats(*balance_to_anchor);
         tbox::plog << "balance_to_anchor edge stats:\n";
         balance_anchor_stats.printNeighborStats(tbox::plog, "BA-> ");
         tbox::plog << "balance_to_anchor:\n";
         balance_to_anchor->recursivePrint(tbox::plog, "BA-> ");

         hier::ConnectorStatistics anchor_balance_stats(*anchor_to_balance);
         tbox::plog << "anchor_to_balance edge stats:\n";
         anchor_balance_stats.printNeighborStats(tbox::plog, "AB-> ");
         tbox::plog << "anchor_to_balance:\n";
         anchor_to_balance->recursivePrint(tbox::plog, "AB-> ");
      }

      {

         const int nnodes = mpi.getSize();
         tbox::RankGroup rank_group_1(mpi);
         if (nnodes == 1) {
            rank_group_1.setMinMax(0, 0);
         } else {
            rank_group_1.setMinMax(0, (nnodes / 2) - 1);
         }

         /*
          * Load balance the unbalanced box_level.
          */
         lb->loadBalanceBoxLevel(
            *balance_box_level,
            balance_to_anchor,
            hierarchy,
            1,
            min_size,
            max_size,
            domain_box_level,
            bad_interval,
            cut_factor,
            rank_group_1);

         balance_to_anchor->assertOverlapCorrectness();
         anchor_to_balance->assertOverlapCorrectness();

         sortNodes(*balance_box_level,
            *anchor_to_balance,
            false,
            true);
      }

      /*
       * Get the balance_to_balance for edge statistics.
       */
      oca.bridge(
         balance_to_balance,
         *balance_to_anchor,
         *anchor_to_balance,
         false);

      {
         /*
          * Output "after" data.
          */
         balance_box_level->cacheGlobalReducedData();
         tbox::plog << "\n\n\nAfter:\n";
         mesh::BalanceUtilities::reduceAndReportLoadBalance(
            std::vector<double>(1, static_cast<double>(balance_box_level->getLocalNumberOfCells())),
            balance_box_level->getMPI());

         hier::BoxLevelStatistics balance_stats(*balance_box_level);
         tbox::plog << "Balance box_level node stats:\n";
         balance_stats.printBoxStats(tbox::plog, "BL-> ");
         tbox::plog << "Balance box_level:\n";
         balance_box_level->recursivePrint(tbox::plog, "BL-> ", 2);

         hier::ConnectorStatistics balance_balance_stats(*balance_to_balance);
         tbox::plog << "balance_to_balance edge stats:\n";
         balance_balance_stats.printNeighborStats(tbox::plog, "BB-> ");
         tbox::plog << "balance_to_balance:\n";
         balance_to_balance->recursivePrint(tbox::plog, "BB-> ");

         hier::ConnectorStatistics balance_anchor_stats(*balance_to_anchor);
         tbox::plog << "balance_to_anchor edge stats:\n";
         balance_anchor_stats.printNeighborStats(tbox::plog, "BA-> ");
         tbox::plog << "balance_to_anchor:\n";
         balance_to_anchor->recursivePrint(tbox::plog, "BA-> ");

         hier::ConnectorStatistics anchor_balance_stats(*anchor_to_balance);
         tbox::plog << "anchor_to_balance edge stats:\n";
         anchor_balance_stats.printNeighborStats(tbox::plog, "AB-> ");
         tbox::plog << "anchor_to_balance:\n";
         anchor_to_balance->recursivePrint(tbox::plog, "AB-> ");

         // Dump summary statistics to output.
         mesh::BalanceUtilities::reduceAndReportLoadBalance(
            std::vector<double>(1, static_cast<double>(balance_box_level->getLocalNumberOfCells())),
            balance_box_level->getMPI(),
            tbox::plog);
      }

      hierarchy->makeNewPatchLevel(
         0,
         anchor_box_level);

      hierarchy->makeNewPatchLevel(
         1,
         balance_box_level);

#ifdef HAVE_HDF5
#if 0
      if (dimval == 2 || dimval == 3) {
         /*
          * Create the VisIt data writer.
          * Write the plot file.
          */
         DerivedVisOwnerData owner_writer;
         const std::string visit_filename = base_name + ".visit";
         appu::VisItDataWriter visit_data_writer(dim,
                                                 "VisIt Writer",
                                                 visit_filename);
         visit_data_writer.registerDerivedPlotQuantity("Owner",
            "SCALAR",
            &owner_writer);
         visit_data_writer.writePlotData(hierarchy, 0);
      }
#endif
#endif

   }

   /*
    * Print input database again to fully show usage.
    */
   plog << "Input database after running..." << std::endl;
   tbox::InputManager::getManager()->getInputDatabase()->printClassData(plog);

   tbox::pout << "\nPASSED:  rank_group" << std::endl;

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

   return fail_count;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void generatePrebalanceByUserShells(
   const std::shared_ptr<tbox::Database>& database,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const hier::IntVector& min_size,
   const hier::IntVector& max_gcw,
   std::shared_ptr<hier::BoxLevel>& balance_box_level,
   const std::shared_ptr<hier::BoxLevel>& anchor_box_level,
   std::shared_ptr<hier::Connector>& anchor_to_balance,
   int tag_level_number)
{

   const tbox::Dimension dim(hierarchy->getDim());
   const int dimval = dim.getValue();

   /*
    * Starting at shell origin, tag cells with centroids
    * at radii[0]<r<radii[1], radii[2]<r<radii[3], and so on.
    */
   std::vector<double> radii;

   std::vector<double> r0(dimval);
   for (int d = 0; d < dimval; ++d) r0[d] = 0;

   std::shared_ptr<tbox::Database> abr_db;
   if (database) {
      if (database->isDouble("r0")) {
         r0 = database->getDoubleVector("r0");
      }
      if (database->isDouble("radii")) {
         radii = database->getDoubleVector("radii");
      }
      abr_db = database->getDatabaseWithDefault("BergerRigoutsos", abr_db);
      TBOX_ASSERT(static_cast<int>(radii.size()) % 2 == 0);
   }

   const int tag_val = 1;

   hier::VariableDatabase* vdb =
      hier::VariableDatabase::getDatabase();
   std::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
         hierarchy->getGridGeometry()));
   TBOX_ASSERT(grid_geometry);

   std::shared_ptr<hier::PatchLevel> tag_level(
      new hier::PatchLevel(
         anchor_box_level,
         grid_geometry,
         vdb->getPatchDescriptor()));
   tag_level->setLevelNumber(tag_level_number);

   std::shared_ptr<pdat::CellVariable<int> > tag_variable(
      new pdat::CellVariable<int>(dim, "TagVariable"));

   std::shared_ptr<hier::VariableContext> default_context(
      vdb->getContext("TagVariable"));

   const int tag_id = vdb->registerVariableAndContext(
         tag_variable,
         default_context,
         hier::IntVector(dim, 0));

   tag_level->allocatePatchData(tag_id);

   const double* xlo = grid_geometry->getXLower();
   const double* h = grid_geometry->getDx();
   for (hier::PatchLevel::iterator pi(tag_level->begin());
        pi != tag_level->end(); ++pi) {
      const std::shared_ptr<hier::Patch>& patch = *pi;
      std::shared_ptr<pdat::CellData<int> > tag_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch->getPatchData(tag_id)));
      TBOX_ASSERT(tag_data);

      tag_data->getArrayData().undefineData();

      pdat::CellData<int>::iterator ciend(pdat::CellGeometry::end(tag_data->getGhostBox()));
      for (pdat::CellData<int>::iterator ci(pdat::CellGeometry::begin(tag_data->getGhostBox()));
           ci != ciend; ++ci) {
         const pdat::CellIndex& idx = *ci;
         double rr = 0;
         std::vector<double> r(dimval);
         for (int d = 0; d < dimval; ++d) {
            r[d] = xlo[d] + h[d] * (idx(d) + 0.5) - r0[d];
            rr += r[d] * r[d];
         }
         rr = sqrt(rr);
         for (int i = 0; i < static_cast<int>(radii.size()); i += 2) {
            if (radii[i] < rr && rr < radii[i + 1]) {
               (*tag_data)(idx) = tag_val;
               break;
            }
         }
      }
   }

   mesh::BergerRigoutsos abr(dim, abr_db);
   abr.useDuplicateMPI(anchor_box_level->getMPI());
   abr.findBoxesContainingTags(
      balance_box_level,
      anchor_to_balance,
      tag_level,
      tag_id,
      tag_val,
      hier::BoxContainer(anchor_box_level->getGlobalBoundingBox(hier::BlockId(0))),
      min_size,
      max_gcw);

   hier::Connector& balance_to_anchor = anchor_to_balance->getTranspose();

   /*
    * The clustering step generated Connectors to/from the temporary
    * tag_level->getBoxLevel(), which is not the same as the
    * anchor BoxLevel.  We need to reset the Connectors to use
    * the anchor_box_level instead.
    */
   anchor_to_balance->setBase(*anchor_box_level);
   anchor_to_balance->setHead(*balance_box_level, true);
   balance_to_anchor.setBase(*balance_box_level);
   balance_to_anchor.setHead(*anchor_box_level, true);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void generatePrebalanceByUserBoxes(
   std::shared_ptr<tbox::Database> database,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const hier::IntVector& min_size,
   const hier::IntVector& max_gcw,
   std::shared_ptr<hier::BoxLevel>& balance_box_level,
   const hier::BoxLevel& anchor_box_level,
   std::shared_ptr<hier::Connector>& anchor_to_balance)
{
   NULL_USE(hierarchy);
   NULL_USE(min_size);

   const tbox::Dimension& dim(hierarchy->getDim());

   std::vector<tbox::DatabaseBox> db_box_vector =
      database->getDatabaseBoxVector("balance_boxes");
   hier::BoxContainer balance_boxes(db_box_vector);
   std::vector<int> initial_owners(1);
   initial_owners[0] = 0;
   initial_owners = database->getIntegerVector("initial_owners");

   balance_box_level.reset(new hier::BoxLevel(hier::IntVector(dim, 1),
         hierarchy->getGridGeometry(),
         anchor_box_level.getMPI()));
   hier::BoxContainer::iterator balance_boxes_itr = balance_boxes.begin();
   for (int i = 0; i < balance_boxes.size(); ++i, ++balance_boxes_itr) {
      const int owner = i % static_cast<int>(initial_owners.size());
      if (owner == balance_box_level->getMPI().getRank()) {
         balance_boxes_itr->setBlockId(hier::BlockId(0));
         balance_box_level->addBox(hier::Box(*balance_boxes_itr,
               hier::LocalId(i), owner));
      }
   }
   hier::OverlapConnectorAlgorithm oca;
   oca.findOverlapsWithTranspose(anchor_to_balance,
      anchor_box_level,
      *balance_box_level,
      max_gcw,
      max_gcw);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void sortNodes(
   hier::BoxLevel& new_box_level,
   hier::Connector& tag_to_new,
   bool sort_by_corners,
   bool sequentialize_global_indices)
{
   const hier::MappingConnectorAlgorithm mca;

   std::shared_ptr<hier::MappingConnector> sorting_map;
   std::shared_ptr<hier::BoxLevel> seq_box_level;
   hier::BoxLevelConnectorUtils dlbg_edge_utils;
   dlbg_edge_utils.makeSortingMap(
      seq_box_level,
      sorting_map,
      new_box_level,
      sort_by_corners,
      sequentialize_global_indices);

   mca.modify(tag_to_new,
      *sorting_map,
      &new_box_level);
}
