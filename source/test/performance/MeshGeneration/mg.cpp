/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test for performance and quality of mesh generation.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include <iomanip>

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/BoxLevelStatistics.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/ConnectorStatistics.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/mesh/CascadePartitioner.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/mesh/TileClustering.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/appu/VisItDataWriter.h"

#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/tbox/BreadthFirstRankTree.h"
#include "SAMRAI/tbox/CenteredRankTree.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/OpenMPUtilities.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include <vector>

#include "test/testlib/SinusoidalFrontGenerator.h"
#include "test/testlib/SphericalShellGenerator.h"
#include "test/testlib/ShrunkenLevelGenerator.h"

using namespace SAMRAI;
using namespace tbox;

/*
 ************************************************************************
 *
 *
 *************************************************************************
 */

void
enforceNesting(
   hier::BoxLevel& L1,
   hier::Connector& L0_to_L1,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int coarser_ln);

void
sortNodes(
   hier::BoxLevel& new_box_level,
   hier::Connector& tag_to_new,
   bool sort_by_corners,
   bool sequentialize_global_indices);

void
refineHead(
   hier::BoxLevel& head,
   hier::Connector& ref_to_head,
   const hier::IntVector& refinement_ratio);

void
outputPostcluster(
   const hier::BoxLevel& cluster,
   const hier::BoxLevel& ref,
   const hier::IntVector& ref_to_cluster_width,
   const std::string& border);

void
outputPrebalance(
   const hier::BoxLevel& pre,
   const hier::BoxLevel& ref,
   const hier::IntVector& pre_width,
   const std::string& border);

void
outputPostbalance(
   const hier::BoxLevel& post,
   const hier::BoxLevel& ref,
   const hier::IntVector& post_width,
   const std::string& border);

std::shared_ptr<mesh::BoxGeneratorStrategy>
createBoxGenerator(
   const std::shared_ptr<tbox::Database>& input_db,
   const std::string& bg_type,
   int ln,
   const tbox::Dimension& dim);

std::shared_ptr<mesh::LoadBalanceStrategy>
createLoadBalancer(
   const std::shared_ptr<tbox::Database>& input_db,
   const std::string& lb_type,
   const std::string& rank_tree_type,
   int ln,
   const tbox::Dimension& dim);

std::shared_ptr<RankTreeStrategy>
getRankTree(
   Database& input_db,
   const std::string& rank_tree_type);

/*!
 * @brief Implementation to tell PatchHierarchy about the request
 * for Connector widths used in enforcing nesting.
 *
 * This is not essential, but we chose to go through the hierarchy to
 * determine how big a Connector width to compute during the level
 * generation.  This more closely resembles what real aplications do.
 * This step is typically done in the mesh generator, and what we are
 * writing here is essentially a mesh generator.
 */
class NestingLevelConnectorWidthRequestor:
   public hier::PatchHierarchy::ConnectorWidthRequestorStrategy
{
public:
   virtual void
   computeRequiredConnectorWidths(
      std::vector<hier::IntVector>& self_connector_widths,
      std::vector<hier::IntVector>& fine_connector_widths,
      const hier::PatchHierarchy& patch_hierarchy) const
   {
      self_connector_widths.clear();
      self_connector_widths.reserve(patch_hierarchy.getMaxNumberOfLevels());
      const hier::IntVector& one = hier::IntVector::getOne((patch_hierarchy.getDim()));
      for (int ln = 0; ln < patch_hierarchy.getMaxNumberOfLevels(); ++ln) {
         self_connector_widths.push_back(
            one * patch_hierarchy.getProperNestingBuffer(ln));
      }
      // fine_connector_widths is same, but doesn't need last level's.
      fine_connector_widths = self_connector_widths;
      fine_connector_widths.pop_back();
   }
};
NestingLevelConnectorWidthRequestor nesting_level_connector_width_requestor;

static std::shared_ptr<tbox::CommGraphWriter> comm_graph_writer;
size_t num_records_written = 0;

/*
 ********************************************************************************
 *
 * Performance testing for mesh generation operations.
 *
 * 1. Build "level 0" from the domain description (input parameter
 * "domain_boxes").  L0 is for doing the test, not for checking load
 * balancer performance.
 *
 * 2. Build "level 1" and write out performance data for balancing it.
 *
 * 3. Build "level 2" and write out performance data for balancing it.
 * The prebalance boxes for L2 are generated by clustering tags on L1.
 * All L1 cells are tagged except for a small margin by the L1 boundary
 * (input parameter "tag_margin".  This configuration tries to mimick
 * real problems where the tags occupy a large portion of the tag
 * level, leading to a greater number of owners for prebalance boxes.
 *
 ********************************************************************************
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

   int error_count = 0;

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
      TimerManager::createManager(input_db->getDatabase("TimerManager"));
      std::shared_ptr<tbox::Timer> t_all =
         tbox::TimerManager::getManager()->getTimer("appu::main::all");
      t_all->start();

      /*
       * Retrieve "Main" section from input database.
       * The main database is used only in main().
       * The base_name variable is a base name for
       * all name strings in this program.
       */

      std::shared_ptr<Database> main_db = input_db->getDatabase("Main");

      const tbox::Dimension
      dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      const hier::IntVector& zero_vec = hier::IntVector::getZero(dim);

      std::string base_name = "unnamed";
      base_name = main_db->getStringWithDefault("base_name", base_name);

      /*
       * Modify basename for this particular run.
       * Add the number of processes and the case name.
       */
      std::string base_name_ext = base_name;
      if (!case_name.empty()) {
         base_name_ext = base_name_ext + '-' + case_name;
      }
      base_name_ext = base_name_ext + '-'
         + tbox::Utilities::nodeToString(mpi.getSize());
      tbox::plog << "Added case name (" << case_name << ") and nprocs ("
                 << mpi.getSize() << ") to base name -> '"
                 << base_name_ext << "'\n";

      if (!case_name.empty()) {
         tbox::plog << "Added case name (" << case_name << ") and nprocs ("
                    << mpi.getSize() << ") to base name -> '"
                    << base_name_ext << "'\n";
      }

      /*
       * Start logging.
       */
      const std::string log_file_name = base_name_ext + ".log";
      bool log_all_nodes = false;
      log_all_nodes = main_db->getBoolWithDefault("log_all_nodes",
            log_all_nodes);
      if (log_all_nodes) {
         PIO::logAllNodes(log_file_name);
      } else {
         PIO::logOnlyNodeZero(log_file_name);
      }

      tbox::plog << "MPI has " << tbox::SAMRAI_MPI::getSAMRAIWorld().getSize()
                 << " processes." << std::endl;
      tbox::plog << "OpenMP version "
                 << TBOX_omp_version << ".\n"
                 << "Running with " << TBOX_omp_get_max_threads() << " threads."
                 << std::endl;

      /*
       * Whether to perform certain steps in mesh generation.
       */

      std::vector<bool> enforce_nesting(1, true);
      if (main_db->isBool("enforce_nesting")) {
         enforce_nesting = main_db->getBoolVector("enforce_nesting");
      }

      std::vector<bool> load_balance(1, true);
      if (main_db->isBool("load_balance")) {
         load_balance = main_db->getBoolVector("load_balance");
      }

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

      std::vector<double> xlo(dim.getValue());
      std::vector<double> xhi(dim.getValue());
      for (int i = 0; i < dim.getValue(); ++i) {
         xlo[i] = 0.0;
         xhi[i] = 1.0;
      }
      if (main_db->isDouble("xlo")) {
         xlo = main_db->getDoubleVector("xlo");
      }
      if (main_db->isDouble("xhi")) {
         xhi = main_db->getDoubleVector("xhi");
      }

      /*
       * Choose the tagging code.
       */
      const std::string mesh_generator_name =
         main_db->getStringWithDefault("mesh_generator_name", "SinusoidalFrontGenerator");
      std::shared_ptr<MeshGenerationStrategy> mesh_gen;
      if (mesh_generator_name == "SinusoidalFrontGenerator") {
         mesh_gen.reset(
            new SinusoidalFrontGenerator(
               "SinusoidalFrontGenerator",
               dim,
               main_db->getDatabaseWithDefault("SinusoidalFrontGenerator",
                  std::shared_ptr<tbox::Database>())));
      } else if (mesh_generator_name == "SphericalShellGenerator") {
         mesh_gen.reset(
            new SphericalShellGenerator(
               "SphericalShellGenerator",
               dim,
               main_db->getDatabaseWithDefault("SphericalShellGenerator",
                  std::shared_ptr<tbox::Database>())));
      } else if (mesh_generator_name == "ShrunkenLevelGenerator") {
         mesh_gen.reset(
            new ShrunkenLevelGenerator(
               "ShrunkenLevelGenerator",
               dim,
               main_db->getDatabaseWithDefault("ShrunkenLevelGenerator",
                  std::shared_ptr<tbox::Database>())));
      } else {
         TBOX_ERROR("Unrecognized MeshGeneratorStrategy " << mesh_generator_name);
      }

      /*
       * If autoscale_base_nprocs is given, take the domain_boxes, xlo and xhi
       * to be the size for the (integer) value of autoscale_base_nprocs.  Scale
       * the problem from there to the number of process running by
       * doubling the size starting with the j direction.
       *
       * The number of processes must be a power of 2 times the value
       * of autoscale_base_nprocs.
       */
      const int autoscale_base_nprocs =
         main_db->getIntegerWithDefault("autoscale_base_nprocs", mpi.getSize());

      mesh_gen->setDomain(domain_boxes, &xlo[0], &xhi[0], autoscale_base_nprocs, mpi);

      hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();

      /*
       * Clustering algorithm.
       */

      std::string box_generator_type =
         main_db->getStringWithDefault("box_generator_type", "BergerRigoutsos");

      /*
       * Create hierarchy.
       */

      tbox::plog << "Building domain with boxes " << domain_boxes.format() << std::endl;
      std::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
         new geom::CartesianGridGeometry(
            "GridGeometry",
            &xlo[0],
            &xhi[0],
            domain_boxes));

      std::shared_ptr<hier::PatchHierarchy> hierarchy(
         new hier::PatchHierarchy(
            "Hierarchy",
            grid_geometry,
            input_db->getDatabase("PatchHierarchy")));

      hierarchy->registerConnectorWidthRequestor(nesting_level_connector_width_requestor);

      mesh_gen->resetHierarchyConfiguration(hierarchy, 0, 1);

      enforce_nesting.resize(hierarchy->getMaxNumberOfLevels(),
         bool(enforce_nesting.back()));

      load_balance.resize(hierarchy->getMaxNumberOfLevels(),
         bool(load_balance.back()));

      const int max_levels = hierarchy->getMaxNumberOfLevels();

      hier::BoxLevel domain_box_level(
         hier::IntVector(dim, 1),
         grid_geometry,
         tbox::SAMRAI_MPI::getSAMRAIWorld(),
         hier::BoxLevel::GLOBALIZED);
      for (hier::BoxContainer::const_iterator bi =
              grid_geometry->getPhysicalDomain().begin();
           bi != grid_geometry->getPhysicalDomain().end(); ++bi) {
         domain_box_level.addBox(*bi);
      }

      /*
       * Set up the patch data for tags.
       */

      std::shared_ptr<pdat::CellVariable<int> > tag_variable(
         new pdat::CellVariable<int>(dim, "ShrinkingLevelTagVariable"));

      std::shared_ptr<hier::VariableContext> default_context =
         vdb->getContext("TagVariable");

      const int tag_data_id = vdb->registerVariableAndContext(
            tag_variable,
            default_context,
            hier::IntVector::getZero(dim));

      /*
       * Set up the load balancers.
       */

      std::vector<std::string> load_balancer_type =
         main_db->getStringVector("load_balancer_type");
      load_balancer_type.reserve(hierarchy->getMaxNumberOfLevels());
      while (static_cast<int>(load_balancer_type.size()) < hierarchy->getMaxNumberOfLevels()) {
         load_balancer_type.push_back(load_balancer_type.back());
      }

      std::string rank_tree_type =
         main_db->getStringWithDefault("rank_tree_type", "CenteredRankTree");

      const bool write_comm_graph = main_db->getBoolWithDefault("write_comm_graph", false);
      if (write_comm_graph) {
         comm_graph_writer.reset(new CommGraphWriter);
         comm_graph_writer->setWriteFullGraph(mpi.getSize() <= 16);
      }

      plog << "Input database after initialization..." << std::endl;
      input_db->printClassData(plog);

      /*
       * Build L0.
       */
      {
         tbox::pout << "\n==================== Generating L0 ====================" << std::endl;

         std::shared_ptr<hier::BoxLevel> L0(
            std::make_shared<hier::BoxLevel>(
               hier::IntVector(dim, 1), grid_geometry));

         hier::BoxContainer L0_boxes(
            grid_geometry->getPhysicalDomain());
         const int boxes_per_proc =
            (L0_boxes.size() + L0->getMPI().getSize() - 1)
            / L0->getMPI().getSize();
         const int my_boxes_start = L0->getMPI().getRank() * boxes_per_proc;
         const int my_boxes_stop =
            tbox::MathUtilities<int>::Min(my_boxes_start + boxes_per_proc,
               L0_boxes.size());
         hier::BoxContainer::iterator L0_boxes_itr = L0_boxes.begin();
         for (int i = 0; i < my_boxes_start; ++i) {
            if (L0_boxes_itr != L0_boxes.end()) {
               ++L0_boxes_itr;
            }
         }
         for (int i = my_boxes_start; i < my_boxes_stop; ++i, ++L0_boxes_itr) {
            L0->addBox(*L0_boxes_itr, hier::BlockId::zero());
         }

         /*
          * Load balance the L0 BoxLevel, using the domain as its L0.
          *
          * This is not a part of the performance test because does not
          * reflect the load balancer use in real apps.  We just neeed a
          * distributed L0 for the real load balancing performance test.
          */
         std::shared_ptr<hier::Connector> domain_to_L0;
         oca.findOverlapsWithTranspose(domain_to_L0,
            domain_box_level,
            *L0,
            hier::IntVector(dim, 2),
            hier::IntVector(dim, 2));
         hier::Connector* L0_to_domain = &domain_to_L0->getTranspose();

         std::shared_ptr<mesh::LoadBalanceStrategy> lb =
            createLoadBalancer(input_db, load_balancer_type[0], rank_tree_type, 0, dim);

         tbox::plog << "\n\tL0 prebalance loads:\n";
         mesh::BalanceUtilities::reduceAndReportLoadBalance(
            std::vector<double>(1, static_cast<double>(L0->getLocalNumberOfCells())),
            L0->getMPI());

         outputPrebalance(*L0, domain_box_level, hierarchy->getRequiredConnectorWidth(0, 0), "L0: ");

         if (load_balance[0]) {
            tbox::pout << "\tPartitioning..." << std::endl;
            tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier();
            lb->loadBalanceBoxLevel(
               *L0,
               L0_to_domain,
               hierarchy,
               0,
               hierarchy->getSmallestPatchSize(0),
               hierarchy->getLargestPatchSize(0),
               domain_box_level,
               hier::IntVector::getOne(dim),
               hier::IntVector::getOne(dim));
         }

         sortNodes(*L0,
            *domain_to_L0,
            false,
            true);

         tbox::plog << "\n\tL0 postbalance loads:\n";
         mesh::BalanceUtilities::reduceAndReportLoadBalance(
            std::vector<double>(1, static_cast<double>(L0->getLocalNumberOfCells())),
            L0->getMPI());

         outputPostbalance(*L0, domain_box_level, hierarchy->getRequiredConnectorWidth(0,
               0), "L0: ");

         if (comm_graph_writer) {
            tbox::pout << "\tWriting communication graph..." << std::endl;
            tbox::plog << "\nCommunication Graph for balancing L0:\n";
            for ( ;
                  num_records_written < comm_graph_writer->getNumberOfRecords();
                  ++num_records_written) {
               comm_graph_writer->writeGraphToTextStream(num_records_written, tbox::plog);
            }
            tbox::plog << "\n";
         }

         L0->cacheGlobalReducedData();

         hierarchy->makeNewPatchLevel(0, L0);
      }

      for (int new_ln = 1; new_ln < max_levels; ++new_ln) {
         /*
          * Build level new_ln.
          */
         tbox::pout << "\n==================== Generating L" << new_ln
                    << " ====================" << std::endl;

         const int tag_ln = new_ln - 1;

         const hier::BoxLevel& Ltag = *hierarchy->getPatchLevel(tag_ln)->getBoxLevel();

         std::shared_ptr<hier::BoxLevel> Lnew;
         std::shared_ptr<hier::Connector> Ltag_to_Lnew;
         std::shared_ptr<hier::Connector> Lnew_to_Lnew;

         // Get the prebalanced Lnew:
         const hier::IntVector required_connector_width =
            hierarchy->getRequiredConnectorWidth(tag_ln, new_ln);
         const hier::IntVector min_size = hier::IntVector::ceilingDivide(
               hierarchy->getSmallestPatchSize(new_ln),
               hierarchy->getRatioToCoarserLevel(new_ln));

         /*
          * Tag cells.
          */
         tbox::pout << "\tTagging..." << std::endl;
         bool exact_tagging = false;
         hierarchy->getPatchLevel(tag_ln)->allocatePatchData(tag_data_id);
         mesh_gen->setTags(exact_tagging, hierarchy, tag_ln, tag_data_id);

         /*
          * Cluster.
          */
         tbox::pout << "\tClustering..." << std::endl;
         std::shared_ptr<mesh::BoxGeneratorStrategy> bg =
            createBoxGenerator(input_db, box_generator_type, new_ln, dim);
         tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier();
         bg->findBoxesContainingTags(
            Lnew,
            Ltag_to_Lnew,
            hierarchy->getPatchLevel(tag_ln),
            tag_data_id,
            1 /* tag_val */,
            hier::BoxContainer(Ltag.getGlobalBoundingBox(hier::BlockId(0))),
            min_size,
            required_connector_width);

         if (Ltag_to_Lnew->getConnectorWidth() != required_connector_width) {
            const hier::Connector& Ltag_to_Ltag =
               Ltag.findConnectorWithTranspose(
                  Ltag,
                  required_connector_width,
                  required_connector_width,
                  hier::CONNECTOR_IMPLICIT_CREATION_RULE);
            hier::OverlapConnectorAlgorithm timed_oca;
            timed_oca.setTimerPrefix(std::string(
                  "apps::fix_zero_width") + tbox::Utilities::intToString(new_ln));
            tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier();
            timed_oca.bridgeWithNesting(
               Ltag_to_Lnew,
               Ltag_to_Ltag,
               hier::Connector(*Ltag_to_Lnew),
               hier::IntVector::getZero(dim),
               hier::IntVector::getZero(dim),
               required_connector_width,
               true);
         }

         outputPostcluster(*Lnew, Ltag, required_connector_width, "Lnew: ");

         if (Lnew->getGlobalNumberOfBoxes() == 0) {
            TBOX_WARNING("Level " << new_ln << " box generator resulted in no boxes.  Stopping.");
            break;
         }

         /*
          * Enforce nesting.
          */
         if (enforce_nesting[new_ln]) {
            enforceNesting(
               *Lnew,
               *Ltag_to_Lnew,
               hierarchy,
               tag_ln);

            if (Lnew->getGlobalNumberOfBoxes() == 0) {
               TBOX_WARNING(
                  "Level " << new_ln << " box generator has no box after proper nesting.");
            }
         }

         if (hierarchy->getRatioToCoarserLevel(new_ln) != zero_vec) {
            refineHead(
               *Lnew,
               *Ltag_to_Lnew,
               hierarchy->getRatioToCoarserLevel(new_ln));
         }

         std::shared_ptr<mesh::LoadBalanceStrategy> lb =
            createLoadBalancer(input_db, load_balancer_type[new_ln], rank_tree_type, new_ln, dim);

         outputPrebalance(*Lnew, Ltag, required_connector_width, "Lnew: ");

         tbox::plog << "\n\tLnew prebalance loads:\n";
         mesh::BalanceUtilities::reduceAndReportLoadBalance(
            std::vector<double>(static_cast<size_t>(1),
               static_cast<double>(Lnew->getLocalNumberOfCells())),
            Lnew->getMPI());

         if (load_balance[new_ln]) {
            tbox::pout << "\tPartitioning..." << std::endl;
            tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier();
            lb->loadBalanceBoxLevel(
               *Lnew,
               &Ltag_to_Lnew->getTranspose(),
               hierarchy,
               tag_ln,
               hierarchy->getSmallestPatchSize(new_ln),
               hierarchy->getLargestPatchSize(new_ln),
               domain_box_level,
               hierarchy->getRatioToCoarserLevel(new_ln),
               hierarchy->getRatioToCoarserLevel(new_ln));
         }

         sortNodes(*Lnew,
            *Ltag_to_Lnew,
            false,
            true);

         outputPostbalance(*Lnew, Ltag, required_connector_width, "Lnew: ");

         tbox::plog << "\n\tLnew postbalance loads:\n";
         mesh::BalanceUtilities::reduceAndReportLoadBalance(
            std::vector<double>(static_cast<size_t>(1),
               static_cast<double>(Lnew->getLocalNumberOfCells())),
            Lnew->getMPI());

         if (comm_graph_writer) {
            tbox::pout << "\tWriting communication graph..." << std::endl;
            tbox::plog << "\nCommunication Graph for balancing Lnew:\n";
            for ( ;
                  num_records_written < comm_graph_writer->getNumberOfRecords();
                  ++num_records_written) {
               comm_graph_writer->writeGraphToTextStream(num_records_written, tbox::plog);
            }
            tbox::plog << "\n\n";
         }

         // Get the Lnew_to_Lnew for edge statistics.
         oca.bridge(
            Lnew_to_Lnew,
            Ltag_to_Lnew->getTranspose(),
            *Ltag_to_Lnew,
            false);

         hierarchy->makeNewPatchLevel(new_ln, Lnew);

      } // end new_ln loop

      tbox::plog << "\n==================== Final hierarchy ====================" << std::endl;
      for (int ln = 0; ln < hierarchy->getNumberOfLevels(); ++ln) {
         tbox::plog << '\n'
                    << "\tL" << ln << " summary:\n"
                    << hierarchy->getPatchLevel(ln)->getBoxLevel()->format("\t\t", 0)
                    << "\tL" << ln << " statistics:\n"
                    << hierarchy->getPatchLevel(ln)->getBoxLevel()->formatStatistics("\t\t");
      }
      tbox::plog << "\n\n";

      bool write_visit =
         main_db->getBoolWithDefault("write_visit", false);
      if (write_visit) {
#ifdef HAVE_HDF5

         if ((dim == tbox::Dimension(2)) || (dim == tbox::Dimension(3))) {
            /*
             * Create the VisIt data writer.
             * Write the plot file.
             */
            const std::string visit_filename = base_name_ext + ".visit";
            appu::VisItDataWriter visit_data_writer(dim,
                                                    "VisIt Writer",
                                                    visit_filename);
            mesh_gen->registerVariablesWithPlotter(visit_data_writer);
            visit_data_writer.writePlotData(hierarchy, 0);
         }
#else
         TBOX_WARNING("main: You set write_visit to TRUE,\n"
            << "but VisIt dumps are not supported due to\n"
            << "not having configured with HDF5.\n");
#endif
      }

      t_all->stop();
      int size = tbox::SAMRAI_MPI::getSAMRAIWorld().getSize();
      if (tbox::SAMRAI_MPI::getSAMRAIWorld().getRank() == 0) {
         std::string timing_file =
            base_name + ".timing" + tbox::Utilities::intToString(size);
         FILE* fp = fopen(timing_file.c_str(), "w");
         fprintf(fp, "%f\n", t_all->getTotalWallclockTime());
         fclose(fp);
      }

   }

   /*
    * Output timer results.
    */
   tbox::TimerManager::getManager()->print(tbox::plog);

   /*
    * Print input database again to fully show usage.
    */
   plog << "Input database after running..." << std::endl;
   tbox::InputManager::getManager()->getInputDatabase()->printClassData(plog);

   if (error_count == 0) {
      tbox::pout << "\nPASSED:  MeshGeneration" << std::endl;
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

   return error_count;
}

/*
 ****************************************************************************
 * Output post-cluster metadata.
 ****************************************************************************
 */
void outputPostcluster(
   const hier::BoxLevel& cluster,
   const hier::BoxLevel& ref,
   const hier::IntVector& ref_to_cluster_width,
   const std::string& border)
{
   cluster.cacheGlobalReducedData();

   const hier::IntVector cluster_to_ref_width =
      hier::Connector::convertHeadWidthToBase(
         cluster.getRefinementRatio(),
         ref.getRefinementRatio(),
         ref_to_cluster_width);

   const hier::Connector& cluster_to_ref = cluster.findConnector(ref,
         ref_to_cluster_width,
         hier::CONNECTOR_CREATE,
         true);

   cluster.cacheGlobalReducedData();
   cluster_to_ref.cacheGlobalReducedData();

   tbox::plog << "\n\n"
              << border << "Cluster summary:\n"
              << cluster.format(border + "\t", 0)
              << border << "Cluster stats:\n"
              << cluster.formatStatistics(border + "\t");

   tbox::plog << '\n'
              << border << "cluster--->ref summary:\n"
              << cluster_to_ref.format(border + "\t", 0)
              << border << "cluster--->ref stats:\n"
              << cluster_to_ref.formatStatistics(border + "\t");
}

/*
 ****************************************************************************
 * Output pre-balance metadata.
 * - pre-balance level
 * - pre--->pre for proximity contrast to post--->post
 * - pre--->ref for proximity contrast to post--->pre
 * - ref--->pre for proximity contrast to ref--->post
 ****************************************************************************
 */
void outputPrebalance(
   const hier::BoxLevel& pre,
   const hier::BoxLevel& ref,
   const hier::IntVector& pre_width,
   const std::string& border)
{
   pre.cacheGlobalReducedData();

   const hier::IntVector ref_width =
      hier::Connector::convertHeadWidthToBase(
         ref.getRefinementRatio(),
         pre.getRefinementRatio(),
         pre_width);

   const hier::Connector& pre_to_pre = pre.findConnector(pre,
         pre_width,
         hier::CONNECTOR_CREATE,
         true);

   tbox::plog << "\n\n"
              << border << "Prebalance summary:\n"
              << pre.format(border + "\t", 0)
              << border << "Prebalance stats:\n"
              << pre.formatStatistics(border + "\t");

   tbox::plog << '\n'
              << border << "pre--->pre summary:\n"
              << pre_to_pre.format(border + "\t", 0)
              << border << "pre--->pre stats:\n"
              << pre_to_pre.formatStatistics(border + "\t");
}

/*
 ****************************************************************************
 * Output post-balance metadata:
 * - post-balance level
 * - post--->post for proximity evaluation
 * - post--->ref for proximity contrast to pre--->ref
 * - ref--->post for proximity contrast to ref--->pre
 ****************************************************************************
 */
void outputPostbalance(
   const hier::BoxLevel& post,
   const hier::BoxLevel& ref,
   const hier::IntVector& post_width,
   const std::string& border)
{
   post.cacheGlobalReducedData();

   const hier::IntVector ref_width =
      hier::Connector::convertHeadWidthToBase(
         ref.getRefinementRatio(),
         post.getRefinementRatio(),
         post_width);

   const hier::Connector& post_to_post = post.findConnector(post,
         post_width,
         hier::CONNECTOR_CREATE,
         true);

   const hier::Connector& post_to_ref = post.findConnectorWithTranspose(ref,
         (post.getRefinementRatio() / ref.getRefinementRatio()) * ref_width,
         ref_width,
         hier::CONNECTOR_CREATE,
         true);

   const hier::Connector& ref_to_post = post_to_ref.getTranspose();

   tbox::plog << "\n\n"
              << border << "Postbalance summary:\n"
              << post.format(border + "\t", 0)
              << border << "Postbalance stats:\n"
              << post.formatStatistics(border + "\t");

   tbox::plog << '\n'
              << border << "post--->post summary:\n"
              << post_to_post.format(border + "\t", 0)
              << border << "post--->post stats:\n"
              << post_to_post.formatStatistics(border + "\t");

   tbox::plog << '\n'
              << border << "post--->ref summary:\n"
              << post_to_ref.format(border + "\t", 0)
              << border << "post--->ref stats:\n"
              << post_to_ref.formatStatistics(border + "\t");

   tbox::plog << '\n'
              << border << "ref--->post summary:\n"
              << ref_to_post.format(border + "\t", 0)
              << border << "ref--->post stats:\n"
              << ref_to_post.formatStatistics(border + "\t");
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
   hier::BoxLevelConnectorUtils dlbg_edge_utils;
   dlbg_edge_utils.setTimerPrefix("apps::sortNodes");

   std::shared_ptr<hier::MappingConnector> sorting_map;
   std::shared_ptr<hier::BoxLevel> seq_box_level;
   dlbg_edge_utils.makeSortingMap(
      seq_box_level,
      sorting_map,
      new_box_level,
      sort_by_corners,
      sequentialize_global_indices);

   hier::MappingConnectorAlgorithm mca;
   mca.setTimerPrefix("apps::sortNodes");

   mca.modify(tag_to_new,
      *sorting_map,
      &new_box_level);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void refineHead(
   hier::BoxLevel& head,
   hier::Connector& ref_to_head,
   const hier::IntVector& refinement_ratio)
{
   hier::Connector& head_to_ref = ref_to_head.getTranspose();

   head.refineBoxes(
      head,
      refinement_ratio,
      head.getRefinementRatio() * refinement_ratio);
   head.finalize();

   hier::IntVector head_to_ref_width =
      refinement_ratio * head_to_ref.getConnectorWidth();
   head_to_ref.setBase(head);
   head_to_ref.setWidth(head_to_ref_width, true);

   ref_to_head.setHead(head, true);
   ref_to_head.refineLocalNeighbors(refinement_ratio);
}

std::shared_ptr<mesh::LoadBalanceStrategy>
createLoadBalancer(
   const std::shared_ptr<tbox::Database>& input_db,
   const std::string& lb_type,
   const std::string& rank_tree_type,
   int ln,
   const tbox::Dimension& dim)
{

   if (lb_type == "TreeLoadBalancer") {

      std::shared_ptr<tbox::RankTreeStrategy> rank_tree = getRankTree(*input_db,
            rank_tree_type);

      std::shared_ptr<mesh::TreeLoadBalancer>
      tree_lb(new mesh::TreeLoadBalancer(
                 dim,
                 std::string("mesh::TreeLoadBalancer") + tbox::Utilities::intToString(ln),
                 input_db->getDatabaseWithDefault("TreeLoadBalancer",
                    std::shared_ptr<tbox::Database>())));
      tree_lb->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());
      tree_lb->setCommGraphWriter(comm_graph_writer);
      return tree_lb;

   } else if (lb_type == "CascadePartitioner") {

      std::shared_ptr<mesh::CascadePartitioner>
      cascade_lb(new mesh::CascadePartitioner(
                    dim,
                    std::string("mesh::CascadePartitioner") + tbox::Utilities::intToString(ln),
                    input_db->getDatabaseWithDefault("CascadePartitioner",
                       std::shared_ptr<tbox::Database>())));
      return cascade_lb;

   } else if (lb_type == "ChopAndPackLoadBalancer") {

      std::shared_ptr<mesh::ChopAndPackLoadBalancer>
      cap_lb(new mesh::ChopAndPackLoadBalancer(
                dim,
                std::string("mesh::ChopAndPackLoadBalancer") + tbox::Utilities::intToString(ln),
                input_db->getDatabaseWithDefault("ChopAndPackLoadBalancer",
                   std::shared_ptr<tbox::Database>())));
      return cap_lb;

   } else {
      TBOX_ERROR(
         "Missing or bad load_balancer specification in Main database.\n"
         << "Specify load_balancer_type = STRING, where STRING can be\n"
         << "\"ChopAndPackLoadBalancer\" or \"TreeLoadBalancer\".");
   }

   return std::shared_ptr<mesh::LoadBalanceStrategy>();
}

std::shared_ptr<mesh::BoxGeneratorStrategy>
createBoxGenerator(
   const std::shared_ptr<tbox::Database>& input_db,
   const std::string& bg_type,
   int ln,
   const tbox::Dimension& dim)
{
   std::shared_ptr<tbox::Database> null_db;

   if (bg_type == "BergerRigoutsos") {

      std::shared_ptr<mesh::BergerRigoutsos>
      berger_rigoutsos(
         new mesh::BergerRigoutsos(
            dim,
            input_db->getDatabaseWithDefault(
               std::string("BergerRigoutsos"), null_db)));
      berger_rigoutsos->useDuplicateMPI(tbox::SAMRAI_MPI::getSAMRAIWorld());
      berger_rigoutsos->setTimerPrefix(
         std::string("mesh::BergerRigoutsos") + tbox::Utilities::intToString(ln));

      return berger_rigoutsos;

   } else if (bg_type == "TileClustering") {

      std::shared_ptr<mesh::TileClustering>
      tile_clustering(
         new mesh::TileClustering(
            dim,
            input_db->getDatabaseWithDefault(
               std::string("TileClustering"), null_db)));
      tile_clustering->setTimerPrefix(
         std::string("mesh::TileClustering") + tbox::Utilities::intToString(ln));

      return tile_clustering;

   } else {
      TBOX_ERROR(
         "Missing or box generator specification in Main database.\n"
         << "Specify load_balancer_type = STRING, where STRING can be\n"
         << "\"BergerRigoutsos\" or \"TileClustering\".");
   }

   return std::shared_ptr<mesh::BoxGeneratorStrategy>();
}

/*
 ****************************************************************************
 * Get the RankTreeStrategy implementation for TreeLoadBalancer
 ****************************************************************************
 */
std::shared_ptr<RankTreeStrategy> getRankTree(
   Database& input_db,
   const std::string& rank_tree_type)
{
   tbox::plog << "Rank tree type is " << rank_tree_type << '\n';

   std::shared_ptr<tbox::RankTreeStrategy> rank_tree;

   if (rank_tree_type == "BalancedDepthFirstTree") {

      BalancedDepthFirstTree * bdfs(new BalancedDepthFirstTree());

      if (input_db.isDatabase("BalancedDepthFirstTree")) {
         std::shared_ptr<tbox::Database> tmp_db = input_db.getDatabase("BalancedDepthFirstTree");
         bool do_left_leaf_switch = tmp_db->getBoolWithDefault("do_left_leaf_switch", true);
         bdfs->setLeftLeafSwitching(do_left_leaf_switch);
      }

      rank_tree.reset(bdfs);

   } else if (rank_tree_type == "CenteredRankTree") {

      CenteredRankTree * crt(new tbox::CenteredRankTree());

      if (input_db.isDatabase("CenteredRankTree")) {
         std::shared_ptr<tbox::Database> tmp_db = input_db.getDatabase("CenteredRankTree");
         bool make_first_rank_the_root = tmp_db->getBoolWithDefault("make_first_rank_the_root",
               true);
         crt->makeFirstRankTheRoot(make_first_rank_the_root);
      }

      rank_tree.reset(crt);

   } else if (rank_tree_type == "BreadthFirstRankTree") {

      BreadthFirstRankTree * dft(new tbox::BreadthFirstRankTree());

      if (input_db.isDatabase("BreadthFirstRankTree")) {
         std::shared_ptr<tbox::Database> tmp_db = input_db.getDatabase("BreadthFirstRankTree");
         const int tree_degree = tmp_db->getIntegerWithDefault("tree_degree", true);
         dft->setTreeDegree(static_cast<unsigned short>(tree_degree));
      }

      rank_tree.reset(dft);

   } else {
      TBOX_ERROR("Unrecognized RankTreeStrategy " << rank_tree_type);
   }

   return rank_tree;
}

/*
 ************************************************************************
 ************************************************************************
 */
void enforceNesting(
   hier::BoxLevel& L1,
   hier::Connector& L0_to_L1,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int coarser_ln)
{
   tbox::pout << "\tEnforcing nesting..." << std::endl;

   const tbox::Dimension dim(hierarchy->getDim());

   const hier::BoxLevel& L0 = L0_to_L1.getBase();

   const size_t cell_count = L1.getGlobalNumberOfCells();

   /*
    * Make L1 nest inside L0 by nesting_width.
    */
   const hier::IntVector nesting_width(dim, hierarchy->getProperNestingBuffer(coarser_ln));
   const hier::IntVector nesting_width_transpose = hier::Connector::convertHeadWidthToBase(
         L0.getRefinementRatio(),
         L1.getRefinementRatio(),
         nesting_width);
   std::shared_ptr<hier::BoxLevel> L1nested;
   std::shared_ptr<hier::MappingConnector> L1_to_L1nested;
   hier::BoxLevelConnectorUtils blcu;
   blcu.computeInternalParts(L1nested,
      L1_to_L1nested,
      L1.findConnectorWithTranspose(L0,
         nesting_width,
         nesting_width_transpose,
         hier::CONNECTOR_CREATE,
         true),
      -nesting_width,
      hierarchy->getGridGeometry()->getDomainSearchTree());
   hier::MappingConnectorAlgorithm mca;
   mca.modify(L0_to_L1,
      *L1_to_L1nested,
      &L1,
      L1nested.get());

   /*
    * Remove overflow nesting.
    */
   blcu.computeInternalParts(L1nested,
      L1_to_L1nested,
      L0_to_L1.getTranspose(),
      hier::IntVector::getZero(dim),
      hierarchy->getGridGeometry()->getDomainSearchTree());
   mca.modify(L0_to_L1,
      *L1_to_L1nested,
      &L1,
      L1nested.get());

   if (cell_count != L1.getGlobalNumberOfCells()) {
      tbox::plog << "\t\tWarning: enforceNesting changed number of cells from " << cell_count
                 << " to " << L1.getGlobalNumberOfCells() << '\n';
   } else {
      tbox::plog << "\t\tenforceNesting left number of cells at " << cell_count << '\n';
   }
}
