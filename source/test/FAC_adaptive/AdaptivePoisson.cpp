/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   AdaptivePoisson class implementation
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/pdat/ArrayDataAccess.h"
#include "patchFcns.h"
#include "AdaptivePoisson.h"
#include "SAMRAI/solv/CellPoissonFACOps.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/InputDatabase.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/geom/CartesianCellDoubleLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianSideDoubleWeightedAverage.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

#include <sstream>
#include <iomanip>
#include <cstring>
#include <stdlib.h>
#include <memory>

#include <cmath>

using namespace SAMRAI;

AdaptivePoisson::AdaptivePoisson(
   const std::string& object_name,
   const tbox::Dimension& dim,
   std::shared_ptr<solv::CellPoissonFACOps>& fac_ops,
   std::shared_ptr<solv::FACPreconditioner>& fac_precond,
   tbox::Database& database,
   /*! Log output stream */ std::ostream* log_stream):
   d_name(object_name),
   d_dim(dim),
   d_fac_ops(fac_ops),
   d_fac_preconditioner(fac_precond),
   d_allocator(tbox::AllocatorDatabase::getDatabase()->getDefaultAllocator()),
   d_context_persistent(new hier::VariableContext("PERSISTENT")),
   d_context_scratch(new hier::VariableContext("SCRATCH")),
   d_diffcoef(new pdat::SideVariable<double>(d_dim, "solution:diffcoef",
                                             hier::IntVector::getOne(d_dim),
                                             d_allocator)),
   d_flux(new pdat::SideVariable<double>(d_dim, "flux",
                                         hier::IntVector::getOne(d_dim),
                                         d_allocator)),
   d_scalar(new pdat::CellVariable<double>(d_dim, "solution:scalar",
                                           d_allocator)),
   d_constant_source(new pdat::CellVariable<double>(d_dim, "poisson source",
                                                    d_allocator)),
   d_ccoef(new pdat::CellVariable<double>(d_dim, "linear source coefficient",
                                          d_allocator)),
   d_rhs(new pdat::CellVariable<double>(d_dim, "linear system rhs",
                                        d_allocator)),
   d_exact(new pdat::CellVariable<double>(d_dim, "solution:exact",
                                          d_allocator)),
   d_resid(new pdat::CellVariable<double>(d_dim, object_name + "residual",
                                          d_allocator)),
   d_weight(new pdat::CellVariable<double>(d_dim, "vector weight",
                                           d_allocator)),
   d_lstream(log_stream),
   d_problem_name("sine"),
   d_sps(object_name + "Poisson solver specifications"),
   d_sine_solution(dim),
   d_gaussian_solution(dim),
   d_multigaussian_solution(dim),
   d_polynomial_solution(dim),
   d_gaussian_diffcoef_solution(dim),
   d_robin_refine_patch(d_dim, object_name + "Refine patch implementation"),
   d_physical_bc_coef(0),
   d_adaption_threshold(0.5),
   d_finest_dbg_plot_ln(database.getIntegerWithDefault("finest_dbg_plot_ln", 99))
{

   /*
    * Register variables with hier::VariableDatabase
    * and get the descriptor indices for those variables.
    * It is not necessary to save the indices obtained from
    * the registration, because they can always be retrieved
    * from the mapVariableAndContextToIndex, but we do it
    * because we refer to the indices often.
    */
   {
      hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

      /*
       * Persistent data.
       */
      d_diffcoef_persistent =
         variable_db->registerVariableAndContext(
            d_diffcoef,
            d_context_persistent,
            hier::IntVector(d_dim, 0));
      d_flux_persistent =
         variable_db->registerVariableAndContext(
            d_flux,
            d_context_persistent,
            hier::IntVector(d_dim, 0));
      d_scalar_persistent =
         variable_db->registerVariableAndContext(
            d_scalar,
            d_context_persistent,
            hier::IntVector(d_dim, 1)   /* ghost cell width is 1 for stencil widths */
            );
      d_constant_source_persistent =
         variable_db->registerVariableAndContext(
            d_constant_source,
            d_context_persistent,
            hier::IntVector(d_dim, 0)
            );
      d_ccoef_persistent =
         variable_db->registerVariableAndContext(
            d_ccoef,
            d_context_persistent,
            hier::IntVector(d_dim, 0)
            );
      d_exact_persistent =
         variable_db->registerVariableAndContext(
            d_exact,
            d_context_persistent,
            hier::IntVector(d_dim, 0)
            );
      d_weight_persistent =
         variable_db->registerVariableAndContext(
            d_weight,
            d_context_persistent,
            hier::IntVector(d_dim, 0)
            );
      /*
       * Scratch data.
       */
      d_rhs_scratch =
         variable_db->registerVariableAndContext(
            d_rhs,
            d_context_scratch,
            hier::IntVector(d_dim, 0)   /* ghost cell width is 0 */
            );
      d_resid_scratch =
         variable_db->registerVariableAndContext(
            d_resid,
            d_context_scratch,
            hier::IntVector(d_dim, 0) /* ghost cell width is 0 */
            );
   }

   /*
    * Experiment with algorithm choices in solv::FACPreconditioner.
    */
   std::string fac_algo = database.getStringWithDefault("fac_algo", "default");
   d_fac_preconditioner->setAlgorithmChoice(fac_algo);

   d_adaption_threshold =
      database.getDoubleWithDefault("adaption_threshold",
         d_adaption_threshold);

   /*
    * Read in the possible solution-specific objects.
    */
   {
      if (database.isDatabase("sine_solution")) {
         std::shared_ptr<tbox::Database> db(
            database.getDatabase("sine_solution"));
         d_sine_solution.setFromDatabase(*db);
      }
      if (database.isDatabase("gaussian_solution")) {
         std::shared_ptr<tbox::Database> db(
            database.getDatabase("gaussian_solution"));
         d_gaussian_solution.setFromDatabase(*db);
      }
      if (database.isDatabase("multigaussian_solution")) {
         std::shared_ptr<tbox::Database> db(
            database.getDatabase("multigaussian_solution"));
         d_multigaussian_solution.setFromDatabase(*db);
      }
      if (database.isDatabase("polynomial_solution")) {
         std::shared_ptr<tbox::Database> db(
            database.getDatabase("polynomial_solution"));
         d_polynomial_solution.setFromDatabase(*db);
      }
      if (database.isDatabase("gaussian_diffcoef_solution")) {
         std::shared_ptr<tbox::Database> db(
            database.getDatabase("gaussian_diffcoef_solution"));
         d_gaussian_diffcoef_solution.setFromDatabase(*db);
      }
   }

   /*
    * We are set up with a choice of problems to solve.
    * The problem name identifies the specific one.
    */

   d_problem_name =
      database.getStringWithDefault("problem_name", d_problem_name);
   if (d_problem_name != "sine"
       && d_problem_name != "gauss"
       && d_problem_name != "multigauss"
       && d_problem_name != "poly"
       && d_problem_name != "gauss-coef"
       ) {
      TBOX_ERROR("Unrecognized problem name " << d_problem_name << "\n");
   }
   if (d_problem_name == "sine") {
      d_sine_solution.setPoissonSpecifications(d_sps,
         d_ccoef_persistent,
         d_diffcoef_persistent);
      d_physical_bc_coef = &d_sine_solution;
   } else if (d_problem_name == "gauss") {
      d_gaussian_solution.setPoissonSpecifications(d_sps,
         d_ccoef_persistent,
         d_diffcoef_persistent);
      d_physical_bc_coef = &d_gaussian_solution;
   } else if (d_problem_name == "multigauss") {
      d_multigaussian_solution.setPoissonSpecifications(d_sps,
         d_ccoef_persistent,
         d_diffcoef_persistent);
      d_physical_bc_coef = &d_multigaussian_solution;
   } else if (d_problem_name == "poly") {
      d_polynomial_solution.setPoissonSpecifications(d_sps,
         d_ccoef_persistent,
         d_diffcoef_persistent);
      d_physical_bc_coef = &d_polynomial_solution;
   } else if (d_problem_name == "gauss-coef") {
      d_gaussian_diffcoef_solution.setPoissonSpecifications(d_sps,
         d_ccoef_persistent,
         d_diffcoef_persistent);
      d_physical_bc_coef = &d_gaussian_diffcoef_solution;
   } else {
      TBOX_ERROR("Unidentified problem name");
   }
   /*
    * Tell ScalarPoissonOperator where to find some of the data
    * we are providing it.
    */
   d_fac_ops->setPoissonSpecifications(d_sps);
   d_fac_ops->setFluxId(-1);

   d_fac_ops->setPhysicalBcCoefObject(d_physical_bc_coef);

   tbox::plog << "Gaussian solution parameters:\n"
              << d_gaussian_solution << "\n\n" << std::endl;
#if 0
   tbox::plog << "Sine solution parameters:\n"
              << d_sine_solution << "\n\n" << std::endl;
   tbox::plog << "Polynomial solution parameters:\n"
              << d_polynomial_solution << "\n\n" << std::endl;
   tbox::plog << "Gaussian diffcoef solution parameters:\n"
              << d_gaussian_diffcoef_solution << "\n\n" << std::endl;
#endif
   tbox::plog << "Problem name is: " << d_problem_name << "\n\n" << std::endl;

   d_fac_ops->setPreconditioner(d_fac_preconditioner.get());
}

void AdaptivePoisson::initializeLevelData(
   /*! Hierarchy to initialize */
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   /*! Level to initialize */
   const int ln,
   const double init_data_time,
   const bool can_be_refined,
   /*! Whether level is being introduced for the first time */
   const bool initial_time,
   /*! Level to copy data from */
   const std::shared_ptr<hier::PatchLevel>& old_level,
   const bool allocate_data)
{
   NULL_USE(init_data_time);
   NULL_USE(can_be_refined);
   NULL_USE(initial_time);

   std::shared_ptr<hier::PatchHierarchy> patch_hierarchy(hierarchy);

   /*
    * Reference the level object with the given index from the hierarchy.
    */
   std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

   /*
    * If instructed, allocate all patch data on the level.
    * Allocate only persistent data.  Scratch data will
    * generally be allocated and deallocated as needed.
    */
   if (allocate_data) {
      if (d_sps.dIsVariable())
         level->allocatePatchData(d_diffcoef_persistent);
      level->allocatePatchData(d_flux_persistent);
      level->allocatePatchData(d_scalar_persistent);
      level->allocatePatchData(d_constant_source_persistent);
      if (d_sps.cIsVariable())
         level->allocatePatchData(d_ccoef_persistent);
      level->allocatePatchData(d_exact_persistent);
      level->allocatePatchData(d_weight_persistent);
   }

   /*
    * Initialize data in all patches in the level.
    */
   for (hier::PatchLevel::iterator pi(level->begin());
        pi != level->end(); ++pi) {

      hier::Patch& patch = **pi;
      hier::Box pbox = patch.getBox();

      std::shared_ptr<pdat::SideData<double> > diffcoef_data(
         SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
            patch.getPatchData(d_diffcoef_persistent)));
      std::shared_ptr<pdat::CellData<double> > exact_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_exact_persistent)));
      std::shared_ptr<pdat::CellData<double> > source_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_constant_source_persistent)));
      TBOX_ASSERT(exact_data);
      TBOX_ASSERT(source_data);

      /* Set source function and exact solution. */
      if (d_problem_name == "sine") {
         d_sine_solution.setGridData(patch,
            *exact_data,
            *source_data);
      } else if (d_problem_name == "gauss") {
         d_gaussian_solution.setGridData(patch,
            *exact_data,
            *source_data);
      } else if (d_problem_name == "multigauss") {
         d_multigaussian_solution.setGridData(patch,
            *exact_data,
            *source_data);
      } else if (d_problem_name == "poly") {
         d_polynomial_solution.setGridData(patch,
            *exact_data,
            *source_data);
      } else if (d_problem_name == "gauss-coef") {
         TBOX_ASSERT(diffcoef_data);
         d_gaussian_diffcoef_solution.setGridData(patch,
            *diffcoef_data,
            *exact_data,
            *source_data);
      } else {
         TBOX_ERROR("Unidentified problem name");
      }

   }

   /*
    * Refine solution data from coarser level and, if provided, old level.
    */
   {
      xfer::RefineAlgorithm refiner;
      std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
            patch_hierarchy->getGridGeometry()));
      TBOX_ASSERT(grid_geometry_);
      geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;
      std::shared_ptr<hier::RefineOperator> accurate_refine_op =
         grid_geometry.
         lookupRefineOperator(d_scalar, "CONSERVATIVE_LINEAR_REFINE");
      TBOX_ASSERT(accurate_refine_op);
      refiner.registerRefine(d_scalar_persistent,
         d_scalar_persistent,
         d_scalar_persistent,
         accurate_refine_op);
      std::shared_ptr<xfer::RefineSchedule> refine_schedule;
      if (ln > 0) {
         /*
          * Include coarser levels in setting data
          */
         refine_schedule =
            refiner.createSchedule(level,
               old_level,
               ln - 1,
               hierarchy,
               &d_robin_refine_patch);
      } else {
         /*
          * There is no coarser level, and source data comes only
          * from old_level, if any.
          */
         if (old_level) {
            refine_schedule =
               refiner.createSchedule(level,
                  old_level,
                  &d_robin_refine_patch);
         }
      }
      if (refine_schedule) {
         d_robin_refine_patch.setCoefImplementation(d_physical_bc_coef);
         d_robin_refine_patch.setTargetDataId(d_scalar_persistent);
         d_robin_refine_patch.setHomogeneousBc(false);
         refine_schedule->fillData(0.0);
         // It is null if this is the bottom level.
      } else {
         math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy, ln, ln);
         hcellmath.setToScalar(d_scalar_persistent, 0.0, false);
      }
      if (0) {
         // begin debug code
         math::HierarchyCellDataOpsReal<double> hcellmath(hierarchy);
         hcellmath.printData(d_scalar_persistent, tbox::pout, false);
         // end debug code
      }
   }

   /* Set vector weight. */
   d_fac_ops->computeVectorWeights(hierarchy, d_weight_persistent);
}

void AdaptivePoisson::resetHierarchyConfiguration(
   /*! New hierarchy */ const std::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
   /*! Coarsest level */ int coarsest_level,
   /*! Finest level */ int finest_level)
{
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);

   d_hierarchy = new_hierarchy;
   /*
    * Recompute or reset internal data tied to the hierarchy,
    * if any.  None at this time.
    */
   /*
    * Log the new hierarchy.
    */
   if (d_lstream) {
      *d_lstream
      << "AdaptivePoisson::resetHierarchyConfiguration\n";
      d_hierarchy->recursivePrint(*d_lstream, "    ", 2);
   }
}

void AdaptivePoisson::applyGradientDetector(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy_,
   const int ln,
   const double error_data_time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation)
{
   NULL_USE(uses_richardson_extrapolation);
   NULL_USE(error_data_time);
   NULL_USE(initial_time);

   if (d_lstream) {
      *d_lstream
      << "AdaptivePoisson(" << d_name << ")::applyGradientDetector"
      << std::endl;
   }
   hier::PatchHierarchy& hierarchy = *hierarchy_;
   hier::PatchLevel& level =
      (hier::PatchLevel &) * hierarchy.getPatchLevel(ln);
   size_t ntag = 0, ntotal = 0;
   double maxestimate = 0;
   for (hier::PatchLevel::iterator pi(level.begin());
        pi != level.end(); ++pi) {
      hier::Patch& patch = **pi;
      std::shared_ptr<hier::PatchData> tag_data(
         patch.getPatchData(tag_index));
      ntotal += patch.getBox().numberCells().getProduct();
      if (!tag_data) {
         TBOX_ERROR(
            "Data index " << tag_index << " does not exist for patch.\n");
      }
      std::shared_ptr<pdat::CellData<int> > tag_cell_data_(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(tag_data));
      TBOX_ASSERT(tag_cell_data_);
      std::shared_ptr<hier::PatchData> soln_data(
         patch.getPatchData(d_scalar_persistent));
      if (!soln_data) {
         TBOX_ERROR("Data index " << d_scalar_persistent
                                  << " does not exist for patch.\n");
      }
      std::shared_ptr<pdat::CellData<double> > soln_cell_data_(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(soln_data));
      TBOX_ASSERT(soln_cell_data_);
      pdat::CellData<double>& soln_cell_data = *soln_cell_data_;
      pdat::CellData<int>& tag_cell_data = *tag_cell_data_;
      pdat::CellData<double> estimate_data(
         patch.getBox(),
         1,
         hier::IntVector(d_dim, 0),
         tbox::AllocatorDatabase::getDatabase()->getTagAllocator());
      computeAdaptionEstimate(estimate_data,
         soln_cell_data);
      tag_cell_data.fill(0);
      hier::Box::iterator iend(patch.getBox().end());
      for (hier::Box::iterator i(patch.getBox().begin()); i != iend; ++i) {
         const pdat::CellIndex cell_index(*i);
         if (maxestimate < estimate_data(cell_index)) maxestimate =
               estimate_data(cell_index);
         if (estimate_data(cell_index) > d_adaption_threshold) {
            tag_cell_data(cell_index) = 1;
            ++ntag;
         }
      }
   }
   tbox::plog << "Adaption threshold is " << d_adaption_threshold << "\n";
   tbox::plog << "Number of cells tagged on level " << ln << " is "
              << ntag << "/" << ntotal << "\n";
   tbox::plog << "Max estimate is " << maxestimate << "\n";
}

#ifdef HAVE_HDF5
int AdaptivePoisson::registerVariablesWithPlotter(
   appu::VisItDataWriter& visit_writer) {

   visit_writer.registerPlotQuantity("Computed solution",
      "SCALAR",
      d_scalar_persistent);
   visit_writer.registerPlotQuantity("Exact solution",
      "SCALAR",
      d_exact_persistent);
   visit_writer.registerPlotQuantity("Poisson source",
      "SCALAR",
      d_constant_source_persistent);
   visit_writer.registerDerivedPlotQuantity("Gradient Function",
      "SCALAR",
      this,
      1.0,
      "CELL");
   visit_writer.registerDerivedPlotQuantity("Patch level number",
      "SCALAR",
      this,
      1.0,
      "CELL");

   std::vector<std::string> expression_keys(1);
   std::vector<std::string> expressions(1);
   std::vector<std::string> expression_types(1);

   {
      expression_keys[0] = "Error";
      expression_types[0] = "scalar";
      expressions[0] = "<Computed solution> - <Exact solution>";
   }

   visit_writer.registerVisItExpressions(expression_keys,
      expressions,
      expression_types);

   return 0;
}
#endif

bool AdaptivePoisson::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_id,
   double simulation_time) const
{
   NULL_USE(depth_id);
   NULL_USE(simulation_time);

   // begin debug code
   // math::HierarchyCellDataOpsReal<double> hcellmath(d_hierarchy);
   // hcellmath.printData( d_exact_persistent, pout, false );
   // end debug code

   if (variable_name == "Gradient Function") {
      std::shared_ptr<pdat::CellData<double> > soln_cell_data_(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_scalar_persistent)));
      TBOX_ASSERT(soln_cell_data_);
      const pdat::CellData<double>& soln_cell_data = *soln_cell_data_;
      pdat::CellData<double> estimate_data(region,
                                           1,
                                           hier::IntVector(d_dim, 0));
      computeAdaptionEstimate(estimate_data,
         soln_cell_data);
      // tbox::plog << "estimate data: " << patch.getBox().size() << "\n";
      // estimate_data.print(region,0,tbox::plog);
      memcpy(buffer, estimate_data.getPointer(), sizeof(double) * region.size());
   } else if (variable_name == "Patch level number") {
      double pln = patch.getPatchLevelNumber();
      for (size_t i = 0; i < region.size(); ++i) buffer[i] = pln;
   } else {
      // Did not register this name.
      TBOX_ERROR(
         "Unregistered variable name '" << variable_name << "' in\n"
                                        << "AdaptivePoisson::packDerivedPatchDataIntoDoubleBuffer");
   }

   // Return TRUE if this patch has derived data on it.
   // FALSE otherwise.
   return true;
}

void AdaptivePoisson::computeAdaptionEstimate(
   pdat::CellData<double>& estimate_data,
   const pdat::CellData<double>& soln_cell_data) const
{
   const int* lower = &estimate_data.getBox().lower()[0];
   const int* upper = &estimate_data.getBox().upper()[0];
   if (d_dim == tbox::Dimension(2)) {
      MDA_AccessConst<double, 2, MDA_OrderColMajor<2> > co =
         pdat::ArrayDataAccess::access<2, double>(soln_cell_data.getArrayData());
      MDA_Access<double, 2, MDA_OrderColMajor<2> > es =
         pdat::ArrayDataAccess::access<2, double>(estimate_data.getArrayData());
      int i, j;
      double estimate, est0, est1, est2, est3, est4, est5;
      for (j = lower[1]; j <= upper[1]; ++j) {
         for (i = lower[0]; i <= upper[0]; ++i) {
            est0 =
               tbox::MathUtilities<double>::Abs(co(i + 1, j) + co(i - 1,
                     j) - 2 * co(i, j));
            est1 =
               tbox::MathUtilities<double>::Abs(co(i, j + 1) + co(i,
                     j - 1) - 2 * co(i, j));
            est2 = 0.5
               * tbox::MathUtilities<double>::Abs(co(i + 1, j
                     + 1) + co(i - 1, j - 1) - 2 * co(i, j));
            est3 = 0.5
               * tbox::MathUtilities<double>::Abs(co(i + 1, j
                     - 1) + co(i - 1, j + 1) - 2 * co(i, j));
            est4 = tbox::MathUtilities<double>::Max(est0, est1);
            est5 = tbox::MathUtilities<double>::Max(est2, est3);
            estimate = tbox::MathUtilities<double>::Max(est4, est5);
            es(i, j) = estimate;
         }
      }
   }
   if (d_dim == tbox::Dimension(3)) {
      MDA_AccessConst<double, 3, MDA_OrderColMajor<3> > co =
         pdat::ArrayDataAccess::access<3, double>(soln_cell_data.getArrayData());
      MDA_Access<double, 3, MDA_OrderColMajor<3> > es =
         pdat::ArrayDataAccess::access<3, double>(estimate_data.getArrayData());
      // math::PatchCellDataOpsReal<double> cops;
      // cops.printData( soln_cell_data_, soln_cell_data_->getGhostBox(), tbox::plog );
      int i, j, k;
      double estimate, est0, est1, est2, est3, est4, est5, est6, est7, est8,
             esta, estb, estc, estd, este, estf, estg;
      for (k = lower[2]; k <= upper[2]; ++k) {
         for (j = lower[1]; j <= upper[1]; ++j) {
            for (i = lower[0]; i <= upper[0]; ++i) {
               est0 =
                  tbox::MathUtilities<double>::Abs(co(i + 1, j, k) + co(i - 1,
                        j,
                        k) - 2 * co(i, j, k));
               est1 =
                  tbox::MathUtilities<double>::Abs(co(i, j + 1, k) + co(i,
                        j - 1,
                        k) - 2 * co(i, j, k));
               est2 =
                  tbox::MathUtilities<double>::Abs(co(i, j, k + 1) + co(i,
                        j,
                        k - 1) - 2 * co(i, j, k));
               est3 = 0.5 * tbox::MathUtilities<double>::Abs(co(i,
                        j + 1,
                        k + 1) + co(i, j - 1, k - 1) - 2 * co(i, j, k));
               est4 = 0.5 * tbox::MathUtilities<double>::Abs(co(i,
                        j + 1,
                        k - 1) + co(i, j - 1, k + 1) - 2 * co(i, j, k));
               est5 = 0.5 * tbox::MathUtilities<double>::Abs(co(i + 1,
                        j,
                        k + 1) + co(i - 1, j, k - 1) - 2 * co(i, j, k));
               est6 = 0.5 * tbox::MathUtilities<double>::Abs(co(i + 1,
                        j,
                        k - 1) + co(i - 1, j, k + 1) - 2 * co(i, j, k));
               est7 = 0.5 * tbox::MathUtilities<double>::Abs(co(i + 1,
                        j + 1,
                        k) + co(i - 1, j - 1, k) - 2 * co(i, j, k));
               est8 = 0.5 * tbox::MathUtilities<double>::Abs(co(i + 1,
                        j - 1,
                        k) + co(i - 1, j + 1, k) - 2 * co(i, j, k));
               esta = tbox::MathUtilities<double>::Max(est0, est1);
               estb = tbox::MathUtilities<double>::Max(est2, est3);
               estc = tbox::MathUtilities<double>::Max(est4, est5);
               estd = tbox::MathUtilities<double>::Max(est6, est7);
               este = tbox::MathUtilities<double>::Max(esta, estb);
               estf = tbox::MathUtilities<double>::Max(estc, estd);
               estg = tbox::MathUtilities<double>::Max(este, estf);
               estimate = tbox::MathUtilities<double>::Max(estg, est8);
               es(i, j, k) = estimate;
            }
         }
      }
   }
}

int AdaptivePoisson::computeError(
   const hier::PatchHierarchy& hierarchy,
   double* l2norm,
   double* linorm,
   std::vector<double>& l2norms,
   std::vector<double>& linorms) const
{

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   int ln;

   /*
    * Compute error on all levels, all patches.
    */
   double diff = 0;
   double l2n = 0, wts = 0, lin = 0;
   /*
    * We give wtsum twice the space required so we can combine
    * the l2norms during the sumReduction, saving a little
    * parallel overhead.
    */
   const int nlevels = hierarchy.getNumberOfLevels();
   std::vector<double> wtsums(2 * nlevels);
   for (ln = nlevels - 1; ln >= 0; --ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy.getPatchLevel(ln));

      double& levelwts(wtsums[ln]);
      double& levell2n(wtsums[ln + nlevels]); // l2n and wts combined in 1 array.
      double& levellin(linorms[ln]);

      levell2n = levellin = levelwts = 0.0;

      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;

         /*
          * Get the patch data.
          */
         std::shared_ptr<pdat::CellData<double> > current_solution(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_scalar_persistent)));
         std::shared_ptr<pdat::CellData<double> > exact_solution(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_exact_persistent)));
         std::shared_ptr<pdat::CellData<double> > weight(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_weight_persistent)));
         TBOX_ASSERT(current_solution);
         TBOX_ASSERT(exact_solution);
         TBOX_ASSERT(weight);

         {
            const int* lower = &current_solution->getBox().lower()[0];
            const int* upper = &current_solution->getBox().upper()[0];
            if (d_dim == tbox::Dimension(2)) {
               MDA_AccessConst<double, 2, MDA_OrderColMajor<2> > ex =
                  pdat::ArrayDataAccess::access<2, double>(
                     exact_solution->getArrayData());
               MDA_AccessConst<double, 2, MDA_OrderColMajor<2> > co =
                  pdat::ArrayDataAccess::access<2, double>(
                     current_solution->getArrayData());
               MDA_AccessConst<double, 2, MDA_OrderColMajor<2> > wt =
                  pdat::ArrayDataAccess::access<2, double>(weight->getArrayData());
               for (int j = lower[1]; j <= upper[1]; ++j) {
                  for (int i = lower[0]; i <= upper[0]; ++i) {
                     /*
                      * Disregard zero weights in error computations
                      * because they are on coarse grids covered by finer grids.
                      */
                     if (wt(i, j) > 0) {
                        diff =
                           tbox::MathUtilities<double>::Abs(co(i, j) - ex(i, j));
                        if (levellin < diff) levellin = diff;
                        levell2n += wt(i, j) * diff * diff;
                        levelwts += wt(i, j);
                     }
                  }
               }
            }
            if (d_dim == tbox::Dimension(3)) {
               MDA_AccessConst<double, 3, MDA_OrderColMajor<3> > ex =
                  pdat::ArrayDataAccess::access<3, double>(
                     exact_solution->getArrayData());
               MDA_AccessConst<double, 3, MDA_OrderColMajor<3> > co =
                  pdat::ArrayDataAccess::access<3, double>(
                     current_solution->getArrayData());
               MDA_AccessConst<double, 3, MDA_OrderColMajor<3> > wt =
                  pdat::ArrayDataAccess::access<3, double>(weight->getArrayData());
               for (int k = lower[2]; k <= upper[2]; ++k) {
                  for (int j = lower[1]; j <= upper[1]; ++j) {
                     for (int i = lower[0]; i <= upper[0]; ++i) {
                        /*
                         * Disregard zero weights in error computations
                         * because they are on coarse grids covered by finer grids.
                         */
                        if (wt(i, j, k) > 0) {
                           diff = tbox::MathUtilities<double>::Abs(co(i,
                                    j,
                                    k) - ex(i, j, k));
                           if (levellin < diff) levellin = diff;
                           levell2n += wt(i, j, k) * diff * diff;
                           levelwts += wt(i, j, k);
                        }
                     }
                  }
               }
            }
         }
      } // end patch loop

   } // end level loop

   if (mpi.getSize() > 1) {
      /*
       * Communicate global data if in parallel.
       * We temporarily combine l2norms and wtsum so we can sumReduction
       * in one shot, saving some parallel overhead.
       */
      if (mpi.getSize() > 1) {
         mpi.AllReduce(&wtsums[0], 2 * nlevels, MPI_SUM);
      }
      if (mpi.getSize() > 1) {
         mpi.AllReduce(&linorms[0], nlevels, MPI_SUM);
      }
   }

   for (ln = 0; ln < nlevels; ++ln) {
      /* Copy l2norm accumulated temporarily in wtsums */
      l2norms[ln] = wtsums[ln + nlevels];
      /* Data for whole hierarchy. */
      l2n += l2norms[ln];
      wts += wtsums[ln];
      lin = linorms[ln] > lin ? linorms[ln] : lin;
      /*
       * Data for level ln.
       * If a level is completely covered by a finer level,
       * wtsums[ln] will legitimately be zero, so that protect
       * it from a zero divide.
       */
      l2norms[ln] =
         tbox::MathUtilities<double>::equalEps(wtsums[ln], 0) ? 0 : sqrt(
            l2norms[ln] / wtsums[ln]);

   }

   if (!tbox::MathUtilities<double>::equalEps(wtsums[ln], 0)) {
      *l2norm = sqrt(l2n / wts);
   } else {
      *l2norm = 0.0;
   }

   *linorm = lin;

   return 0;
}

int AdaptivePoisson::solvePoisson(
   std::shared_ptr<hier::PatchHierarchy> hierarchy,
   std::string initial_u)
{

   const int finest_ln = hierarchy->getFinestLevelNumber();
   const int coarsest_ln = 0;

   /*
    * Allocate scratch data for use in the solve.
    */
   for (int ln = coarsest_ln; ln <= finest_ln; ++ln) {
      hierarchy->getPatchLevel(ln)->allocatePatchData(d_rhs_scratch);
   }

   /*
    * Create vectors x and b for solving Ax=b.
    */
   solv::SAMRAIVectorReal<double>
   x("solution", hierarchy, coarsest_ln, finest_ln),
   b("rhs", hierarchy, coarsest_ln, finest_ln);
   x.addComponent(d_scalar, d_scalar_persistent, d_weight_persistent);
   b.addComponent(d_rhs, d_rhs_scratch, d_weight_persistent);

   /*
    * Fill the rhs vector (by filling the d_rhs_scratch data
    * that forms its component).
    * Fill the boundary condition coefficient data.
    */
   for (int ln = coarsest_ln; ln <= finest_ln; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

      for (hier::PatchLevel::iterator pi(level->begin());
           pi != level->end(); ++pi) {
         const std::shared_ptr<hier::Patch>& patch = *pi;

         const hier::Box& box = patch->getBox();

         std::shared_ptr<pdat::CellData<double> > source_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_constant_source_persistent)));
         std::shared_ptr<pdat::CellData<double> > rhs_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_rhs_scratch)));
         TBOX_ASSERT(source_data);
         TBOX_ASSERT(rhs_data);
         math::PatchCellDataOpsReal<double> cell_ops;
         cell_ops.scale(rhs_data, 1.0, source_data, box);

      }
   }

   // debug:
   if (0) {
      x.setToScalar(1.0);
      double weightsums =
         d_fac_ops->computeResidualNorm(x, finest_ln, coarsest_ln);
      tbox::pout << "weightsums: " << weightsums << std::endl;
   }

   /*
    * Fill the vector x with the initial guess, if one is given.
    * If not given, we assume the initial guess is in place.
    */
   if (!initial_u.empty()) {
      if (initial_u == "random") {
         x.setRandomValues(1.0, 0.0);
      } else {
         x.setToScalar(atof(initial_u.c_str()));
      }
   }

   /*
    * Create the viz data writer for use in debugging.
    */
#ifdef HAVE_HDF5
   d_visit_writer.reset(new appu::VisItDataWriter(d_dim,
         "Internal VisIt Writer",
         "ap-debug.visit"));
   registerVariablesWithPlotter(*d_visit_writer);
#endif

   /*
    * Set up FAC preconditioner object.
    */
   if (d_lstream) {
      d_fac_preconditioner->printClassData(*d_lstream);
   }
   d_fac_preconditioner->initializeSolverState(x, b);

   /*
    * Solve the system.
    */
   d_fac_preconditioner->solveSystem(x, b);
   if (d_lstream) *d_lstream
      << "FAC solve completed with\n"
      << std::setw(30) << "number of iterations: "
      << d_fac_preconditioner->getNumberOfIterations() << "\n"
      << std::setw(30) << "residual norm: "
      << d_fac_preconditioner->getResidualNorm() << "\n"
      ;
   d_fac_preconditioner->deallocateSolverState();

   /*
    * Get data on the solve.
    */
   double avg_convergence_factor, final_convergence_factor;
   d_fac_preconditioner->getConvergenceFactors(avg_convergence_factor,
      final_convergence_factor);
   if (d_lstream) *d_lstream
      << "Final result: \n"
      << std::setw(30) << "average convergence factor: "
      << avg_convergence_factor << "\n"
      << std::setw(30) << "final convergence factor: "
      << final_convergence_factor << "\n"
      ;

   /*
    * Fill in boundary ghosts here to get the correct ghost cells
    * values used to compute the gradient estimator when plotting.
    * We are not sure what state ghost cell values are in after
    * the solver finishes.
    */
   {
      xfer::RefineAlgorithm refiner;
      std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
            hierarchy->getGridGeometry()));
      TBOX_ASSERT(grid_geometry_);
      geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;
      std::shared_ptr<hier::RefineOperator> accurate_refine_op(
         grid_geometry.lookupRefineOperator(d_scalar, "LINEAR_REFINE"));
      TBOX_ASSERT(accurate_refine_op);
      refiner.registerRefine(d_scalar_persistent,
         d_scalar_persistent,
         d_scalar_persistent,
         accurate_refine_op);
      d_robin_refine_patch.setTargetDataId(d_scalar_persistent);
      d_robin_refine_patch.setHomogeneousBc(false);
      d_robin_refine_patch.setCoefImplementation(d_physical_bc_coef);
      std::shared_ptr<xfer::RefineSchedule> refine_schedule;
      for (int ln = coarsest_ln; ln <= finest_ln; ++ln) {
         std::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));
         if (ln > 0) {
            /* Include coarser levels in setting data */
            refine_schedule =
               refiner.createSchedule(level,
                  ln - 1,
                  hierarchy,
                  &d_robin_refine_patch);
         } else {
            /* Exclude coarser levels in setting data */
            refine_schedule =
               refiner.createSchedule(level,
                  &d_robin_refine_patch);
         }
         refine_schedule->fillData(0.0);
      }
   }
   // x.print(plog,false);

   /*
    * Deallocate scratch data.
    */
   for (int ln = coarsest_ln; ln <= finest_ln; ++ln) {
      hierarchy->getPatchLevel(ln)->deallocatePatchData(d_rhs_scratch);
   }

#ifdef HAVE_HDF5
   /*
    * Destroy the viz data writer used for debugging.
    */
   d_visit_writer.reset();
#endif

   return 0;
}
