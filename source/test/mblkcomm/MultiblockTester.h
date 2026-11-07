/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Manager class for patch data communication tests.
 *
 ************************************************************************/

#ifndef included_MultiblockTester
#define included_MultiblockTester

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/geom/GridGeometry.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"

#include <memory>

using namespace SAMRAI;

class PatchMultiblockTestStrategy;

/**
 * Class MultiblockTester serves as a tool to test data communication operations
 * in SAMRAI, such as refining and filling ghost cells.
 *
 * The functions in this class called from main() are:
 * \begin{enumerate}
 *    - [MultiblockTester(...)] constructor which initializes object state and
 *                            creates patch hierarchy and sets initial data.
 *    - [createRefineSchedule(...)] creates communication schedule for
 *                                      refining data to given level.
 *    - [performRefineOperations(...)] refines data to given level.
 * \end{enumerate}
 */

class MultiblockTester:
   public mesh::StandardTagAndInitStrategy,
   public xfer::RefinePatchStrategy,
   public xfer::SingularityPatchStrategy
{
public:
   /**
    * Constructor performs basic setup operations.
    */
   MultiblockTester(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<tbox::Database>& main_input_db,
      std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      PatchMultiblockTestStrategy* strategy,
      const std::string& refine_option = "INTERIOR_FROM_SAME_LEVEL");

   /**
    * Destructor is empty.
    */
   ~MultiblockTester();

   /**
    * Return pointer to patch hierarchy on which communication is tested.
    */
   std::shared_ptr<hier::PatchHierarchy> getPatchHierarchy()
   const
   {
      return d_patch_hierarchy;
   }

   /**
    * Register variable for communication testing.
    *
    * The transfer operator look-up will use the src_variable.
    */
   void
   registerVariable(
      const std::shared_ptr<hier::Variable> src_variable,
      const std::shared_ptr<hier::Variable> dst_variable,
      const hier::IntVector& src_ghosts,
      const hier::IntVector& dst_ghosts,
      const std::shared_ptr<hier::BaseGridGeometry> xfer_geom,
      const std::string& operator_name);

   /**
    * Register variable for communication testing.
    *
    * The transfer operator look-up will use the src_variable.
    */
   void
   registerVariableForReset(
      const std::shared_ptr<hier::Variable> src_variable,
      const std::shared_ptr<hier::Variable> dst_variable,
      const hier::IntVector& src_ghosts,
      const hier::IntVector& dst_ghosts,
      const std::shared_ptr<hier::BaseGridGeometry> xfer_geom,
      const std::string& operator_name);

   /**
    * Create communication schedules for refining data to given level.
    */
   void
   createRefineSchedule(
      const int level_number);
   void
   resetRefineSchedule(
      const int level_number);

   /**
    * Refine data to specified level (or perform interpatch communication
    * on that level).
    */
   void
   performRefineOperations(
      const int level_number);

   /**
    * After communication operations are performed, check results.
    */
   bool
   verifyCommunicationResults() const;

   /**
    * Operations needed by GriddingAlgorithm to construct and
    * initialize levels in patch hierarchy.  These operations are
    * pure virtual in GradientDetectorStrategy.
    */
   void
   initializeLevelData(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      const bool initial_time,
      const std::shared_ptr<hier::PatchLevel>& old_level =
         std::shared_ptr<hier::PatchLevel>(),
      const bool allocate_data = true);

   void
   resetHierarchyConfiguration(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level);

   void
   applyGradientDetector(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation_too);

   /**
    * These routines pass off physicial boundary and pre/postprocess
    * refine operations to patch data test object.  They are
    * pure virtual in RefinePatchStrateg.
    */
   void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double time,
      const hier::IntVector& gcw);

   /*!
    * Set the ghost data at a multiblock singularity.
    */
   void
   fillSingularityBoundaryConditions(
      hier::Patch& patch,
      const hier::PatchLevel& encon_level,
      std::shared_ptr<const hier::Connector> dst_to_encon,
      const hier::Box& fill_box,
      const hier::BoundaryBox& boundary_box,
      const std::shared_ptr<hier::BaseGridGeometry>& grid_geometry);

   hier::IntVector
   getRefineOpStencilWidth(
      const tbox::Dimension& dim) const;

   void
   preprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);

   void
   postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);

   double getLevelDt(
      const std::shared_ptr<hier::PatchLevel>& level,
      const double dt_time,
      const bool initial_time)
   {
      NULL_USE(level);
      NULL_USE(dt_time);
      NULL_USE(initial_time);
      return 0.0;
   }

   /*
    * Construct patch hierarchy and initialize data prior to tests.
    */
   void
   setupHierarchy(
      std::shared_ptr<tbox::Database> main_input_db,
      std::shared_ptr<mesh::StandardTagAndInitialize> cell_tagger);

   /*!
    * @brief Return the dimension of this object.
    */
   const tbox::Dimension& getDim() const
   {
      return d_dim;
   }

private:
   /*
    * Object name for error reporting.
    */
   std::string d_object_name;

   const tbox::Dimension d_dim;

   /*
    * Object supplying operatins for particular patch data test.
    */
   PatchMultiblockTestStrategy* d_data_test_strategy;

   /*
    * String name for refine option; ; i.e., source of interior patch
    * data on refined patches.  Options are "INTERIOR_FROM_SAME_LEVEL"
    * and "INTERIOR_FROM_COARSER_LEVEL".
    */
   std::string d_refine_option;

   /*
    * Patch hierarchy on which tests occur.
    */
   std::shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;

   /*
    * Dummy time stamp for all data operations.
    */
   double d_fake_time;

   /*
    * Dummy cycle for all data operations.
    */
   int d_fake_cycle;

   /*
    * The MultiblockTester uses two variable contexts for each variable.
    * The "source", and "destination" contexts indicate the source
    * and destination patch data for the transfer operation.
    *
    * The "refine_scratch" context is used for managing scratch
    * space during refine operations.
    */
   std::shared_ptr<hier::VariableContext> d_source;
   std::shared_ptr<hier::VariableContext> d_destination;
   std::shared_ptr<hier::VariableContext> d_refine_scratch;

   std::shared_ptr<hier::VariableContext> d_reset_source;
   std::shared_ptr<hier::VariableContext> d_reset_destination;
   std::shared_ptr<hier::VariableContext> d_reset_refine_scratch;

   /*
    * Component selector for allocation/deallocation of variable data.
    */
   hier::ComponentSelector d_patch_data_components;

   xfer::RefineAlgorithm d_reset_refine_algorithm;

   std::shared_ptr<xfer::RefineAlgorithm> d_mblk_refine_alg;

   bool d_is_reset;

   std::vector<std::shared_ptr<xfer::RefineSchedule> > d_refine_schedule;

};

#endif
