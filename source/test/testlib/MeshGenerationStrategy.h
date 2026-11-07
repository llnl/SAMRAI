/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Strategy class for MeshGeneration performance tests.
 *
 ************************************************************************/
#ifndef included_MeshGenerationStrategy
#define included_MeshGenerationStrategy

#include <string>
#include <memory>

/*
 * SAMRAI classes
 */
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Timer.h"


using namespace SAMRAI;

/*!
 * @brief Interface definition for MeshGeneration performance test.
 *
 * This is a combination of mesh::StandardTagAndInitStrategy and
 * appu::VisDerivedDataStrategy, giving no-op implementations to
 * unneeded methods.
 */
class MeshGenerationStrategy:
   public mesh::StandardTagAndInitStrategy,
   public appu::VisDerivedDataStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   MeshGenerationStrategy() {
   }

   virtual ~MeshGenerationStrategy() {
   }

   /*!
    * @brief Set tag on the tag level.
    *
    * @param [o] exact_tagging Set to true if the implementation wants
    * the clustering to match the tags exactly.  Exact clustering
    * match means using clustering efficiency of 1.0.
    */
   virtual void
   setTags(
      bool& exact_tagging,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      int tag_ln,
      int tag_data_id) = 0;

   //@{ @name SAMRAI::mesh::StandardTagAndInitStrategy virtuals

public:
   /*!
    * @brief Set the domain, possibly scaling up the specifications.
    *
    * Take the domain_boxes, xlo and xhi to be the size for the value
    * of autoscale_base_nprocs.  Scale the problem from there to the
    * number of process specified by the argument mpi.
    *
    * The object should also adjust whatever internal data needed to
    * reflect the scaled-up domain.
    *
    * @param [i/o] domain Domain description to be scaled up (or
    * overridden).
    *
    * @param [i/o] xlo Domain lower physical coordinate to be scaled
    * up (or overridden).
    *
    * @param [i/o] xhi Domain upper physical coordinate to be scaled
    * up (or overridden).
    *
    * @param [i] autoscale_base_nprocs Scale up the domain based on
    * the current definition being for this many processes.  Scale the
    * domain up to the number of processes in the given SAMRAI_MPI.
    *
    * @param [i] mpi
    */
   virtual void setDomain(
      hier::BoxContainer & domain,
      double xlo[],
      double xhi[],
      int autoscale_base_nprocs,
      const tbox::SAMRAI_MPI & mpi) = 0;

   /*!
    * @brief Allocate and initialize data for a new level
    * in the patch hierarchy.
    *
    * This is where you implement the code for initialize data on the
    * grid.  Nevermind when it is called or where in the program that
    * happens.  All the information you need to initialize the grid
    * are in the arguments.
    *
    * @see SAMRAI::mesh::StandardTagAndInitStrategy::initializeLevelData()
    */
   void
   initializeLevelData(
      /*! Hierarchy to initialize */
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      /*! Level to initialize */
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      /*! Whether level is being introduced for the first time */
      const bool initial_time,
      /*! Level to copy data from */
      const std::shared_ptr<hier::PatchLevel>& old_level =
         std::shared_ptr<hier::PatchLevel>(),
      /*! Whether data on new patch needs to be allocated */
      const bool allocate_data = true)
   {
      NULL_USE(hierarchy);
      NULL_USE(level_number);
      NULL_USE(init_data_time);
      NULL_USE(can_be_refined);
      NULL_USE(initial_time);
      NULL_USE(old_level);
      NULL_USE(allocate_data);
      TBOX_ERROR("Should not be here");
   }

   virtual void
   resetHierarchyConfiguration(
      /*! New hierarchy */
      const std::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
      /*! Coarsest level */ int coarsest_level,
      /*! Finest level */ int finest_level) = 0;

   void
   applyGradientDetector(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double error_data_time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation)
   {
      NULL_USE(hierarchy);
      NULL_USE(level_number);
      NULL_USE(error_data_time);
      NULL_USE(tag_index);
      NULL_USE(initial_time);
      NULL_USE(uses_richardson_extrapolation);
      TBOX_ERROR("Should not be here");
   }

   //@}

   /*!
    * @brief Compute tag and/or scalar solution on a patch.
    */
   virtual void
   computePatchData(
      const hier::Patch& patch,
      pdat::CellData<double>* uval_data,
      pdat::CellData<int>* tag_data,
      const hier::Box& fill_box) const = 0;

#ifdef HAVE_HDF5
   /*!
    * @brief Tell a VisIt plotter which data to write for this class.
    */
   virtual int
   registerVariablesWithPlotter(
      appu::VisItDataWriter& writer) {
      NULL_USE(writer);
      return 0;
   }
#endif

   virtual bool
   packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& variable_name,
      int depth_index,
      double simulation_time) const
   {
      NULL_USE(buffer);
      NULL_USE(patch);
      NULL_USE(region);
      NULL_USE(variable_name);
      NULL_USE(depth_index);
      NULL_USE(simulation_time);
      TBOX_ERROR("Should not be here");
      return false;
   }

private:
};

#endif  // MeshGenerationStrategy
