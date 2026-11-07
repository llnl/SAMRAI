/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   SinusoidalFrontGenerator class declaration
 *
 ************************************************************************/
#ifndef included_SinusoidalFrontGenerator
#define included_SinusoidalFrontGenerator

#include "MeshGenerationStrategy.h"

#include <string>
#include <memory>

/*
 * SAMRAI classes
 */
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Timer.h"

#include "DerivedVisOwnerData.h"


using namespace SAMRAI;

/*!
 * @brief Class to tag a sinusoidal "front" in given domain.
 *
 * Inputs:
 *
 * init_disp: Initial displacement of the front.
 *
 * period: Period of the front.
 *
 * amplitude: Amplitude of the front.
 *
 * buffer_distance_0, buffer_distance_1, ...:
 * buffer_distance[ln] is the buffer distance when tagging ON
 * level ln.  We tag the fronts and buffer the tags by this amount.
 * Missing buffer distances will use the last values given.
 * Default is zero buffering.
 */
class SinusoidalFrontGenerator:
   public MeshGenerationStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   SinusoidalFrontGenerator(
      /*! Ojbect name */
      const std::string& object_name,
      const tbox::Dimension& dim,
      /*! Input database */
      const std::shared_ptr<tbox::Database>& database = std::shared_ptr<tbox::Database>());

   ~SinusoidalFrontGenerator();

   /*!
    * @brief Set tas on the tag level.
    */
   virtual void
   setTags(
      bool& exact_tagging,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      int tag_ln,
      int tag_data_id);

   //@{ @name SAMRAI::mesh::StandardTagAndInitStrategy virtuals

public:
   /*!
    * @brief Set the domain, possibly scaling up the specifications.
    *
    * Take the domain_boxes, xlo and xhi to be the size for the
    * (integer) value of autoscale_base_nprocs.  Scale the problem
    * from there to the number of process running by doubling the
    * size starting with the j direction.
    *
    * The number of processes must be a power of 2 times the value
    * of autoscale_base_nprocs.
    */
   void setDomain(
      hier::BoxContainer & domain,
      double xlo[],
      double xhi[],
      int autoscale_base_nprocs,
      const tbox::SAMRAI_MPI & mpi);

   //@}

   bool
   packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& variable_name,
      int depth_index,
      double simulation_time) const;

public:
#ifdef HAVE_HDF5
   /*!
    * @brief Tell a VisIt plotter which data to write for this class.
    */
   int
   registerVariablesWithPlotter(
      appu::VisItDataWriter& writer);
#endif

   /*!
    * @brief Compute distance and tag data for a patch.
    *
    * This method is not specific to data on the hierarchy,
    * so it is of more general use.  It does not require the
    * hierarchy.
    */
   void
   computeFrontsData(
      pdat::NodeData<double>* dist_data,
      pdat::CellData<double>* uval_data,
      pdat::CellData<int>* tag_data,
      const hier::Box& fill_box,
      const std::vector<double>& buffer_distance,
      const double xlo[],
      const double dx[]) const;

   /*
    * Compute patch data allocated by this class, on a hierarchy.
    */
   void
   computeHierarchyData(
      hier::PatchHierarchy& hierarchy,
      double time) {
      NULL_USE(hierarchy);
      d_time_shift = time;
   }

   /*!
    * @brief Compute front-dependent data for a patch.
    */
   void
   computePatchData(
      const hier::Patch& patch,
      pdat::CellData<double>* uval_data,
      pdat::CellData<int>* tag_data,
      const hier::Box& fill_box) const;

   /*!
    * @brief Deallocate internally managed patch data on level.
    */
   void
   deallocatePatchData(
      hier::PatchLevel& level) {
      NULL_USE(level);
   }

   /*!
    * @brief Deallocate internally managed patch data on hierarchy.
    */
   void
   deallocatePatchData(
      hier::PatchHierarchy& hierarchy) {
      NULL_USE(hierarchy);
   }

   //@{ @name SAMRAI::mesh::StandardTagAndInitStrategy virtuals

   virtual void
   resetHierarchyConfiguration(
      /*! New hierarchy */
      const std::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
      /*! Coarsest level */ const int coarsest_level,
      /*! Finest level */ const int finest_level);

   /*!
    * @brief Allocate and initialize data for a new level
    * in the patch hierarchy.
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
   }

   virtual void
   applyGradientDetector(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double error_data_time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation);

   //@}

   /*!
    * @brief Set the independent time variable in the front equation.
    */
   void
   setTime(double time) {
      d_time_shift = time;
   }

private:
   std::string d_name;

   const tbox::Dimension d_dim;

   /*!
    * @brief PatchHierarchy for use in implementations of some
    * abstract interfaces that do not specify a hierarch.
    */
   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   /*!
    * @brief Period of sinusoid.
    */
   double d_period[SAMRAI::MAX_DIM_VAL];

   /*!
    * @brief Initial displacement.
    */
   double d_init_disp[SAMRAI::MAX_DIM_VAL];

   /*!
    * @brief Front velocity.
    */
   double d_velocity[SAMRAI::MAX_DIM_VAL];

   /*!
    * @brief Constant time shift to be added to simulation time.
    */
   double d_time_shift;

   /*!
    * @brief Amplitude of sinusoid.
    */
   double d_amplitude;

   /*!
    * @brief Buffer distances for generating tags.
    */
   std::vector<std::vector<double> > d_buffer_distance;

   DerivedVisOwnerData d_vis_owner_data;

   std::shared_ptr<tbox::Timer> t_setup;
   std::shared_ptr<tbox::Timer> t_node_pos;
   std::shared_ptr<tbox::Timer> t_distance;
   std::shared_ptr<tbox::Timer> t_uval;
   std::shared_ptr<tbox::Timer> t_tag_cells;

};

#endif  // included_SinusoidalFrontGenerator
