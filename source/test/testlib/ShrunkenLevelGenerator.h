/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   ShrunkenLevelGenerator class declaration
 *
 ************************************************************************/
#ifndef included_ShrunkenLevelGenerator
#define included_ShrunkenLevelGenerator

#include "MeshGenerationStrategy.h"

#include <string>
#include <memory>

/*
 * SAMRAI classes
 */
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Timer.h"

#include "DerivedVisOwnerData.h"


using namespace SAMRAI;

/*!
 * @brief Class to tag a full level, shrunken by a given IntVector amount.
 *
 * Inputs:
 *
 * shrink_distance_0, shrink_distance_1, ...:
 * shrink_distance[ln] is the shink distance when tagging ON
 * level ln by shrinking the boundaries of level ln.
 */
class ShrunkenLevelGenerator:
   public MeshGenerationStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   ShrunkenLevelGenerator(
      /*! Ojbect name */
      const std::string& object_name,
      const tbox::Dimension& dim,
      /*! Input database */
      const std::shared_ptr<tbox::Database>& database = std::shared_ptr<tbox::Database>());

   ~ShrunkenLevelGenerator();

   /*!
    * @brief Set tags on the tag level.
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

   virtual void
   resetHierarchyConfiguration(
      /*! New hierarchy */
      const std::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
      /*! Coarsest level */ const int coarsest_level,
      /*! Finest level */ const int finest_level);

   //@}

   /*!
    * @brief Compute shell-dependent data for a patch.
    */
   void
   computePatchData(
      const hier::Patch& patch,
      pdat::CellData<double>* uval_data,
      pdat::CellData<int>* tag_data,
      const hier::Box& fill_box) const {
      NULL_USE(patch);
      NULL_USE(uval_data);
      NULL_USE(tag_data);
      NULL_USE(fill_box);
      TBOX_ERROR("Shrunken Level generator doesn't yet support computePatchData.");
   }

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

private:
   /*!
    * @brief Set tags by shrinking the level at its coarse-fine
    * boundary.
    */
   void
   setTagsByShrinkingLevel(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      int tag_ln,
      int tag_data_id,
      const hier::IntVector& shrink_cells,
      const double* shrink_distance);

   std::string d_name;

   const tbox::Dimension d_dim;

   /*!
    * @brief PatchHierarchy for use in implementations of some
    * abstract interfaces that do not specify a hierarch.
    */
   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   /*!
    * @brief Whether to scale up domain by increasing resolution ('r')
    * or by tiling ('t').
    */
   char d_domain_scale_method;

   /*!
    * @brief Shrink distances for generating tags.
    */
   std::vector<std::vector<double> > d_shrink_distance;

   DerivedVisOwnerData d_vis_owner_data;

};

#endif  // included_ShrunkenLevelGenerator
