/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   SphericalShellGenerator class declaration
 *
 ************************************************************************/
#ifndef included_SphericalShellGenerator
#define included_SphericalShellGenerator

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
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Timer.h"

#include "DerivedVisOwnerData.h"


using namespace SAMRAI;

/*!
 * @brief Class to tag spherical shell patternse in given domain.
 *
 * Inputs:
 *
 * radii:
 * Starting at shell origin, tag cells that intersect regions defined by
 * radii[0]<r<radii[1], radii[2]<r<radii[3], radii[2i]<r<radii[2i+1] and so on.
 *
 * buffer_distance_0, buffer_distance_1, ...:
 * buffer_distance[ln] is the buffer distance when tagging ON
 * level ln.  We tag the shells and buffer the tags by this amount.
 * Missing buffer distances will use the last values given.
 * Default is zero buffering.
 */
class SphericalShellGenerator:
   public MeshGenerationStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   SphericalShellGenerator(
      /*! Ojbect name */
      const std::string& object_name,
      const tbox::Dimension& dim,
      /*! Input database */
      const std::shared_ptr<tbox::Database>& database = std::shared_ptr<tbox::Database>());

   ~SphericalShellGenerator();

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
      const hier::Box& fill_box) const;

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
   void
   computeShellsData(
      pdat::CellData<double>* uval_data,
      pdat::CellData<int>* tag_data,
      const hier::Box& fill_box,
      const std::vector<double>& buffer_distance,
      const double xlo[],
      const double dx[]) const;

   std::string d_name;

   const tbox::Dimension d_dim;

   /*!
    * @brief PatchHierarchy for use in implementations of some
    * abstract interfaces that do not specify a hierarch.
    */
   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   /*!
    * @brief Constant time shift to be added to simulation time.
    */
   double d_time_shift;

   /*!
    * @brief Center of shells at time zero.
    */
   double d_init_center[SAMRAI::MAX_DIM_VAL];

   /*!
    * @brief Shell velocity.
    */
   double d_velocity[SAMRAI::MAX_DIM_VAL];

   /*!
    * @brief Radii of shells.
    */
   std::vector<double> d_radii;

   /*!
    * @brief Buffer distances for generating tags.
    */
   std::vector<std::vector<double> > d_buffer_distance;

   DerivedVisOwnerData d_vis_owner_data;

};

#endif  // included_SphericalShellGenerator
