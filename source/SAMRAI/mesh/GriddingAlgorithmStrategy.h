/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   AMR hierarchy generation and regridding routines.
 *
 ************************************************************************/

#ifndef included_mesh_GriddingAlgorithmStrategy
#define included_mesh_GriddingAlgorithmStrategy

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/TagAndInitializeStrategy.h"
#include "SAMRAI/tbox/Serializable.h"

#include <vector>
#include <memory>

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Virtual base class providing interface for gridding
 * algorithms.
 *
 * This class defines operations for building and regridding
 * hierarchies.  The interfaces imply that the implementation of this
 * strategy class works on some hierarchy.
 *
 * @see GriddingAlgorithm
 */

class GriddingAlgorithmStrategy
{
public:
   /*!
    * @brief Constructor
    */
   GriddingAlgorithmStrategy();

   /*!
    * @brief Virtual destructor for GriddingAlgorithmStrategy.
    */
   virtual ~GriddingAlgorithmStrategy();

   /*!
    * @brief Construct the coarsest level in the hierarchy (i.e., level 0).
    *
    * @param level_time Simulation time when level is constructed
    */
   virtual void
   makeCoarsestLevel(
      const double level_time) = 0;

   /*!
    * @brief Attempts to create a new level in the hierarchy finer
    * than the finest level currently residing in the hierarchy.
    *
    * It should select cells for refinement on the finest level and
    * construct a new finest level, if necessary.  If no cells are
    * selected for refinement, no new level should be added to the
    * hierarchy.
    *
    * The boolean argument initial_cycle indicates whether the routine
    * is called at the initial simulation cycle.  If true, this routine
    * is being used to build individual levels during the construction of
    * the AMR hierarchy at the initial simulation time.  If false, the
    * routine is being used to add new levels to the hierarchy at some
    * later point.  In either case, the time value is the current
    * simulation time.  Note that this routine cannot be used to
    * construct the coarsest level in the hierarchy (i.e., level 0).
    * The routine makeCoarsestLevel() above must be used for that
    * purpose.
    *
    * The tag buffer indicates the number of cells by which cells
    * selected for refinement should be buffered before new finer
    * level boxes are constructed.
    *
    * @param tag_buffer Size of buffer around tagged cells that will be
    *                   covered by the fine level.  Must be non-negative.
    * @param initial_cycle Must be true if level_cycle is the initial cycle
    *                      of the simulation, false otherwise
    * @param cycle Current simulation cycle.
    * @param level_time Current simulation time.
    * @param regrid_start_time The simulation time when the regridding
    *                          operation began (this parameter is ignored
    *                          except when using Richardson extrapolation)
    */
   virtual void
   makeFinerLevel(
      const int tag_buffer,
      const bool initial_cycle,
      const int cycle,
      const double level_time,
      const double regrid_start_time = 0.0) = 0;

   /*!
    * @brief Attempt to regrid each level in the PatchHierarchy
    * that is finer than the specified level.
    *
    * The given level number is that of the coarsest level on which
    * cells will be selected for refinement.  In other words, that
    * level is the finest level that will not be subject to a change
    * in its patch configuration during the regridding process.
    * Generally, this routine should be used to alter the pre-existing
    * AMR patch hierarchy based on the need to adapt the computational
    * mesh around some phenomenon of interest.  The routine
    * makeFinerLevel() above should be used to construct the initial
    * hierarchy configuration or to add more than one new level into
    * the hierarchy.  Also, this routine will not reconfigure the
    * patches on level 0 (i.e., the coarsest in the hierarchy).  The
    * routine makeCoarsestLevel() above is provided for that purpose.
    *
    * The boolean level_is_coarsest_to_sync is used for regridding in
    * time-dependent problems.  When true, it indicates that the
    * specified level is the coarsest level to synchronize at the
    * current regrid time before this regridding method is called.
    * This is a pretty idiosyncratic argument but allows some
    * flexibility in the way memory is managed during time-dependent
    * regridding operations.
    *
    * @param level_number Coarsest level on which cells will be tagged for
    *                     refinement
    * @param tag_buffer See tag_buffer in makeFinerLevel.
    *                   There is one value per level in the hierarchy.
    * @param cycle Simulation cycle when regridding occurs
    * @param level_time Simulation time of the level corresponding to level_num
    *                   when regridding occurs
    * @param regrid_start_time The simulation time when the regridding
    *                          operation began on each level (this parameter
    *                          is ignored except when using Richardson
    *                          extrapolation, where it is passed on to the
    *                          estimation methods.)
    * @param level_is_coarsest_to_sync Level is the coarsest to sync
    */
   virtual void
   regridAllFinerLevels(
      const int level_number,
      const std::vector<int>& tag_buffer,
      const int cycle,
      const double level_time,
      const std::vector<double>& regrid_start_time = std::vector<double>(),
      const bool level_is_coarsest_to_sync = true) = 0;

   /*!
    * @brief Return pointer to level gridding strategy data member.
    */
   virtual
   std::shared_ptr<TagAndInitializeStrategy>
   getTagAndInitializeStrategy() const = 0;

   /*!
    * @brief Return pointer to PatchHierarchy data member.
    */
   virtual
   std::shared_ptr<hier::PatchHierarchy>
   getPatchHierarchy() const = 0;

};

}
}
#endif
