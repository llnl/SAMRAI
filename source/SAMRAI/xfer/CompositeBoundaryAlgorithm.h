/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Manages PatchData objects created at coarse-fine stencils
 *
 ************************************************************************/

#ifndef included_xfer_CompositeBoundaryAlgorithm
#define included_xfer_CompositeBoundaryAlgorithm

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/xfer/CompositeBoundarySchedule.h"

#include <memory>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class CompositeBoundaryAlgorithm describes a pattern for
 * a communication and representation of data at coarse-fine boundaries.
 *
 * The purpose of this class is to allow an application code that is operating
 * on a particular Patch from a PatchHierarchy to have a local view into data
 * from the next finer level that exists within a certain stencil width of the
 * Patch's coarse-fine boundaries (if it has any).
 *
 * CompositeBoundaryAlgorithm holds the the desired stencil width for
 * a composite boundary description, and manages the patch data ids for
 * the data that is desired to be accessed at the composity boundaries.
 * The stencil width and data ids should be set once in the life of an
 * instance of this object.  Once its state is set, CompositeBoundaryAlglorithm
 * can repeatedly call createSchedule to create new CompositeBoundarySchedules
 * at different levels and at different times in the run of an application.
 *
 * @see CompositeBoundarySchedule 
 */

class CompositeBoundaryAlgorithm
{
public:

   /*!
    * @brief Constructor sets up basic state
    *
    * The constructor associates the CompositeBoundaryAlgorithm with a hierarchy
    * and a stencil width.
    *
    * The stencil width is interpreted in terms of the resolution of finer
    * levels.  So if the value is 4, that means that code operating on level n
    * will be able to see a coarse-fine stencil of width 4 in the mesh 
    * resolution of level n + 1;
    *
    * @param[in] hierarchy
    * @param[in] stencil_width
    *
    * @pre hierarchy != 0
    * @pre stencil_width >= 1
    */
   CompositeBoundaryAlgorithm(
      std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      int stencil_width);

   /*!
    * @brief Destructor
    */
   ~CompositeBoundaryAlgorithm();

   /*!
    * @brief Add a patch data id to identify data to be managed on the stencils.
    *
    * This may be called multiple times to add as many data ids as desired,
    * but all data ids must be added prior to calling createSchedule().
    *
    * @param[in] data_id
    *
    * @pre data_id >= 0 
    */
   void addDataId(int data_id);

   /*!
    * @brief Create a CompositeBoundarySchedule to communicate and manage
    * data at the coarse-fine boundaries
    *
    * The given level number represents a level from the hierarchy.  The
    * communicate data from the next finer level of resolution.
    *
    * This method should be called to create a new schedule after the level
    * represented by level_num and/or the next finer level have been regridded.
    *
    * @param level_num   level number for a level from hierarchy.
    */
   std::shared_ptr<CompositeBoundarySchedule>
   createSchedule(int level_num);

private:
   CompositeBoundaryAlgorithm(
      const CompositeBoundaryAlgorithm&);                  // not implemented
   CompositeBoundaryAlgorithm&
   operator = (
      const CompositeBoundaryAlgorithm&);                  // not implemented

   /*!
    * @brief Hierarchy where the stencils exist
    */
   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   /*!
    * @brief Width of the stencils.  The width represents a width in
    * terms of fine resolution at coarse-fine boundaries.
    */
   int d_stencil_width;

   /*!
    * @brief Patch data ids for data to be held on the stencil.
    */
   std::set<int> d_data_ids;

   /*!
    * @brief Flag telling that a schedule has been created.  This prevents
    * new data ids from being added.
    */
   bool d_schedule_created;

};

}
}

#endif
