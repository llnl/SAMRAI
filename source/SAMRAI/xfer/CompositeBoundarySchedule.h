/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Manages PatchData objects created at coarse-fine stencils
 *
 ************************************************************************/

#ifndef included_xfer_CompositeBoundarySchedule
#define included_xfer_CompositeBoundarySchedule

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"

#include <memory>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class CompositeBoundarySchedule manages the communication and access
 * of patch data objects that hold data at a composite boundary.
 *
 * The purpose of this class is to allow an application code that is operating
 * on a particular Patch from a PatchHierarchy to have a local view into data
 * from the next finer level that exists within a certain stencil width of the
 * Patch's coarse-fine boundaries (if it has any).
 *
 * The recommended usage is to create a CompositeBoundaryAlgorithm to hold
 * the stencil width and the patch data ids, and use
 * CompositeBoundaryAlgorithm::createSchedule to create instances of this
 * class.  
 */

class CompositeBoundarySchedule
{
public:

   /*!
    * @brief Constructor sets up basic state
    *
    * The constructor associates the CompositeBoundarySchedule with a hierarchy
    * and a stencil width.  It does not manage any data until further
    * setup calls are made.
    *
    * The stencil width is interpreted in terms of the resolution of finer
    * levels.  So if the value is 4, that means that code operating on level n
    * will be able to see a coarse-fine stencil of width 4 in the mesh 
    * resolution of level n + 1;
    *
    * This object will become invalid if either the coarse level or the
    * next finer level is changed by a regrid of the hierarchy.
    *
    * @param[in] hierarchy
    * @param[in] coarse_level_num  Level number of the coarse level at
    *                              coarse-fine boundary
    * @param[in] stencil_width  Width in terms of fine level zones
    * @param[in] data_ids   Patch data ids for data to be communicated 
    *
    * @pre hierarchy != 0
    * @pre stencil_width >= 1
    */
   CompositeBoundarySchedule(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      int coarse_level_num,
      int stencil_width,
      const std::set<int>& data_ids);

   /*!
    * @brief Destructor
    */
   ~CompositeBoundarySchedule();

   /*!
    * @brief Fill data into the stencil that exists at the coarse-fine
    * boundaries of the given level.
    *
    * Data willl be filled from the corresponding fine level from the
    * onto the stencil, where it is then accessible through the method
    * getBoundaryPatchData().
    *
    * This may be called multiple times, to account for changing state of
    * data on the hierarchy.
    *
    * @param[in] fill_time   Time for the filling operation.
    */
   void fillData(double fill_time);

   /*!
    * @brief Get the uncovered boxes for a patch on the schedule's coarse level.
    *
    * The uncovered boxes represent the regions of the patch, if any, that
    * are not overlapped by the mesh of the next finest level.
    *
    * @param patch  A local patch on the the coarse level.  An error will occur
    * the patch is not local or not on the level.
    */
   const hier::BoxContainer&
   getUncoveredBoxes(
      const hier::Patch& patch) const
   {
      std::map<hier::BoxId, hier::BoxContainer>::const_iterator itr =
         d_patch_to_uncovered.find(patch.getBox().getBoxId());

      if (itr == d_patch_to_uncovered.end()) {
         TBOX_ERROR("CompositeBoundarySchedule::getUncoveredBoxes error: Patch " << patch.getBox().getBoxId() << " is not a local patch on the current coarse level." << std::endl);
      }

      return itr->second;
   }

   /*!
    * @brief Get PatchData from the composite boundary stencil for a given
    * data id.
    *
    * Get a vector that holds the stencil PatchData for the parts of the
    * coarse-fine boundary touching the Patch.  The Patch must be a local
    * Patch on the hierarchy, and fillData() must already have been called on
    * this object.
    *
    * If the Patch does not touch any coarse-fine boundaries, the returned
    * vector will be empty.
    *
    * @param patch    a local Patch from the hierarchy
    * @param data_id  patch data id for data stored in the stencil
    */  
   const std::vector<std::shared_ptr<hier::PatchData> >&
   getBoundaryPatchData(
      const hier::Patch& patch,
      int data_id) const;

   /*!
    * @brief Get PatchData from the composite boundary stencil associated with
    * an uncovered box.
    *
    * The returned vector will contain PatchData located at the coarse-fine
    * stencil touching the uncovered box.
    *
    * An error will occur if the uncovered_id is not the id of an uncovered
    * box, or if the data_id is not an id that was passed into the constructor.
    *
    * @param uncovered_id  BoxId for an uncovered box
    * @param data_id       patch data id for the desired data.
    */
   std::vector<std::shared_ptr<hier::PatchData> >
   getDataForUncoveredBox(
      const hier::BoxId& uncovered_id,
      int data_id) const;

private:
   CompositeBoundarySchedule(
      const CompositeBoundarySchedule&);                  // not implemented
   CompositeBoundarySchedule&
   operator = (
      const CompositeBoundarySchedule&);                  // not implemented

   /*!
    * @brief Register the data ids for communication.
    *
    * This must be called prior to calling createStencilForLevel().
    */
   void registerDataIds();

   /*!
    * @brief Create a stencil that exists at the coarse-fine boundaries
    * of the coarse level.
    */
   void createStencilForLevel();

   /*!
    * @brief Set up communication schedule to fill data at composite boundary
    *
    * This sets up the communications that fill the patch data around the
    * coarse-fine boundaries of the coarse level.
    */
   void setUpSchedule();

   /*!
    * @brief Hierarchy where the stencils exist
    */
   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   /*!
    * @brief Levels representing the stencils
    */
   std::shared_ptr<hier::PatchLevel> d_stencil_level;

   /*!
    * @brief Map from BoxId of a hierarchy patch to the BoxIds of its
    * coarse-fine stencil boxes.
    */
   std::map<hier::BoxId, std::set<hier::BoxId> > d_stencil_map;

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
    * @brief Algorithm for communication of data to the stencil.
    */
   RefineAlgorithm d_refine_algorithm;

   /*!
    * @brief Schedules for communication of data from hierarchy to stencils.
    */
   std::shared_ptr<RefineSchedule> d_refine_schedule;

   /*
    * typedefs to simplify nested container declaration.
    */
   typedef std::vector<std::shared_ptr<hier::PatchData> > VectorPatchData;
   typedef std::map<hier::BoxId, VectorPatchData> MapBoxIdPatchData;

   // map BoxId from Patch to container of uncovered boxes
   std::map<hier::BoxId, hier::BoxContainer> d_patch_to_uncovered;

   // map BoxId of uncovered box to container of stencil boxes 
   std::map<hier::BoxId, std::set<hier::BoxId> > d_uncovered_to_stencil;

   // outer map maps a data id to the inner map.  Inner map maps BoxId of
   // a stencil box to the PatchData located on that box.
   std::map<int, std::map<hier::BoxId, std::shared_ptr<hier::PatchData> > >
      d_stencil_to_data;
   
   /*!
    * @brief Container of PatchData on the stencil.
    *
    * The inner nested container is a vector of PatchData which consists of 
    * one PatchData object for each box of the coarse-fine stencil for a single
    * Patch of the hierarchy.
    *
    * The nesting of the containers is organized as:
    *
    * d_data_map[data_id][box_id][stencil_box_index]
    *
    * level_number is the level number of the coarse level at a coarse-fine
    * boundary where the stencil has been created and filled.  data_id
    * is the patch data id of the data that is to be accessed.  box_id is
    * the id of a Patch on the hierarchy.  stencil_box_index is an index into
    * the inner nested vector of PatchData.
    */
   std::map<int, MapBoxIdPatchData> d_data_map;

   /*!
    * @brief Vector of flags telling if stencil has been created for each
    * level of the hierarchy.
    */
   bool d_stencil_created;

   /*!
    * @brief Vector of flags telling if stencil has been filled for each
    * level of the hierarchy.
    */
   bool d_stencil_filled; 

   /*!
    * The level number of the coarse level at the coarse-fine boundary where
    * the composite boundary stencils are created.
    */
   int d_coarse_level_num;
};

}
}

#endif
