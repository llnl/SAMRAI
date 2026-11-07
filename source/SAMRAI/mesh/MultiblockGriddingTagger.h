/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for refining AMR data.
 *
 ************************************************************************/

#ifndef included_mesh_MultiblockGriddingTagger
#define included_mesh_MultiblockGriddingTagger

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/SingularityPatchStrategy.h"

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Class MultiblockGriddingTagger is a concrete implementation
 * of RefinePatchStrategy that is used for boundary filling of tag
 * data.
 *
 * This class is needed for the calls to RefineSchedule in
 * the GriddingAlgorithm.
 *
 * This class implements the interface from SingularityPatchStrategy for
 * fillSingularityBoundaryConditions(), so that boundary conditions for
 * tag data that abuts a singularity can be properly filled.  Also
 * implemented are the interfaces for xfer::RefinePatchStrategy, needed
 * primarily for physical boundary filling.
 *
 * @see GriddingAlgorithm
 * @see xfer::RefineSchedule
 * @see xfer::RefinePatchStrategy
 */

class MultiblockGriddingTagger:
   public xfer::RefinePatchStrategy,
   public xfer::SingularityPatchStrategy
{
public:
   /*!
    * @brief The constructor does nothing interesting.
    */
   MultiblockGriddingTagger();

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~MultiblockGriddingTagger();

   /*!
    * @brief Set the patch data index for tag data.  This routine
    * must be called with a valid cell-centered integer patch data
    * index.
    */
   virtual void
   setScratchTagPatchDataIndex(
      int buf_tag_indx);

   /*!
    * @brief Physical boundary fill
    *
    * Implementation of interface defined in xfer::RefinePatchStrategy.
    * Fills ghost cells of patch data at physical boundaries.
    *
    * @param patch               Patch where data is stored
    * @param fill_time           Simulation time when filling occurs
    * @param ghost_width_to_fill Maximum ghost width of all data to be filled
    */
   virtual void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double fill_time,
      const hier::IntVector& ghost_width_to_fill);

   /*!
    * @brief Set the ghost data at a multiblock singularity.
    *
    * Implementation of interface defined in RefinePatchStrategy.
    * Fills ghost cells of patch data that abut multiblock singularities.
    * The encon_level contains the data from neighboring blocks that also
    * also abut the singularity, and that data from the neighbors
    * is used to fill data on the local patch.
    *
    * @param patch               Local patch containing data to be filled
    * @param encon_level  Level representing enhanced connectivity ghost
    *                     regions
    * @param dst_to_encon  Connector from destination level to encon_level
    * @param fill_box             All ghost data to be filled will be within
    *                             this box
    * @param boundary_box         BoundaryBox object that stores information
    *                             about the type and location of the boundary
    *                             where ghost cells will be filled
    * @param grid_geometry
    *
    * @pre (patch.getDim() == fill_box.getDim()) &&
    *      (patch.getDim() == boundary_box.getDim())
    * @pre !grid_geometry->hasEnhancedConnectivity() || dst_to_encon
    */
   virtual void
   fillSingularityBoundaryConditions(
      hier::Patch& patch,
      const hier::PatchLevel& encon_level,
      std::shared_ptr<const hier::Connector> dst_to_encon,
      const hier::Box& fill_box,
      const hier::BoundaryBox& boundary_box,
      const std::shared_ptr<hier::BaseGridGeometry>& grid_geometry);

   /*!
    * @brief Return maximum stencil width needed for user-defined
    * data interpolation operations.  This is needed to
    * determine the correct interpolation data dependencies.
    *
    * Always returns an IntVector of ones, because that is the maximum
    * stencil needed for the operations in GriddingAlgorithm
    */
   virtual hier::IntVector
   getRefineOpStencilWidth(
      const tbox::Dimension& dim) const;

   /*!
    * Perform user-defined refining operations.  This member function
    * is called before standard refining operations (expressed using
    * concrete subclasses of the hier::RefineOperator base class).
    * The preprocess function must refine data from the scratch components
    * of the coarse patch into the scratch components of the fine patch on the
    * specified fine box region.  Recall that the scratch components are
    * specified in calls to the registerRefine() function in the
    * xfer::RefineAlgorithm class.
    *
    * @param fine        Fine patch containing destination data.
    * @param coarse      Coarse patch containing source data.
    * @param fine_box    Box region on fine patch into which data is refined.
    * @param ratio       Integer vector containing ratio relating index space
    *                    between coarse and fine patches.
    */
   virtual void
   preprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);

   /*!
    * Perform user-defined refining operations.  This member function
    * is called before standard refining operations (expressed using
    * concrete subclasses of the hier::RefineOperator base class).
    * The postprocess function must refine data from the scratch components
    * of the coarse patch into the scratch components of the fine patch on the
    * specified fine box region.  Recall that the scratch components are
    * specified in calls to the registerRefine() function in the
    * xfer::RefineAlgorithm class.
    *
    * @param fine        Fine patch containing destination data.
    * @param coarse      Coarse patch containing source data.
    * @param fine_box    Box region on fine patch into which data is refined.
    * @param ratio       Integer vector containing ratio relating index space
    *                    between coarse and fine patches.
    *
    * @pre (fine.getDim() == coarse.getDim()) &&
    *      (fine.getDim() == fine_box.getDim()) &&
    *      (fine.getDim() == ratio.getDim())
    */
   virtual void
   postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);

private:
   /*
    * Patch data index for
    */
   int d_buf_tag_indx;

};

}
}

#endif
