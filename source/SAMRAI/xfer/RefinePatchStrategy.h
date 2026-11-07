/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for refining AMR data.
 *
 ************************************************************************/

#ifndef included_xfer_RefinePatchStrategy
#define included_xfer_RefinePatchStrategy

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/xfer/RefineClasses.h"

#include <set>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Abstract base class for user-defined patch data refining operations
 * and physical boundary filling operations.
 *
 * RefinePatchStrategy provides an interface for a user to supply methods
 * for application-specific refining of data between levels in an
 * AMR patch hierarchy and filling of physical boundary values.
 * A concrete subclass must define four member functions to perform
 * the following tasks:
 *
 * <ul>
 *    <li> define maximum stencil width for user-defined refine operations
 *    <li> fill physical boundary conditions
 *    <li> preprocess the interpolation
 *    <li> postprocess the interpolation
 * </ul>
 *
 * Note that the preprocess member function is called before standard data
 * refine using RefineOperators and the postprocessor member function is called
 * afterwards.
 *
 * There are two versions of the preprocess and postprocess functions.  The
 * default abstract function only takes a single box.  The user may also
 * over-ride preprocess and postprocess functions that take a box list and
 * process an entire patch at one time.  By default, the box list version
 * loops over all of the boxes in the box list and calls the single box version.
 *
 * For operations on a mesh with block-singularities, implementations of this
 * class must also inherit and implement SingularityPatchStrategy.
 *
 * @see RefineAlgorithm
 * @see RefineSchedule
 */

class RefinePatchStrategy
{
public:
   /*!
    * @brief Get the maximum stencil width over all RefinePatchStrategy objects
    * used in an application.
    *
    * @return The maximum of the return values of calls to
    * getRefineOpStencilWidth() for every RefinePatchStrategy of the
    * given Dimension used in an application.
    *
    * @param[in] dim
    */
   static hier::IntVector
   getMaxRefineOpStencilWidth(
      const tbox::Dimension& dim);

   /*!
    * @brief Constructor.
    *
    * The constructor will register the constructed object with a static
    * set that manages all RefinePatchStrategy objects in an application.
    */
   RefinePatchStrategy();

   /*!
    * @brief Destructor
    */
   virtual ~RefinePatchStrategy();

   /*!
    * @brief Pure virtual function interface for physical boundary filling.
    *
    * Set data in ghost regions at patch boundaries that touch the
    * physical domain boundary.  The specific data values set in physical
    * boundary ghost regions are determined by the boundary conditions needed
    * by the user application.  The patch data components that should be set
    * in this function correspond to the "scratch" components specified in calls
    * to the registerRefine() function in the RefineAlgorithm class.
    *
    * @param[out] patch               Patch on which to fill boundary data.
    * @param[in] fill_time            Simulation time for boundary filling.
    * @param[in] ghost_width_to_fill  Maximum ghost width to fill over
    *                                 all registered scratch components.
    */
   virtual void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double fill_time,
      const hier::IntVector& ghost_width_to_fill) = 0;

   /*!
    * @brief Return maximum stencil width needed for user-defined
    * data refinement operations performed by this object.
    *
    * This is needed to determine the correct interpolatin data dependencies
    * and to ensure that the data has a sufficient amount of ghost width.
    *
    * For any user-defined interpolation operations implemented in the
    * preprocess or postprocess methods, return the maximum stencil needed
    * on a coarse patch to refine data to a fine patch.
    * If your implementation doesn't work with the given dimension, return
    * zero.
    */
   virtual hier::IntVector
   getRefineOpStencilWidth(
      const tbox::Dimension& dim) const = 0;

   /*!
    * @brief Perform user-defined patch data refinement operations.
    *
    * This member function is called before standard refine operations
    * (expressed using concrete subclasses of the RefineOperator base class).
    * The preprocess function must refine data from the scratch components
    * on the coarse patch into the scratch components of the fine patch
    * on the specified fine box region.  Recall that the scratch components
    * are specified in calls to the registerRefine() function in the
    * RefineAlgorithm class.
    *
    * @param[out] fine     Fine patch containing destination data.
    * @param[in] coarse    Coarse patch containing source data.
    * @param[in] fine_box  Box region on fine patch into which data is refined.
    * @param[in] ratio     Refinement ratio between coarse and fine patches.
    */
   virtual void
   preprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio) = 0;

   /*!
    * @brief Perform user-defined patch data refinement operations.
    *
    * This member function is called after standard refine operations
    * (expressed using concrete subclasses of the RefineOperator base class).
    * The preprocess function must refine data from the scratch components
    * on the coarse patch into the scratch components of the fine patch
    * on the specified fine box region.  Recall that the scratch components
    * are specified in calls to the registerRefine() function in the
    * RefineAlgorithm class.
    *
    * @param[out] fine     Fine patch containing destination data.
    * @param[in] coarse    Coarse patch containing source data.
    * @param[in] fine_box  Box region on fine patch into which data is refined.     * @param[in] ratio     Refinement ratio between coarse and fine patches.
    */
   virtual void
   postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio) = 0;

   /*!
    * Perform user-defined patch data refinement operations on a list of boxes.
    *
    * This member function is called before standard refining operations
    * (expressed using concrete subclasses of the RefineOperator base class).
    *
    * The default implementation of this virtual function loops over the
    * box list and calls the preprocessRefine() method for a single box.
    *
    * @param[out] fine   Fine patch containing destination data.
    * @param[in] coarse  Coarse patch containing source data.
    * @param[in] fine_boxes  List of box regions on fine patch into which data
    *                        is refined.
    * @param[in] ratio     Refinement ratio between coarse and fine patches.
    *
    * @pre (fine.getDim() == coarse.getDim()) &&
    *      (fine.getDim() == ratio.getDim())
    */
   virtual void
   preprocessRefineBoxes(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::BoxContainer& fine_boxes,
      const hier::IntVector& ratio)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY3(fine, coarse, ratio);
      for (hier::BoxContainer::const_iterator b = fine_boxes.begin();
           b != fine_boxes.end(); ++b) {
         preprocessRefine(fine, coarse, *b, ratio);
      }
   }

   /*!
    * Perform user-defined patch data refinement operations on a list of boxes.
    *
    * This member function is called after standard refining operations
    * (expressed using concrete subclasses of the RefineOperator base class).
    *
    * The default implementation of this virtual function loops over the
    * box list and calls the postprocessRefine() method for a single box.
    *
    * @param[out] fine   Fine patch containing destination data.
    * @param[in] coarse  Coarse patch containing source data.
    * @param[in] fine_boxes  List of box regions on fine patch into which data
    *                        is refined.
    * @param[in] ratio     Refinement ratio between coarse and fine patches.
    *
    * @pre (fine.getDim() == coarse.getDim()) &&
    *      (fine.getDim() == ratio.getDim())
    */
   virtual void
   postprocessRefineBoxes(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::BoxContainer& fine_boxes,
      const hier::IntVector& ratio)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY3(fine, coarse, ratio);
      for (hier::BoxContainer::const_iterator b = fine_boxes.begin();
           b != fine_boxes.end(); ++b) {
         postprocessRefine(fine, coarse, *b, ratio);
      }
   }

   /*!
    * @brief Perform user-defined patch data pre-refinement operations
    * on a PatchLevel.
    *
    * This member function is called before standard refine operations
    * (expressed using concrete subclasses of the RefineOperator base class).
    *
    * @param[out] fine_level     Fine patch level containing destination data.
    * @param[in] coarse_level    Coarse patch level containing source data.
    * @param[in] coarse_to_fine      Connector coarse to fine
    * @param[in] coarse_to_unfilled  Connector coarse to level representing
    *                                boxes that need to be filled.
    * @param[in] overlaps
    * @param[in] refine_items
    */
   virtual void
   preprocessRefineLevel(
      hier::PatchLevel& fine_level,
      const hier::PatchLevel& coarse_level,
      const hier::Connector& coarse_to_fine,
      const hier::Connector& coarse_to_unfilled,
      const std::vector<std::vector<std::shared_ptr<hier::BoxOverlap> > >& overlaps,
      const RefineClasses::Data** refine_items) {
      NULL_USE(fine_level);
      NULL_USE(coarse_level);
      NULL_USE(coarse_to_fine);
      NULL_USE(coarse_to_unfilled);
      NULL_USE(overlaps);
      NULL_USE(refine_items);
      setNeedRefineSynchronize(false);
   }

   /*!
    * @brief Perform user-defined patch data post-refinement
    * operations on a PatchLevel.
    *
    * This member function is called after standard refine operations
    * (expressed using concrete subclasses of the RefineOperator base class).
    *
    * @param[out] fine_level     Fine patch level containing destination data.
    * @param[in] coarse_level    Coarse patch level containing source data.
    * @param[in] coarse_to_fine      Connector coarse to fine
    * @param[in] coarse_to_unfilled  Connector coarse to level representing
    *                                boxes that need to be filled.
    */
   virtual void
   postprocessRefineLevel(
      hier::PatchLevel& fine_level,
      const hier::PatchLevel& coarse_level,
      const hier::Connector& coarse_to_fine,
      const hier::Connector& coarse_to_unfilled) {
      NULL_USE(fine_level);
      NULL_USE(coarse_level);
      NULL_USE(coarse_to_fine);
      NULL_USE(coarse_to_unfilled);
      setNeedRefineSynchronize(false);
   }

   virtual void
   setPostRefineSyncFlag()
   {
      setNeedRefineSynchronize(true);
   }

   /*!
    * @brief Check flag for if host-device synchronization is needed.
    *
    * Returns current value of the flag while setting the flag back to
    * the default value of true.
    */
   bool
   needSynchronize()
   {
      bool flag = d_need_synchronize;
      d_need_synchronize = true;
      return flag;
   }

protected:

   /*!
    * @brief Set flag indicating if device synchronization is needed after
    * a child class operation.
    *
    * This allows implementations of methods such as preprocessRefine and
    * postprocessRefine to set the flag to false if they have done nothing
    * that requires host-device synchronization and do not need
    * RefineSchedule to call the synchronize routine.
    */
   void
   setNeedRefineSynchronize(bool flag)
   {
      d_need_synchronize = flag;
   }

private:
   /*!
    * @brief Get the set of RefinePatchStrategy objects that have been
    * registered.
    */
   static std::set<RefinePatchStrategy *>&
   getCurrentObjects()
   {
      static std::set<RefinePatchStrategy *> current_objects;
      return current_objects;
   }

   /*!
    * @brief Register the object with a set of all RefinePatchStrategy
    * objects used in an application.
    */
   void
   registerObject()
   {
      std::set<RefinePatchStrategy *>& current_objects =
         RefinePatchStrategy::getCurrentObjects();
      current_objects.insert(this);
   }

   /*!
    * @brief Unregister the object.
    */
   void
   unregisterObject()
   {
      std::set<RefinePatchStrategy *>& current_objects =
         RefinePatchStrategy::getCurrentObjects();
      current_objects.erase(this);
   }

   bool d_need_synchronize = true;

};

}
}

#endif
