/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for coarsening AMR data.
 *
 ************************************************************************/

#ifndef included_xfer_CoarsenPatchStrategy
#define included_xfer_CoarsenPatchStrategy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchLevel.h"

#include <set>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Abstract base class for user-defined patch data coarsening operations.
 *
 * CoarsenPatchStrategy provides an interface for a user to supply
 * methods for application-specific coarsening of data between two levels in
 * an AMR patch hierarchy.  A concrete subclass must define three member
 * functions to perform the following tasks:
 *
 * <ul>
 *    <li> define maximum stencil width for user-defined coarsen operations
 *    <li> preprocess the coarsening
 *    <li> postprocess the coarsening
 * </ul>
 *
 * Note that the preprocess member function is called before standard data
 * coarsening using CoarsenOperators and the postprocess member function is
 * called afterwards.
 *
 * @see CoarsenAlgorithm
 * @see CoarsenSchedule
 */

class CoarsenPatchStrategy
{
public:
   /*!
    * @brief Get the maximum stencil width over all CoarsenPatchStrategy objects
    * used in an application.
    *
    * @return The maximum of the return values of calls to
    * getCoarsenOpStencilWidth() for every CoarsenPatchStrategy of the
    * given Dimension used in an application.
    *
    * @param[in] dim
    */
   static hier::IntVector
   getMaxCoarsenOpStencilWidth(
      const tbox::Dimension& dim);

   /*!
    * @brief Constructor.
    *
    * The constructor will register the constructed object with a static
    * set that manages all CoarsenPatchStrategy objects in an application.
    */
   CoarsenPatchStrategy();

   /*!
    * @brief Destructor.
    */
   virtual ~CoarsenPatchStrategy();

   /*!
    * @brief Return maximum stencil width needed for user-defined
    * data coarsening operations performed by this object.
    *
    * This is needed to determine the correct coarsening data dependencies and
    * to ensure that the data has a sufficient amount of ghost width.
    *
    * For any user-defined coarsening operations implemented in the
    * preprocess or postprocess methods, return the maximum stencil needed
    * on a fine patch to coarsen data to a coarse patch.
    * If your implementation doesn't work with the given dimension, return
    * zero.
    */
   virtual hier::IntVector
   getCoarsenOpStencilWidth(
      const tbox::Dimension& dim) const = 0;

   /*!
    * @brief Perform user-defined patch data coarsening operations.
    *
    * This member function is called before standard coarsening operations
    * (expressed using concrete subclasses of the CoarsenOperator base class).
    * The preprocess function should move data from the source components
    * on the fine patch into the source components on the coarse patch
    * in the specified coarse box region.  Recall that the source components
    * are specified in calls to the registerCoarsen() function in the
    * CoarsenAlgorithm class.
    *
    * @param[out] coarse Coarse patch that will receive coarsened data.
    * @param[in] fine    Fine patch containing source data.
    * @param[in] coarse_box  Box region on coarse patch into which data is
    *                        coarsened.
    * @param[in] ratio   Refinement ratio between coarse and fine patches.
    */
   virtual void
   preprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) = 0;

   /*!
    * @brief Perform user-defined patch data coarsening operations.
    *
    * This member function is called after standard coarsening operations
    * (expressed using concrete subclasses of the CoarsenOperator base class).
    * The postprocess function should move data from the source components on
    * the fine patch into the source components on the coarse patch in the
    * specified coarse box region.  Recall that the source components are
    * specified in calls to the registerCoarsen() function in the
    * CoarsenAlgorithm class.
    *
    * @param[out] coarse  Coarse patch that will receive coarsened data.
    * @param[in] fine     Fine patch containing source data.
    * @param[in] coarse_box  hier::Box region on coarse patch into which data
    *                        is coarsened.
    * @param[in] ratio    Refinement ratio between coarse and fine patches.
    */
   virtual void
   postprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) = 0;


   virtual void
   preprocessCoarsenLevel(
      hier::PatchLevel& coarse_level,
      const hier::PatchLevel& fine_level) {
      NULL_USE(coarse_level);
      NULL_USE(fine_level);
      setNeedCoarsenSynchronize(false);
   }


   virtual void
   postprocessCoarsenLevel(
      hier::PatchLevel& coarse_level,
      const hier::PatchLevel& fine_level) {
      NULL_USE(coarse_level);
      NULL_USE(fine_level);
      setNeedCoarsenSynchronize(false);
   }


   virtual void
   setPostCoarsenSyncFlag()
   {
      setNeedCoarsenSynchronize(true);
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
    * This allows implementations of methods such as preprocessCoarsen and
    * postprocessCoarsen to set the flag to false if they have done nothing
    * that requires host-device synchronization and do not need
    * CoarsenSchedule to call the synchronize routine.
    */
   void
   setNeedCoarsenSynchronize(bool flag)
   {
      d_need_synchronize = flag;
   }

private:
   /*!
    * @brief Get the set of CoarsenPatchStrategy objects that have been
    * registered.
    */
   static std::set<CoarsenPatchStrategy *>&
   getCurrentObjects()
   {
      static std::set<CoarsenPatchStrategy *> current_objects;
      return current_objects;
   }

   /*!
    * @brief Register the object with a set of all CoarsenPatchStrategy
    * objects used in an application.
    */
   void
   registerObject()
   {
      std::set<CoarsenPatchStrategy *>& current_objects =
         CoarsenPatchStrategy::getCurrentObjects();
      current_objects.insert(this);
   }

   bool d_need_synchronize = true;

};

}
}

#endif
