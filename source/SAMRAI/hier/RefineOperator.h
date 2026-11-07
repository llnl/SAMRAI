/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract base class for spatial refinement operators.
 *
 ************************************************************************/

#ifndef included_hier_RefineOperator
#define included_hier_RefineOperator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/OpenMPUtilities.h"

#include <string>
#include <map>

namespace SAMRAI {
namespace hier {

/**
 * Class RefineOperator is an abstract base class for each
 * spatial refinement operator used in the SAMRAI framework.  This class
 * defines the interface between numerical refinement routines and the
 * rest of the framework.  Each concrete refinement operator subclass
 * must provide three four operations:
 *
 *
 *
 * - \b (1) {an implementation of the refinement operation
 *            appropriate for its corresponding patch data type.}
 * - \b (2) {a function that determines whether or not the operator
 *             matches an arbitrary request for a refinement operator.}
 * - \b (3) {a function that returns the stencil width of the operator
 *             (i.e., the number of ghost cells needed by the operator).}
 * - \b (4) {a function that returns an integer stating the priority of the
 *             operator with respect to other refinement operators.}
 *
 *
 *
 * To add a new refinement operator (either for a new patch data type
 * or for a new time refinement routine on an existing type), define
 * the operator by inheriting from this abstract base class.  The operator
 * subclass must implement the refinement operation in the refine()
 * fnction.  The stencil width and operator priority must be returned
 * from the getStencilWidth() and getOperatorPriority() functions,
 * respectively.  Then, the new operator must be added to the
 * operator list for the appropriate transfer geometry object using the
 * BaseGridGeometry::addRefineOperator() function.
 *
 * Since spatial refinement operators usually depend on patch data centering
 * and data type as well as the mesh coordinate system, they are defined
 * in the <EM>geometry</EM> package.
 *
 * @see TransferOperatorRegistry
 */

class RefineOperator
{
public:
   /*!
    * @brief Construct the object with a name to allow the
    * TransferOperatorRegistry class to look up the object using a
    * string.
    *
    * The constructor must be given a name.  The object will be
    * registered under this name with the TransferOperatorRegistry class.
    * The name must be unique, as duplicate names are not allowed.
    */
   RefineOperator(
      const std::string& name);

   /**
    * The virtual destructor for the refinement operator does
    * nothing interesting.
    */
   virtual ~RefineOperator();

   /**
    * Return name std::string identifier of the refinement operation.
    */
   const std::string&
   getOperatorName() const
   {
      return d_name;
   }

   /**
    * Return the priority of this operator relative to other refinement
    * operators.  The SAMRAI transfer routines guarantee that refinement
    * using operators with lower priority values will be performed before
    * those with higher priority.
    */
   virtual int
   getOperatorPriority() const = 0;

   /**
    * Return the stencil width associated with the refinement operator.
    * The SAMRAI transfer routines guarantee that the source patch will
    * contain sufficient ghost cell data surrounding the interior to
    * satisfy the stencil width requirements for each refinement operator.
    * If your implementation doesn't work with the given dimension, return
    * zero.
    */
   virtual IntVector
   getStencilWidth(
      const tbox::Dimension& dim) const = 0;

   /**
    * Refine the source component on the coarse patch to the destination
    * component on the fine patch. The refinement operation is performed
    * on the intersection of the destination patch and the boxes contained
    * in fine_overlap. The coarse patch is guaranteed to contain sufficient
    * data for the stencil width of the refinement operator.
    */
   virtual void
   refine(
      Patch& fine,
      const Patch& coarse,
      const int dst_component,
      const int src_component,
      const BoxOverlap& fine_overlap,
      const IntVector& ratio) const = 0;

   /*!
    * @brief Get the max stencil width of all refine operators.
    *
    * The max stencil width computed from all registered (constructed)
    * refine operators.
    */
   static IntVector
   getMaxRefineOpStencilWidth(
      const tbox::Dimension& dim);

private:
   RefineOperator(
      const RefineOperator&) = delete;
   RefineOperator&
   operator = (
      const RefineOperator&) = delete;

   /*
    * TODO SGS Rich has better way of doing this.
    */

   /*!
    * @brief Associate the given name with this operator.
    *
    * Registering an operator with a name allows that operator
    * to be looked up by name.  Operators registered must have
    * unique names.
    */
   void
   registerInLookupTable(
      const std::string& name);

   /*!
    * @brief Remove the operator with the given name.
    *
    */
   void
   removeFromLookupTable(
      const std::string& name);

   /*!
    * @brief Method registered with ShutdownManager to initialize statics.
    */
   static void
   initializeCallback();

   /*!
    * @brief Method registered with ShutdownManager to cleanup statics.
    */
   static void
   finalizeCallback();

   const std::string d_name;

   static std::multimap<std::string, RefineOperator *> s_lookup_table;
   static TBOX_omp_lock_t l_lookup_table;

   static tbox::StartupShutdownManager::Handler
      s_finalize_handler;

};

}
}

#endif
