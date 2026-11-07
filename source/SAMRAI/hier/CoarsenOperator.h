/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract base class for spatial coarsening operators.
 *
 ************************************************************************/

#ifndef included_hier_CoarsenOperator
#define included_hier_CoarsenOperator

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/OpenMPUtilities.h"

#include <string>
#include <map>
#include <memory>

namespace SAMRAI {
namespace hier {

/**
 * Class CoarsenOperator is an abstract base class for each
 * spatial coarsening operator used in the SAMRAI framework.  This class
 * defines the interface between numerical coarsening routines and the
 * rest of the framework.  Each concrete coarsening operator subclass
 * must provide four operations:
 *
 *
 *
 * - \b (1) {an implementation of the coarsening operation
 *            appropriate for its corresponding patch data type.}
 * - \b (2) {a function that determines whether or not the operator
 *             matches an arbitrary request for a coarsening operator.}
 * - \b (3) {a function that returns the stencil width of the operator
 *             (i.e., the number of ghost cells needed by the operator).}
 * - \b (4) {a function that returns an integer stating the priority of the
 *             operator with respect to other coarsening operators.}
 *
 *
 *
 * To add a new coarsening operator (either for a new patch data type
 * or for a new time coarsening routine on an existing type), define
 * the operator by inheriting from this abstract base class.  The operator
 * subclass must implement the coarsening operation in the coarsen()
 * function.  The stencil width and operator priority must be returned
 * from the getStencilWidth() and getOperatorPriority() functions,
 * respectively.  Then, the new operator must be added to the
 * operator list for the appropriate transfer geometry object using the
 * BaseGridGeometry::addCarsenOperator() function.
 *
 * Since spatial coarsening operators usually depend on patch data centering
 * and data type as well as the mesh coordinate system, they are defined
 * in the <EM>geometry</EM> package.
 *
 * @see TransferOperatorRegistry
 */

class CoarsenOperator
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
   CoarsenOperator(
      const std::string& name);

   /**
    * The virtual destructor for the coarsening operator does
    * nothing interesting.
    */
   virtual ~CoarsenOperator();

   /**
    * Return name std::string identifier of the coarsening operation.
    */
   const std::string&
   getOperatorName() const
   {
      return d_name;
   }

   /**
    * Return the priority of this operator relative to other coarsening
    * operators.  The SAMRAI transfer routines guarantee that coarsening
    * using operators with lower priority values will be performed before
    * those with higher priority.
    */
   virtual int
   getOperatorPriority() const = 0;

   /**
    * Return the stencil width associated with the coarsening operator.
    * The SAMRAI transfer routines guarantee that the source patch will
    * contain sufficient ghost cell data surrounding the interior to
    * satisfy the stencil width requirements for each coarsening operator.
    * If your implementation doesn't work with the given dimension, return
    * zero.
    */
   virtual IntVector
   getStencilWidth(
      const tbox::Dimension& dim) const = 0;

   /**
    * Coarsen the source component on the fine patch to the destination
    * component on the coarse patch. The coarsening operation is performed
    * on the intersection of the destination patch and the coarse box.
    * The fine patch is guaranteed to contain sufficient data for the
    * stencil width of the coarsening operator.
    */
   virtual void
   coarsen(
      Patch& coarse,
      const Patch& fine,
      const int dst_component,
      const int src_component,
      const Box& coarse_box,
      const IntVector& ratio) const = 0;

   /*!
    * @brief Get the max stencil width of all coarsen operators.
    *
    * The max stencil width computed from all registered (constructed)
    * coarsen operators.
    */
   static IntVector
   getMaxCoarsenOpStencilWidth(
      const tbox::Dimension& dim);

private:
   CoarsenOperator(
      const CoarsenOperator&) = delete;
   CoarsenOperator&
   operator = (
      const CoarsenOperator&) = delete;

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
    * @brief Method registered with ShutdownRegister to initialize statics.
    */
   static void
   initializeCallback();

   /*!
    * @brief Method registered with ShutdownRegister to cleanup statics.
    */
   static void
   finalizeCallback();

   const std::string d_name;

   static std::multimap<std::string, CoarsenOperator *> s_lookup_table;
   static TBOX_omp_lock_t l_lookup_table;

   static tbox::StartupShutdownManager::Handler s_finalize_handler;

};

}
}

#endif
