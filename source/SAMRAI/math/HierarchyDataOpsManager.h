/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Singleton manager for hierarchy data operation objects.
 *
 ************************************************************************/

#ifndef included_math_HierarchyDataOpsManager
#define included_math_HierarchyDataOpsManager

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/math/HierarchyDataOpsComplex.h"
#include "SAMRAI/math/HierarchyDataOpsInteger.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/Utilities.h"

#include <vector>
#include <memory>

namespace SAMRAI {
namespace math {

/**
 * Class HierarchyDataOpsManager is a Singleton class that serves as
 * a single point of control for gaining access to operation objects which
 * to treat patch data over multiple levels in an AMR hierarchy. Specifically,
 * this manager returns a pointer to a hierarchy data operation object given
 * a variable and a hierarchy.  Currently, this manager class supports patch
 * data operation objects for cell-, face-, and node-centered data of type
 * double, float, integer and dcomplex.  Although it may be used to manage
 * data operations for multiple patch hierarchies, it supports only one data
 * operations per variable centering and data type per hierarchy.  The manager
 * will create and return a default hierarchy operation object when queried
 * for an operation object if the operation does not already exist.  If a
 * matching operation object exists and is known to the manager, that operation
 * object will be returned.  This manager class may be extended using class
 * inheritance and implementing a manager subclass.  So, for example, if
 * more than one operation object per variable type is needed on a hierarchy,
 * such inheritance should be used to add this capability to the manager.
 *
 * Important note: This manager class is incomplete. Specifically, one cannot
 *                 set an operation object to over-ride the default .  This
 *                 feature will be added in the near future.
 *
 * Important note: If the manager must construct a new hierarchy data
 *                 operation object, the range of levels used in the
 *                 operation object must be set explicitly before the
 *                 operations can be used.
 *
 * See the Design Patterns book by Gamma {\em et al.} for more information
 * about the singleton pattern.
 *
 * @see HierarchyDataOpsComplex
 * @see HierarchyDataOpsInteger
 * @see HierarchyDataOpsReal
 */

class HierarchyDataOpsManager
{
public:
   /**
    * Return a pointer to the single instance of the patch data operation
    * manager.  All access to the HierarchyDataOpsManager object is
    * through the getManager() function.  For example, to obtain a pointer
    * to a hierarchy data operation object appropriate for a variable of type
    * double one makes the following call:
    * HierarchyDataOpsManager::getManager()->getOperationsDouble(
    *                                               var, hierarchy),
    * where ``var'' is a pointer to the variable, and hierarchy is a
    * pointer to the AMR hierarchy.
    *
    * Note that when the manager is accessed for the first time, the
    * Singleton instance is registered with the StartupShutdownManager
    * class which destroys such objects at program completion.  Thus,
    * users of this class do not explicitly allocate or deallocate the
    * Singleton instance.
    */
   static HierarchyDataOpsManager *
   getManager()
   {
      if (!s_pdat_op_manager_instance) {
         s_pdat_op_manager_instance = new HierarchyDataOpsManager();
      }
      return s_pdat_op_manager_instance;
   }

   //@{
   /*!
    * \name Getting operators.
    *
    * If a unique operator object is not requested, and if one already
    * exists for the hierarchy and variable specified, the existing one
    * will be created and returned.  Otherwise, a new one is created.
    * Objects created created for unique requests will not be used later
    * when an equivalent request is made.
    */
   /*!
    * \brief Return pointer to operation object for a double variable
    * on the given hierarchy.
    *
    * @pre variable
    * @pre hierarchy
    * @pre variable->getDim() == hierarchy.getDim()
    * @post returned value not NULL
    */
   virtual std::shared_ptr<HierarchyDataOpsReal<double> >
   getOperationsDouble(
      /*! operation should correspond to this variable */
      const std::shared_ptr<hier::Variable>& variable,
      /*! operation should correspond to this hierarchy */
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      /*! Whether a unique operator is requested */
      bool get_unique = false);

   /*!
    * \brief Return pointer to operation object for a float variable
    * on the given hierarchy.
    *
    * @pre variable
    * @pre hierarchy
    * @pre variable->getDim() == hierarchy.getDim()
    * @post returned value not NULL
    */
   virtual std::shared_ptr<HierarchyDataOpsReal<float> >
   getOperationsFloat(
      /*! operation should correspond to this variable */
      const std::shared_ptr<hier::Variable>& variable,
      /*! operation should correspond to this hierarchy */
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      /*! Whether a unique operator is requested */
      bool get_unique = false);

   /*!
    * \brief Return pointer to operation object for a complex variable
    * on the given hierarchy.
    *
    * @pre variable
    * @pre hierarchy
    * @pre variable->getDim() == hierarchy.getDim()
    * @post returned value not NULL
    */
   virtual std::shared_ptr<HierarchyDataOpsComplex>
   getOperationsComplex(
      /*! operation should correspond to this variable */
      const std::shared_ptr<hier::Variable>& variable,
      /*! operation should correspond to this hierarchy */
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      /*! Whether a unique operator is requested */
      bool get_unique = false);

   /*!
    * \brief Return pointer to operation object for a integer variable
    * on the given hierarchy.
    *
    * @pre variable
    * @pre hierarchy
    * @pre variable->getDim() == hierarchy.getDim()
    * @post returned value not NULL
    */
   virtual std::shared_ptr<HierarchyDataOpsInteger>
   getOperationsInteger(
      /*! operation should correspond to this variable */
      const std::shared_ptr<hier::Variable>& variable,
      /*! operation should correspond to this hierarchy */
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      /*! Whether a unique operator is requested */
      bool get_unique = false);
   //@}

protected:
   /**
    * The constructor for HierarchyDataOpsManager is protected.
    * Consistent with the definition of a Singleton class, only subclasses
    * have access to the constructor for the class.
    */
   HierarchyDataOpsManager();

   /**
    * The destructor for HierarchyDataOpsManager is protected.
    * See the comments for the constructor.
    */
   virtual ~HierarchyDataOpsManager();

   /**
    * Initialize Singleton instance with instance of subclass.  This function
    * is used to make the singleton object unique when inheriting from this
    * base class.
    *
    * @pre !s_pdat_op_manager_instance
    */
   void
   registerSingletonSubclassInstance(
      HierarchyDataOpsManager* subclass_instance)
   {
      if (!s_pdat_op_manager_instance) {
         s_pdat_op_manager_instance = subclass_instance;
      } else {
         TBOX_ERROR("HierarchyDataOpsManager internal error...\n"
            << "Attempting to set Singleton instance to subclass instance,"
            << "\n but Singleton instance already set." << std::endl);
      }
   }

private:
   /**
    * Deallocate the HierarchyDataOpsManager instance.  It is not
    * necessary to explicitly call freeManager() at program termination,
    * since it is automatically called by the StartupShutdownManager class.
    */
   static void
   shutdownCallback()
   {
      if (s_pdat_op_manager_instance) {
         delete s_pdat_op_manager_instance;
      }
      s_pdat_op_manager_instance = 0;
   }

   static HierarchyDataOpsManager* s_pdat_op_manager_instance;

   //@{ \name Operations for data of various types.
   std::vector<std::shared_ptr<HierarchyDataOpsReal<double> > >
   d_cell_ops_double;
   std::vector<std::shared_ptr<HierarchyDataOpsReal<double> > >
   d_face_ops_double;
   std::vector<std::shared_ptr<HierarchyDataOpsReal<double> > >
   d_node_ops_double;
   std::vector<std::shared_ptr<HierarchyDataOpsReal<double> > >
   d_side_ops_double;
   std::vector<std::shared_ptr<HierarchyDataOpsReal<double> > >
   d_edge_ops_double;

   std::vector<std::shared_ptr<HierarchyDataOpsReal<float> > >
   d_cell_ops_float;
   std::vector<std::shared_ptr<HierarchyDataOpsReal<float> > >
   d_face_ops_float;
   std::vector<std::shared_ptr<HierarchyDataOpsReal<float> > >
   d_side_ops_float;
   std::vector<std::shared_ptr<HierarchyDataOpsReal<float> > >
   d_node_ops_float;
   std::vector<std::shared_ptr<HierarchyDataOpsReal<float> > >
   d_edge_ops_float;

   std::vector<std::shared_ptr<HierarchyDataOpsComplex> >
   d_cell_ops_complex;
   std::vector<std::shared_ptr<HierarchyDataOpsComplex> >
   d_face_ops_complex;
   std::vector<std::shared_ptr<HierarchyDataOpsComplex> >
   d_side_ops_complex;
   std::vector<std::shared_ptr<HierarchyDataOpsComplex> >
   d_node_ops_complex;
   std::vector<std::shared_ptr<HierarchyDataOpsComplex> >
   d_edge_ops_complex;

   std::vector<std::shared_ptr<HierarchyDataOpsInteger> >
   d_cell_ops_int;
   std::vector<std::shared_ptr<HierarchyDataOpsInteger> >
   d_face_ops_int;
   std::vector<std::shared_ptr<HierarchyDataOpsInteger> >
   d_side_ops_int;
   std::vector<std::shared_ptr<HierarchyDataOpsInteger> >
   d_node_ops_int;
   std::vector<std::shared_ptr<HierarchyDataOpsInteger> >
   d_edge_ops_int;
   //@}

   static tbox::StartupShutdownManager::Handler
      s_shutdown_handler;

};

}
}

#endif
