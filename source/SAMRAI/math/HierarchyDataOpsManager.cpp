/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Singleton manager for hierarchy data operation objects.
 *
 ************************************************************************/
#include "SAMRAI/math/HierarchyDataOpsManager.h"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/HierarchyFaceDataOpsReal.h"
#include "SAMRAI/math/HierarchyNodeDataOpsReal.h"
#include "SAMRAI/math/HierarchySideDataOpsReal.h"
#include "SAMRAI/math/HierarchyEdgeDataOpsReal.h"

#include "SAMRAI/math/HierarchyCellDataOpsComplex.h"
#include "SAMRAI/math/HierarchyFaceDataOpsComplex.h"
#include "SAMRAI/math/HierarchyNodeDataOpsComplex.h"
#include "SAMRAI/math/HierarchySideDataOpsComplex.h"
#include "SAMRAI/math/HierarchyEdgeDataOpsComplex.h"

#include "SAMRAI/math/HierarchyCellDataOpsInteger.h"
#include "SAMRAI/math/HierarchyFaceDataOpsInteger.h"
#include "SAMRAI/math/HierarchyNodeDataOpsInteger.h"
#include "SAMRAI/math/HierarchySideDataOpsInteger.h"
#include "SAMRAI/math/HierarchyEdgeDataOpsInteger.h"

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

namespace SAMRAI {
namespace math {

/*
 *************************************************************************
 *
 * Static members for Singleton hierarchy operation manager class.
 *
 *************************************************************************
 */

HierarchyDataOpsManager *
HierarchyDataOpsManager::s_pdat_op_manager_instance = 0;

tbox::StartupShutdownManager::Handler
HierarchyDataOpsManager::s_shutdown_handler(
   0,
   0,
   HierarchyDataOpsManager::shutdownCallback,
   0,
   tbox::StartupShutdownManager::priorityHierarchyDataOpsManager);

/*
 *************************************************************************
 *
 * Empty constructor and destructor for hierarchy operation manager.
 *
 *************************************************************************
 */

HierarchyDataOpsManager::HierarchyDataOpsManager()
{
}

HierarchyDataOpsManager::~HierarchyDataOpsManager()
{
}

/*!
 * Return pointer to operation object for a double variable
 * on the given hierarchy.
 *
 * If a unique operator object is not requested, and if one already
 * exists for the hierarchy and variable specified, the existing one
 * will be created and returned.  Otherwise, a new one is created.
 * Objects created created for unique requests will not be used later
 * when an equivalent request is made.
 */

std::shared_ptr<HierarchyDataOpsReal<double> >
HierarchyDataOpsManager::getOperationsDouble(
   const std::shared_ptr<hier::Variable>& variable,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   bool get_unique)
{
   TBOX_ASSERT(variable);
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*variable, *hierarchy);

   const std::shared_ptr<pdat::CellVariable<double> > cellvar(
      std::dynamic_pointer_cast<pdat::CellVariable<double>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::FaceVariable<double> > facevar(
      std::dynamic_pointer_cast<pdat::FaceVariable<double>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::NodeVariable<double> > nodevar(
      std::dynamic_pointer_cast<pdat::NodeVariable<double>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::SideVariable<double> > sidevar(
      std::dynamic_pointer_cast<pdat::SideVariable<double>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::EdgeVariable<double> > edgevar(
      std::dynamic_pointer_cast<pdat::EdgeVariable<double>,
                                  hier::Variable>(variable));

   std::shared_ptr<HierarchyDataOpsReal<double> > ops;

   if (cellvar) {

      if (get_unique) {
         ops.reset(new HierarchyCellDataOpsReal<double>(hierarchy));
      } else {
         const int n = static_cast<int>(d_cell_ops_double.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy !=
                d_cell_ops_double[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_cell_ops_double[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyCellDataOpsReal<double>(hierarchy));
            d_cell_ops_double.resize(n + 1);
            d_cell_ops_double[n] = ops;
         }
      }

   } else if (facevar) {

      if (get_unique) {
         ops.reset(new HierarchyFaceDataOpsReal<double>(hierarchy));
      } else {
         const int n = static_cast<int>(d_face_ops_double.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy !=
                d_face_ops_double[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_face_ops_double[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyFaceDataOpsReal<double>(hierarchy));
            d_face_ops_double.resize(n + 1);
            d_face_ops_double[n] = ops;
         }
      }

   } else if (nodevar) {

      if (get_unique) {
         ops.reset(new HierarchyNodeDataOpsReal<double>(hierarchy));
      } else {
         const int n = static_cast<int>(d_node_ops_double.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy !=
                d_node_ops_double[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_node_ops_double[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyNodeDataOpsReal<double>(hierarchy));
            d_node_ops_double.resize(n + 1);
            d_node_ops_double[n] = ops;
         }
      }

   } else if (sidevar) {

      if (get_unique) {
         ops.reset(new HierarchySideDataOpsReal<double>(hierarchy));
      } else {
         const int n = static_cast<int>(d_side_ops_double.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy !=
                d_side_ops_double[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_side_ops_double[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchySideDataOpsReal<double>(hierarchy));
            d_side_ops_double.resize(n + 1);
            d_side_ops_double[n] = ops;
         }
      }

   } else if (edgevar) {

      if (get_unique) {
         ops.reset(new HierarchyEdgeDataOpsReal<double>(hierarchy));
      } else {
         const int n = static_cast<int>(d_edge_ops_double.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy !=
                d_edge_ops_double[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_edge_ops_double[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyEdgeDataOpsReal<double>(hierarchy));
            d_edge_ops_double.resize(n + 1);
            d_edge_ops_double[n] = ops;
         }
      }

   }

   if (!ops) {
      TBOX_ERROR("HierarchyDataOpsManager internal error...\n"
         << "Operations for variable " << variable->getName()
         << " not defined.");
   }

   return ops;
}

/*!
 * Return pointer to operation object for a float variable
 * on the given hierarchy.
 *
 * If a unique operator object is not requested, and if one already
 * exists for the hierarchy and variable specified, the existing one
 * will be created and returned.  Otherwise, a new one is created.
 * Objects created created for unique requests will not be used later
 * when an equivalent request is made.
 */

std::shared_ptr<HierarchyDataOpsReal<float> >
HierarchyDataOpsManager::getOperationsFloat(
   const std::shared_ptr<hier::Variable>& variable,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   bool get_unique)
{
   TBOX_ASSERT(variable);
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*variable, *hierarchy);

   const std::shared_ptr<pdat::CellVariable<float> > cellvar(
      std::dynamic_pointer_cast<pdat::CellVariable<float>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::FaceVariable<float> > facevar(
      std::dynamic_pointer_cast<pdat::FaceVariable<float>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::NodeVariable<float> > nodevar(
      std::dynamic_pointer_cast<pdat::NodeVariable<float>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::SideVariable<float> > sidevar(
      std::dynamic_pointer_cast<pdat::SideVariable<float>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::EdgeVariable<float> > edgevar(
      std::dynamic_pointer_cast<pdat::EdgeVariable<float>,
                                  hier::Variable>(variable));

   std::shared_ptr<HierarchyDataOpsReal<float> > ops;

   if (cellvar) {

      if (get_unique) {
         ops.reset(new HierarchyCellDataOpsReal<float>(hierarchy));
      } else {
         const int n = static_cast<int>(d_cell_ops_float.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy != d_cell_ops_float[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_cell_ops_float[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyCellDataOpsReal<float>(hierarchy));
            d_cell_ops_float.resize(n + 1);
            d_cell_ops_float[n] = ops;
         }
      }

   } else if (facevar) {

      if (get_unique) {
         ops.reset(new HierarchyFaceDataOpsReal<float>(hierarchy));
      } else {
         const int n = static_cast<int>(d_face_ops_float.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy != d_face_ops_float[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_face_ops_float[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyFaceDataOpsReal<float>(hierarchy));
            d_face_ops_float.resize(n + 1);
            d_face_ops_float[n] = ops;
         }
      }

   } else if (nodevar) {

      if (get_unique) {
         ops.reset(new HierarchyNodeDataOpsReal<float>(hierarchy));
      } else {
         const int n = static_cast<int>(d_node_ops_float.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy != d_node_ops_float[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_node_ops_float[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyNodeDataOpsReal<float>(hierarchy));
            d_node_ops_float.resize(n + 1);
            d_node_ops_float[n] = ops;
         }
      }

   } else if (sidevar) {

      if (get_unique) {
         ops.reset(new HierarchySideDataOpsReal<float>(hierarchy));
      } else {
         const int n = static_cast<int>(d_side_ops_float.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy != d_side_ops_float[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_side_ops_float[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchySideDataOpsReal<float>(hierarchy));
            d_side_ops_float.resize(n + 1);
            d_side_ops_float[n] = ops;
         }
      }

   } else if (edgevar) {

      if (get_unique) {
         ops.reset(new HierarchyEdgeDataOpsReal<float>(hierarchy));
      } else {
         const int n = static_cast<int>(d_edge_ops_float.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy != d_edge_ops_float[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_edge_ops_float[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyEdgeDataOpsReal<float>(hierarchy));
            d_edge_ops_float.resize(n + 1);
            d_edge_ops_float[n] = ops;
         }
      }

   }

   if (!ops) {
      TBOX_ERROR("HierarchyDataOpsManager internal error...\n"
         << "Operations for variable " << variable->getName()
         << " not defined.");
   }

   return ops;
}

/*!
 * Return pointer to operation object for a complex variable
 * on the given hierarchy.
 *
 * If a unique operator object is not requested, and if one already
 * exists for the hierarchy and variable specified, the existing one
 * will be created and returned.  Otherwise, a new one is created.
 * Objects created created for unique requests will not be used later
 * when an equivalent request is made.
 */

std::shared_ptr<HierarchyDataOpsComplex>
HierarchyDataOpsManager::getOperationsComplex(
   const std::shared_ptr<hier::Variable>& variable,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   bool get_unique)
{
   TBOX_ASSERT(variable);
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*variable, *hierarchy);

   const std::shared_ptr<pdat::CellVariable<dcomplex> > cellvar(
      std::dynamic_pointer_cast<pdat::CellVariable<dcomplex>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::FaceVariable<dcomplex> > facevar(
      std::dynamic_pointer_cast<pdat::FaceVariable<dcomplex>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::NodeVariable<dcomplex> > nodevar(
      std::dynamic_pointer_cast<pdat::NodeVariable<dcomplex>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::SideVariable<dcomplex> > sidevar(
      std::dynamic_pointer_cast<pdat::SideVariable<dcomplex>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::EdgeVariable<dcomplex> > edgevar(
      std::dynamic_pointer_cast<pdat::EdgeVariable<dcomplex>,
                                  hier::Variable>(variable));

   std::shared_ptr<HierarchyDataOpsComplex> ops;

   if (cellvar) {

      if (get_unique) {
         ops.reset(new HierarchyCellDataOpsComplex(hierarchy));
      } else {
         const int n = static_cast<int>(d_cell_ops_complex.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy !=
                d_cell_ops_complex[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_cell_ops_complex[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyCellDataOpsComplex(hierarchy));
            d_cell_ops_complex.resize(n + 1);
            d_cell_ops_complex[n] = ops;
         }
      }

   } else if (facevar) {

      if (get_unique) {
         ops.reset(new HierarchyFaceDataOpsComplex(hierarchy));
      } else {
         const int n = static_cast<int>(d_face_ops_complex.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy !=
                d_face_ops_complex[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_face_ops_complex[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyFaceDataOpsComplex(hierarchy));
            d_face_ops_complex.resize(n + 1);
            d_face_ops_complex[n] = ops;
         }
      }

   } else if (nodevar) {

      if (get_unique) {
         ops.reset(new HierarchyNodeDataOpsComplex(hierarchy));
      } else {
         const int n = static_cast<int>(d_node_ops_complex.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy !=
                d_node_ops_complex[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_node_ops_complex[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyNodeDataOpsComplex(hierarchy));
            d_node_ops_complex.resize(n + 1);
            d_node_ops_complex[n] = ops;
         }
      }

   } else if (sidevar) {

      if (get_unique) {
         ops.reset(new HierarchySideDataOpsComplex(hierarchy));
      } else {
         const int n = static_cast<int>(d_side_ops_complex.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy !=
                d_side_ops_complex[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_side_ops_complex[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchySideDataOpsComplex(hierarchy));
            d_side_ops_complex.resize(n + 1);
            d_side_ops_complex[n] = ops;
         }
      }

   } else if (edgevar) {

      if (get_unique) {
         ops.reset(new HierarchyEdgeDataOpsComplex(hierarchy));
      } else {
         const int n = static_cast<int>(d_edge_ops_complex.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy !=
                d_edge_ops_complex[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_edge_ops_complex[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyEdgeDataOpsComplex(hierarchy));
            d_edge_ops_complex.resize(n + 1);
            d_edge_ops_complex[n] = ops;
         }
      }

   }

   if (!ops) {
      TBOX_ERROR("HierarchyDataOpsManager internal error...\n"
         << "Operations for variable " << variable->getName()
         << " not defined.");
   }

   return ops;
}

/*!
 * Return pointer to operation object for an integer variable
 * on the given hierarchy.
 *
 * If a unique operator object is not requested, and if one already
 * exists for the hierarchy and variable specified, the existing one
 * will be created and returned.  Otherwise, a new one is created.
 * Objects created created for unique requests will not be used later
 * when an equivalent request is made.
 */

std::shared_ptr<HierarchyDataOpsInteger>
HierarchyDataOpsManager::getOperationsInteger(
   const std::shared_ptr<hier::Variable>& variable,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   bool get_unique)
{
   TBOX_ASSERT(variable);
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*variable, *hierarchy);

   const std::shared_ptr<pdat::CellVariable<int> > cellvar(
      std::dynamic_pointer_cast<pdat::CellVariable<int>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::FaceVariable<int> > facevar(
      std::dynamic_pointer_cast<pdat::FaceVariable<int>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::NodeVariable<int> > nodevar(
      std::dynamic_pointer_cast<pdat::NodeVariable<int>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::SideVariable<int> > sidevar(
      std::dynamic_pointer_cast<pdat::SideVariable<int>,
                                  hier::Variable>(variable));
   const std::shared_ptr<pdat::EdgeVariable<int> > edgevar(
      std::dynamic_pointer_cast<pdat::EdgeVariable<int>,
                                  hier::Variable>(variable));

   std::shared_ptr<HierarchyDataOpsInteger> ops;

   if (cellvar) {

      if (get_unique) {
         ops.reset(new HierarchyCellDataOpsInteger(hierarchy));
      } else {
         const int n = static_cast<int>(d_cell_ops_int.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy != d_cell_ops_int[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_cell_ops_int[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyCellDataOpsInteger(hierarchy));
            d_cell_ops_int.resize(n + 1);
            d_cell_ops_int[n] = ops;
         }
      }

   } else if (facevar) {

      if (get_unique) {
         ops.reset(new HierarchyFaceDataOpsInteger(hierarchy));
      } else {
         const int n = static_cast<int>(d_face_ops_int.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy != d_face_ops_int[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_face_ops_int[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyFaceDataOpsInteger(hierarchy));
            d_face_ops_int.resize(n + 1);
            d_face_ops_int[n] = ops;
         }
      }

   } else if (nodevar) {

      if (get_unique) {
         ops.reset(new HierarchyNodeDataOpsInteger(hierarchy));
      } else {
         const int n = static_cast<int>(d_node_ops_int.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy != d_node_ops_int[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_node_ops_int[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyNodeDataOpsInteger(hierarchy));
            d_node_ops_int.resize(n + 1);
            d_node_ops_int[n] = ops;
         }
      }

   } else if (sidevar) {

      if (get_unique) {
         ops.reset(new HierarchySideDataOpsInteger(hierarchy));
      } else {
         const int n = static_cast<int>(d_side_ops_int.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy != d_side_ops_int[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_side_ops_int[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchySideDataOpsInteger(hierarchy));
            d_side_ops_int.resize(n + 1);
            d_side_ops_int[n] = ops;
         }
      }

   } else if (edgevar) {

      if (get_unique) {
         ops.reset(new HierarchyEdgeDataOpsInteger(hierarchy));
      } else {
         const int n = static_cast<int>(d_edge_ops_int.size());
         for (int i = 0; i < n && !ops; ++i) {
            if (hierarchy != d_edge_ops_int[i]->getPatchHierarchy()) continue;
            // A compatible operator has been found at i.
            ops = d_edge_ops_int[i];
         }
         if (!ops) {
            // No compatible operator has been found.
            ops.reset(new HierarchyEdgeDataOpsInteger(hierarchy));
            d_edge_ops_int.resize(n + 1);
            d_edge_ops_int[n] = ops;
         }
      }

   }

   if (!ops) {
      TBOX_ERROR("HierarchyDataOpsManager internal error...\n"
         << "Operations for variable " << variable->getName()
         << " not defined.");
   }

   return ops;
}

}
}
