/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Vector class for data on SAMRAI hierarchy.
 *
 ************************************************************************/

#ifndef included_solv_SAMRAIVectorReal_C
#define included_solv_SAMRAIVectorReal_C

#include "SAMRAI/solv/SAMRAIVectorReal.h"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/HierarchyEdgeDataOpsReal.h"
#include "SAMRAI/math/HierarchyFaceDataOpsReal.h"
#include "SAMRAI/math/HierarchyNodeDataOpsReal.h"
#include "SAMRAI/math/HierarchySideDataOpsReal.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <typeinfo>
#include <cfloat>
#include <cmath>

namespace SAMRAI {
namespace solv {

/*
 *************************************************************************
 *
 * Initialize the static operators and counters.
 *
 *************************************************************************
 */

template<class TYPE>
int SAMRAIVectorReal<TYPE>::s_instance_counter[SAMRAI::MAX_DIM_VAL] = { 0 };

template<class TYPE>
std::shared_ptr<math::HierarchyDataOpsReal<TYPE> > SAMRAIVectorReal<TYPE>::
s_cell_ops[SAMRAI::MAX_DIM_VAL];
template<class TYPE>
std::shared_ptr<math::HierarchyDataOpsReal<TYPE> > SAMRAIVectorReal<TYPE>::
s_edge_ops[SAMRAI::MAX_DIM_VAL];
template<class TYPE>
std::shared_ptr<math::HierarchyDataOpsReal<TYPE> > SAMRAIVectorReal<TYPE>::
s_face_ops[SAMRAI::MAX_DIM_VAL];
template<class TYPE>
std::shared_ptr<math::HierarchyDataOpsReal<TYPE> > SAMRAIVectorReal<TYPE>::
s_node_ops[SAMRAI::MAX_DIM_VAL];
template<class TYPE>
std::shared_ptr<math::HierarchyDataOpsReal<TYPE> > SAMRAIVectorReal<TYPE>::
s_side_ops[SAMRAI::MAX_DIM_VAL];

/*
 *************************************************************************
 *
 * The constructor for SAMRAIVectorReal objects initializes
 * vector structure.
 *
 *************************************************************************
 */

template<class TYPE>
SAMRAIVectorReal<TYPE>::SAMRAIVectorReal(
   const std::string& name,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level):
   d_hierarchy(hierarchy),
   d_coarsest_level(coarsest_level),
   d_finest_level(finest_level),
   d_number_components(0)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((coarsest_level >= 0)
      && (finest_level >= coarsest_level)
      && (finest_level <= hierarchy->getFinestLevelNumber()));

   const tbox::Dimension& dim(d_hierarchy->getDim());

   ++SAMRAIVectorReal<TYPE>::s_instance_counter[dim.getValue() - 1];

   if (name.empty()) {
      d_vector_name = "SAMRAIVectorReal";
   } else {
      d_vector_name = name;
   }

   // Set default output stream
   d_output_stream = &tbox::plog;
}

/*
 *************************************************************************
 *
 * Destructor for SAMRAIVectorReal.
 * Component data storage is not deallocated here.
 *
 *************************************************************************
 */
template<class TYPE>
SAMRAIVectorReal<TYPE>::~SAMRAIVectorReal()
{

   const tbox::Dimension& dim(d_hierarchy->getDim());

   --SAMRAIVectorReal<TYPE>::s_instance_counter[dim.getValue() - 1];

   d_number_components = 0;

   d_component_variable.resize(0);
   d_component_data_id.resize(0);
   d_component_operations.resize(0);
   d_control_volume_data_id.resize(0);

   d_variableid_2_vectorcomponent_map.resize(0);

   if (SAMRAIVectorReal<TYPE>::s_instance_counter[dim.getValue() - 1] == 0) {
      if (SAMRAIVectorReal<TYPE>::s_cell_ops[dim.getValue() - 1]) {
         SAMRAIVectorReal<TYPE>::s_cell_ops[dim.getValue() - 1].reset();
      }
      if (SAMRAIVectorReal<TYPE>::s_edge_ops[dim.getValue() - 1]) {
         SAMRAIVectorReal<TYPE>::s_edge_ops[dim.getValue() - 1].reset();
      }
      if (SAMRAIVectorReal<TYPE>::s_face_ops[dim.getValue() - 1]) {
         SAMRAIVectorReal<TYPE>::s_face_ops[dim.getValue() - 1].reset();
      }
      if (SAMRAIVectorReal<TYPE>::s_node_ops[dim.getValue() - 1]) {
         SAMRAIVectorReal<TYPE>::s_node_ops[dim.getValue() - 1].reset();
      }
      if (SAMRAIVectorReal<TYPE>::s_side_ops[dim.getValue() - 1]) {
         SAMRAIVectorReal<TYPE>::s_side_ops[dim.getValue() - 1].reset();
      }
   }

}

/*
 *************************************************************************
 *
 * Set name string identifier for this vector object.
 *
 *************************************************************************
 */

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::setName(
   const std::string& name)
{
   d_vector_name = name;
}

/*
 *************************************************************************
 *
 * Reset vector levels and data operations.
 *
 *************************************************************************
 */

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::resetLevels(
   const int coarsest_level,
   const int finest_level)
{
   d_coarsest_level = coarsest_level;
   d_finest_level = finest_level;
}

/*
 *************************************************************************
 *
 * Create new vector with same structure as this and return new vector.
 *
 *************************************************************************
 */

template<class TYPE>
std::shared_ptr<SAMRAIVectorReal<TYPE> >
SAMRAIVectorReal<TYPE>::cloneVector(
   const std::string& name) const
{

   std::string new_name = (name.empty() ? d_vector_name : name);
   std::shared_ptr<SAMRAIVectorReal<TYPE> > new_vec(
      std::make_shared<SAMRAIVectorReal<TYPE> >(
         new_name,
         d_hierarchy,
         d_coarsest_level,
         d_finest_level));

   new_vec->setNumberOfComponents(d_number_components);

   hier::VariableDatabase* var_db = hier::VariableDatabase::getDatabase();

   for (int i = 0; i < d_number_components; ++i) {

      int new_id =
         var_db->registerClonedPatchDataIndex(d_component_variable[i],
            d_component_data_id[i]);

      new_vec->setComponent(i,
         d_component_variable[i],
         new_id,
         d_control_volume_data_id[i],
         d_component_operations[i]);
   }

   return new_vec;
}

/*
 *************************************************************************
 *
 * Deallocate vector data storage and remove data indices from database.
 *
 *************************************************************************
 */

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::freeVectorComponents()
{
   // deallocate storage for vector components
   deallocateVectorData();

   hier::VariableDatabase* var_db = hier::VariableDatabase::getDatabase();

   // free entries from variable database and return
   // patch descriptor indices
   for (int i = 0; i < d_number_components; ++i) {
      var_db->removePatchDataIndex(d_component_data_id[i]);
   }

   // reset variable state
   d_number_components = 0;

   d_component_variable.resize(0);
   d_component_data_id.resize(0);
   d_component_operations.resize(0);
   d_control_volume_data_id.resize(0);

   d_variableid_2_vectorcomponent_map.resize(0);
}

/*
 *************************************************************************
 *
 * Add new component to vector structure given a variable and the
 * patch descriptor indices for its data and an appropriate
 * control volume.
 *
 *************************************************************************
 */

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::addComponent(
   const std::shared_ptr<hier::Variable>& var,
   const int comp_data_id,
   const int comp_vol_id,
   const std::shared_ptr<math::HierarchyDataOpsReal<TYPE> >& vop)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   hier::VariableDatabase* var_db =
      hier::VariableDatabase::getDatabase();
   std::shared_ptr<hier::PatchDescriptor> patch_descriptor(
      var_db->getPatchDescriptor());
   if (!var_db->checkVariablePatchDataIndexType(var, comp_data_id)) {
      auto varptr = var.get();
      auto factory = patch_descriptor->getPatchDataFactory(comp_data_id).get();
      TBOX_ERROR("Error in SAMRAIVectorReal::addComponent : "
         << "Vector name = " << d_vector_name
         << "\nVariable " << var->getName()
         << " type does not match data type associated with"
         << " comp_data_id patch data index function argument"
         << "\n\t var type = " << typeid(*varptr).name()
         << "\n\t comp_data_id type = "
         << typeid(*factory).name()
         << std::endl);
   }

   if (comp_vol_id >= 0) {
      if (!var_db->checkVariablePatchDataIndexType(var, comp_vol_id)) {
         auto varptr = var.get();
         auto factory = patch_descriptor->getPatchDataFactory(comp_vol_id).get();
         TBOX_ERROR("Error in SAMRAIVectorReal::addComponent : "
            << "Vector name = " << d_vector_name
            << "\nVariable " << var->getName()
            << " type does not match data type associated with"
            << " comp_vol_id patch data index function argument"
            << "\n\t var type = " << typeid(*varptr).name()
            << "\n\t comp_vol_id type = "
            << typeid(*factory).name()
            << std::endl);
      }
   }
#endif

   ++d_number_components;

   d_component_variable.resize(d_number_components);
   d_component_data_id.resize(d_number_components);
   d_component_operations.resize(d_number_components);
   d_control_volume_data_id.resize(d_number_components);

   hier::VariableDatabase::getDatabase()->registerPatchDataIndex(var,
      comp_data_id);

   setComponent(d_number_components - 1,
      var,
      comp_data_id,
      comp_vol_id,
      vop);
}

/*
 *************************************************************************
 *
 * Routines to allocate and deallocate data for all vector components.
 *
 *************************************************************************
 */
template<class TYPE>
void
SAMRAIVectorReal<TYPE>::allocateVectorData(
   const double timestamp)
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (int i = 0; i < d_number_components; ++i) {
         level->allocatePatchData(d_component_data_id[i], timestamp);
      }
   }
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::deallocateVectorData()
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (int i = 0; i < d_number_components; ++i) {
         level->deallocatePatchData(d_component_data_id[i]);
      }
   }
}

/*
 *************************************************************************
 *
 * Print Vector attributes and data.
 *
 *************************************************************************
 */

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::print(
   std::ostream& s,
   bool interior_only) const
{
   s << "\nVector : " << getName() << std::endl;
   s << "coarsest level = " << d_coarsest_level
     << " : finest level = " << d_finest_level << std::endl;
   s << "d_number_components = " << d_number_components << std::endl;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      s << "Printing data components on level " << ln << std::endl;
      for (int i = 0; i < d_number_components; ++i) {
         s << "Vector component index = " << i << std::endl;
         d_component_operations[i]->resetLevels(ln, ln);
         d_component_operations[i]->printData(d_component_data_id[i],
            s,
            interior_only);
      }
   }
}

/*
 *************************************************************************
 *
 * Private member functions to set the number of vector components
 * to set individual components.   These routines are used when cloning
 * vectors and/or adding components.
 *
 *************************************************************************
 */

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::setNumberOfComponents(
   int num_comp)
{
   d_number_components = num_comp;

   d_component_variable.resize(d_number_components);
   d_component_data_id.resize(d_number_components);
   d_component_operations.resize(d_number_components);
   d_control_volume_data_id.resize(d_number_components);
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::setComponent(
   const int comp_id,
   const std::shared_ptr<hier::Variable>& var,
   const int data_id,
   const int vol_id,
   const std::shared_ptr<math::HierarchyDataOpsReal<TYPE> >& vop)
{
   TBOX_ASSERT(comp_id < d_number_components);

   const tbox::Dimension& dim(d_hierarchy->getDim());

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dim, *var);

   d_component_variable[comp_id] = var;
   d_component_data_id[comp_id] = data_id;
   if (!vop) {

      const std::shared_ptr<pdat::CellVariable<TYPE> > cellvar(
         std::dynamic_pointer_cast<pdat::CellVariable<TYPE>, hier::Variable>(
            var));
      const std::shared_ptr<pdat::EdgeVariable<TYPE> > edgevar(
         std::dynamic_pointer_cast<pdat::EdgeVariable<TYPE>, hier::Variable>(
            var));
      const std::shared_ptr<pdat::FaceVariable<TYPE> > facevar(
         std::dynamic_pointer_cast<pdat::FaceVariable<TYPE>, hier::Variable>(
            var));
      const std::shared_ptr<pdat::NodeVariable<TYPE> > nodevar(
         std::dynamic_pointer_cast<pdat::NodeVariable<TYPE>, hier::Variable>(
            var));
      const std::shared_ptr<pdat::SideVariable<TYPE> > sidevar(
         std::dynamic_pointer_cast<pdat::SideVariable<TYPE>, hier::Variable>(
            var));

      if (cellvar) {
         if (!SAMRAIVectorReal<TYPE>::s_cell_ops[dim.getValue() - 1]) {
            SAMRAIVectorReal<TYPE>::s_cell_ops[dim.getValue() - 1].reset(
               new math::HierarchyCellDataOpsReal<TYPE>(d_hierarchy,
                  d_coarsest_level,
                  d_finest_level));
         }
         d_component_operations[comp_id] =
            SAMRAIVectorReal<TYPE>::s_cell_ops[dim.getValue() - 1];
      } else if (edgevar) {
         if (!SAMRAIVectorReal<TYPE>::s_edge_ops[dim.getValue() - 1]) {
            SAMRAIVectorReal<TYPE>::s_edge_ops[dim.getValue() - 1].reset(
               new math::HierarchyEdgeDataOpsReal<TYPE>(d_hierarchy,
                  d_coarsest_level,
                  d_finest_level));
         }
         d_component_operations[comp_id] =
            SAMRAIVectorReal<TYPE>::s_edge_ops[dim.getValue() - 1];
      } else if (facevar) {
         if (!SAMRAIVectorReal<TYPE>::s_face_ops[dim.getValue() - 1]) {
            SAMRAIVectorReal<TYPE>::s_face_ops[dim.getValue() - 1].reset(
               new math::HierarchyFaceDataOpsReal<TYPE>(d_hierarchy,
                  d_coarsest_level,
                  d_finest_level));
         }
         d_component_operations[comp_id] =
            SAMRAIVectorReal<TYPE>::s_face_ops[dim.getValue() - 1];
      } else if (nodevar) {
         if (!SAMRAIVectorReal<TYPE>::s_node_ops[dim.getValue() - 1]) {
            SAMRAIVectorReal<TYPE>::s_node_ops[dim.getValue() - 1].reset(
               new math::HierarchyNodeDataOpsReal<TYPE>(d_hierarchy,
                  d_coarsest_level,
                  d_finest_level));
         }
         d_component_operations[comp_id] =
            SAMRAIVectorReal<TYPE>::s_node_ops[dim.getValue() - 1];
      } else if (sidevar) {
         if (!SAMRAIVectorReal<TYPE>::s_side_ops[dim.getValue() - 1]) {
            SAMRAIVectorReal<TYPE>::s_side_ops[dim.getValue() - 1].reset(
               new math::HierarchySideDataOpsReal<TYPE>(d_hierarchy,
                  d_coarsest_level,
                  d_finest_level));
         }
         d_component_operations[comp_id] =
            SAMRAIVectorReal<TYPE>::s_side_ops[dim.getValue() - 1];
      }
   } else {
      d_component_operations[comp_id] = vop;
   }

   TBOX_ASSERT(d_component_operations[comp_id]);

   d_control_volume_data_id[comp_id] = vol_id;

   int var_id = var->getInstanceIdentifier();

   int oldsize = static_cast<int>(d_variableid_2_vectorcomponent_map.size());
   int newsize = var_id + 1;
   if (oldsize < newsize) {
      newsize = tbox::MathUtilities<int>::Max(
            oldsize + DESCRIPTOR_ID_ARRAY_SCRATCH_SPACE, newsize);
      d_variableid_2_vectorcomponent_map.resize(newsize);
      for (int i = oldsize; i < newsize; ++i) {
         d_variableid_2_vectorcomponent_map[i] = -1;
      }
   }

   TBOX_ASSERT(d_variableid_2_vectorcomponent_map[var_id] == -1);

   d_variableid_2_vectorcomponent_map[var_id] = comp_id;
}

/*
 *************************************************************************
 *
 * The remaining functions are basic vector kernel routines.
 * The operation for each component is performed by its hierarchy
 * data operation object.
 *
 *************************************************************************
 */

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::copyVector(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& src_vec,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->copyData(d_component_data_id[i],
         src_vec->getComponentDescriptorIndex(i),
         interior_only);
   }
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::swapVectors(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& other)
{
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->swapData(d_component_data_id[i],
         other->getComponentDescriptorIndex(i));
   }
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::setToScalar(
   const TYPE& alpha,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->setToScalar(d_component_data_id[i],
         alpha,
         interior_only);
   }
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::scale(
   const TYPE& alpha,
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->scale(d_component_data_id[i],
         alpha,
         x->getComponentDescriptorIndex(i),
         interior_only);
   }
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::addScalar(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
   const TYPE& alpha,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->addScalar(d_component_data_id[i],
         x->getComponentDescriptorIndex(i),
         alpha,
         interior_only);
   }
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::add(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& y,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->add(d_component_data_id[i],
         x->getComponentDescriptorIndex(i),
         y->getComponentDescriptorIndex(i),
         interior_only);
   }
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::subtract(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& y,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->subtract(d_component_data_id[i],
         x->getComponentDescriptorIndex(i),
         y->getComponentDescriptorIndex(i),
         interior_only);
   }
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::multiply(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& y,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->multiply(d_component_data_id[i],
         x->getComponentDescriptorIndex(i),
         y->getComponentDescriptorIndex(i),
         interior_only);
   }
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::divide(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& y,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->divide(d_component_data_id[i],
         x->getComponentDescriptorIndex(i),
         y->getComponentDescriptorIndex(i),
         interior_only);
   }
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::reciprocal(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->reciprocal(d_component_data_id[i],
         x->getComponentDescriptorIndex(i),
         interior_only);
   }
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::linearSum(
   const TYPE& alpha,
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
   const TYPE& beta,
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& y,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->linearSum(d_component_data_id[i],
         alpha,
         x->getComponentDescriptorIndex(i),
         beta,
         y->getComponentDescriptorIndex(i),
         interior_only);
   }
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::axpy(
   const TYPE& alpha,
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& y,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->axpy(d_component_data_id[i],
         alpha,
         x->getComponentDescriptorIndex(i),
         y->getComponentDescriptorIndex(i),
         interior_only);
   }
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::abs(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->abs(d_component_data_id[i],
         x->getComponentDescriptorIndex(i),
         interior_only);
   }
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::setRandomValues(
   const TYPE& width,
   const TYPE& low,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->setRandomValues(d_component_data_id[i],
         width,
         low,
         interior_only);
   }
}

template<class TYPE>
double
SAMRAIVectorReal<TYPE>::L1Norm(
   bool local_only) const
{
   double norm = 0.0;

   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      norm += d_component_operations[i]->L1Norm(d_component_data_id[i],
            d_control_volume_data_id[i],
            local_only);
   }

   return norm;
}

template<class TYPE>
double
SAMRAIVectorReal<TYPE>::L2Norm(
   bool local_only) const
{
   double norm_squared = 0.0;

   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      double comp_norm =
         d_component_operations[i]->L2Norm(d_component_data_id[i],
            d_control_volume_data_id[i],
            local_only);
      norm_squared += comp_norm * comp_norm;
   }

   return sqrt(norm_squared);
}

template<class TYPE>
double
SAMRAIVectorReal<TYPE>::weightedL2Norm(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& wgt) const
{
   double norm_squared = 0.0;

   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      double comp_norm = d_component_operations[i]->weightedL2Norm(
            d_component_data_id[i],
            wgt->getComponentDescriptorIndex(i),
            d_control_volume_data_id[i]);
      norm_squared += comp_norm * comp_norm;
   }

   return sqrt(norm_squared);
}

template<class TYPE>
double
SAMRAIVectorReal<TYPE>::RMSNorm() const
{
   double num = L2Norm();

   double denom = 0.0;
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      if (d_control_volume_data_id[i] < 0) {
         denom += double(d_component_operations[i]->
                         numberOfEntries(d_component_data_id[i], true));
      } else {
         denom += d_component_operations[i]->
            sumControlVolumes(d_component_data_id[i],
               d_control_volume_data_id[i]);
      }
   }

   double norm = 0.0;
   if (denom > 0.0) norm = num / sqrt(denom);
   return norm;
}

template<class TYPE>
double
SAMRAIVectorReal<TYPE>::weightedRMSNorm(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& wgt) const
{
   double num = weightedL2Norm(wgt);

   double denom = 0.0;
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      if (d_control_volume_data_id[i] < 0) {
         denom += double(d_component_operations[i]->
                         numberOfEntries(d_component_data_id[i], true));
      } else {
         denom += d_component_operations[i]->
            sumControlVolumes(d_component_data_id[i],
               d_control_volume_data_id[i]);
      }
   }

   double norm = 0.0;
   if (denom > 0.0) norm = num / sqrt(denom);
   return norm;
}

template<class TYPE>
double
SAMRAIVectorReal<TYPE>::maxNorm(
   bool local_only) const
{
   double norm = 0.0;

   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      norm = tbox::MathUtilities<double>::Max(norm,
            d_component_operations[i]->maxNorm(
               d_component_data_id[i],
               d_control_volume_data_id[i],
               local_only));
   }

   return norm;
}

template<class TYPE>
TYPE
SAMRAIVectorReal<TYPE>::dot(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
   bool local_only) const
{
   TYPE dprod = 0.0;

   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      dprod += d_component_operations[i]->dot(d_component_data_id[i],
            x->getComponentDescriptorIndex(i),
            d_control_volume_data_id[i],
            local_only);
   }

   return dprod;
}

template<class TYPE>
int
SAMRAIVectorReal<TYPE>::computeConstrProdPos(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x) const
{
   int test = 1;

   int i = 0;
   while (test && (i < d_number_components)) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      test = tbox::MathUtilities<int>::Min(test,
            d_component_operations[i]->
            computeConstrProdPos(d_component_data_id[i],
               x->getComponentDescriptorIndex(i),
               d_control_volume_data_id[i]));
      ++i;
   }

   return test;
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::compareToScalar(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
   const TYPE& alpha)
{
   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->
      compareToScalar(d_component_data_id[i],
         x->getComponentDescriptorIndex(i),
         alpha,
         d_control_volume_data_id[i]);
   }
}

template<class TYPE>
int
SAMRAIVectorReal<TYPE>::testReciprocal(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x)
{
   int test = 1;

   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      test = tbox::MathUtilities<int>::Min(test,
            d_component_operations[i]->
            testReciprocal(d_component_data_id[i],
               x->getComponentDescriptorIndex(i),
               d_control_volume_data_id[i]));
   }

   return test;
}

template<class TYPE>
TYPE
SAMRAIVectorReal<TYPE>::maxPointwiseDivide(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& denom) const
{
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());
   TYPE max = 0.0;

   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      TYPE component_max =
         d_component_operations[i]->maxPointwiseDivide(d_component_data_id[i],
            denom->getComponentDescriptorIndex(i),
            true);
      max = tbox::MathUtilities<TYPE>::Max(max, component_max);
   }

   if (mpi.getSize() > 1) {
      mpi.AllReduce(&max, 1, MPI_MAX);
   }
   return max;
}

template<class TYPE>
TYPE
SAMRAIVectorReal<TYPE>::min(
   const bool interior_only) const
{
   TYPE minval = tbox::MathUtilities<TYPE>::getMax();

   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      minval = tbox::MathUtilities<TYPE>::Min(
            minval,
            d_component_operations[i]->min(d_component_data_id[i],
               interior_only));
   }

   return minval;
}

template<class TYPE>
TYPE
SAMRAIVectorReal<TYPE>::max(
   const bool interior_only) const
{
   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();

   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      maxval = tbox::MathUtilities<TYPE>::Max(
            maxval,
            d_component_operations[i]->max(d_component_data_id[i],
               interior_only));
   }

   return maxval;
}

template<class TYPE>
int64_t
SAMRAIVectorReal<TYPE>::getLength(
   const bool interior_only) const
{
   int64_t length = 0;

   for (int i = 0; i < d_number_components; ++i) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      length += d_component_operations[i]->getLength(d_component_data_id[i],
                                                     interior_only);
   }

   return length;
}

template<class TYPE>
void
SAMRAIVectorReal<TYPE>::setOutputStream(
   std::ostream& s)
{
   d_output_stream = &s;
}

template<class TYPE>
std::ostream&
SAMRAIVectorReal<TYPE>::getOutputStream()
{
   return *d_output_stream;
}

template<class TYPE>
const std::string&
SAMRAIVectorReal<TYPE>::getName() const
{
   return d_vector_name;
}

template<class TYPE>
std::shared_ptr<hier::PatchHierarchy>
SAMRAIVectorReal<TYPE>::getPatchHierarchy() const
{
   return d_hierarchy;
}

template<class TYPE>
int
SAMRAIVectorReal<TYPE>::getCoarsestLevelNumber() const
{
   return d_coarsest_level;
}

template<class TYPE>
int
SAMRAIVectorReal<TYPE>::getFinestLevelNumber() const
{
   return d_finest_level;
}

template<class TYPE>
int
SAMRAIVectorReal<TYPE>::getNumberOfComponents() const
{
   return d_number_components;
}

template<class TYPE>
std::shared_ptr<hier::PatchData>
SAMRAIVectorReal<TYPE>::getComponentPatchData(
   const int comp_id,
   const hier::Patch& patch) const
{
   TBOX_ASSERT(comp_id >= 0 && comp_id < d_number_components);
   return patch.getPatchData(d_component_data_id[comp_id]);
}

template<class TYPE>
std::shared_ptr<hier::PatchData>
SAMRAIVectorReal<TYPE>::getComponentPatchData(
   const std::shared_ptr<hier::Variable>& var,
   const hier::Patch& patch) const
{
   TBOX_ASSERT(var);
   TBOX_ASSERT(d_variableid_2_vectorcomponent_map[
         var->getInstanceIdentifier()] >= 0);
   return patch.getPatchData(
      d_component_data_id[
         d_variableid_2_vectorcomponent_map[
            var->getInstanceIdentifier()]]);
}

template<class TYPE>
std::shared_ptr<hier::Variable>
SAMRAIVectorReal<TYPE>::getComponentVariable(
   const int comp_id) const
{
   TBOX_ASSERT(comp_id >= 0 && comp_id < d_number_components);
   return d_component_variable[comp_id];

}

template<class TYPE>
int
SAMRAIVectorReal<TYPE>::getComponentDescriptorIndex(
   const int comp_id) const
{
   TBOX_ASSERT(comp_id >= 0 && comp_id < d_number_components);
   return d_component_data_id[comp_id];
}

template<class TYPE>
int
SAMRAIVectorReal<TYPE>::getControlVolumeIndex(
   const int comp_id) const
{
   TBOX_ASSERT(comp_id >= 0 && comp_id < d_number_components);
   return d_control_volume_data_id[comp_id];
}

}
}
#endif
