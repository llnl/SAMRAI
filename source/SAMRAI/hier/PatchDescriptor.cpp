/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Factory class for patch data objects that live on a patch
 *
 ************************************************************************/
#include "SAMRAI/hier/PatchDescriptor.h"

#include "SAMRAI/tbox/SAMRAIManager.h"

#include <typeinfo>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

const int PatchDescriptor::INDEX_UNDEFINED = -1;

/*
 *************************************************************************
 *
 * The constructor sets the max number of registered components to zero
 * and allocates the factory and name arrays to the fixed length set
 * by the SAMRAIManager utility.  The free list of indices
 * is initialized to the full set of potentially used indices.
 *
 * The destructor clears the free index list and implicitly
 * deallocates the arrays of name strings and factory pointers.
 *
 *************************************************************************
 */

PatchDescriptor::PatchDescriptor():
   d_min_gcw()
{
   const int max_num_patch_data_components_allowed =
      tbox::SAMRAIManager::getMaxNumberPatchDataEntries();
   d_max_number_registered_components = 0;
   d_names.resize(max_num_patch_data_components_allowed);
   d_factories.resize(max_num_patch_data_components_allowed);
   for (int i = 0; i < max_num_patch_data_components_allowed; ++i) {
      d_free_indices.push_back(i);
   }
   for (unsigned short d = 0; d < SAMRAI::MAX_DIM_VAL; ++d) {
      d_min_gcw.push_back(
         IntVector::getZero(tbox::Dimension(static_cast<unsigned short>(d + 1))));
   }
}

PatchDescriptor::~PatchDescriptor()
{
   d_free_indices.clear();
}

/*
 *************************************************************************
 *
 * Add the new factory to the list of patch data factories and assign
 * it an integer index identifier.  Use a free list item if possible.
 *
 *************************************************************************
 */

int
PatchDescriptor::definePatchDataComponent(
   const std::string& name,
   const std::shared_ptr<PatchDataFactory>& factory)
{
   TBOX_ASSERT(!name.empty());
   TBOX_ASSERT(factory);

   if (d_free_indices.empty()) {
      int old_size = static_cast<int>(d_names.size());
      int new_size = old_size +
         tbox::SAMRAIManager::getMaxNumberPatchDataEntries();
      TBOX_ASSERT(new_size > old_size);
      d_names.resize(new_size);
      d_factories.resize(new_size);
      for (int i = old_size; i < new_size; ++i) {
         d_free_indices.push_back(i);
      }
      tbox::SAMRAIManager::setMaxNumberPatchDataEntries(new_size);
   }

   int ret_index = d_free_indices.front();
   d_free_indices.pop_front();
   if (d_max_number_registered_components < ret_index + 1) {
      d_max_number_registered_components = ret_index + 1;
   }
   d_factories[ret_index] = factory;
   d_names[ret_index] = name;

   return ret_index;
}

/*
 *************************************************************************
 *
 * Remove the specified patch data factory index and place the index on
 * the list of free indices.
 *
 *************************************************************************
 */

void
PatchDescriptor::removePatchDataComponent(
   const int id)
{
   if ((id >= 0) && (id < d_max_number_registered_components)) {
      if (!d_names[id].empty()) {
         d_names[id] = std::string();
      }
      if (d_factories[id]) {
         d_factories[id].reset();
         d_free_indices.push_front(id);
      }
   }
}

/*
 *************************************************************************
 *
 * Look up the factory by name; if no matching factory exists, then a
 * pointer to null is returned.  The first matching factory is returned.
 *
 *************************************************************************
 */

std::shared_ptr<PatchDataFactory>
PatchDescriptor::getPatchDataFactory(
   const std::string& name) const
{
   std::shared_ptr<PatchDataFactory> factory;
   const int id = mapNameToIndex(name);
   if (id >= 0) {
      factory = d_factories[id];
   }
   return factory;
}

/*
 *************************************************************************
 *
 * Search the factory list for a match and return the associated
 * factory.  If no match exists, return a negative identifier.
 *
 *************************************************************************
 */

int
PatchDescriptor::mapNameToIndex(
   const std::string& name) const
{
   int ret_index = INDEX_UNDEFINED;
   int id = 0;
   while ((ret_index == INDEX_UNDEFINED) &&
          (id < d_max_number_registered_components)) {
      if (name == d_names[id]) {
         ret_index = id;
      }
      ++id;
   }
   return ret_index;
}

/*
 *************************************************************************
 *
 * Print index, name, and factory data for the patch descriptor.
 *
 *************************************************************************
 */

void
PatchDescriptor::printClassData(
   std::ostream& stream) const
{
   stream << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
          << std::endl;
   stream << "Printing PatchDescriptor state ..." << std::endl;
   stream << "this = " << (PatchDescriptor *)this << std::endl;
   stream << "d_max_number_registered_components = "
          << d_max_number_registered_components << std::endl;
   stream << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
          << std::endl;
   for (int i = 0; i < d_max_number_registered_components; ++i) {
      stream << "Patch Data Index=" << i << std::endl;
      if (d_factories[i]) {
         auto& f = *d_factories[i];
         stream << "   Patch Data Factory Name = "
                << d_names[i] << std::endl;
         stream << "   Patch Data Factory = "
                << typeid(f).name() << std::endl;
      } else {
         stream << "   Patch Data Factory = NULL" << std::endl;
      }
   }
   stream << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
          << std::endl;
}

/*
 *************************************************************************
 * Return the maximum ghost cell width across all factories and the
 * user-specified minimum value.
 *************************************************************************
 */

IntVector
PatchDescriptor::getMaxGhostWidth(
   const tbox::Dimension& dim) const
{
   IntVector max_gcw(d_min_gcw[dim.getValue() - 1]);
   for (int i = 0; i < d_max_number_registered_components; ++i) {
      if (d_factories[i] && (d_factories[i]->getDim() == dim)) {
         max_gcw.max(d_factories[i]->getGhostCellWidth());
      }
   }
   return max_gcw;
}

/*
 *************************************************************************
 * Return the dimension of the data for the given data_id.
 *************************************************************************
 */

tbox::Dimension
PatchDescriptor::getPatchDataDim(
   int patch_id) const
{
   return d_factories[patch_id]->getDim();
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
