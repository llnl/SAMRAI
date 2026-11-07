/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   tbox
 *
 ************************************************************************/

#ifndef included_hier_ProcessorMapping
#define included_hier_ProcessorMapping

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/Utilities.h"

#include <vector>

namespace SAMRAI {
namespace hier {

/**
 * Class ProcessorMapping represents the processor assignments of
 * patches to processors.  It makes sure that all processor assignments
 * are in the range from 0 through NODES-1 and answers whether a particular
 * assignment is local to the processor.
 */

class ProcessorMapping
{
public:
   /**
    * Create a default processor mapping array with 0 elements.  Before
    * the mapping can be used, its size should be set using the function
    * setMappingSize() and each element of the mapping should be set
    * by setProcessorAssignment().
    */
   ProcessorMapping();

   /**
    * Create a processor mapping array with enough space for n elements.
    * All elements of the mapping are initialized to processor zero, but
    * they should be set by setProcessorAssignment() later.
    */
   explicit ProcessorMapping(
      const int n);

   /**
    * Create a new processor mapping and copy the processor assignments
    * from the argument.
    */
   ProcessorMapping(
      const ProcessorMapping& mapping);

   /**
    * Create a new processor mapping and get processor assignments
    * from the the std::vector<int> argument.
    */
   explicit ProcessorMapping(
      const std::vector<int>& mapping);

   /**
    * The destructor simply releases the storage for the mapping.
    */
   ~ProcessorMapping() = default;

   /**
    * Resize the mapping so that it has n elements.  Before it can be
    * used, each element should be set using setProcessorAssignment().
    */
   void
   setMappingSize(
      const size_t n);

   /**
    * Sets the number of boxes to n.
    * IMPORTANT NOTE: This method should only be used for
    * testing purposes.  Under normal circumstances, the number of
    * boxes is set by a call to tbox::SAMRAI_MPI::getNodes() and should NOT
    * be changed.
    */
   void
   setNumberNodes(
      const int n)
   {
      d_nodes = n;
   }

   /**
    * Return the processor assignment for the specified patch index.
    *
    * @pre (i >= 0) && (i < getProcessorMapping().size())
    */
   int
   getProcessorAssignment(
      const int i) const
   {
      TBOX_ASSERT((i >= 0) && (i < static_cast<int>(d_mapping.size())));
      return d_mapping[i];
   }

   /**
    * Set the processor assignment (second argument) for the specified
    * patch index (first argument).
    *
    * @pre (i >= 0) && (i < getProcessorMapping().size())
    * @pre (p >= 0) && (p < d_nodes)
    */
   void
   setProcessorAssignment(
      const int i,
      const int p)
   {
      TBOX_ASSERT((i >= 0) && (i < static_cast<int>(d_mapping.size())));
      TBOX_ASSERT((p >= 0) && (p < d_nodes));
      d_mapping[i] = p % d_nodes;
   }

   /**
    * Return an std::vector<int> of the processor mappings.
    */
   const std::vector<int>&
   getProcessorMapping() const
   {
      return d_mapping;
   }

   /**
    * Sets the processor mappings from an std::vector<int>.  Remaps the
    * processors so that patches are not accidentally mapped to
    * non-existent boxes.
    */
   void
   setProcessorMapping(
      const std::vector<int>& mapping);

   /**
    * Return the number of local indices (that is, those indices mapped to
    * the local processor).
    */
   int
   getNumberOfLocalIndices() const
   {
      computeLocalIndices();
      return d_local_id_count;
   }

   /**
    * Return a vector containing the local indices (that is,
    * those indices mapped to the local processor).
    */
   const std::vector<int>&
   getLocalIndices() const
   {
      computeLocalIndices();
      return d_local_indices;
   }

   /**
    * Return the total number of indices in the mapping array.
    */
   int
   getSizeOfMappingArray() const
   {
      return static_cast<int>(d_mapping.size());
   }

   /**
    * Check whether the specified index is a local index (that is, mapped
    * to the local processor).
    *
    * @pre (i >= 0) && (i < getProcessorMapping().size())
    */
   bool
   isMappingLocal(
      const int i) const
   {
      TBOX_ASSERT((i >= 0) && (i < static_cast<int>(d_mapping.size())));
      return d_mapping[i] == d_my_rank;
   }

private:
   /**
    * Fills in the array d_local_indices, and sets d_local_id_count.
    */
   void
   computeLocalIndices() const;

   ProcessorMapping&
   operator = (
      const ProcessorMapping&) = delete;

   int d_my_rank;
   int d_nodes;
   std::vector<int> d_mapping;
   mutable int d_local_id_count;
   mutable std::vector<int> d_local_indices;
};

}
}

#endif
