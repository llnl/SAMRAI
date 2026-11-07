/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   tbox
 *
 ************************************************************************/

#include "SAMRAI/hier/ProcessorMapping.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

namespace SAMRAI {
namespace hier {

ProcessorMapping::ProcessorMapping():
   d_my_rank(-1),
   d_nodes(-1),
   d_mapping(0),
   d_local_id_count(-1)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   d_my_rank = mpi.getRank();
   d_nodes = mpi.getSize();
}

ProcessorMapping::ProcessorMapping(
   const int n):
   d_my_rank(-1),
   d_nodes(-1),
   d_mapping(n),
   d_local_id_count(-1)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   d_my_rank = mpi.getRank();
   d_nodes = mpi.getSize();
   for (int i = 0; i < n; ++i) {
      d_mapping[i] = 0;
   }
}

ProcessorMapping::ProcessorMapping(
   const ProcessorMapping& mapping):
   d_my_rank(-1),
   d_nodes(-1),
   d_mapping(mapping.d_mapping.size()),
   d_local_id_count(-1)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   d_my_rank = mpi.getRank();
   d_nodes = mpi.getSize();
   const int n = static_cast<int>(d_mapping.size());
   for (int i = 0; i < n; ++i) {
      d_mapping[i] = mapping.d_mapping[i];
   }
}

ProcessorMapping::ProcessorMapping(
   const std::vector<int>& mapping):
   d_my_rank(-1),
   d_nodes(-1),
   d_local_id_count(-1)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   d_my_rank = mpi.getRank();
   d_nodes = mpi.getSize();
   setProcessorMapping(mapping);
}

void
ProcessorMapping::setMappingSize(
   const size_t n)
{
   d_mapping.resize(n);

   for (size_t i = 0; i < n; ++i) {
      d_mapping[i] = 0;
   }
   d_local_id_count = -1;
}

void
ProcessorMapping::setProcessorMapping(
   const std::vector<int>& mapping)
{
   d_mapping.resize(mapping.size());

   for (int i = 0; i < static_cast<int>(d_mapping.size()); ++i) {
      //  (mapping[i] % d_nodes) keeps patches from being assigned
      //  non-existent processors.
      setProcessorAssignment(i, mapping[i] % d_nodes);
   }
   d_local_id_count = -1;
}

void
ProcessorMapping::computeLocalIndices() const
{
   if (d_local_id_count != -1) {
      return;
   }

   /*
    * first, count the number of local indices,
    * so we can set the array size.
    */
   const int n = static_cast<int>(d_mapping.size());
   d_local_id_count = 0;

   for (int i = 0; i < n; ++i) {
      if (d_mapping[i] == d_my_rank) {
         ++d_local_id_count;
      }
   }

   /*
    * second, resize the array and fill in the data
    */
   d_local_indices.resize(d_local_id_count);
   int idx = 0;
   for (int i = 0; i < n; ++i) {
      if (d_mapping[i] == d_my_rank) {
         d_local_indices[idx++] = i;
      }
   }
}

}
}
