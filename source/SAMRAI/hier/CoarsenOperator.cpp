/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract base class for spatial coarsening operators.
 *
 ************************************************************************/
#include "SAMRAI/hier/CoarsenOperator.h"

#include "SAMRAI/tbox/StartupShutdownManager.h"

#include "SAMRAI/tbox/OpenMPUtilities.h"

namespace SAMRAI {
namespace hier {

std::multimap<std::string, CoarsenOperator *> CoarsenOperator::s_lookup_table;
TBOX_omp_lock_t CoarsenOperator::l_lookup_table;

tbox::StartupShutdownManager::Handler
CoarsenOperator::s_finalize_handler(
   CoarsenOperator::initializeCallback,
   0,
   0,
   CoarsenOperator::finalizeCallback,
   tbox::StartupShutdownManager::priorityList);

CoarsenOperator::CoarsenOperator(
   const std::string& name):
   d_name(name)
{
   registerInLookupTable(name);
}

CoarsenOperator::~CoarsenOperator()
{
   removeFromLookupTable(d_name);
}

void
CoarsenOperator::registerInLookupTable(
   const std::string& name)
{
   TBOX_omp_set_lock(&l_lookup_table);
   s_lookup_table.insert(
      std::pair<std::string, CoarsenOperator *>(name, this));
   TBOX_omp_unset_lock(&l_lookup_table);
}

void
CoarsenOperator::removeFromLookupTable(
   const std::string& name)
{
   /*
    * The lookup table might be empty if static CoarsenOperator's are used
    * in which case the table will have been cleared before the statics
    * are destroyed.
    */
   TBOX_omp_set_lock(&l_lookup_table);
   if (!s_lookup_table.empty()) {
      std::multimap<std::string, CoarsenOperator *>::iterator mi =
         s_lookup_table.find(name);
      TBOX_ASSERT(mi != s_lookup_table.end());
      while (mi->first == name && mi->second != this) {
         ++mi;
         TBOX_ASSERT(mi != s_lookup_table.end());
      }
      TBOX_ASSERT(mi->first == name);
      TBOX_ASSERT(mi->second == this);
      mi->second = 0;
      s_lookup_table.erase(mi);
   }
   TBOX_omp_unset_lock(&l_lookup_table);
}
/*
 *************************************************************************
 * Compute the max coarsen stencil width from all constructed
 * coarsen operators.
 *************************************************************************
 */
IntVector
CoarsenOperator::getMaxCoarsenOpStencilWidth(
   const tbox::Dimension& dim)
{
   IntVector max_width(dim, 0);

   TBOX_omp_set_lock(&l_lookup_table);
   for (std::multimap<std::string, CoarsenOperator *>::const_iterator
        mi = s_lookup_table.begin(); mi != s_lookup_table.end(); ++mi) {
      const CoarsenOperator* op = mi->second;
      max_width.max(op->getStencilWidth(dim));
   }
   TBOX_omp_unset_lock(&l_lookup_table);

   return max_width;
}

/*
 *************************************************************************
 *************************************************************************
 */
void
CoarsenOperator::initializeCallback()
{
   TBOX_omp_init_lock(&l_lookup_table);
}

/*
 *************************************************************************
 *************************************************************************
 */
void
CoarsenOperator::finalizeCallback()
{
   s_lookup_table.clear();
   TBOX_omp_destroy_lock(&l_lookup_table);
}

}
}
