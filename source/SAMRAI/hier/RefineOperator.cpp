/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract base class for spatial refinement operators.
 *
 ************************************************************************/
#include "SAMRAI/hier/RefineOperator.h"

#include "SAMRAI/tbox/StartupShutdownManager.h"

#include "SAMRAI/tbox/OpenMPUtilities.h"

namespace SAMRAI {
namespace hier {

std::multimap<std::string, RefineOperator *> RefineOperator::s_lookup_table;
TBOX_omp_lock_t RefineOperator::l_lookup_table;

tbox::StartupShutdownManager::Handler
RefineOperator::s_finalize_handler(
   RefineOperator::initializeCallback,
   0,
   0,
   RefineOperator::finalizeCallback,
   tbox::StartupShutdownManager::priorityList);

RefineOperator::RefineOperator(
   const std::string& name):
   d_name(name)
{
   registerInLookupTable(name);
}

RefineOperator::~RefineOperator()
{
   removeFromLookupTable(d_name);
}

void
RefineOperator::registerInLookupTable(
   const std::string& name)
{
   TBOX_omp_set_lock(&l_lookup_table);
   s_lookup_table.insert(
      std::pair<std::string, RefineOperator *>(name, this));
   TBOX_omp_unset_lock(&l_lookup_table);
}

void
RefineOperator::removeFromLookupTable(
   const std::string& name)
{
   /*
    * The lookup table might be empty if static RefineOperator's are used
    * in which case the table will have been cleared before the statics
    * are destroyed.
    */
   TBOX_omp_set_lock(&l_lookup_table);
   if (!s_lookup_table.empty()) {
      std::multimap<std::string, RefineOperator *>::iterator mi =
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
 * Compute the max refine stencil width from all constructed
 * refine operators.
 *************************************************************************
 */
IntVector
RefineOperator::getMaxRefineOpStencilWidth(
   const tbox::Dimension& dim)
{
   IntVector max_width(dim, 0);

   TBOX_omp_set_lock(&l_lookup_table);
   for (std::multimap<std::string, RefineOperator *>::const_iterator
        mi = s_lookup_table.begin(); mi != s_lookup_table.end(); ++mi) {
      const RefineOperator* op = mi->second;
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
RefineOperator::initializeCallback()
{
   TBOX_omp_init_lock(&l_lookup_table);
}

/*
 *************************************************************************
 *************************************************************************
 */
void
RefineOperator::finalizeCallback()
{
   s_lookup_table.clear();
   TBOX_omp_destroy_lock(&l_lookup_table);
}

}
}
