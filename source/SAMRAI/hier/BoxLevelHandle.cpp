/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Registry of BoxLevelHandles incident from a common BoxLevel.
 *
 ************************************************************************/
#include "SAMRAI/hier/BoxLevelHandle.h"
#include "SAMRAI/hier/BoxLevel.h"

namespace SAMRAI {
namespace hier {

/*
 ************************************************************************
 ************************************************************************
 */
BoxLevelHandle::BoxLevelHandle(
   const BoxLevel* box_level):
   d_box_level(box_level)
{
}

/*
 ************************************************************************
 ************************************************************************
 */
BoxLevelHandle::~BoxLevelHandle()
{
   detachMyBoxLevel();
}

/*
 ************************************************************************
 ************************************************************************
 */
const BoxLevel&
BoxLevelHandle::getBoxLevel() const
{
   if (d_box_level == nullptr) {
      TBOX_ERROR(
         "BoxLevelHandle::getBoxLevel Attempted to access a BoxLevel\n"
         << "that has been detached from its handle.  Detachment happens\n"
         << "when the BoxLevel changes in a way that can invalidate\n"
         << "Connector data.  Therefore, Connectors should not attempt\n"
         << "to access the BoxLevel using a detatched handle.");
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   // Sanity check: the BoxLevel's handle should be this handle.
   if (d_box_level->getBoxLevelHandle().get() != this) {
      TBOX_ERROR("Library error in BoxLevelHandle::getBoxLevel");
   }
#endif
   return *d_box_level;
}

}
}
