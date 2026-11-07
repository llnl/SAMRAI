/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated operations for real node data on multiple levels.
 *
 ************************************************************************/

#ifndef included_math_HierarchyNodeDataOpsReal_C
#define included_math_HierarchyNodeDataOpsReal_C

#include "SAMRAI/math/HierarchyNodeDataOpsReal.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/pdat/NodeDataFactory.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include <typeinfo>
#include <stdlib.h>
#include <float.h>
#include <math.h>

namespace SAMRAI {
namespace math {

template<class TYPE>
HierarchyNodeDataOpsReal<TYPE>::HierarchyNodeDataOpsReal(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level):
   HierarchyDataOpsReal<TYPE>(),
   d_hierarchy(hierarchy)
{
   TBOX_ASSERT(hierarchy);

   if ((coarsest_level < 0) || (finest_level < 0)) {
      if (d_hierarchy->getNumberOfLevels() == 0) {
         d_coarsest_level = coarsest_level;
         d_finest_level = finest_level;
      } else {
         resetLevels(0, d_hierarchy->getFinestLevelNumber());
      }
   } else {
      resetLevels(coarsest_level, finest_level);
   }
}

template<class TYPE>
HierarchyNodeDataOpsReal<TYPE>::~HierarchyNodeDataOpsReal()
{
}

/*
 *************************************************************************
 *
 * Routines to set the hierarchy and level informtation.
 *
 *************************************************************************
 */

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::setPatchHierarchy(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
   TBOX_ASSERT(hierarchy);

   d_hierarchy = hierarchy;
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::resetLevels(
   const int coarsest_level,
   const int finest_level)
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((coarsest_level >= 0)
      && (finest_level >= coarsest_level)
      && (finest_level <= d_hierarchy->getFinestLevelNumber()));

   d_coarsest_level = coarsest_level;
   d_finest_level = finest_level;

   d_nonoverlapping_node_boxes.resize(d_finest_level + 1);

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      hier::BoxContainer node_boxes = level->getBoxes();
      for (hier::BoxContainer::iterator i = node_boxes.begin();
           i != node_boxes.end(); ++i) {
         *i = pdat::NodeGeometry::toNodeBox(*i);
      }
      hier::BoxUtilities::makeNonOverlappingBoxContainers(
         d_nonoverlapping_node_boxes[ln],
         node_boxes);
   }
}

template<class TYPE>
const std::shared_ptr<hier::PatchHierarchy>
HierarchyNodeDataOpsReal<TYPE>::getPatchHierarchy() const
{
   return d_hierarchy;
}

/*
 *************************************************************************
 *
 * Basic generic operations.
 *
 *************************************************************************
 */

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::copyData(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src_id)));

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.copyData(dst, src, box);
      }
   }
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::swapData(
   const int data1_id,
   const int data2_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_hierarchy);
   std::shared_ptr<pdat::NodeDataFactory<TYPE> > d1fact(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeDataFactory<TYPE>, hier::PatchDataFactory>(
         d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data1_id)));
   TBOX_ASSERT(d1fact);
   std::shared_ptr<pdat::NodeDataFactory<TYPE> > d2fact(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeDataFactory<TYPE>, hier::PatchDataFactory>(
         d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data2_id)));
   TBOX_ASSERT(d2fact);
   TBOX_ASSERT(d1fact->getDepth() == d2fact->getDepth());
   TBOX_ASSERT(d1fact->getGhostCellWidth() == d2fact->getGhostCellWidth());
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         d_patch_ops.swapData(p, data1_id, data2_id);
      }
   }
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::printData(
   const int data_id,
   std::ostream& s,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   auto factory = d_hierarchy->getPatchDescriptor()->
      getPatchDataFactory(data_id).get();
   s << "Patch descriptor id = " << data_id << std::endl;
   s << "Factory = " << typeid(*factory).name()
     << std::endl;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      s << "Level number = " << ln << std::endl;
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(data_id)));

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.printData(d, box, s);
      }
   }
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::setToScalar(
   const int data_id,
   const TYPE& alpha,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(data_id)));

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.setToScalar(d, alpha, box);
      }
   }
}

/*
 *************************************************************************
 *
 * Basic generic arithmetic operations.
 *
 *************************************************************************
 */

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::scale(
   const int dst_id,
   const TYPE& alpha,
   const int src_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src_id)));

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.scale(dst, alpha, src, box);
      }
   }
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::addScalar(
   const int dst_id,
   const int src_id,
   const TYPE& alpha,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src_id)));

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.addScalar(dst, src, alpha, box);
      }
   }
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::add(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src1(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src1_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src2(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src2_id)));

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src1);
         TBOX_ASSERT(src2);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.add(dst, src1, src2, box);
      }
   }
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::subtract(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src1(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src1_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src2(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src2_id)));

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src1);
         TBOX_ASSERT(src2);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.subtract(dst, src1, src2, box);
      }
   }
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::multiply(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src1(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src1_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src2(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src2_id)));

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src1);
         TBOX_ASSERT(src2);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.multiply(dst, src1, src2, box);
      }
   }
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::divide(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src1(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src1_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src2(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src2_id)));

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src1);
         TBOX_ASSERT(src2);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.divide(dst, src1, src2, box);
      }
   }
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::reciprocal(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src_id)));

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.reciprocal(dst, src, box);
      }
   }
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::linearSum(
   const int dst_id,
   const TYPE& alpha,
   const int src1_id,
   const TYPE& beta,
   const int src2_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src1(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src1_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src2(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src2_id)));

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src1);
         TBOX_ASSERT(src2);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.linearSum(dst, alpha, src1, beta, src2, box);
      }
   }
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::axpy(
   const int dst_id,
   const TYPE& alpha,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src1(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src1_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src2(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src2_id)));

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src1);
         TBOX_ASSERT(src2);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.axpy(dst, alpha, src1, src2, box);
      }
   }
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::axmy(
   const int dst_id,
   const TYPE& alpha,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src1(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src1_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src2(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src2_id)));

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src1);
         TBOX_ASSERT(src2);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.axmy(dst, alpha, src1, src2, box);
      }
   }
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::abs(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src_id)));

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.abs(dst, src, box);
      }
   }
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::setRandomValues(
   const int data_id,
   const TYPE& width,
   const TYPE& low,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > data(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(data_id)));

         TBOX_ASSERT(data);

         hier::Box box = (interior_only ? p->getBox() : data->getGhostBox());

         d_patch_ops.setRandomValues(data, width, low, box);
      }
   }
}

/*
 *************************************************************************
 *
 * Generic norm and order operations.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
HierarchyNodeDataOpsReal<TYPE>::numberOfEntries(
   const int data_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   size_t entries = 0;

   if (interior_only) {

      std::shared_ptr<pdat::NodeDataFactory<TYPE> > dfact(
         SAMRAI_SHARED_PTR_CAST<pdat::NodeDataFactory<TYPE>, hier::PatchDataFactory>(
            d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data_id)));

      TBOX_ASSERT(dfact);

      for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
         std::shared_ptr<hier::PatchLevel> level(
            d_hierarchy->getPatchLevel(ln));
         const int npatches = level->getNumberOfPatches();

         TBOX_ASSERT(npatches ==
            static_cast<int>(d_nonoverlapping_node_boxes[ln].size()));

         for (int il = 0; il < npatches; ++il) {
            hier::BoxContainer::const_iterator lb =
               ((d_nonoverlapping_node_boxes[ln])[il]).begin();
            for ( ; lb != ((d_nonoverlapping_node_boxes[ln])[il]).end(); ++lb) {
               entries += lb->size();
            }
         }
      }

      entries *= dfact->getDepth();

   } else {

      for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
         std::shared_ptr<hier::PatchLevel> level(
            d_hierarchy->getPatchLevel(ln));
         for (hier::PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            std::shared_ptr<pdat::NodeData<TYPE> > d(
               SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
                  (*ip)->getPatchData(data_id)));

            TBOX_ASSERT(d);

            entries += d_patch_ops.numberOfEntries(d, d->getGhostBox());
         }
      }

      unsigned long int global_entries = entries;
      if (mpi.getSize() > 1) {
         mpi.Allreduce(&entries, &global_entries, 1, MPI_UNSIGNED_LONG, MPI_SUM);
      }
      entries = global_entries;

   }

   return entries;
}

template<class TYPE>
double
HierarchyNodeDataOpsReal<TYPE>::sumControlVolumes(
   const int data_id,
   const int vol_id) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   double sum = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > data(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(data_id)));
         std::shared_ptr<pdat::NodeData<double> > cv(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
               p->getPatchData(vol_id)));

         TBOX_ASSERT(data);
         TBOX_ASSERT(cv);

         hier::Box box = cv->getGhostBox();

         sum += d_patch_ops.sumControlVolumes(data, cv, box);
      }
   }

   double global_sum = sum;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM);
   }
   return global_sum;
}

template<class TYPE>
double
HierarchyNodeDataOpsReal<TYPE>::L1Norm(
   const int data_id,
   const int vol_id,
   bool local_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   double norm = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(data_id)));
         std::shared_ptr<hier::PatchData> pd;

         TBOX_ASSERT(d);

         hier::Box box = p->getBox();
         if (vol_id >= 0) {

            box = d->getGhostBox();
            pd = p->getPatchData(vol_id);
         }

         std::shared_ptr<pdat::NodeData<double> > cv(
            std::dynamic_pointer_cast<pdat::NodeData<double>,
                                        hier::PatchData>(pd));
         norm += d_patch_ops.L1Norm(d, box, cv);
      }
   }

   if (!local_only) {
      double global_norm = norm;
      if (mpi.getSize() > 1) {
         mpi.Allreduce(&norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM);
      }
      norm = global_norm;
   }
   return norm;
}

template<class TYPE>
double
HierarchyNodeDataOpsReal<TYPE>::L2Norm(
   const int data_id,
   const int vol_id,
   bool local_only) const
{
   double norm_squared = HierarchyNodeDataOpsReal<TYPE>::dot(data_id,
         data_id,
         vol_id,
         local_only);

   return sqrt(norm_squared);
}

template<class TYPE>
double
HierarchyNodeDataOpsReal<TYPE>::weightedL2Norm(
   const int data_id,
   const int wgt_id,
   const int vol_id) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   double norm_squared = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(data_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > w(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(wgt_id)));
         std::shared_ptr<hier::PatchData> pd;

         TBOX_ASSERT(d);
         TBOX_ASSERT(w);

         hier::Box box = p->getBox();
         if (vol_id >= 0) {

            box = d->getGhostBox();
            pd = p->getPatchData(vol_id);
         }

         std::shared_ptr<pdat::NodeData<double> > cv(
            std::dynamic_pointer_cast<pdat::NodeData<double>,
                                        hier::PatchData>(pd));
         double pnorm = d_patch_ops.weightedL2Norm(d, w, box, cv);

         norm_squared += pnorm * pnorm;
      }
   }

   double global_norm_squared = norm_squared;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&norm_squared, &global_norm_squared, 1, MPI_DOUBLE, MPI_SUM);
   }
   return sqrt(global_norm_squared);
}

template<class TYPE>
double
HierarchyNodeDataOpsReal<TYPE>::RMSNorm(
   const int data_id,
   const int vol_id) const
{
   double l2_norm = L2Norm(data_id, vol_id);

   double volume = ((vol_id < 0) ? (double)numberOfEntries(data_id, true)
                    : sumControlVolumes(data_id, vol_id));

   double rms_norm = l2_norm / sqrt(volume);
   return rms_norm;
}

template<class TYPE>
double
HierarchyNodeDataOpsReal<TYPE>::weightedRMSNorm(
   const int data_id,
   const int wgt_id,
   const int vol_id) const
{

   double l2_norm = weightedL2Norm(data_id, wgt_id, vol_id);

   double volume = ((vol_id < 0) ? (double)numberOfEntries(data_id, true)
                    : sumControlVolumes(data_id, vol_id));

   double rms_norm = l2_norm / sqrt(volume);
   return rms_norm;
}

template<class TYPE>
double
HierarchyNodeDataOpsReal<TYPE>::maxNorm(
   const int data_id,
   const int vol_id,
   bool local_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   double norm = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(data_id)));
         std::shared_ptr<hier::PatchData> pd;

         TBOX_ASSERT(d);

         hier::Box box = p->getBox();
         if (vol_id >= 0) {

            box = d->getGhostBox();
            pd = p->getPatchData(vol_id);
         }

         std::shared_ptr<pdat::NodeData<double> > cv(
            std::dynamic_pointer_cast<pdat::NodeData<double>,
                                        hier::PatchData>(pd));
         norm = tbox::MathUtilities<double>::Max(norm,
               d_patch_ops.maxNorm(d, box, cv));
      }
   }

   if (!local_only) {
      double global_norm = norm;
      if (mpi.getSize() > 1) {
         mpi.Allreduce(&norm, &global_norm, 1, MPI_DOUBLE, MPI_MAX);
      }
      norm = global_norm;
   }
   return norm;
}

template<class TYPE>
TYPE
HierarchyNodeDataOpsReal<TYPE>::dot(
   const int data1_id,
   const int data2_id,
   const int vol_id,
   bool local_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   TYPE dprod = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > d1(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(data1_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > d2(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(data2_id)));
         std::shared_ptr<hier::PatchData> pd;

         TBOX_ASSERT(d1);
         TBOX_ASSERT(d2);

         hier::Box box = p->getBox();
         if (vol_id >= 0) {

            box = d1->getGhostBox();
            pd = p->getPatchData(vol_id);
         }

         std::shared_ptr<pdat::NodeData<double> > cv(
            std::dynamic_pointer_cast<pdat::NodeData<double>,
                                        hier::PatchData>(pd));
         dprod += d_patch_ops.dot(d1, d2, box, cv);
      }
   }

   if (!local_only) {
      if (mpi.getSize() > 1) {
         mpi.AllReduce(&dprod, 1, MPI_SUM);
      }
   }
   return dprod;
}

template<class TYPE>
TYPE
HierarchyNodeDataOpsReal<TYPE>::integral(
   const int data_id,
   const int vol_id) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   TYPE local_integral = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > data(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(data_id)));
         std::shared_ptr<pdat::NodeData<double> > vol(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
               p->getPatchData(vol_id)));

         TBOX_ASSERT(data);
         TBOX_ASSERT(vol);

         hier::Box box = data->getGhostBox();

         local_integral += d_patch_ops.integral(data, box, vol);
      }
   }

   TYPE global_integral = local_integral;
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&global_integral, 1, MPI_SUM);
   }
   return global_integral;
}

/*
 *************************************************************************
 *
 * Generic miscellaneous operations for real data.
 *
 *************************************************************************
 */

template<class TYPE>
int
HierarchyNodeDataOpsReal<TYPE>::computeConstrProdPos(
   const int data1_id,
   const int data2_id,
   const int vol_id) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   int test = 1;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > d1(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(data1_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > d2(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(data2_id)));
         std::shared_ptr<hier::PatchData> pd;

         TBOX_ASSERT(d1);
         TBOX_ASSERT(d2);

         hier::Box box = p->getBox();
         if (vol_id >= 0) {

            box = d1->getGhostBox();
            pd = p->getPatchData(vol_id);
         }

         std::shared_ptr<pdat::NodeData<double> > cv(
            std::dynamic_pointer_cast<pdat::NodeData<double>,
                                        hier::PatchData>(pd));
         test = tbox::MathUtilities<int>::Min(test,
               d_patch_ops.computeConstrProdPos(d1, d2, box, cv));
      }
   }

   int global_test = test;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&test, &global_test, 1, MPI_INT, MPI_MIN);
   }
   return global_test;
}

template<class TYPE>
void
HierarchyNodeDataOpsReal<TYPE>::compareToScalar(
   const int dst_id,
   const int src_id,
   const TYPE& alpha,
   const int vol_id) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src_id)));
         std::shared_ptr<hier::PatchData> pd;

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src);

         hier::Box box = p->getBox();
         if (vol_id >= 0) {

            box = dst->getGhostBox();
            pd = p->getPatchData(vol_id);
         }

         std::shared_ptr<pdat::NodeData<double> > cv(
            std::dynamic_pointer_cast<pdat::NodeData<double>,
                                        hier::PatchData>(pd));
         d_patch_ops.compareToScalar(dst, src, alpha, box, cv);
      }
   }
}

template<class TYPE>
int
HierarchyNodeDataOpsReal<TYPE>::testReciprocal(
   const int dst_id,
   const int src_id,
   const int vol_id) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   int test = 1;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > src(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(src_id)));
         std::shared_ptr<hier::PatchData> pd;

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src);

         hier::Box box = p->getBox();
         if (vol_id >= 0) {

            box = dst->getGhostBox();
            pd = p->getPatchData(vol_id);
         }

         std::shared_ptr<pdat::NodeData<double> > cv(
            std::dynamic_pointer_cast<pdat::NodeData<double>,
                                        hier::PatchData>(pd));
         test = tbox::MathUtilities<int>::Min(test,
               d_patch_ops.testReciprocal(dst, src, box, cv));
      }
   }

   int global_test = test;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&test, &global_test, 1, MPI_INT, MPI_MIN);
   }
   return global_test;
}

template<class TYPE>
TYPE
HierarchyNodeDataOpsReal<TYPE>::maxPointwiseDivide(
   const int numer_id,
   const int denom_id,
   bool local_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   TYPE max = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > numer(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(numer_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > denom(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(denom_id)));

         TBOX_ASSERT(numer);
         TBOX_ASSERT(denom);

         hier::Box box = p->getBox();

         max = tbox::MathUtilities<TYPE>::Max(max,
               d_patch_ops.maxPointwiseDivide(numer, denom, box));
      }
   }

   if (!local_only) {
      if (mpi.getSize() > 1) {
         mpi.AllReduce(&max, 1, MPI_MAX);
      }
   }
   return max;
}

template<class TYPE>
TYPE
HierarchyNodeDataOpsReal<TYPE>::minPointwiseDivide(
   const int numer_id,
   const int denom_id,
   bool local_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   TYPE min = tbox::MathUtilities<TYPE>::getMax();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > numer(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(numer_id)));
         std::shared_ptr<pdat::NodeData<TYPE> > denom(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(denom_id)));

         TBOX_ASSERT(numer);
         TBOX_ASSERT(denom);

         hier::Box box = p->getBox();

         min = tbox::MathUtilities<TYPE>::Min(min,
               d_patch_ops.minPointwiseDivide(numer, denom, box));
      }
   }

   if (!local_only) {
      if (mpi.getSize() > 1) {
         mpi.AllReduce(&min, 1, MPI_MIN);
      }
   }
   return min;
}

template<class TYPE>
TYPE
HierarchyNodeDataOpsReal<TYPE>::min(
   const int data_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   TYPE minval = tbox::MathUtilities<TYPE>::getMax();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(data_id)));

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         minval = tbox::MathUtilities<TYPE>::Min(minval, d_patch_ops.min(d, box));
      }
   }

   TYPE global_min = minval;
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&global_min, 1, MPI_MIN);
   }
   return global_min;
}

template<class TYPE>
TYPE
HierarchyNodeDataOpsReal<TYPE>::max(
   const int data_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(data_id)));

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         maxval = tbox::MathUtilities<TYPE>::Max(maxval, d_patch_ops.max(d, box));
      }
   }

   TYPE global_max = maxval;
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&global_max, 1, MPI_MAX);
   }
   return global_max;
}

template<class TYPE>
int64_t HierarchyNodeDataOpsReal<TYPE>::getLength(
   const int data_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   int64_t length = 0;
   hier::Box data_box(d_hierarchy->getDim());

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::NodeData<TYPE> > data(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
               p->getPatchData(data_id)));

         TBOX_ASSERT(data);

         if (interior_only) {
            data_box = pdat::NodeGeometry::toNodeBox(data->getBox());
         } else {
            data_box = data->getArrayData().getBox();
         }
         length += static_cast<int64_t>(data_box.size() * data->getDepth());
      }
   }

   int64_t global_length = length;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&length, &global_length, 1, MPI_INT64_T, MPI_SUM);
   }
   return global_length;
}


}
}
#endif
