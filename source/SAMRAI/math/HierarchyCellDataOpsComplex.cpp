/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex cell data on multiple levels.
 *
 ************************************************************************/
#include "SAMRAI/math/HierarchyCellDataOpsComplex.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include <typeinfo>
#include <stdlib.h>
#include <float.h>
#include <math.h>

namespace SAMRAI {
namespace math {

HierarchyCellDataOpsComplex::HierarchyCellDataOpsComplex(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level):
   HierarchyDataOpsComplex(),
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

HierarchyCellDataOpsComplex::~HierarchyCellDataOpsComplex()
{
}

/*
 *************************************************************************
 *
 * Routines to set the hierarchy and level information.
 *
 *************************************************************************
 */

void
HierarchyCellDataOpsComplex::setPatchHierarchy(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
   TBOX_ASSERT(hierarchy);

   d_hierarchy = hierarchy;
}

void
HierarchyCellDataOpsComplex::resetLevels(
   const int coarsest_level,
   const int finest_level)
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((coarsest_level >= 0)
      && (finest_level >= coarsest_level)
      && (finest_level <= d_hierarchy->getFinestLevelNumber()));

   d_coarsest_level = coarsest_level;
   d_finest_level = finest_level;
}

const std::shared_ptr<hier::PatchHierarchy>
HierarchyCellDataOpsComplex::getPatchHierarchy() const
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

void
HierarchyCellDataOpsComplex::copyData(
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

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > s(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src_id)));

         TBOX_ASSERT(d);
         TBOX_ASSERT(s);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.copyData(d, s, box);
      }
   }
}

void
HierarchyCellDataOpsComplex::swapData(
   const int data1_id,
   const int data2_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   std::shared_ptr<pdat::CellDataFactory<dcomplex> > d1fact(
      SAMRAI_SHARED_PTR_CAST<pdat::CellDataFactory<dcomplex>, hier::PatchDataFactory>(
         d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data1_id)));
   TBOX_ASSERT(d1fact);
   std::shared_ptr<pdat::CellDataFactory<dcomplex> > d2fact(
      SAMRAI_SHARED_PTR_CAST<pdat::CellDataFactory<dcomplex>, hier::PatchDataFactory>(
         d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data2_id)));
   TBOX_ASSERT(d2fact);
   TBOX_ASSERT(d1fact->getDepth() == d2fact->getDepth());
   TBOX_ASSERT(d1fact->getGhostCellWidth() == d2fact->getGhostCellWidth());
#endif

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

         d_patch_ops.swapData(p, data1_id, data2_id);
      }
   }
}

void
HierarchyCellDataOpsComplex::printData(
   const int data_id,
   std::ostream& s,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   auto& pdf = *d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data_id);
   s << "Patch descriptor id = " << data_id << std::endl;
   s << "Factory = " << typeid(pdf).name() << std::endl;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      s << "Level number = " << ln << std::endl;
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(data_id)));

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.printData(d, box, s);
      }
   }
}

void
HierarchyCellDataOpsComplex::setToScalar(
   const int data_id,
   const dcomplex& alpha,
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

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
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

void
HierarchyCellDataOpsComplex::scale(
   const int dst_id,
   const dcomplex& alpha,
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

         std::shared_ptr<pdat::CellData<dcomplex> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > src(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src_id)));

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.scale(dst, alpha, src, box);
      }
   }
}

void
HierarchyCellDataOpsComplex::addScalar(
   const int dst_id,
   const int src_id,
   const dcomplex& alpha,
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

         std::shared_ptr<pdat::CellData<dcomplex> > dst(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > src(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src_id)));

         TBOX_ASSERT(dst);
         TBOX_ASSERT(src);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.addScalar(dst, src, alpha, box);
      }
   }
}

void
HierarchyCellDataOpsComplex::add(
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

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > s1(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src1_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > s2(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src2_id)));

         TBOX_ASSERT(d);
         TBOX_ASSERT(s1);
         TBOX_ASSERT(s2);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.add(d, s1, s2, box);
      }
   }
}

void
HierarchyCellDataOpsComplex::subtract(
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

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > s1(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src1_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > s2(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src2_id)));

         TBOX_ASSERT(d);
         TBOX_ASSERT(s1);
         TBOX_ASSERT(s2);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.subtract(d, s1, s2, box);
      }
   }
}

void
HierarchyCellDataOpsComplex::multiply(
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

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > s1(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src1_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > s2(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src2_id)));

         TBOX_ASSERT(d);
         TBOX_ASSERT(s1);
         TBOX_ASSERT(s2);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.multiply(d, s1, s2, box);
      }
   }
}

void
HierarchyCellDataOpsComplex::divide(
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

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > s1(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src1_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > s2(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src2_id)));

         TBOX_ASSERT(d);
         TBOX_ASSERT(s1);
         TBOX_ASSERT(s2);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.divide(d, s1, s2, box);
      }
   }
}

void
HierarchyCellDataOpsComplex::reciprocal(
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

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > src(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src_id)));

         TBOX_ASSERT(d);
         TBOX_ASSERT(src);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.reciprocal(d, src, box);
      }
   }
}

void
HierarchyCellDataOpsComplex::linearSum(
   const int dst_id,
   const dcomplex& alpha,
   const int src1_id,
   const dcomplex& beta,
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

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > s1(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src1_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > s2(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src2_id)));

         TBOX_ASSERT(d);
         TBOX_ASSERT(s1);
         TBOX_ASSERT(s2);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.linearSum(d, alpha, s1, beta, s2, box);
      }
   }
}

void
HierarchyCellDataOpsComplex::axpy(
   const int dst_id,
   const dcomplex& alpha,
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

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > s1(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src1_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > s2(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src2_id)));

         TBOX_ASSERT(d);
         TBOX_ASSERT(s1);
         TBOX_ASSERT(s2);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.axpy(d, alpha, s1, s2, box);
      }
   }
}

void
HierarchyCellDataOpsComplex::axmy(
   const int dst_id,
   const dcomplex& alpha,
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

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > s1(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src1_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > s2(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src2_id)));

         TBOX_ASSERT(d);
         TBOX_ASSERT(s1);
         TBOX_ASSERT(s2);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.axmy(d, alpha, s1, s2, box);
      }
   }
}

void
HierarchyCellDataOpsComplex::abs(
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

         std::shared_ptr<pdat::CellData<double> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               p->getPatchData(dst_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > src(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(src_id)));

         TBOX_ASSERT(d);
         TBOX_ASSERT(src);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.abs(d, src, box);
      }
   }
}

void
HierarchyCellDataOpsComplex::setRandomValues(
   const int data_id,
   const dcomplex& width,
   const dcomplex& low,
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

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(data_id)));

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.setRandomValues(d, width, low, box);
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

size_t
HierarchyCellDataOpsComplex::numberOfEntries(
   const int data_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   size_t entries = 0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(data_id)));

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         if (d) {
            entries += d_patch_ops.numberOfEntries(d, box);
         }
      }
   }

   unsigned long int global_entries = entries;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&entries, &global_entries, 1, MPI_UNSIGNED_LONG, MPI_SUM);
   }
   return global_entries;
}

double
HierarchyCellDataOpsComplex::sumControlVolumes(
   const int data_id,
   const int vol_id) const
{
   TBOX_ASSERT(vol_id >= 0);
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

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(data_id)));
         std::shared_ptr<pdat::CellData<double> > cv(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               p->getPatchData(vol_id)));

         TBOX_ASSERT(d);
         TBOX_ASSERT(cv);

         hier::Box box = cv->getGhostBox();

         sum += d_patch_ops.sumControlVolumes(d, cv, box);
      }
   }

   double global_sum = sum;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM);
   }
   return global_sum;
}

double
HierarchyCellDataOpsComplex::L1Norm(
   const int data_id,
   const int vol_id) const
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

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(data_id)));
         std::shared_ptr<hier::PatchData> pd;

         TBOX_ASSERT(d);

         hier::Box box = p->getBox();
         if (vol_id >= 0) {

            box = d->getGhostBox();
            pd = p->getPatchData(vol_id);
         }

         std::shared_ptr<pdat::CellData<double> > cv(
            std::dynamic_pointer_cast<pdat::CellData<double>,
                                        hier::PatchData>(pd));
         norm += d_patch_ops.L1Norm(d, box, cv);
      }
   }

   double global_norm = norm;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM);
   }
   return global_norm;
}

double
HierarchyCellDataOpsComplex::L2Norm(
   const int data_id,
   const int vol_id) const
{
   dcomplex dotprod = HierarchyCellDataOpsComplex::dot(data_id,
         data_id,
         vol_id);

   return sqrt(real(dotprod));
}

double
HierarchyCellDataOpsComplex::weightedL2Norm(
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

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(data_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > w(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(wgt_id)));
         std::shared_ptr<hier::PatchData> pd;

         TBOX_ASSERT(d);
         TBOX_ASSERT(w);

         hier::Box box = p->getBox();
         if (vol_id >= 0) {

            box = d->getGhostBox();
            pd = p->getPatchData(vol_id);
         }

         std::shared_ptr<pdat::CellData<double> > cv(
            std::dynamic_pointer_cast<pdat::CellData<double>,
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

double
HierarchyCellDataOpsComplex::RMSNorm(
   const int data_id,
   const int vol_id) const
{
   double l2_norm = L2Norm(data_id, vol_id);

   double volume = ((vol_id < 0) ? (double)numberOfEntries(data_id, true)
                    : sumControlVolumes(data_id, vol_id));

   double rms_norm = l2_norm / sqrt(volume);
   return rms_norm;
}

double
HierarchyCellDataOpsComplex::weightedRMSNorm(
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

double
HierarchyCellDataOpsComplex::maxNorm(
   const int data_id,
   const int vol_id) const
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

         std::shared_ptr<pdat::CellData<dcomplex> > d(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(data_id)));
         std::shared_ptr<hier::PatchData> pd;

         TBOX_ASSERT(d);

         hier::Box box = p->getBox();
         if (vol_id >= 0) {

            box = d->getGhostBox();
            pd = p->getPatchData(vol_id);
         }

         std::shared_ptr<pdat::CellData<double> > cv(
            std::dynamic_pointer_cast<pdat::CellData<double>,
                                        hier::PatchData>(pd));
         norm = tbox::MathUtilities<double>::Max(norm,
               d_patch_ops.maxNorm(d, box, cv));
      }
   }

   double global_norm = norm;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&norm, &global_norm, 1, MPI_DOUBLE, MPI_MAX);
   }
   return global_norm;
}

dcomplex
HierarchyCellDataOpsComplex::dot(
   const int data1_id,
   const int data2_id,
   const int vol_id) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   dcomplex dprod = dcomplex(0.0, 0.0);

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::CellData<dcomplex> > d1(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(data1_id)));
         std::shared_ptr<pdat::CellData<dcomplex> > d2(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(data2_id)));
         std::shared_ptr<hier::PatchData> pd;

         TBOX_ASSERT(d1);
         TBOX_ASSERT(d2);

         hier::Box box = p->getBox();
         if (vol_id >= 0) {

            box = d1->getGhostBox();
            pd = p->getPatchData(vol_id);
         }

         std::shared_ptr<pdat::CellData<double> > cv(
            std::dynamic_pointer_cast<pdat::CellData<double>,
                                        hier::PatchData>(pd));
         dprod += d_patch_ops.dot(d1, d2, box, cv);
      }
   }

   if (mpi.getSize() > 1) {
      // It should be possible to do this with a single Allreduce and a
      // datatype of MPI_C_DOUBLE_COMPLEX.  However, while recent versions of
      // openmpi define this datatype their implementations of Allreduce do not
      // recognize it.
      double real_part = dprod.real();
      double imag_part = dprod.imag();
      double global_real_part;
      double global_imag_part;
      mpi.Allreduce(&real_part, &global_real_part, 1, MPI_DOUBLE, MPI_SUM);
      mpi.Allreduce(&imag_part, &global_imag_part, 1, MPI_DOUBLE, MPI_SUM);
      dcomplex global_dot(global_real_part, global_imag_part);
      return global_dot;
   } else {
      return dprod;
   }
}

dcomplex
HierarchyCellDataOpsComplex::integral(
   const int data_id,
   const int vol_id) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   dcomplex local_integral = dcomplex(0.0, 0.0);

   for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::CellData<dcomplex> > data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               p->getPatchData(data_id)));
         std::shared_ptr<pdat::CellData<double> > vol(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               p->getPatchData(vol_id)));

         TBOX_ASSERT(data);
         TBOX_ASSERT(vol);

         hier::Box box = data->getGhostBox();

         local_integral += d_patch_ops.integral(data, box, vol);
      }
   }

   if (mpi.getSize() > 1) {
      // It should be possible to do this with a single Allreduce and a
      // datatype of MPI_C_DOUBLE_COMPLEX.  However, while recent versions of
      // openmpi define this datatype their implementations of Allreduce do not
      // recognize it.
      double real_part = local_integral.real();
      double imag_part = local_integral.imag();
      double global_real_part;
      double global_imag_part;
      mpi.Allreduce(&real_part, &global_real_part, 1, MPI_DOUBLE, MPI_SUM);
      mpi.Allreduce(&imag_part, &global_imag_part, 1, MPI_DOUBLE, MPI_SUM);
      dcomplex global_integral(global_real_part, global_imag_part);
      return global_integral;
   } else {
      return local_integral;
   }
}

}
}
