/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Misc patch functions used in FAC solver tests.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/pdat/ArrayDataAccess.h"
#include "QuarticFcn.h"
#include "SinusoidFcn.h"

#include "setArrayData.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/OutersideData.h"

#include <memory>

using namespace SAMRAI;

/*!
 * \file
 * \brief AMR-unaware functions to operate on a given single patch,
 * to support FAC Poisson solve.
 */

/*!
 * \brief Scale pdat::ArrayData.
 */
void scaleArrayData(
   pdat::ArrayData<double>& ad,
   double scale)
{
   if (ad.getDim() == tbox::Dimension(2)) {
      MDA_Access<double, 2,
                 MDA_OrderColMajor<2> > t4 =
         pdat::ArrayDataAccess::access<2, double>(ad);
      setArrayDataToScaled(t4,
         &ad.getBox().lower()[0],
         &ad.getBox().upper()[0],
         scale);
   } else if (ad.getDim() == tbox::Dimension(3)) {
      MDA_Access<double, 3,
                 MDA_OrderColMajor<3> > t4 =
         pdat::ArrayDataAccess::access<3, double>(ad);
      setArrayDataToScaled(t4,
         &ad.getBox().lower()[0],
         &ad.getBox().upper()[0],
         scale);
   }
}

/*!
 * \brief Set pdat::ArrayData to a constant.
 */
void setArrayDataToConstant(
   pdat::ArrayData<double>& ad,
   const geom::CartesianPatchGeometry& patch_geom,
   double value)
{
   if (ad.getDim() == tbox::Dimension(2)) {
      MDA_Access<double, 2,
                 MDA_OrderColMajor<2> > t4 =
         pdat::ArrayDataAccess::access<2, double>(ad);
      setArrayDataToConstant(t4,
         &ad.getBox().lower()[0],
         &ad.getBox().upper()[0],
         patch_geom.getXLower(),
         patch_geom.getXUpper(),
         patch_geom.getDx(),
         value);
   } else if (ad.getDim() == tbox::Dimension(3)) {
      MDA_Access<double, 3,
                 MDA_OrderColMajor<3> > t4 =
         pdat::ArrayDataAccess::access<3, double>(ad);
      setArrayDataToConstant(t4,
         &ad.getBox().lower()[0],
         &ad.getBox().upper()[0],
         patch_geom.getXLower(),
         patch_geom.getXUpper(),
         patch_geom.getDx(),
         value);
   }
}

/*!
 * \brief Set pdat::ArrayData to the x coordinate.
 */
void setArrayDataTo(
   pdat::ArrayData<double>& ad,
   const geom::CartesianPatchGeometry& patch_geom)
{
   if (ad.getDim() == tbox::Dimension(2)) {
      MDA_Access<double, 2,
                 MDA_OrderColMajor<2> > t4 =
         pdat::ArrayDataAccess::access<2, double>(ad);
      setArrayDataTo(t4,
         &ad.getBox().lower()[0],
         &ad.getBox().upper()[0],
         patch_geom.getXLower(),
         patch_geom.getXUpper(),
         patch_geom.getDx());
   } else if (ad.getDim() == tbox::Dimension(3)) {
      MDA_Access<double, 3,
                 MDA_OrderColMajor<3> > t4 =
         pdat::ArrayDataAccess::access<3, double>(ad);
      setArrayDataTo(t4,
         &ad.getBox().lower()[0],
         &ad.getBox().upper()[0],
         patch_geom.getXLower(),
         patch_geom.getXUpper(),
         patch_geom.getDx());
   }
}

/*!
 * \brief Set pdat::CellData to a sinusoid function.
 */
void setCellDataToSinusoid(
   pdat::CellData<double>& cd,
   const hier::Patch& patch,
   const SinusoidFcn& fcn)
{
   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   if (cd.getDim() == tbox::Dimension(2)) {
      MDA_Access<double, 2, MDA_OrderColMajor<2> >
      t4 = pdat::ArrayDataAccess::access<2, double>(cd.getArrayData());
      setArrayDataToSinusoid(t4,
         &cd.getGhostBox().lower()[0],
         &cd.getGhostBox().upper()[0],
         &cd.getBox().lower()[0],
         patch_geom->getXLower(),
         patch_geom->getDx(),
         fcn);
   } else if (cd.getDim() == tbox::Dimension(3)) {
      MDA_Access<double, 3, MDA_OrderColMajor<3> >
      t4 = pdat::ArrayDataAccess::access<3, double>(cd.getArrayData());
      setArrayDataToSinusoid(t4,
         &cd.getGhostBox().lower()[0],
         &cd.getGhostBox().upper()[0],
         &cd.getBox().lower()[0],
         patch_geom->getXLower(),
         patch_geom->getDx(),
         fcn);
   }
}

/*!
 * \brief Set pdat::ArrayData to Michael's exact solution.
 */
void setArrayDataToPerniceExact(
   pdat::ArrayData<double>& ad,
   const geom::CartesianPatchGeometry& patch_geom)
{
   if (ad.getDim() == tbox::Dimension(2)) {
      MDA_Access<double, 2,
                 MDA_OrderColMajor<2> > t4 =
         pdat::ArrayDataAccess::access<2, double>(ad);
      setArrayDataToPerniceExact(t4,
         &ad.getBox().lower()[0],
         &ad.getBox().upper()[0],
         patch_geom.getXLower(),
         patch_geom.getXUpper(),
         patch_geom.getDx());
   } else if (ad.getDim() == tbox::Dimension(3)) {
      MDA_Access<double, 3,
                 MDA_OrderColMajor<3> > t4 =
         pdat::ArrayDataAccess::access<3, double>(ad);
      setArrayDataToPerniceExact(t4,
         &ad.getBox().lower()[0],
         &ad.getBox().upper()[0],
         patch_geom.getXLower(),
         patch_geom.getXUpper(),
         patch_geom.getDx());
   }
}

/*!
 * \brief Set pdat::ArrayData to Michael's source function.
 */
void setArrayDataToPerniceSource(
   pdat::ArrayData<double>& ad,
   const geom::CartesianPatchGeometry& patch_geom)
{
   if (ad.getDim() == tbox::Dimension(2)) {
      MDA_Access<double, 2,
                 MDA_OrderColMajor<2> > t4 =
         pdat::ArrayDataAccess::access<2, double>(ad);
      setArrayDataToPerniceSource(t4,
         &ad.getBox().lower()[0],
         &ad.getBox().upper()[0],
         patch_geom.getXLower(),
         patch_geom.getXUpper(),
         patch_geom.getDx());
   } else if (ad.getDim() == tbox::Dimension(3)) {
      MDA_Access<double, 3,
                 MDA_OrderColMajor<3> > t4 =
         pdat::ArrayDataAccess::access<3, double>(ad);
      setArrayDataToPerniceSource(t4,
         &ad.getBox().lower()[0],
         &ad.getBox().upper()[0],
         patch_geom.getXLower(),
         patch_geom.getXUpper(),
         patch_geom.getDx());
   }
}

/*!
 * \brief Set pdat::ArrayData to a quartic function.
 */
void setCellDataToQuartic(
   pdat::CellData<double>& cd,
   const hier::Patch& patch,
   const QuarticFcn& fcn)
{
   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   if (cd.getDim() == tbox::Dimension(2)) {
      MDA_Access<double, 2, MDA_OrderColMajor<2> >
      t4 = pdat::ArrayDataAccess::access<2, double>(cd.getArrayData());
      setArrayDataToQuartic(t4,
         &cd.getGhostBox().lower()[0],
         &cd.getGhostBox().upper()[0],
         &cd.getBox().lower()[0],
         patch_geom->getXLower(),
         patch_geom->getDx(),
         fcn);
   } else if (cd.getDim() == tbox::Dimension(3)) {
      MDA_Access<double, 3, MDA_OrderColMajor<3> >
      t4 = pdat::ArrayDataAccess::access<3, double>(cd.getArrayData());
      setArrayDataToQuartic(t4,
         &cd.getGhostBox().lower()[0],
         &cd.getGhostBox().upper()[0],
         &cd.getBox().lower()[0],
         patch_geom->getXLower(),
         patch_geom->getDx(),
         fcn);
   }
}
