/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Misc patch functions used in FAC solver tests.
 *
 ************************************************************************/
#include "SinusoidFcn.h"
#include "QuarticFcn.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/OutersideData.h"

using namespace SAMRAI;

void
scaleArrayData(
   pdat::ArrayData<double>& ad,
   double scale);

void
setArrayDataToConstant(
   pdat::ArrayData<double>& ad,
   const geom::CartesianPatchGeometry& patch_geom,
   double value);

void
setArrayDataTo(
   pdat::ArrayData<double>& ad,
   const geom::CartesianPatchGeometry& patch_geom);

void
setCellDataToSinusoid(
   pdat::CellData<double>& cd,
   const hier::Patch& patch,
   const SinusoidFcn& fcn);
/*!
 * \brief Set pdat::ArrayData to Michael's exact solution.
 */
void
setArrayDataToPerniceExact(
   pdat::ArrayData<double>& ad,
   const geom::CartesianPatchGeometry& patch_geom);
/*!
 * \brief Set pdat::ArrayData to Michael's source function.
 */
void
setArrayDataToPerniceSource(
   pdat::ArrayData<double>& ad,
   const geom::CartesianPatchGeometry& patch_geom);

void
setCellDataToQuartic(
   pdat::CellData<double>& cd,
   const hier::Patch& patch,
   const QuarticFcn& fcn);
