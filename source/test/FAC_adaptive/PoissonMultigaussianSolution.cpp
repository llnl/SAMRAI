/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   PoissonMultigaussianSolution class implementation
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/ArrayDataAccess.h"
#include "patchFcns.h"
#include "PoissonMultigaussianSolution.h"
#include STL_SSTREAM_HEADER_FILE

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include "SAMRAI/tbox/Utilities.h"

using namespace SAMRAI;

PoissonMultigaussianSolution::PoissonMultigaussianSolution(
   const tbox::Dimension& dim)
#ifndef PACKAGE
   :d_dim(dim)
#endif
{
}

PoissonMultigaussianSolution::PoissonMultigaussianSolution(
   const std::string& object_name
   ,
   const tbox::Dimension& dim
   ,
   tbox::Database& database
   , /*! Standard output stream */
   std::ostream* out_stream
   , /*! Log output stream */
   std::ostream* log_stream)
#ifndef PACKAGE
   :d_dim(dim)
#endif
{
   NULL_USE(object_name);
   NULL_USE(out_stream);
   NULL_USE(log_stream);

   setFromDatabase(database);
}

PoissonMultigaussianSolution::~PoissonMultigaussianSolution()
{
}

void PoissonMultigaussianSolution::setFromDatabase(
   tbox::Database& database)
{
   std::string singlegauss = "GaussianFcnControl_"
      + tbox::Utilities::intToString(static_cast<int>(d_gauss_size));
   if (!database.isString(singlegauss)) {
      TBOX_ERROR(
         "You must have at least " << singlegauss << " defined in the\n"
         "database for PoissonMultigaussianSolution.\n");
   }
   GaussianFcn gauss(d_dim);
   do {
      std::string istr = database.getString(singlegauss);
      std::istringstream ist(istr);
      ist >> gauss;
      d_gauss.push_back(gauss);
      singlegauss = "GaussianFcnControl_"
         + tbox::Utilities::intToString(static_cast<int>(d_gauss_size));
   } while (database.isString(singlegauss));
}

void PoissonMultigaussianSolution::setPoissonSpecifications(
   solv::PoissonSpecifications& sps,
   int C_patch_data_id,
   int D_patch_data_id) const
{
   NULL_USE(C_patch_data_id);
   NULL_USE(D_patch_data_id);

   sps.setDConstant(1.0);
   sps.setCZero();
}

double PoissonMultigaussianSolution::exactFcn(
   double x,
   double y) const {
   double rval = 0;
   d_gauss_const_iterator i;
   for (i = d_gauss_begin; i != d_gauss_end; ++i) {
      rval += (*i)(x, y);
   }
   return rval;
}

double PoissonMultigaussianSolution::exactFcn(
   double x,
   double y,
   double z) const {
   double rval = 0;
   d_gauss_const_iterator i;
   for (i = d_gauss_begin; i != d_gauss_end; ++i) {
      rval += (*i)(x, y, z);
   }
   return rval;
}

double PoissonMultigaussianSolution::sourceFcn(
   double x,
   double y) const {
   double rval = 0;
   d_gauss_const_iterator i;
   for (i = d_gauss_begin; i != d_gauss_end; ++i) {
      const GaussianFcn& gauss = *i;
      double gauss_ctr[SAMRAI::MAX_DIM_VAL];
      gauss.getCenter(gauss_ctr);
      double tval;
      tval = 4 * gauss.getLambda() * ((x - gauss_ctr[0]) * (x - gauss_ctr[0])
                                      + (y - gauss_ctr[1]) * (y - gauss_ctr[1])
                                      );
      tval += 2 * d_dim.getValue();
      tval *= gauss(x, y) * gauss.getLambda();
      rval += tval;
   }
   return rval;
}

double PoissonMultigaussianSolution::sourceFcn(
   double x,
   double y,
   double z) const {
   double rval = 0;
   d_gauss_const_iterator i;
   for (i = d_gauss_begin; i != d_gauss_end; ++i) {
      const GaussianFcn& gauss = *i;
      double gauss_ctr[SAMRAI::MAX_DIM_VAL];
      gauss.getCenter(gauss_ctr);
      double tval;
      tval = 4 * gauss.getLambda() * ((x - gauss_ctr[0]) * (x - gauss_ctr[0])
                                      + (y - gauss_ctr[1]) * (y - gauss_ctr[1])
                                      + (z - gauss_ctr[2]) * (z - gauss_ctr[2])
                                      );
      tval += 2 * d_dim.getValue();
      tval *= gauss(x, y, z) * gauss.getLambda();
      rval += tval;
   }
   return rval;
}

void PoissonMultigaussianSolution::setGridData(
   hier::Patch& patch,
   pdat::CellData<double>& exact_data,
   pdat::CellData<double>& source_data)
{
   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);

   const double* h = patch_geom->getDx();
   const double* xl = patch_geom->getXLower();
   const int* il = &patch.getBox().lower()[0];

   {
      /* Set cell-centered data. */
      double sl[SAMRAI::MAX_DIM_VAL]; // Like XLower, except for cell.
      int j;
      for (j = 0; j < d_dim.getValue(); ++j) {
         sl[j] = xl[j] + 0.5 * h[j];
      }
      pdat::CellData<double>::iterator iter(pdat::CellGeometry::begin(patch.getBox()));
      pdat::CellData<double>::iterator iterend(pdat::CellGeometry::end(patch.getBox()));
      if (d_dim == tbox::Dimension(2)) {
         double x, y;
         for ( ; iter != iterend; ++iter) {
            const pdat::CellIndex& index = *iter;
            x = sl[0] + (index[0] - il[0]) * h[0];
            y = sl[1] + (index[1] - il[1]) * h[1];
            exact_data(index) = exactFcn(x, y);
            source_data(index) = sourceFcn(x, y);
         }
      } else if (d_dim == tbox::Dimension(3)) {
         double x, y, z;
         for ( ; iter != iterend; ++iter) {
            const pdat::CellIndex& index = *iter;
            x = sl[0] + (index[0] - il[0]) * h[0];
            y = sl[1] + (index[1] - il[1]) * h[1];
            z = sl[2] + (index[2] - il[2]) * h[2];
            exact_data(index) = exactFcn(x, y, z);
            source_data(index) = sourceFcn(x, y, z);
         }
      }
   }
}       // End patch loop.

std::ostream& operator << (
   std::ostream& os,
   const PoissonMultigaussianSolution& r) {
   os << "{\n";
   for (unsigned int i = 0; i < r.d_gauss_size; ++i) {
      os << "GaussianFcnControl_" << i << " " << r.d_gauss[i] << "\n";
   }
   os << "}\n";
   return os;
}

void PoissonMultigaussianSolution::setBcCoefs(
   const std::shared_ptr<pdat::ArrayData<double> >& acoef_data,
   const std::shared_ptr<pdat::ArrayData<double> >& bcoef_data,
   const std::shared_ptr<pdat::ArrayData<double> >& gcoef_data,
   const std::shared_ptr<hier::Variable>& variable,
   const hier::Patch& patch,
   const hier::BoundaryBox& bdry_box,
   const double fill_time) const
{
   NULL_USE(variable);
   NULL_USE(fill_time);

   if (!acoef_data && !gcoef_data) {
      return;
   }

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   /*
    * Set to an inhomogeneous Dirichlet boundary condition.
    */
   hier::Box patch_box(patch.getBox());

   const double* xlo = patch_geom->getXLower();
   const double* xup = patch_geom->getXUpper();
   const double* dx = patch_geom->getDx();

   if (bdry_box.getBoundaryType() != 1) {
      // Must be a face boundary.
      TBOX_ERROR("Bad boundary type in\n"
         << "PoissonMultigaussianSolution::setBcCoefs \n");
   }
   const hier::Box& box = bdry_box.getBox();
   hier::Index lower = box.lower();
   hier::Index upper = box.upper();

   if (d_dim == tbox::Dimension(2)) {
      hier::Box::iterator boxit(acoef_data ?
                                acoef_data->getBox().begin() :
                                gcoef_data->getBox().begin());
      hier::Box::iterator boxitend(acoef_data ?
                                   acoef_data->getBox().end() :
                                   gcoef_data->getBox().end());
      int i, j;
      double x, y;
      switch (bdry_box.getLocationIndex()) {
         case 0:
            // min i edge
            x = xlo[0];
            for ( ; boxit != boxitend; ++boxit) {
               j = (*boxit)[1];
               y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
               if (acoef_data) (*acoef_data)(*boxit, 0) = 1.0;
               if (bcoef_data) (*bcoef_data)(*boxit, 0) = 0.0;
               if (gcoef_data) (*gcoef_data)(*boxit, 0) = exactFcn(x, y);
            }
            break;
         case 1:
            // max i edge
            x = xup[0];
            for ( ; boxit != boxitend; ++boxit) {
               j = (*boxit)[1];
               y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
               if (acoef_data) (*acoef_data)(*boxit, 0) = 1.0;
               if (bcoef_data) (*bcoef_data)(*boxit, 0) = 0.0;
               if (gcoef_data) (*gcoef_data)(*boxit, 0) = exactFcn(x, y);
            }
            break;
         case 2:
            // min j edge
            y = xlo[1];
            for ( ; boxit != boxitend; ++boxit) {
               i = (*boxit)[0];
               x = xlo[1] + dx[1] * (i - patch_box.lower()[1] + 0.5);
               if (acoef_data) (*acoef_data)(*boxit, 0) = 1.0;
               if (bcoef_data) (*bcoef_data)(*boxit, 0) = 0.0;
               if (gcoef_data) (*gcoef_data)(*boxit, 0) = exactFcn(x, y);
            }
            break;
         case 3:
            // max j edge
            y = xup[1];
            for ( ; boxit != boxitend; ++boxit) {
               i = (*boxit)[0];
               x = xlo[1] + dx[1] * (i - patch_box.lower()[1] + 0.5);
               if (acoef_data) (*acoef_data)(*boxit, 0) = 1.0;
               if (bcoef_data) (*bcoef_data)(*boxit, 0) = 0.0;
               if (gcoef_data) (*gcoef_data)(*boxit, 0) = exactFcn(x, y);
            }
            break;
         default:
            TBOX_ERROR("Invalid location index in\n"
            << "PoissonMultigaussianSolution::setBcCoefs");
      }
   }

   if (d_dim == tbox::Dimension(3)) {
      MDA_Access<double, 3, MDA_OrderColMajor<3> > a_array, b_array, g_array;
      if (acoef_data) a_array = pdat::ArrayDataAccess::access<3, double>(
               *acoef_data);
      if (bcoef_data) b_array = pdat::ArrayDataAccess::access<3, double>(
               *bcoef_data);
      if (gcoef_data) g_array = pdat::ArrayDataAccess::access<3, double>(
               *gcoef_data);
      int i, j, k, ibeg, iend, jbeg, jend, kbeg, kend;
      double x, y, z;
      switch (bdry_box.getLocationIndex()) {
         case 0:
            // min i side
            jbeg = box.lower()[1];
            jend = box.upper()[1];
            kbeg = box.lower()[2];
            kend = box.upper()[2];
            i = box.lower()[0] + 1;
            x = xlo[0];
            for (k = kbeg; k <= kend; ++k) {
               z = xlo[2] + dx[2] * (k - patch_box.lower()[2] + 0.5);
               for (j = jbeg; j <= jend; ++j) {
                  y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
                  if (a_array) a_array(i, j, k) = 1.0;
                  if (b_array) b_array(i, j, k) = 0.0;
                  if (g_array) g_array(i, j, k) = exactFcn(x, y, z);
               }
            }
            break;
         case 1:
            // max i side
            jbeg = box.lower()[1];
            jend = box.upper()[1];
            kbeg = box.lower()[2];
            kend = box.upper()[2];
            i = box.upper()[0];
            x = xup[0];
            for (k = kbeg; k <= kend; ++k) {
               z = xlo[2] + dx[2] * (k - patch_box.lower()[2] + 0.5);
               for (j = jbeg; j <= jend; ++j) {
                  y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
                  if (a_array) a_array(i, j, k) = 1.0;
                  if (b_array) b_array(i, j, k) = 0.0;
                  if (g_array) g_array(i, j, k) = exactFcn(x, y, z);
               }
            }
            break;
         case 2:
            // min j side
            ibeg = box.lower()[0];
            iend = box.upper()[0];
            kbeg = box.lower()[2];
            kend = box.upper()[2];
            j = box.lower()[1] + 1;
            y = xlo[1];
            for (k = kbeg; k <= kend; ++k) {
               z = xlo[2] + dx[2] * (k - patch_box.lower()[2] + 0.5);
               for (i = ibeg; i <= iend; ++i) {
                  x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                  if (a_array) a_array(i, j, k) = 1.0;
                  if (b_array) b_array(i, j, k) = 0.0;
                  if (g_array) g_array(i, j, k) = exactFcn(x, y, z);
               }
            }
            break;
         case 3:
            // max j side
            ibeg = box.lower()[0];
            iend = box.upper()[0];
            kbeg = box.lower()[2];
            kend = box.upper()[2];
            j = box.upper()[1];
            y = xup[1];
            for (k = kbeg; k <= kend; ++k) {
               z = xlo[2] + dx[2] * (k - patch_box.lower()[2] + 0.5);
               for (i = ibeg; i <= iend; ++i) {
                  x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                  if (a_array) a_array(i, j, k) = 1.0;
                  if (b_array) b_array(i, j, k) = 0.0;
                  if (g_array) g_array(i, j, k) = exactFcn(x, y, z);
               }
            }
            break;
         case 4:
            // min k side
            ibeg = box.lower()[0];
            iend = box.upper()[0];
            jbeg = box.lower()[1];
            jend = box.upper()[1];
            k = box.lower()[2] + 1;
            z = xlo[2];
            for (j = jbeg; j <= jend; ++j) {
               y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
               for (i = ibeg; i <= iend; ++i) {
                  x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                  if (a_array) a_array(i, j, k) = 1.0;
                  if (b_array) b_array(i, j, k) = 0.0;
                  if (g_array) g_array(i, j, k) = exactFcn(x, y, z);
               }
            }
            break;
         case 5:
            // max k side
            ibeg = box.lower()[0];
            iend = box.upper()[0];
            jbeg = box.lower()[1];
            jend = box.upper()[1];
            k = box.upper()[2];
            z = xup[2];
            for (j = jbeg; j <= jend; ++j) {
               y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
               for (i = ibeg; i <= iend; ++i) {
                  x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                  if (a_array) a_array(i, j, k) = 1.0;
                  if (b_array) b_array(i, j, k) = 0.0;
                  if (g_array) g_array(i, j, k) = exactFcn(x, y, z);
               }
            }
            break;
         default:
            TBOX_ERROR("Invalid location index in\n"
            << "PoissonMultigaussianSolution::setBcCoefs");
      }
   }
}

/*
 ***********************************************************************
 * This class uses analytical boundary condition, so it can
 * an unlimited number of extensions past the corner of a patch.
 ***********************************************************************
 */
hier::IntVector PoissonMultigaussianSolution::numberOfExtensionsFillable()
const
{
   return hier::IntVector(tbox::Dimension(d_dim), 1000);
}
