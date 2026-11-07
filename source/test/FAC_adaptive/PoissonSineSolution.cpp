/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   PoissonSineSolution class implementation
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/pdat/ArrayDataAccess.h"
#include "patchFcns.h"
#include "PoissonSineSolution.h"
#include STL_SSTREAM_HEADER_FILE

#include "SAMRAI/geom/CartesianPatchGeometry.h"

using namespace SAMRAI;

PoissonSineSolution::PoissonSineSolution(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_linear_coef(0.0),
   d_exact(dim)
{
   int i;
   for (i = 0; i < 2 * d_dim.getValue(); ++i) {
      d_neumann_location[i] = false;
   }
}

PoissonSineSolution::PoissonSineSolution(
   const std::string& object_name,
   const tbox::Dimension& dim,
   tbox::Database& database,
   std::ostream* out_stream,
   std::ostream* log_stream):
   d_dim(dim),
   d_linear_coef(0.0),
   d_exact(dim)
{
   NULL_USE(object_name);
   NULL_USE(out_stream);
   NULL_USE(log_stream);

   int i;
   for (i = 0; i < 2 * d_dim.getValue(); ++i) {
      d_neumann_location[i] = false;
   }
   setFromDatabase(database);
}

PoissonSineSolution::~PoissonSineSolution()
{
}

void PoissonSineSolution::setFromDatabase(
   tbox::Database& database)
{
   std::string istr = database.getStringWithDefault("SinusoidFcnControl", "{}");
   std::istringstream ist(istr);
   ist >> d_exact;
   if (database.isBool("neumann_locations")) {
      std::vector<bool> neumann_locations =
         database.getBoolVector("neumann_locations");
      if (static_cast<int>(neumann_locations.size()) > 2 * d_dim.getValue()) {
         TBOX_ERROR(
            "'neumann_locations' should have at most " << 2 * d_dim.getValue()
                                                       << " entries in " << d_dim << "D.\n");
      }
      int i;
      for (i = 0; i < static_cast<int>(neumann_locations.size()); ++i) {
         d_neumann_location[i] = neumann_locations[i];
      }
   }
   d_linear_coef = database.getDoubleWithDefault("linear_coef",
         d_linear_coef);
}

void PoissonSineSolution::setNeumannLocation(
   int location_index,
   bool flag)
{
   TBOX_ASSERT(location_index < 2 * d_dim.getValue());
   d_neumann_location[location_index] = flag;
}

void PoissonSineSolution::setPoissonSpecifications(
   solv::PoissonSpecifications& sps,
   int C_patch_data_id,
   int D_patch_data_id) const
{
   NULL_USE(C_patch_data_id);
   NULL_USE(D_patch_data_id);

   sps.setDConstant(1.0);
   sps.setCConstant(d_linear_coef);
}

void PoissonSineSolution::setGridData(
   hier::Patch& patch,
   pdat::CellData<double>& exact_data,
   pdat::CellData<double>& source_data)
{
   /* Set source function and exact solution. */
   /*
    * For the forcing function
    * -( 1 + (nx^2 + ny^2) pi^2 ) * sin(pi*(nx*x+px))*sin(pi*(ny*y+py))
    * the exact solution is sin(pi*(nx*x+px))*sin(pi*(ny*y+py))
    */
   setCellDataToSinusoid(exact_data,
      patch,
      d_exact);
   setCellDataToSinusoid(source_data,
      patch,
      d_exact);
   double npi[SAMRAI::MAX_DIM_VAL],
          ppi[SAMRAI::MAX_DIM_VAL];
   d_exact.getWaveNumbers(npi);
   d_exact.getPhaseAngles(ppi);
   double source_scale = 0.0;
   if (d_dim == tbox::Dimension(2)) {
      source_scale = d_linear_coef
         - ((npi[0] * npi[0] + npi[1] * npi[1]) * M_PI * M_PI);
   }
   if (d_dim == tbox::Dimension(3)) {
      source_scale = d_linear_coef
         - ((npi[0] * npi[0] + npi[1] * npi[1] + npi[2] * npi[2]) * M_PI * M_PI);
   }
   scaleArrayData(source_data.getArrayData(), source_scale);
}       // End patch loop.

std::ostream& operator << (
   std::ostream& os,
   const PoissonSineSolution& r) {
   os << r.d_exact << "\n";
   return os;
}

void PoissonSineSolution::setBcCoefs(
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

   if (bdry_box.getBoundaryType() != 1) {
      // Must be a face boundary.
      TBOX_ERROR("Bad boundary type in\n"
         << "PoissonSineSolution::setBcCoefs \n");
   }

   const int location_index = bdry_box.getLocationIndex();

   // a is either 0 (Neumann) or 1 (Dirichlet).
   if (acoef_data)
      acoef_data->fill(d_neumann_location[location_index] ? 0.0 : 1.0, 0);
   if (bcoef_data)
      bcoef_data->fill(d_neumann_location[location_index] ? 1.0 : 0.0, 0);

   /*
    * Get geometry information needed to compute coordinates
    * of side centers.
    */
   hier::Box patch_box(patch.getBox());
   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const double* xlo = patch_geom->getXLower();
   const double* xup = patch_geom->getXUpper();
   const double* dx = patch_geom->getDx();
   const hier::Box& box = bdry_box.getBox();
   hier::Index lower = box.lower();
   hier::Index upper = box.upper();

   if (gcoef_data) {
      if (d_dim == tbox::Dimension(2)) {
         hier::Box::iterator boxit(gcoef_data->getBox().begin());
         hier::Box::iterator boxitend(gcoef_data->getBox().end());
         int i, j;
         double x, y;
         switch (location_index) {
            case 0:
               // min i edge
               x = xlo[0];
               if (d_neumann_location[location_index]) {
                  SinusoidFcn slope = d_exact.differentiate(1, 0);
                  if (gcoef_data) {
                     for ( ; boxit != boxitend; ++boxit) {
                        j = (*boxit)[1];
                        y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
                        (*gcoef_data)(*boxit, 0) = -slope(x, y);
                     }
                  }
               } else {
                  if (gcoef_data) {
                     for ( ; boxit != boxitend; ++boxit) {
                        j = (*boxit)[1];
                        y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
                        (*gcoef_data)(*boxit, 0) = d_exact(x, y);
                     }
                  }
               }
               break;
            case 1:
               // max i edge
               x = xup[0];
               if (d_neumann_location[location_index]) {
                  SinusoidFcn slope = d_exact.differentiate(1, 0);
                  if (gcoef_data) {
                     for ( ; boxit != boxitend; ++boxit) {
                        j = (*boxit)[1];
                        y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
                        (*gcoef_data)(*boxit, 0) = slope(x, y);
                     }
                  }
               } else {
                  if (gcoef_data) {
                     for ( ; boxit != boxitend; ++boxit) {
                        j = (*boxit)[1];
                        y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
                        (*gcoef_data)(*boxit, 0) = d_exact(x, y);
                     }
                  }
               }
               break;
            case 2:
               // min j edge
               y = xlo[1];
               if (d_neumann_location[location_index]) {
                  SinusoidFcn slope = d_exact.differentiate(0, 1);
                  if (gcoef_data) {
                     for ( ; boxit != boxitend; ++boxit) {
                        i = (*boxit)[0];
                        x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                        if (gcoef_data) {
                           (*gcoef_data)(*boxit, 0) = -slope(x, y);
                        }
                     }
                  }
               } else {
                  if (gcoef_data) {
                     for ( ; boxit != boxitend; ++boxit) {
                        i = (*boxit)[0];
                        x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                        (*gcoef_data)(*boxit, 0) = d_exact(x, y);
                     }
                  }
               }
               break;
            case 3:
               // max j edge
               y = xup[1];
               if (d_neumann_location[location_index]) {
                  SinusoidFcn slope = d_exact.differentiate(0, 1);
                  if (gcoef_data) {
                     for ( ; boxit != boxitend; ++boxit) {
                        i = (*boxit)[0];
                        x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                        (*gcoef_data)(*boxit, 0) = slope(x, y);
                     }
                  }
               } else {
                  if (gcoef_data) {
                     for ( ; boxit != boxitend; ++boxit) {
                        i = (*boxit)[0];
                        x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                        (*gcoef_data)(*boxit, 0) = d_exact(x, y);
                     }
                  }
               }
               break;
            default:
               TBOX_ERROR("Invalid location index in\n"
               << "PoissonSineSolution>::setBcCoefs");
         }
      }

      if (d_dim == tbox::Dimension(3)) {
         MDA_Access<double, 3, MDA_OrderColMajor<3> > g_array;
         if (gcoef_data) {
            g_array = pdat::ArrayDataAccess::access<3, double>(*gcoef_data);
         }
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
               if (d_neumann_location[location_index]) {
                  SinusoidFcn slope = d_exact.differentiate(1, 0, 0);
                  for (k = kbeg; k <= kend; ++k) {
                     z = xlo[2] + dx[2] * (k - patch_box.lower()[2] + 0.5);
                     for (j = jbeg; j <= jend; ++j) {
                        y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
                        if (g_array) g_array(i, j, k) = -slope(x, y, z);
                     }
                  }
               } else {
                  for (k = kbeg; k <= kend; ++k) {
                     z = xlo[2] + dx[2] * (k - patch_box.lower()[2] + 0.5);
                     for (j = jbeg; j <= jend; ++j) {
                        y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
                        if (g_array) g_array(i, j, k) = d_exact(x, y, z);
                     }
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
               if (d_neumann_location[location_index]) {
                  SinusoidFcn slope = d_exact.differentiate(1, 0, 0);
                  for (k = kbeg; k <= kend; ++k) {
                     z = xlo[2] + dx[2] * (k - patch_box.lower()[2] + 0.5);
                     for (j = jbeg; j <= jend; ++j) {
                        y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
                        if (g_array) g_array(i, j, k) = slope(x, y, z);
                     }
                  }
               } else {
                  for (k = kbeg; k <= kend; ++k) {
                     z = xlo[2] + dx[2] * (k - patch_box.lower()[2] + 0.5);
                     for (j = jbeg; j <= jend; ++j) {
                        y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
                        if (g_array) g_array(i, j, k) = d_exact(x, y, z);
                     }
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
               if (d_neumann_location[location_index]) {
                  SinusoidFcn slope = d_exact.differentiate(0, 1, 0);
                  for (k = kbeg; k <= kend; ++k) {
                     z = xlo[2] + dx[2] * (k - patch_box.lower()[2] + 0.5);
                     for (i = ibeg; i <= iend; ++i) {
                        x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                        if (g_array) g_array(i, j, k) = -slope(x, y, z);
                     }
                  }
               } else {
                  for (k = kbeg; k <= kend; ++k) {
                     z = xlo[2] + dx[2] * (k - patch_box.lower()[2] + 0.5);
                     for (i = ibeg; i <= iend; ++i) {
                        x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                        if (g_array) g_array(i, j, k) = d_exact(x, y, z);
                     }
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
               if (d_neumann_location[location_index]) {
                  SinusoidFcn slope = d_exact.differentiate(0, 1, 0);
                  for (k = kbeg; k <= kend; ++k) {
                     z = xlo[2] + dx[2] * (k - patch_box.lower()[2] + 0.5);
                     for (i = ibeg; i <= iend; ++i) {
                        x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                        if (g_array) g_array(i, j, k) = slope(x, y, z);
                     }
                  }
               } else {
                  for (k = kbeg; k <= kend; ++k) {
                     z = xlo[2] + dx[2] * (k - patch_box.lower()[2] + 0.5);
                     for (i = ibeg; i <= iend; ++i) {
                        x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                        if (g_array) g_array(i, j, k) = d_exact(x, y, z);
                     }
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
               if (d_neumann_location[location_index]) {
                  SinusoidFcn slope = d_exact.differentiate(0, 0, 1);
                  for (j = jbeg; j <= jend; ++j) {
                     y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
                     for (i = ibeg; i <= iend; ++i) {
                        x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                        if (g_array) g_array(i, j, k) = -slope(x, y, z);
                     }
                  }
               } else {
                  for (j = jbeg; j <= jend; ++j) {
                     y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
                     for (i = ibeg; i <= iend; ++i) {
                        x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                        if (g_array) g_array(i, j, k) = d_exact(x, y, z);
                     }
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
               if (d_neumann_location[location_index]) {
                  SinusoidFcn slope = d_exact.differentiate(0, 0, 1);
                  for (j = jbeg; j <= jend; ++j) {
                     y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
                     for (i = ibeg; i <= iend; ++i) {
                        x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                        if (g_array) g_array(i, j, k) = slope(x, y, z);
                     }
                  }
               } else {
                  for (j = jbeg; j <= jend; ++j) {
                     y = xlo[1] + dx[1] * (j - patch_box.lower()[1] + 0.5);
                     for (i = ibeg; i <= iend; ++i) {
                        x = xlo[0] + dx[0] * (i - patch_box.lower()[0] + 0.5);
                        if (g_array) g_array(i, j, k) = d_exact(x, y, z);
                     }
                  }
               }
               break;
            default:
               TBOX_ERROR("Invalid location index in\n"
               << "PoissonSineSolution::setBcCoefs");
         }
      }
   }
}

/*
 ***********************************************************************
 * This class uses analytical boundary condition, so it can
 * an unlimited number of extensions past the corner of a patch.
 ***********************************************************************
 */
hier::IntVector PoissonSineSolution::numberOfExtensionsFillable() const
{
   return hier::IntVector(d_dim, 1000);
}
