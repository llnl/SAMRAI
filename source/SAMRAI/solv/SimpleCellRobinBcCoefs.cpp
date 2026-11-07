/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Level solver for diffusion-like elliptic problems.
 *
 ************************************************************************/

#include "SAMRAI/solv/SimpleCellRobinBcCoefs.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/OuterfaceData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace solv {

/*
 *************************************************************************
 *
 * Construct an unitialized boundary specification.
 *
 *************************************************************************
 */

SimpleCellRobinBcCoefs::SimpleCellRobinBcCoefs(
   const tbox::Dimension& dim,
   const std::string& object_name):
   d_dim(dim),
   d_object_name(object_name),
   d_ln_min(-1),
   d_ln_max(-1),
   d_flux_id(-1),
   d_flag_id(-1),
   d_diffusion_coef_id(-1)
{
   t_set_bc_coefs = tbox::TimerManager::getManager()->
      getTimer("solv::SimpleCellRobinBcCoefs::setBcCoefs()");
}

SimpleCellRobinBcCoefs::~SimpleCellRobinBcCoefs()
{
}

void
SimpleCellRobinBcCoefs::setHierarchy(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int ln_min,
   const int ln_max)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *hierarchy);

   d_hierarchy = hierarchy;
   d_ln_min = ln_min;
   d_ln_max = ln_max;
   d_dirichlet_data.clear();
   d_dirichlet_data_pos.clear();

   if (d_ln_min == -1) {
      d_ln_min = 0;
   }
   if (d_ln_max == -1) {
      d_ln_min = d_hierarchy->getFinestLevelNumber();
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_ln_min < 0 || d_ln_max < 0 || d_ln_min > d_ln_max) {
      TBOX_ERROR(d_object_name
         << ": Bad range of levels in setHierarchy().\n");
   }
#endif
}

void
SimpleCellRobinBcCoefs::setBoundaries(
   const std::string& boundary_type,
   const int fluxes,
   const int flags,
   int* bdry_types)
{

   int k;

   if (boundary_type == "Dirichlet") {
      d_flux_id = -1;
      d_flag_id = -1;
      for (k = 0; k < 2 * d_dim.getValue(); ++k) {
         d_bdry_types[k] = DIRICHLET;
      }
   } else if (boundary_type == "Neumann") {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (fluxes < 0) {
         TBOX_ERROR(
            d_object_name << ": bad flux patch data index ("
                          << fluxes
                          << ") for Neumann boundary condition.\n");
      }
#endif
      for (k = 0; k < 2 * d_dim.getValue(); ++k) {
         d_bdry_types[k] = NEUMANN;
      }
      d_flux_id = fluxes;
      d_flag_id = -1;
   } else if (boundary_type == "Mixed") {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (fluxes < 0) {
         TBOX_ERROR(
            d_object_name << ": bad flux patch data index ("
                          << fluxes
                          << ") for Mixed boundary condition.\n");
      }
      if (flags < 0) {
         TBOX_ERROR(
            d_object_name << ": bad flag patch data index ("
                          << flags
                          << ") for Mixed boundary condition.\n");
      }
#endif
      d_flux_id = fluxes;
      d_flag_id = flags;
      if (bdry_types != 0) {
         for (k = 0; k < 2 * d_dim.getValue(); ++k) {
            d_bdry_types[k] = bdry_types[k];
         }
      } else {
         for (k = 0; k < 2 * d_dim.getValue(); ++k) {
            d_bdry_types[k] = MIXED;
         }
      }
   } else {
      TBOX_ERROR(
         d_object_name << ": Non-existing case of\n"
                       << "boundary_type in SimpleCellRobinBcCoefs::setBoundaries()");
   }

}

/*
 ************************************************************************
 *
 * Set the bc coefficients according to information received in the
 * call to setBoundaries.
 *
 * Do a lot of error checking before hand.
 *
 * For Dirichlet, we need the Dirichlet boundary
 * values stored in the ghost cells of the solution.
 *
 * For Neumann bc, we need the flux data and the
 * diffusion coefficient to determine the required
 * normal gradient of the solution.  However, the
 * diffusion coefficient is assumed to be 1.0 if
 * left unspecified.
 *
 * For mixed bc, we need the flag stating whether
 * Dirichlet or Neumann at any face, in addition to
 * Dirichlet and Neumann data.
 *
 ************************************************************************
 */

void
SimpleCellRobinBcCoefs::setBcCoefs(
   const std::shared_ptr<pdat::ArrayData<double> >& acoef_data,
   const std::shared_ptr<pdat::ArrayData<double> >& bcoef_data,
   const std::shared_ptr<pdat::ArrayData<double> >& gcoef_data,
   const std::shared_ptr<hier::Variable>& variable,
   const hier::Patch& patch,
   const hier::BoundaryBox& bdry_box,
   double fill_time) const
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY2(d_dim, patch, bdry_box);

#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   if (acoef_data) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *acoef_data);
   }
   if (bcoef_data) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *bcoef_data);
   }
   if (gcoef_data) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *gcoef_data);
   }
#endif

   NULL_USE(variable);
   NULL_USE(fill_time);

   t_set_bc_coefs->start();

   const int ln = patch.getPatchLevelNumber();
   const hier::GlobalId& global_id = patch.getGlobalId();
   const int location_index = bdry_box.getLocationIndex();

   std::shared_ptr<hier::PatchData> patch_data;
   std::shared_ptr<pdat::OuterfaceData<double> > flux_data_ptr;
   std::shared_ptr<pdat::SideData<double> > diffcoef_data_ptr;
   std::shared_ptr<pdat::OuterfaceData<int> > flag_data_ptr;

#ifdef DEBUG_CHECK_ASSERTIONS
   if (gcoef_data) {
      if (d_bdry_types[location_index] == DIRICHLET
          || d_bdry_types[location_index] == MIXED) {
         /*
          * For Dirichlet and mixed boundaries, we use cached data
          * to get the Dirichlet value.  Data specific to the
          * d_hierarchy has been cached by cacheDirichletData().
          * This implementation can only set Dirichlet coefficients
          * when the patch is in the correct hierarchy.
          */
         if (!patch.inHierarchy()) {
            TBOX_ERROR(
               d_object_name << ": patch is not in any hierarchy.\n"
                             << "SimpleCellRobinBcCoefs can only set\n"
                             << "boundary coefficients for patches in\n"
                             << "the same hierarchy as cached\n"
                             << "Dirichlet coefficients.");
         }
         std::shared_ptr<hier::PatchLevel> level(
            d_hierarchy->getPatchLevel(ln));
         if (!level->getPatch(global_id)->getBox().isSpatiallyEqual(patch.getBox())) {
            TBOX_ERROR(
               d_object_name << ": patch is not in the hierarchy\n"
                             << "of cached boundary data.\n"
                             << "SimpleCellRobinBcCoefs can only set\n"
                             << "boundary coefficients for patches in\n"
                             << "the same hierarchy as cached\n"
                             << "Dirichlet coefficients.");
         }
      }
      if (d_bdry_types[location_index] == NEUMANN
          || d_bdry_types[location_index] == MIXED) {
         patch_data = patch.getPatchData(d_flux_id);
         if (!patch_data) {
            TBOX_ERROR(d_object_name << ": Flux data (patch data id = "
                                     << d_flux_id << ") does not exist.");
         }
         flux_data_ptr =
            SAMRAI_SHARED_PTR_CAST<pdat::OuterfaceData<double>, hier::PatchData>(patch_data);
         TBOX_ASSERT(flux_data_ptr);
         if (d_diffusion_coef_id != -1) {
            patch_data = patch.getPatchData(d_diffusion_coef_id);
            if (!patch_data) {
               TBOX_ERROR(d_object_name << ": Diffusion coefficient data\n"
                  "(patch data id = " << d_diffusion_coef_id
                                        << ") does not exist.");
            }
            diffcoef_data_ptr =
               SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(patch_data);
            TBOX_ASSERT(diffcoef_data_ptr);
         }
      }
   }
   if (acoef_data) {
      if (d_bdry_types[location_index] == MIXED) {
         patch_data = patch.getPatchData(d_flag_id);
         if (!patch_data) {
            TBOX_ERROR(d_object_name << ": Flags data (patch data id = "
                                     << d_flag_id << ") does not exist.");
         }
         flag_data_ptr = SAMRAI_SHARED_PTR_CAST<pdat::OuterfaceData<int>, hier::PatchData>(
               patch.getPatchData(d_flag_id));
         TBOX_ASSERT(flag_data_ptr);
      }
   }
#endif

   int bn;

   if (d_bdry_types[location_index] == DIRICHLET) {

      if (acoef_data) {
         acoef_data->fill(1.0);
      }
      if (bcoef_data) {
         bcoef_data->fill(0.0);
      }

      if (gcoef_data) {

         std::shared_ptr<geom::CartesianPatchGeometry> pg(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch.getPatchGeometry()));

         TBOX_ASSERT(pg);

         const std::vector<hier::BoundaryBox>& codim1_boxes =
            pg->getCodimensionBoundaries(1);
         /*
          * Search for cached boundary box containing current boundary box.
          */
         for (bn = 0; bn < static_cast<int>(codim1_boxes.size()); ++bn) {
            const hier::BoundaryBox& cdb = codim1_boxes[bn];
            if (bdry_box.getLocationIndex() == cdb.getLocationIndex()
                && bdry_box.getBox().lower() >= cdb.getBox().lower()
                && bdry_box.getBox().upper() <= cdb.getBox().upper()
                ) break;
         }
#ifdef DEBUG_CHECK_ASSERTIONS
         if (bn == static_cast<int>(codim1_boxes.size())) {
            TBOX_ERROR(
               d_object_name << " cannot find cached Dirichlet data.\n"
                             << "This is most likely caused by not calling\n"
                             << "SimpleCellRobinBcCoefs::cacheDirichletData()\n"
                             << "after the hierarchy changed.\n");
         }
#endif
         hier::BoxId box_id(global_id);

         std::map<hier::BoxId, int> foo = d_dirichlet_data_pos[ln];
         int position = foo[box_id] + bn;
         gcoef_data->copy(*d_dirichlet_data[position],
            d_dirichlet_data[position]->getBox(),
            hier::IntVector::getZero(d_dim));

      }
   } else if (d_bdry_types[location_index] == NEUMANN) {

      if (acoef_data) {
         acoef_data->fill(0.0);
      }
      if (bcoef_data) {
         bcoef_data->fill(1.0);
      }

      if (gcoef_data) {
         flux_data_ptr =
            SAMRAI_SHARED_PTR_CAST<pdat::OuterfaceData<double>, hier::PatchData>(
               patch.getPatchData(d_flux_id));
         TBOX_ASSERT(flux_data_ptr);
         pdat::OuterfaceData<double>& flux_data(*flux_data_ptr);
         const int axis = location_index / 2;
         const int face = location_index % 2;
         pdat::ArrayData<double>& g = *gcoef_data;
         pdat::ArrayDataIterator ai(g.getBox(), true);
         pdat::ArrayDataIterator aiend(g.getBox(), false);
         hier::Index offset_to_inside(d_dim, 0);
         if (face != 0) offset_to_inside(axis) = -1;
         if (d_diffusion_coef_id == -1) {
            for ( ; ai != aiend; ++ai) {
               pdat::FaceIndex fi(*ai + offset_to_inside, axis, face);
               g(*ai, 0) = flux_data(fi, face) / d_diffusion_coef_constant;
            }
         } else {
            diffcoef_data_ptr =
               SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                  patch.getPatchData(d_diffusion_coef_id));
            TBOX_ASSERT(diffcoef_data_ptr);
            const pdat::ArrayData<double>& diffcoef_array_data =
               diffcoef_data_ptr->getArrayData(axis);
            for ( ; ai != aiend; ++ai) {
               pdat::FaceIndex fi(*ai + offset_to_inside, axis, face);
               g(*ai, 0) = flux_data(fi, face) / diffcoef_array_data(*ai, 0);
            }
         }
      }

   } else if (d_bdry_types[location_index] == MIXED) {

      const int axis = location_index / 2;
      const int face = location_index % 2;
      flag_data_ptr = SAMRAI_SHARED_PTR_CAST<pdat::OuterfaceData<int>, hier::PatchData>(
            patch.getPatchData(d_flag_id));
      TBOX_ASSERT(flag_data_ptr);
      pdat::OuterfaceData<int>& flag_data(*flag_data_ptr);
      hier::Index offset_to_inside(d_dim, 0);
      if (face != 0) offset_to_inside(axis) = -1;

      if (acoef_data) {
         pdat::ArrayData<double>& a = *acoef_data;
         pdat::ArrayDataIterator ai(a.getBox(), true);
         pdat::ArrayDataIterator aiend(a.getBox(), false);
         for ( ; ai != aiend; ++ai) {
            pdat::FaceIndex fi(*ai + offset_to_inside, axis, face);
            if (flag_data(fi, face) == 0) {
               a(*ai, 0) = 1.0;
            } else {
               a(*ai, 0) = 0.0;
            }
         }
      }

      if (bcoef_data) {
         pdat::ArrayData<double>& b = *bcoef_data;
         pdat::ArrayDataIterator bi(b.getBox(), true);
         pdat::ArrayDataIterator biend(b.getBox(), false);
         for ( ; bi != biend; ++bi) {
            pdat::FaceIndex fi(*bi + offset_to_inside, axis, face);
            if (flag_data(fi, face) == 0) {
               b(*bi, 0) = 0.0;
            } else {
               b(*bi, 0) = 1.0;
            }
         }
      }

      if (gcoef_data) {
         std::shared_ptr<geom::CartesianPatchGeometry> pg(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch.getPatchGeometry()));

         TBOX_ASSERT(pg);

         const std::vector<hier::BoundaryBox>& codim1_boxes =
            pg->getCodimensionBoundaries(1);
         /*
          * Search for cached boundary box containing current boundary box.
          */
         for (bn = 0; bn < static_cast<int>(codim1_boxes.size()); ++bn) {
            const hier::BoundaryBox& cdb = codim1_boxes[bn];
            if (bdry_box.getLocationIndex() == cdb.getLocationIndex()
                && bdry_box.getBox().lower() >= cdb.getBox().lower()
                && bdry_box.getBox().upper() <= cdb.getBox().upper()
                ) break;
         }
#ifdef DEBUG_CHECK_ASSERTIONS
         if (bn == static_cast<int>(codim1_boxes.size())) {
            TBOX_ERROR(
               d_object_name << " cannot find cached Dirichlet data.\n"
                             << "This is most likely caused by not calling\n"
                             << "SimpleCellRobinBcCoefs::cacheDirichletData() after the\n"
                             << "hierarchy changed.\n");
         }
#endif
         hier::BoxId box_id(global_id);

         std::map<hier::BoxId, int> foo = d_dirichlet_data_pos[ln];
         int position = foo[box_id] + bn;

         const pdat::ArrayData<double>& dirichlet_array_data =
            *d_dirichlet_data[position];
         diffcoef_data_ptr =
            SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
               patch.getPatchData(d_diffusion_coef_id));
         TBOX_ASSERT(diffcoef_data_ptr);
         pdat::ArrayData<double>& g = *gcoef_data;
         pdat::OuterfaceData<double>& flux_data(*flux_data_ptr);
         pdat::ArrayDataIterator ai(g.getBox(), true);
         pdat::ArrayDataIterator aiend(g.getBox(), false);
         for ( ; ai != aiend; ++ai) {
            pdat::FaceIndex fi(*ai + offset_to_inside, axis, face);
            if (flag_data(fi, face) == 0) {
               g(*ai, 0) = dirichlet_array_data(*ai, 0);
            } else {
               pdat::FaceIndex fi2(*ai + offset_to_inside, axis, face);
               if (d_diffusion_coef_id == -1) {
                  g(*ai, 0) = flux_data(fi2, face) / d_diffusion_coef_constant;
               } else {
                  g(*ai, 0) = flux_data(fi2, face)
                     / diffcoef_data_ptr->getArrayData(axis) (*ai, 0);
               }
            }
         }
      }

   }

   t_set_bc_coefs->stop();
}

/*
 ***********************************************************************
 * This class cannot set coefficients for boundary boxes that extend
 * past the patch in the direction parallel to the boundary,
 * because it relies on data, such as pdat::OutersideData<TYPE>,
 * that does not extend.
 ***********************************************************************
 */
hier::IntVector
SimpleCellRobinBcCoefs::numberOfExtensionsFillable() const
{
   return hier::IntVector::getZero(d_dim);
}

/*
 ************************************************************************
 *
 * Copy and save cell-centered Dirichlet data in ghost cells.
 * For each boundary box in the hierarchy, we create a pdat::ArrayData<TYPE>
 * object indexed on the side indices corresponding to boundary boxes.
 * The ghost-cell-centered data is shifted to the side indices and
 * saved in the pdat::ArrayData<TYPE> objects.
 *
 * First, loop through the hierarchy to compute how many
 * pdat::ArrayData<TYPE> objects we need and the position of each one.
 *
 * Second, allocate the pdat::ArrayData<TYPE> objects.
 *
 * Third, loop through the hierarchy again to allocate the data in each
 * pdat::ArrayData<TYPE> object and cache the ghost data.
 *
 * The position of the appropriate boundary box bn of patch pn
 * of level ln is d_dirichlet_data_pos[ln][pn]+bn
 *
 ************************************************************************
 */
void
SimpleCellRobinBcCoefs::cacheDirichletData(
   int dirichlet_data_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!d_hierarchy) {
      TBOX_ERROR(
         d_object_name << ": hierarchy has not been set.\n"
                       << "use setHierarchy() to set the hierarchy before\n"
                       << "caching boundary ghost cell data.\n");
   }
#endif

   tbox::ResourceAllocator allocator =
      tbox::AllocatorDatabase::getDatabase()->getDefaultAllocator();

   d_dirichlet_data.clear();
   d_dirichlet_data_pos.clear();
   int ln, bn, position, n_reqd_boxes = 0;
   d_dirichlet_data_pos.resize(d_ln_max + 1);
   for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
      hier::PatchLevel& level = (hier::PatchLevel &)
         * d_hierarchy->getPatchLevel(ln);
      hier::PatchLevel::iterator pi(level.begin());
      for ( ; pi != level.end(); ++pi) {
         hier::Patch& patch = **pi;
         const hier::GlobalId& global_id = patch.getGlobalId();
         hier::BoxId box_id(global_id);
         std::shared_ptr<geom::CartesianPatchGeometry> pg(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch.getPatchGeometry()));

         TBOX_ASSERT(pg);

         const std::vector<hier::BoundaryBox>& codim1_boxes =
            pg->getCodimensionBoundaries(1);
         d_dirichlet_data_pos[ln][box_id] = n_reqd_boxes;
         n_reqd_boxes += static_cast<int>(codim1_boxes.size());
      }
   }
   d_dirichlet_data.resize(n_reqd_boxes);
   for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
      hier::PatchLevel& level = (hier::PatchLevel &)
         * d_hierarchy->getPatchLevel(ln);
      hier::PatchLevel::iterator pi(level.begin());
      for ( ; pi != level.end(); ++pi) {
         hier::Patch& patch = **pi;
         const hier::GlobalId& global_id = patch.getGlobalId();
         hier::BoxId box_id(global_id);
         std::shared_ptr<pdat::CellData<double> > cell_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(dirichlet_data_id)));
         std::shared_ptr<geom::CartesianPatchGeometry> pg(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                       hier::PatchGeometry>(patch.getPatchGeometry()));

         TBOX_ASSERT(cell_data);
         TBOX_ASSERT(pg);

         const std::vector<hier::BoundaryBox>& codim1_boxes =
            pg->getCodimensionBoundaries(1);
         for (bn = 0; bn < static_cast<int>(codim1_boxes.size()); ++bn) {
            const hier::BoundaryBox& bdry_box = codim1_boxes[bn];
            position = d_dirichlet_data_pos[ln][box_id] + bn;
            hier::Box databox = makeSideBoundaryBox(bdry_box);
            d_dirichlet_data[position].reset(
               new pdat::ArrayData<double>(databox, 1, allocator));
            pdat::ArrayData<double>& array_data = *d_dirichlet_data[position];
            hier::IntVector shift_amount(d_dim, 0);
            const int location_index = bdry_box.getLocationIndex();
            if (location_index % 2 == 0) shift_amount[location_index / 2] = 1;
            array_data.copy(cell_data->getArrayData(),
               databox,
               shift_amount);
         }
      }
   }
}

/*
 ************************************************************************
 *
 * Reverse action of cacheDirichletData by copying cached data back
 * into the ghost cells.
 *
 * The cached data is not deallocated.
 *
 ************************************************************************
 */
void
SimpleCellRobinBcCoefs::restoreDirichletData(
   int dirichlet_data_id)
{
   if (d_dirichlet_data_pos.empty()) {
      TBOX_ERROR(d_object_name << ".restoreDirichletData(): Dirichlet\n"
                               << "data has not been set.\n");
   }
   int ln, bn, position;
   for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
      hier::PatchLevel& level = (hier::PatchLevel &)
         * d_hierarchy->getPatchLevel(ln);
      hier::PatchLevel::iterator pi(level.begin());
      for ( ; pi != level.end(); ++pi) {
         hier::Patch& patch = **pi;
         const hier::GlobalId& global_id = patch.getGlobalId();
         hier::BoxId box_id(global_id);
         std::shared_ptr<pdat::CellData<double> > cell_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(dirichlet_data_id)));
         std::shared_ptr<geom::CartesianPatchGeometry> pg(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch.getPatchGeometry()));

         TBOX_ASSERT(cell_data);
         TBOX_ASSERT(pg);

         const std::vector<hier::BoundaryBox>& codim1_boxes =
            pg->getCodimensionBoundaries(1);
         for (bn = 0; bn < static_cast<int>(codim1_boxes.size()); ++bn) {
            const hier::BoundaryBox& bdry_box = codim1_boxes[bn];
            position = d_dirichlet_data_pos[ln][box_id] + bn;
            pdat::ArrayData<double>& array_data = *d_dirichlet_data[position];
            hier::IntVector shift_amount(d_dim, 0);
            const int location_index = bdry_box.getLocationIndex();
            hier::Box dst_box = array_data.getBox();
            if (location_index % 2 == 0) {
               shift_amount[location_index / 2] = -1;
               dst_box.shift(static_cast<tbox::Dimension::dir_t>(location_index / 2), -1);
            }
            cell_data->getArrayData().copy(array_data,
               dst_box,
               shift_amount);
         }
      }
   }
}

/*
 ************************************************************************
 * Make surface box on boundary using standard boundary box
 ************************************************************************
 */

hier::Box
SimpleCellRobinBcCoefs::makeSideBoundaryBox(
   const hier::BoundaryBox& boundary_box) const
{
   if (boundary_box.getBoundaryType() != 1) {
      TBOX_ERROR(
         d_object_name
         << ": CartesianRobinBcHelper::makeSideBoundaryBox called with\n"
         << "improper boundary box\n"
         << "for " << d_object_name);
   }
   hier::Box face_indices = boundary_box.getBox();
   int location_index = boundary_box.getLocationIndex();
   if (location_index % 2 == 0) {
      /*
       * On the min index side, the face indices are one higher
       * than the boundary cell indices, in the direction normal
       * to the boundary.
       */
      face_indices.shift(static_cast<tbox::Dimension::dir_t>(location_index / 2), 1);
   }
   return face_indices;
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
