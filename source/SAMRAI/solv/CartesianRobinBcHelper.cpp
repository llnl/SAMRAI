/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Robin boundary condition support on cartesian grids.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include IOMANIP_HEADER_FILE

#include "SAMRAI/solv/CartesianRobinBcHelper.h"

#include <vector>

extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

void SAMRAI_F77_FUNC(settype1cells2d, SETTYPE1CELLS2D) (
   double* data,
   const int& difirst, const int& dilast,
   const int& djfirst, const int& djlast,
   const double* a, const double* b, const double* g,
   const int& ifirst, const int& ilast,
   const int& jfirst, const int& jlast,
   const int& ibeg, const int& iend,
   const int& jbeg, const int& jend,
   const int& face, const int& ghos, const int& inte, const int& location,
   const double& h, const int& zerog);
void SAMRAI_F77_FUNC(settype2cells2d, SETTYPE2CELLS2D) (
   double* data,
   const int& difirst, const int& dilast,
   const int& djfirst, const int& djlast,
   const int* lower, const int* upper, const int& location);
void SAMRAI_F77_FUNC(settype1cells3d, SETTYPE1CELLS3D) (
   double* data,
   const int& difirst, const int& dilast,
   const int& djfirst, const int& djlast,
   const int& dkfirst, const int& dklast,
   const double* a, const double* b, const double* g,
   const int& ifirst, const int& ilast,
   const int& jfirst, const int& jlast,
   const int& kfirst, const int& klast,
   const int& ibeg, const int& iend,
   const int& jbeg, const int& jend,
   const int& kbeg, const int& kend,
   const int& face, const int& ghos, const int& inte, const int& location,
   const double& h, const int& zerog);
void SAMRAI_F77_FUNC(settype2cells3d, SETTYPE2CELLS3D) (
   double* data,
   const int& difirst, const int& dilast,
   const int& djfirst, const int& djlast,
   const int& dkfirst, const int& dklast,
   const int* lower, const int* upper, const int& location);
void SAMRAI_F77_FUNC(settype3cells3d, SETTYPE3CELLS3D) (
   double* data,
   const int& difirst, const int& dilast,
   const int& djfirst, const int& djlast,
   const int& dkfirst, const int& dklast,
   const int* lower, const int* upper, const int& location);
}

namespace SAMRAI {
namespace solv {

/*
 ************************************************************************
 * Constructor
 ************************************************************************
 */

CartesianRobinBcHelper::CartesianRobinBcHelper(
   const tbox::Dimension& dim,
   std::string object_name,
   RobinBcCoefStrategy* coef_strategy):
   xfer::RefinePatchStrategy(),
   d_object_name(object_name),
   d_dim(dim),
   d_coef_strategy(0),
   d_target_data_id(-1),
   d_homogeneous_bc(false)
{

   NULL_USE(coef_strategy);

   d_allocator =
      tbox::AllocatorDatabase::getDatabase()->getDefaultAllocator();

   t_set_boundary_values_in_cells = tbox::TimerManager::getManager()->
      getTimer("solv::CartesianRobinBcHelper::setBoundaryValuesInCells()");
   t_use_set_bc_coefs = tbox::TimerManager::getManager()->
      getTimer(
         "solv::CartesianRobinBcHelper::setBoundaryValuesInCells()_setBcCoefs");
}

/*
 ************************************************************************
 * Destructor
 ************************************************************************
 */

CartesianRobinBcHelper::~CartesianRobinBcHelper() {
}

/*
 ************************************************************************
 *     Set physical boundary conditions in cells.
 ************************************************************************
 */

void
CartesianRobinBcHelper::setBoundaryValuesInCells(
   hier::Patch& patch,
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill,
   int target_data_id,
   bool homogeneous_bc) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(patch, ghost_width_to_fill);

   NULL_USE(fill_time);

   const tbox::Dimension& dim(patch.getDim());

   t_set_boundary_values_in_cells->start();

#ifdef DEBUG_CHECK_ASSERTIONS
   if (!d_coef_strategy) {
      TBOX_ERROR(d_object_name << ": coefficient strategy is not set.\n"
                               << "Use setCoefImplementation() to set it.\n");
   }
#endif

   if (patch.getDim() == tbox::Dimension(1)) {
      TBOX_ERROR(d_object_name << ": dim = 1 not supported");
   }
   math::PatchCellDataOpsReal<double> cops;

   /*
    * Get info on the data.
    */
   hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();
   std::shared_ptr<hier::Variable> variable_ptr;
   vdb->mapIndexToVariable(target_data_id, variable_ptr);
   if (!variable_ptr) {
      TBOX_ERROR(d_object_name << ": No variable for index "
                               << target_data_id);
   }
   std::shared_ptr<pdat::CellVariable<double> > cell_variable_ptr(
      SAMRAI_SHARED_PTR_CAST<pdat::CellVariable<double>, hier::Variable>(variable_ptr));
   TBOX_ASSERT(cell_variable_ptr);

   /*
    * Get the data.
    */
   std::shared_ptr<hier::PatchData> data_ptr(
      patch.getPatchData(target_data_id));
   if (!data_ptr) {
      TBOX_ERROR(d_object_name << ": No data for index " << target_data_id);
   }
   std::shared_ptr<pdat::CellData<double> > cell_data_ptr(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(data_ptr));
   TBOX_ASSERT(cell_data_ptr);
   pdat::CellData<double>& data = *cell_data_ptr;

   const hier::IntVector& ghost_cells =
      cell_data_ptr->getGhostCellWidth();
   hier::IntVector gcw_to_fill = hier::IntVector::min(ghost_cells,
         ghost_width_to_fill);
   if (!(gcw_to_fill == hier::IntVector::getZero(dim))) {
      /*
       * Given a and g in a*u + (1-a)*un = g,
       * where un is the derivative in the outward normal direction,
       * and ui (the value of u in the first interior cell),
       * we compute the value on the outer face
       * uf = ...
       * and the normal derivative on the outer face
       * un = ...
       * and the uo (the value in the first ghost cell)
       * uo = ...
       */
      const hier::Box& patch_box(patch.getBox());

      /*
       * These definitions can go in the next block.
       * They are kept her for debugging.
       */
      std::shared_ptr<geom::CartesianPatchGeometry> pg(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));

      TBOX_ASSERT(pg);

      const std::vector<hier::BoundaryBox>& codim1_boxes =
         pg->getCodimensionBoundaries(1);

      const int n_codim1_boxes = static_cast<int>(codim1_boxes.size());

      const hier::Box& ghost_box = data.getGhostBox();
      const double* h = pg->getDx();
      const int num_coefs(homogeneous_bc ? 1 : 2);
      const int zerog = num_coefs == 1;

      for (int n = 0; n < n_codim1_boxes; ++n) {

         const int location_index = codim1_boxes[n].getLocationIndex();
         const int normal_dir = location_index / 2;
         if (!gcw_to_fill(normal_dir)) {
            // Zero ghost width to fill for this boundary box.
            continue;
         }
         hier::IntVector extension_amount(d_dim, 1);
         extension_amount(normal_dir) = 0;
         const hier::BoundaryBox boundary_box =
            d_coef_strategy->numberOfExtensionsFillable() >= extension_amount ?
            trimBoundaryBox(codim1_boxes[n], ghost_box) :
            trimBoundaryBox(codim1_boxes[n], patch_box);
         const hier::Index& lower = boundary_box.getBox().lower();
         const hier::Index& upper = boundary_box.getBox().upper();
         const hier::Box coefbox = makeFaceBoundaryBox(boundary_box);
         std::shared_ptr<pdat::ArrayData<double> > acoef_data(
            std::make_shared<pdat::ArrayData<double> >(coefbox, 1, d_allocator));
         std::shared_ptr<pdat::ArrayData<double> > bcoef_data(
            std::make_shared<pdat::ArrayData<double> >(coefbox, 1, d_allocator));
         std::shared_ptr<pdat::ArrayData<double> > gcoef_data(
            homogeneous_bc ? 0 :
            new pdat::ArrayData<double>(coefbox, 1, d_allocator));
         t_use_set_bc_coefs->start();
         d_coef_strategy->setBcCoefs(acoef_data,
            bcoef_data,
            gcoef_data,
            variable_ptr,
            patch,
            boundary_box,
            fill_time);
         t_use_set_bc_coefs->stop();

         int igho, ifac, iint, ibeg, iend;
         double dx;
         int jgho, jfac, jint, jbeg, jend;
         double dy;
         int kgho, kfac, kint, kbeg, kend;
         double dz;

         if (d_dim == tbox::Dimension(2)) {
            switch (location_index) {
               case 0:
                  // min i edge
                  dx = h[0];
                  igho = lower[0]; // Lower and upper are the same.
                  ifac = igho + 1;
                  iint = igho + 1;
                  jbeg = lower[1];
                  jend = upper[1];
                  SAMRAI_F77_FUNC(settype1cells2d, SETTYPE1CELLS2D) (data.getPointer(0),
                  ghost_box.lower()[0], ghost_box.upper()[0],
                  ghost_box.lower()[1], ghost_box.upper()[1],
                  acoef_data->getPointer(),
                  bcoef_data->getPointer(),
                  gcoef_data ? gcoef_data->getPointer() : 0,
                  coefbox.lower()[0], coefbox.upper()[0],
                  coefbox.lower()[1], coefbox.upper()[1],
                  igho, igho, jbeg, jend,
                  ifac, igho, iint, location_index, dx, zerog
                  );
                  break;
               case 1:
                  // max i edge
                  dx = h[0];
                  igho = lower[0]; // Lower and upper are the same.
                  ifac = igho;
                  iint = igho - 1;
                  jbeg = lower[1];
                  jend = upper[1];
                  SAMRAI_F77_FUNC(settype1cells2d, SETTYPE1CELLS2D) (data.getPointer(0),
                  ghost_box.lower()[0], ghost_box.upper()[0],
                  ghost_box.lower()[1], ghost_box.upper()[1],
                  acoef_data->getPointer(),
                  bcoef_data->getPointer(),
                  gcoef_data ? gcoef_data->getPointer() : 0,
                  coefbox.lower()[0], coefbox.upper()[0],
                  coefbox.lower()[1], coefbox.upper()[1],
                  igho, igho, jbeg, jend,
                  ifac, igho, iint, location_index, dx, zerog
                  );
                  break;
               case 2:
                  // min j edge
                  dy = h[1];
                  jgho = lower[1]; // Lower and upper are the same.
                  jfac = jgho + 1;
                  jint = jgho + 1;
                  ibeg = lower[0];
                  iend = upper[0];
                  SAMRAI_F77_FUNC(settype1cells2d, SETTYPE1CELLS2D) (data.getPointer(0),
                  ghost_box.lower()[0], ghost_box.upper()[0],
                  ghost_box.lower()[1], ghost_box.upper()[1],
                  acoef_data->getPointer(),
                  bcoef_data->getPointer(),
                  gcoef_data ? gcoef_data->getPointer() : 0,
                  coefbox.lower()[0], coefbox.upper()[0],
                  coefbox.lower()[1], coefbox.upper()[1],
                  ibeg, iend, jgho, jgho,
                  jfac, jgho, jint, location_index, dy, zerog
                  );
                  break;
               case 3:
                  // max j edge
                  dy = h[1];
                  jgho = lower[1]; // Lower and upper are the same.
                  jfac = jgho;
                  jint = jgho - 1;
                  ibeg = lower[0];
                  iend = upper[0];
                  SAMRAI_F77_FUNC(settype1cells2d, SETTYPE1CELLS2D) (data.getPointer(0),
                  ghost_box.lower()[0], ghost_box.upper()[0],
                  ghost_box.lower()[1], ghost_box.upper()[1],
                  acoef_data->getPointer(),
                  bcoef_data->getPointer(),
                  gcoef_data ? gcoef_data->getPointer() : 0,
                  coefbox.lower()[0], coefbox.upper()[0],
                  coefbox.lower()[1], coefbox.upper()[1],
                  ibeg, iend, jgho, jgho,
                  jfac, jgho, jint, location_index, dy, zerog
                  );
                  break;
               default:
                  TBOX_ERROR(d_object_name << ": Invalid location index ("
                                           << location_index << ") in\n"
                                           << "setBoundaryValuesInCells");
            }
         } else if (d_dim == tbox::Dimension(3)) {
            switch (location_index) {
               case 0:
                  // min i face
                  dx = h[0];
                  igho = lower[0]; // Lower and upper are the same.
                  ifac = igho + 1;
                  iint = igho + 1;
                  jbeg = lower[1];
                  jend = upper[1];
                  kbeg = lower[2];
                  kend = upper[2];
                  SAMRAI_F77_FUNC(settype1cells3d, SETTYPE1CELLS3D) (data.getPointer(0),
                  ghost_box.lower()[0], ghost_box.upper()[0],
                  ghost_box.lower()[1], ghost_box.upper()[1],
                  ghost_box.lower()[2], ghost_box.upper()[2],
                  acoef_data->getPointer(),
                  bcoef_data->getPointer(),
                  gcoef_data ? gcoef_data->getPointer() : 0,
                  coefbox.lower()[0], coefbox.upper()[0],
                  coefbox.lower()[1], coefbox.upper()[1],
                  coefbox.lower()[2], coefbox.upper()[2],
                  igho, igho, jbeg, jend, kbeg, kend,
                  ifac, igho, iint, location_index, dx, zerog
                  );
                  break;
               case 1:
                  // max i face
                  dx = h[0];
                  igho = lower[0]; // Lower and upper are the same.
                  ifac = igho;
                  iint = igho - 1;
                  jbeg = lower[1];
                  jend = upper[1];
                  kbeg = lower[2];
                  kend = upper[2];
                  SAMRAI_F77_FUNC(settype1cells3d, SETTYPE1CELLS3D) (data.getPointer(0),
                  ghost_box.lower()[0], ghost_box.upper()[0],
                  ghost_box.lower()[1], ghost_box.upper()[1],
                  ghost_box.lower()[2], ghost_box.upper()[2],
                  acoef_data->getPointer(),
                  bcoef_data->getPointer(),
                  gcoef_data ? gcoef_data->getPointer() : 0,
                  coefbox.lower()[0], coefbox.upper()[0],
                  coefbox.lower()[1], coefbox.upper()[1],
                  coefbox.lower()[2], coefbox.upper()[2],
                  igho, igho, jbeg, jend, kbeg, kend,
                  ifac, igho, iint, location_index, dx, zerog
                  );
                  break;
               case 2:
                  // min j face
                  dy = h[1];
                  jgho = lower[1]; // Lower and upper are the same.
                  jfac = jgho + 1;
                  jint = jgho + 1;
                  ibeg = lower[0];
                  iend = upper[0];
                  kbeg = lower[2];
                  kend = upper[2];
                  SAMRAI_F77_FUNC(settype1cells3d, SETTYPE1CELLS3D) (data.getPointer(0),
                  ghost_box.lower()[0], ghost_box.upper()[0],
                  ghost_box.lower()[1], ghost_box.upper()[1],
                  ghost_box.lower()[2], ghost_box.upper()[2],
                  acoef_data->getPointer(),
                  bcoef_data->getPointer(),
                  gcoef_data ? gcoef_data->getPointer() : 0,
                  coefbox.lower()[0], coefbox.upper()[0],
                  coefbox.lower()[1], coefbox.upper()[1],
                  coefbox.lower()[2], coefbox.upper()[2],
                  ibeg, iend, jgho, jgho, kbeg, kend,
                  jfac, jgho, jint, location_index, dy, zerog
                  );
                  break;
               case 3:
                  // max j face
                  dy = h[1];
                  jgho = lower[1]; // Lower and upper are the same.
                  jfac = jgho;
                  jint = jgho - 1;
                  ibeg = lower[0];
                  iend = upper[0];
                  kbeg = lower[2];
                  kend = upper[2];
                  SAMRAI_F77_FUNC(settype1cells3d, SETTYPE1CELLS3D) (data.getPointer(0),
                  ghost_box.lower()[0], ghost_box.upper()[0],
                  ghost_box.lower()[1], ghost_box.upper()[1],
                  ghost_box.lower()[2], ghost_box.upper()[2],
                  acoef_data->getPointer(),
                  bcoef_data->getPointer(),
                  gcoef_data ? gcoef_data->getPointer() : 0,
                  coefbox.lower()[0], coefbox.upper()[0],
                  coefbox.lower()[1], coefbox.upper()[1],
                  coefbox.lower()[2], coefbox.upper()[2],
                  ibeg, iend, jgho, jgho, kbeg, kend,
                  jfac, jgho, jint, location_index, dy, zerog
                  );
                  break;
               case 4:
                  // min k face
                  dz = h[2];
                  kgho = lower[2]; // Lower and upper are the same.
                  kfac = kgho + 1;
                  kint = kgho + 1;
                  ibeg = lower[0];
                  iend = upper[0];
                  jbeg = lower[1];
                  jend = upper[1];
                  SAMRAI_F77_FUNC(settype1cells3d, SETTYPE1CELLS3D) (data.getPointer(0),
                  ghost_box.lower()[0], ghost_box.upper()[0],
                  ghost_box.lower()[1], ghost_box.upper()[1],
                  ghost_box.lower()[2], ghost_box.upper()[2],
                  acoef_data->getPointer(),
                  bcoef_data->getPointer(),
                  gcoef_data ? gcoef_data->getPointer() : 0,
                  coefbox.lower()[0], coefbox.upper()[0],
                  coefbox.lower()[1], coefbox.upper()[1],
                  coefbox.lower()[2], coefbox.upper()[2],
                  ibeg, iend, jbeg, jend, kgho, kgho,
                  kfac, kgho, kint, location_index, dz, zerog
                  );
                  break;
               case 5:
                  // max k face
                  dz = h[2];
                  kgho = lower[2]; // Lower and upper are the same.
                  kfac = kgho;
                  kint = kgho - 1;
                  ibeg = lower[0];
                  iend = upper[0];
                  jbeg = lower[1];
                  jend = upper[1];
                  SAMRAI_F77_FUNC(settype1cells3d, SETTYPE1CELLS3D) (data.getPointer(0),
                  ghost_box.lower()[0], ghost_box.upper()[0],
                  ghost_box.lower()[1], ghost_box.upper()[1],
                  ghost_box.lower()[2], ghost_box.upper()[2],
                  acoef_data->getPointer(),
                  bcoef_data->getPointer(),
                  gcoef_data ? gcoef_data->getPointer() : 0,
                  coefbox.lower()[0], coefbox.upper()[0],
                  coefbox.lower()[1], coefbox.upper()[1],
                  coefbox.lower()[2], coefbox.upper()[2],
                  ibeg, iend, jbeg, jend, kgho, kgho,
                  kfac, kgho, kint, location_index, dz, zerog
                  );
                  break;
               default:
                  TBOX_ERROR(d_object_name << ": Invalid location index ("
                                           << location_index << ") in\n"
                                           << "setBoundaryValuesInCells");
            }
         }
      }

      /*
       * Now that the surface boundary boxes have been set,
       * the rest of this function set the lower-dimensional
       * boundary boxes.  Users may not need to have these
       * set, but refiners may.
       */

      if (d_dim == tbox::Dimension(2)) {
         /*
          * The node boundary conditions are set from a linear interpolation
          * through the nearest interior cell and the two nearest edge values.
          * This data may be used by refinement operators to do interpolation.
          */

         const std::vector<hier::BoundaryBox>& node_boxes =
            pg->getNodeBoundaries();
         const int n_node_boxes = static_cast<int>(node_boxes.size());
         for (int n = 0; n < n_node_boxes; ++n) {
            const hier::BoundaryBox& bb = node_boxes[n];
            TBOX_ASSERT(bb.getBoundaryType() == 2);        // Must be a node boundary.
            const hier::Box& bb_box = bb.getBox();
            const hier::Index& lower = bb_box.lower();
            const hier::Index& upper = bb_box.upper();
            const int location_index = bb.getLocationIndex();
            SAMRAI_F77_FUNC(settype2cells2d, SETTYPE2CELLS2D) (data.getPointer(0),
               ghost_box.lower()[0], ghost_box.upper()[0],
               ghost_box.lower()[1], ghost_box.upper()[1],
               &lower[0], &upper[0], location_index);
         }
      } else if (d_dim == tbox::Dimension(3)) {
         /*
          * The edge boundary conditions are set from a linear interpolation
          * through the nearest interior cell and the two nearest side values.
          * This data may be used by refinement operators to do interpolation.
          */
         const std::vector<hier::BoundaryBox>& edge_boxes =
            pg->getEdgeBoundaries();
         const int n_edge_boxes = static_cast<int>(edge_boxes.size());
         for (int n = 0; n < n_edge_boxes; ++n) {
            const int location_index = edge_boxes[n].getLocationIndex();
            const int edge_dir = 2 - (location_index / 4);
            hier::IntVector extension_amount(d_dim, 0);
            extension_amount(edge_dir) = 1;
            const hier::BoundaryBox boundary_box =
               d_coef_strategy->numberOfExtensionsFillable() >= extension_amount ?
               trimBoundaryBox(edge_boxes[n], ghost_box) :
               trimBoundaryBox(edge_boxes[n], patch_box);
            TBOX_ASSERT(boundary_box.getBoundaryType() == 2);
            const hier::Index& lower = boundary_box.getBox().lower();
            const hier::Index& upper = boundary_box.getBox().upper();
            SAMRAI_F77_FUNC(settype2cells3d, SETTYPE2CELLS3D) (data.getPointer(0),
               ghost_box.lower()[0], ghost_box.upper()[0],
               ghost_box.lower()[1], ghost_box.upper()[1],
               ghost_box.lower()[2], ghost_box.upper()[2],
               &lower[0], &upper[0], location_index);
         }

         /*
          * The node boundary conditions are set from a linear interpolation
          * through the nearest interior cell and the three nearest edge values.
          * This data may be used by refinement operators to do interpolation.
          */
         const std::vector<hier::BoundaryBox>& node_boxes =
            pg->getNodeBoundaries();
         const int n_node_boxes = static_cast<int>(node_boxes.size());
         for (int n = 0; n < n_node_boxes; ++n) {
            const hier::BoundaryBox& bb = node_boxes[n];
            TBOX_ASSERT(bb.getBoundaryType() == 3); // Must be an node boundary.
            const hier::Box& bb_box = bb.getBox();
            const hier::Index& lower = bb_box.lower();
            const hier::Index& upper = bb_box.upper();
            TBOX_ASSERT(lower == upper);
            const int location_index = bb.getLocationIndex();
            SAMRAI_F77_FUNC(settype3cells3d, SETTYPE3CELLS3D) (data.getPointer(0),
               ghost_box.lower()[0], ghost_box.upper()[0],
               ghost_box.lower()[1], ghost_box.upper()[1],
               ghost_box.lower()[2], ghost_box.upper()[2],
               &lower[0], &upper[0], location_index);
         }
      } else {
         TBOX_ERROR("CartesianRobinBcHelper::setBoundaryValuesInCells error ..."
            << "\n not implemented for dim>3" << std::endl);
      }
   }

   t_set_boundary_values_in_cells->stop();
}

/*
 ************************************************************************
 * Set physical boundary conditions in cells, for all patches in a
 * given level.
 ************************************************************************
 */

void
CartesianRobinBcHelper::setBoundaryValuesInCells(
   hier::PatchLevel& level,
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill,
   int target_data_id,
   bool homogeneous_bc) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(level, ghost_width_to_fill);

   for (hier::PatchLevel::iterator p(level.begin());
        p != level.end(); ++p) {
      const std::shared_ptr<hier::Patch>& patch = *p;
      setBoundaryValuesInCells(*patch,
         fill_time,
         ghost_width_to_fill,
         target_data_id,
         homogeneous_bc);
   }
}

/*
 ************************************************************************
 * Set physical boundary conditions at nodes.
 ************************************************************************
 */

void
CartesianRobinBcHelper::setBoundaryValuesAtNodes(
   hier::Patch& patch,
   const double fill_time,
   int target_data_id,
   bool homogeneous_bc) const
{
   NULL_USE(patch);
   NULL_USE(fill_time);
   NULL_USE(target_data_id);
   NULL_USE(homogeneous_bc);

   TBOX_ERROR(
      d_object_name << ": Using incomplete implementation"
                    << "CartesianRobinBcHelper::setBoundaryValuesAtNodes"
                    << "is not implemented because there is not a need for it (yet)"
                    << std::endl);
}

/*
 ***********************************************************************
 *
 *  Virtual functions or xfer::RefinePatchStrategy.
 *
 ***********************************************************************
 */

void
CartesianRobinBcHelper::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(patch, ghost_width_to_fill);

   setBoundaryValuesInCells(patch,
      fill_time,
      ghost_width_to_fill,
      d_target_data_id,
      d_homogeneous_bc);
}

hier::IntVector
CartesianRobinBcHelper::getRefineOpStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getZero(dim);
}

void
CartesianRobinBcHelper::preprocessRefineBoxes(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::BoxContainer& fine_boxes,
   const hier::IntVector& ratio)
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_boxes);
   NULL_USE(ratio);
}
void
CartesianRobinBcHelper::preprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}
void
CartesianRobinBcHelper::postprocessRefineBoxes(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::BoxContainer& fine_box,
   const hier::IntVector& ratio)
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}
void
CartesianRobinBcHelper::postprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_boxes,
   const hier::IntVector& ratio)
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_boxes);
   NULL_USE(ratio);
}

/*
 ************************************************************************
 * Trim a boundary box so it does not stick out past the corners of a
 * given box.  This removes the extension parallel to the boundary,
 * past the corner of the limit box.
 ************************************************************************
 */

hier::BoundaryBox
CartesianRobinBcHelper::trimBoundaryBox(
   const hier::BoundaryBox& boundary_box,
   const hier::Box& limit_box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(boundary_box, limit_box);

   if (boundary_box.getBoundaryType() == d_dim.getValue()) {
      // This is a node boundary box and cannot be trimmed anymore.
      return boundary_box;
   }

   const hier::Box& bbox = boundary_box.getBox();
   const hier::Index& plo = limit_box.lower();
   const hier::Index& pup = limit_box.upper();
   const hier::Index& blo = bbox.lower();
   const hier::Index& bup = bbox.upper();
   hier::Index newlo(d_dim), newup(d_dim);
   int key_direction;
   int d;
   switch (boundary_box.getBoundaryType()) {
      case 2:
         key_direction = 2 - (boundary_box.getLocationIndex() / 4);
         for (d = 0; d < d_dim.getValue(); ++d) {
            if (d == key_direction) {
               newlo(d) = tbox::MathUtilities<int>::Max(blo(d), plo(d));
               newup(d) = tbox::MathUtilities<int>::Min(bup(d), pup(d));
            } else {
               newlo(d) = blo(d);
               newup(d) = bup(d);
            }
         }
         break;
      case 1:
         key_direction = boundary_box.getLocationIndex() / 2;
         /*
          * Loop through directions.
          * Preserve box size in direction normal to boundary.
          * Trim box size in direction transverse to boundary.
          */
         for (d = 0; d < d_dim.getValue(); ++d) {
            if (d == key_direction) {
               newlo(d) = blo(d);
               newup(d) = bup(d);
            } else {
               // Min side.  Use max between boundary and patch boxes.
               newlo(d) = tbox::MathUtilities<int>::Max(blo(d), plo(d));
               // Max side.  Use min between boundary and patch boxes.
               newup(d) = tbox::MathUtilities<int>::Min(bup(d), pup(d));
            }
         }
         break;
   }
   const hier::Box newbox(newlo, newup, boundary_box.getBox().getBlockId());
   const hier::BoundaryBox newbbox(newbox,
                                   boundary_box.getBoundaryType(),
                                   boundary_box.getLocationIndex());
   return newbbox;
}

/*
 ************************************************************************
 * Make surface box on boundary using standard boundary box
 ************************************************************************
 */

hier::Box
CartesianRobinBcHelper::makeFaceBoundaryBox(
   const hier::BoundaryBox& boundary_box) const
{
   if (boundary_box.getBoundaryType() != 1) {
      TBOX_ERROR(d_object_name << ": makeFaceBoundaryBox called with\n"
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

/*
 ************************************************************************
 * Make node box on boundary using standard boundary box
 ************************************************************************
 */

hier::Box
CartesianRobinBcHelper::makeNodeBoundaryBox(
   const hier::BoundaryBox& boundary_box) const
{
   if (boundary_box.getBoundaryType() != 1) {
      TBOX_ERROR(d_object_name << ": makeNodeBoundaryBox called with\n"
                               << "improper boundary box\n"
                               << "for " << d_object_name);
   }
   hier::Box node_indices = boundary_box.getBox();
   int location_index = boundary_box.getLocationIndex();
   if (location_index % 2 == 0) {
      /*
       * On the min index side, the node indices are one higher
       * than the boundary cell indices, in the direction normal
       * to the boundary.
       */
      node_indices.shift(static_cast<tbox::Dimension::dir_t>(location_index / 2), 1);
   }
   /*
    * The node indices range one higher than the cell indices,
    * in the directions parallel to the boundary.
    */
   hier::IntVector parallel_growth(d_dim, 1);
   parallel_growth(location_index / 2) = 0;
   node_indices.growUpper(parallel_growth);
   return node_indices;
}

}
}
