/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Robin boundary condition support on cartesian grids.
 *
 ************************************************************************/
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/math/ArrayDataBasicOps.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include IOMANIP_HEADER_FILE

#include "SAMRAI/solv/GhostCellRobinBcCoefs.h"

namespace SAMRAI {
namespace solv {

/*
 ************************************************************************
 * Constructor
 ************************************************************************
 */

GhostCellRobinBcCoefs::GhostCellRobinBcCoefs(
   const tbox::Dimension& dim,
   std::string object_name):
   d_object_name(object_name),
   d_dim(dim),
   d_ghost_data_id(-1),
   d_extensions_fillable(dim)
{

   t_set_bc_coefs = tbox::TimerManager::getManager()->
      getTimer("solv::GhostCellRobinBcCoefs::setBcCoefs()");
}

/*
 ************************************************************************
 * Destructor
 ************************************************************************
 */

GhostCellRobinBcCoefs::~GhostCellRobinBcCoefs()
{
}

/*
 ************************************************************************
 * Set the index of the data providing ghost cell values
 ************************************************************************
 */

void
GhostCellRobinBcCoefs::setGhostDataId(
   int ghost_data_id,
   hier::IntVector extensions_fillable)
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, extensions_fillable);

   d_ghost_data_id = ghost_data_id;
   d_extensions_fillable = extensions_fillable;
   /*
    * Check for correctness of data index.
    * Unfortunately, the ghost width is not provided by the
    * variable database, so we cannot check that also.
    */
   if (d_ghost_data_id != -1) {
      hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();
      std::shared_ptr<hier::Variable> variable_ptr;
      vdb->mapIndexToVariable(ghost_data_id, variable_ptr);
      if (!variable_ptr) {
         TBOX_ERROR(d_object_name << ": hier::Index " << ghost_data_id
                                  << " does not correspond to any variable.");
      }
      std::shared_ptr<pdat::CellVariable<double> > cell_variable_ptr(
         SAMRAI_SHARED_PTR_CAST<pdat::CellVariable<double>, hier::Variable>(variable_ptr));
      TBOX_ASSERT(cell_variable_ptr);
   }
}

/*
 ************************************************************************
 * Set the bc coefficients reflect the value at the ghost cell centers.
 * The factor 1.0/(1+0.5*h) appears in a and g.  This factor comes
 * from a linear approximation of the data through the patch boundary,
 * going through the centers of the first interior and ghost cells
 * and having the specified values there.
 ************************************************************************
 */

void
GhostCellRobinBcCoefs::setBcCoefs(
   const std::shared_ptr<pdat::ArrayData<double> >& acoef_data,
   const std::shared_ptr<pdat::ArrayData<double> >& bcoef_data,
   const std::shared_ptr<pdat::ArrayData<double> >& gcoef_data,
   const std::shared_ptr<hier::Variable>& variable,
   const hier::Patch& patch,
   const hier::BoundaryBox& bdry_box,
   double fill_time) const
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(d_dim, *variable, patch, bdry_box);

   NULL_USE(variable);
   NULL_USE(fill_time);

   t_set_bc_coefs->start();

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));

   TBOX_ASSERT(patch_geom);

   const int norm_dir = bdry_box.getLocationIndex() / 2;
   const double* dx = patch_geom->getDx();
   const double h = dx[norm_dir];

   /*
    * Set acoef_data to 1.0/(1+0.5*h) uniformly.  This value
    * corresponds to the fact that the solution is fixed at
    * the ghost cell centers.  bcoef_data is 1-acoef_data.
    */
   if (acoef_data) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *acoef_data);

      acoef_data->fill(1.0 / (1 + 0.5 * h));
   }
   if (bcoef_data) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *bcoef_data);

      bcoef_data->fill(0.5 * h / (1 + 0.5 * h));
   }

   if (gcoef_data) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *gcoef_data);

      if (d_ghost_data_id == -1) {
         TBOX_ERROR(d_object_name << ": Coefficient g requested without\n"
                                  << "having valid ghost data id.\n");
      }

      /*
       * Fill in gcoef_data with data from d_ghost_data_id.
       * The data is first looked for in a pdat::OutersideData<TYPE> object
       * and a pdat::CellData<TYPE> object in that order.  Data from the
       * first place with allocated storage is used.
       */
      std::shared_ptr<hier::PatchData> patch_data(
         patch.getPatchData(d_ghost_data_id));
      if (!patch_data) {
         TBOX_ERROR(d_object_name << ": hier::Patch data for index "
                                  << d_ghost_data_id << " does not exist.");
      }
      std::shared_ptr<pdat::CellData<double> > cell_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(patch_data));

      TBOX_ASSERT(cell_data);

      const int location_index = bdry_box.getLocationIndex();
      const hier::IntVector& gw = cell_data->getGhostCellWidth();
      if (gw[norm_dir] < 1) {
         TBOX_ERROR(
            d_object_name << ": hier::Patch data for index "
                          << d_ghost_data_id
                          << " has zero ghost width.");
      }
      const pdat::ArrayData<double>& cell_array_data =
         cell_data->getArrayData();
      hier::IntVector shift_amount(d_dim, 0);
      if (location_index % 2 == 0) shift_amount[location_index / 2] = 1;
      gcoef_data->copy(cell_array_data,
         makeSideBoundaryBox(bdry_box),
         shift_amount);
      math::ArrayDataBasicOps<double> aops;
      /*
       * To convert from the value at the ghost cell center
       * to the coefficient g, we must scale the data by
       * 1/(1+h/2), according to our linear approximation
       * of the data at the patch boundary.
       */
      aops.scale(*gcoef_data, 1.0 / (1 + 0.5 * h), *gcoef_data,
         makeSideBoundaryBox(bdry_box));

   }

   t_set_bc_coefs->stop();
}

/*
 ***********************************************************************
 * This class can only set coeficients for boundary boxes that extend
 * no more than what the data it uses provides.
 ***********************************************************************
 */
hier::IntVector
GhostCellRobinBcCoefs::numberOfExtensionsFillable() const
{
   return d_extensions_fillable;
}

/*
 ************************************************************************
 * Make surface box on boundary using standard boundary box
 ************************************************************************
 */

hier::Box
GhostCellRobinBcCoefs::makeSideBoundaryBox(
   const hier::BoundaryBox& boundary_box) const
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, boundary_box);

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
