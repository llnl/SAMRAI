/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Simple Cartesian grid geometry for an AMR hierarchy.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace geom {


/*
 *************************************************************************
 *
 * Constructor for CartesianPatchGeometry allocates and sets
 * patch coordinate system information.
 *
 *************************************************************************
 */
CartesianPatchGeometry::CartesianPatchGeometry(
   const hier::IntVector& ratio_to_level_zero,
   const TwoDimBool& touches_regular_bdry,
   const hier::BlockId& block_id,
   const double* dx,
   const double* x_lo,
   const double* x_up):
   hier::PatchGeometry(ratio_to_level_zero,
                       touches_regular_bdry,
                       block_id)
{
   TBOX_ASSERT(dx != 0);
   TBOX_ASSERT(x_lo != 0);
   TBOX_ASSERT(x_up != 0);

   const tbox::Dimension& dim(ratio_to_level_zero.getDim());

   for (int id = 0; id < dim.getValue(); ++id) {
      d_dx[id] = dx[id];
      d_x_lo[id] = x_lo[id];
      d_x_up[id] = x_up[id];
   }
}

/*
 *************************************************************************
 *
 * Put a descripton of this geometry into the Blueprint format for
 * coordinates.
 *
 *************************************************************************
 */
void
CartesianPatchGeometry::putBlueprintCoords(
   const std::shared_ptr<tbox::Database>& coords_db,
   const hier::Box& box) const
{
   const tbox::Dimension& dim(box.getDim());

   coords_db->putString("type", "uniform");
   std::shared_ptr<tbox::Database> dims_db(
      coords_db->putDatabase("dims"));
   std::shared_ptr<tbox::Database> origin_db(
      coords_db->putDatabase("origin"));
   std::shared_ptr<tbox::Database> spacing_db(
      coords_db->putDatabase("spacing"));

   hier::IntVector box_size(box.numberCells());

   dims_db->putInteger("i", box_size[0]);
   origin_db->putDouble("x", d_x_lo[0]);
   spacing_db->putDouble("dx", d_dx[0]);

   if (dim.getValue() > 1) {
      dims_db->putInteger("j", box_size[1]);
      origin_db->putDouble("y", d_x_lo[1]);
      spacing_db->putDouble("dy", d_dx[1]);
   }

   if (dim.getValue() > 2) {
      dims_db->putInteger("k", box_size[2]);
      origin_db->putDouble("z", d_x_lo[2]);
      spacing_db->putDouble("dz", d_dx[2]);
   }

}

/*
 *************************************************************************
 *
 * Destructor for CartesianPatchGeometry deallocates dx array.
 *
 *************************************************************************
 */
CartesianPatchGeometry::~CartesianPatchGeometry()
{
}

/*
 *************************************************************************
 *
 * Print CartesianPatchGeometry class data.
 *
 *************************************************************************
 */
void
CartesianPatchGeometry::printClassData(
   std::ostream& os) const
{
   const tbox::Dimension& dim(getRatio().getDim());

   os << "Printing CartesianPatchGeometry data: this = "
      << (CartesianPatchGeometry *)this << std::endl;
   os << "x_lo = ";
   for (int id1 = 0; id1 < dim.getValue(); ++id1) {
      os << d_x_lo[id1] << "   ";
   }
   os << std::endl;
   os << "x_up = ";
   for (int id2 = 0; id2 < dim.getValue(); ++id2) {
      os << d_x_up[id2] << "   ";
   }
   os << std::endl;
   os << "dx = ";
   for (int id3 = 0; id3 < dim.getValue(); ++id3) {
      os << d_dx[id3] << "   ";
   }
   os << std::endl;

   hier::PatchGeometry::printClassData(os);
}

}
}
