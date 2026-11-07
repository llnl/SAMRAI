/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   set geometry for multiblock domain
 *
 ************************************************************************/

#include "MblkGeometry.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <cmath>
#include <vector>

#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)

#define POLY3(i, j, k, imin, jmin, kmin, nx, nxny) \
   ((i - imin) + (j - jmin) * (nx) + (k - kmin) * (nxny))

#define MBLK_GEOM_NUM_CHARS 128

/*
 *************************************************************************
 *
 * This class creates mapped multi-block grid geometrys
 *
 *************************************************************************
 */
MblkGeometry::MblkGeometry(
   const std::string& object_name,
   const tbox::Dimension& dim,
   std::shared_ptr<tbox::Database> input_db,
   const size_t nblocks):
   d_dim(dim)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(input_db);

   d_object_name = object_name;
   //tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   d_nblocks = nblocks;

   getFromInput(input_db);

}

/*
 *************************************************************************
 *
 * Empty destructor.
 *
 *************************************************************************
 */
MblkGeometry::~MblkGeometry()
{
}

/*
 *************************************************************************
 *
 * Return the geometry (CARTESIAN, WEDGE, TRILINEAR or SPHERICAL_SHELL)
 *
 *************************************************************************
 */
std::string MblkGeometry::getGeometryType()
{
   return d_geom_problem;
}

/*
 *************************************************************************
 *
 * Return the user-specified refine boxes for a particular block/level
 * number.  If no boxes exist, it returns false.  If they do exist, it
 * returns true and sets the refine boxes argument.
 *
 *************************************************************************
 */
bool MblkGeometry::getRefineBoxes(
   hier::BoxContainer& refine_boxes,
   const hier::BlockId::block_t block_number,
   const int level_number)
{
   bool boxes_exist = false;
   if (block_number < d_refine_boxes.size()) {
      if (level_number < static_cast<int>(d_refine_boxes[level_number].size())) {
         boxes_exist = true;
         refine_boxes = d_refine_boxes[block_number][level_number];
      }
   }
   return boxes_exist;
}

/*
 *************************************************************************
 *
 * Read data members from input.  All values set from restart can be
 * overridden by values in the input database.
 *
 *************************************************************************
 */
void MblkGeometry::getFromInput(
   std::shared_ptr<tbox::Database> input_db)
{
   TBOX_ASSERT(input_db);

   std::shared_ptr<tbox::Database> db(
      input_db->getDatabase("MblkGeometry"));

   d_geom_problem = db->getString("problem_type");

   bool found = false;
   int i;
   char block_name[128];
   double temp_domain[SAMRAI::MAX_DIM_VAL];

   //
   // Cartesian geometry
   //
   if (d_geom_problem == "CARTESIAN") {

      std::shared_ptr<tbox::Database> cart_db(
         db->getDatabase("CartesianGeometry"));

      d_cart_xlo.resize(d_nblocks);
      d_cart_xhi.resize(d_nblocks);

      for (hier::BlockId::block_t nb = 0; nb < d_nblocks; ++nb) {

         // xlo
         snprintf(block_name, MBLK_GEOM_NUM_CHARS, "domain_xlo_%d", nb);
         if (!cart_db->keyExists(block_name)) {
            TBOX_ERROR(d_object_name << ":  "
                                     << "Key data `" << block_name
                                     << "' domain_xlo for block " << nb
                                     << " not found in input." << std::endl);
         }
         d_cart_xlo[nb].resize(d_dim.getValue());
         cart_db->getDoubleArray(block_name, temp_domain, d_dim.getValue());
         for (i = 0; i < d_dim.getValue(); ++i) {
            d_cart_xlo[nb][i] = temp_domain[i];
         }

         // xhi
         snprintf(block_name, MBLK_GEOM_NUM_CHARS, "domain_xhi_%d", nb);
         if (!cart_db->keyExists(block_name)) {
            TBOX_ERROR(d_object_name << ":  "
                                     << "Key data `" << block_name
                                     << "' domain_xhi for block " << nb
                                     << " not found in input." << std::endl);
         }
         d_cart_xhi[nb].resize(d_dim.getValue());
         cart_db->getDoubleArray(block_name, temp_domain, d_dim.getValue());
         for (i = 0; i < d_dim.getValue(); ++i) {
            d_cart_xhi[nb][i] = temp_domain[i];
         }

      }
      found = true;
   }

   //
   // Wedge geometry
   //
   if (d_geom_problem == "WEDGE") {

      std::shared_ptr<tbox::Database> wedge_db(
         db->getDatabase("WedgeGeometry"));

      d_wedge_rmin.resize(d_nblocks);
      d_wedge_rmax.resize(d_nblocks);

      for (hier::BlockId::block_t nb = 0; nb < d_nblocks; ++nb) {

         // rmin
         snprintf(block_name, MBLK_GEOM_NUM_CHARS, "rmin_%d", nb);
         if (!wedge_db->keyExists(block_name)) {
            TBOX_ERROR(d_object_name << ":  "
                                     << "Key data `" << block_name
                                     << "' rmin for block " << nb
                                     << " not found in input." << std::endl);
         }

         d_wedge_rmin[nb] = wedge_db->getDouble(block_name);

         // rmax
         snprintf(block_name, MBLK_GEOM_NUM_CHARS, "rmax_%d", nb);
         if (!wedge_db->keyExists(block_name)) {
            TBOX_ERROR(d_object_name << ":  "
                                     << "Key data `" << block_name
                                     << "' rmax for block " << nb
                                     << " not found in input." << std::endl);
         }

         d_wedge_rmax[nb] = wedge_db->getDouble(block_name);

         // theta min/max
         d_wedge_thmin = wedge_db->getDouble("thmin");
         d_wedge_thmax = wedge_db->getDouble("thmax");

         // Z min/max
         d_wedge_zmin = wedge_db->getDouble("zmin");
         d_wedge_zmax = wedge_db->getDouble("zmax");

      }
      found = true;
   }

   //
   // Trilinear geometry, just read in the base blocks
   //
   if (d_geom_problem == "TRILINEAR") {

      std::shared_ptr<tbox::Database> tri_db(
         db->getDatabase("TrilinearGeometry"));

      d_tri_mesh_filename = tri_db->getString("mesh_filename");

      FILE* fid = fopen(d_tri_mesh_filename.c_str(), "rb");
      int rtest = 0;

      // --------- number of blocks
      rtest = static_cast<int>(fread(&d_tri_nblocks, sizeof(int), 1, fid));
      TBOX_ASSERT(rtest == 1);
      //fscanf( fid, "%d\n", &d_tri_nblocks );

      d_tri_nxp = new int[d_tri_nblocks];
      d_tri_nyp = new int[d_tri_nblocks];
      d_tri_nzp = new int[d_tri_nblocks];

      d_tri_nbr = new hier::BlockId::block_t *[d_tri_nblocks];
      d_tri_x = new double *[d_tri_nblocks];
      d_tri_y = new double *[d_tri_nblocks];
      d_tri_z = new double *[d_tri_nblocks];

      d_tri_node_size = new int[d_tri_nblocks];

      for (int ib = 0; ib < d_tri_nblocks; ++ib) {

         // --------- the size of each block
         rtest = static_cast<int>(fread(&d_tri_nxp[ib], sizeof(int), 1, fid));
         TBOX_ASSERT(rtest == 1);
         rtest = static_cast<int>(fread(&d_tri_nyp[ib], sizeof(int), 1, fid));
         TBOX_ASSERT(rtest == 1);
         rtest = static_cast<int>(fread(&d_tri_nzp[ib], sizeof(int), 1, fid));
         TBOX_ASSERT(rtest == 1);
         //fscanf( fid, "%d %d %d\n",
         //      &d_tri_nxp[ib],
         //      &d_tri_nyp[ib],
         //      &d_tri_nzp[ib] );

         int nsize = d_tri_nxp[ib] * d_tri_nyp[ib] * d_tri_nzp[ib];
         d_tri_node_size[ib] = nsize;
         d_tri_nbr[ib] = new hier::BlockId::block_t[6];
         d_tri_x[ib] = new double[nsize];
         d_tri_y[ib] = new double[nsize];
         d_tri_z[ib] = new double[nsize];

         // --------- the neighbors of each block
         rtest = static_cast<int>(fread(&d_tri_nbr[ib][0], sizeof(int), 6, fid));
         TBOX_ASSERT(rtest == 6);
         //fscanf( fid, "%d %d %d %d %d %d\n",
         //      &d_tri_nbr[ib][0],
         //      &d_tri_nbr[ib][1],
         //      &d_tri_nbr[ib][2],
         //      &d_tri_nbr[ib][3],
         //      &d_tri_nbr[ib][4],
         //      &d_tri_nbr[ib][5] );

         // --------- the mesh positions
         rtest = static_cast<int>(fread(&d_tri_x[ib][0], sizeof(double), nsize, fid));
         TBOX_ASSERT(rtest == nsize);
         rtest = static_cast<int>(fread(&d_tri_y[ib][0], sizeof(double), nsize, fid));
         TBOX_ASSERT(rtest == nsize);
         rtest = static_cast<int>(fread(&d_tri_z[ib][0], sizeof(double), nsize, fid));
         TBOX_ASSERT(rtest == nsize);

         //for ( int ii = 0 ; ii < d_tri_node_size[ib]; ++ii ) {
         //  fscanf( fid, "%20.12e\n", &d_tri_x[ib][ii] );
         //}

         //for ( int ii = 0 ; ii < d_tri_node_size[ib]; ++ii ) {
         //   fscanf( fid, "%20.12e\n", &d_tri_y[ib][ii] );
         //}

         //for ( int ii = 0 ; ii < d_tri_node_size[ib]; ++ii ) {
         //   fscanf( fid, "%20.12e\n", &d_tri_z[ib][ii] );
         //}

      }
      NULL_USE(rtest);

      fclose(fid);
      found = true;
   }

   //
   // Spherical shell (works only in 3d)
   //
   if (d_geom_problem == "SPHERICAL_SHELL") {

      if (d_dim < tbox::Dimension(3)) {
         TBOX_ERROR(d_object_name << ": The " << d_geom_problem
                                  << "only works in 3D." << std::endl);
      }

      std::shared_ptr<tbox::Database> sshell_db(
         db->getDatabase("ShellGeometry"));

      d_sshell_rmin = sshell_db->getDouble("rmin");
      d_sshell_rmax = sshell_db->getDouble("rmax");
      d_sshell_type = sshell_db->getString("shell_type");
      // types are: SOLID, OCTANT

      if (d_sshell_type == "SOLID") {
         double lpi = 3.14159265358979325873406851315;

         d_sangle_degrees = sshell_db->getDouble("solid_angle_degrees");
         d_sangle_thmin =
            -(lpi / 180.0) * tbox::MathUtilities<double>::Abs(d_sangle_degrees);
         // we only need the minimum of the angle coverage in radians
      }

      if (d_sshell_type == "OCTANT") {
         d_tag_velocity = 0.;
         if (sshell_db->keyExists("tag_velocity")) {
            d_tag_velocity = sshell_db->getDouble("tag_velocity");
         }
         if (sshell_db->keyExists("tag_width")) {
            d_tag_width = sshell_db->getDouble("tag_width");
         }
      }
      found = true;
   }

   if (!found) {
      TBOX_ERROR(
         d_object_name << ": Could not identify problem.\n"
                       << d_geom_problem << " is not a valid input."
                       << std::endl);
   }

   // ------------------------------

   /*
    * Block refine boxes.  These are supplied via input of the
    * form:
    *
    *    refine_boxes_<block number>_<level number>
    *
    * For example,
    *
    *    refine_boxes_2_0 = [(0,0),(9,9)],...
    *
    * would specify the refinement region on block 2, level 0.
    *
    */
   d_refine_boxes.resize(d_nblocks);
   for (hier::BlockId::block_t nb = 0; nb < d_nblocks; ++nb) {

      // see what the max number of levels is
      int max_ln = 0;
      int ln;
      for (ln = 0; ln < 10; ++ln) {
         snprintf(block_name, MBLK_GEOM_NUM_CHARS, "refine_boxes_%d_%d", nb, ln);
         if (db->keyExists(block_name)) {
            ++max_ln;
         }
      }
      d_refine_boxes[nb].resize(max_ln);

      for (ln = 0; ln < max_ln; ++ln) {
         snprintf(block_name, MBLK_GEOM_NUM_CHARS, "refine_boxes_%d_%d", nb, ln);
         if (db->keyExists(block_name)) {
            std::vector<tbox::DatabaseBox> db_box_vector =
               db->getDatabaseBoxVector(block_name);
            d_refine_boxes[nb][ln] = db_box_vector;
         } else {
            TBOX_ERROR(
               d_object_name << ": input entry `"
                             << block_name << "' does not exist."
                             << std::endl);
         }
      }
   }

}

/*
 *************************************************************************
 *
 * build the array of local blocks
 *
 *************************************************************************
 */

void MblkGeometry::buildLocalBlocks(
   const hier::Box& pbox,                                     // the patch box
   const hier::Box& domain,                                   // the block box
   const hier::BlockId::block_t block_number,
   hier::BlockId::block_t* dom_local_blocks)                                     // this returns the blocks neighboring this patch
{
   // by default single block simulation
   for (int i = 0; i < 6; ++i) {
      dom_local_blocks[i] = block_number;
   }

   if (d_geom_problem == "TRILINEAR") {
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast = pbox.upper();
      if (ifirst(0) == domain.lower(0)) dom_local_blocks[0] =
            d_tri_nbr[block_number][0];
      if (ifirst(1) == domain.lower(1)) dom_local_blocks[1] =
            d_tri_nbr[block_number][1];
      if (ifirst(2) == domain.lower(2)) dom_local_blocks[2] =
            d_tri_nbr[block_number][2];

      if (ilast(0) == domain.upper(0)) dom_local_blocks[3] =
            d_tri_nbr[block_number][3];
      if (ilast(1) == domain.upper(1)) dom_local_blocks[4] =
            d_tri_nbr[block_number][4];
      if (ilast(2) == domain.upper(2)) dom_local_blocks[5] =
            d_tri_nbr[block_number][5];
   }

   if (d_geom_problem == "SPHERICAL_SHELL") {
      if (d_sshell_type == "OCTANT") {

         //
         // process the block neighbors for this grid
         // using the order, this could be improved
         //
         // 0 imin, 1 jmin, 2 kmin
         // 3 imax, 4 jmax, 5 kmax
         //
         //
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast = pbox.upper();

         int jmin = domain.lower(1);
         int jmax = domain.upper(1);
         int kmax = domain.upper(2);

         if (block_number == 0) {
            if (ilast(1) == jmax) dom_local_blocks[4] = 1;
            if (ilast(2) == kmax) dom_local_blocks[5] = 2;
         }

         if (block_number == 1) {
            if (ifirst(1) == jmin) dom_local_blocks[1] = 0;
            if (ilast(2) == kmax) dom_local_blocks[5] = 2;
         }

         if (block_number == 2) {
            if (ilast(1) == jmax) dom_local_blocks[4] = 0;
            if (ilast(2) == kmax) dom_local_blocks[5] = 1;
         }
      }
   }

}

/*
 *************************************************************************
 *
 * Build grid on patch for supplied inputs for different goemetries
 *
 *************************************************************************
 */

void MblkGeometry::buildGridOnPatch(
   const hier::Patch& patch,
   const hier::Box& domain,
   const int xyz_id,
   const hier::BlockId::block_t block_number,
   hier::BlockId::block_t* dom_local_blocks)                                    // this returns the blocks neighboring this patch
{
   buildLocalBlocks(patch.getBox(), domain, block_number, dom_local_blocks);

   if (d_geom_problem == "CARTESIAN") {
      buildCartesianGridOnPatch(patch,
         domain,
         xyz_id);
   }

   if (d_geom_problem == "WEDGE") {
      buildWedgeGridOnPatch(patch,
         domain,
         xyz_id,
         block_number);
   }

   if (d_geom_problem == "TRILINEAR") {
      buildTrilinearGridOnPatch(patch,
         domain,
         xyz_id,
         block_number);
   }

   if (d_geom_problem == "SPHERICAL_SHELL") {
      buildSShellGridOnPatch(patch,
         domain,
         xyz_id,
         block_number);
   }
}

/*
 *************************************************************************
 *
 * Build the Cartesian grid
 *
 *************************************************************************
 */

void MblkGeometry::buildCartesianGridOnPatch(
   const hier::Patch& patch,
   const hier::Box& domain,
   const int xyz_id)
{
   //
   // compute dx
   //
   hier::Index lower(domain.lower());
   hier::Index upper(domain.upper());
   hier::Index diff(upper - lower + hier::Index(lower.getDim(), 1));

   double dx[SAMRAI::MAX_DIM_VAL];
   double xlo[SAMRAI::MAX_DIM_VAL];
   for (int i = 0; i < d_dim.getValue(); ++i) {
      dx[i] = (d_cart_xhi[0][i] - d_cart_xlo[0][i]) / (double)diff(i);
      xlo[i] = d_cart_xlo[0][i];
   }

   //
   // get the coordinates array information
   //
   std::shared_ptr<pdat::NodeData<double> > xyz(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
         patch.getPatchData(xyz_id)));

   TBOX_ASSERT(xyz);

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();
   hier::IntVector nghost_cells = xyz->getGhostCellWidth();

   int nd_imin = ifirst(0) - nghost_cells(0);
   int nd_imax = ilast(0) + 1 + nghost_cells(0);
   int nd_jmin = ifirst(1) - nghost_cells(1);
   int nd_jmax = ilast(1) + 1 + nghost_cells(1);
   int nd_kmin = ifirst(2) - nghost_cells(2);
   int nd_kmax = ilast(2) + 1 + nghost_cells(2);
   int nd_nx = nd_imax - nd_imin + 1;
   int nd_ny = nd_jmax - nd_jmin + 1;
   int nd_nxny = nd_nx * nd_ny;

   double* x = xyz->getPointer(0);
   double* y = xyz->getPointer(1);
   double* z = xyz->getPointer(2);

   //
   // ----------- set the nodal positions
   //
   for (int k = nd_kmin; k <= nd_kmax; ++k) {
      for (int j = nd_jmin; j <= nd_jmax; ++j) {
         for (int i = nd_imin; i <= nd_imax; ++i) {

            int ind = POLY3(i, j, k, nd_imin, nd_jmin, nd_kmin, nd_nx, nd_nxny);
            x[ind] = xlo[0] + i * dx[0];
            y[ind] = xlo[1] + j * dx[1];
            if (d_dim > tbox::Dimension(2)) {
               z[ind] = xlo[2] + k * dx[2];
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Build the Wedge grid
 *
 *************************************************************************
 */

void MblkGeometry::buildWedgeGridOnPatch(
   const hier::Patch& patch,
   const hier::Box& domain,
   const int xyz_id,
   const hier::BlockId::block_t block_number)
{
   //
   // Set dx (dr, dth, dz) for the level
   //
   double dx[SAMRAI::MAX_DIM_VAL];

   double nr = (domain.upper(0) - domain.lower(0) + 1);
   double nth = (domain.upper(1) - domain.lower(1) + 1);
   double nz = (domain.upper(2) - domain.lower(2) + 1);
   dx[0] = (d_wedge_rmax[0] - d_wedge_rmin[0]) / nr;
   dx[1] = (d_wedge_thmax - d_wedge_thmin) / nth;
   if (d_dim > tbox::Dimension(2)) {
      dx[2] = (d_wedge_zmax - d_wedge_zmin) / nz;
   } else {
      dx[2] = 0.0;
   }

   std::shared_ptr<pdat::NodeData<double> > xyz(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
         patch.getPatchData(xyz_id)));

   TBOX_ASSERT(xyz);

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();
   hier::IntVector nghost_cells = xyz->getGhostCellWidth();

   int nd_imin = ifirst(0) - nghost_cells(0);
   int nd_imax = ilast(0) + 1 + nghost_cells(0);
   int nd_jmin = ifirst(1) - nghost_cells(1);
   int nd_jmax = ilast(1) + 1 + nghost_cells(1);
   int nd_kmin = ifirst(2) - nghost_cells(2);
   int nd_kmax = ilast(2) + 1 + nghost_cells(2);
   int nd_nx = nd_imax - nd_imin + 1;
   int nd_ny = nd_jmax - nd_jmin + 1;
   int nd_nxny = nd_nx * nd_ny;

   double* x = xyz->getPointer(0);
   double* y = xyz->getPointer(1);
   double* z = xyz->getPointer(2);

   //
   // ----------- set the wedge nodal positions
   //

   for (int k = nd_kmin; k <= nd_kmax; ++k) {
      for (int j = nd_jmin; j <= nd_jmax; ++j) {
         for (int i = nd_imin; i <= nd_imax; ++i) {

            int ind = POLY3(i, j, k, nd_imin, nd_jmin, nd_kmin, nd_nx, nd_nxny);

            double r = d_wedge_rmin[block_number] + dx[0] * (i);
            double th = d_wedge_thmin + dx[1] * (j);

            double xx = r * cos(th);
            double yy = r * sin(th);
            double zz = d_wedge_zmin + dx[2] * (k);

            x[ind] = xx;
            y[ind] = yy;
            z[ind] = zz;
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Build the Trilinear grid
 *
 *************************************************************************
 */

void MblkGeometry::buildTrilinearGridOnPatch(
   const hier::Patch& patch,
   const hier::Box& domain,
   const int xyz_id,
   const hier::BlockId::block_t block_number)
{
   //
   // Set dx (dr, dth, dz) for the level
   //
   double nx = (domain.upper(0) - domain.lower(0) + 1);
   double ny = (domain.upper(1) - domain.lower(1) + 1);
   double nz = (domain.upper(2) - domain.lower(2) + 1);

   std::shared_ptr<pdat::NodeData<double> > xyz(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
         patch.getPatchData(xyz_id)));

   TBOX_ASSERT(xyz);

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();
   hier::IntVector nghost_cells = xyz->getGhostCellWidth();

   int nd_imin = ifirst(0) - nghost_cells(0);
   int nd_imax = ilast(0) + 1 + nghost_cells(0);
   int nd_jmin = ifirst(1) - nghost_cells(1);
   int nd_jmax = ilast(1) + 1 + nghost_cells(1);
   int nd_kmin = ifirst(2) - nghost_cells(2);
   int nd_kmax = ilast(2) + 1 + nghost_cells(2);
   int nd_nx = nd_imax - nd_imin + 1;
   int nd_ny = nd_jmax - nd_jmin + 1;
   int nd_nxny = nd_nx * nd_ny;

   double* x = xyz->getPointer(0);
   double* y = xyz->getPointer(1);
   double* z = xyz->getPointer(2);

   double* bx = d_tri_x[block_number];
   double* by = d_tri_y[block_number];
   double* bz = d_tri_z[block_number];

   int mxp = d_tri_nxp[block_number]; // nodal base mesh size
   int myp = d_tri_nyp[block_number];
   int mzp = d_tri_nzp[block_number];
   int mxpmyp = mxp * myp;

   int mx = mxp - 1; // zonal base mesh size
   int my = myp - 1;
   int mz = mzp - 1;

   double rx = nx / mx;
   double ry = ny / my;
   double rz = nz / mz;

   //
   // ----------- compute the nodal tri-linear interpolation
   //

   for (int k = nd_kmin; k <= nd_kmax; ++k) {
      for (int j = nd_jmin; j <= nd_jmax; ++j) {
         for (int i = nd_imin; i <= nd_imax; ++i) {

            int ind = POLY3(i, j, k, nd_imin, nd_jmin, nd_kmin, nd_nx, nd_nxny);

            double ric = ((double)i) / rx;
            double rjc = ((double)j) / ry;
            double rkc = ((double)k) / rz;

            int ic = (int)MAX(0, MIN(mx - 1, ric));
            int jc = (int)MAX(0, MIN(my - 1, rjc));
            int kc = (int)MAX(0, MIN(mz - 1, rkc));

            double xi = ric - ic;
            double eta = rjc - jc;
            double zeta = rkc - kc;

            int n1 = POLY3(ic, jc, kc, 0, 0, 0, mxp, mxpmyp);
            int n2 = POLY3(ic + 1, jc, kc, 0, 0, 0, mxp, mxpmyp);
            int n3 = POLY3(ic, jc + 1, kc, 0, 0, 0, mxp, mxpmyp);
            int n4 = POLY3(ic + 1, jc + 1, kc, 0, 0, 0, mxp, mxpmyp);
            int n5 = POLY3(ic, jc, kc + 1, 0, 0, 0, mxp, mxpmyp);
            int n6 = POLY3(ic + 1, jc, kc + 1, 0, 0, 0, mxp, mxpmyp);
            int n7 = POLY3(ic, jc + 1, kc + 1, 0, 0, 0, mxp, mxpmyp);
            int n8 = POLY3(ic + 1, jc + 1, kc + 1, 0, 0, 0, mxp, mxpmyp);

            double xx = (bx[n1] * (1 - xi) * (1 - eta) * (1 - zeta)
                         + bx[n2] * (xi) * (1 - eta) * (1 - zeta)
                         + bx[n3] * (1 - xi) * (eta) * (1 - zeta)
                         + bx[n4] * (xi) * (eta) * (1 - zeta)
                         + bx[n5] * (1 - xi) * (1 - eta) * (zeta)
                         + bx[n6] * (xi) * (1 - eta) * (zeta)
                         + bx[n7] * (1 - xi) * (eta) * (zeta)
                         + bx[n8] * (xi) * (eta) * (zeta));

            double yy = (by[n1] * (1 - xi) * (1 - eta) * (1 - zeta)
                         + by[n2] * (xi) * (1 - eta) * (1 - zeta)
                         + by[n3] * (1 - xi) * (eta) * (1 - zeta)
                         + by[n4] * (xi) * (eta) * (1 - zeta)
                         + by[n5] * (1 - xi) * (1 - eta) * (zeta)
                         + by[n6] * (xi) * (1 - eta) * (zeta)
                         + by[n7] * (1 - xi) * (eta) * (zeta)
                         + by[n8] * (xi) * (eta) * (zeta));

            double zz = (bz[n1] * (1 - xi) * (1 - eta) * (1 - zeta)
                         + bz[n2] * (xi) * (1 - eta) * (1 - zeta)
                         + bz[n3] * (1 - xi) * (eta) * (1 - zeta)
                         + bz[n4] * (xi) * (eta) * (1 - zeta)
                         + bz[n5] * (1 - xi) * (1 - eta) * (zeta)
                         + bz[n6] * (xi) * (1 - eta) * (zeta)
                         + bz[n7] * (1 - xi) * (eta) * (zeta)
                         + bz[n8] * (xi) * (eta) * (zeta));

            x[ind] = xx;
            y[ind] = yy;
            z[ind] = zz;
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Build the spherical shell grid
 *
 *************************************************************************
 */

void MblkGeometry::buildSShellGridOnPatch(
   const hier::Patch& patch,
   const hier::Box& domain,
   const int xyz_id,
   const hier::BlockId::block_t block_number)
{

   bool xyz_allocated = patch.checkAllocated(xyz_id);
   if (!xyz_allocated) {
      TBOX_ERROR("xyz data not allocated" << std::endl);
      //patch.allocatePatchData(xyz_id);
   }

   std::shared_ptr<pdat::NodeData<double> > xyz(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
         patch.getPatchData(xyz_id)));

   TBOX_ASSERT(xyz);

   if (d_dim == tbox::Dimension(3)) {

      const hier::Index ifirst = patch.getBox().lower();
      const hier::Index ilast = patch.getBox().upper();
      hier::IntVector nghost_cells = xyz->getGhostCellWidth();

      int nd_imin = ifirst(0) - nghost_cells(0);
      int nd_imax = ilast(0) + 1 + nghost_cells(0);
      int nd_jmin = ifirst(1) - nghost_cells(1);
      int nd_jmax = ilast(1) + 1 + nghost_cells(1);
      int nd_kmin = ifirst(2) - nghost_cells(2);
      int nd_kmax = ilast(2) + 1 + nghost_cells(2);
      int nd_nx = nd_imax - nd_imin + 1;
      int nd_ny = nd_jmax - nd_jmin + 1;
      int nd_nxny = nd_nx * nd_ny;

      double* x = xyz->getPointer(0);
      double* y = xyz->getPointer(1);
      double* z = xyz->getPointer(2);

      bool found = false;

      int nrad = (domain.upper(0) - domain.lower(0) + 1);
      int nth = (domain.upper(1) - domain.lower(1) + 1);
      int nphi = (domain.upper(2) - domain.lower(2) + 1);

      /*
       * If its a solid shell, its a single block and dx = dr, dth, dphi
       */
      if (d_sshell_type == "SOLID") {

         double dx[3];
         dx[0] = (d_sshell_rmax - d_sshell_rmin) / (double)nrad;
         dx[1] =
            2.0 * tbox::MathUtilities<double>::Abs(d_sangle_thmin)
            / (double)nth;
         dx[2] =
            2.0 * tbox::MathUtilities<double>::Abs(d_sangle_thmin)
            / (double)nphi;

         //
         // step in a radial direction in x and set y and z appropriately
         // for a solid angle we go -th to th and -phi to phi
         //
         for (int k = nd_kmin; k <= nd_kmax; ++k) {
            for (int j = nd_jmin; j <= nd_jmax; ++j) {

               double theta = d_sangle_thmin + j * dx[1]; // dx used for dth
               double phi = d_sangle_thmin + k * dx[2];

               double xface = cos(theta) * cos(phi);
               double yface = sin(theta) * cos(phi);
               double zface = sin(phi);

               for (int i = nd_imin; i <= nd_imax; ++i) {

                  int ind = POLY3(i,
                        j,
                        k,
                        nd_imin,
                        nd_jmin,
                        nd_kmin,
                        nd_nx,
                        nd_nxny);

                  double r = d_sshell_rmin + dx[0] * (i);

                  double xx = r * xface;
                  double yy = r * yface;
                  double zz = r * zface;

                  x[ind] = xx;
                  y[ind] = yy;
                  z[ind] = zz;
               }
            }
         }

         found = true;
      }

      /*
       * If its an octant problem, then its got multiple (three) blocks
       */
      if (d_sshell_type == "OCTANT") {

         double drad = (d_sshell_rmax - d_sshell_rmin) / nrad;

         //
         // as in the solid angle we go along a radial direction in
         // x setting y and z appropriately, but here we have logic for
         // the block we are in.  This is contained in the dispOctant.m
         // matlab code.
         //
         for (int k = nd_kmin; k <= nd_kmax; ++k) {
            for (int j = nd_jmin; j <= nd_jmax; ++j) {

               //
               // compute the position on the unit sphere for our radial line
               //
               double xface, yface, zface;
               computeUnitSphereOctant(block_number, nth, j, k,
                  &xface, &yface, &zface);

               for (int i = nd_imin; i <= nd_imax; ++i) {
                  int ind = POLY3(i,
                        j,
                        k,
                        nd_imin,
                        nd_jmin,
                        nd_kmin,
                        nd_nx,
                        nd_nxny);

                  double r = d_sshell_rmin + drad * (i);

                  double xx = r * xface;
                  double yy = r * yface;
                  double zz = r * zface;

                  x[ind] = xx;
                  y[ind] = yy;
                  z[ind] = zz;
               }
            }
         }
         found = true;
      }

      if (!found) {
         TBOX_ERROR(
            d_object_name << ": "
                          << "spherical shell nodal positions for "
                          << d_sshell_type
                          << " not found" << std::endl);
      }

   }

}

/*
 *************************************************************************
 *
 * For a given j, k, compute the position on the unit sphere with
 * the supplied block number.  "nth" is the number of cells in the theta
 * direction (it should be the same for all blocks).  The three faces
 * "xface", "yface", "zface" are the different block faces.  The
 * code that performs this operation is in Dave's dispOctant.m matlab
 * code.
 *
 *************************************************************************
 */
void MblkGeometry::computeUnitSphereOctant(
   hier::BlockId::block_t nblock,
   int nth,
   int j,
   int k,
   double* xface,
   double* yface,
   double* zface)
{
   static int jmn[3] = { 0, 0, 2 }; // matrix of rotations
   static int jmx[3] = { 1, 1, 1 };
   static int kmn[3] = { 0, 1, 0 };
   static int kmx[3] = { 2, 2, 2 };

   //
   // coefficients are of the form
   //
   // jto = j0 + jslp*jfrom + k0 + kslp*kfrom
   // kto = j0 + jslp*jfrom + k0 + kslp*kfrom

   static int jcoef[3][3][4] = // how to convert j from one block to another
                               // [from][to][4]
   { { { 0, 1, 0, 0 },         // [0][0]
       { -nth, 1, 0, 0 },      // [0][1]
       { 0, 0, 0, 1 } },         // [0][2]

     { { nth, 1, 0, 0 },         // [1][0]
       { 0, 1, 0, 0 },         // [1][1]
       { 0, 0, nth, -1 } },      // [1][2]

     { { 0, 0, 2 * nth, -1 },    // [2][0]
       { 0, 0, 1, 0 },         // [2][1]
       { nth, -1, 0, 0 } } };    // [2][2]

   static int kcoef[3][3][4] = // the k conversion
/*    {   0,  0,    0,  1,   // [0][0]
*      0,  0,    0,  1,   // [0][1]
*      0,  0,2*nth, -1,   // [0][2]
*
*      0,  0,    0,  1,   // [1][0]
*      0,  0,    0,  1,   // [1][1]
*      0,  0,2*nth, -1,   // [1][2]
*
*      0,  1,    0,  0,   // [2][0]
*  2*nth, -1,    0,  0,   // [2][1]
*      0,  0,    0,  1 }; // [2][2] */
   { { { 0, 0, 0, 1 },          // [0][0]
       { 0, 0, 0, 1 },         // [0][1]
       { 0, 0, 2 * nth, -1 } },  // [0][2]

     { { 0, 0, 0, 1 },          // [1][0]
       { 0, 0, 0, 1 },         // [1][1]
       { 0, 0, 2 * nth, -1 } },  // [1][2]

     { { 0, 1, 0, 0 },          // [2][0]
       { 2 * nth, -1, 0, 0 },  // [2][1]
       { 0, 0, 0, 1 } } };       // [2][2]

   static double pio4 = 0.25 * 3.14159265358979325873406851315;

   //
   // decide which position in the unit sphere we break out on
   //
   // nblock = 1, xface = 2, yface, = 3, zface

   int tb = static_cast<int>(nblock);
   int tj = j;
   int tk = k;

   int sb = tb;
   int ttj = tj;
   int ttk = tk;
   if (j < 0) {
      tb = jmn[sb];
      tj = jcoef[sb][tb][0] + jcoef[sb][tb][1] * ttj + jcoef[sb][tb][2]
         + jcoef[sb][tb][3] * ttk;
      tk = kcoef[sb][tb][0] + kcoef[sb][tb][1] * ttj + kcoef[sb][tb][2]
         + kcoef[sb][tb][3] * ttk;
   }

   sb = tb;
   ttj = tj;
   ttk = tk;
   if (j > nth) {
      tb = jmx[sb];
      tj = jcoef[sb][tb][0] + jcoef[sb][tb][1] * ttj + jcoef[sb][tb][2]
         + jcoef[sb][tb][3] * ttk;
      tk = kcoef[sb][tb][0] + kcoef[sb][tb][1] * ttj + kcoef[sb][tb][2]
         + kcoef[sb][tb][3] * ttk;
   }

   sb = tb;
   ttj = tj;
   ttk = tk;
   if (k < 0) {   // this catches corner terms now
      tb = kmn[sb];
      tj = jcoef[sb][tb][0] + jcoef[sb][tb][1] * ttj + jcoef[sb][tb][2]
         + jcoef[sb][tb][3] * ttk;
      tk = kcoef[sb][tb][0] + kcoef[sb][tb][1] * ttj + kcoef[sb][tb][2]
         + kcoef[sb][tb][3] * ttk;
   }

   sb = tb;
   ttj = tj;
   ttk = tk;
   if (k > nth) {
      tb = kmx[sb];
      tj = jcoef[sb][tb][0] + jcoef[sb][tb][1] * ttj + jcoef[sb][tb][2]
         + jcoef[sb][tb][3] * ttk;
      tk = kcoef[sb][tb][0] + kcoef[sb][tb][1] * ttj + kcoef[sb][tb][2]
         + kcoef[sb][tb][3] * ttk;
   }

   //
   // once we know where we are in our block, we can return our position on the unit sphere
   //
   if (tb == 0) {   // we go along an x face  (block 0)
      double lj = tj / (double)nth;
      double lk = tk / (double)nth;

      double sz = sin(pio4 * lj);
      double cz = cos(pio4 * lj);
      double sy = sin(pio4 * lk);
      double cy = cos(pio4 * lk);

      double den = sqrt(1. - sy * sy * sz * sz);

      *xface = cy * cz / den;
      *yface = cy * sz / den;
      *zface = sy * cz / den;
   } else if (tb == 1) { // a y face (block 1)
      double li = 1 - tj / (double)nth;
      double lk = tk / (double)nth;

      double sx = sin(pio4 * lk);
      double cx = cos(pio4 * lk);
      double sz = sin(pio4 * li);
      double cz = cos(pio4 * li);

      double den = sqrt(1. - sx * sx * sz * sz);

      *xface = cx * sz / den;
      *yface = cx * cz / den;
      *zface = sx * cz / den;
   } else { // a z face (block 2)
      double li = tj / (double)nth; // 1 - tj;
      double lj = tk / (double)nth; // 1 - tk;

      double sx = sin(pio4 * lj);
      double cx = cos(pio4 * lj);
      double sy = sin(pio4 * li);
      double cy = cos(pio4 * li);

      double den = sqrt(1. - sx * sx * sy * sy);

      *xface = cx * sy / den;
      *yface = sx * cy / den;
      *zface = cx * cy / den;
   }

}

/*
 *************************************************************************
 *
 * Compute the rotations for the particular block number.  The
 * "local_blocks" argument is basically saying who is the ne
 *
 *************************************************************************
 */
void computeBlocksOctant(
   const hier::Box& bb,
   int local_blocks[6],
   int nblock,
   int nth)
{
   const hier::Index ifirst = bb.lower();
   const hier::Index ilast = bb.upper();

   local_blocks[0] = nblock; // imin stays local
   local_blocks[3] = nblock; // jmin stays local

   static int jmn[3] = { 0, 0, 2 }; // matrix of rotations
   static int jmx[3] = { 1, 1, 1 };
   static int kmn[3] = { 0, 1, 0 };
   static int kmx[3] = { 2, 2, 2 };

   //
   // bounds of the patch zones go from 0 to nth-1
   //
   if (ifirst(1) <= 0) local_blocks[1] = jmn[nblock];
   if (ifirst(2) <= 0) local_blocks[2] = kmn[nblock];

   if (ilast(1) >= nth - 1) local_blocks[4] = jmx[nblock];
   if (ilast(2) >= nth - 1) local_blocks[5] = kmx[nblock];
}
