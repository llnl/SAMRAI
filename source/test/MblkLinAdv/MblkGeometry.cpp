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

#define POLY3(i, j, k, imin, jmin, kmin, nx, nxny) \
   ((i - imin) + (j - jmin) * (nx) + (k - kmin) * (nxny))

#define MBLK_GEOM_NUM_CHARS 128

/*
 *************************************************************************
 *
 * This class creates the mapped multi-block grid geometry used
 * for calculations in the MblkLinAdv code.
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

   d_metrics_set.resize(10);
   for (int i = 0; i < 10; ++i) {
      d_metrics_set[i] = false;
   }

   /*
    * Initialize object with data read from given input/restart databases.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
//      getFromRestart();  // ADD
   }
   getFromInput(input_db, is_from_restart);

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
 * Return the geometry (CARTESIAN, WEDGE, or SPHERICAL_SHELL)
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
 * Tag cells for spherical octant problem
 *
 *************************************************************************
 */
void MblkGeometry::tagOctantCells(
   hier::Patch& patch,
   const int xyz_id,
   std::shared_ptr<pdat::CellData<int> >& temp_tags,
   const double regrid_time,
   const int refine_tag_val)
{
   TBOX_ASSERT(d_geom_problem == "SPHERICAL_SHELL" &&
      d_sshell_type == "OCTANT");
   TBOX_ASSERT(temp_tags);

   std::shared_ptr<pdat::NodeData<double> > xyz(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
         patch.getPatchData(xyz_id)));
   TBOX_ASSERT(xyz);

   if (d_dim == tbox::Dimension(3)) {
      /*
       * Tag in X direction only
       */
      double xtag_loc_lo = d_sshell_rmin
         + (regrid_time * d_tag_velocity) - (0.5 * d_tag_width);
      double xtag_loc_hi = d_sshell_rmin
         + (regrid_time * d_tag_velocity) + (0.5 * d_tag_width);

      hier::Box pbox = patch.getBox();
      for (int k = pbox.lower(2); k <= pbox.upper(2) + 1; ++k) {
         for (int j = pbox.lower(1); j <= pbox.upper(1) + 1; ++j) {
            for (int i = pbox.lower(0); i <= pbox.upper(0) + 1; ++i) {
               hier::Index ic(i, j, k);
               pdat::NodeIndex node(ic, pdat::NodeIndex::LLL);
               hier::Index icm1(i - 1, j - 1, k - 1);
               pdat::CellIndex cell(icm1);

               double node_x_loc = (*xyz)(node, 0);

               if ((node_x_loc > xtag_loc_lo) &&
                   (node_x_loc < xtag_loc_hi)) {
                  (*temp_tags)(cell) = refine_tag_val;
               }
            }
         }
      }
   }

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
   std::shared_ptr<tbox::Database> input_db,
   bool is_from_restart)
{
   TBOX_ASSERT(input_db);

   NULL_USE(is_from_restart);

   std::shared_ptr<tbox::Database> db(input_db->getDatabase("MblkGeometry"));

   d_geom_problem = db->getString("problem_type");

   bool found = false;
   int i;
   char block_name[MBLK_GEOM_NUM_CHARS];
   double temp_domain[SAMRAI::MAX_DIM_VAL];

   /*
    * Cartesian geometry
    */
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

   /*
    * Wedge geometry
    */
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

         if (d_dim == tbox::Dimension(3)) {
            // Z min/max
            d_wedge_zmin = wedge_db->getDouble("zmin");
            d_wedge_zmax = wedge_db->getDouble("zmax");
         }

      }
      found = true;
   }

   /*
    * Spherical shell
    */
   if (d_geom_problem == "SPHERICAL_SHELL") {

      /*
       * This case only works in 3 dimensions
       */
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

   /*
    * Block rotation
    */
   d_block_rotation.resize(d_nblocks);
   for (hier::BlockId::block_t nb = 0; nb < d_nblocks; ++nb) {
      d_block_rotation[nb] = 0;
      snprintf(block_name, MBLK_GEOM_NUM_CHARS, "rotation_%d", nb);
      if (db->keyExists(block_name)) {
         d_block_rotation[nb] = db->getInteger(block_name);
      }
   }

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
 * Build grid on patch for supplied inputs for different goemetries
 *
 *************************************************************************
 */

void MblkGeometry::buildGridOnPatch(
   const hier::Patch& patch,
   const hier::Box& domain,
   const int xyz_id,
   const int level_number,
   const hier::BlockId::block_t block_number)
{

   if (d_geom_problem == "CARTESIAN") {
      if (!d_metrics_set[level_number]) {
         setCartesianMetrics(domain,
            level_number);
      }
      buildCartesianGridOnPatch(patch,
         xyz_id,
         level_number,
         block_number);
   }

   if (d_geom_problem == "WEDGE") {
      if (!d_metrics_set[level_number]) {
         setWedgeMetrics(domain,
            level_number);
      }
      buildWedgeGridOnPatch(patch,
         xyz_id,
         level_number,
         block_number);
   }

   if (d_geom_problem == "SPHERICAL_SHELL") {
      if (!d_metrics_set[level_number]) {
         setSShellMetrics(domain,
            level_number);
      }
      buildSShellGridOnPatch(patch,
         domain,
         xyz_id,
         level_number,
         block_number);
   }

}

/*
 *************************************************************************
 *
 * Access the stored dx
 *
 *************************************************************************
 */

void MblkGeometry::getDx(
   const hier::Box& domain,
   const int level_number,
   double* dx)

{
   if (d_geom_problem == "CARTESIAN") {
      if (!d_metrics_set[level_number]) {
         setCartesianMetrics(domain,
            level_number);
      }
   }

   if (d_geom_problem == "WEDGE") {
      if (!d_metrics_set[level_number]) {
         setWedgeMetrics(domain,
            level_number);
      }
   }

   if (d_geom_problem == "SPHERICAL_SHELL") {
      if (!d_metrics_set[level_number]) {
         setSShellMetrics(domain,
            level_number);
      }
   }

   getDx(level_number,
      dx);

}

void MblkGeometry::getDx(
   const int level_number,
   double* dx)

{
   if (!d_metrics_set[level_number]) {
      TBOX_ERROR(
         d_object_name << ":metrics have not been set.\n"
                       << "Use the alternative 'getDx()' method call that "
                       << "takes in the domain." << std::endl);
   }

   for (int i = 0; i < d_dim.getValue(); ++i) {
      dx[i] = d_dx[level_number][i];
   }

}

/*
 *************************************************************************
 *
 * Access the block rotation
 *
 *************************************************************************
 */

int MblkGeometry::getBlockRotation(
   const hier::BlockId::block_t block_number)

{
   return d_block_rotation[block_number];
}

/*
 *************************************************************************
 *
 * Set the Cartesian metrics (dx)
 *
 *************************************************************************
 */

void MblkGeometry::setCartesianMetrics(
   const hier::Box& domain,
   const int level_number)
{
   hier::Index lower(domain.lower());
   hier::Index upper(domain.upper());
   hier::Index diff(upper - lower + hier::Index(lower.getDim(), 1));

   if (static_cast<int>(d_dx.size()) < (level_number + 1)) {
      d_dx.resize(level_number + 1);
      d_dx[level_number].resize(d_dim.getValue());
   }

   /*
    * Compute dx from first grid geometry block only.  Its assumed
    * to be uniform across the multiple blocks.
    */
   for (int i = 0; i < d_dim.getValue(); ++i) {
      d_dx[level_number][i] =
         (d_cart_xhi[0][i] - d_cart_xlo[0][i]) / (double)diff(i);
   }

   d_metrics_set[level_number] = true;

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
   const int xyz_id,
   const int level_number,
   const hier::BlockId::block_t block_number)
{

   std::shared_ptr<pdat::NodeData<double> > xyz(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
         patch.getPatchData(xyz_id)));

   TBOX_ASSERT(xyz);

   pdat::NodeIterator niend(pdat::NodeGeometry::end(patch.getBox()));
   for (pdat::NodeIterator ni(pdat::NodeGeometry::begin(patch.getBox()));
        ni != niend; ++ni) {
      pdat::NodeIndex node = *ni;
      if (d_block_rotation[block_number] == 0) {

         (*xyz)(node, 0) =
            d_cart_xlo[block_number][0] + node(0) * d_dx[level_number][0];
         (*xyz)(node, 1) =
            d_cart_xlo[block_number][1] + node(1) * d_dx[level_number][1];
         if (d_dim == tbox::Dimension(3)) {
            (*xyz)(node, 2) =
               d_cart_xlo[block_number][2] + node(2) * d_dx[level_number][2];
         }
      }
      if (d_block_rotation[block_number] == 1) { // I sideways, J down

         (*xyz)(node, 0) =
            d_cart_xlo[block_number][0] - node(0) * d_dx[level_number][0];
         (*xyz)(node, 1) =
            d_cart_xlo[block_number][1] + node(1) * d_dx[level_number][1];
         if (d_dim == tbox::Dimension(3)) {
            (*xyz)(node, 2) =
               d_cart_xlo[block_number][2] + node(2) * d_dx[level_number][2];
         }
      }

   }

#if 1
   tbox::plog << "xyz locations...."
              << "\tblock: " << block_number
              << "\tlevel: " << level_number
              << "\tpatch: " << patch.getLocalId()
              << std::endl;

   xyz->print(patch.getBox(), tbox::plog);
#endif

}

/*
 *************************************************************************
 *
 * Set the wedge metrics (dx)
 *
 *************************************************************************
 */

void MblkGeometry::setWedgeMetrics(
   const hier::Box& domain,
   const int level_number)
{

   //
   // Set dx (dr, dth, dz) for the level
   //
   d_dx.resize(level_number + 1);
   d_dx[level_number].resize(d_dim.getValue());

   double nr = (domain.upper(0) - domain.lower(0) + 1);
   double nth = (domain.upper(1) - domain.lower(1) + 1);
   d_dx[level_number][0] = (d_wedge_rmax[0] - d_wedge_rmin[0]) / nr;
   d_dx[level_number][1] = (d_wedge_thmax - d_wedge_thmin) / nth;

   if (d_dim == tbox::Dimension(3)) {
      double nz = (domain.upper(2) - domain.lower(2) + 1);
      d_dx[level_number][2] = (d_wedge_zmax - d_wedge_zmin) / nz;
   }

   d_metrics_set[level_number] = true;
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
   const int xyz_id,
   const int level_number,
   const hier::BlockId::block_t block_number)
{

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
   int nd_nx = nd_imax - nd_imin + 1;
   int nd_ny = nd_jmax - nd_jmin + 1;
   //int nd_nz   = nd_kmax - nd_kmin + 1;
   int nd_nxny = nd_nx * nd_ny;
   //int nd_nel  = nd_nx*nd_ny*nd_nz;

   double dx[SAMRAI::MAX_DIM_VAL];
   dx[0] = d_dx[level_number][0];
   dx[1] = d_dx[level_number][1];

   double* x = xyz->getPointer(0);
   double* y = xyz->getPointer(1);

   int nd_kmin;
   int nd_kmax;
   dx[2] = d_dx[level_number][2];
   double* z = 0;
   if (d_dim == tbox::Dimension(3)) {
      nd_kmin = ifirst(2) - nghost_cells(2);
      nd_kmax = ilast(2) + 1 + nghost_cells(2);
      dx[2] = d_dx[level_number][2];
      z = xyz->getPointer(2);
   } else {
      nd_kmin = 0;
      nd_kmax = 0;
   }

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

            x[ind] = xx;
            y[ind] = yy;

            if (d_dim == tbox::Dimension(3)) {
               double zz = d_wedge_zmin + dx[2] * (k);
               z[ind] = zz;
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Set the spherical shell metrics
 *
 *************************************************************************
 */

void MblkGeometry::setSShellMetrics(
   const hier::Box& domain,
   const int level_number)
{

   //
   // Set dx (drad, dth, dphi) for the level
   //
   d_dx.resize(level_number + 1);
   d_dx[level_number].resize(d_dim.getValue());

   double nrad = (domain.upper(0) - domain.lower(0) + 1);
   double nth = (domain.upper(1) - domain.lower(1) + 1);
   double nphi = 0;
   if (d_dim == tbox::Dimension(3)) {
      nphi = (domain.upper(2) - domain.lower(2) + 1);
   }

   /*
    * If its a solid shell, its a single block and dx = dr, dth, dphi
    */
   if (d_sshell_type == "SOLID") {

      d_dx[level_number][0] = (d_sshell_rmax - d_sshell_rmin) / nrad;
      d_dx[level_number][1] =
         2.0 * tbox::MathUtilities<double>::Abs(d_sangle_thmin) / nth;
      if (d_dim == tbox::Dimension(3)) {
         d_dx[level_number][2] =
            2.0 * tbox::MathUtilities<double>::Abs(d_sangle_thmin) / nphi;
      }
   } else {
      d_dx[level_number][0] = 0.0001;
      d_dx[level_number][1] = 0.0001;
      if (d_dim == tbox::Dimension(3)) {
         d_dx[level_number][2] = 0.0001;
      }
   }

   /*
    * If its an OCTANT shell, then everything is set in the
    * computeUnitSphereOctant() method so all we do here is allocate
    * space for d_dx.
    */
   d_metrics_set[level_number] = true;
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
   const int level_number,
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

      //int imin = ifirst(0);
      //int imax = ilast(0)  + 1;
      //int jmin = ifirst(1);
      //int jmax = ilast(1)  + 1;
      //int kmin = ifirst(2);
      //int kmax = ilast(2)  + 1;
      //int nx   = imax - imin + 1;
      //int ny   = jmax - jmin + 1;
      //int nxny = nx*ny;

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

         d_dx[level_number][0] = (d_sshell_rmax - d_sshell_rmin) / (double)nrad;
         d_dx[level_number][1] =
            2.0 * tbox::MathUtilities<double>::Abs(d_sangle_thmin)
            / (double)nth;
         d_dx[level_number][2] =
            2.0 * tbox::MathUtilities<double>::Abs(d_sangle_thmin)
            / (double)nphi;

         //
         // step in a radial direction in x and set y and z appropriately
         // for a solid angle we go -th to th and -phi to phi
         //
         for (int k = nd_kmin; k <= nd_kmax; ++k) {
            for (int j = nd_jmin; j <= nd_jmax; ++j) {

               double theta = d_sangle_thmin + j * d_dx[level_number][1]; // dx used for dth
               double phi = d_sangle_thmin + k * d_dx[level_number][2];

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

                  double r = d_sshell_rmin + d_dx[level_number][0] * (i);

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
       { -nth, 1, 0, 0 },    // [0][1]
       { 0, 0, 0, 1 } },       // [0][2]

     { { nth, 1, 0, 0 },      // [1][0]
       { 0, 1, 0, 0 },       // [1][1]
       { 0, 0, nth, -1 } },    // [1][2]

     { { 0, 0, 2 * nth, -1 }, // [2][0]
       { 0, 0, 1, 0 },       // [2][1]
       { nth, -1, 0, 0 } } };  // [2][2]

   static int kcoef[3][3][4] = // the k conversion
   { { { 0, 0, 0, 1 },        // [0][0]
       { 0, 0, 0, 1 },       // [0][1]
       { 0, 0, 2 * nth, -1 } }, // [0][2]

     { { 0, 0, 0, 1 },        // [1][0]
       { 0, 0, 0, 1 },       // [1][1]
       { 0, 0, 2 * nth, -1 } }, // [1][2]

     { { 0, 1, 0, 0 },        // [2][0]
       { 2 * nth, -1, 0, 0 }, // [2][1]
       { 0, 0, 0, 1 } } };    // [2][2]

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
