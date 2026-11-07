/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   set geometry for multiblock domain
 *
 ************************************************************************/

#ifndef included_MblkGeometryXD
#define included_MblkGeometryXD

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Serializable.h"

#include <memory>

using namespace SAMRAI;

/*!
 * This class creates the mapped multi-block grid geometry used
 * for calculations in the MblkLinAdv code.  The supported grid types
 * include Cartesian, Wedge, and Spherical shell.  The spherical shell
 * case is a full multi-block grid with 3 blocks.
 */
class MblkGeometry
{
public:
   /*!
    * Reads geometry information from the "MblkGeometry" input file
    * entry.
    */
   MblkGeometry(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<tbox::Database>& input_db,
      std::shared_ptr<hier::BaseGridGeometry>& grid_geom);

   ~MblkGeometry();

   /*!
    * Return the geometry type (CARTESIAN, WEDGE, or SPHERICAL_SHELL)
    */
   std::string
   getGeometryType();

   /*!
    * Return the user-specified refine boxes, given a block and
    * level number
    */
   bool
   getRefineBoxes(
      hier::BoxContainer& refine_boxes,
      const int block_number,
      const int level_number);

   /*!
    * Build mapped grid on patch.  The method defers the actual grid
    * construction to private members, depending on the geometry
    * choice in input.
    */
   void
   buildGridOnPatch(
      const hier::Patch& patch,
      const hier::Box& domain,
      const int xyz_id,
      const int level_number,
      const int block_number);

   /*!
    * Access the stored dx
    */

   void
   getDx(
      const hier::Box& domain,
      const int level_number,
      const int block_number,
      double* dx);

   void
   getDx(
      const int level_number,
      const int block_number,
      double* dx);

   /*!
    * Access the block rotation
    */
   int
   getBlockRotation(
      const int block_number);

   /*!
    * Tag cells for the octant problem.
    */
   void
   tagOctantCells(
      hier::Patch& patch,
      const int xyz_id,
      std::shared_ptr<pdat::CellData<int> >& temp_tags,
      const double regrid_time,
      const int refine_tag_val);

private:
   /*
    * Read data members from input.
    */
   void
   getFromInput(
      std::shared_ptr<tbox::Database> input_db,
      bool is_from_restart);

   /*
    * Cartesian grid construction.
    */
   void
   setCartesianMetrics(
      const hier::Box& domain,
      const int level_number,
      const int block_number);

   void
   buildCartesianGridOnPatch(
      const hier::Patch& patch,
      const int xyz_id,
      const int level_number,
      const int block_number);

   /*
    * Wedge grid construction.
    */
   void
   setWedgeMetrics(
      const hier::Box& domain,
      const int level_number);

   void
   buildWedgeGridOnPatch(
      const hier::Patch& patch,
      const int xyz_id,
      const int level_number,
      const int block_number);

   /*
    * Spherical shell grid construction
    */
   void
   setSShellMetrics(
      const hier::Box& domain,
      const int level_number);

   void
   buildSShellGridOnPatch(
      const hier::Patch& patch,
      const hier::Box& domain,
      const int xyz_id,
      const int level_number,
      const int block_number);

   /*
    * For the spherical shell construction, i always points in the r direction
    * and j,k are points on the shell.  Given a certain j,k compute the
    * unit sphere locations for block face (actual xyz is computed
    * by x = r*xface, y = r*yface, z = r*zface.  Note that the size
    * in the theta direction (nth) should be the same for each block.
    */
   void
   computeUnitSphereOctant(
      int nblock,
      int nth,
      int j,
      int k,
      double* xface,
      double* yface,
      double* zface);

   /*
    * Geometry type.  Choices are CARTESIAN, WEDGE, SPHERICAL_SHELL
    */
   std::string d_geom_problem;
   std::string d_object_name;

   std::shared_ptr<hier::BaseGridGeometry> d_grid_geom;

   const tbox::Dimension d_dim;

   /*
    * The number of blocks and the set of skelton grid geometries that make
    * up a multiblock mesh.
    */
   int d_nblocks;
   std::vector<std::vector<bool> > d_metrics_set;

   /*
    * The grid spacing.  For cartesian, d_dx = (dx,dy,dz).  For wedge,
    * d_dx = (dr, dth, dz). For spherical shell, d_dx = (dr, dth, dphi)
    */
   std::vector<std::vector<std::vector<double> > > d_dx;

   /*
    * Cartesian inputs
    */
   std::vector<std::vector<double> > d_cart_xlo;
   std::vector<std::vector<double> > d_cart_xhi;

   /*
    * Wedge inputs
    */
   std::vector<double> d_wedge_rmin;
   std::vector<double> d_wedge_rmax;
   double d_wedge_thmin;
   double d_wedge_thmax;
   double d_wedge_zmin;
   double d_wedge_zmax;

   /*
    * Shell inputs
    */
   double d_sshell_rmin;
   double d_sshell_rmax;

   // options are SOLID, OCTANT
   std::string d_sshell_type;

   /*
    * For tagging in the spherical octant case
    */
   double d_tag_velocity;
   double d_tag_width;

   // if SOLID, need to read in these...
   double d_sangle_degrees;
   double d_sangle_thmin;

   /*
    * Specify block rotation.
    */
   std::vector<int> d_block_rotation;

   /*
    * Refine boxes for different blocks/levels
    */
   std::vector<std::vector<hier::BoxContainer> > d_refine_boxes;

};

#endif // included_MblkGeometry
