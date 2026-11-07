/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Simple Cartesian grid geometry for an AMR hierarchy.
 *
 ************************************************************************/

#ifndef included_geom_CartesianGridGeometry
#define included_geom_CartesianGridGeometry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/geom/GridGeometry.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Serializable.h"

#include <string>
#include <memory>

namespace SAMRAI {
namespace geom {

/**
 * Class CartesianGridGeometry provides simple Cartesian mesh geometry
 * management on an AMR hierarchy.  The mesh coordinates on each hierarchy
 * level are limited to mesh increments specified as DIM-tuple
 * (dx[0],...,dx[DIM-1]) and spatial coordinates of the lower and upper
 * corners of the smallest parallelepiped bounding the entire computational
 * domain.  The mesh increments on each level are defined with respect to
 * the coarsest hierarchy level and multiplying those values by the proper
 * refinement ratio.  This class sets geometry information on each patch in
 * an AMR hierarchy.  This class is derived from the GridGeometry base class.
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b>
 *    - \b    domain_boxes
 *       an array of boxes representing the index space for the entire domain
 *       on the coarsest mesh level; i.e., level zero.
 *
 *    - \b    periodic_dimension
 *       An array of integer values (expected number of values is equal to
 *       the spatial dimension of the mesh) representing the directions in
 *       which the physical domain is periodic.  A non-zero value indicates
 *       that the direction is periodic.  A zero value indicates that the
 *       direction is not periodic.
 *
 *    - \b    x_lo
 *       values representing the spatial coordinates of the lower corner of the
 *       physical domain.
 *
 *    - \b    x_up
 *       values representing the spatial coordinates of the upper corner of the
 *       physical domain.
 *
 * No values read in from a restart database may be overridden by input
 * database values.
 *
 * <b> Details: </b> <br>
 * <table>
 *   <tr>
 *     <th>parameter</th>
 *     <th>type</th>
 *     <th>default</th>
 *     <th>range</th>
 *     <th>opt/req</th>
 *     <th>behavior on restart</th>
 *   </tr>
 *   <tr>
 *     <td>domain_boxes</td>
 *     <td>array of DatabaseBoxes</td>
 *     <td>none</td>
 *     <td>all Boxes must be non-empty</td>
 *     <td>req</td>
 *     <td>May not be modified by input db on restart</td>
 *   </tr>
 *   <tr>
 *     <td>periodic_dimension</td>
 *     <td>int[]</td>
 *     <td>all values 0</td>
 *     <td>any int</td>
 *     <td>opt</td>
 *     <td>May not be modified by input db on restart</td>
 *   </tr>
 *   <tr>
 *     <td>x_lo</td>
 *     <td>double[]</td>
 *     <td>none</td>
 *     <td>each x_lo < corresponding x_up</td>
 *     <td>req</td>
 *     <td>May not be modified by input db on restart</td>
 *   </tr>
 *   <tr>
 *     <td>x_up</td>
 *     <td>double[]</td>
 *     <td>none</td>
 *     <td>each x_up > corresponding x_lo</td>
 *     <td>req</td>
 *     <td>May not be modified by input db on restart</td>
 *   </tr>
 * </table>
 *
 * A sample input file for a two-dimensional problem might look like:
 *
 * @code
 *    domain_boxes = [(0,0) , (49,39)]
 *    x_lo = 0.0 , 0.0
 *    x_up = 50.0 , 40.0
 *    periodic_dimension = 0, 1  // periodic in y only
 * @endcode
 *
 * This generates a two-dimensional rectangular domain periodic in the
 * y-direction, and having 50 cells in the x-direction and 40 cells in
 * the y-direction, with the cell size 1 unit in each direction.
 *
 * @see GridGeometry
 */

class CartesianGridGeometry:
   public GridGeometry
{
   typedef hier::PatchGeometry::TwoDimBool TwoDimBool;

public:
   /**
    * Constructor for CartesianGridGeometry initializes data
    * members based on parameters read from the specified input database
    * or from the restart database corresponding to the specified
    * object name.
    *
    * @pre !object_name.empty()
    * @pre input_db
    */
   CartesianGridGeometry(
      const tbox::Dimension& dim,
      const std::string& object_name,
      const std::shared_ptr<tbox::Database>& input_db);

   /**
    * Constructor for CartesianGridGeometry sets data members
    * based on arguments.
    *
    * @pre !object_name.empty()
    * @pre domain.size() > 0
    * @pre x_lo != 0
    * @pre x_up != 0
    */
   CartesianGridGeometry(
      const std::string& object_name,
      const double* x_lo,
      const double* x_up,
      hier::BoxContainer& domain);

   /*!
    * @brief Construct a new coarsened/refined CartesianGridGeometry object
    * with the supplied domain.
    *
    * This method is intended to be called only by std::make_shared from the
    * make[Coarsened, Refined]GridGeometry methods to make a coarsened or
    * refined version of a given CartesianGridGeometry.
    *
    * @param[in] object_name The same name as the uncoarsened/unrefined grid
    *            geometry.
    * @param[in] x_lo The same lower corner as the uncoarsened/unrefined grid
    *            geometry.
    * @param[in] x_up The same upper corner as the uncoarsened/unrefined grid
    *            geometry.
    * @param[in] domain The coarsened/refined domain.
    * @param[in] op_reg The same operator registry as the uncoarsened/unrefined
    *            grid geometry.
    *
    * @pre !object_name.empty()
    * @pre domain.size() > 0
    * @pre x_lo != 0
    * @pre x_up != 0
    */
   CartesianGridGeometry(
      const std::string& object_name,
      const double* x_lo,
      const double* x_up,
      hier::BoxContainer& domain,
      const std::shared_ptr<hier::TransferOperatorRegistry>& op_reg);

   /**
    * Destructor for CartesianGridGeometry deallocates
    * data describing grid geometry and unregisters the object with
    * the restart manager.
    */
   virtual ~CartesianGridGeometry();

   /**
    * Create and return a pointer to a refined version of this Cartesian grid
    * geometry object.
    *
    * @pre !fine_geom_name.empty()
    * @pre fine_geom_name != getObjectName()
    * @pre refine_ratio > hier::IntVector::getZero(getDim())
    */
   std::shared_ptr<hier::BaseGridGeometry>
   makeRefinedGridGeometry(
      const std::string& fine_geom_name,
      const hier::IntVector& refine_ratio) const;

   /**
    * Create and return a pointer to a coarsened version of this Cartesian grid
    * geometry object.
    *
    * @pre !coarse_geom_name.empty()
    * @pre coarse_geom_name != getObjectName()
    * @pre coarsen_ratio > hier::IntVector::getZero(getDim())
    */
   std::shared_ptr<hier::BaseGridGeometry>
   makeCoarsenedGridGeometry(
      const std::string& coarse_geom_name,
      const hier::IntVector& coarsen_ratio) const;

   /*
    * Compute grid data for patch and assign new geom_CartesianPatchGeometry
    * object to patch.
    *
    * @pre (getDim() == patch.getDim()) &&
    *      (getDim() == ratio_to_level_zero.getDim())
    * @pre ratio_to_level_zero != hier::IntVector::getZero(getDim())
    */
   void
   setGeometryDataOnPatch(
      hier::Patch& patch,
      const hier::IntVector& ratio_to_level_zero,
      const TwoDimBool& touches_regular_bdry) const;

   /**
    * Set data members for this CartesianGridGeometry object.
    *
    * @pre x_lo != 0
    * @pre x_up != 0
    */
   void
   setGeometryData(
      const double* x_lo,
      const double* x_up,
      const hier::BoxContainer& domain);

   /**
    * Return const pointer to dx array for reference level in hierarchy.
    */
   const double *
   getDx() const
   {
      return d_dx;
   }

   /**
    * Return const pointer to lower spatial coordinate for reference
    * level in hierarchy.
    */
   const double *
   getXLower() const
   {
      return d_x_lo;
   }

   /**
    * Return const pointer to upper spatial coordinate for reference
    * level in hierarchy.
    */
   const double *
   getXUpper() const
   {
      return d_x_up;
   }

   /**
    * Print class data representation.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

   /**
    * Writes the state of the CartesianGridGeometry object to the
    * restart database.
    *
    * @pre restart_db
    */
   virtual void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

protected:
   /*!
    * @brief Build operators appropriate for a CartesianGridGeometry.
    */
   virtual void
   buildOperators();

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int GEOM_CARTESIAN_GRID_GEOMETRY_VERSION;

   /*
    * Reads in domain_boxes, x_lo, and x_up from the input database.
    * Data is read from input only if the simulation is not from restart.
    * Otherwise, all values specified in the input database are ignored.
    *
    * @pre is_from_restart || input_db
    */
   void
   getFromInput(
      const std::shared_ptr<tbox::Database>& input_db,
      bool is_from_restart);

   /*
    * Read object state from the restart file and initialize class data
    * members.  The database from which the restart data is read is
    * determined by the object_name specified in the constructor.
    *
    * Unrecoverable Errors:
    *
    *    -The database corresponding to object_name is not found
    *     in the restart file.
    *
    *    -The class version number and restart version number do not
    *     match.
    *
    */
   void
   getFromRestart();

   double d_dx[SAMRAI::MAX_DIM_VAL];     // mesh increments for level 0.
   double d_x_lo[SAMRAI::MAX_DIM_VAL];   // spatial coordinates of lower corner
   // (i.e., box corner) of problem domain.
   double d_x_up[SAMRAI::MAX_DIM_VAL];   // spatial coordinates of upper corner
   // (i.e., box corner) of problem domain.

   hier::Box d_domain_box;           // smallest box covering coarsest level
                                     // (i.e., reference level) index space.

};

}
}

#endif
