/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Base class for geometry management in AMR hierarchy
 *
 ************************************************************************/

#ifndef included_geom_GridGeometry
#define included_geom_GridGeometry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <memory>

namespace SAMRAI {
namespace geom {

/*!
 * @brief Class GridGeometry manages the index space that determines the
 * extent of the coarse-level domain of a SAMRAI hierarchy.
 *
 * A GridGeometry object can be directly constructed with a consistent
 * state, or it may be used as a base class to derive child classes that manage
 * particular grid types (%e.g., Cartesian, cylindrical, etc.).  Direct
 * construction of a GridGeometry object is recommended for multiblock
 * problems and other problems where the physical locations of mesh
 * coordinates are managed by user code.
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b>
 *    - @b    num_blocks
 *       specifies the number of blocks in the mesh configuration.
 *
 *    - @b    domain_boxes_
 *       For each block, an array of boxes representing the index space for the
 *       entire domain within a block on the coarsest mesh level; i.e., level
 *       zero.  The key must have an integer value as a suffix
 *       (domain_boxes_0, domain_boxes_1, etc.), and there must be an entry
 *       for every block from 0 to num_blocks-1.
 *
 *    - @b    periodic_dimension
 *       An array of integer values (expected number of values is equal to
 *       the spatial dimension of the mesh) representing the directions in
 *       which the physical domain is periodic.  A non-zero value indicates
 *       that the direction is periodic.  A zero value indicates that the
 *       direction is not periodic.  This key should only be used when the
 *       number of blocks is 1 (a single block mesh), as periodic boundaries
 *       are not supported for multiblock meshes.
 *
 *    - @b    Singularity
 *       When there is a reduced or enhanced connectivity singularity, this key
 *       must be used to identify which blocks touch the singularity and the
 *       position of the singularity in relation to each block's index space.
 *       The key for this entry must include a unique trailing integer, and the
 *       integers for the full set of Singularity keys must be a continuous
 *       sequence beginning with 0.
 *
 *    - @b    BlockNeighbors
 *       For multiblock grids, a BlockNeighbors entry must be given for every
 *       pair of blocks that touch each other in any way.  Like Singularity,
 *       each entry must have a trailing integer beginning with 0.
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
 *     <td>num_blocks</td>
 *     <td>int</td>
 *     <td>1</td>
 *     <td>>=1</td>
 *     <td>opt</td>
 *     <td>May not be modified by input db on restart</td>
 *   </tr>
 *   <tr>
 *     <td>domain_boxes_N</td>
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
 *     <td>Singularity</td>
 *     <td>see Multiblock.pdf for subentries</td>
 *     <td>none</td>
 *     <td>see Multiblock.pdf for subentries</td>
 *     <td>opt</td>
 *     <td>May not be modified by input db on restart</td>
 *   </tr>
 *   <tr>
 *     <td>BlockNeighbors</td>
 *     <td>see Multiblock.pdf for subentries</td>
 *     <td>none</td>
 *     <td>see Multiblock.pdf for subentries</td>
 *     <td>opt</td>
 *     <td>May not be modified by input db on restart</td>
 *   </tr>
 * </table>
 *
 * A description of the input format for Singularity* and BlockNeighbors*
 * is included in the Multiblock.pdf document in the docs/userdocs
 * directory of the SAMRAI distribution.
 *
 * @see hier::BaseGridGeometry
 */

class GridGeometry:
   public hier::BaseGridGeometry
{
public:
   /*!
    * @brief Construct a new GridGeometry object and initialize from
    * input.
    *
    * This constructor for GridGeometry initializes data members
    * based on parameters read from the specified input database.
    *
    * This constructor is intended for use when directly constructing a
    * GridGeometry without using a derived child class.  The object will
    * contain all index space grid information for a mesh, but nothing about
    * the physical coordinates of the mesh.
    *
    * @param[in]  dim
    * @param[in]  object_name
    * @param[in]  input_db
    * @param[in]  allow_multiblock set to false if called by inherently single
    *             block derived class such as CartesianGridGeometry
    */
   GridGeometry(
      const tbox::Dimension& dim,
      const std::string& object_name,
      const std::shared_ptr<tbox::Database>& input_db,
      bool allow_multiblock = true);

   /*!
    * @brief Construct a new GridGeometry object based on arguments.
    *
    * This constructor creates a new GridGeometry object based on the
    * arguments, rather than relying on input or restart data.
    *
    * @param[in]  object_name
    * @param[in]  domain      Each element of the array describes the index
    *                         space for a block.
    */
   GridGeometry(
      const std::string& object_name,
      hier::BoxContainer& domain);

   /*!
    * @brief Construct a new coarsened/refined GridGeometry object with the
    * supplied domain.
    *
    * This method is intended to be called only by std::make_shared from the
    * make[Coarsened, Refined]GridGeometry methods to make a coarsened or
    * refined version of a given GridGeometry.
    *
    * @param[in]  object_name The same name as the uncoarsened/unrefined grid
    *                         geometry.
    * @param[in]  domain The coarsened/refined domain.
    * @param[in]  op_reg The same operator registry as the
    *                    uncoarsened/unrefined grid geometry.
    */
   GridGeometry(
      const std::string& object_name,
      hier::BoxContainer& domain,
      const std::shared_ptr<hier::TransferOperatorRegistry>& op_reg);

   /*!
    * @brief Virtual destructor
    */
   virtual ~GridGeometry();

   /*!
    * @brief Create a pointer to a refined version of this grid geometry
    *        object.
    *
    * Virtual method -- should be overridden in specialized grid geometry
    * classes
    *
    * @param[in]     fine_geom_name std::string name of the geometry object
    * @param[in]     refine_ratio the refinement ratio.
    *
    * @return The pointer to the grid geometry object.
    *
    * @pre !fine_geom_name.empty()
    * @pre fine_geom_name != getObjectName()
    * @pre getDim() == refine_ratio.getDim()
    * @pre refine_ratio > hier::IntVector::getZero(getDim())
    */
   virtual std::shared_ptr<hier::BaseGridGeometry>
   makeRefinedGridGeometry(
      const std::string& fine_geom_name,
      const hier::IntVector& refine_ratio) const;

   /*!
    * @brief Create a pointer to a coarsened version of this grid geometry
    *        object.
    *
    * Virtual method -- should be overridden in specialized grid geometry
    * classes
    *
    * @param[in]     coarse_geom_name std::string name of the geometry object
    * @param[in]     coarsen_ratio the coasening ratio
    *
    * @return The pointer to a coarsened version of this grid geometry object.
    *
    * @pre !coarse_geom_name.empty()
    * @pre coarse_geom_name != getObjectName()
    * @pre getDim() == coarsen_ratio.getDim()
    * @pre coarsen_ratio > hier::IntVector::getZero(getDim())
    */
   virtual std::shared_ptr<hier::BaseGridGeometry>
   makeCoarsenedGridGeometry(
      const std::string& coarse_geom_name,
      const hier::IntVector& coarsen_ratio) const;

protected:
   /*!
    * @brief Build operators appropriate for a GridGeometry.
    */
   virtual void
   buildOperators();
};

}
}

#endif
