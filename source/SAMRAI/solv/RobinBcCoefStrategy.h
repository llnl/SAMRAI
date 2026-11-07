/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Robin boundary condition problem-dependent interfaces
 *
 ************************************************************************/
#ifndef included_solv_RobinBcCoefStrategy
#define included_solv_RobinBcCoefStrategy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Patch.h"

#include <memory>

namespace SAMRAI {
namespace solv {

/*!
 * @brief Interface for specifying Robin boundary conditions.
 *
 * The Robin boundary conditions are specified in terms of
 * the coefficients @f$ \alpha @f$, @f$ \beta @f$ and @f$ \gamma @f$
 * in the Robin formula
 * @f[  \alpha u + \beta u_n = \gamma @f]
 * applied on the boundary with outward normal n.
 *
 * This class specifies the interfaces for communicating the
 * boundary condition coefficients.
 */
class RobinBcCoefStrategy
{

public:
   /*!
    * @brief Constructor
    */
   RobinBcCoefStrategy();

   /*!
    * @brief Destructor.
    */
   virtual ~RobinBcCoefStrategy();

   //@{

   /*!
    * @name Functions to set boundary condition coefficients a and g.
    */

   /*!
    * @brief User-supplied function to fill arrays of Robin boundary
    * condition coefficients at a patch boundary.
    *
    * This class specifies the Robin boundary condition coefficients
    * at discrete locations on the patch boundary.
    * Though these locations are defined by boundary box object,
    * they do not necessarily coincide with
    * the centers of the cells referred to by those boxes.
    * These locations typically coincide with the nodes
    * or face centers which do lie on the patch boundary.
    * Accordingly, you use this function to provide the
    * boundary coefficients at those locations by filling an array
    * at indices corresponding to those locations.
    *
    * When setting the values of the boundary condition coefficients
    * it is useful to note that for any cell (i,j,k),
    * the indices of the sides, edges and nodes are easily determined.
    * The index on the lower
    * side of the cell is the same as the index of the cell, whereas
    * the index on the upper side of the cell has the next higher
    * value.  In 2D, the cell and its surrounding nodes and faces
    * has the following indices:
    * @verbatim
    *
    *       (i,j+1)----(i,j+1)---(i+1,j+1)
    *          |                     |
    *          |                     |
    *          |                     |
    *          |                     |
    *        (i,j)      (i,j)     (i+1,j)
    *          |                     |
    *          |                     |
    *          |                     |
    *          |                     |
    *        (i,j)------(i,j)-----(i+1,j)
    *
    * @endverbatim
    * Once this is understood, translation between the index
    * in the boundary box index space to the index of things
    * on the boundary is simple.
    *
    * The boundary condition coefficients should be placed
    * in the pdat::ArrayData<TYPE> objects, @c acoef_data and @c gcoef_data
    * (see argument list), which are dimensioned to contain the indices
    * of the points alligned with @c variable and lying on the
    * the boundary defined by @c bdry_box.
    *
    * This function is only used with type-1 boundary boxes,
    * such as faces in 3D.
    * Other types of boundaries do not have a well-defined
    * surface normal.
    *
    * The parameter @c variable is passed through to tell
    * the implementation of this function what variable
    * to set the coefficients for.  You may wish to ignore
    * it if your implementation is intended for a specific
    * variable.
    *
    * @param acoef_data boundary coefficient data.
    *        The array will have been defined to include index range
    *        for corresponding to the boundary box @c bdry_box
    *        and appropriate for the alignment of the given variable.
    *        If this is a null pointer, then the calling function
    *        is not interested in a, and you can disregard it.
    * @param bcoef_data boundary coefficient data.
    *        This array is exactly like @c acoef_data,
    *        except that it is to be filled with the b coefficient.
    * @param gcoef_data boundary coefficient data.
    *        This array is exactly like @c acoef_data,
    *        except that it is to be filled with the g coefficient.
    * @param variable variable to set the coefficients for.
    *        If implemented for multiple variables, this parameter
    *        can be used to determine which variable's coefficients
    *        are being sought.
    * @param patch patch requiring bc coefficients
    * @param bdry_box boundary box showing where on the boundary
    *        the coefficient data is needed.
    * @param fill_time solution time corresponding to filling, for use
    *        when coefficients are time-dependent.
    */
   virtual void
   setBcCoefs(
      const std::shared_ptr<pdat::ArrayData<double> >& acoef_data,
      const std::shared_ptr<pdat::ArrayData<double> >& bcoef_data,
      const std::shared_ptr<pdat::ArrayData<double> >& gcoef_data,
      const std::shared_ptr<hier::Variable>& variable,
      const hier::Patch& patch,
      const hier::BoundaryBox& bdry_box,
      const double fill_time = 0.0) const = 0;

   /*
    * @brief Return how many cells past the edge or corner of the
    * patch the object can fill.
    *
    * The "extension" used here is the number of cells that
    * a boundary box extends past the patch in the direction
    * parallel to the boundary.
    *
    * Note that the inability to fill the sufficient
    * number of cells past the edge or corner of the patch
    * may preclude the child class from being used in data
    * refinement operations that require the extra data,
    * such as linear refinement.
    *
    * The boundary box that setBcCoefs() is required to fill
    * should not extend past the limits returned by this
    * function.
    */
   virtual hier::IntVector
   numberOfExtensionsFillable() const = 0;

   //@}

};

}
}

#endif  // included_solv_RobinBcCoefStrategy
