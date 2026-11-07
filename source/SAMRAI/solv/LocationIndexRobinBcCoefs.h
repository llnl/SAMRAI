/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Robin boundary condition problem-dependent interfaces
 *
 ************************************************************************/
#ifndef included_solv_LocationIndexRobinBcCoefs
#define included_solv_LocationIndexRobinBcCoefs

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/solv/RobinBcCoefStrategy.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Utilities.h"

#include <memory>

namespace SAMRAI {
namespace solv {

/*!
 * @brief A prefabricated Robin boundary condition coefficients
 * for coefficients that are entirely specified by the boundary
 * box location index.
 *
 * This implementation of the strategy
 * class RobinBcCoefStrategy may be used when your
 * Robin boundary condition coefficients are completely determined
 * by the location index of the boundary box.
 *
 * Before this class is used in to provide the boundary condition
 * coefficients, you must specify what boundary conditions to
 * associate with what location index.  Methods for specifying
 * these are setBoundaryValue(), setBoundarySlope() and
 * setRawCoefficients().  The first two are for Dirichlet
 * and Neumann boundary conditions, respectively.  If the boundary
 * condition is the more general Robin boundary condition,
 * the third function should be used to set the coefficients
 * a and g directly (see RobinBcCoefStrategy) for the
 * meanings of a and g.
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions:</b>
 *    - \b    boundary_N
 *       You can specify the boundary conditions for any location index
 *       through the input database.  One line is required for each
 *       location index.  The input parameters are "boundary_N", where
 *       N is the index of the location.  Each parameter must be
 *       a std::string array so that all boundary types can be accomodated
 *       the same way.  The first std::string must be one of "value",
 *       "slope" or "coefficients".  If the std::string is "value" or "slope"
 *       the next std::string is the value you want to set.  If the first
 *       std::string is "coefficients", the next two strings specify the
 *       values of a and b.
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
 *     <td>boundary_N</td>
 *     <td>array of strings</td>
 *     <td>none</td>
 *     <td>"value", "gval"; "slope", "gval"; "coefficients", "aVal", "bVal"</td>
 *     <td>req</td>
 *     <td>Not written to restart.  Value in input db used.</td>
 * </table>
 *
 */
class LocationIndexRobinBcCoefs:
   public RobinBcCoefStrategy
{

public:
   /*!
    * @brief Constructor using input database.
    *
    * @pre input_db
    */
   LocationIndexRobinBcCoefs(
      const tbox::Dimension& dim,
      const std::string& object_name,
      const std::shared_ptr<tbox::Database>& input_db =
         std::shared_ptr<tbox::Database>());

   /*!
    * @brief Destructor.
    */
   virtual ~LocationIndexRobinBcCoefs();

   /*!
    * @brief Function to fill arrays of Robin boundary
    * condition coefficients at a patch boundary.
    *
    * This implementation of the virtual function
    * RobinBcCoefStrategy::setBcCoefs()
    * fills the coefficient arrays with constant values
    * set according to the location index of the boundary box.
    *
    * @param acoef_data boundary coefficient data.
    *        This is defined to include index range for
    *        the boundary faces on the boundary box @c bdry_box.
    *        If this is a null pointer, then the calling function
    *        is not interested in a, and you can disregard it.
    * @param bcoef_data boundary coefficient data.
    *        This is defined to include index range for
    *        the boundary faces on the boundary box @c bdry_box.
    * @param gcoef_data boundary coefficient data.
    *        This is defined to include index range for
    *        the boundary faces on the boundary box @c bdry_box.
    * @param variable variable to set the coefficients for.
    * @param patch patch requiring bc coefficients
    * @param bdry_box boundary box showing where on the boundary
    *        the coefficient data is needed.
    * @param fill_time Solution time corresponding to filling,
    *        for use when coefficients are time-dependent.
    *
    * @pre (d_dim == patch.getDim()) && (d_dim == bdry_box.getDim())
    * @pre (!acoef_data || (d_dim == acoef_data.getDim())
    * @pre (!bcoef_data || (d_dim == acoef_data.betDim())
    * @pre (!gcoef_data || (d_dim == acoef_data.getDim())
    * @pre (bdry_box.getLocationIndex() >= 0) &&
    *      (bdry_box.getLocationIndex() < 2 * d_dim.getValue())
    */
   void
   setBcCoefs(
      const std::shared_ptr<pdat::ArrayData<double> >& acoef_data,
      const std::shared_ptr<pdat::ArrayData<double> >& bcoef_data,
      const std::shared_ptr<pdat::ArrayData<double> >& gcoef_data,
      const std::shared_ptr<hier::Variable>& variable,
      const hier::Patch& patch,
      const hier::BoundaryBox& bdry_box,
      double fill_time = 0.0) const;

   /*
    * @brief Return how many cells past the edge or corner of the
    * patch the object can fill.
    */
   hier::IntVector
   numberOfExtensionsFillable() const;

   /*!
    * @brief Set the boundary value at a given location index.
    *
    * @param location_index Set coefficients for this index.
    * @param value Boundary value at @c location_index.
    *
    * @pre (location_index >= 0) && (location_index < 2 * d_dim.getValue())
    */
   void
   setBoundaryValue(
      int location_index,
      double value)
   {
      if (location_index < 0 || location_index >= 2 * d_dim.getValue()) {
         TBOX_ERROR("Location index in " << d_dim.getValue() << "D must be\n"
                                         << "in [0," << 2 * d_dim.getValue() - 1 << ").\n");
      }
      d_a_map[location_index] = 1.0;
      d_b_map[location_index] = 0.0;
      d_g_map[location_index] = value;
   }

   /*!
    * @brief Set the boundary slope at a given location index.
    *
    * @param location_index Set coefficients for this index.
    * @param slope Boundary slope at @c location_index.
    *
    * @pre (location_index >= 0) && (location_index < 2 * d_dim.getValue())
    */
   void
   setBoundarySlope(
      int location_index,
      double slope)
   {
      if (location_index < 0 || location_index >= 2 * d_dim.getValue()) {
         TBOX_ERROR("Location index in " << d_dim.getValue() << "D must be\n"
                                         << "in [0," << 2 * d_dim.getValue() - 1 << ").\n");
      }
      d_a_map[location_index] = 0.0;
      d_b_map[location_index] = 1.0;
      d_g_map[location_index] = slope;
   }

   /*!
    * @brief Set the values of coefficients a and g at a
    * given location index.
    *
    * See RobinBcCoefStrategy for the definitions
    * of coefficients a and g.
    *
    * If the boundary condition is neither Dirichlet nor
    * Neumann (a general Robin boundary condition), use
    * this function to set the values of the bc coefficients.
    *
    * @param location_index Set coefficients for this index.
    * @param a Value of coefficient a at given location index.
    * @param b Value of coefficient b at given location index.
    * @param g Value of coefficient g at given location index.
    *
    * @pre (location_index >= 0) && (location_index < 2 * d_dim.getValue())
    */
   void
   setRawCoefficients(
      int location_index,
      double a,
      double b,
      double g)
   {
      if (location_index < 0 || location_index >= 2 * d_dim.getValue()) {
         TBOX_ERROR("Location index in " << d_dim.getValue() << "D must be\n"
                                         << "in [0," << 2 * d_dim.getValue() - 1 << ").\n");
      }
      d_a_map[location_index] = a;
      d_b_map[location_index] = b;
      d_g_map[location_index] = g;
   }

   /*!
    * @brief Access coefficients.
    *
    * @param location_index Set coefficients for this index.
    * @param a Value of coefficient a at given location index.
    * @param b Value of coefficient b at given location index.
    * @param g Value of coefficient g at given location index.
    *
    * @pre (location_index >= 0) && (location_index < 2 * d_dim.getValue())
    */
   void
   getCoefficients(
      int location_index,
      double& a,
      double& b,
      double& g) const
   {
      if (location_index < 0 || location_index >= 2 * d_dim.getValue()) {
         TBOX_ERROR("Location index in " << d_dim.getValue() << "D must be\n"
                                         << "in [0," << 2 * d_dim.getValue() - 1 << ").\n");
      }
      a = d_a_map[location_index];
      b = d_b_map[location_index];
      g = d_g_map[location_index];
   }

   /**
    * @brief Get the name of this object.
    *
    * @return The name of this object.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

   /*!
    * @brief Assignment operator.
    */
   LocationIndexRobinBcCoefs&
   operator = (
      const LocationIndexRobinBcCoefs& r);

private:
   /*
    * @brief Set state from input database.
    *
    * See the class description for the parameters that can be set
    * from a database.
    *
    * @param input_db Input database.
    *
    * @pre input_db
    */
   void
   getFromInput(
      const std::shared_ptr<tbox::Database>& input_db);

   /*
    * @brief Object dimension
    */
   const tbox::Dimension d_dim;

   /*
    * @brief Object name.
    */
   std::string d_object_name;

   /*
    * @brief Mapping for a coefficient.
    */
   double d_a_map[2 * SAMRAI::MAX_DIM_VAL];

   /*
    * @brief Mapping for b coefficient.
    */
   double d_b_map[2 * SAMRAI::MAX_DIM_VAL];
   /*
    * @brief Mapping for g coefficient.
    */
   double d_g_map[2 * SAMRAI::MAX_DIM_VAL];

};

}
}

#endif
