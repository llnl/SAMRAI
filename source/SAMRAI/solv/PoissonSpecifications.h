/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Specifications for the scalar Poisson equation
 *
 ************************************************************************/
#ifndef included_solv_PoissonSpecifications
#define included_solv_PoissonSpecifications

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Utilities.h"

#include <string>

namespace SAMRAI {
namespace solv {

/*!
 * @brief Light class holding specifications for cell-centered
 * implementation of the scalar Poisson equation.
 *
 * The scalar Poisson equation is
 * @f$ \nabla ( D \nabla u ) + C u = f @f$,
 * where C is a scalar field, D is the diffusion coefficient.
 * and u and f are scalar quantities.
 *
 * This class describes the things you can set: C, D.
 *
 * Note that the storage and alignment of u, f, C and D depend
 * on the implementation of the solver.  For example, if the
 * solver is cell centered, u, f and C are cell-centered while
 * D is side-centered.
 */

class PoissonSpecifications
{
public:
   /*!
    * @brief Constructor.
    *
    * Sets the specifications to their default state:
    * - C is zero
    * - D is uniformly 1
    *
    * @param object_name Name of object.
    */
   explicit PoissonSpecifications(
      const std::string& object_name);

   /*!
    * @brief Copy constructor.
    */
   PoissonSpecifications(
      const std::string& object_name,
      const PoissonSpecifications& r);

   /*!
    * @brief Destructor (does nothing).
    */
   virtual ~PoissonSpecifications();

   /*!
    * @brief Assignment operator
    *
    * Assign everything except name.
    */
   PoissonSpecifications&
   operator = (
      const PoissonSpecifications& r)
   {
      d_D_id = r.d_D_id;
      d_D_constant = r.d_D_constant;
      d_C_zero = r.d_C_zero;
      d_C_id = r.d_C_id;
      d_C_constant = r.d_C_constant;
      return *this;
   }

   /*!
    * @brief Print out class data.
    */
   virtual void
   printClassData(
      std::ostream& stream) const;

   //@{
   //! @name Functions for setting and getting D

   /*!
    * @brief Set the patch data index for variable D.
    *
    * In addition, disregard any previous value
    * specified by setDConstant().
    *
    * @pre id >= 0
    */
   void
   setDPatchDataId(
      int id)
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (id < 0) {
         TBOX_ERROR(d_object_name << ": Invalid patch data id.\n");
      }
#endif
      d_D_id = id;
      d_D_constant = 0.0;
   }

   /*!
    * @brief Set the constant value variable D.
    *
    * In addition, disregard any previous patch data index
    * specified by setDPatchDataId().
    */
   void
   setDConstant(
      double constant)
   {
      d_D_id = -1;
      d_D_constant = constant;
   }

   /*!
    * @brief Whether D is variable (described by a patch data id).
    *
    * @return True if D is variable, described by the patch data
    *         id given in setCPatchDataId().
    */
   bool
   dIsVariable() const
   {
      return d_D_id != -1;
   }

   /*!
    * @brief Whether D is constant.
    *
    * @return True if D is constant, as specified by setCConstant().
    */
   bool
   dIsConstant() const
   {
      return d_D_id == -1;
   }

   /*!
    * @brief Get D's patch data id
    *
    * Error if D is not represented by a patch data id.
    *
    * @return D's id
    *
    * @pre d_D_id != -1
    */
   int
   getDPatchDataId() const
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_D_id == -1) {
         TBOX_ERROR(d_object_name << ": D not prepresented by a patch data.\n");
      }
#endif
      return d_D_id;
   }

   /*!
    * @brief Get D constant value
    *
    * Error if D is not represented by a constant.
    *
    * @return D's constant value
    *
    * @pre d_D_id != -1
    */
   double
   getDConstant() const
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_D_id != -1) {
         TBOX_ERROR(d_object_name << ": D not prepresented by a constant.\n");
      }
#endif
      return d_D_constant;
   }

   //@}

   //@{
   //! @name Functions for setting and getting C

   /*!
    * @brief Set the patch data index for C.
    *
    * In addition, disregard any previous values
    * specified by setCConstant() or setCZero().
    *
    * @pre id >= 0
    */
   void
   setCPatchDataId(
      int id)
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (id < 0) {
         TBOX_ERROR(d_object_name << ": Invalid patch data id.\n");
      }
#endif
      d_C_zero = false;
      d_C_id = id;
      d_C_constant = 0.0;
   }

   /*!
    * @brief Set C to a constant.
    *
    * In addition, disregard any previous value
    * specified by setCPatchDataId() or setCZero().
    *
    * If you want to set C to zero, use setCZero() instead.
    * This allows solvers to take advantage of fact C is absent.
    */
   void
   setCConstant(
      double constant)
   {
      d_C_zero = false;
      d_C_id = -1;
      d_C_constant = constant;
   }

   /*!
    * @brief Set the value of C to zero.
    *
    * In addition, disregard any previous patch data index
    * specified by setCPatchDataId() and any previous constant
    * specified by setCConstant().
    */
   void
   setCZero()
   {
      d_C_zero = true;
      d_C_id = -1;
      d_C_constant = 0.0;
   }

   /*!
    * @brief Whether C is variable (described by a patch data id).
    *
    * @return True if C is variable, described by the patch data
    *         id given in setCPatchDataId().
    */
   bool
   cIsVariable() const
   {
      return d_C_id != -1;
   }

   /*!
    * @brief Whether C is zero.
    *
    * As it pertains to what this function returns,
    * C is zero @em only by calling setCZero().
    * Calling setCConstant() does @em not make C zero,
    * even if you pass in the value of zero.
    *
    * @return True if C is exactly zero, as set by setCZero().
    */
   bool
   cIsZero() const
   {
      return d_C_zero;
   }

   /*!
    * @brief Whether C is constant.
    *
    * As it pertains to what this function returns,
    * C is constant @em only by calling setCConstant().
    * Calling setCZero() does @em not make C a constant.
    *
    * @return True if C is constant, as specified by setCConstant().
    */
   bool
   cIsConstant() const
   {
      return !d_C_zero && (d_C_id == -1);
   }

   /*!
    * @brief Get C's patch data id
    *
    * Error if C is not represented by a patch data id.
    *
    * @return C's patch data id
    *
    * @pre d_C_id != -1
    */
   int
   getCPatchDataId() const
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_C_id == -1) {
         TBOX_ERROR(d_object_name << ": C not prepresented by a an index.\n");
      }
#endif
      return d_C_id;
   }

   /*!
    * @brief Get C as a constant value.
    *
    * Error if C is not represented by a constant.
    *
    * @return C's constant value
    *
    * @pre (d_C_id == -1) && !d_C_zero
    */
   double
   getCConstant() const
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_C_id != -1 || d_C_zero) {
         TBOX_ERROR(d_object_name << ": C is not prepresented by a constant.\n");
      }
#endif
      return d_C_constant;
   }

   //@}

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

private:
   /*!
    * @brief Object name.
    */
   std::string d_object_name;

   int d_D_id;
   double d_D_constant;

   bool d_C_zero;
   int d_C_id;
   double d_C_constant;

};

} // namespace SAMRAI
}

#endif
