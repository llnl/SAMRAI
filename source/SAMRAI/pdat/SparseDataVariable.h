/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   SparseDataVariable
 *
 ************************************************************************/
#ifndef included_pdat_SparseDataVariable
#define included_pdat_SparseDataVariable

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Variable.h"

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Variable class used to define sparse data on a set of box indices
 * defined by a template parameter defining the box geometry.
 *
 * In reality, this class can only be used for cell-centered and node-centered
 * box geoemtries because of the internal calling semantics of the underlying
 * overlap operations.
 *
 * @see SparseData
 * @see SparseDataFactory
 * @see Variable
 */
template<typename BOX_GEOMETRY>
class SparseDataVariable:public hier::Variable
{
public:
   /*!
    * @brief Create a sparse data variable object with the specified name.
    *
    * The creation of the variable creates a "default" SparseDataFactory
    * with ghost width set to zero.
    *
    * @param [in] dim
    * @param [in] name
    */
   SparseDataVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      const std::vector<std::string>& dbl_attributes,
      const std::vector<std::string>& int_attributes);

   /*!
    * @brief Destructor
    */
   ~SparseDataVariable();

   /*!
    * @brief Returns true since this class can be used for either cell-centered
    *        or node-centered index spaces and this covers all cases.
    */
   bool
   fineBoundaryRepresentsVariable() const {
      return true;
   }

   /*!
    * @brief Returns true since this class can be used for either cell-centered
    *        or node-centered index spaces and this covers all cases.
    */
   bool
   dataLivesOnPatchBorder() const {
      return true;
   }

private:
   /*
    * copy c'tor and assignment operator are private to prevent the
    * compiler from generating a default.
    */
   SparseDataVariable(
      const SparseDataVariable& rhs);

   SparseDataVariable&
   operator = (
      const SparseDataVariable& rhs);

}; // end class SparseDataVariable.

} // end namespace pdat.
} // end namespace SAMRAI

#include "SAMRAI/pdat/SparseDataVariable.cpp"

#endif
