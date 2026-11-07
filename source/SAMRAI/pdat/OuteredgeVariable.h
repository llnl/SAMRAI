/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Variable class for defining outeredge centered variables
 *
 ************************************************************************/

#ifndef included_pdat_OuteredgeVariable
#define included_pdat_OuteredgeVariable

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/Complex.h"

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class OuteredgeVariable<TYPE> is a templated variable class
 * used to define edge-centered data quantities only on patch boundaries.
 * It is a subclass of hier::Variable and is templated on the type
 * of the underlying data (e.g., double, int, bool, etc.).
 *
 * Note that the data layout in the outeredge data arrays matches the corresponding
 * array sections provided by the edge data implementation.  See header file for
 * the OuteredgeData<TYPE> class for a more detailed description of the data layout.
 *
 * @see EdgeData
 * @see OuteredgeData
 * @see OuteredgeDataFactory
 * @see hier::Variable
 */

template<class TYPE>
class OuteredgeVariable:public hier::Variable
{
public:
   /*!
    * @brief Create an outeredge variable object having properties
    * specified by the name and depth (i.e., number of data values
    * at each index location).  The default depth is one.
    *
    * Note that The ghost cell width for all outeredge data is currently
    * fixed at zero; this may be changed in the future if needed.
    */
   OuteredgeVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      int depth = 1);

   /*!
    * @brief Virtual destructor for outeredge variable objects.
    */
   virtual ~OuteredgeVariable();

   /*!
    * @brief Return a boolean true value indicating that fine patch
    * values take precedence on coarse-fine interfaces.
    */
   bool fineBoundaryRepresentsVariable() const {
      return true;
   }

   /*!
    * @brief Return true indicating that outeredge data
    * exists on the patch boundary.
    */
   bool dataLivesOnPatchBorder() const {
      return true;
   }

   /*!
    * @brief Return the the depth (number of components).
    */
   int
   getDepth() const;

private:
   // neither of the following functions are implemented
   OuteredgeVariable(
      const OuteredgeVariable&);
   OuteredgeVariable&
   operator = (
      const OuteredgeVariable&);

};

}
}

#include "SAMRAI/pdat/OuteredgeVariable.cpp"

#endif
