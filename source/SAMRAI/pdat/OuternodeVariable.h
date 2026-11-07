/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Variable class for defining outernode centered variables
 *
 ************************************************************************/

#ifndef included_pdat_OuternodeVariable
#define included_pdat_OuternodeVariable

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/hier/Variable.h"

#include <string>

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class OuternodeVariable<TYPE> is a templated variable class
 * used to define node-centered data quantities only on patch boundaries.
 * It is a subclass of hier::Variable and is templated on the type
 * of the underlying data (e.g., double, int, bool, etc.).
 *
 * Note that the data layout in the outernode data arrays matches the corresponding
 * array sections provided by the node data implementation.  See header file for
 * the OuternodeData<TYPE> class for a more detailed description of the data layout.
 *
 * @see NodeData<TYPE>
 * @see OuternodeData<TYPE>
 * @see OuternodeDataFactory<TYPE>
 * @see hier::Variable
 */

template<class TYPE>
class OuternodeVariable:public hier::Variable
{
public:
   /*!
    * @brief Create an outernode variable object having properties
    * specified by the name and depth (i.e., number of data values
    * at each index location).  The default depth is one.
    *
    * Note that The ghost cell width for all outernode data is currently
    * fixed at zero; this may be changed in the future if needed.
    */
   OuternodeVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      int depth = 1);

   /*!
    * @brief Virtual destructor for outernode variable objects.
    */
   virtual ~OuternodeVariable();

   /*!
    * @brief Return a boolean true value indicating that fine patch
    * values take precedence on coarse-fine interfaces.
    */
   bool fineBoundaryRepresentsVariable() const {
      return true;
   }

   /*!
    * @brief Return true indicating that outernode data
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
   OuternodeVariable(
      const OuternodeVariable&);
   OuternodeVariable&
   operator = (
      const OuternodeVariable&);

};

}
}

#include "SAMRAI/pdat/OuternodeVariable.cpp"

#endif
