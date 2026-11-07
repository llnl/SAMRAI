/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_OutersideVariable
#define included_pdat_OutersideVariable

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#include <string>

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class OutersideVariable<TYPE> is a templated variable class
 * used to define side-centered data quantities only on patch boundaries.
 * It is a subclass of hier::Variable and is templated on the type
 * of the underlying data (e.g., double, int, bool, etc.).
 *
 * Note that the data layout in the outerside data arrays matches the corresponding
 * array sections provided by the side data implementation. See header file for
 * the OutersideData<TYPE> class for a more detailed description of the data layout.
 *
 * IMPORTANT: The class OuterfaceVariable<TYPE> and associated "outerface
 * data" classes define the same storage as this outerside variable class,
 * except that the individual array indices are permuted in the outerface
 * data type.
 *
 * @see SideData
 * @see OutersideData
 * @see OutersideDataFactory
 * @see hier::Variable
 */

template<class TYPE>
class OutersideVariable:public hier::Variable
{
public:
   /*!
    * @brief Create an outerside variable object having properties
    * specified by the name and depth (i.e., number of data values
    * at each index location).  The default depth is one.
    *
    * Note that The ghost cell width for all outerside data is currently
    * fixed at zero; this may be changed in the future if needed.
    */
   OutersideVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      int depth = 1);

   /*!
    * @brief Constructor that also includes an Umpire allocator for
    * allocations of the underlying data.
    */
   OutersideVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      tbox::ResourceAllocator allocator,
      int depth = 1);

   /*!
    * @brief Virtual destructor for outerside variable objects.
    */
   virtual ~OutersideVariable();

   /*!
    * @brief Return a boolean true value indicating that fine patch
    * values take precedence on coarse-fine interfaces.
    */
   bool fineBoundaryRepresentsVariable() const {
      return true;
   }

   /*!
    * @brief Return true indicating that outerside data
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
   OutersideVariable(
      const OutersideVariable&);
   OutersideVariable&
   operator = (
      const OutersideVariable&);

};

}
}

#include "SAMRAI/pdat/OutersideVariable.cpp"

#endif
