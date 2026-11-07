/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_CellVariable
#define included_pdat_CellVariable

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#include <string>

namespace SAMRAI {
namespace pdat {

/*!
 * Class CellVariable<TYPE> is a templated variable class used to define
 * cell-centered quantities on an AMR mesh.   It is a subclass of
 * hier::Variable and is templated on the type of the underlying data
 * (e.g., double, int, bool, etc.).
 *
 * See header file for CellData<TYPE> class for a more detailed
 * description of the data layout.
 *
 * @see CellData
 * @see CellDataFactory
 * @see hier::Variable
 */

template<class TYPE>
class CellVariable:public hier::Variable
{
public:
   /*!
    * @brief Create a cell-centered variable object with the given name and
    * depth (i.e., number of data values at each cell index location).
    * A default depth of one is provided.
    */
   CellVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      int depth = 1);

   /*!
    * @brief Create a cell-centered variable object with the given name,
    * allocator, and depth (i.e., number of data values at each cell index
    * location).  A default depth of one is provided.
    */
   CellVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      tbox::ResourceAllocator allocator,
      int depth = 1);

   /*!
    * @brief Virtual destructor for cell variable objects.
    */
   virtual ~CellVariable();

   /*!
    * @brief Return true indicating that cell data quantities will always
    * be treated as though fine values take precedence on coarse-fine
    * interfaces.  Note that this is really artificial since the cell
    * data index space matches the cell-centered index space for AMR
    * patches.  However, some value must be supplied for communication
    * operations.
    */
   bool
   fineBoundaryRepresentsVariable() const;

   /*!
    * @brief Return false indicating that cell data on a patch interior
    * does not exist on the patch boundary.
    */
   bool
   dataLivesOnPatchBorder() const;

   /*!
    * @brief Return the the depth (number of components).
    */
   int
   getDepth() const;

private:
   // Unimplemented copy constructor
   CellVariable(
      const CellVariable&);

   // Unimplemented assignment operator
   CellVariable&
   operator = (
      const CellVariable&);
};

}
}

#include "SAMRAI/pdat/CellVariable.cpp"

#endif
