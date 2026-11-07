/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   IndexVariable
 *
 ************************************************************************/

#ifndef included_pdat_IndexVariable
#define included_pdat_IndexVariable

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Variable.h"

#include <string>

namespace SAMRAI {
namespace pdat {

/**
 * Class IndexVariable<TYPE,BOX_GEOMETRY> is a templated variable
 * class used to define quantities that exist on an irregular
 * cell-centered index set.  The template parameter TYPE defines the
 * storage at each index location.  For example, this class is used to
 * represent embedded boundary features as a regular patch data type
 * using the BoundaryCell class as the template type.  The template
 * parameter BOX_GEOMETRY allows IndexVariables to be instantiated
 * with a provided centering and geometry in index space via a
 * BoxGeometry (e.g. CellGeometry, NodeGeometry).
 *
 * Please consult the README file in the index data source directory for
 * instructions on using this class to provide other irregular index set
 * types.
 *
 * @see IndexData
 * @see IndexDataFactory
 * @see Variable
 */

template<class TYPE, class BOX_GEOMETRY>
class IndexVariable:public hier::Variable
{
public:
   /**
    * Create an index variable object with the specified name.
    */
   IndexVariable(
      const tbox::Dimension& dim,
      const std::string& name);

   /**
    * Virtual destructor for index variable objects.
    */
   virtual ~IndexVariable();

   /**
    * Return true so that the index data quantities will always be treated as cell-
    * centered quantities as far as communication is concerned.  Note that this is
    * really artificial since the cell data index space matches the cell-centered
    * index space for AMR patches.  Thus, cell data does not live on patch borders
    * and so there is no ambiguity reagrding coarse-fine interface values.
    */
   bool fineBoundaryRepresentsVariable() const {
      return true;
   }

   /**
    * Return false since the index data index space matches the cell-centered
    * index space for AMR patches.  Thus, index data does not live on patch borders.
    */
   bool dataLivesOnPatchBorder() const {
      return false;
   }

private:
   // Unimplemented copy constructor
   IndexVariable(
      const IndexVariable&);

   // Unimplemented assignment operator
   IndexVariable&
   operator = (
      const IndexVariable&);

};

}
}

#include "SAMRAI/pdat/IndexVariable.cpp"

#endif
