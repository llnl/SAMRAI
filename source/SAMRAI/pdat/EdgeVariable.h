/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_EdgeVariable
#define included_pdat_EdgeVariable

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/Complex.h"

#include <string>

namespace SAMRAI {
namespace pdat {

/*!
 * Class EdgeVariable<TYPE> is a templated variable class used to define
 * edge-centered quantities on an AMR mesh.   It is a subclass of
 * hier::Variable and is templated on the type of the underlying data
 * (e.g., double, int, bool, etc.).
 *
 * See header file for EdgeData<TYPE> class for a more detailed
 * description of the data layout.
 *
 * @see EdgeData
 * @see EdgeDataFactory
 * @see hier::Variable
 */

template<class TYPE>
class EdgeVariable:public hier::Variable
{
public:
   /*!
    * @brief Create an edge-centered variable object with the given name and
    * depth (i.e., number of data values at each edge index location).
    * A default depth of one is provided.   The fine boundary representation
    * boolean argument indicates which values (either coarse or fine) take
    * precedence at coarse-fine mesh boundaries during coarsen and refine
    * operations.  The default is that fine data values take precedence
    * on coarse-fine interfaces.
    */
   EdgeVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      int depth = 1,
      bool fine_boundary_represents_var = true);

   /*!
    * @brief Virtual destructor for edge variable objects.
    */
   virtual ~EdgeVariable();

   /*!
    * @brief Return boolean indicating which edge data values (coarse
    * or fine) take precedence at coarse-fine mesh interfaces.  The
    * value is set in the constructor.
    */
   bool fineBoundaryRepresentsVariable() const
   {
      return d_fine_boundary_represents_var;
   }

   /*!
    * @brief Return true indicating that edge data on a patch interior
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
   bool d_fine_boundary_represents_var;

   // Unimplemented copy constructor
   EdgeVariable(
      const EdgeVariable&);

   // Unimplemented assignment operator
   EdgeVariable&
   operator = (
      const EdgeVariable&);

};

}
}

#include "SAMRAI/pdat/EdgeVariable.cpp"

#endif
