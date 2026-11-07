/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_NodeVariable
#define included_pdat_NodeVariable

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#include <string>

namespace SAMRAI {
namespace pdat {

/*!
 * Class NodeVariable<TYPE> is a templated variable class used to define
 * node-centered quantities on an AMR mesh.   It is a subclass of
 * hier::Variable and is templated on the type of the underlying data
 * (e.g., double, int, bool, etc.).
 *
 * See header file for NodeData<TYPE> class for a more detailed
 * description of the data layout.
 *
 * @see NodeData
 * @see NodeDataFactory
 * @see hier::Variable
 */

template<class TYPE>
class NodeVariable:public hier::Variable
{
public:
   /*!
    * @brief Create a node-centered variable object with the given name and
    * depth (i.e., number of data values at each node index location).
    * A default depth of one is provided.  The fine boundary representation
    * boolean argument indicates which values (either coarse or fine) take
    * precedence at coarse-fine mesh boundaries during coarsen and refine
    * operations.  The default is that fine data values take precedence
    * on coarse-fine interfaces.
    */
   NodeVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      int depth = 1,
      bool fine_boundary_represents_var = true);

   /*!
    * @brief Constructor that also includes an Umpire allocator for
    * allocations of the underlying data.
    */
   NodeVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      tbox::ResourceAllocator allocator,
      int depth = 1,
      bool fine_boundary_represents_var = true);

   /*!
    * @brief Virtual destructor for node variable objects.
    */
   virtual ~NodeVariable();

   /*!
    * @brief Return boolean indicating which node data values (coarse
    * or fine) take precedence at coarse-fine mesh interfaces.  The
    * value is set in the constructor.
    */
   bool fineBoundaryRepresentsVariable() const
   {
      return d_fine_boundary_represents_var;
   }

   /*!
    * @brief Return true indicating that node data on a patch interior
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
   NodeVariable(
      const NodeVariable&);

   // Unimplemented assignment operator
   NodeVariable&
   operator = (
      const NodeVariable&);
};

}
}

#include "SAMRAI/pdat/NodeVariable.cpp"

#endif
