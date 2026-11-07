/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_FaceVariable
#define included_pdat_FaceVariable

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#include <string>

namespace SAMRAI {
namespace pdat {

/*!
 * Class FaceVariable<TYPE> is a templated variable class used to define
 * face-centered quantities on an AMR mesh.   It is a subclass of
 * hier::Variable and is templated on the type of the underlying data
 * (e.g., double, int, bool, etc.).
 *
 * Note that the indices in the face data arrays are permuted so that
 * the leading index in each array corresponds to the associated face
 * normal coordinate direction. See header file for FaceData<TYPE> class
 * for a more detailed description of the data layout.
 *
 * IMPORTANT: The class SideVariable<TYPE> and associated "side data" classes
 * define the same storage as this face variable class, except that the
 * individual array indices are not permuted in the side data type.
 *
 * @see FaceData
 * @see FaceDataFactory
 * @see FaceGeometry
 * @see hier::Variable
 */

template<class TYPE>
class FaceVariable:public hier::Variable
{
public:
   /*!
    * @brief Create an face-centered variable object with the given name and
    * depth (i.e., number of data values at each edge index location).
    * A default depth of one is provided.   The fine boundary representation
    * boolean argument indicates which values (either coarse or fine) take
    * precedence at coarse-fine mesh boundaries during coarsen and refine
    * operations.  The default is that fine data values take precedence
    * on coarse-fine interfaces.
    */
   FaceVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      int depth = 1,
      bool fine_boundary_represents_var = true);

   /*!
    * @brief Constructor that also includes an Umpire allocator for
    * allocations of the underlying data.
    */
   FaceVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      tbox::ResourceAllocator allocator,
      int depth = 1,
      bool fine_boundary_represents_var = true);

   /*!
    * @brief Virtual destructor for face variable objects.
    */
   virtual ~FaceVariable();

   /*!
    * @brief Return boolean indicating which face data values (coarse
    * or fine) take precedence at coarse-fine mesh interfaces.  The
    * value is set in the constructor.
    */
   bool fineBoundaryRepresentsVariable() const
   {
      return d_fine_boundary_represents_var;
   }

   /*!
    * @brief Return true indicating that face data on a patch interior
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
   FaceVariable(
      const FaceVariable&);

   // Unimplemented assignment operator
   FaceVariable&
   operator = (
      const FaceVariable&);
};

}
}

#include "SAMRAI/pdat/FaceVariable.cpp"

#endif
