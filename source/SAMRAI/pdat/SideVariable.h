/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_SideVariable
#define included_pdat_SideVariable

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#include <string>

namespace SAMRAI {
namespace pdat {

/*!
 * Class SideVariable<TYPE> is a templated variable class used to define
 * side-centered quantities on an AMR mesh.   It is a subclass of
 * hier::Variable and is templated on the type of the underlying data
 * (e.g., double, int, bool, etc.).  See header file for SideData<TYPE> class
 * for a more detailed description of the data layout.
 *
 * Note that it is possible to create a side data object for managing
 * data at cell sides associated with a single coordinate direction only.
 * See the constructor for more information.
 *
 * IMPORTANT: The class FaceVariable<TPYE> and associated "face data" classes
 * define the same storage as this side variable class, except that the
 * individual array indices are permuted in the face data type.
 *
 * @see SideData
 * @see SideDataFactory
 * @see SideGeometry
 * @see hier::Variable
 */

template<class TYPE>
class SideVariable:public hier::Variable
{
public:
   /*!
    * @brief Create an side-centered variable object with the given name and
    * depth (i.e., number of data values at each edge index location).
    *
    * A default depth of one is provided.   The fine boundary representation
    * boolean argument indicates which values (either coarse or fine) take
    * precedence at coarse-fine mesh boundaries during coarsen and refine
    * operations.  The default is that fine data values take precedence
    * on coarse-fine interfaces.
    */
   SideVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      const hier::IntVector& directions,
      int depth = 1,
      bool fine_boundary_represents_var = true);

   /*!
    * @brief Constructor that will assume a directions vector of 1 in
    * all coordinate directions.
    *
    * Constructs a SideVariable the same as the first constructor, except
    * there is no directions argument, meaning that data is expected to
    * exist on the sides in all coordinate directions.
    */
   SideVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      int depth = 1,
      bool fine_boundary_represents_var = true);

   /*!
    * @brief Constructor with directions vector and  an Umpire allocator for
    * allocations of the underlying data.
    */
   SideVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      const hier::IntVector& directions,
      tbox::ResourceAllocator allocator,
      int depth = 1,
      bool fine_boundary_represents_var = true);

   /*!
    * @brief Constructor that asumes all coordinate directions are used and
    * also includes an Umpire allocator for allocations of the underlying
    * data.
    */
   SideVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      tbox::ResourceAllocator allocator,
      int depth = 1,
      bool fine_boundary_represents_var = true);

   /*!
    * @brief Virtual destructor for side variable objects.
    */
   virtual ~SideVariable();

   /*!
    * @brief Return constant reference to vector describing which coordinate
    * directions have data associated with this side variable object.
    *
    * A vector entry of zero indicates that there is no data array
    * allocated for the corresponding coordinate direction for side data
    * created via this side variable object.  A non-zero value indicates
    * that a valid data array will be allocated for that coordinate
    * direction.
    */
   const hier::IntVector&
   getDirectionVector() const;

   /*!
    * @brief Return boolean indicating which side data values (coarse
    * or fine) take precedence at coarse-fine mesh interfaces.  The
    * value is set in the constructor.
    */
   bool fineBoundaryRepresentsVariable() const
   {
      return d_fine_boundary_represents_var;
   }

   /*!
    * @brief Return true indicating that side data on a patch interior
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
   hier::IntVector d_directions;

   // Unimplemented copy constructor
   SideVariable(
      const SideVariable&);

   // Unimplemented assignment operator
   SideVariable&
   operator = (
      const SideVariable&);

};

}
}

#include "SAMRAI/pdat/SideVariable.cpp"

#endif
