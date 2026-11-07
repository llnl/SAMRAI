/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Class for managing tanssformations between index spaces in
 *                an AMR hierarchy.
 *
 ************************************************************************/

#ifndef included_hier_Transformation
#define included_hier_Transformation

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

namespace SAMRAI {
namespace hier {

class Box;
class Patch;

/*!
 * @brief Class Transformation represents generalized transformations between
 * coordinate systems or within coordinate systems.
 *
 * Transformation is used to provide a way to handle communication of data
 * between patches that lie on different coordinate systems (such as different
 * blocks in a multiblock mesh) or communication between opposite sides of
 * a periodic mesh.  This class is able to transform the index space of a
 * source patch to the index space of a destination patch (or vice versa) so
 * that the overlap between data on the patches can be computed.
 *
 * Transformation can also handle the trival case where the source and
 * destination index spaces need no transformation.
 *
 * The Transformation is represented by a RotationIdentifier and an offset
 * IntVector.  When transforming a Box, this class rotates the Box around the
 * origin according to RotationIdentifier, and then shifts the Box by the
 * offset.  When dealing with periodic conditions, there rotation is always
 * the zero value NO_ROTATE and only the offset shift is invoked.  For
 * trivial transformations, the rotation is NO_ROTATE and the IntVector is
 * a zero vector.
 *
 * @see BoxOverlap
 * @see BoxGeometry
 * @see VariableFillPattern
 */

class Transformation
{
public:
   /*!
    * @brief Type that identifies the rotation relationship between two blocks.
    */
   enum RotationIdentifier {
      NO_ROTATE = 0,
      IUP = 0,
      IDOWN = 1,
      IUP_JUP = 0,
      JUP_IDOWN = 1,
      IDOWN_JDOWN = 2,
      JDOWN_IUP = 3,
      IUP_JUP_KUP = 0,
      KUP_IUP_JUP = 1,
      JUP_KUP_IUP = 2,
      IDOWN_KUP_JUP = 3,
      KUP_JUP_IDOWN = 4,
      JUP_IDOWN_KUP = 5,
      KDOWN_JUP_IUP = 6,
      IUP_KDOWN_JUP = 7,
      JUP_IUP_KDOWN = 8,
      KDOWN_IDOWN_JUP = 9,
      IDOWN_JUP_KDOWN = 10,
      JUP_KDOWN_IDOWN = 11,
      JDOWN_IUP_KUP = 12,
      IUP_KUP_JDOWN = 13,
      KUP_JDOWN_IUP = 14,
      JDOWN_KUP_IDOWN = 15,
      IDOWN_JDOWN_KUP = 16,
      KUP_IDOWN_JDOWN = 17,
      JDOWN_KDOWN_IUP = 18,
      KDOWN_IUP_JDOWN = 19,
      IUP_JDOWN_KDOWN = 20,
      JDOWN_IDOWN_KDOWN = 21,
      KDOWN_JDOWN_IDOWN = 22,
      IDOWN_KDOWN_JDOWN = 23
   };

   /*!
    * @brief Constructor to set rotation and offset
    *
    * If begin_block and end_block have different values, all Boxes
    * transformed by an object constructed with this constructor will have
    * their BlockIds changed from begin_block to end_block.
    *
    * @param[in] rotation  specifies rotation, if any, between blocks
    * @param[in] offset    offset to be applied after rotation
    * @param[in] begin_block  Block before the transformation
    * @param[in] end_block    Block after the transformation
    */
   Transformation(
      const RotationIdentifier rotation,
      const IntVector& offset,
      const BlockId& begin_block,
      const BlockId& end_block);

   /*!
    * @brief Constructor that sets only offset
    *
    * The RotationIdentifier is set to NO_ROTATE.  This constructor should
    * be used in cases where the calling code knows that the Transformation
    * is trivial or periodic.  Pass in a zero IntVector to construct a trivial
    * Transformation.
    *
    * Any Boxes transformed by an object constructed with this constructor
    * will not have their BlockIds changed.
    *
    * @param[in]  offset
    */
   explicit Transformation(
      const IntVector& offset);

   /*!
    * @brief Copy constructor
    *
    * @param[in] copy_trans
    */
   Transformation(
      const Transformation& copy_trans) = default;

   Transformation(Transformation&& other) = default;

   /*!
    * @brief Destructor
    */
   ~Transformation() = default;

   /*!
    * @brief Get the rotation
    */
   constexpr RotationIdentifier
   getRotation() const noexcept
   {
      return d_rotation;
   }

   /*!
    * @brief Get the offset
    */
   const IntVector&
   getOffset() const noexcept
   {
      return d_offset;
   }

   /*!
    * @brief Transform the Box in the way defined by this object
    *
    * @param[in,out] box  The Box will be transformed
    *
    * @pre (box.getBlockId() == getBeginBlock()) ||
    *      (getBeginBlock() == BlockId::invalidId())
    */
   void
   transform(
      Box& box) const;

   /*!
    * @brief Apply the inverse of this object's transformation to the Box
    *
    * If transform() and inverseTransform() are called consecutively on
    * a Box, the Box will end up in its original state.
    *
    * @param[in,out] box
    *
    * @pre (box.getBlockId() == getBeginBlock()) ||
    *      (getBeginBlock() == BlockId::invalidId())
    */
   void
   inverseTransform(
      Box& box) const;

   /*!
    * @brief Get a Tranformation object that defines the inverse of this
    * tranformation.
    */
   Transformation
   getInverseTransformation() const;

   /*!
    * @brief Assignment operator
    */
   Transformation&
   operator = (
      const Transformation& rhs) = default;

   Transformation&
   operator = (
      Transformation&& rhs) = default;

   /*!
    * @brief Get the BlockId for the Box before transformation.
    */
   constexpr const BlockId&
   getBeginBlock() const noexcept
   {
      return d_begin_block;
   }

   /*!
    * @brief Get the BlockId for the Box after transformation.
    */
   constexpr const BlockId&
   getEndBlock() const noexcept
   {
      return d_end_block;
   }

   /*!
    * @brief Map a string-based identifier of a rotation operation to a
    * RotationIdentifier value.
    *
    * The rotation_string array must have each entry intended to be used
    * to create a value of the enumerated type RotationIdentifier.  For
    * example, to create the value JUP_IDOWN_KUP, the strings in the array
    * should be "J_UP", "I_DOWN", and "K_UP", in that order.
    * See the comments for the RotationIdentifier definition for the
    * explanation of the meaning of the RotationIdentifier values.
    *
    * A run-time error will occur if the strings do not match the format
    * needed to create a valid Transformation::RotationIdentifier value.
    *
    * @return RotationIdentifier determined by the strings.
    *
    * @param[in] rotation_string
    * @param[in] dim
    *
    * @pre rotation_string.size() == dim.getValue()
    * @pre (dim.getValue() == 1) || (dim.getValue() == 2) ||
    *      (dim.getValue() == 3)
    */
   static RotationIdentifier
   getRotationIdentifier(
      const std::vector<std::string>& rotation_string,
      const tbox::Dimension& dim);

   /*!
    * @brief static method to get a reverse rotation identifier
    *
    * A rotation identifier signifies a specific rotation of an index space.
    * For each rotation there is another rotation that rotates in the exact
    * opposite manner.  This routine returns the identifier of the reverse
    * rotation corresponding to the given rotation.
    *
    * @return  The RotationIdentifier opposite the input parameter
    *
    * @param[in] rotation Rotation for which the reverse rotation is sought
    * @param[in] dim      Dimension being used
    *
    * @pre (dim.getValue() == 1) || (dim.getValue() == 2) ||
    *      (dim.getValue() == 3)
    */
   static RotationIdentifier
   getReverseRotationIdentifier(
      const RotationIdentifier rotation,
      const tbox::Dimension& dim);

   /*!
    * @brief static method to get a reverse shift.
    *
    * Given a rotation and shift that define the relationship between two
    * neighboring blocks, get a reverse shift that, combined with the
    * reverse rotation from getReverseRotationIdentifier, can be used
    * to reverse the effect of the original rotation and shift.
    *
    * @param[out] back_shift
    * @param[in] shift
    * @param[in] rotation
    *
    * @pre back_shift.getDim() == shift.getDim()
    * @pre (dim.getValue() == 1) || (dim.getValue() == 2) ||
    *      (dim.getValue() == 3)
    */
   static void
   calculateReverseShift(
      IntVector& back_shift,
      const IntVector& shift,
      const RotationIdentifier rotation);

   /*!
    * @brief rotate a cell centered index from one index space to another
    *
    * The parameter index is an int pointer with points to an array of
    * int data, length DIM.  It signifies an ijk location in a cell centered
    * index space.  According to the rotation number, the location will be
    * rotated around the origin, with the new values overwriting the original
    * values in the array pointed to by index.
    *
    * @param index array identifying a cell centered point in index space
    * @param dim Dimension of the index and the hierarchy where it is located
    * @param rotation    identifier of the rotation that will be applied
    *                        to index
    *
    * @pre (dim.getValue() == 1) || (dim.getValue() == 2) ||
    *      (dim.getValue() == 3)
    */
   static void
   rotateIndex(
      int* index,
      const tbox::Dimension& dim,
      const RotationIdentifier rotation);

   /*!
    * @brief rotate a cell centered index from one index space to another
    *
    * According to the rotation number, the location of the given cell centered
    * Index will be rotated around the origin, overwriting the original value
    * of the Index.
    *
    * @param index a cell centered point in index space
    * @param rotation    identifier of the rotation that will be applied
    *                        to index
    */
   static void
   rotateIndex(
      Index& index,
      const RotationIdentifier rotation)
   {
      rotateIndex(&index[0], index.getDim(), rotation);
   }

   static void
   setOrientationVector(
      std::vector<int>& orientation,
      RotationIdentifier rotation);

private:
   /*!
    * @brief private routine to rotate a cell centered index around an axis
    *
    * In 3D, rotation of a cell centered index about the origin is decomposed
    * into a series of rotations about an axis.  This function performs one
    * such rotation.
    *
    * @param dim space associated rotation performed in
    * @param index array identifying a cell centered point in index space
    * @param axis axis around which index will be rotated
    * @param num_rotations number of 90-degree rotations around the axis
    *
    * @pre (dim.getValue() != 3) || (axis < dim.getValue())
    */
   static void
   rotateAboutAxis(
      const tbox::Dimension& dim,
      int* index,
      const int axis,
      const int num_rotations);

   /*
    * Unimplemented default constructor.
    */
   Transformation() = delete;

   RotationIdentifier d_rotation;
   IntVector d_offset;

   BlockId d_begin_block;
   BlockId d_end_block;

};

}
}

#endif
