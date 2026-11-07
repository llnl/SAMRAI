/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Class for managing transformations between index spaces in
 *                an AMR hierarchy.
 *
 ************************************************************************/
#include "SAMRAI/hier/Transformation.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

namespace SAMRAI {
namespace hier {

/*
 * ************************************************************************
 *
 * Constructors
 *
 * ************************************************************************
 */

Transformation::Transformation(
   RotationIdentifier rotation,
   const IntVector& src_offset,
   const BlockId& begin_block,
   const BlockId& end_block):
   d_rotation(rotation),
   d_offset(src_offset),
   d_begin_block(begin_block),
   d_end_block(end_block)
{
}

Transformation::Transformation(
   const IntVector& src_offset):
   d_rotation(NO_ROTATE),
   d_offset(src_offset),
   d_begin_block(BlockId::invalidId()),
   d_end_block(BlockId::invalidId())
{
}

/*
 * ************************************************************************
 *
 * Transform a box
 *
 * ************************************************************************
 */
void
Transformation::transform(Box& box) const
{
   TBOX_ASSERT(box.getBlockId() == d_begin_block ||
      d_begin_block == BlockId::invalidId());
   box.rotate(d_rotation);
   box.shift(d_offset);
   if (d_begin_block != d_end_block) {
      box.setBlockId(d_end_block);
   }
}

/*
 * ************************************************************************
 *
 * Do inverse tranformation on a box
 *
 * ************************************************************************
 */
void
Transformation::inverseTransform(
   Box& box) const
{
   TBOX_ASSERT(box.getBlockId() == d_end_block ||
      d_end_block == BlockId::invalidId());
   IntVector reverse_offset(d_offset.getDim());
   calculateReverseShift(reverse_offset, d_offset, d_rotation);

   box.rotate(getReverseRotationIdentifier(d_rotation, d_offset.getDim()));
   box.shift(reverse_offset);
   if (d_begin_block != d_end_block) {
      box.setBlockId(d_begin_block);
   }
}

/*
 * ************************************************************************
 *
 * Get Transformation object that is the inverse of 'this'.
 *
 * ************************************************************************
 */
Transformation
Transformation::getInverseTransformation() const
{
   const tbox::Dimension& dim = d_offset.getDim();
   IntVector inv_offset(dim);
   calculateReverseShift(inv_offset, d_offset, d_rotation);

   RotationIdentifier inv_rotate =
      getReverseRotationIdentifier(d_rotation, dim);

   return Transformation(inv_rotate,
      inv_offset,
      d_end_block,
      d_begin_block);
}

/*
 * ************************************************************************
 *
 * Get a RotationIdentifier value associated with given string input
 *
 * ************************************************************************
 */

Transformation::RotationIdentifier
Transformation::getRotationIdentifier(
   const std::vector<std::string>& rotation_string,
   const tbox::Dimension& dim)
{
   TBOX_ASSERT(static_cast<int>(rotation_string.size()) == dim.getValue());

   RotationIdentifier id = NO_ROTATE;
   bool is_error = false;

   if (dim.getValue() == 1) {
      if (rotation_string[0] == "I_UP") {
         id = IUP; //0;
      } else if (rotation_string[0] == "I_DOWN") {
         id = IDOWN; //1;
      } else {
         is_error = true;
      }
      if (is_error) {
         TBOX_ERROR("Rotation_input " << rotation_string[0]
                                      << " is invalid.\n");
      }
   } else if (dim.getValue() == 2) {
      if (rotation_string[0] == "I_UP") {
         if (rotation_string[1] == "J_UP") {
            id = IUP_JUP; //0;
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "I_DOWN") {
         if (rotation_string[1] == "J_DOWN") {
            id = IDOWN_JDOWN; //2;
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "J_UP") {
         if (rotation_string[1] == "I_DOWN") {
            id = JUP_IDOWN; //1;
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "J_DOWN") {
         if (rotation_string[1] == "I_UP") {
            id = JDOWN_IUP; //3;
         } else {
            is_error = true;
         }
      } else {
         is_error = true;
      }
      if (is_error) {
         TBOX_ERROR("Transformation::getRotationIdentifier "
            << rotation_string[0] << " "
            << rotation_string[1] << " "
            << " is invalid.\n");
      }

   } else if (dim.getValue() == 3) {
      if (rotation_string[0] == "I_UP") {
         if (rotation_string[1] == "J_UP") {
            if (rotation_string[2] == "K_UP") {
               id = IUP_JUP_KUP; //0;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "J_DOWN") {
            if (rotation_string[2] == "K_DOWN") {
               id = IUP_JDOWN_KDOWN; //20;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_UP") {
            if (rotation_string[2] == "J_DOWN") {
               id = IUP_KUP_JDOWN; //13;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_DOWN") {
            if (rotation_string[2] == "J_UP") {
               id = IUP_KDOWN_JUP; //7;
            } else {
               is_error = true;
            }
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "I_DOWN") {
         if (rotation_string[1] == "J_UP") {
            if (rotation_string[2] == "K_DOWN") {
               id = IDOWN_JUP_KDOWN; //10;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "J_DOWN") {
            if (rotation_string[2] == "K_UP") {
               id = IDOWN_JDOWN_KUP; //16;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_UP") {
            if (rotation_string[2] == "J_UP") {
               id = IDOWN_KUP_JUP; //3;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_DOWN") {
            if (rotation_string[2] == "J_DOWN") {
               id = IDOWN_KDOWN_JDOWN; //23;
            } else {
               is_error = true;
            }
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "J_UP") {
         if (rotation_string[1] == "I_UP") {
            if (rotation_string[2] == "K_DOWN") {
               id = JUP_IUP_KDOWN; //8;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "I_DOWN") {
            if (rotation_string[2] == "K_UP") {
               id = JUP_IDOWN_KUP; //5;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_UP") {
            if (rotation_string[2] == "I_UP") {
               id = JUP_KUP_IUP; //2;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_DOWN") {
            if (rotation_string[2] == "I_DOWN") {
               id = JUP_KDOWN_IDOWN; //11;
            } else {
               is_error = true;
            }
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "J_DOWN") {
         if (rotation_string[1] == "I_UP") {
            if (rotation_string[2] == "K_UP") {
               id = JDOWN_IUP_KUP; //12;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "I_DOWN") {
            if (rotation_string[2] == "K_DOWN") {
               id = JDOWN_IDOWN_KDOWN; //21;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_UP") {
            if (rotation_string[2] == "I_DOWN") {
               id = JDOWN_KUP_IDOWN; //15;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_DOWN") {
            if (rotation_string[2] == "I_UP") {
               id = JDOWN_KDOWN_IUP; //18;
            } else {
               is_error = true;
            }
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "K_UP") {
         if (rotation_string[1] == "I_UP") {
            if (rotation_string[2] == "J_UP") {
               id = KUP_IUP_JUP; //1;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "I_DOWN") {
            if (rotation_string[2] == "J_DOWN") {
               id = KUP_IDOWN_JDOWN; //17;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "J_UP") {
            if (rotation_string[2] == "I_DOWN") {
               id = KUP_JUP_IDOWN; //4;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "J_DOWN") {
            if (rotation_string[2] == "I_UP") {
               id = KUP_JDOWN_IUP; //14;
            } else {
               is_error = true;
            }
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "K_DOWN") {
         if (rotation_string[1] == "I_UP") {
            if (rotation_string[2] == "J_DOWN") {
               id = KDOWN_IUP_JDOWN; //19;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "I_DOWN") {
            if (rotation_string[2] == "J_UP") {
               id = KDOWN_IDOWN_JUP; //9;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "J_UP") {
            if (rotation_string[2] == "I_UP") {
               id = KDOWN_JUP_IUP; //6;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "J_DOWN") {
            if (rotation_string[2] == "I_DOWN") {
               id = KDOWN_JDOWN_IDOWN; //22;
            } else {
               is_error = true;
            }
         } else {
            is_error = true;
         }
      } else {
         is_error = true;
      }

      if (is_error) {
         TBOX_ERROR("Transformation::getRotationIdentifier "
            << rotation_string[0] << " "
            << rotation_string[1] << " "
            << rotation_string[2]
            << " is invalid.\n");
      }
   } else {
      TBOX_ERROR(
         "Transformation::getRotationIdentifier : DIM > 3 not implemented");
   }

   return id;
}

/*
 * ************************************************************************
 *
 * Get a RotationIdentifier value for the reverse of the given rotation.
 *
 * ************************************************************************
 */

Transformation::RotationIdentifier
Transformation::getReverseRotationIdentifier(
   const RotationIdentifier rotation,
   const tbox::Dimension& dim)
{
   RotationIdentifier reverse_id = (RotationIdentifier)0;

   if (rotation == NO_ROTATE) {
      reverse_id = rotation;
   } else if (dim.getValue() == 1) {
      reverse_id = rotation;
   } else if (dim.getValue() == 2) {
      reverse_id = (RotationIdentifier)((4 - (int)rotation) % 4);
   } else if (dim.getValue() == 3) {
      switch (rotation) {

         case IUP_JUP_KUP:
            reverse_id = IUP_JUP_KUP;
            break;

         case KUP_IUP_JUP:
            reverse_id = JUP_KUP_IUP;
            break;

         case JUP_KUP_IUP:
            reverse_id = KUP_IUP_JUP;
            break;

         case IDOWN_KUP_JUP:
            reverse_id = IDOWN_KUP_JUP;
            break;

         case KUP_JUP_IDOWN:
            reverse_id = KDOWN_JUP_IUP;
            break;

         case JUP_IDOWN_KUP:
            reverse_id = JDOWN_IUP_KUP;
            break;

         case KDOWN_JUP_IUP:
            reverse_id = KUP_JUP_IDOWN;
            break;

         case IUP_KDOWN_JUP:
            reverse_id = IUP_KUP_JDOWN;
            break;

         case JUP_IUP_KDOWN:
            reverse_id = JUP_IUP_KDOWN;
            break;

         case KDOWN_IDOWN_JUP:
            reverse_id = JDOWN_KUP_IDOWN;
            break;

         case IDOWN_JUP_KDOWN:
            reverse_id = IDOWN_JUP_KDOWN;
            break;

         case JUP_KDOWN_IDOWN:
            reverse_id = KDOWN_IUP_JDOWN;
            break;

         case JDOWN_IUP_KUP:
            reverse_id = JUP_IDOWN_KUP;
            break;

         case IUP_KUP_JDOWN:
            reverse_id = IUP_KDOWN_JUP;
            break;

         case KUP_JDOWN_IUP:
            reverse_id = KUP_JDOWN_IUP;
            break;

         case JDOWN_KUP_IDOWN:
            reverse_id = KDOWN_IDOWN_JUP;
            break;

         case IDOWN_JDOWN_KUP:
            reverse_id = IDOWN_JDOWN_KUP;
            break;

         case KUP_IDOWN_JDOWN:
            reverse_id = JDOWN_KDOWN_IUP;
            break;

         case JDOWN_KDOWN_IUP:
            reverse_id = KUP_IDOWN_JDOWN;
            break;

         case KDOWN_IUP_JDOWN:
            reverse_id = JUP_KDOWN_IDOWN;
            break;

         case IUP_JDOWN_KDOWN:
            reverse_id = IUP_JDOWN_KDOWN;
            break;

         case JDOWN_IDOWN_KDOWN:
            reverse_id = JDOWN_IDOWN_KDOWN;
            break;

         case KDOWN_JDOWN_IDOWN:
            reverse_id = KDOWN_JDOWN_IDOWN;
            break;

         case IDOWN_KDOWN_JDOWN:
            reverse_id = IDOWN_KDOWN_JDOWN;
            break;

         default:
            TBOX_ERROR(
            "Transformation::getReverseRotationIdentifier error...\n"
            << " Invalid RotationIdentifier value given" << std::endl);

            reverse_id = IUP_JUP_KUP;
            break;
      }
   } else {
      TBOX_ERROR(
         "Transformation::getReverseRotationIdentifier : DIM > 3 with rotation not implemented");
   }

   return reverse_id;
}

/*
 * ************************************************************************
 *
 * Get the reverse of the given shift
 *
 * ************************************************************************
 */

void
Transformation::calculateReverseShift(
   IntVector& back_shift,
   const IntVector& shift,
   const RotationIdentifier rotation)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(back_shift, shift);

   const tbox::Dimension& dim(back_shift.getDim());

   if (rotation == NO_ROTATE) {
      back_shift = -shift;
   } else if (dim.getValue() == 1) {
      if (rotation == IUP) {
         back_shift = -shift;
      } else if (rotation == IDOWN) {
         back_shift = shift;
      } else {
         TBOX_ERROR("Transformation::calculateReverseShift error...\n"
            << " Invalid RotationIdentifier value given" << std::endl);
      }
   } else if (dim.getValue() == 2) {

      if (rotation == IUP_JUP) {
         back_shift = -shift;
      } else if (rotation == JUP_IDOWN) {
         back_shift(0) = shift(1);
         back_shift(1) = -shift(0);
      } else if (rotation == IDOWN_JDOWN) {
         back_shift(0) = shift(0);
         back_shift(1) = shift(1);
      } else if (rotation == JDOWN_IUP) {
         back_shift(0) = -shift(1);
         back_shift(1) = shift(0);
      } else {
         TBOX_ERROR("Transformation::calculateReverseShift error...\n"
            << " Invalid RotationIdentifier value given" << std::endl);
      }

   } else if (dim.getValue() == 3) {

      RotationIdentifier back_rotation =
         getReverseRotationIdentifier(rotation, dim);

      if (back_rotation == IUP_JUP_KUP) {
         back_shift = -shift;
      } else if (back_rotation == KUP_IUP_JUP) {
         back_shift(0) = -shift(2);
         back_shift(1) = -shift(0);
         back_shift(2) = -shift(1);
      } else if (back_rotation == JUP_KUP_IUP) {
         back_shift(0) = -shift(1);
         back_shift(1) = -shift(2);
         back_shift(2) = -shift(0);
      } else if (back_rotation == IDOWN_KUP_JUP) {
         back_shift(0) = shift(0);
         back_shift(1) = -shift(2);
         back_shift(2) = -shift(1);
      } else if (back_rotation == KUP_JUP_IDOWN) {
         back_shift(0) = -shift(2);
         back_shift(1) = -shift(1);
         back_shift(2) = shift(0);
      } else if (back_rotation == JUP_IDOWN_KUP) {
         back_shift(0) = -shift(1);
         back_shift(1) = shift(0);
         back_shift(2) = -shift(2);
      } else if (back_rotation == KDOWN_JUP_IUP) {
         back_shift(0) = shift(2);
         back_shift(1) = -shift(1);
         back_shift(2) = -shift(0);
      } else if (back_rotation == IUP_KDOWN_JUP) {
         back_shift(0) = -shift(0);
         back_shift(1) = shift(2);
         back_shift(2) = -shift(1);
      } else if (back_rotation == JUP_IUP_KDOWN) {
         back_shift(0) = -shift(1);
         back_shift(1) = -shift(0);
         back_shift(2) = shift(2);
      } else if (back_rotation == KDOWN_IDOWN_JUP) {
         back_shift(0) = shift(2);
         back_shift(1) = shift(0);
         back_shift(2) = -shift(1);
      } else if (back_rotation == IDOWN_JUP_KDOWN) {
         back_shift(0) = shift(0);
         back_shift(1) = -shift(1);
         back_shift(2) = shift(2);
      } else if (back_rotation == JUP_KDOWN_IDOWN) {
         back_shift(0) = -shift(1);
         back_shift(1) = shift(2);
         back_shift(2) = shift(0);
      } else if (back_rotation == JDOWN_IUP_KUP) {
         back_shift(0) = shift(1);
         back_shift(1) = -shift(0);
         back_shift(2) = -shift(2);
      } else if (back_rotation == IUP_KUP_JDOWN) {
         back_shift(0) = -shift(0);
         back_shift(1) = -shift(2);
         back_shift(2) = shift(1);
      } else if (back_rotation == KUP_JDOWN_IUP) {
         back_shift(0) = -shift(2);
         back_shift(1) = shift(1);
         back_shift(2) = -shift(0);
      } else if (back_rotation == JDOWN_KUP_IDOWN) {
         back_shift(0) = shift(1);
         back_shift(1) = -shift(2);
         back_shift(2) = shift(0);
      } else if (back_rotation == IDOWN_JDOWN_KUP) {
         back_shift(0) = shift(0);
         back_shift(1) = shift(1);
         back_shift(2) = -shift(2);
      } else if (back_rotation == KUP_IDOWN_JDOWN) {
         back_shift(0) = -shift(2);
         back_shift(1) = shift(0);
         back_shift(2) = shift(1);
      } else if (back_rotation == JDOWN_KDOWN_IUP) {
         back_shift(0) = shift(1);
         back_shift(1) = shift(2);
         back_shift(2) = -shift(0);
      } else if (back_rotation == KDOWN_IUP_JDOWN) {
         back_shift(0) = shift(2);
         back_shift(1) = -shift(0);
         back_shift(2) = shift(1);
      } else if (back_rotation == IUP_JDOWN_KDOWN) {
         back_shift(0) = -shift(0);
         back_shift(1) = shift(1);
         back_shift(2) = shift(2);
      } else if (back_rotation == JDOWN_IDOWN_KDOWN) {
         back_shift(0) = shift(1);
         back_shift(1) = shift(0);
         back_shift(2) = shift(2);
      } else if (back_rotation == KDOWN_JDOWN_IDOWN) {
         back_shift(0) = shift(2);
         back_shift(1) = shift(1);
         back_shift(2) = shift(0);
      } else if (back_rotation == IDOWN_KDOWN_JDOWN) {
         back_shift(0) = shift(0);
         back_shift(1) = shift(2);
         back_shift(2) = shift(1);
      } else {
         TBOX_ERROR("Transformation::calculateReverseShift error...\n"
            << " Invalid RotationIdentifier value given" << std::endl);
      }
   } else {
      TBOX_ERROR("Transformation::calculateReverseShift : DIM > 3 with rotation not implemented");
   }
}

/*
 *************************************************************************
 *
 * rotate an index around the origin.
 *
 *************************************************************************
 */

void
Transformation::rotateIndex(
   int* index,
   const tbox::Dimension& dim,
   const RotationIdentifier rotation)
{
   if (dim.getValue() == 1) {
      if (rotation == IUP) {
         return;
      } else if (rotation == IDOWN) {
         index[0] = -index[0] - 1;
      }
   } else if (dim.getValue() == 2) {
      int num_rotations = static_cast<int>(rotation);

      for (int j = 0; j < num_rotations; ++j) {
         int tmp_in[2];
         tmp_in[0] = index[0];
         tmp_in[1] = index[1];

         index[0] = tmp_in[1];
         index[1] = -tmp_in[0] - 1;
      }
   } else if (dim.getValue() == 3) {
      if (rotation == IUP_JUP_KUP) {
         return;
      } else if (rotation == KUP_IUP_JUP) {
         rotateAboutAxis(dim, index, 0, 3);
         rotateAboutAxis(dim, index, 2, 3);
      } else if (rotation == JUP_KUP_IUP) {
         rotateAboutAxis(dim, index, 1, 1);
         rotateAboutAxis(dim, index, 2, 1);
      } else if (rotation == IDOWN_KUP_JUP) {
         rotateAboutAxis(dim, index, 1, 2);
         rotateAboutAxis(dim, index, 0, 3);
      } else if (rotation == KUP_JUP_IDOWN) {
         rotateAboutAxis(dim, index, 1, 3);
      } else if (rotation == JUP_IDOWN_KUP) {
         rotateAboutAxis(dim, index, 2, 1);
      } else if (rotation == KDOWN_JUP_IUP) {
         rotateAboutAxis(dim, index, 1, 1);
      } else if (rotation == IUP_KDOWN_JUP) {
         rotateAboutAxis(dim, index, 0, 3);
      } else if (rotation == JUP_IUP_KDOWN) {
         rotateAboutAxis(dim, index, 0, 2);
         rotateAboutAxis(dim, index, 2, 3);
      } else if (rotation == KDOWN_IDOWN_JUP) {
         rotateAboutAxis(dim, index, 0, 3);
         rotateAboutAxis(dim, index, 2, 1);
      } else if (rotation == IDOWN_JUP_KDOWN) {
         rotateAboutAxis(dim, index, 1, 2);
      } else if (rotation == JUP_KDOWN_IDOWN) {
         rotateAboutAxis(dim, index, 0, 3);
         rotateAboutAxis(dim, index, 1, 3);
      } else if (rotation == JDOWN_IUP_KUP) {
         rotateAboutAxis(dim, index, 2, 3);
      } else if (rotation == IUP_KUP_JDOWN) {
         rotateAboutAxis(dim, index, 0, 1);
      } else if (rotation == KUP_JDOWN_IUP) {
         rotateAboutAxis(dim, index, 0, 2);
         rotateAboutAxis(dim, index, 1, 1);
      } else if (rotation == JDOWN_KUP_IDOWN) {
         rotateAboutAxis(dim, index, 0, 1);
         rotateAboutAxis(dim, index, 1, 3);
      } else if (rotation == IDOWN_JDOWN_KUP) {
         rotateAboutAxis(dim, index, 0, 2);
         rotateAboutAxis(dim, index, 1, 2);
      } else if (rotation == KUP_IDOWN_JDOWN) {
         rotateAboutAxis(dim, index, 0, 1);
         rotateAboutAxis(dim, index, 2, 1);
      } else if (rotation == JDOWN_KDOWN_IUP) {
         rotateAboutAxis(dim, index, 0, 3);
         rotateAboutAxis(dim, index, 1, 1);
      } else if (rotation == KDOWN_IUP_JDOWN) {
         rotateAboutAxis(dim, index, 0, 1);
         rotateAboutAxis(dim, index, 2, 3);
      } else if (rotation == IUP_JDOWN_KDOWN) {
         rotateAboutAxis(dim, index, 0, 2);
      } else if (rotation == JDOWN_IDOWN_KDOWN) {
         rotateAboutAxis(dim, index, 0, 2);
         rotateAboutAxis(dim, index, 2, 1);
      } else if (rotation == KDOWN_JDOWN_IDOWN) {
         rotateAboutAxis(dim, index, 0, 2);
         rotateAboutAxis(dim, index, 1, 3);
      } else if (rotation == IDOWN_KDOWN_JDOWN) {
         rotateAboutAxis(dim, index, 1, 2);
         rotateAboutAxis(dim, index, 0, 1);
      }
   } else {
      TBOX_ERROR("Transformation::rotateIndex : DIM > 3 not implemented");
   }

}

/*
 *************************************************************************
 *
 * Private routine to rotate an index about an axis.
 *
 *************************************************************************
 */

void
Transformation::rotateAboutAxis(
   const tbox::Dimension& dim,
   int* index,
   const int axis,
   const int num_rotations)
{
   if (dim.getValue() == 3) {
      TBOX_ASSERT(axis < dim.getValue());

      const int a = (axis + 1) % dim.getValue();
      const int b = (axis + 2) % dim.getValue();

      for (int j = 0; j < num_rotations; ++j) {
         int tmp_in[3] = { index[0], index[1], index[2] };
         index[a] = tmp_in[b];
         index[b] = -tmp_in[a] - 1;
      }
   }
}

void
Transformation::setOrientationVector(
   std::vector<int>& orientation,
   RotationIdentifier rotation)
{
    orientation.resize(3);
    if (rotation == NO_ROTATE) {
        orientation[0] = 1;
        orientation[1] = 2;
        orientation[2] = 3;
    } else if (rotation == IUP_JUP) {
        orientation[0] = 1;
        orientation[1] = 2;
        orientation[2] = 3;
    } else if (rotation == JUP_IDOWN) {
        orientation[0] = 2;
        orientation[1] = -1;
        orientation[2] = 3;
    } else if (rotation == IDOWN_JDOWN) {
        orientation[0] = -1;
        orientation[1] = -2;
        orientation[2] = 3;
    } else if (rotation == JDOWN_IUP) {
        orientation[0] = -2;
        orientation[1] = 1;
        orientation[2] = 3;
    } else if (rotation == IUP_JUP_KUP) {
        orientation[0] = 1;
        orientation[1] = 2;
        orientation[2] = 3;
    } else if (rotation == KUP_IUP_JUP) {
        orientation[0] = 3;
        orientation[1] = 1;
        orientation[2] = 2;
    } else if (rotation == JUP_KUP_IUP) {
        orientation[0] = 2;
        orientation[1] = 3;
        orientation[2] = 1;
    } else if (rotation == IDOWN_KUP_JUP) {
        orientation[0] = -1;
        orientation[1] = 3;
        orientation[2] = 2;
    } else if (rotation == KUP_JUP_IDOWN) {
        orientation[0] = 3;
        orientation[1] = 2;
        orientation[2] = -1;
    } else if (rotation == JUP_IDOWN_KUP) {
        orientation[0] = 2;
        orientation[1] = -1;
        orientation[2] = 3;
    } else if (rotation == KDOWN_JUP_IUP) {
        orientation[0] = -3;
        orientation[1] = 2;
        orientation[2] = 1;
    } else if (rotation == IUP_KDOWN_JUP) {
        orientation[0] = 1;
        orientation[1] = -3;
        orientation[2] = 2;
    } else if (rotation == JUP_IUP_KDOWN) {
        orientation[0] = 2;
        orientation[1] = 1;
        orientation[2] = -3;
    } else if (rotation == KDOWN_IDOWN_JUP) {
        orientation[0] = -3;
        orientation[1] = -1;
        orientation[2] = 2;
    } else if (rotation == IDOWN_JUP_KDOWN) {
        orientation[0] = -1;
        orientation[1] = 2;
        orientation[2] = -3;
    } else if (rotation == JUP_KDOWN_IDOWN) {
        orientation[0] = 2;
        orientation[1] = -3;
        orientation[2] = -1;
    } else if (rotation == JDOWN_IUP_KUP) {
        orientation[0] = -2;
        orientation[1] = 1;
        orientation[2] = 3;
    } else if (rotation == IUP_KUP_JDOWN) {
        orientation[0] = 1;
        orientation[1] = 3;
        orientation[2] = -2;
    } else if (rotation == KUP_JDOWN_IUP) {
        orientation[0] = 3;
        orientation[1] = -2;
        orientation[2] = 1;
    } else if (rotation == JDOWN_KUP_IDOWN) {
        orientation[0] = -2;
        orientation[1] = 3;
        orientation[2] = -1;
    } else if (rotation == IDOWN_JDOWN_KUP) {
        orientation[0] = -1;
        orientation[1] = -2;
        orientation[2] = 3;
    } else if (rotation == KUP_IDOWN_JDOWN) {
        orientation[0] = 3;
        orientation[1] = -1;
        orientation[2] = -2;
    } else if (rotation == JDOWN_KDOWN_IUP) {
        orientation[0] = -2;
        orientation[1] = -3;
        orientation[2] = 1;
    } else if (rotation == KDOWN_IUP_JDOWN) {
        orientation[0] = -3;
        orientation[1] = 1;
        orientation[2] = -2;
    } else if (rotation == IUP_JDOWN_KDOWN) {
        orientation[0] = 1;
        orientation[1] = -2;
        orientation[2] = -3;
    } else if (rotation == JDOWN_IDOWN_KDOWN) {
        orientation[0] = -2;
        orientation[1] = -1;
        orientation[2] = -3;
    } else if (rotation == KDOWN_JDOWN_IDOWN) {
        orientation[0] = -3;
        orientation[1] = -2;
        orientation[2] = -1;
    } else if (rotation == IDOWN_KDOWN_JDOWN) {
        orientation[0] = -1;
        orientation[1] = -3;
        orientation[2] = -2;
    }
}


}
}
