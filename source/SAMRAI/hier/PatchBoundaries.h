/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Describes boundaries for a patch
 *
 ************************************************************************/

#ifndef included_hier_PatchBoundaries
#define included_hier_PatchBoundaries

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/tbox/Utilities.h"

#include <map>
#include <vector>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

/*!
 * @brief Class PatchBoundaries is a container class for storing
 * BoundaryBox objects for a single patch.
 *
 * @see BoundaryBox
 */

class PatchBoundaries
{
public:
   /*!
    * @brief Explicit constructor requires dimension argument.
    *
    * @param[in] dim
    */
   explicit PatchBoundaries(
      const tbox::Dimension& dim);

   /*!
    * @brief Copy constructor.
    *
    * @param[in] r  Patchboundaries object to be copied in constructor.
    */
   PatchBoundaries(
      const PatchBoundaries& r) = default;

   PatchBoundaries(
      PatchBoundaries&& r) = default;

   /*!
    * @brief Assignment operator.
    *
    * @param[in] r  Patchboundaries object to be copied in assignment.
    */
   PatchBoundaries&
   operator = (
      const PatchBoundaries& r)
   {
      for (unsigned int d = 0; d < getDim().getValue(); ++d) {
         d_array_of_bboxes[d] = r.d_array_of_bboxes[d];
      }
      return *this;
   }

   /*!
    * @brief Vector access operator.
    *
    * @param[in] i  Vector index.
    *
    * @pre i < getDim().getValue()
    */
   std::vector<BoundaryBox>&
   operator [] (
      unsigned int i)
   {
      TBOX_ASSERT(i < getDim().getValue());
      return d_array_of_bboxes[i];
   }

   /*!
    * @brief Const Vector access operator.
    *
    * @param[in] i  Vector index.
    *
    * @pre i < getDim().getValue()
    */
   const std::vector<BoundaryBox>&
   operator [] (
      unsigned int i) const
   {
      TBOX_ASSERT(i < getDim().getValue());
      return d_array_of_bboxes[i];
   }

   /*!
    * @brief Get copy of the internal vectors.
    *
    * @return  Copy of the internal vectors.
    */
   std::vector<std::vector<BoundaryBox> >&
   getVectors()
   {
      return d_array_of_bboxes;
   }

   /*!
    * @brief Get const copy of the internal vectors.
    *
    * @return  Const copy of the internal vectors.
    */
   const std::vector<std::vector<BoundaryBox> >&
   getVectors() const
   {
      return d_array_of_bboxes;
   }

   const tbox::Dimension&
   getDim() const
   {
      return d_dim;
   }

private:
   /*
    * Unimplemented default constructor.
    */
   PatchBoundaries() = delete;

   /*!
    * @brief Dimension of the object.
    */
   const tbox::Dimension d_dim;

   /*
    * @brief Internal arrays of BoundaryBox
    */
   std::vector<std::vector<BoundaryBox> > d_array_of_bboxes;
};

} // SAMRAI namespace
} // hier namespace

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
