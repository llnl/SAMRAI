/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated class providing RAJA indexing into data
 *
 ************************************************************************/

#ifndef included_pdat_ArrayView
#define included_pdat_ArrayView

#include "SAMRAI/SAMRAI_config.h"

#if defined(HAVE_RAJA)

#include "RAJA/RAJA.hpp"

namespace SAMRAI {
namespace pdat {

namespace detail {

struct layout_traits {
   using Layout1d = RAJA::OffsetLayout<1, RAJA::Index_type>;
   using Layout2d = RAJA::OffsetLayout<2, RAJA::Index_type>;
   using Layout3d = RAJA::OffsetLayout<3, RAJA::Index_type>;
};

} // namespace detail

/*!
 * @brief ArrayView<DIM,TYPE> is a templated struct that provides an
 * indexing interface into the arrays held by class ArrayData<TYPE> for
 * use within RAJA loops as defined by the hier::parallel_for_all loops
 * provided in the file ForAll.h
 *
 * This can be used with any of the standard PatchData implementations
 * that depend on ArrayData in the pdat component of SAMRAI: CellData,
 * FaceData, etc.
 *
 * Usage example
 * \verbatim
 *
 * CellData<double> old_data; // Assume this is allocated and initialized
 * CellData<double> new_data; // Assume this is allocated
 *
 * auto old_array = old_data.getConstView<3>();
 * auto new_array = new_data.getView<3>();
 * const hier::Box& box = new_data.getBox();
 *
 * hier::parallel_for_all(box, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
 *    new_array(i, j, k) = old_array(i, j, k);
 * });
 *
 * \endverbatim
 */

template<int DIM, class TYPE>
struct ArrayView {};

template<class TYPE>
struct ArrayView<1, TYPE> : public RAJA::View<TYPE, detail::layout_traits::Layout1d>
{
   using Layout = detail::layout_traits::Layout1d;

   ArrayView(TYPE* data, const hier::Box& box) :
      RAJA::View<TYPE, Layout>(
         data,
         RAJA::make_permuted_offset_layout(
            std::array<RAJA::Index_type, 1>{ {box.lower()[0]} },
            std::array<RAJA::Index_type, 1>{ {box.upper()[0]+1} },
            RAJA::as_array<RAJA::PERM_I>::get())){}
};

template<class TYPE>
struct ArrayView<2, TYPE> : public RAJA::View<TYPE, detail::layout_traits::Layout2d>
{
   using Layout = detail::layout_traits::Layout2d;

   SAMRAI_INLINE ArrayView(TYPE* data, const hier::Box& box) :
      RAJA::View<TYPE, Layout>(
         data,
         RAJA::make_permuted_offset_layout(
            std::array<RAJA::Index_type, 2>{ {box.lower()[0], box.lower()[1]} },
            std::array<RAJA::Index_type, 2>{ {box.upper()[0]+1, box.upper()[1]+1} },
            RAJA::as_array<RAJA::PERM_JI>::get())){}
};

template<class TYPE>
struct ArrayView<3, TYPE> : public RAJA::View<TYPE, detail::layout_traits::Layout3d>
{
   using Layout = detail::layout_traits::Layout3d;

   SAMRAI_INLINE ArrayView(TYPE* data, const hier::Box& box) :
      RAJA::View<TYPE, Layout>(
         data,
         RAJA::make_permuted_offset_layout(
            std::array<RAJA::Index_type, 3>{ {box.lower()[0], box.lower()[1], box.lower()[2]} },
            std::array<RAJA::Index_type, 3>{ {box.upper()[0]+1, box.upper()[1]+1, box.upper()[2]+1} },
            RAJA::as_array<RAJA::PERM_KJI>::get())){};
};

template<class TYPE>
struct ArrayView<1, const TYPE> : public RAJA::View<const TYPE, detail::layout_traits::Layout1d>
{
   using Layout = detail::layout_traits::Layout1d;

   ArrayView(const TYPE* data, const hier::Box& box) :
      RAJA::View<const TYPE, Layout>(
         data,
         RAJA::make_permuted_offset_layout(
            std::array<RAJA::Index_type, 1>{ {box.lower()[0]} },
            std::array<RAJA::Index_type, 1>{ {box.upper()[0]+1} },
            RAJA::as_array<RAJA::PERM_I>::get())){}
};

template<class TYPE>
struct ArrayView<2, const TYPE> : public RAJA::View<const TYPE, detail::layout_traits::Layout2d>
{
   using Layout = detail::layout_traits::Layout2d;

   SAMRAI_INLINE ArrayView(const TYPE* data, const hier::Box& box) :
      RAJA::View<const TYPE, Layout>(
         data,
         RAJA::make_permuted_offset_layout(
            std::array<RAJA::Index_type, 2>{ {box.lower()[0], box.lower()[1]} },
            std::array<RAJA::Index_type, 2>{ {box.upper()[0]+1, box.upper()[1]+1} },
            RAJA::as_array<RAJA::PERM_JI>::get())){}
};


template<class TYPE>
struct ArrayView<3, const TYPE> : public RAJA::View<const TYPE, detail::layout_traits::Layout3d>
{
   using Layout = detail::layout_traits::Layout3d;

   SAMRAI_INLINE ArrayView(const TYPE* data, const hier::Box& box) :
      RAJA::View<const TYPE, Layout>(
         data,
         RAJA::make_permuted_offset_layout(
            std::array<RAJA::Index_type, 3>{ {box.lower()[0], box.lower()[1], box.lower()[2]} },
            std::array<RAJA::Index_type, 3>{ {box.upper()[0]+1, box.upper()[1]+1, box.upper()[2]+1} },
            RAJA::as_array<RAJA::PERM_KJI>::get())){};
};

}
}

#endif // HAVE_RAJA

#endif
