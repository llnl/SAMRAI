/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Class to record statistics during program execution.
 *
 ************************************************************************/

#ifndef included_hier_ForAll
#define included_hier_ForAll

#include "SAMRAI/SAMRAI_config.h"

#if defined(HAVE_RAJA)
#include "RAJA/RAJA.hpp"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/ExecutionPolicy.h"
#include "SAMRAI/tbox/KernelFuser.h"

#include <type_traits>
#include <tuple>
#include <cstdlib>  // for std::size_t

namespace SAMRAI
{
namespace hier
{

/*!
 * Three looping structures for use with RAJA policies are provided;
 *
 * parallel_for_all() uses a default parallel policy based on the
 * configuration used to build SAMRAI.  If SAMRAI is confugured with CUDA,
 * the loop kernel will be launched on the GPU.  If cuda is not used
 * then the loop will run on the CPU.
 *
 * host_parallel_for_all() is provided to define loops that always run
 * on the host CPU even when configured with CUDA.  This allows RAJA loop
 * kernels to be implemented to operate on data that the calling code has
 * chosen to keep on the CPU.
 * 
 * Both parallel_for_all() and host_parallel_for_all() are overloaded to
 * provide interfaces that loop over a mesh defined by a hier::Box of
 * a Dimension from 1 to 3.  To use it, provide the Box and a lambda function
 * that specifies the integer loop indices.  The number of indices must
 * equal the Dimension of the Box.  The SAMRAI_HOST_DEVICE macro is used to
 * set the default policy for parallel_for_all().  SAMRAI_HOST_DEVICE is
 * not recommended to be used in host_parallel_for_all().
 * 
 * \verbatim
 *
 * hier::parallel_for_all(box, [=] SAMRAI_HOST_DEVICE(int i, int j, int k)
 * {
 *    // loop interior
 * });
 *
 * hier::host_parallel_for_all(box, [=] (int i, int j, int k)
 * {
 *    // loop interior
 * });
 *
 * \endverbatim
 *
 * for_all<policy>() uses a template to allow the calling code to choose
 * a custom RAJA policy.  Usage is similar, as it also loops over a hier::Box.
 *
 * \verbatim
 *
 * hier::for_all<policy>(box, [=] (int i, int j, int k)
 * {
 *    // loop interior
 * });
 *
 * \endverbatim
 *
 */

namespace detail
{

template <typename T>
struct function_traits : function_traits<decltype(&T::operator())> {
};

// free function
template <typename R, typename... Args>
struct function_traits<R(Args...)> {
   using result_type = R;
   using argument_types = std::tuple<Args...>;
   enum { argument_count = sizeof...(Args) };
};

// pointer to function
template <typename R, typename... Args>
struct function_traits<R (*)(Args...)> {
   using result_type = R;
   using argument_types = std::tuple<Args...>;
   enum { argument_count = sizeof...(Args) };
};

// member function
template <typename T, typename R, typename... Args>
struct function_traits<R (T::*)(Args...)> {
   using result_type = R;
   using argument_types = std::tuple<Args...>;
   enum { argument_count = sizeof...(Args) };
};

// const member function
template <typename T, typename R, typename... Args>
struct function_traits<R (T::*)(Args...) const> {
   using result_type = R;
   using argument_types = std::tuple<Args...>;
   enum { argument_count = sizeof...(Args) };
};

inline RAJA::RangeSegment make_range(const hier::Index& ifirst, const hier::Index& ilast, std::size_t index)
{
   return RAJA::RangeSegment(ifirst(index), ilast(index) + 1);
}

template <int ArgumentCount>
struct for_all {
};

template <>
struct for_all<1> {
   template <typename Policy, typename LoopBody,
             typename std::enable_if<std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<typename tbox::detail::policy_traits<Policy>::Policy1d>(
          RAJA::make_tuple(make_range(ifirst, ilast, 0)),
          body);
   }

   template <typename Policy, typename LoopBody,
             typename std::enable_if<!std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<Policy>(
          RAJA::make_tuple(make_range(ifirst, ilast, 0)),
          body);
   }

   template <typename Policy, typename LoopBody,
             typename std::enable_if<std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(tbox::KernelFuser* fuser, const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      if (fuser == nullptr) {
         RAJA::kernel<typename tbox::detail::policy_traits<Policy>::Policy1d>(
             RAJA::make_tuple(make_range(ifirst, ilast, 0)),
             body);
      } else {
         fuser->enqueue(ifirst(0), ilast(0), body);
      }
   }

   template <typename Policy, typename LoopBody,
             typename std::enable_if<!std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(tbox::KernelFuser* fuser, const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      if (fuser == nullptr) {
         RAJA::kernel<Policy>(
             RAJA::make_tuple(make_range(ifirst, ilast, 0)),
             body);
      } else {
         fuser->enqueue(ifirst(0), ilast(0), body);
      }
   }
};


template <>
struct for_all<2> {
   template <typename Policy, typename LoopBody,
             typename std::enable_if<std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<typename tbox::detail::policy_traits<Policy>::Policy2d>(
          RAJA::make_tuple(make_range(ifirst, ilast, 0),
                           make_range(ifirst, ilast, 1)),
          body);
   }

   template <typename Policy, typename LoopBody,
             typename std::enable_if<!std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<Policy>(
          RAJA::make_tuple(make_range(ifirst, ilast, 0),
                           make_range(ifirst, ilast, 1)),
          body);
   }

   template <typename Policy, typename LoopBody,
             typename std::enable_if<std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(tbox::KernelFuser* fuser, const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      NULL_USE(fuser);
      RAJA::kernel<typename tbox::detail::policy_traits<Policy>::Policy2d>(
          RAJA::make_tuple(make_range(ifirst, ilast, 0),
                           make_range(ifirst, ilast, 1)),
          body);
   }

   template <typename Policy, typename LoopBody,
             typename std::enable_if<!std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(tbox::KernelFuser* fuser, const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      NULL_USE(fuser);
      RAJA::kernel<Policy>(
          RAJA::make_tuple(make_range(ifirst, ilast, 0),
                           make_range(ifirst, ilast, 1)),
          body);
   }
};

template <>
struct for_all<3> {
   template <typename Policy, typename LoopBody,
             typename std::enable_if<std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<typename tbox::detail::policy_traits<Policy>::Policy3d>(
          RAJA::make_tuple(make_range(ifirst, ilast, 0),
                           make_range(ifirst, ilast, 1),
                           make_range(ifirst, ilast, 2)),
          body);
   }

   template <typename Policy, typename LoopBody,
             typename std::enable_if<!std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      RAJA::kernel<Policy>(
          RAJA::make_tuple(make_range(ifirst, ilast, 0),
                           make_range(ifirst, ilast, 1),
                           make_range(ifirst, ilast, 2)),
          body);
   }

   template <typename Policy, typename LoopBody,
             typename std::enable_if<std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(tbox::KernelFuser* fuser, const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      NULL_USE(fuser);
      RAJA::kernel<typename tbox::detail::policy_traits<Policy>::Policy3d>(
          RAJA::make_tuple(make_range(ifirst, ilast, 0),
                           make_range(ifirst, ilast, 1),
                           make_range(ifirst, ilast, 2)),
          body);
   }

   template <typename Policy, typename LoopBody,
             typename std::enable_if<!std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(tbox::KernelFuser* fuser, const hier::Index& ifirst, const hier::Index& ilast, LoopBody body)
   {
      NULL_USE(fuser);
      RAJA::kernel<Policy>(
          RAJA::make_tuple(make_range(ifirst, ilast, 0),
                           make_range(ifirst, ilast, 1),
                           make_range(ifirst, ilast, 2)),
          body);
   }
};

}  // namespace detail

// does NOT include end
template <typename Policy, typename LoopBody,
          typename std::enable_if<std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
inline void for_all(int begin, int end, LoopBody body)
{
   RAJA::forall<typename tbox::detail::policy_traits<Policy>::Policy>(RAJA::RangeSegment(begin, end), body);
}

template <typename Policy, typename LoopBody,
          typename std::enable_if<std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
inline void for_all(tbox::KernelFuser* fuser, int begin, int end, LoopBody body)
{
   if (fuser == nullptr) {
      RAJA::forall<typename tbox::detail::policy_traits<Policy>::Policy>(RAJA::RangeSegment(begin, end), body);
   } else {
      fuser->enqueue(begin, end, body);
   }
}

template <typename Policy, typename LoopBody,
          typename std::enable_if<!std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
inline void for_all(int begin, int end, LoopBody body)
{
   RAJA::forall<Policy>(RAJA::RangeSegment(begin, end), body);
}

template <typename Policy, typename LoopBody,
          typename std::enable_if<!std::is_base_of<tbox::policy::base, Policy>::value, int>::type = 0>
inline void for_all(tbox::KernelFuser* fuser, int begin, int end, LoopBody body)
{
   if (fuser == nullptr) {
      RAJA::forall<Policy>(RAJA::RangeSegment(begin, end), body);
   } else {
      fuser->enqueue(begin, end, body);
   }
}

// does NOT include end
template <typename LoopBody>
inline void parallel_for_all(int begin, int end, LoopBody body)
{
   for_all<tbox::policy::parallel>(begin, end, body);
}

template <typename LoopBody>
inline void parallel_for_all(tbox::KernelFuser* fuser, int begin, int end, LoopBody body)
{
   if (fuser == nullptr) {
      for_all<tbox::policy::parallel>(begin, end, body);
   } else {
      for_all<tbox::policy::parallel>(fuser, begin, end, body);
   }
}

template <typename LoopBody>
inline void host_parallel_for_all(int begin, int end, LoopBody body)
{
   for_all<tbox::policy::host_parallel>(begin, end, body);
}

template <typename Policy, typename LoopBody>
inline void for_all(const hier::Box& box, const int dim, LoopBody body)
{
   for_all<Policy>(box.lower()(dim), box.upper()(dim) + 1, body);
}


template <typename Policy, typename LoopBody>
inline void for_all(tbox::KernelFuser* fuser, const hier::Box& box, const int dim, LoopBody body)
{
   for_all<Policy>(fuser, box.lower()(dim), box.upper()(dim) + 1, body);
}

template <typename LoopBody>
inline void parallel_for_all(const hier::Box& box, const int dim, LoopBody body)
{
   for_all<tbox::policy::parallel>(box.lower()(dim), box.upper()(dim) + 1, body);
}

template <typename LoopBody>
inline void parallel_for_all(tbox::KernelFuser* fuser, const hier::Box& box, const int dim, LoopBody body)
{
   for_all<tbox::policy::parallel>(fuser, box.lower()(dim), box.upper()(dim) + 1, body);
}

template <typename LoopBody>
inline void host_parallel_for_all(const hier::Box& box, const int dim, LoopBody body)
{
   for_all<tbox::policy::host_parallel>(box.lower()(dim), box.upper()(dim) + 1, body);
}

template <typename Policy, typename LoopBody>
inline void for_all(const hier::Box& box, LoopBody body)
{
   constexpr int arg_count = detail::function_traits<LoopBody>::argument_count;
   detail::for_all<arg_count>::template eval<Policy>(box.lower(), box.upper(), body);
}

template <typename Policy, typename LoopBody>
inline void for_all(tbox::KernelFuser* fuser, const hier::Box& box, LoopBody body)
{
   if (fuser == nullptr) {
      for_all<Policy,LoopBody>(box, body);
   } else {
      constexpr int arg_count = detail::function_traits<LoopBody>::argument_count;
      detail::for_all<arg_count>::template eval<Policy>(fuser, box.lower(), box.upper(), body);
   }
}

template <typename LoopBody>
inline void parallel_for_all(const hier::Box& box, LoopBody body)
{
   for_all<tbox::policy::parallel>(box, body);
}

template <typename LoopBody>
inline void parallel_for_all(tbox::KernelFuser* fuser, const hier::Box& box, LoopBody body)
{
   for_all<tbox::policy::parallel>(fuser, box, body);
}

template <typename LoopBody>
inline void host_parallel_for_all(const hier::Box& box, LoopBody body)
{
   for_all<tbox::policy::host_parallel>(box, body);
}


}  // namespace hier
}  // namespace SAMRAI

#endif

#endif  // included_hier_ForAll
