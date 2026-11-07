/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   RAJA collective operations
 *
 ************************************************************************/

#ifndef included_tbox_Collectives
#define included_tbox_Collectives

#if defined(HAVE_RAJA)

#include "RAJA/RAJA.hpp"

#include "SAMRAI/tbox/ExecutionPolicy.h"
#include "SAMRAI/tbox/GPUUtilities.h"

namespace SAMRAI {
namespace tbox {

template <typename policy>
inline void
synchronize() {}

#if defined(HAVE_CUDA)
template<>
inline void
synchronize<policy::parallel>()
{
   RAJA::synchronize<RAJA::cuda_synchronize>();
}
#elif defined(HAVE_HIP)
template<>
inline void
synchronize<policy::parallel>()
{
   RAJA::synchronize<RAJA::hip_synchronize>();
}
#endif

inline void
parallel_synchronize() {
   GPUUtilities::parallel_synchronize();
}

// Reductions see https://raja.readthedocs.io/en/master/feature/reduction.html for these options
enum class Reduction {
   Sum,
   Min,
   Max,
   MinLoc,
   MaxLoc
};

template<typename Policy, Reduction R, typename TYPE = double>
struct reduction_variable;

template<typename Policy, typename TYPE>
struct reduction_variable<Policy, Reduction::Sum, TYPE> {
   using type = RAJA::ReduceSum<typename detail::policy_traits<Policy>::ReductionPolicy, TYPE>;
};

template<typename Policy, typename TYPE>
struct reduction_variable<Policy, Reduction::Min, TYPE> {
   using type = RAJA::ReduceMin<typename detail::policy_traits<Policy>::ReductionPolicy, TYPE>;
};

template<typename Policy, typename TYPE>
struct reduction_variable<Policy, Reduction::Max, TYPE> {
   using type = RAJA::ReduceMax<typename detail::policy_traits<Policy>::ReductionPolicy, TYPE>;
};

template<typename Policy, typename TYPE>
struct reduction_variable<Policy, Reduction::MinLoc, TYPE> {
   using type = RAJA::ReduceMinLoc<typename detail::policy_traits<Policy>::ReductionPolicy, TYPE>;
};

template<typename Policy, typename TYPE>
struct reduction_variable<Policy, Reduction::MaxLoc, TYPE> {
   using type = RAJA::ReduceMaxLoc<typename detail::policy_traits<Policy>::ReductionPolicy, TYPE>;
};

template<Reduction R, typename TYPE = double>
using parallel_reduction_variable = reduction_variable<tbox::policy::parallel, R, TYPE>;

template<typename Policy, Reduction R, typename TYPE = double>
using reduction_variable_t = typename reduction_variable<Policy, R, TYPE>::type;

template<Reduction R, typename TYPE = double>
using parallel_reduction_variable_t = typename parallel_reduction_variable<R, TYPE>::type;

}
}
#endif
#endif
