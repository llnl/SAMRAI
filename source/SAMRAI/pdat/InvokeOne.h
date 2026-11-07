/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Copy operation on single array data elements templated on data type
 *
 ************************************************************************/

#ifndef included_pdat_InvokeOne
#define included_pdat_InvokeOne

template<typename F1, typename F2, typename F3, typename F4, typename ... Args>
auto invokeOneOfFour(F1 f1, F2 /*f2*/, F3 /*f3*/, F4 /*f4*/, Args && ... args) -> decltype(f1(args...))
{
    return f1(args...);
}

template<typename F1, typename F2, typename F3, typename F4, typename ... Args>
auto invokeOneOfFour(F1 /*f1*/, F2 f2, F3 /*f3*/, F4 /*f4*/, Args && ... args) -> decltype(f2(args...))
{
    return f2(args...);
}

template<typename F1, typename F2, typename F3, typename F4, typename ... Args>
auto invokeOneOfFour(F1 /*f1*/, F2 /*f2*/, F3 f3, F4 /*f4*/, Args && ... args) -> decltype(f3(args...))
{
    return f3(args...);
}

template<typename F1, typename F2, typename F3, typename F4, typename ... Args>
auto invokeOneOfFour(F1 /*f1*/, F2 /*f2*/, F3 /*f3*/, F4 f4, Args && ... args) -> decltype(f4(args...))
{
    return f4(args...);
}


template<typename F1, typename F2, typename F3, typename ... Args>
auto invokeOneOfThree(F1 f1, F2 /*f2*/, F3 /*f3*/, Args && ... args) -> decltype(f1(args...))
{
    return f1(args...);
}

template<typename F1, typename F2, typename F3, typename ... Args>
auto invokeOneOfThree(F1 /*f1*/, F2 f2, F3 /*f3*/, Args && ... args) -> decltype(f2(args...))
{
    return f2(args...);
}

template<typename F1, typename F2, typename F3, typename ... Args>
auto invokeOneOfThree(F1 /*f1*/, F2 /*f2*/, F3 f3, Args && ... args) -> decltype(f3(args...))
{
    return f3(args...);
}

namespace SAMRAI {
namespace pdat {


}
}


#endif
