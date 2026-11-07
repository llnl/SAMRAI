/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:
 *
 ************************************************************************/

#ifndef included_pdat_ArrayDataAccess_C
#define included_pdat_ArrayDataAccess_C

#include "SAMRAI/pdat/ArrayDataAccess.h"

namespace SAMRAI {
namespace pdat {

template<int DIM>
MDA_Access<double, DIM, MDA_OrderColMajor<DIM> >
access(
   ArrayData<double>& array_data,
   int depth)
{
   MDA_Access<double, DIM, MDA_OrderColMajor<DIM> > r(
      array_data.getPointer(depth),
      &array_data.getBox().lower()[0],
      &array_data.getBox().upper()[0]);
   return r;
}

template<int DIM>
const MDA_AccessConst<double, DIM, MDA_OrderColMajor<DIM> >
access(
   const ArrayData<double>& array_data,
   int depth)
{
   MDA_AccessConst<double, DIM, MDA_OrderColMajor<DIM> > r(
      array_data.getPointer(depth),
      &array_data.getBox().lower()[0],
      &array_data.getBox().upper()[0]);
   return r;
}

}
}
#endif
