/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Light weight array class.
 *
 ************************************************************************/
#ifndef included_ArrayDataAccess_h
#define included_ArrayDataAccess_h

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/MDA_Access.h"

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Utility for wrapping data from ArrayData class in an
 * MDA_Access or MDA_AccessConst object.
 *
 * This, and other classes in the folder array_access are meant
 * only for certain SAMRAI development efforts.  Please do not
 * use in user code.
 */
class ArrayDataAccess
{

public:
   /*!
    * @brief Create an MDA_Access version of an ArrayData object.
    */
   template<int DIM, class TYPE>
   static MDA_Access<TYPE, DIM, MDA_OrderColMajor<DIM> >
   access(
      ArrayData<TYPE>& array_data,
      int depth = 0) {
      return MDA_Access<TYPE, DIM, MDA_OrderColMajor<DIM> >(
                array_data.getPointer(depth),
                &array_data.getBox().lower()[0],
                &array_data.getBox().upper()[0]);
   }

   /*!
    * @brief Create an MDA_AccessConst version of a const ArrayData object.
    */
   template<int DIM, class TYPE>
   static MDA_AccessConst<TYPE, DIM, MDA_OrderColMajor<DIM> >
   access(
      const ArrayData<TYPE>& array_data,
      int depth = 0) {
      return MDA_AccessConst<TYPE, DIM, MDA_OrderColMajor<DIM> >(
                array_data.getPointer(depth),
                &array_data.getBox().lower()[0],
                &array_data.getBox().upper()[0]);
   }

};

}
}

#include "SAMRAI/pdat/ArrayDataAccess.cpp"

#endif
