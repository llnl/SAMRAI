/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   $Description
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/pdat/MDA_Access.h"

#ifndef LACKS_EXPLICIT_TEMPLATE_INSTANTIATION

template class MDA_Access<double, 2, MDA_OrderColMajor<2> >;
template class MDA_Access<int, 2, MDA_OrderColMajor<2> >;

template class MDA_Access<double, 3, MDA_OrderColMajor<3> >;
template class MDA_Access<int, 3, MDA_OrderColMajor<3> >;

#endif
