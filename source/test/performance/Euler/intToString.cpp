/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   $Description
 *
 ************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "intToString.h"

#include <string>

#include <sstream>
#include <iomanip>


string intToString(
   int i,
   int min_length)
{
   ostringstream co;
   co << setw(min_length) << setfill('0') << i;
   return string(co.str());
}
