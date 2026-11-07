/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Specifications for the scalar Poisson equation
 *
 ************************************************************************/
#include "SAMRAI/solv/PoissonSpecifications.h"

namespace SAMRAI {
namespace solv {

/*
 *******************************************************************
 * Default constructor
 *******************************************************************
 */

PoissonSpecifications::PoissonSpecifications(
   const std::string& object_name):d_object_name(object_name),
   d_D_id(-1),
   d_D_constant(1.0),
   d_C_zero(true),
   d_C_id(-1),
   d_C_constant(0.0)
{
}

/*
 *******************************************************************
 * Copy constructor
 *******************************************************************
 */

PoissonSpecifications::PoissonSpecifications(
   const std::string& object_name,
   const PoissonSpecifications& r):d_object_name(object_name),
   d_D_id(r.d_D_id),
   d_D_constant(r.d_D_constant),
   d_C_zero(r.d_C_zero),
   d_C_id(r.d_C_id),
   d_C_constant(r.d_C_constant)
{
}

/*
 *******************************************************************
 * Destructor (does nothing).
 *******************************************************************
 */
PoissonSpecifications::~PoissonSpecifications()
{
}

void
PoissonSpecifications::printClassData(
   std::ostream& stream) const
{
   stream << "PoissonSpecifications " << d_object_name << "\n"
          << "   D is ";
   if (d_D_id != -1) {
      stream << "variable with patch id " << d_D_id << "\n";
   } else {
      stream << "constant with value " << d_D_constant << "\n";
   }
   stream << "   C is ";
   if (d_C_zero) {
      stream << "zero\n";
   } else if (d_C_id != -1) {
      stream << "variable with patch id " << d_C_id << "\n";
   } else {
      stream << "constant with value " << d_C_constant << "\n";
   }
}

} // namespace solv
} // namespace SAMRAI
