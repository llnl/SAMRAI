/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface for the AMR Index object
 *
 ************************************************************************/
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

namespace SAMRAI {
namespace hier {

std::array<std::optional<Index>, SAMRAI::MAX_DIM_VAL> Index::s_zeros{};
std::array<std::optional<Index>, SAMRAI::MAX_DIM_VAL> Index::s_ones{};
std::array<std::optional<Index>, SAMRAI::MAX_DIM_VAL> Index::s_mins{};
std::array<std::optional<Index>, SAMRAI::MAX_DIM_VAL> Index::s_maxs{};

tbox::StartupShutdownManager::Handler
Index::s_initialize_finalize_handler(
   Index::initializeCallback,
   0,
   0,
   Index::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

Index::Index(
   const tbox::Dimension& dim) noexcept :
   d_dim(dim),
   d_index{}
{
#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (unsigned int i = 0; i < SAMRAI::MAX_DIM_VAL; ++i) {
      d_index[i] = tbox::MathUtilities<int>::getMin();
   }
#endif
}

Index::Index(
   const tbox::Dimension& dim,
   int value) noexcept :
   d_dim(dim),
   d_index{}
{
   for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
      d_index[i] = value;
   }

#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (unsigned int i = d_dim.getValue(); i < SAMRAI::MAX_DIM_VAL; ++i) {
      d_index[i] = tbox::MathUtilities<int>::getMin();
   }
#endif
}

Index::Index(
   const int i,
   const int j) noexcept :
   d_dim(2),
   d_index{}
{
   TBOX_DIM_ASSERT(tbox::Dimension::getMaxDimension() >= tbox::Dimension(2));

   d_index[0] = i;
   if (SAMRAI::MAX_DIM_VAL > 1) {
      d_index[1] = j;
   }
}

Index::Index(
   const int i,
   const int j,
   const int k) noexcept :
   d_dim(3),
   d_index{}
{
   TBOX_DIM_ASSERT(tbox::Dimension::getMaxDimension() >= tbox::Dimension(3));

   d_index[0] = i;
   if (SAMRAI::MAX_DIM_VAL > 1) {
      d_index[1] = j;
   }

   if (SAMRAI::MAX_DIM_VAL > 2) {
      d_index[2] = k;
   }

}

Index::Index(
   const std::vector<int>& a) noexcept :
   d_dim(static_cast<unsigned short>(a.size())),
   d_index{}
{
   TBOX_ASSERT(a.size() > 0 && a.size() <= SAMRAI::MAX_DIM_VAL);
   for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
      d_index[i] = a[i];
   }

#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (unsigned int i = d_dim.getValue(); i < SAMRAI::MAX_DIM_VAL; ++i) {
      d_index[i] = tbox::MathUtilities<int>::getMin();
   }
#endif
}

Index::Index(
   const tbox::Dimension& dim,
   const int array[]) noexcept :
   d_dim(dim),
   d_index{}
{
   for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
      d_index[i] = array[i];
   }
}

Index::Index(
   const Index& rhs) noexcept :
   d_dim(rhs.d_dim),
   d_index{}
{
   for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
      d_index[i] = rhs.d_index[i];
   }
}

Index::Index(
   const IntVector& rhs) noexcept :
   d_dim(rhs.getDim()),
   d_index{}
{
   TBOX_ASSERT(rhs.getNumBlocks() == 1);
   for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
      d_index[i] = rhs[i];
   }
}

void
Index::initializeCallback()
{
   for (unsigned short d = 0; d < SAMRAI::MAX_DIM_VAL; ++d) {
      tbox::Dimension dim(d + 1);
      s_zeros[d].emplace(dim, 0);
      s_ones[d].emplace(dim, 1);
      s_mins[d].emplace(dim, tbox::MathUtilities<int>::getMin());
      s_maxs[d].emplace(dim, tbox::MathUtilities<int>::getMax());
   }
}

void
Index::finalizeCallback()
{
   for (unsigned short d = 0; d < SAMRAI::MAX_DIM_VAL; ++d) {
      s_zeros[d].reset();
      s_ones[d].reset();
      s_mins[d].reset();
      s_maxs[d].reset();
   }
}

std::istream&
operator >> (
   std::istream& s,
   Index& rhs)
{
   while (s.get() != '(') ;

   for (int i = 0; i < rhs.getDim().getValue(); ++i) {
      s >> rhs(i);
      if (i < rhs.getDim().getValue() - 1)
         while (s.get() != ',') ;
   }

   while (s.get() != ')') ;

   return s;
}

std::ostream& operator << (
   std::ostream& s,
   const Index& rhs)
{
   s << '(';

   for (int i = 0; i < rhs.getDim().getValue(); ++i) {
      s << rhs(i);
      if (i < rhs.getDim().getValue() - 1)
         s << ",";
   }
   s << ')';

   return s;
}


}
}
