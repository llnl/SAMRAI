/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A n-dimensional integer vector
 *
 ************************************************************************/
#include "SAMRAI/hier/IntVector.h"

#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"


namespace SAMRAI {
namespace hier {

std::array<std::optional<IntVector>, SAMRAI::MAX_DIM_VAL> IntVector::s_zeros{};
std::array<std::optional<IntVector>, SAMRAI::MAX_DIM_VAL> IntVector::s_ones{};

tbox::StartupShutdownManager::Handler
IntVector::s_initialize_finalize_handler(
   IntVector::initializeCallback,
   0,
   0,
   IntVector::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);


/*
 * *************************************************************************
 * Constructors
 * *************************************************************************
 */

IntVector::IntVector(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_num_blocks(1),
   d_vector(dim.getValue())
{
#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
      d_vector[i] = tbox::MathUtilities<int>::getMin();
   }
#endif
}

IntVector::IntVector(
   size_t num_blocks,
   const tbox::Dimension& dim):
   d_dim(dim),
   d_num_blocks(num_blocks),
   d_vector(dim.getValue()*num_blocks)
{
   TBOX_ASSERT(num_blocks >=1);
#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (unsigned int i = 0; i < num_blocks*dim.getValue(); ++i) {
      d_vector[i] = tbox::MathUtilities<int>::getMin();
   }
#endif
}

IntVector::IntVector(
   const tbox::Dimension& dim,
   int value,
   size_t num_blocks):
   d_dim(dim),
   d_num_blocks(num_blocks),
   d_vector(dim.getValue()*num_blocks, value)
{
   TBOX_ASSERT(num_blocks >=1);
}

IntVector::IntVector(
   const std::vector<int>& vec,
   size_t num_blocks):
   d_dim(static_cast<unsigned short>(vec.size())),
   d_num_blocks(num_blocks),
   d_vector(vec.size()*num_blocks)
{
   TBOX_ASSERT(vec.size() >= 1);
   for (BlockId::block_t b = 0; b < num_blocks; ++b) {
      unsigned int offset = b*d_dim.getValue();
      for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
         d_vector[offset + i] = vec[i];
      }
   }
}

IntVector::IntVector(
   const tbox::Dimension& dim,
   const int array[],
   size_t num_blocks):
   d_dim(dim),
   d_num_blocks(num_blocks),
   d_vector(dim.getValue()*num_blocks)
{
   for (BlockId::block_t b = 0; b < num_blocks; ++b) {
      unsigned int offset = b*d_dim.getValue();
      for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
         d_vector[offset + i] = array[i];
      }
   }
}

IntVector::IntVector(
   const IntVector& rhs,
   size_t num_blocks):
   d_dim(rhs.getDim()),
   d_num_blocks(num_blocks),
   d_vector(rhs.getDim().getValue() * num_blocks)
{
   TBOX_ASSERT(d_num_blocks >= 1);
   TBOX_ASSERT(rhs.d_num_blocks == d_num_blocks || rhs.d_num_blocks == 1); 
   if (rhs.d_num_blocks == 1 && d_num_blocks != 1) {
      for (BlockId::block_t b = 0; b < d_num_blocks; ++b) {
         unsigned int offset = b*d_dim.getValue();
         for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
            d_vector[offset + i] = rhs.d_vector[i];
         }
      }
   } else {
      d_vector = rhs.d_vector;
   }
}

IntVector::IntVector(
   const Index& rhs,
   size_t num_blocks):
   d_dim(rhs.getDim()),
   d_num_blocks(num_blocks),
   d_vector(rhs.getDim().getValue() * num_blocks)
{
   TBOX_ASSERT(d_num_blocks >= 1);
   for (BlockId::block_t b = 0; b < num_blocks; ++b) {
      unsigned int offset = b*d_dim.getValue();
      for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
         d_vector[offset + i] = rhs[i];
      } 
   }
}

/*
 * *************************************************************************
 * Assignment
 * *************************************************************************
 */
IntVector&
IntVector::operator = (
   const Index& rhs)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
   if (d_num_blocks != 1) {
      d_num_blocks = 1;
      d_vector.resize(d_dim.getValue());
   }

   for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
      d_vector[i] = rhs[i];
   }
   return *this;
}

/*
 * *************************************************************************
 * Streaming I/O
 * *************************************************************************
 */
std::istream&
operator >> (
   std::istream& s,
   IntVector& rhs)
{
   for (BlockId::block_t b = 0; b < rhs.getNumBlocks(); ++b) {
      while (s.get() != '(') ;
      for (unsigned int i = 0; i < rhs.getDim().getValue(); ++i) {
         s >> rhs(b,i);
         if (static_cast<int>(i) < rhs.getDim().getValue() - 1)
            while (s.get() != ',') ;
      }
      while (s.get() != ')') ;
   }

   return s;
}

std::ostream& operator << (
   std::ostream& s,
   const IntVector& rhs)
{

   for (BlockId::block_t b = 0; b < rhs.getNumBlocks(); ++b) {
      s << '(';
      for (unsigned int i = 0; i < rhs.getDim().getValue(); ++i) {
         s << rhs(b,i);
         if (static_cast<int>(i) < rhs.getDim().getValue() - 1)
            s << ",";
      }
      s << ')';
   }

   return s;
}

/*
 * *************************************************************************
 * Write/read for restart
 * *************************************************************************
 */
void
IntVector::putToRestart(
   tbox::Database& restart_db,
   const std::string& name) const
{
   std::shared_ptr<tbox::Database> intvec_db =
      restart_db.putDatabase(name);
   intvec_db->putInteger("d_num_blocks", static_cast<int>(d_num_blocks));
   intvec_db->putIntegerVector("d_vector",
                               d_vector);

}

void
IntVector::getFromRestart(
   tbox::Database& restart_db,
   const std::string& name)
{
   std::shared_ptr<tbox::Database> intvec_db =
      restart_db.getDatabase(name);

   d_num_blocks = static_cast<size_t>(intvec_db->getInteger("d_num_blocks"));
   d_vector = intvec_db->getIntegerVector("d_vector"); 

   TBOX_ASSERT(d_num_blocks * d_dim.getValue() == d_vector.size());

}

/*
 *************************************************************************
 * Sort the values of the given IntVector from smallest to largest value.
 *************************************************************************
 */
void
IntVector::sortIntVector(
   const IntVector& values)
{
   TBOX_ASSERT(d_dim == values.getDim());
   d_num_blocks = values.getNumBlocks();
   d_vector.resize(values.d_vector.size());

   for (BlockId::block_t b = 0; b < d_num_blocks; ++b ) {
      unsigned int offset = b*d_dim.getValue();
      for (unsigned int d = 0; d < d_dim.getValue(); ++d) {
         d_vector[offset + d] = static_cast<int>(d);
      }
      for (unsigned int d0 = 0;
           d0 < static_cast<unsigned int>(d_dim.getValue() - 1); ++d0) {
         for (unsigned int d1 = d0 + 1; d1 < d_dim.getValue(); ++d1) {
            unsigned int v0 = static_cast<unsigned int>(d_vector[offset + d0]);
            unsigned int v1 = static_cast<unsigned int>(d_vector[offset + d1]);
            if (values(b, v0) > values(b, v1)) {
               int tmp_d = d_vector[offset + d0];
               d_vector[offset + d0] = d_vector[offset + d1];
               d_vector[offset + d1] = tmp_d;
            }
         }
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      for (unsigned int d = 0;
           d < static_cast<unsigned int>(d_dim.getValue() - 1); ++d) {
         unsigned int v0 = static_cast<unsigned int>(d_vector[offset + d]);
         unsigned int v1 = static_cast<unsigned int>(d_vector[offset + d + 1]);
         TBOX_ASSERT(values(b, v0) <= values(b, v1));
      }
#endif
   }
}

/*
 * *************************************************************************
 * Callback routines
 * *************************************************************************
 */
void
IntVector::initializeCallback()
{
   for (unsigned short d = 0; d < SAMRAI::MAX_DIM_VAL; ++d) {
      tbox::Dimension dim(d + 1);
      s_zeros[d].emplace(dim, 0);
      s_ones[d].emplace(dim, 1);
   }
}

void
IntVector::finalizeCallback()
{
   for (int d = 0; d < SAMRAI::MAX_DIM_VAL; ++d) {
      s_zeros[d].reset();
      s_ones[d].reset();
   }
}

}
}
