/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface for using efficient communication tree.
 *
 ************************************************************************/
#ifndef included_tbox_RankTreeStrategy
#define included_tbox_RankTreeStrategy

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/RankGroup.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Strategy pattern for getting information on MPI ranks
 * aranged on a user-defined tree.
 *
 * The purpose of a rank tree is to facilitate "collective-like"
 * communications where information is communicated up or down the
 * tree along its edges.  It gives the parent rank, children ranks,
 * etc., for a node in the tree.
 *
 * To improve communication performance, the tree should map well to
 * the communication hardware.  That is, processes on each ends of an
 * edge should be close to each other on the communication network.
 * Implementations of this strategy should consider the network
 * topology when deciding how to lay out the process on the tree.
 */
class RankTreeStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   RankTreeStrategy();

   /*!
    * @brief Destructor.
    */
   virtual ~RankTreeStrategy();

   /*!
    * @brief Set up the tree.
    *
    * Set up the tree for the processors in the given RankGroup.
    * Prepare to provide tree data for the given rank.
    *
    * @param[in] rank_group
    *
    * @param[in] my_rank The rank whose parent and children are
    * sought, usually the local process.
    */
   virtual void
   setupTree(
      const RankGroup& rank_group,
      int my_rank) = 0;

   /*!
    * @brief Access the rank used in initialization.
    */
   virtual int
   getRank() const = 0;

   /*!
    * @brief Access the parent rank.
    *
    * @returns parent rank or getInvalidRank() if is root of the tree.
    */
   virtual int
   getParentRank() const = 0;

   /*!
    * @brief Access a child rank.
    *
    * @param child_number If outside [0, getNumberOfChildren()),
    * the invalid rankd should be returned.
    *
    * @return Child's rank, or getInvalidRank() if child_number
    * is outside [0, getNumberOfChildren()).
    */
   virtual int
   getChildRank(
      unsigned int child_number) const = 0;

   /*!
    * @brief Return the number of children.
    */
   virtual unsigned int
   getNumberOfChildren() const = 0;

   /*!
    * @brief Return the child number.
    *
    * Return child number or invalidChildNumber() if is root of the
    * tree.
    */
   virtual unsigned int
   getChildNumber() const = 0;

   /*!
    * @brief Return the degree of the tree (the maximum number of
    * children each node may have).
    */
   virtual unsigned int
   getDegree() const = 0;

   /*!
    * @brief Return the generation number.
    *
    * The generation number starts at zero, at the root, and increases
    * by one each generation.
    */
   virtual unsigned int
   getGenerationNumber() const = 0;

   /*!
    * @brief Return the rank of the root of the tree.
    */
   virtual int
   getRootRank() const = 0;

   /*!
    * @brief Return whether the local process is the root of its tree.
    */
   bool isRoot() const {
      return getRank() == getRootRank();
   }

   /*!
    * @brief What this class considers an invalid rank.
    *
    * When a parent or child does not exist, this value is returned for
    * the rank.
    */
   static int getInvalidRank() {
      return s_invalid_rank;
   }

   /*!
    * @brief What this class considers an invalid child number.
    *
    * The root of the tree is not anyone's child, so it should return
    * this for getChildNumber();.
    */
   static unsigned int getInvalidChildNumber() {
      return s_invalid_child_number;
   }

private:
   // Unimplemented copy constructor.
   RankTreeStrategy(
      const RankTreeStrategy& other);

   // Unimplemented assignment opperator.
   RankTreeStrategy&
   operator = (
      const RankTreeStrategy& rhs);

   static const int s_invalid_rank;
   static const unsigned int s_invalid_child_number;

};

}
}

#endif  // included_tbox_RankTreeStrategy
