/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Binary tree of Boxes for overlap searches.
 *
 ************************************************************************/

#ifndef included_hier_BoxTree
#define included_hier_BoxTree

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/tbox/Timer.h"

#include <vector>
#include <list>
#include <memory>

namespace SAMRAI {
namespace hier {

class BoxContainer;

/*!
 * @brief Utility sorting Boxes into tree-like form for finding
 * box overlaps.
 *
 * This class recursively splits a set of Boxes into tree-like
 * form and stores them for fast searches.  The recursive
 * splitting stops when the number of boxes in a leaf node of the tree
 * is less than a minimum number specified in the constructor.
 *
 * All boxes in a BoxTree must exist in the same index space.
 * This means that they must all have the same BlockId value.
 *
 * Overlap searches are done by
 * - hasOverlap()
 * - findOverlapBoxes()
 *
 * Except for two static methods and a destructor needed by shared_ptr,
 * the entire interface is private.
 */

class BoxTree
{
   friend class MultiblockBoxTree;

public:

   /*!
    * @brief Print statistics on number of constructor calls, tree
    * builds, tree searches, etc.
    *
    * This method is for developers to analyze performance.
    */
   static void
   printStatistics(
      const tbox::Dimension& dim);

   /*!
    * @brief Reset statistics on number of constructor calls, tree
    * builds, tree searches, etc.
    *
    * This method is for developers to analyze performance.
    */
   static void
   resetStatistics(
      const tbox::Dimension& dim);

   /*!
    * @brief Destructor.
    */
   ~BoxTree() = default;

private:

   BoxTree(
      const std::list<const Box *> boxes,
      int min_number = 10);

   /*!
    * @brief Constructs a BoxTree from set of Boxes.
    *
    * @param[in] dim
    *
    * @param[in] boxes
    *
    * @param[in] min_number Split up sets of boxes while the number of
    * boxes in a subset is greater than this value.  Setting to a
    * larger value tends to make tree building faster but tree
    * searching slower, and vice versa.  @b Default: 10
    *
    * @pre for each box in boxes, !box.empty()
    * @pre each box in boxes has a valid, identical BlockId
    */
   BoxTree(
      const tbox::Dimension& dim,
      const BoxContainer& boxes,
      int min_number = 10);

   /*!
    * @brief Constructor building an uninitialized object.
    *
    * Private as it is used only internally to create child trees.
    * The object can be initialized using generateTree().
    *
    * @param[in] dim
    */
   explicit BoxTree(
      const tbox::Dimension& dim);

   /*!
    * Default constructor is unimplemented and should not be used.
    */
   BoxTree() = delete;

   /*!
    * @brief Reset to uninitialized state.
    *
    * The dimension of boxes in the tree cannot be changed.
    *
    * Uninitialized trees can be initialized using generateTree().
    */
   void
   clear()
   {
      d_bounding_box.setEmpty();
      d_left_child.reset();
      d_right_child.reset();
      d_boxes.clear();
      d_center_child.reset();
   }

   /*!
    * @brief Check whether the tree has been initialized.
    *
    * Uninitialized trees can be initialized using generateTree().
    */
   bool
   isInitialized() const
   {
      return !d_bounding_box.empty();
   }

   //@{

   //! @name Access to state data

   /*!
    * @brief Return the dimension of the boxes in the tree.
    */
   const tbox::Dimension&
   getDim() const
   {
      return d_dim;
   }

   const BlockId&
   getBlockId() const
   {
      return d_block_id;
   }

   //@}

   //@{

   //! @name Overlap checks

   /*!
    * @brief Whether the given box has an overlap with Boxes in the
    * tree.
    *
    * @param[in] box The box is assumed to be in same index space as
    * those in the tree.
    *
    * @pre getDim() == box.getDim()
    */
   bool
   hasOverlap(
      const Box& box) const;

   /*!
    * @brief Find all boxes that overlap the given \b box.
    *
    * Analogous to findOverlapBoxes returning a BoxContainer
    * but avoids the copies.  If the returned overlap_boxes are used
    * in a context in which the BoxTree is constant there is no point
    * in incurring the cost of copying the tree's Boxes.  Just return
    * a vector of their addresses.
    *
    * @param[out] overlap_boxes Pointers to Boxes that overlap
    * with box.
    *
    * @param[in] box the specified box whose overlaps are requested.
    *
    * @param[in] recursive_call Disable logging of information except at first
    * call.
    *
    * @pre getDim() == box.getDim()
    * @pre box.getBlockId() == getBlockId()
    */
   void
   findOverlapBoxes(
      std::vector<const Box *>& overlap_boxes,
      const Box& box,
      bool recursive_call = false) const;

   /*!
    * @brief Find all boxes that overlap the given \b box.
    *
    * To avoid unneeded work, the output @b overlap_boxes container
    * is not emptied.  Overlapping Boxes are simply added.
    *
    * @param[out] overlap_boxes Boxes that overlap with box.  The ordered/
    * unordered state of this container is not changed from its state at
    * entry of this method
    *
    * @param[in] box the specified box whose overlaps are requested.
    *
    * @param[in] recursive_call Disable logging of information except at first
    * call.
    *
    * @pre getDim() == box.getDim()
    * @pre box.getBlockId() == getBlockId()
    */
   void
   findOverlapBoxes(
      BoxContainer& overlap_boxes,
      const Box& box,
      bool recursive_call = false) const;

   //@}

   /*!
    * @brief Private recursive function for generating the search tree.
    *
    * d_boxes is changed in the process (for efficiency reasons).
    * On output it will contain any boxes that are not assigned to a child
    * tree.
    *
    * The object is not cleared in this method.  If the object has
    * been initialized, it should be cleared before calling this
    * method.  @see clear().
    *
    * @param min_number  @b Default: 10
    */
   void
   privateGenerateTree(
      int min_number = 10);

   /*!
    * @brief Set up the child branches.
    *
    * This method is called after splitting the Boxes into the
    * left_boxes and right_boxes, with boxes straddling
    * the divider stored in d_boxes.  It generates
    * d_left_child, d_right_child and, if needed, d_center_child.
    *
    * @param[in] min_number
    *
    * @param[in,out] left_boxes
    *
    * @param[in,out] right_boxes
    */
   void
   setupChildren(
      const int min_number,
      std::list<const Box*>& left_boxes,
      std::list<const Box*>& right_boxes);

   /*!
    * @brief Set up static class members.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * @brief Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   /*!
    * @brief Dimension corresponds to the dimension of boxes in the
    * tree.
    */
   const tbox::Dimension d_dim;

   /*!
    * @brief Bounding box of all the Boxes in this tree.
    */
   Box d_bounding_box;

   /*!
    * @brief BlockId
    */
   BlockId d_block_id;

   /*!
    * std::shared_ptrs to familial boxes.
    */
   std::shared_ptr<BoxTree> d_left_child;
   std::shared_ptr<BoxTree> d_right_child;

   /*!
    * @brief A tree for Boxes that are not given to the left or
    * right children.
    */
   std::shared_ptr<BoxTree> d_center_child;

   /*!
    * @brief Boxes that are contained within the physical domain
    * that this tree represents.  When we have a small number of boxes
    * that do not warant the overhead of a child tree, the boxes go here.
    */
   std::list<const Box*> d_boxes;

   /*!
    * @brief Dimension along which the input box triples are
    * partitioned.
    */
   tbox::Dimension::dir_t d_partition_dir;

   /*
    * Timers are static to keep the objects light-weight.
    */
   static std::shared_ptr<tbox::Timer> t_build_tree[SAMRAI::MAX_DIM_VAL];
   static std::shared_ptr<tbox::Timer> t_search[SAMRAI::MAX_DIM_VAL];

   static unsigned int s_num_build[SAMRAI::MAX_DIM_VAL];
   static unsigned int s_num_generate[SAMRAI::MAX_DIM_VAL];
   static unsigned int s_num_duplicate[SAMRAI::MAX_DIM_VAL];
   static unsigned int s_num_search[SAMRAI::MAX_DIM_VAL];
   static unsigned int s_num_sorted_box[SAMRAI::MAX_DIM_VAL];
   static unsigned int s_num_found_box[SAMRAI::MAX_DIM_VAL];
   static unsigned int s_max_sorted_box[SAMRAI::MAX_DIM_VAL];
   static unsigned int s_max_found_box[SAMRAI::MAX_DIM_VAL];
   static unsigned int s_max_lin_search[SAMRAI::MAX_DIM_VAL];

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#endif
