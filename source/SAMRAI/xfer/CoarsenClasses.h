/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Simple structure for managing coarsening data in equivalence classes.
 *
 ************************************************************************/

#ifndef included_xfer_CoarsenClasses
#define included_xfer_CoarsenClasses

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/VariableFillPattern.h"

#include <iostream>
#include <list>
#include <vector>
#include <memory>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Maintain a collection of coarsen items and organize them
 * into equivalence classes.
 *
 * CoarsenClasses is used by the CoarsenSchedule and CoarsenAlgorithm
 * classes to manage coarsen data items that describe coarsening of
 * patch data between two levels in an AMR hierarchy.  Specifically, this
 * class organizes these items into equivalence clases, so that items are
 * grouped together if they have the same data communication dependencies.
 * See documentation for the method itemsAreEquivalent() for definition
 * of equivalence.
 */

class CoarsenClasses
{
public:
   /*!
    * @brief Nested class used to describe a coarsening operation
    * between patch data components on an AMR hierarchy.
    */
   class Data
   {
public:
      /*!
       * @brief Destination patch data component
       */
      int d_dst;

      /*!
       * @brief Source patch data component
       */
      int d_src;

      /*!
       * @brief Boolean flag that is set to true when it is desired that fine
       * data values have priority over coarse data values when data exists
       * at the same location in the mesh on levels with different resolutions.
       */
      bool d_fine_bdry_reps_var;

      /*!
       * If the coarsen operation requires data to be coarsened from the fine
       * level's ghost regions onto the coarse data representing the same
       * mesh space, then this IntVector tells how wide the ghost region to
       * coarsen will be.  It is represented in terms of the coarse level's
       * index space.
       */
      hier::IntVector d_gcw_to_coarsen;

      /*!
       * @brief Coarsening operator.
       */
      std::shared_ptr<hier::CoarsenOperator> d_opcoarsen;

      /*!
       * @brief Index of equivalence class where this item belongs.  All
       * items of the same equivalence class will have the same value.
       */
      int d_class_index;

      /*!
       * @brief An array index telling where this item sits in an array of
       * coarsen items.
       */
      int d_tag;

      /*!
       * @brief VariableFillPattern that can restrict the stencil of the data
       * coarsened by the CoarsenSchedule.
       */
      std::shared_ptr<VariableFillPattern> d_var_fill_pattern;

      /*!
       * @brief Constructor.
       *
       * @param[in] dim Dimension.
       */
      explicit Data(
         tbox::Dimension dim);

private:
      Data();  //not implemented
   };

   /*!
    * @brief The default constructor creates an empty array of coarsen classes.
    */
   CoarsenClasses();

   /*!
    * @brief The destructor destroys the coarsen data items owned
    * by this object.
    */
   ~CoarsenClasses();

   /*!
    * Return number of equivalence classes maintained by this object.
    */
   int
   getNumberOfEquivalenceClasses() const
   {
      return static_cast<int>(d_equivalence_class_indices.size());
   }

   /*!
    * @brief Return total number of coarsen items that have been registered
    * and stored in the CoarsenClasses object
    */
   int
   getNumberOfCoarsenItems() const
   {
      return d_num_coarsen_items;
   }

   /*!
    * @brief Get representative item for a given equivalence class index.
    *
    * @return Given an index of an existing equivalence class, one item
    * from that class is returned.
    *
    * @param[in] equiv_class_index
    *
    * @pre (equiv_class_index >= 0) &&
    *      (equiv_class_index < getNumberOfEquivalenceClasses())
    */
   const CoarsenClasses::Data&
   getClassRepresentative(
      int equiv_class_index) const
   {
      TBOX_ASSERT((equiv_class_index >= 0) &&
         (equiv_class_index < getNumberOfEquivalenceClasses()));
      return d_coarsen_classes_data_items[
                d_equivalence_class_indices[equiv_class_index].front()];
   }

   /*!
    * @brief Get a coarsen item from the array of all coarsen items held by
    * this object.
    *
    * The internal storage of the coarsen items held by this class is not
    * controlled by the user, so this method is intended for use when looping
    * over all of the items, from 0 to getNumberOfCoarsenItems()-1, or when
    * looping over the integers in the List obtained from getIterator().
    *
    * @return A coarsen classes data object identified by an integer id.
    *
    * @param[in] coarsen_item_array_id
    */
   CoarsenClasses::Data&
   getCoarsenItem(
      const int coarsen_item_array_id)
   {
      return d_coarsen_classes_data_items[coarsen_item_array_id];
   }

   /*!
    * @brief Return an iterator for the list of array ids corresponding to the
    * equivalence class with the given integer index.
    *
    * The number of quivalence classes can be determined via the
    * getNumberOfEquivalenceClasses() member function.  Valid integer
    * arguments are from 0 to getNumberOfEquivalenceClasses()-1.
    *
    * @note The list should not be modified through this iterator.
    *
    * @return The iterator iterates over a list of integers which are array
    * ids that can be passed into getCoarsenItem().  The array ids in a
    * single list all correspond to coarsen items in a single equivalence
    * class.
    *
    * @param[in] equiv_class_index
    *
    * @pre (equiv_class_index >= 0) &&
    *      (equiv_class_index < getNumberOfEquivalenceClasses())
    */
   std::list<int>::iterator
   getIterator(
      int equiv_class_index)
   {
      TBOX_ASSERT((equiv_class_index >= 0) &&
         (equiv_class_index < getNumberOfEquivalenceClasses()));
      return d_equivalence_class_indices[equiv_class_index].begin();
   }

   /*!
    * @brief Return an iterator for the list of array ids corresponding to the
    * equivalence class with the given integer index.
    *
    * The number of quivalence classes can be determined via the
    * getNumberOfEquivalenceClasses() member function.  Valid integer
    * arguments are from 0 to getNumberOfEquivalenceClasses()-1.
    *
    * @note The list should not be modified through this iterator.
    *
    * @return The iterator iterates over a list of integers which are array
    * ids that can be passed into getCoarsenItem().  The array ids in a
    * single list all correspond to coarsen items in a single equivalence
    * class.
    *
    * @param[in] equiv_class_index
    *
    * @pre (equiv_class_index >= 0) &&
    *      (equiv_class_index < getNumberOfEquivalenceClasses())
    */
   std::list<int>::iterator
   getIteratorEnd(
      int equiv_class_index)
   {
      TBOX_ASSERT((equiv_class_index >= 0) &&
         (equiv_class_index < getNumberOfEquivalenceClasses()));
      return d_equivalence_class_indices[equiv_class_index].end();
   }

   /*!
    * @brief Given a CoarsenClasses::Data object, insert it into the proper
    * equivalence class.
    *
    * If the item belongs in an existing equivalence class, it will be added
    * to that class. Otherwise, a new equivalence class will be created for
    * this item.  The integer class index in the data item will set to the
    * index of the equivalence class into which it is inserted.
    *
    * If a null patch descriptor argument is passed (or ommitted), the
    * descriptor associated with the variable database Singleton object will be
    * used.
    *
    * @param[in,out] data
    * @param[in] descriptor
    *
    * @pre itemIsValid(data, descriptor)
    */
   void
   insertEquivalenceClassItem(
      CoarsenClasses::Data& data,
      const std::shared_ptr<hier::PatchDescriptor>& descriptor =
         std::shared_ptr<hier::PatchDescriptor>());

   /*!
    * @brief Check coarsen data item for validity.
    *
    * A coarsen data item is invalid if any of its patch data ids are
    * negative, or if its source data does not have sufficient ghost width
    * for the stencil of the coarsen operator, or if copying from
    * the source to the destination data is not valid.
    *
    * An error will occur with a descriptive message if the item is invalid.
    *
    * If a null patch descriptor argument is passed (or ommitted), the
    * descriptor associated with the variable database Singleton object will
    * be used.
    *
    * @return True if the item is valid; else false.
    *
    * @param[in] data_item
    * @param[in] descriptor
    */
   bool
   itemIsValid(
      const CoarsenClasses::Data& data_item,
      const std::shared_ptr<hier::PatchDescriptor>& descriptor =
         std::shared_ptr<hier::PatchDescriptor>()) const;

   /*!
    * @brief Compare CoarsenClasses object with another CoarsenClasses object;
    *        return true if they match, else false.
    *
    * This method checks whether all equivalence classes match between this
    * CoarsenClasses object and the argument object.  To match, the number of
    * equivalence classes held by the objects must be the same and each
    * equivalence class in this object must match the class with the same
    * equivalence class number in the argument object.  Two classes match if
    * they have the same number of items and their representative items are
    * equialvent as defined by the method itemsAreEquivalent().
    *
    * If a null patch descriptor argument is passed (or ommitted), the
    * descriptor associated with the variable database Singleton object will
    * be used.
    *
    * @return true if test_classes matches this object.
    *
    * @param[in] test_classes  CoarsenClasses object to compare with this.
    * @param[in] descriptor
    */
   bool
   classesMatch(
      const std::shared_ptr<CoarsenClasses>& test_classes,
      const std::shared_ptr<hier::PatchDescriptor>& descriptor =
         std::shared_ptr<hier::PatchDescriptor>()) const;

   /*!
    * @brief Compare CoarsenClasses::Data objects for equivalence;
    *        return true if equivalent, else false.
    *
    * Two CoarsenClasses::Data objects are equivalent if and only if
    * the following conditions hold:
    *
    * <ul>
    *    <li> Each corresponding patch data component (d_dst and d_src)
    *         must have the same patch data type and ghost cell width.
    *    <li> d_fine_bdry_reps_var flag must have the same value for
    *         each object.
    *    <li> The coarsening operator ptr d_opcoarsen must be null or non-null
    *         for both objects.  If non-null, both operators must have the
    *         same stencil width.
    *    <li> The same variable fill pattern is used.
    * </ul>
    *
    * If a null patch descriptor argument is passed (or ommitted), the
    * descriptor associated with the variable database Singleton object will
    * be used.
    *
    * @return true if test_classes matches this object.
    *
    * @param[in] data1  CoarsenClasses::Data object to compare.
    * @param[in] data2  CoarsenClasses::Data object to compare.
    * @param[in] descriptor
    */
   bool
   itemsAreEquivalent(
      const CoarsenClasses::Data& data1,
      const CoarsenClasses::Data& data2,
      const std::shared_ptr<hier::PatchDescriptor>& descriptor =
         std::shared_ptr<hier::PatchDescriptor>()) const;

   /*!
    * @brief Get the size that has been allocated for the array storing coarsen
    * items.
    *
    * Note that this is not necessarily the same as the number of registered
    * coarsen items, which can be retrieved using getNumberOfCoarsenItems().
    * The coarsen item array is allocated to a default size and grown when
    * necessary or when increaseCoarsenItemArraySize() is called.
    */
   int
   getCoarsenItemArraySize() const
   {
      return static_cast<int>(d_coarsen_classes_data_items.size());
   }

   /*!
    * @brief Increase the allocated size of the array storing coarsen items.
    *
    * This should be used in cases where there is a large number of coarsen
    * items being registered with the CoarsenAlgorithm, to avoid frequent
    * resizing of the array.  If the size argument is less than the current
    * allocated size of the array, then the size of the array is not changed.
    *
    * @param[in] size
    * @param[in] dim
    */
   void
   increaseCoarsenItemArraySize(
      const int size,
      const tbox::Dimension& dim)
   {
      if (size > static_cast<int>(d_coarsen_classes_data_items.size())) {
         d_coarsen_classes_data_items.resize(size, Data(dim));
      }
   }

   /*!
    * @brief Print data for all coarsen items to the specified output stream.
    *
    * @param[out] stream
    */
   void
   printClassData(
      std::ostream& stream) const;

   /*!
    * @brief Print single coarsen item to the specified output stream.
    *
    * @param[out] stream
    * @param[in] data
    */
   void
   printCoarsenItem(
      std::ostream& stream,
      const CoarsenClasses::Data& data) const;

private:
   CoarsenClasses(
      const CoarsenClasses&);                   // not implemented
   void
   operator = (
      const CoarsenClasses&);                     // not implemented

   /*!
    * @brief Check two patch data items (with given descriptor indices)
    * to see whether they match.
    *
    * Two patch data items match if the are of the same patch data type and
    * have the same ghost width.
    *
    * @return true if consistent; false otherwise.
    *
    * @param[in] item_id1
    * @param[in] item_id2
    * @param[in] pd  descriptor
    */
   bool
   patchDataMatch(
      int item_id1,
      int item_id2,
      const std::shared_ptr<hier::PatchDescriptor>& pd) const;

   /*!
    * @brief Determine the equivalence class index of given
    * CoarsenClasses::Data object.
    *
    * The coarsen data item is compared with existing equivalence classes to
    * determine if it can be a member of any of them.
    *
    * @return If the item matches an existing equivalence class the integer
    * index for that equivalence class is returned.  Otherwise -1 is
    * returned.
    *
    * @param[in] data
    * @param[in] descriptor
    */
   int
   getEquivalenceClassIndex(
      const CoarsenClasses::Data& data,
      const std::shared_ptr<hier::PatchDescriptor>& descriptor =
         std::shared_ptr<hier::PatchDescriptor>()) const;

   /*!
    * The default length of the coarsen item array.
    */
   static int s_default_coarsen_item_array_size;

   /*!
    * The array of coarsen items.
    */
   std::vector<CoarsenClasses::Data> d_coarsen_classes_data_items;

   /*!
    * The array managing equivalence classes.  Each element of the array
    * represents one equivalence class.  Each List holds integers identifying
    * which items are part of an equivalence class.  The integers index into
    * the array d_coarsen_classes_data_items.
    */
   std::vector<std::list<int> > d_equivalence_class_indices;

   /*!
    * The number of coarsen items that have been registered.
    */
   int d_num_coarsen_items;

};

}
}

#endif
