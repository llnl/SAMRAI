/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Simple structure for managing refinement data in equivalence classes.
 *
 ************************************************************************/

#ifndef included_xfer_RefineClasses
#define included_xfer_RefineClasses

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/TimeInterpolateOperator.h"
#include "SAMRAI/xfer/VariableFillPattern.h"

#include <iostream>
#include <list>
#include <vector>
#include <memory>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Maintain a collection of refine items and organize them
 * into equivalence classes.
 *
 * RefineClasses is used by the RefineSchedule and RefineAlgorithm
 * classes to manage refinement data items that describe communication
 * of patch data on an AMR hierarchy.  Specifically, this class organizes
 * these items into equivalence clases, so that items are grouped
 * together if they have the same data communication dependencies.
 * See documentation for the method itemsAreEquivalent() for definition
 * of equivalence.
 */

class RefineClasses
{
public:
   /*!
    * @brief Data structure used to describe a refinement operation
    * between patch data components on an AMR hierarchy.
    */
   struct Data {
      /*!
       * @brief Destination patch data component
       */
      int d_dst;

      /*!
       * @brief Source patch data component
       */
      int d_src;

      /*!
       * @brief Patch data component for source data at the old time in
       * a time interpolation operation.
       */
      int d_src_told;

      /*!
       * @brief Patch data component for source data at the new time in
       * a time interpolation operation.
       */
      int d_src_tnew;

      /*!
       * @brief Scratch patch data component
       */
      int d_scratch;

      std::vector<int> d_work;

      /*!
       * @brief Boolean flag that is set to true when it is desired that fine
       * data values have priority over coarse data values when data exists
       * at the same location in the mesh on levels with different resolutions
       * (e.g., nodes that are coincident on consecutive mesh levels).
       */
      bool d_fine_bdry_reps_var;

      /*!
       * @brief Boolean flag telling if this item uses time interpolation.
       */
      bool d_time_interpolate;

      /*!
       * @brief Refinement operator
       */
      std::shared_ptr<hier::RefineOperator> d_oprefine;

      /*!
       * @brief Time interpolation operator
       */
      std::shared_ptr<hier::TimeInterpolateOperator> d_optime;

      /*!
       * @brief Index of equivalence class where this item belongs.  All
       * items of the same equivalence class will have the same value.
       */
      int d_class_index;

      /*!
       * @brief An array index telling where this item sits in an array of
       * refine items.
       */
      int d_tag;

      /*!
       * @brief VariableFillPattern that can restrict the stencil of the data
       * filled by the RefineSchedule.
       */
      std::shared_ptr<VariableFillPattern> d_var_fill_pattern;
   };

   /*!
    * @brief The constructor creates an empty array of refine classes.
    */
   RefineClasses();

   /*!
    * @brief The destructor destroys the refinement data items owned
    * by this object.
    */
   ~RefineClasses();

   /*!
    * @brief Return number of equivalence classes maintained by this object.
    */
   int
   getNumberOfEquivalenceClasses() const
   {
      return static_cast<int>(d_equivalence_class_indices.size());
   }

   /*!
    * @brief Return total number of refine items that have been registered and
    * stored in the RefineClasses object
    */
   int
   getNumberOfRefineItems() const
   {
      return static_cast<int>(d_refine_classes_data_items.size());
   }

   /*!
    * @brief Get representative item for a given equivalence class index.
    *
    * @return Given index of an existing equivalence class, one item
    * from that class is returned.
    *
    * @param[in] equiv_class_index
    *
    * @pre (equiv_class_index >= 0) &&
    *      (equiv_class_index < getNumberOfEquivalenceClasses())
    */
   const RefineClasses::Data&
   getClassRepresentative(
      int equiv_class_index) const
   {
      TBOX_ASSERT((equiv_class_index >= 0) &&
         (equiv_class_index < getNumberOfEquivalenceClasses()));
      return d_refine_classes_data_items[
                d_equivalence_class_indices[equiv_class_index].front()];
   }

   /*!
    * @brief Get a refine item from the array of all refine items held by
    * this object.
    *
    * The internal storage of the refine items held by this class is not
    * controlled by the user, so this method is intended for use when looping
    * over all of the items, from 0 to getNumberOfRefineItems()-1, or when
    * looping over the integers in the List obtained from getIterator().
    *
    * @return A refine classes data object identified by an integer id.
    *
    * @param[in] refine_item_array_id
    */
   RefineClasses::Data&
   getRefineItem(
      const int refine_item_array_id)
   {
      return d_refine_classes_data_items[refine_item_array_id];
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
    * ids that can be passed into getRefineItem().  The array ids in a
    * single list all correspond to refine items in a single equivalence
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
    * ids that can be passed into getRefineItem().  The array ids in a
    * single list all correspond to refine items in a single equivalence
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
    * @brief Given a RefineClasses::Data object, insert it into the proper
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
    * @param[in,out] data_item
    * @param[in] descriptor
    *
    * @pre itemIsValid(data, descriptor)
    */
   void
   insertEquivalenceClassItem(
      RefineClasses::Data& data_item,
      const std::shared_ptr<hier::PatchDescriptor>& descriptor =
         std::shared_ptr<hier::PatchDescriptor>());

   /*!
    * @brief Check refine data item for validity.
    *
    * A refine data item is invalid if any of its patch data integer ids are
    * negative, or if its scratch data does not have sufficient ghost width
    * for the stencil of the refine operator or the fill pattern, or if
    * it is not a valid operation to copy from source data to scratch data or
    * from scratch data to destination data, or when time interpolation is used
    * the old and new time patch data enties are either undefined or their
    * types do not match the source data type.
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
      const RefineClasses::Data& data_item,
      const std::shared_ptr<hier::PatchDescriptor>& descriptor =
         std::shared_ptr<hier::PatchDescriptor>()) const;

   /*!
    * @brief Compare RefineClasses object with another RefineClasses object;
    *        return true if they match, else false.
    *
    * This method checks whether all equivalence classes match between this
    * RefineClasses object and the argument object.  To match, the number of
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
    * @param[in] test_classes  RefineClasses object to compare with this.
    * @param[in] descriptor
    */
   bool
   classesMatch(
      const std::shared_ptr<RefineClasses>& test_classes,
      const std::shared_ptr<hier::PatchDescriptor>& descriptor =
         std::shared_ptr<hier::PatchDescriptor>()) const;

   /*!
    * @brief Compare RefineClasses::Data objects for equivalence;
    *        return true if equivalent, else false.
    *
    * Two RefineClasses::Data objects are equivalent if and only if
    * the following conditions hold:
    *
    * <ul>
    *    <li> Each corresponding patch data component (d_dst, d_src, etc.)
    *         must have the same patch data type and ghost cell width.
    *    <li> The d_time_interpolate flag must be true or false for both
    *         objects. If true, each corresponding patch data component
    *         (d_src_told, d_src_tnew) must have the same data type and
    *         ghost cell width.  Also, the time interpolate operators
    *         must be the same type.
    *    <li> d_fine_bdry_reps_var flag must have the same value for
    *         each object.
    *    <li> The refinement operator ptr d_oprefine must be null or non-null
    *         for both objects.  If non-null, both operators must have the
    *         same stencil width.
    *    <li> The type of the d_var_fill_pattern must be the same for both
    *         objects.
    * </ul>
    *
    * If a null patch descriptor argument is passed (or ommitted), the
    * descriptor associated with the variable database Singleton object will
    * be used.
    *
    * @return true if test_classes matches this object.
    *
    * @param[in] data1  RefineClasses::Data object to compare.
    * @param[in] data2  RefineClasses::Data object to compare.
    * @param[in] descriptor
    */
   bool
   itemsAreEquivalent(
      const RefineClasses::Data& data1,
      const RefineClasses::Data& data2,
      const std::shared_ptr<hier::PatchDescriptor>& descriptor =
         std::shared_ptr<hier::PatchDescriptor>()) const;

   /*!
    * @brief Increase the allocated size of the array storing refine items.
    *
    * This should be used in cases where there is a large number of refine
    * items being registered with the RefineAlgorithm, to avoid frequent
    * resizing of the array.  If the size argument is less than the current
    * allocated size of the array, then the size of the array is not changed.
    *
    * @param[in] size
    */
   void
   increaseRefineItemArraySize(
      const int size)
   {
      if (size > static_cast<int>(d_refine_classes_data_items.size())) {
         d_refine_classes_data_items.resize(size);
      }
   }

   /*!
    * @brief Print data for all refine items to the specified output stream.
    *
    * @param[out] stream
    */
   void
   printClassData(
      std::ostream& stream) const;

   /*!
    * @brief Print single refine item to the specified output stream.
    *
    * @param[out] stream
    * @param[in] data
    */
   void
   printRefineItem(
      std::ostream& stream,
      const RefineClasses::Data& data) const;

private:
   RefineClasses(
      const RefineClasses&);            // not implemented
   RefineClasses&
   operator = (
      const RefineClasses&);                     // not implemented

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
    * @brief Determine the equivalence class index of given RefineClasses::Data
    * object.
    *
    * The refine data item is compared with existing equivalence classes to
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
      const RefineClasses::Data& data,
      const std::shared_ptr<hier::PatchDescriptor>& descriptor) const;

   /*!
    * The default length of the refine item array.
    */
   static int s_default_refine_item_array_size;

   /*!
    * The array of refine items.
    */
   std::vector<Data> d_refine_classes_data_items;

   /*!
    * The array managing equivalence classes.  Each element of the array
    * represents one equivalence class.  Each List holds integers identifying
    * which items are part of an equivalence class.  The integers index into
    * the array d_refine_classes_data_items.
    */
   std::vector<std::list<int> > d_equivalence_class_indices;

};

}
}

#endif
