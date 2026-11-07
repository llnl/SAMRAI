/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   pdat
 *
 ************************************************************************/

#ifndef included_pdat_SparseData_h
#define included_pdat_SparseData_h

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/pdat/IntegerAttributeId.h"
#include "SAMRAI/pdat/DoubleAttributeId.h"

#include <list>
#include <string>
#include <vector>
#include <unordered_map>

namespace SAMRAI {
namespace pdat {

// forward declarations
template<typename BOX_GEOMETRY>
class SparseDataIterator;
template<typename BOX_GEOMETRY>
class SparseDataAttributeIterator;

template<typename BOX_GEOMETRY>
std::ostream&
operator << (
   std::ostream& out,
   SparseDataIterator<BOX_GEOMETRY>& sparse_data_iterator);

template<typename BOX_GEOMETRY>
std::ostream&
operator << (
   std::ostream& out,
   SparseDataAttributeIterator<BOX_GEOMETRY>& attr_iterator);

/*!
 * @brief SparseData is used for storing sparse data.
 *
 * SparseData is a PatchData type that allows storage of a collection
 * of attributes at an Index location.  Each SparseData object is initialized
 * with the attributes that it will support (unique to each instance).
 * Each Index can have any number of attribute collections associated with it.
 *
 * Attributes names are registered when the object is created.  Once created,
 * no additional attributes names can be registered for that instance.  The
 * names have a unique associated identifier which allows for operator[]
 * access to the individual attribute.
 *
 * Access to the Indexed elements is done through the SparseData iterator.
 * Clients should use the SparseData<BOX_GEOMETRY>::iterator interface.
 *
 * Access to an Index's attribute list is done through the AttributeIterator.
 * Individual elements within an attribute are then accessed through the
 * operator[].
 *
 * <code>
 * SparseData<BOX_GEOMETRY>::iterator index_iter = sparse_data.begin();
 * SparseData<BOX_GEOMETRY>::iterator index_iter_end = sparse_data.end();
 *
 * for (; index_iter != index_iter_end; ++index_iter) {
 *    // operator over all Index elements in the SparseData collection.
 * }
 *
 * SparseData<BOX_GEOMETRY::AttributeIterator attr_list_iter =
 *   sparse_data.begin(idx);
 * SparseData<BOX_GEOMETRY::AttributeIterator attr_list_iter_end =
 *   sparse_data.end(idx);
 *
 * for ( ; attr_list_iter != attr_list_iter_end; ++attr_list_iter) {
 *    // do something with attr_list_iter[ dbl_attr_id ];
 *    // do somethign with attr_list_iter[ int_attr_id ];
 * }
 * </code>
 *
 * You can also create the AttributeIterator directly as such:
 *
 * <code>
 * SparseData<BOX_GEOMETRY>::AttributeIterator attributes(sparse_data, index);
 * </code>
 *
 * Access to the individual elements is done through the iterator's operator[].
 *
 * For copy operations (between similar SparseData (PatchData) types,
 * see copy() and copy2() methods.
 *
 * To remove a single element from the SparseData list, you must use
 * the iterator interface to erase it.
 *
 * To erase all elements in the list, use the clear() method.
 *
 * Since SparseData is derived from hier::PatchData, its interface conforms
 * to the standard interface that PatchData defines.
 *
 * TEMPLATE PARAMETERS
 *
 * The BOX_GEOMETRY template parameter defines the geometry.  It must have
 * a nested class name Overlap that implements the following two methods:
 *
 * <li> getSourceOffset
 * <li> getDestinationBoxList
 *
 * As with all other PatchData types, SparseData objects are created by the
 * SparseDataFactory.
 *
 * @see hier::PatchData
 * @see SparseDataFactory
 * @see SparseDataVariable
 */
template<typename BOX_GEOMETRY>
class SparseData:public hier::PatchData
{
private:
   class Attributes;
   struct index_hash;

   template <class T>
   static
   void hash_combine(std::size_t& seed, const T& v);

   template <class T>
   void to_lower(T& input);

public:
   /*!
    * @brief Iterator access through SparseData<BOX_GEOMETRY>::iterator
    */
   typedef SparseDataIterator<BOX_GEOMETRY> iterator;

   /*!
    * @brief AttributeIterator access through
    * SparseData<BOX_GEOMETRY>::AttributeIterator
    */
   typedef SparseDataAttributeIterator<BOX_GEOMETRY> AttributeIterator;

   /*!
    * @brief Construct a SparseData object.
    *
    * SparseData objects are constructed as normal PatchData objects,
    * and adds the addition of two string vectors which specify the names
    * of the attributes for both double values and integer values that
    * will be stored in Index attributes in this instance.  Once the
    * instance is constructed, no other attribute (names) can be added.
    *
    * @param [in] box Describes the interior of the index space
    * @param [in] ghosts Describes the ghost nodes in each coordinate
    *                    direction.
    * @param [in] dbl_attributes The double (named) attributes
    * @param [in] int_attributes The integer (named) attributes
    *
    * @pre box.getDim() == ghosts.getDim()
    */
   SparseData(
      const hier::Box& box,
      const hier::IntVector& ghosts,
      const std::vector<std::string>& dbl_names,
      const std::vector<std::string>& int_names);

   /*!
    * @brief Destructor
    *
    * Destruction will clear the list and all associated memory.
    */
   virtual ~SparseData();

   /*!
    * @brief A fast copy between source and destination.
    *
    * All data is copied from the source into the destination where there is
    * overlap in the index space.
    *
    * @param [in] src The source PatchData from which to copy.
    *
    * @pre getDim() == ghosts.getDim()
    * @pre dynamic_cast<const SparseData<BOX_GEOMETRY> *>(&src) != 0
    */
   void
   copy(
      const hier::PatchData& src);

   /*!
    * @brief A fast copy between source and destination.
    *
    * All data is copied from the source into the destination where there is
    * overlap in the index space.  This copy does not change the state of
    * <tt>this</tt>
    *
    * @param [in] dst The destination PatchData
    *
    * @pre getDim() == dst.getDim()
    */
   void
   copy2(
      hier::PatchData& dst) const;

   /*!
    * @brief Copy data from the source using the designated overlap descriptor.
    *
    * @param [in] src The source PatchData from which to copy.
    * @param [in] overlap The overlap description
    *
    * @pre getDim() == ghosts.getDim()
    * @pre dynamic_cast<const SparseData<BOX_GEOMETRY> *>(&src) != 0
    * @pre dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap) != 0
    */
   void
   copy(
      const hier::PatchData& src,
      const hier::BoxOverlap& overlap);

   /*!
    * @brief Copy data to the destination using the overlap descriptor.
    *
    * All data is copied from <tt>this</tt> to the destination, without
    * changing the state of <tt>this</tt>.
    * @param [in] dst The destination PatchData
    * @param [in] overlap The overlap description
    *
    * @pre getDim() == dst.getDim()
    */
   void
   copy2(
      hier::PatchData& dst,
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Return true if the patch data object can estimate the
    * stream size required to fit its data using only index
    * space information (i.e., a box); else false.
    */
   bool
   canEstimateStreamSizeFromBox() const;

   /*!
    * @brief Calculate the number of bytes needed to stream the data.
    *
    * @param [in] overlap
    *
    * @pre dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap) != 0
    */
   size_t
   getDataStreamSize(
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Pack data residing on the specified index space.
    *
    * @param [out] stream The message stream
    * @param [in] overlap
    *
    * @pre dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap)
    */
   void
   packStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Unpack data from the message stream
    *
    * @param [in,out] stream The message stream
    * @param [in] overlap
    *
    * @pre dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap) != 0
    */
   void
   unpackStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap);

   /*!
    * @brief Obtain the specific information for this SparseData object
    *        from the restart Database.
    *
    * The parent PatchData part of this class will call this function
    * to get the parts of it specified within this concrete child.  It
    * will check the version number of this instance to ensure that
    * it is valid.
    *
    * @pre restart_db
    */
   void
   getFromRestart(
      const std::shared_ptr<tbox::Database>& restart_db);

   /*!
    * @brief Write out the specialized data to the restart Database.
    *
    * A version number is written out as well, in order to maintain
    * validity across runs.
    *
    * @pre restart_db
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /*!
    * @brief Returns the attribute ID associated with the named
    *        attribute.
    *
    * The attribute ID is unique for this instance.
    *
    * ASSERTIONS:
    *    An assertion will occur if double attributes have not been
    * registered.
    *
    * ERRORS: A TBOX_ERROR occurs if the attribute does not
    * exist in the list of registered attributes.
    */
   const DoubleAttributeId
   getDblAttributeId(
      const std::string& attribute) const;

   /*!
    * @brief Returns the attribute ID associated with the named
    *        attribute.
    *
    * The attribute ID is unique for this instance.
    *
    * ASSERTIONS:
    *    An assertion will occur if integer attributes have not been
    * registered.
    *
    * ERRORS: A TBOX_ERROR occurs if the attribute does not
    * exist in the list of registered attributes.
    */
   const IntegerAttributeId
   getIntAttributeId(
      const std::string& attribute) const;

   /*!
    * @brief Returns true if there are no elements within the SparseData.
    *
    * ASSERTIONS:  This must be a valid (created) SparseData object.
    */
   bool
   empty();

   /*!
    * @brief Registers a new set of Index-Attributes items to this
    *        SparseData.
    *
    * ASSERTIONS:
    *    This must be a valid/initialized object.
    *    The attribute value passed must exist in the set of registered
    * attributes.
    *
    * @return The iterator pointing to the index item just added
    */
   iterator
   registerIndex(
      const hier::Index& index);

   /*!
    * @brief Remove this Index and its associated attributes
    * from the object.
    *
    * Remove will invalidate iterators. The proper way to use remove in
    * a loop with iterators is as follows:
    *
    * <code>
    * while (iter != sparseDataObject->end()) {
    *    if ( <condition to test> ) {
    *       sparseDataObject->remove(iter);
    *    } else {
    *       ++iter;
    *    }
    * }
    * </code>
    *
    * "remove" will automatically increment the iterator such that the
    * next test of the loop will occur correctly.
    */
   void
   remove(
      iterator& iter);

   /*!
    * @brief Erases all elements within the object.
    *
    * NOTE:  If this is the last reference to the items added with
    * addItem, the underlying item will be cleaned up as well.
    */
   void
   clear();

   /*!
    * @brief Provides the number of Index elements contained within this
    * SparseData.
    *
    * ASSERTIONS:  This must be a valid (created) SparseData object.
    */
   int
   size();

   /*!
    * @brief Returns true if the id is valid.
    */
   bool
   isValid(
      const DoubleAttributeId& id) const;

   /*!
    * @brief Returns true if the id is valid.
    */
   bool
   isValid(
      const IntegerAttributeId& id) const;

   /*!
    * @brief Provides an iterator to the first SparseData item.
    */
   iterator
   begin();

   /*!
    * @brief Provides an iterator to indicate the end of the SparseData items.
    *
    * This iterator is a special value, and should never be considered to
    * contain any valid data.
    */
   iterator
   end();

   /*!
    * @brief Provides an iterator to the first attribute item in
    * the list associated with the Index.
    *
    */
   AttributeIterator
   begin(
      const hier::Index& index);

   /*!
    * @brief Provides an iterator to the last attribute item in
    * the list associated with the Index.
    */
   AttributeIterator
   end(
      const hier::Index& index);

   /*!
    * @brief Print the names of the attributes and their ID
    */
   void
   printNames(
      std::ostream& out) const;

   /*!
    * @brief Print the Index-Attributes for this sparse data object.
    */
   void
   printAttributes(
      std::ostream& out) const;

   /*!
    * @brief equality operator
    */
   bool
   operator == (
      const SparseData<BOX_GEOMETRY>& rhs) const;
   bool
   operator != (
      const SparseData<BOX_GEOMETRY>& rhs) const;

   // friends.
   friend class SparseDataIterator<BOX_GEOMETRY>;
   friend class SparseDataAttributeIterator<BOX_GEOMETRY>;
private:
   // Internal typedefs
   typedef std::list<Attributes> AttributeList;
   typedef std::unordered_map<
      hier::Index, AttributeList, index_hash> IndexMap;
   typedef std::unordered_map<
      std::string, DoubleAttributeId> DoubleAttrNameMap;
   typedef std::unordered_map<
      std::string, IntegerAttributeId> IntAttrNameMap;

   /*
    * Copy c'tor and assignment operator are private to prevent
    * the compiler from generating a default
    */
   SparseData(
      const SparseData& rhs);
   SparseData&
   operator = (
      const SparseData& rhs);

   /*
    * The dimension
    */
   const tbox::Dimension d_dim;

   /*
    * The map of index to attributes.  See typedef above
    */
   IndexMap d_index_to_attribute_map;

   /*
    * Registered name to ID maps.  See typedef above
    */
   DoubleAttrNameMap d_dbl_names;
   IntAttrNameMap d_int_names;

   // cached values.
   int d_dbl_attr_size;
   int d_int_attr_size;

   /*
    * Unique version number for this patch data for restart purposes.
    */
   static const int PDAT_SPARSEDATA_VERSION;

   /*
    * Invalid attributes ID.
    */
   static const int INVALID_ID;

   /**********************************************************************
    * Private methods for this class
    *********************************************************************/

   /*
    * Get the attribute list for this index.
    *
    * ASSERTION: Dimensions of the <tt>index</tt> and this must be the same.
    */
   //std::list<Attributes>& _get(const hier::Index& index) const;
   AttributeList&
   _get(
      const hier::Index& index) const;

   /*
    * Add this item to the index to attribute map
    */
   void
   _add(
      const typename IndexMap::const_iterator& item_to_add);

   /*
    * iterate through the index elements in this object, and if
    * the index is contained within the box, remove it from the
    * d_list.
    */
   void
   _removeInsideBox(
      const hier::Box& box);

   /*
    * The index hash function for adding elements to the std::unordered_map's
    * buckets.
    */
   struct index_hash {
      std::size_t operator () (const hier::Index& index) const;
   };

   /*
    * Internal attribute class.  Clients of SparseData will never use this
    * directly.
    *
    * This is intentionally kept as simple an implementation as possible.
    * We will add complexity if it becomes necessary.
    */
   class Attributes
   {
public:
      /**********************************************************************
      * Constructors, assignment, destructor
      **********************************************************************/
      Attributes(
         const int dsize,
         const int isize);
      Attributes(
         const Attributes& other);
      Attributes&
      operator = (
         const Attributes& rhs);
      ~Attributes();

      /**********************************************************************
      * modifiers
      **********************************************************************/
      void
      add(
         const double * dvals,
         const int * ivals);
      void
      add(
         const std::vector<double>& dvals,
         const std::vector<int>& ivals);

      /**********************************************************************
      * non-modifying operations
      **********************************************************************/
      const double *
      getDoubleAttributes() const;
      const int *
      getIntAttributes() const;

      /**********************************************************************
      * operators (modifying and non-modifying
      **********************************************************************/
      double&
      operator [] (
         const DoubleAttributeId& id);
      const double&
      operator [] (
         const DoubleAttributeId& id) const;
      int&
      operator [] (
         const IntegerAttributeId& id);
      const int&
      operator [] (
         const IntegerAttributeId& id) const;

      bool
      operator == (
         const Attributes& rhs) const;

      /**********************************************************************
      * output
      **********************************************************************/
      void
      printAttributes(
         std::ostream& out) const;
private:
      // internal typedefs
      typedef std::vector<double>::const_iterator dbl_iterator;
      typedef std::vector<int>::const_iterator int_iterator;

      // data members
      std::vector<double> d_dbl_attrs;
      std::vector<int> d_int_attrs;
   };  // end class Attributes

};

/**********************************************************************
 * Iterate over all sparse data items.
 *********************************************************************/
template<typename BOX_GEOMETRY>
class SparseDataIterator
{
public:
   template<typename T>
   friend
   std::ostream&
   operator << (
      std::ostream& out,
      SparseDataIterator<T>& sparse_data_iterator);

   /*!
    * @brief Constructs a default SparseDataIterator
    */
   SparseDataIterator();

   /*!
    * @brief Constructs a SparseDataIterator from a SparseData<BOX_GEOMETRY>
    * object.
    *
    * @param [in] sparse_data the SparseData oject
    */
   explicit SparseDataIterator(
      SparseData<BOX_GEOMETRY>& sparse_data);

   /*!
    * @brief Constructs a SparseDataIterator from a SparseData<BOX_GEOMETRY>
    * object.
    *
    * @param [in] sparse_data the SparseData oject
    */
   explicit SparseDataIterator(
      SparseData<BOX_GEOMETRY> * sparse_data);

   /*!
    * @brief Copy constructor
    *
    * @param [in] other
    */
   SparseDataIterator(
      const SparseDataIterator& other);

   /*!
    * @brief Assignment operator
    *
    * @param [in] rhs
    */
   SparseDataIterator&
   operator = (
      const SparseDataIterator& rhs);

   /*!
    * @brief  Destructor
    */
   ~SparseDataIterator();

   /*!
    * @brief Equality operator
    *
    * @param[in] rhs
    */
   bool
   operator == (
      const SparseDataIterator& rhs) const;

   /*!
    * @brief Inequality operator
    *
    * @param [in] rhs
    */
   bool
   operator != (
      const SparseDataIterator& rhs) const;

   /*!
    * @brief pre-increment operator
    */
   SparseDataIterator&
   operator ++ ();

   /*!
    * @brief post-increment operator
    */
   SparseDataIterator
   operator ++ (
      int);

   /*!
    * @brief returns the index of this iterator
    */
   const hier::Index&
   getIndex() const;

   /*!
    * @brief Insert a collection of attributes into the Attribute list
    * of the current Index element.
    *
    * NOTE:  The order of the elements is assumed to be the order of the
    * attribute names provided during registration of the Index element.
    *
    * @param [in] dvals The array of double attribute values
    * @param [in] ivals The array of integer attribute values
    */
   void
   insert(
      const double* dvals,
      const int* ivals);

   /*!
    * @brief Insert a collection of attributes into the Attribute list
    * of the current Index element.
    *
    * NOTE:  The order of the elements is assumed to be the order of the
    * attribute names provided during registration of the Index element.
    *
    * @param [in] dvals The vector of double attribute values
    * @param [in] ivals The vector of integer attribute values
    */
   void
   insert(
      const std::vector<double>& dvals,
      const std::vector<int>& ivals);

   /*!
    * @brief Checks equality of the contents of this interator against
    * the <tt>other</tt> iterator.
    *
    * @param [in] other
    * @return True if the contents are the same.
    */
   bool
   equals(
      const SparseDataIterator& other) const;

   /*!
    * @brief Move the attributes of this iterator to the given Index's list.
    *
    * PRECONDITIONS:  The Index parameter given must exist in this sparse data
    * object.
    *
    * POSTCONDITIONS: This iterator's list and Index no longer exist after the
    * move.
    *
    * ASSERTIONS: If the parameter Index does not exist, an assertion will be
    * triggered.
    *
    * @param [in] index The hier::Index to move this iterator's list elements
    * into.
    *
    */
   void
   move(
      const hier::Index& toIndex);

private:
   friend class SparseData<BOX_GEOMETRY>;

   typedef typename SparseData<BOX_GEOMETRY>::AttributeList AttributeList;
   /**********************************************************************
   * Data members
   **********************************************************************/
   SparseData<BOX_GEOMETRY>* d_data;
   typename SparseData<BOX_GEOMETRY>::IndexMap::iterator d_iterator;

   /**********************************************************************
   * private  methods for internal use only since they expose
   * implementation.
   **********************************************************************/
   SparseDataIterator(
      SparseData<BOX_GEOMETRY>& sparse_data,
      typename SparseData<BOX_GEOMETRY>::IndexMap::iterator iterator);

   void
   _insert(
      const typename SparseData<BOX_GEOMETRY>::Attributes& attributes);

   /*!
    * @brief prints the contents of this Iterator (Index + all attributes).
    *
    * @param [out] out The output stream
    */
   void
   printIterator(
      std::ostream& out) const;
};

/*!
 * @brief Iterator for attributes in an attribute list associated with an
 * Index.
 *
 * Each Index in a SparseData element has zero or more attributes collections
 * associated with it.  Each iterator element consists of a object representing
 * this collection.  The Iterator acts as both an iterator and an individual
 * element (located at the index associated with it).
 */
template<typename BOX_GEOMETRY>
class SparseDataAttributeIterator
{
public:
   template<typename T>
   friend
   std::ostream&
   operator << (
      std::ostream& out,
      SparseDataAttributeIterator<T>& attr_iterator);

   /*!
    * @brief Constructor
    *
    * @param [in] sparse_data
    * @param [in] index
    */
   SparseDataAttributeIterator(
      const SparseData<BOX_GEOMETRY>& sparse_data,
      const hier::Index& index);

   /*!
    * @brief Copy constructor
    *
    * @param [in] other
    */
   SparseDataAttributeIterator(
      const SparseDataAttributeIterator& other);

   /*!
    * @brief Iterator equality.
    *
    * @param [in] rhs
    */
   bool
   operator == (
      const SparseDataAttributeIterator& rhs) const;

   /*!
    * @brief Iterator inequality.
    *
    * @param [in] rhs
    */
   bool
   operator != (
      const SparseDataAttributeIterator& rhs) const;

   /*!
    * @brief pre-increment operator
    */
   SparseDataAttributeIterator&
   operator ++ ();

   /*!
    * @brief post-increment operator
    */
   SparseDataAttributeIterator
   operator ++ (
      int);

   /*!
    * @brief double attribute element access operator
    */
   double&
   operator [] (
      const DoubleAttributeId& id);

   /*!
    * @brief const double attribute element access operator
    */
   const double&
   operator [] (
      const DoubleAttributeId& id) const;

   /*!
    * @brief integer attribute element access operator
    */
   int&
   operator [] (
      const IntegerAttributeId& id);

   /*!
    * @brief const integer attribute element access operator
    */
   const int&
   operator [] (
      const IntegerAttributeId& id) const;

private:
   friend class SparseData<BOX_GEOMETRY>;
   typedef typename SparseData<BOX_GEOMETRY>::Attributes Attributes;
   typedef typename SparseData<BOX_GEOMETRY>::AttributeList AttributeList;

   SparseDataAttributeIterator(
      const AttributeList& attributes,
      const typename AttributeList::iterator& iterator);

   AttributeList d_list;
   typename AttributeList::iterator d_list_iterator;

   /*!
    * @brief print the attributes.  Called from the ostream<< operator.
    */
   void
   printAttribute(
      std::ostream& out) const;
};

} // namespace pdat
} // namespace SAMRAI

#include "SAMRAI/pdat/SparseData.cpp"

#endif // included_pdat_SparseData
