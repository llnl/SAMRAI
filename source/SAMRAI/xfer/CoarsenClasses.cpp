/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Simple structure for managing coarsening data in equivalence classes.
 *
 ************************************************************************/
#include "SAMRAI/xfer/CoarsenClasses.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/VariableDatabase.h"

#include <typeinfo>

namespace SAMRAI {
namespace xfer {

int CoarsenClasses::s_default_coarsen_item_array_size = 20;

/*
 *************************************************************************
 *
 * Constructor sets boolean for filling coarse data and creates new
 * array of equivalence classes.
 *
 *************************************************************************
 */

CoarsenClasses::CoarsenClasses():
   d_coarsen_classes_data_items(),
   d_num_coarsen_items(0)
{
}

/*
 *************************************************************************
 *
 * The destructor implicitly deletes the item storage associated with
 * the equivalence classes (and also the coarsen algorithm).
 *
 *************************************************************************
 */

CoarsenClasses::~CoarsenClasses()
{
}

/*
 *************************************************************************
 *
 * Insert a data item into the proper equivalence class.
 *
 *************************************************************************
 */

void
CoarsenClasses::insertEquivalenceClassItem(
   CoarsenClasses::Data& data,
   const std::shared_ptr<hier::PatchDescriptor>& descriptor)
{

   if (!itemIsValid(data, descriptor)) {
      tbox::perr << "Invalid coarsen class data passed to "
                 << "CoarsenClasses::insertEquivalenceClassItem\n";
      printCoarsenItem(tbox::perr, data);
      TBOX_ERROR("Check entries..." << std::endl);
   } else {

      int eq_index = getEquivalenceClassIndex(data, descriptor);

      if (eq_index < 0) {
         eq_index = static_cast<int>(d_equivalence_class_indices.size());
         d_equivalence_class_indices.resize(eq_index + 1);
      }

      data.d_class_index = eq_index;

      if (d_num_coarsen_items >=
          static_cast<int>(d_coarsen_classes_data_items.size())) {
         d_coarsen_classes_data_items.resize(
            d_num_coarsen_items + s_default_coarsen_item_array_size,
            Data(data.d_gcw_to_coarsen.getDim()));
      }

      d_coarsen_classes_data_items[d_num_coarsen_items] = data;

      d_equivalence_class_indices[eq_index].push_back(d_num_coarsen_items);

      ++d_num_coarsen_items;
   }

}

/*
 *************************************************************************
 *
 * Check for valid patch data ids, patch data types, and that source and
 * destination data entries have sufficient ghost cells to satisfy the
 * coarsen operator and necessary copy operations.  If so, return true;
 * else return false.  A descriptive error message is sent to TBOX_ERROR
 * when a problem appears.  If a null patch descriptor argument is
 * passed, the descriptor associated with the variable database
 * Singleton object will be used.
 *
 *************************************************************************
 */

bool
CoarsenClasses::itemIsValid(
   const CoarsenClasses::Data& data_item,
   const std::shared_ptr<hier::PatchDescriptor>& descriptor) const
{

   bool item_good = true;

   std::shared_ptr<hier::PatchDescriptor> pd(descriptor);
   if (!pd) {
      pd = hier::VariableDatabase::getDatabase()->getPatchDescriptor();
   }
   const tbox::Dimension& dim = pd->getPatchDataFactory(data_item.d_dst)->getDim();

   const int dst_id = data_item.d_dst;
   const int src_id = data_item.d_src;

   if (dst_id < 0) {
      item_good = false;
      TBOX_ERROR("Bad data given to CoarsenClasses...\n"
         << "`Destination' patch data id invalid (< 0!)" << std::endl);
   }
   if (item_good && (src_id < 0)) {
      item_good = false;
      TBOX_ERROR("Bad data given to CoarsenClasses...\n"
         << "`Source' patch data id invalid (< 0!)" << std::endl);
   }

   std::shared_ptr<hier::PatchDataFactory> dfact(
      pd->getPatchDataFactory(dst_id));
   std::shared_ptr<hier::PatchDataFactory> sfact(
      pd->getPatchDataFactory(src_id));

   if (item_good && !(sfact->validCopyTo(dfact))) {
      item_good = false;
      TBOX_ERROR("Bad data given to CoarsenClasses...\n"
         << "It is not a valid operation to copy from `Source' patch data \n"
         << pd->mapIndexToName(src_id) << " to `Destination' patch data "
         << pd->mapIndexToName(dst_id) << std::endl);
   }

   std::shared_ptr<hier::CoarsenOperator> coarsop(data_item.d_opcoarsen);
   if (item_good && coarsop) {
      if (coarsop->getStencilWidth(dim) > sfact->getGhostCellWidth()) {
         item_good = false;
         TBOX_ERROR("Bad data given to CoarsenClasses...\n"
            << "Coarsen operator " << coarsop->getOperatorName()
            << "\nhas larger stencil width than ghost cell width"
            << "of `Source' patch data" << pd->mapIndexToName(src_id)
            << "\noperator stencil width = " << coarsop->getStencilWidth(dim)
            << "\n`Source'  ghost width = "
            << sfact->getGhostCellWidth()
            << std::endl);
      }
   }

   return item_good;

}

/*
 *************************************************************************
 *
 * Compare the equivalence classes in this coarsen classes object against
 * those in the argument coarsen classes object.  Return true if both
 * object contain the same number of classes and the classes with the
 * class number match.  Two equivalence classes match if their
 * representatives are equivalent as defined by the method
 * itemsAreEquivalent().
 *
 * If a null patch descriptor argument is passed, the
 * descriptor associated with the variable database Singleton object
 * will be used.
 *
 *************************************************************************
 */

bool
CoarsenClasses::classesMatch(
   const std::shared_ptr<CoarsenClasses>& test_classes,
   const std::shared_ptr<hier::PatchDescriptor>& descriptor) const
{
   NULL_USE(descriptor);

   bool items_match = true;

   if (getNumberOfEquivalenceClasses() !=
       test_classes->getNumberOfEquivalenceClasses()) {

      items_match = false;

   } else {

      int eq_index = 0;
      while (items_match && eq_index < getNumberOfEquivalenceClasses()) {

         if (d_equivalence_class_indices[eq_index].size() !=
             test_classes->d_equivalence_class_indices[eq_index].size()) {

            items_match = false;

         } else {

            const CoarsenClasses::Data& my_item =
               getClassRepresentative(eq_index);
            const CoarsenClasses::Data& test_item =
               test_classes->getClassRepresentative(eq_index);

            items_match = itemsAreEquivalent(my_item, test_item);

         } // if number of items in equivalence class match

         ++eq_index;

      } // while equivalence classes match

   } // else number of equivalence classes do not match

   return items_match;

}

/*
 *************************************************************************
 *
 * Return true if data items are equivalent; false otherwise.
 * This routine defines coarsen item equivalence.
 *
 *************************************************************************
 */

bool
CoarsenClasses::itemsAreEquivalent(
   const CoarsenClasses::Data& data1,
   const CoarsenClasses::Data& data2,
   const std::shared_ptr<hier::PatchDescriptor>& descriptor) const
{
   bool equivalent = true;

   std::shared_ptr<hier::PatchDescriptor> pd(descriptor);
   if (!pd) {
      pd = hier::VariableDatabase::getDatabase()->getPatchDescriptor();
   }
   const tbox::Dimension& dim = pd->getPatchDataFactory(data1.d_dst)->getDim();

   equivalent = patchDataMatch(data1.d_dst, data2.d_dst, pd);

   equivalent &= patchDataMatch(data1.d_src, data2.d_src, pd);

   equivalent &= (data1.d_fine_bdry_reps_var == data2.d_fine_bdry_reps_var);

   equivalent &= (data1.d_gcw_to_coarsen == data2.d_gcw_to_coarsen);

   equivalent &= (!data1.d_opcoarsen == !data2.d_opcoarsen);
   if (equivalent && data1.d_opcoarsen) {
      equivalent &= (data1.d_opcoarsen->getStencilWidth(dim) ==
                     data2.d_opcoarsen->getStencilWidth(dim));
   }

   equivalent &= (!data1.d_var_fill_pattern ==
                  !data2.d_var_fill_pattern);
   if (equivalent && data1.d_var_fill_pattern) {
      auto& d1 = *(data1.d_var_fill_pattern);
      auto& d2 = *(data2.d_var_fill_pattern);
      equivalent &= (typeid(d1) == typeid(d2));
   }

   return equivalent;
}

/*
 *************************************************************************
 *
 * Print the data in the coarsen item lists to the specified stream.
 *
 *************************************************************************
 */

void
CoarsenClasses::printClassData(
   std::ostream& stream) const
{
   stream << "CoarsenClasses::printClassData()\n";
   stream << "--------------------------------------\n";
   for (int i = 0; i < static_cast<int>(d_equivalence_class_indices.size()); ++i) {
      stream << "EQUIVALENCE CLASS # " << i << std::endl;
      int j = 0;
      const std::list<int>& indices = d_equivalence_class_indices[i];
      for (std::list<int>::const_iterator li(indices.begin());
           li != indices.end(); ++li) {

         stream << "Item # " << j << std::endl;
         stream << "-----------------------------\n";

         printCoarsenItem(stream, d_coarsen_classes_data_items[*li]);

         ++j;
      }
      stream << std::endl;
   }

}

void
CoarsenClasses::printCoarsenItem(
   std::ostream& stream,
   const CoarsenClasses::Data& data) const
{
   stream << "\n";
   stream << "desination component:   "
          << data.d_dst << std::endl;
   stream << "source component:       "
          << data.d_src << std::endl;
   stream << "fine boundary represents variable:       "
          << data.d_fine_bdry_reps_var << std::endl;
   stream << "gcw to coarsen:       "
          << data.d_gcw_to_coarsen << std::endl;
   stream << "tag:       "
          << data.d_tag << std::endl;

   if (!data.d_opcoarsen) {
      stream << "NULL coarsening operator" << std::endl;
   } else {
      auto& d = *data.d_opcoarsen;
      stream << "coarsen operator name:          "
             << typeid(d).name()
             << std::endl;
      stream << "operator priority:      "
             << data.d_opcoarsen->getOperatorPriority()
             << std::endl;
      stream << "operator stencil width: "
             << data.d_opcoarsen->getStencilWidth(
         hier::VariableDatabase::getDatabase()->getPatchDescriptor()->getPatchDataDim(data.d_dst))
             << std::endl;
   }
   stream << std::endl;
}

/*
 *************************************************************************
 *
 * Private member function to determine whether two patch data items
 * match (are same type and have same ghost width).
 *
 *************************************************************************
 */

bool
CoarsenClasses::patchDataMatch(
   int item_id1,
   int item_id2,
   const std::shared_ptr<hier::PatchDescriptor>& pd) const
{

   bool items_match = ((item_id1 >= 0) && (item_id2 >= 0));

   if (items_match) {

      std::shared_ptr<hier::PatchDataFactory> pdf1(
         pd->getPatchDataFactory(item_id1));
      std::shared_ptr<hier::PatchDataFactory> pdf2(
         pd->getPatchDataFactory(item_id2));

      auto& p1 = *pdf1;
      auto& p2 = *pdf2;
      items_match = (typeid(p1) == typeid(p2));

      if (items_match) {
         items_match = (pdf1->getGhostCellWidth() ==
                        pdf2->getGhostCellWidth());
      }

   }

   return items_match;

}

/*
 *************************************************************************
 *
 * Private member function to determine equivalence class for
 * given data item.  Return value of -1 indicates no match found; else
 * return value is index of match.
 *
 *************************************************************************
 */

int
CoarsenClasses::getEquivalenceClassIndex(
   const CoarsenClasses::Data& data,
   const std::shared_ptr<hier::PatchDescriptor>& descriptor) const
{
   NULL_USE(descriptor);

   int eq_index = -1;

   bool class_found = false;
   int check_index = 0;
   while (!class_found && check_index < getNumberOfEquivalenceClasses()) {

      const CoarsenClasses::Data& class_rep =
         getClassRepresentative(check_index);

      class_found = itemsAreEquivalent(data, class_rep);

      if (class_found) {
         eq_index = check_index;
      }

      ++check_index;
   }

   return eq_index;

}

/*
 *************************************************************************
 * Constructor
 *************************************************************************
 */
CoarsenClasses::Data::Data(
   tbox::Dimension dim):
   d_gcw_to_coarsen(dim)
{
}

}
}
