/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Simple structure for managing refinement data in equivalence classes.
 *
 ************************************************************************/
#include <typeinfo>

#include "SAMRAI/xfer/RefineClasses.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace xfer {

int RefineClasses::s_default_refine_item_array_size = 20;

/*
 *************************************************************************
 *
 * Default constructor.
 *
 *************************************************************************
 */

RefineClasses::RefineClasses()
{
   d_refine_classes_data_items.reserve(s_default_refine_item_array_size);
}

/*
 *************************************************************************
 *
 * The destructor implicitly deletes the item storage associated with
 * the equivalence classes (and also the refine algorithm).
 *
 *************************************************************************
 */

RefineClasses::~RefineClasses()
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
RefineClasses::insertEquivalenceClassItem(
   RefineClasses::Data& data,
   const std::shared_ptr<hier::PatchDescriptor>& descriptor)
{

   if (!itemIsValid(data, descriptor)) {
      tbox::perr << "Invalid refine class data passed to "
                 << "RefineClasses::insertEquivalenceClassItem\n";
      printRefineItem(tbox::perr, data);
      TBOX_ERROR("Check entries..." << std::endl);
   } else {

      int eq_index = getEquivalenceClassIndex(data, descriptor);

      if (eq_index < 0) {
         eq_index = static_cast<int>(d_equivalence_class_indices.size());
         d_equivalence_class_indices.resize(eq_index + 1);
      }

      data.d_class_index = eq_index;

      d_equivalence_class_indices[eq_index].push_back(static_cast<int>(d_refine_classes_data_items.
                                                                       size()));
      d_refine_classes_data_items.push_back(data);

   }

}

/*
 *************************************************************************
 *
 * Check for valid patch data ids, patch data types, and that scratch
 * data entry has at least as many ghost cells as destination data entry
 * and stencil width of operator.  If so, return true; else return false.
 * A descriptive error message is sent to TBOX_ERROR when a problem
 * appears.  If a null patch descriptor argument is passed, the
 * descriptor associated with the variable database Singleton object
 * will be used.
 *
 *************************************************************************
 */

bool
RefineClasses::itemIsValid(
   const RefineClasses::Data& data_item,
   const std::shared_ptr<hier::PatchDescriptor>& descriptor) const
{

   bool item_good = true;

   std::shared_ptr<hier::PatchDescriptor> pd(descriptor);
   if (!pd) {
      pd = hier::VariableDatabase::getDatabase()->getPatchDescriptor();
   }

   const int dst_id = data_item.d_dst;
   const int src_id = data_item.d_src;
   const int scratch_id = data_item.d_scratch;

   if (dst_id < 0) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses...\n"
         << "`Destination' patch data id invalid (< 0!)" << std::endl);
   }
   if (item_good && (src_id < 0)) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses...\n"
         << "`Source' patch data id invalid (< 0!)" << std::endl);
   }
   if (item_good && (scratch_id < 0)) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses...\n"
         << "`Scratch' patch data id invalid (< 0!)" << std::endl);
   }

   const std::vector<int>& work_ids = data_item.d_work;
   if (item_good && !work_ids.empty()) {
      for (std::vector<int>::const_iterator itr = work_ids.begin();
           itr != work_ids.end(); ++itr) {
         if (item_good && (*itr < 0)) {
            item_good = false;
            TBOX_ERROR("Bad data given to RefineClasses...\n"
               << "`Work' patch data id invalid (< 0!)" << std::endl);
         }
      }
   }

   std::shared_ptr<hier::PatchDataFactory> dst_fact(
      pd->getPatchDataFactory(dst_id));
   std::shared_ptr<hier::PatchDataFactory> src_fact(
      pd->getPatchDataFactory(src_id));
   std::shared_ptr<hier::PatchDataFactory> scratch_fact(
      pd->getPatchDataFactory(scratch_id));

   const tbox::Dimension& dim = dst_fact->getDim();

   if (item_good && !(src_fact->validCopyTo(scratch_fact))) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses...\n"
         << "It is not a valid operation to copy from `Source' patch data \n"
         << pd->mapIndexToName(src_id) << " to `Scratch' patch data "
         << pd->mapIndexToName(scratch_id) << std::endl);
   }

   if (item_good && !(scratch_fact->validCopyTo(dst_fact))) {
      item_good = false;
      pd->mapIndexToName(scratch_id);
      pd->mapIndexToName(dst_id);
      TBOX_ERROR("Bad data given to RefineClasses...\n"
         << "It is not a valid operation to copy from `Scratch' patch data \n"
         << pd->mapIndexToName(scratch_id) << " to `Destination' patch data "
         << pd->mapIndexToName(dst_id) << std::endl);
   }

   const hier::IntVector& scratch_gcw = scratch_fact->getGhostCellWidth();

   if (item_good && (dst_fact->getGhostCellWidth() > scratch_gcw)) {
      item_good = false;
      TBOX_ERROR("Bad data given to RefineClasses...\n"
         << "`Destination' patch data " << pd->mapIndexToName(dst_id)
         << " has a larger ghost cell width than \n"
         << "`Scratch' patch data " << pd->mapIndexToName(scratch_id)
         << "\n`Destination' ghost width = "
         << dst_fact->getGhostCellWidth()
         << "\n`Scratch' ghost width = " << scratch_gcw << std::endl);
   }

   std::shared_ptr<hier::RefineOperator> refop(data_item.d_oprefine);
   if (item_good && refop) {
      if (refop->getStencilWidth(dim) > scratch_gcw) {
         item_good = false;
         TBOX_ERROR("Bad data given to RefineClasses...\n"
            << "Refine operator " << refop->getOperatorName()
            << "\nhas larger stencil width than ghost cell width"
            << "of `Scratch' patch data" << pd->mapIndexToName(scratch_id)
            << "\noperator stencil width = " << refop->getStencilWidth(scratch_gcw.getDim())
            << "\n`Scratch'  ghost width = " << scratch_gcw << std::endl);
      }
   }

   std::shared_ptr<VariableFillPattern> fill_pattern(
      data_item.d_var_fill_pattern);
   if (item_good && fill_pattern) {
      if (fill_pattern->getPatternName() != "BOX_GEOMETRY_FILL_PATTERN") {
         if (fill_pattern->getStencilWidth() > scratch_gcw) {
            item_good = false;
            TBOX_ERROR("Bad data given to RefineClasses...\n"
               << "VariableFillPattern " << fill_pattern->getPatternName()
               << "\nhas larger stencil width than ghost cell width"
               << "of `Scratch' patch data" << pd->mapIndexToName(
                  scratch_id)
               << "\noperator stencil width = "
               << fill_pattern->getStencilWidth()
               << "\n`Scratch'  ghost width = " << scratch_gcw << std::endl);
         }
      }
   }

   if (item_good && data_item.d_time_interpolate) {
      const int src_told_id = data_item.d_src_told;
      const int src_tnew_id = data_item.d_src_tnew;

      if (src_told_id < 0) {
         item_good = false;
         TBOX_ERROR("Bad data given to RefineClasses...\n"
            << "`Source old' patch data id invalid (< 0!),\n"
            << "yet a request has made to time interpolate" << std::endl);
      }
      if (item_good && src_tnew_id < 0) {
         item_good = false;
         TBOX_ERROR("Bad data given to RefineClasses...\n"
            << "`Source new' patch data id invalid (< 0!),\n"
            << "yet a request has made to time interpolate with them"
            << std::endl);
      }

      std::shared_ptr<hier::PatchDataFactory> src_told_fact(
         pd->getPatchDataFactory(src_told_id));
      std::shared_ptr<hier::PatchDataFactory> src_tnew_fact(
         pd->getPatchDataFactory(src_tnew_id));

      auto& told_fact = *src_told_fact;
      auto& fact      = *src_fact;
      if (item_good && typeid(told_fact) != typeid(fact)) {
         item_good = false;
         TBOX_ERROR("Bad data given to RefineClasses...\n"
            << "`Source' patch data " << pd->mapIndexToName(src_id)
            << " and `Source old' patch data "
            << pd->mapIndexToName(src_told_id)
            << " have different patch data types, yet a request has"
            << "\n been made to time interpolate with them" << std::endl);
      }

      auto& tnew_fact = *src_tnew_fact;
      if (item_good && typeid(tnew_fact) != typeid(fact)) {
         item_good = false;
         TBOX_ERROR("Bad data given to RefineClasses...\n"
            << "`Source' patch data " << pd->mapIndexToName(src_id)
            << " and `Source new' patch data "
            << pd->mapIndexToName(src_tnew_id)
            << " have different patch data types, yet a request has"
            << "\n been made to time interpolate with them" << std::endl);
      }

   }

   return item_good;

}

/*
 *************************************************************************
 *
 * Compare the equivalence classes in this refine classes object against
 * those in the argument refine classes object.  Return true if both
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
RefineClasses::classesMatch(
   const std::shared_ptr<RefineClasses>& test_classes,
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

            const RefineClasses::Data& my_item =
               getClassRepresentative(eq_index);
            const RefineClasses::Data& test_item =
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
 * This routine defines refine item equivalence.
 *
 *************************************************************************
 */

bool
RefineClasses::itemsAreEquivalent(
   const RefineClasses::Data& data1,
   const RefineClasses::Data& data2,
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

   equivalent &= patchDataMatch(data1.d_scratch, data2.d_scratch, pd);

   equivalent &= (data1.d_work.size() == data2.d_work.size());

   if (equivalent && !data1.d_work.empty()) {
      std::vector<int>::const_iterator itr1 = data1.d_work.begin();
      std::vector<int>::const_iterator itr2 = data2.d_work.begin();
      for ( ; itr1 != data1.d_work.end() && itr2 != data1.d_work.end();
           ++itr1, ++itr2 ) {
         equivalent &= patchDataMatch(*itr1, *itr2, pd);
      }
   }

   equivalent &= (data1.d_time_interpolate == data2.d_time_interpolate);
   if (equivalent && data1.d_time_interpolate) {
      equivalent &= patchDataMatch(data1.d_src_told, data2.d_src_told, pd);
      equivalent &= patchDataMatch(data1.d_src_tnew, data2.d_src_tnew, pd);
   }

   equivalent &= (data1.d_fine_bdry_reps_var == data2.d_fine_bdry_reps_var);

   equivalent &= (!data1.d_oprefine == !data2.d_oprefine);
   if (equivalent && data1.d_oprefine) {
      equivalent &= (data1.d_oprefine->getStencilWidth(dim) ==
                     data2.d_oprefine->getStencilWidth(dim));
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
 * Print the data in the refine item lists to the specified stream.
 *
 *************************************************************************
 */

void
RefineClasses::printClassData(
   std::ostream& stream) const
{
   stream << "RefineClasses::printClassData()\n";
   stream << "--------------------------------------\n";
   for (int i = 0; i < static_cast<int>(d_equivalence_class_indices.size()); ++i) {
      stream << "EQUIVALENCE CLASS # " << i << std::endl;
      int j = 0;
      const std::list<int>& indices = d_equivalence_class_indices[i];
      for (std::list<int>::const_iterator li(indices.begin());
           li != indices.end(); ++li) {

         stream << "Item # " << j << std::endl;
         stream << "-----------------------------\n";

         printRefineItem(stream, d_refine_classes_data_items[*li]);

         ++j;
      }
      stream << std::endl;
   }

}

void
RefineClasses::printRefineItem(
   std::ostream& stream,
   const RefineClasses::Data& data) const
{
   stream << "\n";
   stream << "desination component:   "
          << data.d_dst << std::endl;
   stream << "source component:       "
          << data.d_src << std::endl;
   stream << "scratch component:      "
          << data.d_scratch << std::endl;
   stream << "fine boundary represents variable:      "
          << data.d_fine_bdry_reps_var << std::endl;
   stream << "tag:      "
          << data.d_tag << std::endl;

   if (!data.d_work.empty()) {
      const std::vector<int>& work_ids = data.d_work;
      for (std::vector<int>::const_iterator itr = work_ids.begin();
           itr != work_ids.end(); ++itr) {
         stream << "work component:      "
                << *itr << std::endl;
      }
   }

   if (!data.d_oprefine) {
      stream << "NULL refining operator" << std::endl;
   } else {
      auto& oprefine = *data.d_oprefine;
      stream << "refine operator name:          "
             << typeid(oprefine).name()
             << std::endl;
      stream << "operator priority:      "
             << data.d_oprefine->getOperatorPriority()
             << std::endl;
      stream << "operator stencil width: "
             << data.d_oprefine->getStencilWidth(
         hier::VariableDatabase::getDatabase()->getPatchDescriptor()->getPatchDataDim(data.d_dst))
             << std::endl;
   }
   if (!data.d_time_interpolate) {
      stream << "time interpolate is false" << std::endl;
   } else {
      auto& optime = *data.d_optime;
      stream << "old source component:   "
             << data.d_src_told << std::endl;
      stream << "new source component:   "
             << data.d_src_tnew << std::endl;
      stream << "time interpolation operator name:          "
             << typeid(optime).name()
             << std::endl;
   }
   if (!data.d_var_fill_pattern) {
      stream << "var fill pattern is null" << std::endl;
   } else {
      auto& d = *data.d_var_fill_pattern;
      stream << "var fill pattern name:          "
             << typeid(d).name()
             << std::endl;
   }
   stream << std::endl;
}

/*
 *************************************************************************
 *
 * Private member function to determine whether two patch data items
 * match (are same type and have same ghost width.)
 *
 *************************************************************************
 */

bool
RefineClasses::patchDataMatch(
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
 * Private member function to determine equivalence class index for
 * given data item.  Return value of -1 indicates no match found; else
 * return value is index of match.
 *
 *************************************************************************
 */

int
RefineClasses::getEquivalenceClassIndex(
   const RefineClasses::Data& data,
   const std::shared_ptr<hier::PatchDescriptor>& descriptor) const
{
   NULL_USE(descriptor);

   int eq_index = -1;

   bool class_found = false;
   int check_index = 0;
   while (!class_found && check_index < getNumberOfEquivalenceClasses()) {

      const RefineClasses::Data& class_rep =
         getClassRepresentative(check_index);

      class_found = itemsAreEquivalent(data, class_rep);

      if (class_found) {
         eq_index = check_index;
      }

      ++check_index;
   }

   return eq_index;

}

}
}
