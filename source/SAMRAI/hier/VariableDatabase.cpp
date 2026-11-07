/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Manager class for variables used in a SAMRAI application.
 *
 ************************************************************************/
#include "SAMRAI/hier/VariableDatabase.h"

#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

 /*
  * Static data members used to control access to and destruction of
  * singleton variable database instance.
  */
class VariableDatabaseInstance{
public:
  static VariableDatabaseInstance &instance(){
    static VariableDatabaseInstance self;
    return self;
  }
  VariableDatabase * get(){
    return s_variable_database_instance;
  }
  void set(VariableDatabase * ptr){
    s_variable_database_instance = ptr;
  }
private:
  VariableDatabase * s_variable_database_instance = nullptr;
};

const int VariableDatabase::s_context_array_alloc_size(10);
const int VariableDatabase::s_variable_array_alloc_size(100);
const int VariableDatabase::s_descriptor_array_alloc_size(200);

tbox::StartupShutdownManager::Handler
VariableDatabase::s_shutdown_handler(
   0,
   0,
   VariableDatabase::shutdownCallback,
   0,
   tbox::StartupShutdownManager::priorityVariableDatabase);

/*
 *************************************************************************
 *
 * Static database member functions.
 *
 *************************************************************************
 */

VariableDatabase *
VariableDatabase::getDatabase()
{
   if (!VariableDatabaseInstance::instance().get()) {
      VariableDatabaseInstance::instance().set(new VariableDatabase());
   }
   return VariableDatabaseInstance::instance().get();
}

void
VariableDatabase::shutdownCallback()
{
   VariableDatabase * s_variable_database_instance = VariableDatabaseInstance::instance().get();
   if (s_variable_database_instance) {
      delete s_variable_database_instance;
   }
   VariableDatabaseInstance::instance().set(nullptr);
}

/*
 *************************************************************************
 *
 * Protected VariableDatabase constructor, destructor, and function to
 * register Singleton subclass instance for inheritance.
 *
 *************************************************************************
 */

VariableDatabase::VariableDatabase():
   d_patch_descriptor(std::make_shared<PatchDescriptor>())
{
   d_max_variable_id = idUndefined();
   d_max_context_id = idUndefined();
   d_max_descriptor_id = idUndefined();
   d_num_registered_patch_data_ids = 0;

   d_internal_SAMRAI_context = getContext("Internal_SAMRAI_Variable");
}

void
VariableDatabase::registerSingletonSubclassInstance(
   VariableDatabase* subclass_instance)
{
   if (!VariableDatabaseInstance::instance().get()) {
      VariableDatabaseInstance::instance().set(subclass_instance);
   } else {
      TBOX_ERROR("hier::VariableDatabase internal error...\n"
         << "Attemptng to set Singleton instance to subclass instance,"
         << "\n but Singleton instance already set." << std::endl);
   }
}

/*
 *************************************************************************
 *
 * Accessory functions to retrieve data members.
 *
 *************************************************************************
 */

std::shared_ptr<PatchDescriptor>
VariableDatabase::getPatchDescriptor() const
{
   return d_patch_descriptor;
}

int
VariableDatabase::getNumberOfRegisteredPatchDataIndices() const
{
   return d_num_registered_patch_data_ids;
}

int
VariableDatabase::getNumberOfRegisteredVariableContexts() const
{
   // currently, we do not allow removal of variable contexts
   // so this suffices
   return d_max_context_id + 1;
}

/*
 *************************************************************************
 *
 * Return the context in the database with the given name, or add a
 * context to the database with that name if no such context exists.
 *
 *************************************************************************
 */

std::shared_ptr<VariableContext>
VariableDatabase::getContext(
   const std::string& name)
{
   std::shared_ptr<VariableContext> context;

   if (!name.empty()) {

      int ctxt_id = getContextId_Private(name);

      if (ctxt_id == idUndefined()) {
         context.reset(new VariableContext(name));
         addContext_Private(context);
      } else {
         context = d_contexts[ctxt_id];
      }

   }

   return context;
}

/*
 *************************************************************************
 *
 * Return true if context with given name exists in database;
 * otherwise return false.
 *
 *************************************************************************
 */

bool
VariableDatabase::checkContextExists(
   const std::string& name) const
{
   int ctxt_id = getContextId_Private(name);

   return ctxt_id != idUndefined();
}

/*
 *************************************************************************
 *
 * Add user-defined variable to database if it doesn't already exist in
 * the database.
 *
 *************************************************************************
 */

void
VariableDatabase::addVariable(
   const std::shared_ptr<Variable>& variable)
{
   TBOX_ASSERT(variable);

   const bool user_variable = true;
   bool variable_added = addVariable_Private(variable, user_variable);

   if (!variable_added) {
      TBOX_ERROR("hier::VariableDatabase::addVariable() error...\n"
         << "Attempt to add variable with duplicate name " << variable->getName()
         << " to database is not allowed.\n"
         << "Another variable with this name already exists in database."
         << std::endl);
   }

}

/*
 *************************************************************************
 *
 * Return variable in database with given name.  If no such variable
 * resides in database, return a null pointer.
 *
 *************************************************************************
 */

std::shared_ptr<Variable>
VariableDatabase::getVariable(
   const std::string& name) const
{
   std::shared_ptr<Variable> variable;

   int var_id = getVariableId(name);

   if (var_id != idUndefined()) {
      variable = d_variables[var_id];
   }

   return variable;
}

/*
 *************************************************************************
 *
 * Return true if variable with given name exists in database.
 * Otherwise, return false.
 *
 *************************************************************************
 */

bool
VariableDatabase::checkVariableExists(
   const std::string& name) const
{
   int var_id = getVariableId(name);

   return var_id != idUndefined();
}

/*
 *************************************************************************
 *
 * Create new patch data index index by cloning factory for variable
 * at the old index and return index of new factory.   Note that the
 * function checkVariablePatchDataIndex() checks type of variable
 * against given patch data index.   If these types match, then add
 * variable and new patch data index to database.  If the types do not
 * match, the program will abort with an error message in the private
 * routine checkVariablePatchDataIndex().
 *
 *************************************************************************
 */

int
VariableDatabase::registerClonedPatchDataIndex(
   const std::shared_ptr<Variable>& variable,
   int old_id)
{
   TBOX_ASSERT(variable);

   int new_id = idUndefined();

   if (checkVariablePatchDataIndex(variable, old_id)) {

      std::string old_name = d_patch_descriptor->mapIndexToName(old_id);
      std::string old_id_string(tbox::Utilities::intToString(old_id, 4));

      std::string new_name;
      if (old_name.find("-clone_of_id=") == std::string::npos) {
         new_name = old_name + "-clone_of_id=" + old_id_string;
      } else {
         std::string::size_type last_dash = old_name.rfind("=");
         new_name = old_name.substr(0, last_dash + 1) + old_id_string;
      }

      new_id = d_patch_descriptor->definePatchDataComponent(
            new_name,
            d_patch_descriptor->getPatchDataFactory(old_id)->
            cloneFactory(d_patch_descriptor->getPatchDataFactory(old_id)->
               getGhostCellWidth()));

      const bool user_variable = true;
      addVariablePatchDataIndexPairToDatabase_Private(variable,
         new_id,
         user_variable);

   } else {

      auto& pdf = *(d_patch_descriptor->getPatchDataFactory(old_id));

      TBOX_ERROR("hier::VariableDatabase::registerClonedPatchDataIndex()"
         << "  error...\n"
         << "Variable with name " << variable->getName()
         << "\n does not match type at descriptor index = " << old_id
         << "\n That type is " << typeid(pdf).name()
         << std::endl);
   }

   return new_id;
}

/*
 *************************************************************************
 *
 * Add patch data index and variable pair to the database.  Note
 * that the function checkVariablePatchDataIndex() checks type of
 * variable against given patch data index.  If the types do not match,
 * the program will abort with an error message in the private routine
 * checkVariablePatchDataIndex().   If the input index is undefined,
 * we clone the default variable factory and add this new index to the
 * database.  In any case, the index of the index-variable pair that
 * is added to the database is returned.
 *
 *************************************************************************
 */

int
VariableDatabase::registerPatchDataIndex(
   const std::shared_ptr<Variable>& variable,
   int data_id)
{
   TBOX_ASSERT(variable);

   int new_id = data_id;

   if (new_id == idUndefined()) {

      new_id = d_patch_descriptor->definePatchDataComponent(
            variable->getName(),
            variable->getPatchDataFactory()->cloneFactory(
               variable->getPatchDataFactory()->getGhostCellWidth()));

      const bool user_variable = true;
      addVariablePatchDataIndexPairToDatabase_Private(variable,
         new_id,
         user_variable);

   } else {

      if (checkVariablePatchDataIndex(variable, new_id)) {

         const bool user_variable = true;
         addVariablePatchDataIndexPairToDatabase_Private(variable,
            new_id,
            user_variable);

      } else {

         auto& pdf = *(d_patch_descriptor->getPatchDataFactory(data_id));
         TBOX_ERROR("hier::VariableDatabase::registerPatchDataIndex()"
            << "  error...\n"
            << "Variable with name " << variable->getName()
            << "\n does not match type at patch data index = " << new_id
            << "\n That type is " << typeid(pdf).name()
            << std::endl);

      }

   }

   return new_id;
}

/*
 *************************************************************************
 *
 * Remove the given patch data index from the database.  Also, clear
 * the index from the patch descriptor if the index is in the database.
 *
 *************************************************************************
 */

void
VariableDatabase::removePatchDataIndex(
   int data_id)
{

   if ((data_id >= 0) && (data_id <= d_max_descriptor_id)) {

      std::shared_ptr<Variable> variable(d_index2variable_map[data_id]);

      if (variable) {

         std::vector<int>& indx_array =
            d_variable_context2index_map[variable->getInstanceIdentifier()];
         int array_size = static_cast<int>(indx_array.size());
         for (int i = 0; i < array_size; ++i) {
            if (indx_array[i] == data_id) {
               indx_array[i] = idUndefined();
               break;
            }
         }

         d_patch_descriptor->removePatchDataComponent(data_id);

         if (d_index2variable_map[data_id]) {
            --d_num_registered_patch_data_ids;
         }

         d_index2variable_map[data_id].reset();
         if (data_id == d_max_descriptor_id) {
            for (int id = d_max_descriptor_id; id >= 0; --id) {
               if (!d_index2variable_map[id]) {
                  --d_max_descriptor_id;
               } else {
                  break;
               }
            }
         }

      }

   }

}

/*
 *************************************************************************
 *
 * Return true if the given variable is mapped to the given patch data
 * index.  Otherwise, return false.
 *
 *************************************************************************
 */

bool
VariableDatabase::checkVariablePatchDataIndex(
   const std::shared_ptr<Variable>& variable,
   int data_id) const
{
   TBOX_ASSERT(variable);
   TBOX_ASSERT(data_id >= 0 &&
      data_id < d_patch_descriptor->getMaxNumberRegisteredComponents());

   bool ret_value = false;

   std::shared_ptr<Variable> test_variable;

   if ((data_id >= 0) && (data_id <= d_max_descriptor_id)) {
      test_variable = d_index2variable_map[data_id];
   }

   if (test_variable) {

      ret_value = (variable.get() == test_variable.get());

   }

   return ret_value;
}

/*
 *************************************************************************
 *
 * Return true if the type of the variable matches the type of the
 * patch data at the given patch data index.  Otherwise, return false.
 *
 *************************************************************************
 */

bool
VariableDatabase::checkVariablePatchDataIndexType(
   const std::shared_ptr<Variable>& variable,
   int data_id) const
{
   TBOX_ASSERT(variable);
   TBOX_ASSERT(data_id >= 0 &&
      data_id < d_patch_descriptor->getMaxNumberRegisteredComponents());

   bool ret_value = false;

   if (d_patch_descriptor->getPatchDataFactory(data_id)) {

      std::shared_ptr<PatchDataFactory> dfact(
         d_patch_descriptor->getPatchDataFactory(data_id));

      auto& pdf = *(variable->getPatchDataFactory());
      auto& df  = *dfact;
      if (dfact && (typeid(pdf) == typeid(df))) {
         ret_value = true;
      }

   }

   return ret_value;
}

/*
 *************************************************************************
 *
 * Register variable-context pair with the database and return
 * patch daya index corresponding to this pair and given ghost width.
 *
 *************************************************************************
 */

int
VariableDatabase::registerVariableAndContext(
   const std::shared_ptr<Variable>& variable,
   const std::shared_ptr<VariableContext>& context,
   const IntVector& ghosts)
{
   TBOX_ASSERT(variable);
   TBOX_ASSERT(context);
   TBOX_ASSERT(ghosts.min() >= 0);

   bool user_variable = true;
   return registerVariableAndContext_Private(variable,
      context,
      ghosts,
      user_variable);

}

/*
 *************************************************************************
 *
 * Return patch data index that is mapped to given variable-context
 * pair.  If variable-context pair does not exist in database, return
 * an undefined patch data index of idUndefined().
 *
 *************************************************************************
 */

int
VariableDatabase::mapVariableAndContextToIndex(
   const std::shared_ptr<Variable>& variable,
   const std::shared_ptr<VariableContext>& context) const
{
   TBOX_ASSERT(variable);
   TBOX_ASSERT(context);

   int index = idUndefined();

   int var_id = variable->getInstanceIdentifier();
   int ctxt_id = context->getIndex();

   if ((var_id <= d_max_variable_id) &&
       (ctxt_id < static_cast<int>(d_variable_context2index_map[var_id].size()))) {

      index = d_variable_context2index_map[var_id][ctxt_id];

   }

   return index;

}

/*
 *************************************************************************
 *
 * Return true if given patch data index is mapped to some variable
 * in the database and set the variable pointer to that variable.
 * Otherwise, return false and set the variable pointer to null.
 *
 *************************************************************************
 */

bool
VariableDatabase::mapIndexToVariable(
   const int index,
   std::shared_ptr<Variable>& variable) const
{
   variable.reset();

   if ((index >= 0) && (index <= d_max_descriptor_id)) {
      variable = d_index2variable_map[index];
   }

   return variable.get();
}

/*
 *************************************************************************
 *
 * Return true if specified index is mapped to some variable-context
 * pair in the database and set the variable and context pointers
 * appropriately.  Otherwise, return false and set the pointers to null.
 *
 *************************************************************************
 */

bool
VariableDatabase::mapIndexToVariableAndContext(
   const int index,
   std::shared_ptr<Variable>& variable,
   std::shared_ptr<VariableContext>& context) const
{
   bool found = false;

   variable.reset();
   context.reset();

   if ((index >= 0) && (index <= d_max_descriptor_id)) {

      variable = d_index2variable_map[index];

      if (variable) {

         const std::vector<int>& var_indx_array =
            d_variable_context2index_map[variable->getInstanceIdentifier()];
         int arr_size = static_cast<int>(var_indx_array.size());
         for (int i = 0; i < arr_size; ++i) {
            if (var_indx_array[i] == index) {
               found = true;
               context = d_contexts[i];
               break;
            }
         }

      }

   }

   return found;

}

/*
 *************************************************************************
 *
 * Print all context, variable, and patch data index data
 * contained in database to given output stream.
 *
 *************************************************************************
 */

void
VariableDatabase::printClassData(
   std::ostream& os,
   bool print_only_user_defined_variables) const
{
   int i;
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << std::endl;
   os << "Printing hier::VariableDatabase information...";
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << std::endl;
   os << "Variable Contexts registered with database:";
   for (i = 0; i <= d_max_context_id; ++i) {
      os << "\nContext id = " << i;
      if (d_contexts[i]) {
         os << " : Context name = " << d_contexts[i]->getName();
      } else {
         os << " : NOT IN DATABASE";
      }
   }
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << std::endl << std::flush;
   os << "Variables registered with database:";
   for (i = 0; i <= d_max_variable_id; ++i) {
      os << "\nVariable instance = " << i;
      if (d_variables[i]) {
         os << "\n";
         if (!print_only_user_defined_variables ||
             (print_only_user_defined_variables &&
              d_is_user_variable[i])) {
            auto& v = *(d_variables[i]);
            os << "   Variable name = " << d_variables[i]->getName();
            os << "\n   Variable type = " << typeid(v).name();
         } else {
            os << "   internal SAMRAI variable";
         }
      } else {
         os << " : NOT IN DATABASE";
      }
   }
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << std::endl << std::flush;
   os << "Variable-Context pairs mapping to Patch Data Indices in database:";
   for (i = 0; i <= d_max_variable_id; ++i) {
      if (d_variables[i]) {
         if (!print_only_user_defined_variables ||
             (print_only_user_defined_variables &&
              d_is_user_variable[i])) {
            os << "\nVariable name = " << d_variables[i]->getName();
            int nctxts =
               static_cast<int>(d_variable_context2index_map[i].size());
            if (nctxts > 0) {
               for (int j = 0; j < nctxts; ++j) {
                  if (d_variable_context2index_map[i][j] != idUndefined()) {
                     os << "\n   context id = " << j << ", name = "
                        << d_contexts[j]->getName()
                        << " :  patch data id = "
                        << d_variable_context2index_map[i][j];
                  } else {
                     os << "\n   context id = " << j
                        << " UNDEFINED for this variable";
                  }
               }
            } else {
               os << "\n   --- No contexts defined ---";
            }
         }
      }
   }
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << std::endl << std::flush;
   os << "Mapping from Patch Data Indices to Variables:";
   for (i = 0; i <= d_max_descriptor_id; ++i) {
      os << "\nPatch data id = " << i << " -- ";
      if (!d_index2variable_map[i]) {
         os << "UNDEFINED in database";
      } else {
         int vid = d_index2variable_map[i]->getInstanceIdentifier();
         if (!print_only_user_defined_variables ||
             (print_only_user_defined_variables &&
              d_is_user_variable[vid])) {
            os << "data factory name = "
               << d_patch_descriptor->mapIndexToName(i);
         } else {
            os << "internal SAMRAI patch data";
         }
      }
   }
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << std::endl << std::flush;
   os
   << "Printing contents of patch descriptor for comparison to database..."
   << std::endl;
   d_patch_descriptor->printClassData(os);
   os << std::flush;
   os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      << std::endl << std::flush;
}

/*
 *************************************************************************
 *
 * Register internal SAMRAI variable with database using internal
 * SAMRAI variable context.
 *
 * If the variable is already registered with the database as a user-
 * defined variable, then an unrecoverable error results.  This avoids
 * potential naming conflicts with user variables.
 *
 *************************************************************************
 */

int
VariableDatabase::registerInternalSAMRAIVariable(
   const std::shared_ptr<Variable>& variable,
   const IntVector& ghosts)
{
   TBOX_ASSERT(variable);
   TBOX_ASSERT(ghosts.min() >= 0);

   int data_id = idUndefined();

   int var_id = variable->getInstanceIdentifier();
   if (var_id <= d_max_variable_id) {

      if (d_variables[var_id] &&
          d_is_user_variable[var_id]) {
         TBOX_ERROR(
            "hier::VariableDatabase::registerInternalSAMRAIVariable error...\n"
            << "Attempt to register internal SAMRAI variable named "
            << variable->getName() << " with database,\n"
            << "But, that variable is already registered with the database"
            << " as a user-defined variable."
            << std::endl);
      }

   }

   bool user_variable = false;
   data_id = registerVariableAndContext_Private(variable,
         d_internal_SAMRAI_context,
         ghosts,
         user_variable);

   return data_id;

}

/*
 *************************************************************************
 *
 * Remove the given patch data index from the database if it has been
 * generated as an internal SAMRAI variable patch data index.  Also,
 * clear the index from the patch descriptor.
 *
 *************************************************************************
 */

void
VariableDatabase::removeInternalSAMRAIVariablePatchDataIndex(
   int data_id)
{
   if ((data_id >= 0) && (data_id <= d_max_descriptor_id)) {

      std::shared_ptr<Variable> variable(d_index2variable_map[data_id]);

      if (variable &&
          !d_is_user_variable[variable->getInstanceIdentifier()]) {
         removePatchDataIndex(data_id);
      }

   }
}

/*
 *************************************************************************
 *
 * Protected member function to get variable id by string name.
 *
 *************************************************************************
 */

int
VariableDatabase::getVariableId(
   const std::string& name) const
{
   int ret_id = idUndefined();

   if (!name.empty()) {
      for (int i = 0; i <= d_max_variable_id; ++i) {
         if (d_variables[i] &&
             (d_variables[i]->getName() == name)) {
            ret_id = i;
            break;
         }
      }
   }

   return ret_id;
}

/*
 *************************************************************************
 *
 * Private member functions to add contexts to database and to look up
 * context by string name.
 *
 *************************************************************************
 */

int
VariableDatabase::getContextId_Private(
   const std::string& name) const
{
   int ret_id = idUndefined();

   if (!name.empty()) {
      for (int i = 0; i <= d_max_context_id; ++i) {
         if (d_contexts[i] &&
             (d_contexts[i]->getName() == name)) {
            ret_id = i;
            break;
         }
      }
   }

   return ret_id;
}

void
VariableDatabase::addContext_Private(
   const std::shared_ptr<VariableContext>& context)
{
   int new_id = context->getIndex();
   int oldsize = static_cast<int>(d_contexts.size());
   int newsize = new_id + 1;
   if (oldsize < newsize) {
      newsize =
         tbox::MathUtilities<int>::Max(oldsize + s_context_array_alloc_size,
            newsize);
      d_contexts.resize(newsize);
   }
   d_contexts[new_id] = context;
   d_max_context_id = tbox::MathUtilities<int>::Max(d_max_context_id, new_id);
}

/*
 *************************************************************************
 *
 * Private member functions to add mapping from data index to variable
 * to the database.  Note that no error checking is done.
 *
 *************************************************************************
 */

void
VariableDatabase::addVariablePatchDataIndexPairToDatabase_Private(
   const std::shared_ptr<Variable>& variable,
   int data_id,
   bool user_variable)
{
   bool variable_added = addVariable_Private(variable, user_variable);

   if (!variable_added) {
      TBOX_ERROR("Internal hier::VariableDatabase error...\n"
         << "Attempt to add variable with duplicate name " << variable->getName()
         << " to database is not allowed.\n"
         << "Another variable with this name already exists in database."
         << std::endl);
   }

   int oldsize = static_cast<int>(d_index2variable_map.size());
   if (data_id >= oldsize) {
      d_index2variable_map.resize(
         tbox::MathUtilities<int>::Max(oldsize + s_descriptor_array_alloc_size,
            data_id + 1));
   }

   if (!d_index2variable_map[data_id] &&
       variable) {
      ++d_num_registered_patch_data_ids;
   }

   d_index2variable_map[data_id] = variable;
   d_max_descriptor_id =
      tbox::MathUtilities<int>::Max(d_max_descriptor_id, data_id);
}

void
VariableDatabase::removeVariable(
   const std::string& name)
{
   /*
    * 1. find the variable in d_variables by looking up variable id by
    *    given variable name.
    *
    * if valid id:
    * 2. unregister and/or unmap it from collection of user-defined
    *    variables in database.
    * 3. remove each patch data id associated with variable
    *    from PatchDescriptor
    * 4. remove variable from collection of variables held by database
    *    (i.e., d_variables array).
    * 5. reset max variable instance identifier if necessary
    * 6. return
    */
   int var_id = getVariableId(name);
   // if we have a valid variable id, then we'll unregister, unmap and
   // remove the variable.
   if (var_id != idUndefined()) {
      d_is_user_variable[var_id] = false;

      std::vector<int>& index_array = d_variable_context2index_map[var_id];
      for (int context_id = 0;
           context_id < static_cast<int>(index_array.size()); ++context_id) {
         if (index_array[context_id] != idUndefined()) {
            int desc_id = index_array[context_id];
            removePatchDataIndex(desc_id);
         }
      }
      // We cannot erase the item from the list, because the list's index is
      // assumed to be the instance identifier.  So, we just set this item to
      // undefined.
      d_variables[var_id].reset();
      if (var_id == d_max_variable_id) {
         --d_max_variable_id;
      }
   }
}

/*
 *************************************************************************
 *
 * Add variable to database if it doesn't already exist in the database.
 * If variable already exists in the database, do nothing.  Note that
 * we check ensure that no two distinct user-defined variables can exist
 * in the database with the same name.
 *
 *************************************************************************
 */

bool
VariableDatabase::addVariable_Private(
   const std::shared_ptr<Variable>& variable,
   bool user_variable)
{
   bool ret_value = true;

   int var_id = variable->getInstanceIdentifier();
   bool var_found = false;
   bool grow_array = false;

   if (var_id < static_cast<int>(d_variables.size())) {
      var_found = d_variables[var_id].get();
   } else {
      grow_array = true;
   }

   if (!var_found) {

      if (getVariableId(variable->getName()) != idUndefined()) {
         ret_value = false;
      }

      if (ret_value) {

         if (grow_array) {
            const int newsize =
               tbox::MathUtilities<int>::Max(static_cast<int>(d_variables.size())
                  + s_variable_array_alloc_size,
                  var_id + 1);
            d_variables.resize(newsize);
            d_variable_context2index_map.resize(newsize);

            const int oldsize = static_cast<int>(d_is_user_variable.size());
            d_is_user_variable.resize(newsize);
            for (int i = oldsize; i < newsize; ++i) {
               d_is_user_variable[i] = false;
            }
         }

         d_variables[var_id] = variable;
         d_is_user_variable[var_id] = user_variable;
         d_max_variable_id =
            tbox::MathUtilities<int>::Max(d_max_variable_id, var_id);

      }

   } // if !var_found

   return ret_value;

}

/*
 *************************************************************************
 *
 * Private member function to register variable-context pair with the
 * database and return patch data index corresponding to this pair and
 * given ghost width. The steps are:
 *
 * (1) Check whether variable-context pair maps to a valid patch data
 *     index in the database.  If it does, then check whether the
 *     index is null in the patch descriptor.  If it is, then we will
 *     create a new patch data index.  If the index is not null in
 *     the patch descriptor, the we check to see if the ghost width of
 *     that patch data index matches that in the argument list.  If the
 *     ghost width does not match, then we report an error and abort.
 *
 * (2) If we find a matching patch data index in step 1, we are done.
 *     We return the index.
 *
 * (3) If we need to create a new patch data index, do the following:
 *
 *     (3a) Create a new patch data factory, add it to the patch
 *          descriptor, and record the index.
 *
 *     (3b) We add the context to the database, if not already there.
 *
 *     (3c) We add the variable, and index to variable map to the
 *          database, if not already there.
 *
 *     (3d) We add the variable-context to index map to the database.
 *
 * (4) In the end, we return the patch data index for the
 *     variable-context pair.
 *
 *************************************************************************
 */

int
VariableDatabase::registerVariableAndContext_Private(
   const std::shared_ptr<Variable>& variable,
   const std::shared_ptr<VariableContext>& context,
   const IntVector& ghosts,
   bool user_variable)
{

   static std::string separator = "##";

   int desc_id = idUndefined();

   bool make_new_factory = true;
   int context_id = context->getIndex();
   int variable_id = variable->getInstanceIdentifier();

   // Check for valid variable_id, and get the
   // associated context to index map if valid.
   if (variable_id <= d_max_variable_id) {
      std::vector<int>& test_indx_array =
         d_variable_context2index_map[variable_id];

      // Check for valid context id and get the patch
      // descriptor id if valid.
      if (context_id < static_cast<int>(test_indx_array.size())) {
         desc_id = test_indx_array[context_id];

         // Check the descriptor id. If valid, get the associated
         // PatchDataFactory instance.
         if (desc_id != idUndefined()) {
            std::shared_ptr<PatchDataFactory> factory(
               d_patch_descriptor->getPatchDataFactory(desc_id));

            // Ensure the factory is not null and that the ghost
            // cells are the same as what we passed in.  If the ghost
            // cells aren't the same, we'll report an error and abort.
            if (factory &&
                (factory->getGhostCellWidth() != ghosts)) {
               TBOX_ERROR("hier::VariableDatabase::registerVariableAndContext"
                  << " error ...\n" << "Attempting to to register variable "
                  << variable->getName()
                  << " and context " << context->getName()
                  << " with ghost width = " << ghosts
                  << "\n This variable-context pair is already "
                  << "registered with a different ghost width. " << std::endl);
            } else {
               // reset the boolean flag if necessary
               if (factory) {
                  make_new_factory = false;
               }
            }
         }     // if (desc_id != idUndefined())
      }     // if (context_id < test_indx_array.size())
   }     // if (variable_id <= d_max_variable_id)

   // Create the new factory if necessary
   if (make_new_factory) {

      std::shared_ptr<PatchDataFactory> new_factory(
         variable->getPatchDataFactory()->cloneFactory(ghosts));

      std::string tmp(variable->getName());
      tmp += separator;
      tmp += context->getName();
      desc_id = d_patch_descriptor->definePatchDataComponent(tmp, new_factory);

      addContext_Private(context);

      addVariablePatchDataIndexPairToDatabase_Private(variable,
         desc_id,
         user_variable);

      std::vector<int>& var_indx_array =
         d_variable_context2index_map[variable_id];
      int oldsize = static_cast<int>(var_indx_array.size());
      int newsize = context_id + 1;
      if (oldsize < newsize) {
         var_indx_array.resize(newsize);
         for (int i = oldsize; i < newsize; ++i) {
            var_indx_array[i] = idUndefined();
         }
      }
      var_indx_array[context_id] = desc_id;

   }  // if (make_new_factory)

   return desc_id;

}

}
}
