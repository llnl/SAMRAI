/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Singleton database class for managing variables and contexts.
 *
 ************************************************************************/

#ifndef included_hier_VariableDatabase
#define included_hier_VariableDatabase

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"

#include <string>
#include <iostream>
#include <vector>
#include <memory>

namespace SAMRAI {
namespace hier {

/*
 *  TODO: Should we be more explicit in this class documentation to describe
 *  how variables of different dimensions are handled?  For example, this
 *  places constraints on several methods (i.e., variable dimension must
 *  match the ghost width dimension, etc.).  The Array containers that we use
 *  do not know about the dimension associated with the objects they contain.
 *  Should we manage things more explicitly for each dimension?
 */
/*!
 * @brief Generates and manages mappings between
 * Variables and/or (Variable,VariableContext) pairs and patch data indices.
 *
 * This class is a Singleton, and serves as a globally accessible lookup
 * table for the bookkeeping information related to variable storage on a
 * SAMRAI patch hierarchy.  The database may be used to manage the mapping
 * between a Variable or (Variable,VariableContext) pair and the patch data
 * index associated with the data for the Variable or
 * (Variable,VariableContext) pair.
 *
 * Typically, numerical routines or solution algorithms manage and access
 * data in terms of Variables and/or (Variable,VariableContext) pairs
 * (e.g., "OLD", "NEW").
 *
 * @par Definitions
 *
 * <b><em> %Patch data index:</em></b> SAMRAI uses a patch data index value (an
 * integer) to lookup patch data information in the patch descriptor
 * object that is shared by all patches in a SAMR hierarchy.
 *
 * <b><em> %Variable:</em></b>  A storage management class.  Each
 * concrete implementation (cell-centered, node-centered, etc.) provides
 * the appropriate storage management for that type.
 *
 * <b><em>%VariableContext:</em></b>  A way to consistently address the state
 * of a variable in multiple contexts.  For example, an application may need
 * to manage data for "OLD" pressure values as well as "NEW" pressure
 * values.  The context allows the user to have a single variable in multiple
 * contexts, without having to create multiple variables.  A context is
 * optional and a user could create multiple variables; e.g., "pressure_old"
 * and "pressure_new".
 *
 * When the association between Variable and VariableContext is needed,
 * the registration operations are used to create this association.
 * Each variable and context pair maps to a unique integer
 * patch data index.  The integer indices are the same for every patch in the
 * hierarchy so data objects are accessed identically for all patches.
 *
 * The following example shows the mapping between a PatchData indices
 * and a Variable/VariableContext pair.
 *
 * @code
 *    PatchData Index   Variable   Context
 *    ================  ========   =======
 *    index = 0         uval       CURRENT
 *    index = 1         flux       CURRENT
 *    index = 2         flux       SCRATCH
 * @endcode
 *
 * @par Usage
 * The typical way to use the database is as follows:
 *
 * <ol>
 *   <li>Add Variable to the database
 *   <li>Get or create the VariableContext via the getContext() function
 *   <li>Register the VariableContext pair in the database
 * </ol>
 *
 *  @code
 *
 *    // other setup including dimension and IntVectors specifying ghosts
 *    ...
 *
 *    // Get the singleton variable database.
 *    VariableDatabase* var_db = VariableDatabase.getDatabase();
 *
 *    // Create the context (or access it if it exists already).
 *    std::shared_ptr<VariableContext> current =
 *       var_db->getContext("CURRENT");
 *
 *    // Define the variable.
 *    std::shared_ptr<pdat::FaceVariable<double> > flux(
 *       new pdat::FaceVariable<double>(dim, "flux", 1));
 *
 *    // Register the Variable/VariableContext pair.
 *    // "ghosts" is an IntVector defining the data ghost width for this
 *    // usage context.  Note that other contexts may have different ghost
 *    // widthsi for the same variable.
 *    const int flux_current_id =
 *       var_db->registerVariableAndContext(flux, current, ghosts);
 *
 *  @endcode
 *
 * @par
 * A Variable can also be added to the database using the
 * addVariable() function.  This function does not register a mapping
 * between a Variable and VariableContext, rather it is used when
 * there are multiple objects needing to access a single Variable.
 * This eliminates issues with shared objects, order of creation, etc.
 * For example, two different algorithms in an
 * application may need to access the current value for pressure.
 * Each of these algorithms may use the database to access the
 * patch data index of this shared state.  This is not needed if
 * registerVariableAndContext() is used.
 *
 * @par Alternate Usage
 * The database can be used to maintain a mapping between a single patch data
 * index and a Variable.  This type of mapping is constructed using either the
 * registerPatchDataIndex() function or registerClonedPatchDataIndex() function.
 *
 * <ul>
 *   <li> Add variable and index pair to the database via the
 *       addVariablePatchDataIndex() function.   Or, clone an existing
 *       variable-index pair via registerClonedPatchDataIndex().
 *       Either of these functions will check to make sure that the type
 *       of the variable matches the given patch data index.
 *   <li> Get variable for the index via the mapIndexToVariable() function.
 * </ul>
 *
 * @par
 * The database is also used in SAMRAI to manage the reading/writing of
 * patch data objects to/from restart files.
 *
 * @note
 * <ul>
 *    <li> Variable names are unique in the database
 *    <li> VariableContext names are unique in the database
 *    <li> VariableContext should be generated solely through the use of the
 *         getContext() function to avoid unexpected results; e.g.,
 *         multiple context objects with the same name.
 *    <li> A Variable can be removed from the database using removeVariable().
 *    <li> The database can free PatchData indices using removepatchDataIndex()
 *         function.
 *    <li> A VariableContext persists once created and cannot be removed from
 *         the database.
 *    <li> The database is reset (all data wiped) during shutdown.  It is
 *         started from a clean state during a new startup.
 * </ul>
 *
 *
 * @see PatchDescriptor
 * @see VariableContext
 * @see Variable
 * @see Patch
 */

class VariableDatabase
{
public:
   /*!
    * @brief Return a pointer to the singleton variable database instance.
    *
    * @note
    * When the database is accessed for the first time, the
    * Singleton instance is registered with the StartupShutdownManager
    * class which destroys such objects at program completion.  Thus,
    * this class is not explicitly allocated or deallocated.
    *
    * @return  Bare pointer to variable database instance.
    */
   static VariableDatabase *
   getDatabase();

   /*!
    * @brief Get the number of patch data indices registered with the database.
    */
   virtual int
   getNumberOfRegisteredPatchDataIndices() const;

   /*!
    * @brief Get number of variable contexts registered with the database.
    */
   virtual int
   getNumberOfRegisteredVariableContexts() const;

   /*!
    * @brief Get the patch descriptor managed by the database.
    *
    * This descriptor is shared by all patches in the hierarchy.
    */
   virtual std::shared_ptr<PatchDescriptor>
   getPatchDescriptor() const;

   /*!
    * @brief Creates or returns a VariableContext object with the given name.
    *
    * If a context exists in the database with the given name, it is
    * returned.  Otherwise, a new context is created and returned.
    *
    * @note
    * It is impossible to add two distinct variable context
    * objects to the database with the same name.
    *
    * The name must not be empty; a null pointer will be return in this case.
    *
    * @param[in] context_name
    *
    * @return  Variable context.
    */
   virtual std::shared_ptr<VariableContext>
   getContext(
      const std::string& context_name);

   /*!
    * @brief Check whether context with given name exists in the database.
    *
    * @param[in] context_name
    */
   virtual bool
   checkContextExists(
      const std::string& context_name) const;

   /*!
    * @brief Add the Variable to the database.
    *
    * If the given Variable already exists nothing is done; the same
    * Variable may be added multiple times.  If a Variable exists with
    * the same name identifier but is a different Variable object, an
    * error will be logged and the program will abort.  This prevents
    * multiple Variables being associated with the same name.
    *
    * @param[in] variable std::shared_ptr to variable
    *
    * @pre variable
    */
   virtual void
   addVariable(
      const std::shared_ptr<Variable>& variable);

   /*!
    * @brief Remove the Variable from the database identified by @c name.
    *
    * @param[in] variable_name
    */
   virtual void
   removeVariable(
      const std::string& variable_name);

   /*!
    * @brief Get variable in database with given name string identifier.
    *
    * @param[in] variable_name
    *
    * @return  Variable in the database with given name.
    *          If no such variable exists, a null pointer is returned.
    */
   virtual std::shared_ptr<Variable>
   getVariable(
      const std::string& variable_name) const;

   /*!
    * @brief Check whether a variable with given name exists in the database.
    *
    * @param[in] variable_name
    *
    * @return  True if variable with name exists in database;
    *          otherwise, false.
    */
   virtual bool
   checkVariableExists(
      const std::string& variable_name) const;

   /*!
    * @brief Create and register a new patch data index by cloning the data
    * factory associated with the Variable at the index provided.
    *
    * The new index and variable pair is added to the VariableDatabase.
    * A variable-patch data index pair generated using this function cannot
    * be looked up using a VariableContext. If the @c old_id is invalid
    * or undefined, or does not map to patch data of the same type as
    * @c variable, the program will abort with an error message.
    *
    * @note
    * This function does not deallocate any patch data storage associated
    * with the new patch data index.
    *
    * @param[in]  variable std::shared_ptr to @c Variable.  If the variable
    *             is unknown to the database, then an invalid patch data index
    *             (< 0) will be returned
    * @param[in]  old_id Integer patch data index currently associated with
    *             variable. If this value is not a valid patch data index
    *             (< 0) or does not map to patch data matching the
    *             type of the given variable, the program will abort with an
    *             error message.
    *
    * @return New integer patch data index. If new patch data not added,
    *         return value is an invalid (undefined) patch data index (< 0).
    *
    * @pre variable
    * @pre checkVariablePatchDataIndex(variable, old_id)
    */
   virtual int
   registerClonedPatchDataIndex(
      const std::shared_ptr<Variable>& variable,
      int old_id);

   /*!
    * @brief Add given patch data index and variable pair to the database.
    *
    * This registration function is primarily intended for Variable objects
    * (i.e., DO NOT use for internal SAMRAI variables) that are not
    * associated with a VariableContext and for which a patch data index
    * is already known.
    *
    * @par [Default Case]
    * If the index is unspecified, the default variable factory is cloned
    * and the Variable and new index are added to the database.
    * In this case, the patch data will have the default
    * ghost associated with the given @c variable (defined
    * by the patch data factory it generates).
    *
    * @note
    * <ul>
    *   <li> variable-patch data index pair generated with this function
    *        cannot be looked up using a VariableContext.
    *   <li> This function does not allocate any patch data storage associated
    *        with the integer index.
    *   <li> This function must not be used by SAMRAI developers for
    *        creating patch data indices for internal SAMRAI variables.  The
    *        routine registerInternalSAMRAIVariable() must be used for that
    *        case.
    * </ul>
    * @param[in]  variable std::shared_ptr to Variable
    * @param[in]  data_id  Optional integer patch data index to be added
    *                      (along with variable) to the database.  If the value
    *                      is unspecified (default case), the default variable
    *                      patch data factory is used to generate a new factory.
    *                      If the value is provided and does not map to patch
    *                      data matching the type of the given variable, the
    *                      program will abort with an error message.
    *
    * @return New integer patch data index.  If new patch data index not
    *         added, return value is an invalid patch data index (< 0).
    *
    * @pre variable
    * @pre data_id == idUndefined() || checkVariablePatchDataIndex(variable, data_id)
    */
   virtual int
   registerPatchDataIndex(
      const std::shared_ptr<Variable>& variable,
      int data_id = idUndefined());

   /*!
    * @brief Remove the patch data index from the VariableDatabase if it exists.
    *
    * @par Side Effects
    * This function also removes the given index from the patch descriptor
    * and any mapping between the index and a variable from the
    * database.
    *
    * @note
    * This function does not deallocate any patch data storage associated
    * with the integer index.
    *
    * @param[in]  data_id  Integer patch data index to be removed from
    *                  the database
    */
   virtual void
   removePatchDataIndex(
      int data_id);

   /*!
    * @brief Check whether the given variable is mapped to the given patch data
    * index in the database.
    *
    * @param[in]  variable  std::shared_ptr to variable
    * @param[in]  data_id   Integer patch data index
    *
    * @return  Boolean true if the variable is mapped the given patch
    * data index; false otherwise.
    *
    * @pre variable
    * @pre (data_id >= 0) &&
    *      (data_id < getPatchDescriptor()->getMaxNumberRegisteredComponents())
    */
   virtual bool
   checkVariablePatchDataIndex(
      const std::shared_ptr<Variable>& variable,
      int data_id) const;

   /*!
    * @brief Check whether the given variable matches the patch data type
    * associated with the given patch data index in the database.
    *
    * @param[in]  variable  std::shared_ptr to variable
    * @param[in] data_id   Integer patch data index
    *
    * @return  Boolean true if the type of the variable matches the type of
    *          the patch data at the given patch data index; false otherwise.
    *
    * @pre variable
    * @pre (data_id >= 0) &&
    *      (data_id < getPatchDescriptor()->getMaxNumberRegisteredComponents())
    */
   virtual bool
   checkVariablePatchDataIndexType(
      const std::shared_ptr<Variable>& variable,
      int data_id) const;

   /*!
    * @brief Register variable and context pair along with the ghost
    * width for the patch data mapped to the (Variable, VariableContext)
    * pair with the variable database.
    *
    * @par
    * Typically, this function will generate a new patch data
    * index for the variable and ghost width and add the
    * variable-context pair and index to the database.  If the
    * variable-context pair is already mapped to some patch data index
    * in the database, and the given ghost width matches that of the
    * patch data, then that index will be returned and the function will
    * do nothing.  However, if the variable-context pair is already mapped
    * to some patch data index with a different ghost width, the program
    * will abort with a descriptive error message.
    *
    * @par
    * If either the variable or the context is unknown to the database
    * prior to calling this routine, both items will be added to the
    * database, if possible.  The constraints for the getContext() and
    * addVariable() routines apply.
    *
    * @note
    * It is an error to map a (Variable,VariableContext) pair plus data
    * index having a different ghost width than that passed in the
    * argument list.  The program will abort with a descriptive message
    * in this case.
    *
    *
    * @param[in]  variable  std::shared_ptr to variable
    * @param[in] context    std::shared_ptr to variable context
    * @param[in] ghosts     Optional ghost width for patch data associated
    *                       with variable-context pair.
    *
    * @return Integer patch data index of variable-context pair in database.
    *
    * @pre variable
    * @pre context
    * @pre ghosts.min() >= 0
    */
   virtual int
   registerVariableAndContext(
      const std::shared_ptr<Variable>& variable,
      const std::shared_ptr<VariableContext>& context,
      const IntVector& ghosts // NOTE: old default (zero ghost width)
                              // does not work since dimension of
                              // variable and IntVector must match.
      );

   /*!
    * @brief Map variable-context pair in database to patch data index.
    *
    * If there is no such pair in the database (either the variable does
    * not exist, the context does not exist, or the pair has not been
    * registered), then an invalid patch data index (i.e., < 0) is returned.
    *
    * @note
    * For this function to operate as expected, the database mapping
    * information must have been generated using the
    * registerVariableAndContext() function.  If the variable was
    * registered without a variable context, then the patch data index
    * associated with the variable will not be returned.  See the other
    * map...() functions declared in this class.
    *
    * @param[in]  variable  std::shared_ptr to variable
    * @param[in]  context   std::shared_ptr to variable context
    *
    * @return Integer patch data index of variable-context pair in database.
    *         If the variable-context pair was not registered with the
    *         database, then an invalid data index (< 0) will be returned.
    *
    * @pre variable
    * @pre context
    */
   virtual int
   mapVariableAndContextToIndex(
      const std::shared_ptr<Variable>& variable,
      const std::shared_ptr<VariableContext>& context) const;

   /*!
    * @brief Map patch data index to variable associated with the data, if
    * possible, and set the variable pointer to the variable in the database.
    *
    * @param[in]   index  Integer patch data index
    * @param[out]  variable  std::shared_ptr to variable that maps to patch
    *                    data index in database.  If there is no index in the
    *                    database matching the index input value, then the
    *                    variable pointer is set to null.
    *
    * @return  Boolean true if patch data index maps to variable in the
    *          database; otherwise false.
    */
   virtual bool
   mapIndexToVariable(
      const int index,
      std::shared_ptr<Variable>& variable) const;

   /*!
    * @brief Map patch data index to variable-context pair associated with
    * the data, if possible, and set the variable and context pointers to
    * the corresponding database entries.
    *
    * @note
    * For this function to operate as expected, the database
    * mapping information must have been generated using the
    * registerVariableAndContext() function.  If the variable was
    * registered without a variable context, then the variable and
    * variable context returned may not be what is expected by the
    * user; e.g., they may be associated with internal SAMRAI
    * variables.
    *
    * @param[in]   index patch data index
    * @param[out]   variable std::shared_ptr to variable set to matching
    *          variable in database.  If no match is found, it is set to null.
    * @param[out]   context  std::shared_ptr to variable context set to
    *          matching variable context in database. If no match is found, it
    *          is set to null.
    *
    * @return  Boolean true if patch data index maps to variable-context
    *          pair in the database; otherwise false.
    */
   virtual bool
   mapIndexToVariableAndContext(
      const int index,
      std::shared_ptr<Variable>& variable,
      std::shared_ptr<VariableContext>& context) const;

   /*!
    * @brief Print variable, context, and patch descriptor information
    * contained in the database to the specified output stream.
    *
    * @param[in] os  Optional output stream.  If not given, tbox::plog is used.
    * @param[in] print_only_user_defined_variables Optional boolean value
    *        indicating whether to print information for all variables
    *        in database or only those that are associated with user-
    *        defined quantities; i.e., not internal SAMRAI variables.
    *        The default is true, indicating that only user-defined
    *        variable information will be printed.
    */
   virtual void
   printClassData(
      std::ostream& os = tbox::plog,
      bool print_only_user_defined_variables = true) const;

   /*!
    * @brief Register internal SAMRAI variable and ghost width
    * with the variable database.
    *
    * This function will generate a new patch data index for the variable
    * and ghost width unless the variable is already mapped to some
    * patch data index in the database with a different ghost width
    * or as a user-defined variable.  If the variable is unknown to the
    * database prior to calling this routine, it will be added to the
    * database.
    *
    * @note
    * This routine is intended for managing internal SAMRAI
    * work variables that are typically unseen by users.  It should not be
    * called by users for registering variables, or within SAMRAI for
    * registering any user-defined variables with the variable database.
    * This function enforces the same constraints on variable registration
    * that are applied for registering user-defined variables; e.g., using
    * the routine registerVariableAndContext().  Thus, it is the
    * responsibility of SAMRAI developers to avoid naming and ghost
    * width conflicts with user-defined variables or other internal
    * SAMRAI variables.
    *
    *
    * @param[in]  variable  std::shared_ptr to variable
    * @param[in] ghosts     Ghost width for patch data associated with the
    *                       variable
    * @return Integer patch data index of variable-ghost width pair
    * in database.
    *
    * @pre variable
    * @pre ghosts.min() >= 0
    */
   virtual int
   registerInternalSAMRAIVariable(
      const std::shared_ptr<Variable>& variable,
      const IntVector& ghosts);

   /*!
    * @brief Remove the given index from the variable database if it exists in
    * the database and is associated with an internal SAMRAI variable registered
    * with the function registerInternalSAMRAIVariable().
    *
    * @par
    * Also, remove the given index from the patch descriptor and
    * remove any mapping between the index and a variable from the
    * variable database.
    *
    * @note
    * This function does not deallocate
    * any patch data storage associated with the integer index.
    *
    * @note
    * This routine is intended for managing internal SAMRAI
    * work variables that are typically unseen by users.  It should not be
    * called by users for removing patch data indices, or within SAMRAI for
    * removing any patch data indices associated with user-defined variables.
    *
    * @note
    * The given index will not be removed if is not associated with
    * an internal SAMRAI variable in the variable database; i.e., a
    * user-defined variable.
    *
    * @param[in]  data_id  Integer patch data identifier to be removed from
    *                  the database.
    */
   virtual void
   removeInternalSAMRAIVariablePatchDataIndex(
      int data_id);

protected:
   /**
    * @brief The constructor for VariableDatabase is protected.
    * Consistent with the definition of a Singleton class, only the
    * database object has access to the constructor for the class.
    *
    * The constructor initializes the state of database contents.
    */
   VariableDatabase();

   /**
    * @brief The destructor for VariableDatabase is protected. See the
    * comments for the constructor.
    *
    * The destructor deallocates database contents.
    */
   virtual ~VariableDatabase() = default;

   /**
    * @brief Return integer value used to indicate undefined variable or
    * context identifier.  This routine is protected to allow subclasses
    * to be consistent with this database class.
    */
   static int
   idUndefined()
   {
      return -1;
   }

   /**
    * @brief Return integer identifier for first variable found matching given
    * string name identifier (care must be used to avoid adding variables to
    * database with same name to insure correct behavior), or return an
    * undefined integer index if no such variable exists in the database.
    */
   int
   getVariableId(
      const std::string& name) const;

   /**
    * @brief Initialize Singleton instance with instance of subclass.  This
    * function is used to make the singleton object unique when inheriting
    * from this base class.
    *
    * @pre !s_variable_database_instance
    */
   void
   registerSingletonSubclassInstance(
      VariableDatabase* subclass_instance);

private:
   /*!
    * @brief Deallocate the Singleton VariableDatabase instance.
    *
    * It is not necessary to call this function at program termination,
    * since it is automatically called by the StartupShutdownManager class.
    */
   static void
   shutdownCallback();

   /*
    * Private member functions to search for indices of contexts
    * given a string identifier, and to add a new context
    * to database.
    */
   int
   getContextId_Private(
      const std::string& name) const;

   void
   addContext_Private(
      const std::shared_ptr<VariableContext>& context);

   /*
    * Private member function to add variable to database (either
    * user-defined or internal SAMRAI variable depending on boolean
    * argument).  Boolean return value is true if variable is either
    * added to variable database or already found in variable database.
    * Boolean return value is false only when a user variable has
    * the same name string identifier of another variable already in
    * the database.  In this case, the variable is not added to the database.
    */
   bool
   addVariable_Private(
      const std::shared_ptr<Variable>& variable,
      bool user_variable);

   /*
    * Private member function to add variable-patch data index mapping
    * to the database.  Boolean indicates whether variable is a user-defined
    * variable or an internal SAMRAI work variable.
    */
   void
   addVariablePatchDataIndexPairToDatabase_Private(
      const std::shared_ptr<Variable>& variable,
      int data_id,
      bool user_variable);

   /*
    * Private member function to add variable-context pair to the database.
    * Boolean indicates whether variable is a user-defined variable or
    * an internal SAMRAI work variable.
    */
   int
   registerVariableAndContext_Private(
      const std::shared_ptr<Variable>& variable,
      const std::shared_ptr<VariableContext>& context,
      const IntVector& ghosts,
      bool user_variable);

   /*
    * Static data members used to control allocation of arrays.
    */
   static const int s_context_array_alloc_size;
   static const int s_variable_array_alloc_size;
   static const int s_descriptor_array_alloc_size;

   /*
    * Data members that store variable, context, patch data index information.
    */
   std::shared_ptr<PatchDescriptor> d_patch_descriptor;

   std::shared_ptr<VariableContext> d_internal_SAMRAI_context;

   int d_num_registered_patch_data_ids;

   /*
    * Vector of VariableContext pointers is indexed as
    * d_contexts[ <context id> ]
    */
   int d_max_context_id;
   std::vector<std::shared_ptr<VariableContext> > d_contexts;

   /*
    * Vector of Variable pointers is indexed as d_variables[ <variable id> ]
    */
   int d_max_variable_id;
   std::vector<std::shared_ptr<Variable> > d_variables;

   /*
    * Vector of VariableContext to patch descriptor indices is indexed as
    * d_variable_context2index_map[ <context id> ]
    */
   std::vector<std::vector<int> > d_variable_context2index_map;

   /*
    * Vector of patch descriptor indices to Variables is indexed as
    * d_index2variable_map[ <descriptor id> ]
    */
   int d_max_descriptor_id;
   std::vector<std::shared_ptr<Variable> > d_index2variable_map;

   /*
    * Vector of user variable booleans is indexed as
    * d_is_user_variable[ <variable id> ]
    */
   std::vector<bool> d_is_user_variable;

   static tbox::StartupShutdownManager::Handler
      s_shutdown_handler;
};

}
}

#endif
