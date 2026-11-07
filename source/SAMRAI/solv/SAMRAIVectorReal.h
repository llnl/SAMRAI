/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Vector class for real data on SAMRAI hierarchy.
 *
 ************************************************************************/

#ifndef included_solv_SAMRAIVectorReal
#define included_solv_SAMRAIVectorReal

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/PIO.h"

#include <string>
#include <iostream>
#include <memory>

namespace SAMRAI {
namespace solv {

/**
 * Class SAMRAIVectorReal allows a collection of patch data types
 * (with double or float data, but not both) defined over a SAMRAI hierarchy
 * to be manipulated as though they are all part of a single vector.
 * Specifically, this class provides a set of common vector operations
 * to manipulate all of the data components as a whole.  The most obvious
 * use of this class is in SAMRAI applications that use solver libraries,
 * such as KINSOL, CVODE, or PETSc.  Specific vector objects that can be
 * used with these packages are defined elsewhere in SAMRAI.  However, all
 * these vactor interfaces are built using this vector class.
 *
 * This class defines a vector to be any collection of patch data objects
 * (either cell-, edge-, face-, node, or side-centered, or any combination
 * of these) defined over a set of patch levels in an AMR hierarchy.  All of
 * the data objects must have the same underlying data type, either double
 * or float.  The vector structure is composed by adding individual
 * variable quantities to the vector after it is constructed.  When a
 * component is added, a weighting or "control volume" component
 * (having the same type as the vector component) may also be added
 * to the vector.  These weights are used to define the contribution of
 * each vector entry to summing operations such as norms and dot products.
 * For example, the weights can be used to mask out coarse level vector
 * data entries in cells that are covered by fine cells when the coarse data
 * are not actually part of the solution vector.  The weights can also be
 * used to map the vector operations to grid-based operations that define
 * control volume weights.   It is important to note that the centering
 * of each control volume component must match that of the vector component
 * with which it is associated.
 *
 * Typical usage of this vector class is as follows:
 *
 * - @b (1) Construct a vector instance by specifying the patch hierarchy
 *            and range of levels over which the vector is defined.
 *            The levels must exist in the hierarchy before the vector
 *            can be used or an assertion will result.  However, a vector
 *            may be created before the levels exist.  The range of levels
 *            can be reset (such as after remeshing) by calling the
 *            resetLevels() function.
 * - @b (2) Register each data component with the vector by providing
 *            the variable and its storage location (i.e., patch data
 *            index), and the control volume index if needed.
 *            See the addComponent() functions.
 * - @b (3) Manipulate data using vector operations.
 *
 *
 * Before the vector operations can be used, the storage for each of its
 * components must be allocated.  Storage allocation is only possible
 * through a vector object after all component variables are added to
 * the vector (using the addComponent() function).
 * Then, the allocateVectorData() function will allocate storage for all
 * components when called.  Alternatively, patch data objects (corresponding
 * to the variables and vector patch data indices) may be
 * explicitly created elsewhere.  However, depending on the circumstance, this
 * second alternative may be more confusing and require more bookkeeping on
 * the user's part.  See the documentation accompanying the addComponent()
 * function for more information.
 *
 * @see math::HierarchyDataOpsReal
 */

template<class TYPE>
class SAMRAIVectorReal
{
public:
   /**
    * Constructor for SAMRAIVectorReal class is used to construct
    * each unique vector within an application.  That is, each vector
    * that is used to represent a unique set of variable quantities is
    * considered unique.  This constructor is used to create a solution
    * vector for an application or solver algorithm.  The cloneVector()
    * function is provided to generate copies of a given vector.  For example,
    * the clone function may be used by a solver to generate copies of
    * the vector as needed; e.g., in a Krylov subspace method like GMRES.
    *
    * Before the vector may be used, data components must be added to it using
    * the adddComponent(0 function.  Also, this constructor does not allocate
    * storage for vector data.  This is usually done after all components are
    * added.  The allocateVectorData() function is used for this
    * purpose.  Otherwise, existing patch data quantities can be
    * added as vector components.  In any case, storage for all
    * components must be allocated before the vector can be used.
    *
    * The range levels can be reset at any time (e.g., if the level
    * configuration changes by re-meshing), by calling the resetLevels() member
    * function.
    *
    * Although an empty std::string may be passed as the vector name, it is
    * recommended that a descriptive name be used to facilitate debugging
    * and error reporting.
    *
    * By default the vector component information and data will be sent to
    * the "plog" output stream when the print() function is called.  This
    * stream can be changed at any time via the setOutputStream() function.
    *
    * @pre hierarchy
    * @pre (coarsest_level >= 0) && (finest_level >= coarsest_level) &&
    *      (finest_level <= hierarchy->getFinestLevelNumber())
    */
   SAMRAIVectorReal(
      const std::string& name,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level);

   /**
    * Virtual destructor for SAMRAIVectorReal class.  The destructor
    * destroys all vector component information.  However, the destructor
    * does not deallocate the vector component storage, nor does it return the
    * vector patch data indices to the patch descriptor free list.
    * The freeVectorComponents() function is provided for this task.  The
    * reason for this is that an application may create a vector based on
    * some pre-existing patch data objects that must live beyond the
    * destruction of the vector object.
    */
   virtual ~SAMRAIVectorReal();

   /**
    * Set std::string identifier for this vector object.
    */
   void
   setName(
      const std::string& name);

   /**
    * Set output stream for vector object.  When the print() function is
    * called, all vector data will be sent to the given output stream.
    */
   void
   setOutputStream(
      std::ostream& s);

   /**
    * Return reference to the output stream used by this vector object.
    * This function is primarily used by classes which define interfaces
    * between this vector class and vector kernels defined by other
    * packages.  Specifically, SAMRAI vectors and package-specific wrappers
    * for those vectors may all access the same output stream.
    */
   std::ostream&
   getOutputStream();

   /**
    * Reset range of patch levels over which vector is defined.  This
    * function resets the data operations for all vector components.
    * Note that the levels must exist in the hierarchy when this function
    * is called or a non-recoverable assertion will result.
    */
   void
   resetLevels(
      const int coarsest_level,
      const int finest_level);

   /**
    * Return std::string identifier for this vector object.
    */
   const std::string&
   getName() const;

   /**
    * Return pointer to patch hierarchy associated with the vector.
    */
   std::shared_ptr<hier::PatchHierarchy>
   getPatchHierarchy() const;

   /**
    * Return integer number of coarsest hierarchy level for vector.
    */
   int
   getCoarsestLevelNumber() const;

   /**
    * Return integer number of finest hierarchy level for vector.
    */
   int
   getFinestLevelNumber() const;

   /**
    * Return integer number of patch data components in vector.
    */
   int
   getNumberOfComponents() const;

   /**
    * Return patch data object for given vector component index.
    *
    * @pre (comp_id >= 0) && (comp_id < getNumberOfComponents())
    */
   std::shared_ptr<hier::PatchData>
   getComponentPatchData(
      const int comp_id,
      const hier::Patch& patch) const;

   /**
    * Return patch data object associated with given variable.
    *
    * @pre var
    * @pre d_variableid_2_vectorcomponent_map[var->getInstanceIdentifier()] >= 0
    */
   std::shared_ptr<hier::PatchData>
   getComponentPatchData(
      const std::shared_ptr<hier::Variable>& var,
      const hier::Patch& patch) const;

   /**
    * Return pointer to variable for specified vector component.
    *
    * @pre (component >= 0) && (component < getNumberOfComponents())
    */
   std::shared_ptr<hier::Variable>
   getComponentVariable(
      const int component) const;

   /**
    * Return patch data index for specified vector component.
    *
    * @pre (component >= 0) && (component < getNumberOfComponents())
    */
   int
   getComponentDescriptorIndex(
      const int component) const;

   /**
    * Return patch data index of control volume data for vector component.
    *
    * @pre (component >= 0) && (component < getNumberOfComponents())
    */
   int
   getControlVolumeIndex(
      const int component) const;

   /**
    * Clone this vector object and return a pointer to the vector copy
    * (i.e., a new vector).  Each patch data component in the new vector
    * will match the corresponding component in this vector object.  However,
    * the data for the components of the new vector will be assigned to
    * different patch data indices than the original.  In short,
    * the cloned vector will have an identical structure to the original,
    * but its data storage will be distinct.  Before the new vector object
    * can be used, its data must be allocated explicitly.
    *
    * Note that this function maps the variables associated with the new
    * vector to the new vector component data indices (i.e., patch
    * data indices) in the variable database.  Thus the mapping
    * between variables and patch data for the new vector can be obtained
    * from the variable database if needed.
    *
    * If an empty std::string is passed in, the name of this vector object
    * is used for the new vector.
    */
   std::shared_ptr<SAMRAIVectorReal<TYPE> >
   cloneVector(
      const std::string& name) const;

   /**
    * Destroy the storage corresponding to the vector components and
    * free the associated patch data entries from the variable database
    * (which will also clear the indices from the patch descriptor).
    */
   void
   freeVectorComponents();

   /**
    * Add a new variable and patch data component to this vector.
    * The integer values passed in represent the patch data indices
    * for the vector component data and the component control volume data.
    * If the control volume patch data index is not specified
    * (i.e., control_vol_id < 0), no weighting will be applied in
    * vector operations associated with the component.  This routine
    * also accepts a hierarchy data operation object for the component
    * should the user want to provide a special set of such operations.
    * If left unspecified (nearly all cases), the standard operations
    * for the given variable type are used.
    *
    * Note that this function maps the variable to the component data
    * index (i.e., patch data index) in the variable database.
    * Thus, the mapping between the variable and its patch data for the
    * vector can be obtained from the variable database if needed.
    *
    * @pre hier::VariableDatabase::getDatabase()->checkVariablePatchDataIndexType(var, comp_data_id)
    * @pre (comp_vol_id < 0) || (hier::VariableDatabase::getDatabase()->checkVariablePatchDataIndexType(var, comp_vol_id))
    */
   void
   addComponent(
      const std::shared_ptr<hier::Variable>& var,
      const int comp_data_id,
      const int control_vol_id = -1,
      const std::shared_ptr<math::HierarchyDataOpsReal<TYPE> >& vop =
         std::shared_ptr<math::HierarchyDataOpsReal<TYPE> >());

   /**
    * Allocate data storage for all components of this vector object.
    *
    * @pre getPatchHierarchy
    * @pre (getCoarsestLevelNumber() >= 0) &&
    *      (getFinestLevelNumber() >= getCoarsestLevelNumber()) &&
    *      (getFinestLevelNumber( <= d_hierarchy->getFinestLevelNumber())
    */
   void
   allocateVectorData(
      const double timestamp = 0.0);

   /**
    * Deallocate data storage for all components of this vector object.
    * Note that this routine will not free the associated data
    * indices in the patch descriptor.  See freeVectorComponents() function.
    *
    * @pre getPatchHierarchy
    * @pre (getCoarsestLevelNumber() >= 0) &&
    *      (getFinestLevelNumber() >= getCoarsestLevelNumber()) &&
    *      (getFinestLevelNumber( <= d_hierarchy->getFinestLevelNumber())
    */
   void
   deallocateVectorData();

   /**
    * Print component information and data for this vector object.
    */
   void
   print(
      std::ostream& s = tbox::plog,
      const bool interior_only = true) const;

   /**
    * Copy data from source vector components to components of this vector.
    */
   void
   copyVector(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& src_vec,
      const bool interior_only = true);

   /**
    * Swap data components (i.e. storage) between this vector object and
    * argument vector.
    */
   void
   swapVectors(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& other);

   //@{
   /*!
    * @name Vector arithmetic functions
    */

   /**
    * Set all components of this vector to given scalar value.
    */
   void
   setToScalar(
      const TYPE& alpha,
      const bool interior_only = true);

   /**
    * Set this vector to src vector multiplied by given scalar.
    */
   void
   scale(
      const TYPE& alpha,
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
      const bool interior_only = true);

   /**
    * Set this vector to sum of given vector and scalar.
    */
   void
   addScalar(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
      const TYPE& alpha,
      const bool interior_only = true);

   /**
    * Set this vector to sum of two given vectors.
    */
   void
   add(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& y,
      const bool interior_only = true);

   /**
    * Set this vector to difference of two given vectors (i.e., x - y).
    */
   void
   subtract(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& y,
      const bool interior_only = true);

   /**
    * Set each entry in this vector to product of corresponding entries in
    * input vectors.
    */
   void
   multiply(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& y,
      const bool interior_only = true);

   /**
    * Set each entry in this vector to ratio of corresponding entries in
    * input vectors (i.e., x / y).  No check for division by zero.
    */
   void
   divide(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& y,
      const bool interior_only = true);

   /**
    * Set each entry of this vector to reciprocal of corresponding entry
    * in input vector.  No check is made for division by zero.
    */
   void
   reciprocal(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
      const bool interior_only = true);

   /**
    * Set this vector to the linear sum @f$ \alpha x + @beta y @f$ , where
    * @f$ \alpha, @beta @f$  are scalars and @f$ x, y @f$  are vectors.
    */
   void
   linearSum(
      const TYPE& alpha,
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
      const TYPE& beta,
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& y,
      const bool interior_only = true);

   /**
    * Set this vector to the sum @f$ \alpha x + y @f$ , where @f$ \alpha @f$  is a scalar
    * and @f$ x, y @f$  are vectors.
    */
   void
   axpy(
      const TYPE& alpha,
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& y,
      const bool interior_only = true);

   /**
    * Set each entry of this vector to absolute values of corresponding
    * entry in input vector.
    */
   void
   abs(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
      const bool interior_only = true);

   /**
    * Return the minimum data entry in this vector.  Note that this routine
    * returns a global min over all vector components and makes no adjustment
    * for coarser level vector data that may be masked out by the existence
    * of underlying fine values.  In particular, the control volumes are not
    * used in this operation.  This may change based on user needs.
    */
   TYPE
   min(
      const bool interior_only = true) const;

   /**
    * Return the maximum entry of this vector. Note that this routine
    * returns a global max over all vector components and makes no adjustment
    * for coarser level vector data that may be masked out by the existence
    * of underlying fine values.  In particular, the control volumes are not
    * used in this operation.  This may change based on user needs.
    */
   TYPE
   max(
      const bool interior_only = true) const;

   /**
    * Set data in this vector to random values.
    */
   void
   setRandomValues(
      const TYPE& width,
      const TYPE& low,
      const bool interior_only = true);

   /**
    * Return discrete @f$ L_1 @f$ -norm of this vector using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum @f$ \sum_i ( \| data_i \| cvol_i ) @f$ .  If the
    * control volume is not defined for a component, the contribution is
    * @f$ \sum_i ( \| data_i \| ) @f$  for that data component.  Thus, to have
    * a consistent norm calculation all components must have control
    * volumes, or no control volumes should be used at all.
    */
   double
   L1Norm(
      bool local_only = false) const;

   /**
    * Return discrete @f$ L_2 @f$ -norm of this vector using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum @f$ \sqrt{ \sum_i ( (data_i)^2 cvol_i ) } @f$ .
    * If the control volume is not defined for a component, the contribution
    * is @f$ \sqrt{ \sum_i ( (data_i)^2 ) } @f$  for that data component.
    * Thus, to have a consistent norm calculation all components must
    * have control volumes, or no control volumes should be used at all.
    */
   double
   L2Norm(
      bool local_only = false) const;

   /**
    * Return discrete weighted @f$ L_2 @f$ -norm of this vector using the control
    * volume to weight the contribution of the data and weight entries to
    * the sum.  That is, the return value is the sum @f$ \sqrt{ \sum_i (
    * (data_i * weight_i)^2 cvol_i ) } @f$ .  If the control volume is not defined
    * for a component, the contribution is
    * @f$ \sqrt{ \sum_i ( (data_i * weight_i)^2 ) } @f$  for that data component.
    * Thus, to have a consistent norm calculation all components must
    * have control volumes, or no control volumes should be used at all.
    */
   double
   weightedL2Norm(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& wgt) const;

   /**
    * Return discrete root mean squared norm of this vector.  If control
    * volumes are defined for all components, the return value is the
    * @f$ L_2 @f$ -norm divided by the square root of the sum of the control volumes.
    * If the control volume is not defined for a component, its contribution
    * to the norm corresponds to its @f$ L_2 @f$ -norm divided by the square root
    * of the number of data entries.  Thus, to have a consistent norm
    * calculation all components must have control volumes, or no
    * control volumes should be used at all.
    */
   double
   RMSNorm() const;

   /**
    * Return discrete weighted root mean squared norm of this vector.  If
    * control volumes are defined for all components, the return value is the
    * weighted @f$ L_2 @f$ -norm divided by the square root of the sum of the
    * control volumes.  If the control volume is not defined for a component,
    * its contribution to the norm corresponds to its weighted @f$ L_2 @f$ -norm
    * divided by the square root of the number of data entries.  Thus, to
    * have a consistent norm calculation all components must have
    * control volumes, or no control volumes should be used at all.
    */
   double
   weightedRMSNorm(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& wgt) const;

   /**
    * Return the @f$ \max @f$ -norm of this vector.  If control volumes are defined
    * for all components, the return value is the max norm over all data
    * values where the control volumes are non-zero.  If the control volume
    * is not defined for a component, its contribution to the norm will
    * take a max over all of its data values.  Thus, to have a consistent
    * norm calculation all components must have control volumes,
    * or no control volumes should be used at all.
    */
   double
   maxNorm(
      bool local_only = false) const;

   /**
    * Return the dot product of this vector with the argument vector.
    * If control volumes are defined for all components, the return value
    * is a weighted sum involving all data values where the control volumes
    * are non-zero.  If the control volume is not defined for a component,
    * its contribution to the sum will involve all of its data values.
    * Thus, to have a consistent dot product calculation all components must
    * have control volumes, or no control volumes should be used at all.
    */
   TYPE
   dot(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
      bool local_only = false) const;

   /**
    * Return 1 if @f$ \|x_i\| > 0 @f$  and @f$ w_i * x_i \leq 0 @f$ , for any @f$ i @f$  in
    * the set of vector data indices, where @f$ cvol_i > 0 @f$ .  Here, @f$ w_i @f$  is
    * a data entry in this vector.  Otherwise, return 0.  If the control
    * volume is undefined for a component, all data values for the component
    * are considered in the test.  Thus, to have a consistent test all
    * components must have control volumes, or no control volumes
    * should be used at all.
    */
   int
   computeConstrProdPos(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x) const;

   /**
    * Wherever @f$ cvol_i > 0 @f$  in the set of vector data indices, set @f$ w_i = 1 @f$
    * if @f$ \|x_i\| > \alpha @f$ , and @f$ w_i = 0 @f$  otherwise.  Here, @f$ w_i @f$  is a data
    * entry in this vector.  If the control volume is undefined for a
    * component, all data values for the component are involved in the
    * comparison.  Thus, to have a consistent comparison all components
    * must have control volumes, or no control volumes should be used at all.
    */
   void
   compareToScalar(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x,
      const TYPE& alpha);

   /**
    * Wherever @f$ cvol_i > 0 @f$  in the set of vector data indices, set
    * @f$ w_i = 1/x_i @f$  if @f$ x_i \neq 0 @f$ , and @f$ w_i = 0 @f$  otherwise.  Here, @f$ w_i @f$
    * is a data entry in this vector.  If the control volume is undefined for a
    * component, all data values for the component are involved in the
    * operation.  Thus, to have a consistent operation all components
    * must have control volumes, or no control volumes should be used at all.
    */
   int
   testReciprocal(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& x);

   /*!
    * @brief Compute max of "conditional" quotients of two arrays.
    *
    * Return the maximum of pointwise "conditional" quotients of the numerator
    * and denominator.
    *
    * The "conditional" quotient is defined as |numerator/denominator|
    * if the denominator is nonzero.  Otherwise, it is defined as
    * |numerator|.
    *
    * @b Note: This method is currently intended to support the
    * PETSc-2.1.6 vector wrapper only.  Please do not use it!
    */
   TYPE
   maxPointwiseDivide(
      const std::shared_ptr<SAMRAIVectorReal<TYPE> >& denom) const;

   /*!
    * @brief Get the length of this vector
    *
    * @param interior_only   If true, only count the elements in interior
    *                        of patches.
    * @return The length (number of elements in the underlying data)
    */
   int64_t getLength(const bool interior_only = true) const;

   //@}

private:
   /*
    * Static integer constant.
    */
   static const int DESCRIPTOR_ID_ARRAY_SCRATCH_SPACE = 10;

   // The following are not implemented
   SAMRAIVectorReal(
      const SAMRAIVectorReal&);
   SAMRAIVectorReal&
   operator = (
      const SAMRAIVectorReal&);

   /*
    * Private member function to set number of vector components.  This
    * is used during cloning.
    */
   void
   setNumberOfComponents(
      int num_comp);

   /*
    * Private function to set attributes for the specified vector component.
    * The component will be associated with the given variable and patch
    * data index information.  This function is specialized for
    * either double or float vector types.  This function is called from
    * addComponent() and clonevector().
    *
    * @pre comp_id < getNumberOfComponents()
    * @pre getPatchHierarchy()->getDim() == var->getDim()
    *
    * @post d_component_operations[comp_id]
    */
#ifdef _MSC_VER
   std::shared_ptr<math::HierarchyDataOpsReal<TYPE> > _bug_in_msvc;
#endif
   void
   setComponent(
      const int comp_id,
      const std::shared_ptr<hier::Variable>& var,
      const int data_id,
      const int control_vol_id = -1,
      const std::shared_ptr<math::HierarchyDataOpsReal<TYPE> >& vop =
         std::shared_ptr<math::HierarchyDataOpsReal<TYPE> >());

   static int s_instance_counter[SAMRAI::MAX_DIM_VAL];

   // shared data operations for variaous array-based types...
   static std::shared_ptr<math::HierarchyDataOpsReal<TYPE> >
   s_cell_ops[SAMRAI::MAX_DIM_VAL];
   static std::shared_ptr<math::HierarchyDataOpsReal<TYPE> >
   s_edge_ops[SAMRAI::MAX_DIM_VAL];
   static std::shared_ptr<math::HierarchyDataOpsReal<TYPE> >
   s_face_ops[SAMRAI::MAX_DIM_VAL];
   static std::shared_ptr<math::HierarchyDataOpsReal<TYPE> >
   s_node_ops[SAMRAI::MAX_DIM_VAL];
   static std::shared_ptr<math::HierarchyDataOpsReal<TYPE> >
   s_side_ops[SAMRAI::MAX_DIM_VAL];

   std::string d_vector_name;

   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;
   int d_coarsest_level;
   int d_finest_level;

   int d_number_components;

   // arrays for component information whose size is the number of components
   std::vector<std::shared_ptr<hier::Variable> > d_component_variable;
   std::vector<int> d_component_data_id;
   std::vector<std::shared_ptr<math::HierarchyDataOpsReal<TYPE> > >
   d_component_operations;
   std::vector<int> d_control_volume_data_id;

   // map from variable instance id to vector component index:
   // size = largest instance id over all variables in vector.
   std::vector<int> d_variableid_2_vectorcomponent_map;

   // output stream for vector data
   std::ostream* d_output_stream;

};

}
}

#include "SAMRAI/solv/SAMRAIVectorReal.cpp"

#endif
