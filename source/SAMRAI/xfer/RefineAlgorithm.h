/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Refine algorithm for data transfer between AMR levels
 *
 ************************************************************************/

#ifndef included_xfer_RefineAlgorithm
#define included_xfer_RefineAlgorithm

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/RefineClasses.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/hier/TimeInterpolateOperator.h"
#include "SAMRAI/xfer/VariableFillPattern.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"

#include <memory>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class RefineAlgorithm encapsulates the AMR communication
 * pattern to refine data to, copy data to, or fill physical boundary data on
 * any destination patch level.
 *
 * The basic procedure for moving data follows three steps:
 *
 * <ul>
 *    <li> interpolate data (spatial and possibly temporal) from coarser levels
 *    <li> copy data from the same level of refinement
 *    <li> fill physical boundary conditions regions
 * </ul>
 *
 * Each data communication procedure generally consists of three parts: an
 * algorithm, a schedule, and a patch strategy.  The algorithm describes
 * the patch data components and time and space interpolation operations,
 * but is independent of the configuration of the patches in an AMR hierarchy.
 * Patch data items and their associated spatial and time interpolation
 * operators are registered with an instantiation of this algorithm class.
 * To generate the communication dependencies for a particular patch hierarchy
 * configuration, the algorithm creates a refine schedule based on the state
 * of a given hierarchy and the information in the algorithm.  The schedule
 * can then perform the communication operations that move data to the
 * destination patch level using the associated operators.  User-defined
 * operations (such as filling physical boundaries and special interpolation
 * procedures) are provided through a refine patch strategy object.
 * User-defined interpolation operations can be written using the interfaces
 * in RefinePatchStrategy for preprocessRefine() and postProcessRefine().
 *
 * In general, source data is copied into the designated scratch data for
 * temporary processing.  The scratch space must contain sufficient ghost
 * cells to accommodate the stencil width of the given interpolation operators
 * and any physical boundary data that must be filled.  The scratch storage is
 * copied into the destination data space at the end of the communication
 * process.  Thus, copy operations between source, scratch, and destination
 * patch data objects must be defined.
 *
 * In general, the destination and scratch data components may be the same
 * (assuming that the scratch component has a sufficient ghost cells width).
 * The source and scratch components SHOULD NOT be the same, since the
 * interiors of the source space could be changed by the use of the scratch
 * data as temporary work space.
 *
 * It is the user's responsibility to register valid operations with the
 * refine algorithm so that the data communication can occur.  In particular,
 * communication operations (e.g., data refinement, data copy, etc.) are
 * performed in the order that items are registered for refinement with
 * a refine algorithm object.  Thus, order of registration must repect any
 * dependencies among patch data communicated.  Also, users who use
 * the preprocessRefine() and postProcessRefine() operations in the patch
 * strategy object must make sure that all data that is needed in those
 * operations are registered with the RefineAlgorithm using registerRefine()
 * whether or not the data is to be refined.
 *
 * Typical usage of a refine algorithm to perform inter-patch communication
 * on an AMR hierarchy involves four steps:
 *
 * <ul>
 *    <li> Construct a refine algorithm object.
 *    <li> Register refine operations with the refine algorithm.  Using the
 *         registerRefine() methods(s), one provides source and destination
 *         patch data information, as well as time and space interpolation
 *         operators as needed.  Two registerRefine() methods appear in this
 *         class; one supports time interpolation, one does not.
 *    <li> After all operations are registered with the algorithm, one
 *         creates a communication schedule using one of the createSchedule()
 *         methods.  These methods are distinguished by the resulting data
 *         communication pattern (e.g., interpatch communication on a single
 *         level, between two different levels at the same grid resolution,
 *         interpolation of data between different AMR hierarchy levels, etc.)
 *         Note that when creating a communication schedule, a concrete
 *         instance of a RefinePatchStrategy object may be required
 *         to supply physical boundary conditions as well as user-defined
 *         spatial data interpolation operations.
 *    <li> Invoke the fillData() method in the communication schedule to
 *         perform the data transfers.
 * </ul>
 *
 * Note that each refine schedule created by a refine algorithm remains valid
 * as long as the patches involved in the communication process do not change;
 * thus, they can be used for multiple data communication cycles.
 *
 * @see RefineSchedule
 * @see RefinePatchStrategy
 * @see RefineClasses
 */

class RefineAlgorithm
{
public:
   /*!
    * @brief Construct a refinement algorithm and initialize its basic state.
    *
    * Refinement operations must be registered with this algorithm
    * before it can do anything useful.  See the registerRefine() routines
    * for details
    */
   RefineAlgorithm();

   /*!
    * @brief The destructor releases all internal storage.
    */
   ~RefineAlgorithm();

   /*!
    * @brief Register a refine operation with the refine algorithm object.
    *
    * This method does not support time interpolation.  Data values will be
    * moved from the source data to the destination data using scratch
    * data as a temporary work space.  The scratch data must have
    * sufficient ghost cells to cover the required operator stencil width and
    * any needed physical boundary ghost cells.
    *
    * @param[in] dst   Patch data index filled on the destination level.
    * @param[in] src   Patch data index for source data.
    * @param[in] scratch   Patch data index for temporary work space data.
    * @param[in] oprefine  Refinement operator.  This may be a
    *                      null pointer.  In this case, refinement must be
    *                      handled by the refine patch strategy member
    *                      functions.  See the comments for
    *                      RefinePatchStrategy::preprocessRefine() and
    *                      RefinePatchStrategy::postprocessRefine().
    * @param[in] var_fill_pattern  std::shared_ptr to the variable fill
    *                              pattern, which can be used to restrict the
    *                              filling of data to a specific stencil.  If
    *                              the NULL default is used, then class
    *                              BoxGeometryVariableFillPattern will be used
    *                              internally.
    *
    * @pre !d_schedule_created
    */
   void
   registerRefine(
      const int dst,
      const int src,
      const int scratch,
      const std::shared_ptr<hier::RefineOperator>& oprefine,
      const std::shared_ptr<VariableFillPattern>& var_fill_pattern =
         std::shared_ptr<VariableFillPattern>(),
      const std::vector<int>& work_ids = std::vector<int>());

   /*!
    * @brief Register a refine operation with the refine algorithm object.
    *
    * This method supports time interpolation.  Time interpolation will take
    * place between the old and new source data components on coarser levels.
    * On the destination level, data will be moved from the source data to
    * the destination data using scratch data as a temporary work
    * space.  The scratch data must have sufficient ghost cells to cover
    * the required operator stencil width and any needed physical boundary
    * ghost cells.
    *
    * @param[in] dst   Patch data index filled on the destination level.
    * @param[in] src   Patch data index for source data.
    * @param[in] src_told  Patch data index for old data used in time
    *                      interpolation.
    * @param[in] src_tnew  Patch data index for new data used in time
    *                      interpolation.
    * @param[in] scratch   Patch data index for temporary work space data.
    * @param[in] oprefine  Refinement operator.  This may be a
    *                      null pointer.  In this case, refinement must be
    *                      handled by the refine patch strategy member
    *                      functions.  See the comments for
    *                      RefinePatchStrategy::preprocessRefine() and
    *                      RefinePatchStrategy::postprocessRefine().
    * @param[in] optime    Time interpolation operator.  This pointer may
    *                      not be null.
    * @param[in] var_fill_pattern  std::shared_ptr to the variable fill
    *                              pattern, which can be used to restrict the
    *                              filling of data to a specific stencil.  If
    *                              the NULL default is used, then class
    *                              BoxGeometryVariableFillPattern will be used
    *                              internally.
    *
    * @pre optime
    * @pre !d_schedule_created
    */
   void
   registerRefine(
      const int dst,
      const int src,
      const int src_told,
      const int src_tnew,
      const int scratch,
      const std::shared_ptr<hier::RefineOperator>& oprefine,
      const std::shared_ptr<hier::TimeInterpolateOperator>& optime,
      const std::shared_ptr<VariableFillPattern>& var_fill_pattern =
         std::shared_ptr<VariableFillPattern>(),
      const std::vector<int>& work_ids = std::vector<int>());

   /*!
    * @brief Create a communication schedule for communicating data within a
    * single level.
    *
    * The schedule will communicate data from the interiors of the source
    * data into the interior and ghosts of the destination data on the
    * the same level where those sources and destinations overlap.
    *
    * Neither time nor spatial interpolation is performed.
    *
    * Note that the schedule remains valid as long as the level does not
    * change; thus, it can be used for multiple data communication cycles
    * cycles.
    *
    * @return std::shared_ptr to refine schedule that performs the data
    *         transfers.
    *
    * @param[in] level           Level on which communication occurs.  This
    *                            pointer cannot be null.
    * @param[in] patch_strategy  Optional pointer to a refine patch strategy
    *                            that provides user-defined physical boundary
    *                            filling operations.  If this patch strategy is
    *                            null (default state), then no physical
    *                            boundary filling is performed.
    * @param[in] transaction_factory Optional std::shared_ptr to a refine
    *                                transaction factory that creates data
    *                                transactions for the schedule.  If this
    *                                pointer is null (default state), then a
    *                                StandardRefineTransactionFactory object
    *                                will be used.
    *
    * @pre level
    */
   std::shared_ptr<RefineSchedule>
   createSchedule(
      const std::shared_ptr<hier::PatchLevel>& level,
      RefinePatchStrategy* patch_strategy = 0,
      const std::shared_ptr<RefineTransactionFactory>& transaction_factory =
         std::shared_ptr<RefineTransactionFactory>());

   /*
    * @brief Same as the above, except with fill_pattern specified.
    *
    * @param[in] fill_pattern  Indicates which parts of the destination level
    *                          to fill.  See RefineSchedule for available
    *                          patterns.
    *
    * @pre level
    */
   std::shared_ptr<RefineSchedule>
   createSchedule(
      const std::shared_ptr<PatchLevelFillPattern>& fill_pattern,
      const std::shared_ptr<hier::PatchLevel>& level,
      RefinePatchStrategy* patch_strategy = 0,
      const std::shared_ptr<RefineTransactionFactory>& transaction_factory =
         std::shared_ptr<RefineTransactionFactory>());

   /*!
    * @brief Create a communication schedule that communicates data between
    * two levels that are at the same level of resolution.
    *
    * Data will be communicated from the interiors of the source data on the
    * source level into the interior and ghosts of the destination data on a
    * destination level where those sources and destinations overlap.
    *
    * Note that both levels must reside in the same AMR hierarchy index space,
    * or in index spaces that represent the same level of mesh refinement.
    * No spatial interpolation is performed.
    *
    * In certain rare cases it may be desired to use this schedule to perform
    * time interpolation, in which case the use_time_interpolation optional
    * argument should be set to true.
    *
    * Note that the schedule remains valid as long as the levels do not
    * change; thus, it can be used for multiple data communication cycles.
    *
    * @return std::shared_ptr to refine schedule that performs the data
    *         transfers.
    *
    * @param[in] dst_level       std::shared_ptr to destination level; cannot
    *                            be null.
    * @param[in] src_level       std::shared_ptr to source level; cannot be
    *                            null.
    * @param[in] patch_strategy  std::shared_ptr to a refine patch strategy
    *                            that provides user-defined physical boundary
    *                            filling operations.  If this patch strategy is
    *                            null (default state), then no physical
    *                            boundary filling is performed.
    * @param[in]  use_time_interpolation  Flag to create the schedule with
    *                                     the ability to perform time
    *                                     interpolation.
    * @param[in] transaction_factory  std::shared_ptr to a refine transaction
    *                                 factory that creates data transactions
    *                                 for the schedule.  If this pointer is
    *                                 null (default state), then a
    *                                 StandardRefineTransactionFactory object
    *                                 will be used.
    *
    * @pre dst_level
    * @pre src_level
    * @pre dst_level->getDim() == src_level->getDim()
    */
   std::shared_ptr<RefineSchedule>
   createSchedule(
      const std::shared_ptr<hier::PatchLevel>& dst_level,
      const std::shared_ptr<hier::PatchLevel>& src_level,
      RefinePatchStrategy* patch_strategy = 0,
      bool use_time_interpolation = false,
      const std::shared_ptr<RefineTransactionFactory>& transaction_factory =
         std::shared_ptr<RefineTransactionFactory>());

   /*!
    * @brief Same as the above, except with fill_pattern specified.
    *
    * @param[in] fill_pattern  Indicates which parts of the destination level
    *                          to fill.  See RefineSchedule for available
    *                          patterns.
    * @param[in] dst_level
    * @param[in] src_level
    * @param patch_strategy
    * @param[in] use_time_interpolation
    * @param[in] transaction_factory
    *
    * @pre dst_level
    * @pre src_level
    * @pre dst_level->getDim() == src_level->getDim()
    */
   std::shared_ptr<RefineSchedule>
   createSchedule(
      const std::shared_ptr<PatchLevelFillPattern>& fill_pattern,
      const std::shared_ptr<hier::PatchLevel>& dst_level,
      const std::shared_ptr<hier::PatchLevel>& src_level,
      RefinePatchStrategy* patch_strategy = 0,
      bool use_time_interpolation = false,
      const std::shared_ptr<RefineTransactionFactory>& transaction_factory =
         std::shared_ptr<RefineTransactionFactory>());

   /*!
    * @brief Create a communication schedule that communicates data within a
    * single level and interpolates data from coarser hierarchy levels
    * where needed.
    *
    * Data will be communicated from the interiors of the source data on
    * on the given level to the interiors and ghosts of destination data on
    * the same level where those sources and destinations overlap.  Where they
    * do not overlap, data will be interpolated from source data on coarser
    * levels in the patch hierarchy.
    *
    * Data is time interpolated between old and new sources on coarser levels
    * when and where time interpolation is needed and copied from the source
    * components on the patch level into the destination components otherwise.
    *
    * In certain rare cases in may be necessary to perform time interpolation
    * between old and new sources on the given patch level.  In this case
    * the optional argument use_time_interpolation should be set to true.
    * Regardless of the value of this argument, time interpolation on
    * coarser levels will always occur whenever needed.
    *
    * Note that the next coarser level number must correspond to a level
    * in the hierarchy that represents a region of coarser index space than
    * the destination level.
    *
    * Note that the schedule remains valid as long as the levels involved
    * in its creation do not change; thus, it can be used for multiple
    * data communication cycles.
    *
    * @return std::shared_ptr to refine schedule that performs the data
    *         transfers.
    *
    * @param[in] level         std::shared_ptr to destination level; cannot
    *                          be null.
    * @param[in] next_coarser_level  Level number of next coarser patch level
    *                                in the patch hierarchy relative to the
    *                                destination level.  Note that when the
    *                                destination level has number zero (i.e.,
    *                                the coarsest level), this value should
    *                                value should be < 0.
    * @param[in] hierarchy     std::shared_ptr to patch hierarchy from which
    *                          data to fill level should come.  This pointer
    *                          may be null only when the next_coarser_level is
    *                          < 0.
    * @param[in] patch_strategy  std::shared_ptr to a refine patch strategy
    *                            that provides user-defined physical boundary
    *                            filling operations and user-defined spatial
    *                            interpolation operations.  If this patch
    *                            strategy is null (default state), then no
    *                            physical boundary filling or user-defined
    *                            interpolation is performed.  Note that this
    *                            may cause problems if the interpolation
    *                            stencils require physical boundary data on
    *                            the coarser levels.
    * @param[in] use_time_interpolation  Boolean flag to create the schedule
    *                                    the ability to perform time
    *                                    interpolation on the destination
    *                                    level.  Default is no time
    *                                    interpolation (false).
    * @param[in] transaction_factory  std::shared_ptr to a refine transaction
    *                                 factory that creates data transactions
    *                                 for the schedule.  If this pointer is
    *                                 null (default state), then a
    *                                 StandardRefineTransactionFactory object
    *                                 will be used.
    *
    * @pre level
    * @pre (next_coarser_level == -1) || hierarchy
    * @pre !hierarchy || (level->getDim() == hierarchy->getDim())
    */
   std::shared_ptr<RefineSchedule>
   createSchedule(
      const std::shared_ptr<hier::PatchLevel>& level,
      const int next_coarser_level,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      RefinePatchStrategy* patch_strategy = 0,
      bool use_time_interpolation = false,
      const std::shared_ptr<RefineTransactionFactory>& transaction_factory =
         std::shared_ptr<RefineTransactionFactory>());

   /*!
    * @brief Same as the above, except with fill_pattern specified.
    *
    * @param[in] fill_pattern  Indicates which parts of the destination level
    *                          to fill.  See RefineSchedule for available
    *                          patterns.
    * @param[in] level
    * @param[in] next_coarser_level
    * @param[in] hierarchy
    * @param patch_strategy
    * @param[in] use_time_interpolation
    * @param[in] transaction_factory
    *
    * @pre level
    * @pre (next_coarser_level == -1) || hierarchy
    * @pre !hierarchy || (level->getDim() == hierarchy->getDim())
    */
   std::shared_ptr<RefineSchedule>
   createSchedule(
      const std::shared_ptr<PatchLevelFillPattern>& fill_pattern,
      const std::shared_ptr<hier::PatchLevel>& level,
      const int next_coarser_level,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      RefinePatchStrategy* patch_strategy = 0,
      bool use_time_interpolation = false,
      const std::shared_ptr<RefineTransactionFactory>& transaction_factory =
         std::shared_ptr<RefineTransactionFactory>());

   /*!
    * @brief Create a communication schedule that communicates data from a
    * source level to a destination level and interpolates data from coarser
    * hierarchy levels where needed.
    *
    * Data will be communicated from the interiors of the source data on
    * on the source level to the interiors and ghosts of destination data on
    * the destination level where those sources and destinations overlap.
    * Where they do not overlap, data will be interpolated from source data on
    * coarser levels in the patch hierarchy.
    *
    * Data is time interpolated between old and new sources on coarser levels
    * when and where time interpolation is needed and copied from the source
    * components on the source level into the destination components otherwise.
    *
    * This form of schedule construction is typically used after regridding
    * (where the source level is the patch level being replaced by the
    * destination level in the patch hierarchy) or when the data on destination
    * patch level is to be overwritten by data interpolated from coarser levels
    * in the patch hierarchy.  In the first case, data on the destination level
    * will be copied from the source level in regions where those two levels
    * overlap and filled with interpolated values from the hierarchy elsewhere.
    * In the latter case, the source level pointer may be null.  Then, data on
    * the destination level will be filled using interpolated data from coarser
    * hierarchy levels.
    *
    * In certain rare cases in may be desired to perform time interpolation
    * between old and new sources onto the destination level.  In this case
    * the optional argument use_time_interpolation should be set to true.
    * Regardless of the value of this argument, time interpolation on
    * coarser levels will always occur whenever needed.
    *
    * Note that when the source level pointer is non-null, the index spaces of
    * the source and destination levels must be aligned with one another.
    *
    * Note that the schedule remains valid as long as the levels involved
    * in its creation do not change; thus, it can be used for multiple
    * data communication cycles.
    *
    * @return std::shared_ptr to refine schedule that performs the data
    *         transfers.
    *
    * @param[in] dst_level      std::shared_ptr to destination level; cannot
    *                           be null.
    * @param[in] src_level      std::shared_ptr to source level. This pointer
    *                           may be null.  In this case, data on the
    *                           destination level will be filled only using
    *                           interpolated data from coarser hierarchy
    *                           levels.  When this pointer is not null, the
    *                           source level must live in the same AMR
    *                           hierarchy index space as the destination level.
    * @param[in] next_coarser_level  Level number of next coarser patch level
    *                                in the patch hierarchy relative to the
    *                                destination level.  Note that when the
    *                                destination level has number zero (i.e.,
    *                                the coarsest level), this value should
    *                                value should be < 0.
    * @param[in] hierarchy    std::shared_ptr to patch hierarchy from which
    *                         data to fill level should come.  This pointer may
    *                         be null only when the next_coarser_level is < 0.
    * @param[in] patch_strategy  std::shared_ptr to a refine patch strategy
    *                            that provides user-defined physical boundary
    *                            filling operations and user-defined spatial
    *                            interpolation operations.  If this patch
    *                            strategy is null (default state), then no
    *                            physical boundary filling or user-defined
    *                            interpolation is performed.  Note that this
    *                            may cause problems if the interpolation
    *                            stencils require physical boundary data on
    *                            the coarser levels.
    * @param[in] use_time_interpolation  Boolean flag to create the schedule
    *                                    the ability to perform time
    *                                    interpolation on the destination
    *                                    level.  Default is no time
    *                                    interpolation (false).
    * @param[in] transaction_factory  std::shared_ptr to a refine transaction
    *                                 factory that creates data transactions
    *                                 for the schedule.  If this pointer is
    *                                 null (default state), then a
    *                                 StandardRefineTransactionFactory object
    *                                 will be used.
    *
    * @pre dst_level
    * @pre (next_coarser_level == -1) || hierarchy
    * @pre !src_level || (dst_level->getDim() == src_level->getDim())
    * @pre !hierarchy || (dst_level->getDim() == hierarchy->getDim())
    */
   std::shared_ptr<RefineSchedule>
   createSchedule(
      const std::shared_ptr<hier::PatchLevel>& dst_level,
      const std::shared_ptr<hier::PatchLevel>& src_level,
      const int next_coarser_level,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      RefinePatchStrategy* patch_strategy = 0,
      bool use_time_interpolation = false,
      const std::shared_ptr<RefineTransactionFactory>& transaction_factory =
         std::shared_ptr<RefineTransactionFactory>());

   /*!
    * @brief Same as the above, except with fill_pattern specified.
    *
    * @param[in] fill_pattern  Indicates which parts of the destination level
    *                          to fill.  See RefineSchedule for available
    *                          patterns.
    * @param[in] dst_level
    * @param[in] src_level
    * @param[in] next_coarser_level
    * @param[in] hierarchy
    * @param patch_strategy
    * @param[in] use_time_interpolation
    * @param[in] transaction_factory
    *
    * @pre dst_level
    * @pre (next_coarser_level == -1) || hierarchy
    * @pre !src_level || (dst_level->getDim() == src_level->getDim())
    * @pre !hierarchy || (dst_level->getDim() == hierarchy->getDim())
    */
   std::shared_ptr<RefineSchedule>
   createSchedule(
      const std::shared_ptr<PatchLevelFillPattern>& fill_pattern,
      const std::shared_ptr<hier::PatchLevel>& dst_level,
      const std::shared_ptr<hier::PatchLevel>& src_level,
      const int next_coarser_level,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      RefinePatchStrategy* patch_strategy = 0,
      bool use_time_interpolation = false,
      const std::shared_ptr<RefineTransactionFactory>& transaction_factory =
         std::shared_ptr<RefineTransactionFactory>());

   /*!
    * @brief Given a previously-generated refine schedule, check for
    * consistency with this refine algorithm object to see whether a call to
    * resetSchedule() is a valid operation.
    *
    * Consistency means that the number of operations registered must be the
    * same and the source and destination patch data items and operators must
    * have identical characteristics (i.e., data centering, ghost cell widths,
    * stencil requirements, etc.).  However, the specific source, destination
    * patch data ids and refine operators can be different.  The specific
    * time interpolation operators and variable fill patterns must be
    * the same.  See RefineClasses::classesMatch() for more details.
    *
    * @return true if schedule reset is valid; false otherwise.
    *
    * @param[in] schedule  std::shared_ptr to refine schedule, which cannot
    *                      be null.
    *
    * @pre schedule
    */
   bool
   checkConsistency(
      const std::shared_ptr<RefineSchedule>& schedule) const;

   /*!
    * @brief Given a previously-generated refine schedule, reconfigure it to
    * peform the communication operations registered with THIS refine
    * algorithm object.
    *
    * That is, the schedule will be transformed so that it will function as
    * though this refine algorithm created it.  Note that the set of
    * operations registered with this refine algorithm must be essentially
    * the same as those registered with the refine algorithm that created the
    * schedule originally, and this is enforced using a call to
    * checkConsistency().  An error will result if the schedule is not
    * consistent with this RefineAlgorithm object according to the criteria
    * in checkConsistency().
    *
    * @param[in,out] schedule  std::shared_ptr to refine schedule, which
    *                          cannot be null.
    *
    * @pre schedule
    * @pre d_refine_classes->classesMatch(schedule->getEquivalenceClasses())
    */
   void
   resetSchedule(
      const std::shared_ptr<RefineSchedule>& schedule) const;

   /*!
    * @brief Return the refine equivalence classes used in the algorithm.
    */
   const std::shared_ptr<RefineClasses>&
   getEquivalenceClasses() const
   {
      return d_refine_classes;
   }

   /*!
    * @brief Set the pointer to the refine equivalence classes to be equal
    * to the given argument.
    *
    * @param[in] refine_classes A pointer to refine equivalence classes
    */
   void
   setEquivalenceClasses(
      const std::shared_ptr<RefineClasses>& refine_classes)
   {
      d_refine_classes = refine_classes;
   }

   /*!
    * @brief Print the refine algorithm state to the specified data stream.
    *
    * @param[out] stream Output data stream.
    */
   void
   printClassData(
      std::ostream& stream) const;

private:
   RefineAlgorithm(
      const RefineAlgorithm&);                  // not implemented
   RefineAlgorithm&
   operator = (
      const RefineAlgorithm&);                  // not implemented

   /*!
    * RefineClasses object holds all of the registered refine items
    */
   std::shared_ptr<RefineClasses> d_refine_classes;

   /*!
    * Tells if any schedule has yet been created using this object.
    */
   bool d_schedule_created;

};

}
}

#endif
