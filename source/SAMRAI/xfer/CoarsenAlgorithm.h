/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Coarsening algorithm for data transfer between AMR levels
 *
 ************************************************************************/

#ifndef included_xfer_CoarsenAlgorithm
#define included_xfer_CoarsenAlgorithm

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/CoarsenClasses.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/CoarsenPatchStrategy.h"
#include "SAMRAI/xfer/CoarsenTransactionFactory.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iostream>
#include <memory>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class CoarsenAlgorithm encapsulates the AMR communication
 * pattern to coarsen data from a finer level to any coarser level.
 *
 * Most often, data is coarsened from the interiors of source patch
 * components on the source patch level into interiors of destination
 * patch components on the destination level.  If the coarsening operators
 * require ghost cells on a source component, then sufficient ghost cell
 * storage must be provided by the source patch data component, and those
 * ghost cells must be filled before calling the data coarsening routines.
 *
 * Communication algorithms generally consist of three parts: an algorithm,
 * a schedule, and a patch strategy.  The algorithm describes the communication
 * between patch data items but is independent of the configuration of the
 * AMR hierarchy.  PatchData items and their associated coarsening
 * operators are registered with the algorithm.  To generate the communication
 * dependencies for a particular hierarchy configuration, the algorithm
 * generates a schedule based on the current hierarchy configuration.  This
 * schedule then performs the communication based on the registered data types
 * and their associated operators.  User-defined coarsening operations can
 * be written using the interfaces in CoarsenPatchStrategy for
 * preprocessCoarsen() and postProcessCoarsen().
 *
 * The source patch data space is used during processing to store temporary
 * data.  Thus, the user-defined coarsening operators should operate on the
 * source space by using the patch data with those indices.
 *
 * It is the user's responsibility to register valid operations with the
 * coarsen algorithm so that the data communication can occur.  In particular,
 * communication operations (e.g., data coarsening, data copy, etc.) are
 * performed in the order that items are registered for coarsening with
 * a coarsen algorithm object.  Thus, order of registration must repect any
 * dependencies among patch data communicated.  Also, users who use
 * the preprocessCoarsen() and postProcessCoarsen() operations in the patch
 * strategy object must make sure that all data that is needed in those
 * operations are registered with the CoarsenAlgorithm using registerCoarsen()
 * whether or not the data is to be coarsened.
 *
 * Typical usage of a coarsen algorithm to perform data coarsening
 * on an AMR hierarchy involves four steps:
 *
 * <ul>
 *    <li> Construct a coarsen algorithm object.
 *    <li> Register coarsen operations with the coarsen algorithm.  Using the
 *         registerCoarsen() methods(s), one provides source and destination
 *         patch data information, as well as spatial coarsening operators
 *         as needed.
 *    <li> After all operations are registered with the algorithm, one
 *         creates a communication schedule using the createSchedule()
 *         method.  This method identifies the source (fine) and destination
 *         (coarse) patch levels for data coarsening.  Note that when creating
 *         a communication schedule, a concrete instance of a
 *         CoarsenPatchStrategy object may be required to supply
 *         user-defined spatial data coarsening operations.
 *    <li> Invoke the coarsenData() method in the communication schedule to
 *         perform the data transfers.
 * </ul>
 *
 * Note that each coarsen schedule created by a coarsen algorithm remains valid
 * as long as the levels involved in the communication process do not change;
 * thus, they can be used for multiple data communication cycles.
 *
 * @see CoarsenSchedule
 * @see CoarsenPatchStrategy
 * @see CoarsenClasses
 */

class CoarsenAlgorithm
{
public:
   /*!
    * @brief Construct a coarsening algorithm and initialize its basic state.
    *
    * Coarsening operations must be registered with this algorithm
    * before it can do anything useful.  See the registerCoarsen() routine
    * for details
    *
    * @param[in] dim  Dimension
    * @param[in] fill_coarse_data  boolean flag indicating whether pre-existing
    *                              coarse level data is needed for the data
    *                              coarsening operations.  By default this
    *                              argument is false.  If a true value is
    *                              value is provided, then source data will be
    *                              filled on a temporary coarse patch level
    *                              (copied from the actual coarse level source
    *                              data) for use in coarsening operations
    *                              registered with this algorithm.  This option
    *                              should only be used by those who
    *                              specifically require this behavior and who
    *                              know how to properly process the patch
    *                              data on coarse and fine patch levels during
    *                              the coarsening process.
    */
   explicit CoarsenAlgorithm(
      const tbox::Dimension& dim,
      bool fill_coarse_data = false);

   /*!
    * @brief The destructor releases all internal storage.
    */
   ~CoarsenAlgorithm();

   /*!
    * @brief Register a coarsening operation with the coarsening algorithm.
    *
    * Data from the interiors of the source component on a source (fine) patch
    * level will be coarsened into the source component of a temporary (coarse)
    * patch level and then copied into the destination component on the
    * destination (coarse) patch level.  If the coarsening operator requires
    * data in ghost cells outside of the patch interiors (i.e., a non-zero
    * stencil width), then those ghost cells must exist in the source patch
    * data component and the ghost cells must be filled with valid data on the
    * source level before a call to invoke the communication schedule.  Note
    * that the source and destination components may be the same.
    *
    * Some special circumstances require that data be coarsened from the
    * ghost cell regions of a finer level and the resulting coarsened data
    * should be copied to the destination patch level.  When this is the case,
    * the optional integer vector argument should be set to the cell width, in
    * the destination (coarser) level index space, of the region around the
    * fine level where this coarsening should occur.  For example, if the
    * coarser level needs data in a region two (coarse) cells wide around the
    * boundary of the finer level, then the gcw_to_coarsen should be set to a
    * vector with all entries set to two.  Moreover, if the ratio
    * between coarse and fine levels is four in this case, then the source
    * patch data is required to have at least eight ghost cells.
    *
    * @param[in] dst       Patch data index filled on destination level.
    * @param[in] src       Patch data index coarsened from the source level.
    * @param[in] opcoarsen Coarsening operator.  This may be a null pointer.
    *                  If null, coarsening must be handled by the coarsen
    *                  patch strategy member functions.  See the comments for
    *                  CoarsenPatchStrategy::preprocessCoarsen() and
    *                  CoarsenPatchStrategy::postprocessCoarsen().
    * @param[in] gcw_to_coarsen Ghost cell width to be used when data should
    *                           be coarsened from ghost cell regions of the
    *                           source (finer) level into the destination
    *                           (coarser) level. If coarsening from fine ghost
    *                           cell regions is not desired, then it should be
    *                           a zero IntVector.  If this argument is nonzero,
    *                           its value should be the cell width, in the
    *                           destination (coarser) index space, of the
    *                           region around the fine level where this
    *                           coarsening should occur.  This argument should
    *                           only be made nonzero by those who specifically
    *                           require this special behavior and know how to
    *                           properly process the patch data on coarse and
    *                           fine patch levels during the coarsening
    *                           process.
    * @param[in] var_fill_pattern std::shared_ptr to the variable fill
    *                             pattern, which controls how box overlaps are
    *                             constructed.  If the NULL default is used,
    *                             then class BoxGeometryVariableFillPattern
    *                             will be used internally.
    *
    * @pre !d_schedule_created
    */
   void
   registerCoarsen(
      const int dst,
      const int src,
      const std::shared_ptr<hier::CoarsenOperator>& opcoarsen,
      const hier::IntVector& gcw_to_coarsen,
      const std::shared_ptr<VariableFillPattern>& var_fill_pattern =
         std::shared_ptr<VariableFillPattern>());

   /*!
    * @brief Register a coarsening operation with the coarsening algorithm.
    *
    * This will do all of the same things as the above registerCoarsen(),
    * except it does not have the gcw_to_coarsen parameter.
    */
   void
   registerCoarsen(
      const int dst,
      const int src,
      const std::shared_ptr<hier::CoarsenOperator>& opcoarsen,
      const std::shared_ptr<VariableFillPattern>& var_fill_pattern =
         std::shared_ptr<VariableFillPattern>())
   {
      registerCoarsen(dst, src, opcoarsen,
         hier::IntVector::getZero(d_dim), var_fill_pattern);
   }

   /*!
    * @brief Create a communication schedule to coarsen data from the given
    * fine patch level to the given coarse patch level.
    *
    * This communication schedule may then be executed to perform
    * the data transfers.  This schedule creation procedure assumes that
    * the coarse level represents a region of coarser index space than the
    * fine level.  To avoid potentially erroneous behavior, the coarse level
    * domain should cover the domain of the fine level.
    *
    * Note that the schedule remains valid as long as the levels do not
    * change; thus, it can be used for multiple data communication cycles.
    *
    * @return std::shared_ptr to coarsen schedule that performs the data
    *         transfers.
    *
    * @param[in] crse_level     std::shared_ptr to coarse (destination) level.
    * @param[in] fine_level     std::shared_ptr to fine (source) level.
    * @param[in] coarsen_strategy std::shared_ptr to a coarsen patch strategy
    *                           that provides user-defined coarsen operations.
    *                           If this patch strategy is null (default state),
    *                           then no user-defined coarsen operations will be
    *                           performed.
    * @param[in] transaction_factory Optional std::shared_ptr to a coarsen
    *                                transaction factory that creates data
    *                                transactions for the schedule.  If this
    *                                pointer is null default state), then a
    *                                StandardCoarsenTransactionFactory object
    *                                will be used.
    *
    * @pre crse_level && fine_level
    * @pre (getDim() == crse_level->getDim()) &&
    *      (getDim() == fine_level->getDim())
    */
   std::shared_ptr<CoarsenSchedule>
   createSchedule(
      const std::shared_ptr<hier::PatchLevel>& crse_level,
      const std::shared_ptr<hier::PatchLevel>& fine_level,
      CoarsenPatchStrategy* coarsen_strategy = 0,
      const std::shared_ptr<CoarsenTransactionFactory>& transaction_factory =
         std::shared_ptr<CoarsenTransactionFactory>());

   /*!
    * @brief Given a previously-generated coarsen schedule, check for
    * consistency with this coarsen algorithm object to see whether a call to
    * resetSchedule is a valid operation.
    *
    * Consistency means that the number of operations registered must be the
    * same and the source and destination patch data items and operators must
    * have identical characteristics (i.e., data centering, ghost cell widths,
    * stencil requirements, etc.).  However, the specific source, destination
    * patch data ids and operators can be different.
    *
    * @return true if schedule reset is valid; false otherwise.
    *
    * @param[in] schedule  std::shared_ptr to coarsen schedule, which cannot
    *                      be null.
    *
    * @pre schedule
    */
   bool
   checkConsistency(
      const std::shared_ptr<CoarsenSchedule>& schedule) const
   {
      TBOX_ASSERT(schedule);
      return d_coarsen_classes->classesMatch(schedule->getEquivalenceClasses());
   }

   /*!
    * @brief Given a previously-generated coarsen schedule, reconfigure it to
    * peform the communication operations registered with this coarsen
    * algorithm object.
    *
    * That is, the schedule will be transformed so that it will function as
    * though this coarsen algorithm created it.  Note that the set of
    * operations registered with this coarsen algorithm must be essentially
    * the same as those registered with the coarsen algorithm that created the
    * schedule originally, and this is enforced using a call to
    * checkConsistency().
    *
    * @param[in,out] schedule  std::shared_ptr to coarsen schedule, which
    *                          cannot be null.
    *
    * @pre schedule
    * @pre d_coarsen_classes->classesMatch(schedule->getEquivalenceClasses())
    */
   void
   resetSchedule(
      const std::shared_ptr<CoarsenSchedule>& schedule) const;

   /*!
    * @brief Print the coarsen algorithm state to the specified data stream.
    *
    * @param[out] stream Output data stream.
    */
   void
   printClassData(
      std::ostream& stream) const;

   /*!
    * @brief Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const
   {
      return d_dim;
   }

private:
   CoarsenAlgorithm(
      const CoarsenAlgorithm&);               // not implemented
   CoarsenAlgorithm&
   operator = (
      const CoarsenAlgorithm&);                   // not implemented

   /*!
    * @brief Dimension of the object.
    */
   const tbox::Dimension d_dim;

   /*!
    * CoarsenClasses object holds all of the registered coarsen items.
    */
   std::shared_ptr<CoarsenClasses> d_coarsen_classes;

   /*!
    * Tells if special behavior to pre-fill the temporary coarse level with
    * existing coarse data values is turned on.
    */
   bool d_fill_coarse_data;

   /*!
    * Tells if any schedule has yet been created using this object.
    */
   bool d_schedule_created;
};

}
}

#endif
