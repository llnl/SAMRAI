/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Coarsening schedule for data transfer between AMR levels
 *
 ************************************************************************/

#ifndef included_xfer_CoarsenSchedule
#define included_xfer_CoarsenSchedule

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Schedule.h"
#include "SAMRAI/tbox/ScheduleOpsStrategy.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/xfer/CoarsenClasses.h"
#include "SAMRAI/xfer/CoarsenPatchStrategy.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/CoarsenTransactionFactory.h"

#include <iostream>
#include <memory>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class CoarsenSchedule performs the communication operations
 * to coarsen data from a finer level to a coarser level.
 *
 * Data is typically coarsened from the interiors of source patch
 * components on the source patch level into interiors of destination
 * patch components on the destination level.  Variations are possible
 * for special situations; see the CoarsenAlgorithm class header for
 * more information.  Generally, the source patch data must contain
 * sufficient ghost cells to satisfy the coarsening operators
 * involved.  If a coarsen operator has a non-zero ghost cell width,
 * then the source ghost cells must be filled before the coarsen
 * schedule is executed.  The communication schedule is executed by
 * calling member function coarsenData().
 *
 * Each schedule object is typically created by a coarsen algorithm
 * and represents communication dependencies for a particular
 * configuration of the AMR hierarchy.  The communication schedule is
 * only valid for that particular configuration and must be
 * regenerated when the AMR patch hierarchy changes.  As long as the
 * patch levels involved in the creation of the schedule remain
 * unchanged, the schedule may be used for multiple communication
 * cycles.  For more information about creating refine schedules, see
 * the CoarsenAlgorithm header file.
 *
 * NOTE: Algorithmic variations are available by calling the static method
 *       CoarsenSchedule::setScheduleGenerationMethod(), which
 *       sets the option for all instances of the class.
 *
 * @see CoarsenAlgorithm
 * @see CoarsenPatchStrategy
 * @see CoarsenClasses
 */

class CoarsenSchedule
{
public:
   /*!
    * @brief Constructor
    *
    * Creates a coarsen schedule that coarsens data from source patch data
    * components on the fine level into the destination patch data components
    * on the coarse level.
    *
    * In general, this constructor is called by a CoarsenAlgorithm object.
    * For possible variations on data coarsening, see the CoarsenAlgorithm
    * class header information.
    *
    * If the coarsening operators require data from ghost cells, then the
    * associated source patch data components must have a sufficient ghost
    * cell width and and they must be filled with valid data before calling
    * coarsenData().
    *
    * @param[in] crse_level      std::shared_ptr to coarse (destination)
    *                            patch level.
    * @param[in] fine_level      std::shared_ptr to fine (source) patch level.
    * @param[in] coarsen_classes std::shared_ptr to structure containing
    *                            patch data and operator information.  In
    *                            general, this is constructed by the calling
    *                            CoarsenAlgorithm object.
    * @param[in] transaction_factory  std::shared_ptr to a factory object
    *                                 that will create data transactions.
    * @param[in] patch_strategy  std::shared_ptr to a coarsen patch strategy
    *                            object that provides user-defined coarsen
    *                            operations.  This pointer may be null, in
    *                            which case no user-defined coarsen operations
    *                            will be performed.
    * @param[in] fill_coarse_data  Boolean indicating whether coarse data
    *                              should be filled before coarsening
    *                              operations are done.
    *
    * @pre crse_level
    * @pre fine_level
    * @pre coarsen_classes
    * @pre transaction_factory
    * @pre crse_level->getDim() == fine_level->getDim()
    */
   CoarsenSchedule(
      const std::shared_ptr<hier::PatchLevel>& crse_level,
      const std::shared_ptr<hier::PatchLevel>& fine_level,
      const std::shared_ptr<CoarsenClasses>& coarsen_classes,
      const std::shared_ptr<CoarsenTransactionFactory>& transaction_factory,
      CoarsenPatchStrategy* patch_strategy,
      bool fill_coarse_data);

   /*!
    * @brief The destructor for the schedule releases all internal
    * storage.
    */
   ~CoarsenSchedule();

   /*!
    * @brief Read static data from input database.
    */
   void
   getFromInput();

   /*!
    * @brief Reset this coarsen schedule to perform data transfers asssociated
    * items in the given CoarsenClasses argument.
    *
    * The schedule will be changed to operate on data given by the
    * coarsen_classes argument rather than the data it has previously been
    * set to operate on.
    *
    * @param[in] coarsen_classes  std::shared_ptr to structure containing
    *                             patch data and operator information.  In
    *                             general, this is constructed by the calling
    *                             CoarsenAlgorithm object.  This pointer must
    *                             be non-null.
    *
    * @pre coarsen_classes
    */
   void
   reset(
      const std::shared_ptr<CoarsenClasses>& coarsen_classes);

   /*!
    * @brief Execute the stored communication schedule and perform the data
    * movement.
    */
   void
   coarsenData() const;

   /*!
    * @brief Return the coarsen equivalence classes used in the schedule.
    */
   const std::shared_ptr<CoarsenClasses>&
   getEquivalenceClasses() const
   {
      return d_coarsen_classes;
   }

   /*!
    * @brief Set whether to unpack messages in a deterministic order.
    *
    * By default message unpacking is ordered by receive time, which
    * is not deterministic.  If your results are dependent on unpack
    * ordering and you want deterministic results, set this flag to
    * true.
    *
    * @param [in] flag
    */
   void
   setDeterministicUnpackOrderingFlag(
      bool flag);

   /*!
    * @brief Static function to set box intersection algorithm to use during
    * schedule construction for all CoarsenSchedule objects.
    *
    * If this method is not called, the default will be used.  If an invalid
    * string is passed, an unrecoverable error will result.
    *
    * @param[in] method  string identifying box intersection method.  Valid
    *                    choices are:  "DLBG" (default case),
    *                    and "ORIG_NSQUARED".   More details can be found below
    *                    in the comments for the generateSchedule() routine.
    *
    * @pre (method == "ORIG_NSQUARED") || (method == "DLBG")
    */
   static void
   setScheduleGenerationMethod(
      const std::string& method);

   /*!
    * @brief Set a pointer to a ScheduleOpsStrategy object
    *
    * @param strategy   Pointer to a concrete instance of ScheduleOpsStrategy
    */
   void setScheduleOpsStrategy(tbox::ScheduleOpsStrategy* strategy);

   /*!
    * @brief Print the coarsen schedule state to the specified data stream.
    *
    * @param[out] stream Output data stream.
    */
   void
   printClassData(
      std::ostream& stream) const;

private:
   CoarsenSchedule(
      const CoarsenSchedule&);              // not implemented
   CoarsenSchedule&
   operator = (
      const CoarsenSchedule&);              // not implemented

   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   //! @brief Mapping from a (potentially remote) Box to a set of neighbors.
   typedef std::map<hier::Box, hier::BoxContainer, hier::Box::id_less> FullNeighborhoodSet;

   /*!
    * @brief Main schedule generation routine for moving data from temporary
    * level to destination level.
    *
    * This function passes control to one of the algorithmic variations
    * of schedule generation based on the which method of generation is
    * selected.
    *
    * The resulting communication schedule will move source patch data from a
    * temporary coarse level (i.e., coarsened version of fine level) into the
    * destination patch data of the destination (coarse) level.
    *
    * The generateSchedule() routine invokes various versions of the schedule
    * generation process implemented in the similarly named routines below
    * based on the chosen schedule generation method. The different options
    * will not change the result of the application but may improve its
    * performance, especially for large numbers of processors.  Note that the
    * algorithm choice may be changed by calling the
    * setScheduleGenerationMethod() routine.
    *
    * The possibilities are as follows:
    *
    * <ul>
    *    <li>   if setScheduleGenerationMethod("DLBG") is called use
    *           generateScheduleDLBG() to generate the schedule.
    *           NOTE: THIS IS THE DEFAULT OPTION.
    *
    *    <li>   if setScheduleGenerationMethod("ORIG_NSQUARED") is called use
    *           generateScheduleNSquared() to generate the schedule.
    * </ul>
    *
    * @pre (s_schedule_generation_method == "ORIG_NSQUARED") ||
    *      (s_schedule_generation_method == "DLBG")
    */
   void
   generateSchedule();

   /*!
    * @brief Generate schedule using N^2 algorithms to determing box
    * intersections.
    *
    * This uses the original SAMRAI implementation which has a global view
    * of all distributed patches and checks every box against every other.
    */
   void
   generateScheduleNSquared();

   /*!
    * @brief Generate schedule using distributed data.
    *
    * This uses the DLBG distributed algortion to determine which source
    * patches contributed data to each destination patch and to compute
    * unfilled_boxes.
    */
   void
   generateScheduleDLBG();

   /*!
    * @brief Generate a temporary coarse level by coarsening the fine level.
    *
    * Note that this function does not allocate patch data storage.
    */
   void
   generateTemporaryLevel();

   /*!
    * @brief Set up refine algorithms to transfer coarsened data and to fill
    * temporary coarse level before coarsening operations, if needed.
    *
    * This is used when fill_coarse_data is set to true in the constructor.
    * The associated schedules are set in the generateSchedule() routine.
    */
   void
   setupRefineAlgorithm();

   /*!
    * @brief Coarsen source patch data from the fine patch level into the
    * source patch data on the coarse temporary patch level.
    *
    * @param[in] patch_strategy  Provides interface for user-defined functions
    *                            that may be used for coarsening.  Can be
    *                            null if no user-defined functions are needed.
    */
   void
   coarsenSourceData(
      CoarsenPatchStrategy* patch_strategy) const;

   /*!
    * @brief Calculate the maximum ghost cell width to grow boxes to check
    * for overlaps.
    */
   hier::IntVector
   getMaxGhostsToGrow() const;

   /*!
    * @brief Construct schedule transactions that communicate or copy coarsened
    * data from temporary coarse level to the destination level.
    *
    * @param[in] dst_level      The destination level for the schedule
    * @param[in] dst_box        Owned by a Patch on the destination level
    * @param[in] src_level      The temporary coarse level that will have
    *                           coarsened data
    * @param[in] src_box        Owned by a Patch on the temporary coarse level
    *
    * @pre dst_level
    * @pre src_level
    * @pre (d_crse_level.getDim() == dst_level->getDim()) &&
    *      (d_crse_level.getDim() == src_level->getDim()) &&
    *      (d_crse_level.getDim() == dst_box.getDim()) &&
    *      (d_crse_level.getDim() == src_box.getDim())
    */
   void
   constructScheduleTransactions(
      const std::shared_ptr<hier::PatchLevel>& dst_level,
      const hier::Box& dst_box,
      const std::shared_ptr<hier::PatchLevel>& src_level,
      const hier::Box& src_box);

   /*!
    * @brief Restructure the neighborhood sets from a src_to_dst Connector
    * so they can be used in schedule generation.
    *
    * First, this puts the neighborhood set data in src_to_dst into dst-major
    * order so the src owners can easily loop through the dst-src edges in the
    * same order that dst owners see them.  Transactions must have the same
    * order on the sending and receiving processors.
    *
    * Section, it shifts periodic image dst boxes back to the zero-shift
    * position, and applies a similar shift to src boxes so that the
    * overlap is unchanged.  The constructScheduleTransactions method requires
    * all shifts to be absorbed in the src box.
    *
    * The restructured neighboorhood sets are added to the output parameter.
    *
    * @param[out] full_inverted_edges
    * @param[in]  src_to_dst
    */
   void
   restructureNeighborhoodSetsByDstNodes(
      FullNeighborhoodSet& full_inverted_edges,
      const hier::Connector& src_to_dst) const;

   /*!
    * @brief Utility function to set up local copies of coarsen items.
    *
    * An array of coarsen data items obtained from the CoarsenClasses object
    * is stored locally here to facilitate interaction with transactions.
    *
    * @param[in] coarsen_classes
    */
   void
   setCoarsenItems(
      const std::shared_ptr<CoarsenClasses>& coarsen_classes);

   /*!
    * @brief Utility function to clear local copies of coarsen items.
    */
   void
   clearCoarsenItems();

   /*!
    * @brief Utility function to check coarsen items to see whether their
    * patch data components have sufficient ghost width.
    *
    * Specifically, the destination data must have a ghost cell width at least
    * as large as the coarsen item's d_gcw_to_coarsen data member.  The source
    * data must have a ghost cell width at least as large as d_gcw_to_coarsen
    * refined to the source (finer) level index space.  Also, if their are
    * any user-defined coarsening operations provided through the
    * CoarsenPatchStrategy interface, the source ghost cell width must be
    * at least as large as the stencil of those user-defined operations.
    *
    * If any of the tested ghost cell widths are insufficient, an error
    * will occur with a descriptive message.
    */
   void
   initialCheckCoarsenClassItems() const;

   /*!
    * @brief Selects algorithm used to generate communication schedule.
    */
   static std::string s_schedule_generation_method;

   /*!
    * @brief Shared debug checking flag.
    */
   static bool s_extra_debug;

   /*!
    * @brief Flag indicating if any RefineSchedule has read the input database
    * for static data.
    */
   static bool s_read_static_input;

   /*!
    * @brief Structures that store coarsen data items.
    */
   std::shared_ptr<CoarsenClasses> d_coarsen_classes;

   /*!
    * @brief number of coarsen data items
    */
   size_t d_number_coarsen_items;

   /*!
    * @brief used as array to store copy of coarsen data items.
    */
   const CoarsenClasses::Data** d_coarsen_items;

   /*!
    * @brief Cached pointers to the coarse, fine, and temporary patch levels.
    */
   std::shared_ptr<hier::PatchLevel> d_crse_level;
   std::shared_ptr<hier::PatchLevel> d_fine_level;
   std::shared_ptr<hier::PatchLevel> d_temp_crse_level;

   /*!
    * @brief Connector from coarse box_level to temporary
    * (coarsened fine) box_level.
    */
   std::shared_ptr<hier::Connector> d_coarse_to_temp;

   /*!
    * @brief Object supporting interface to user-defined spatial data
    * coarsening operations.
    */
   CoarsenPatchStrategy* d_coarsen_patch_strategy;

   /*!
    * @brief Factory object used to create data transactions when schedule is
    * constructed.
    */
   std::shared_ptr<CoarsenTransactionFactory> d_transaction_factory;

   /*!
    * @brief Cached ratio between source (fine) level and destination (coarse)
    * level.
    */
   hier::IntVector d_ratio_between_levels;

   /*!
    * @brief Source patch data indices for rapid data allocation/deallocation.
    */
   hier::ComponentSelector d_sources;

   /*!
    * @brief Level-to-level communication schedule between the temporary coarse
    * level and (actual) destination level.
    */
   std::shared_ptr<tbox::Schedule> d_schedule;

   /*!
    * @brief Boolean indicating whether source data on the coarse temporary
    * level must be filled before coarsening operations (see comments for class
    * constructor).
    */
   bool d_fill_coarse_data;

   /*!
    * @brief Algorithm used to set up schedule to fill temporary level if
    * d_fill_coarse_data is true.
    */
   std::shared_ptr<RefineAlgorithm> d_precoarsen_refine_algorithm;

   /*!
    * @brief Schedule used to fill temporary level if d_fill_coarse_data is
    * true.
    */
   std::shared_ptr<RefineSchedule> d_precoarsen_refine_schedule;

   //@{

   /*!
    * @brief Flag that turns on barrier calls for use in performance analysis.
    */
   static bool s_barrier_and_time;

   /*!
    * @name Timer objects for performance measurement.
    */
   static std::shared_ptr<tbox::Timer> t_coarsen_schedule;
   static std::shared_ptr<tbox::Timer> t_coarsen_data;
   static std::shared_ptr<tbox::Timer> t_gen_sched_n_squared;
   static std::shared_ptr<tbox::Timer> t_gen_sched_dlbg;
   static std::shared_ptr<tbox::Timer> t_coarse_data_fill;

   //*}

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#endif
