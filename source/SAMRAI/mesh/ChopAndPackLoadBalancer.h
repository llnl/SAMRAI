/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Load balance routines for uniform and non-uniform workloads.
 *
 ************************************************************************/

#ifndef included_mesh_ChopAndPackLoadBalancer
#define included_mesh_ChopAndPackLoadBalancer

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/ProcessorMapping.h"
#include "SAMRAI/mesh/LoadBalanceStrategy.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"

#include <memory>

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Class ChopAndPackLoadBalancer provides load balancing routines for
 * AMR hierarchy levels based on either uniform or non-uniform
 * workload estimates.
 *
 * This class is derived from the abstract base class
 * LoadBalanceStrategy; thus, it is a concrete implementation of
 * the load balance Strategy pattern interface.
 *
 * Load balancing operations, whether based on uniform or non-uniform
 * workloads, can be specified for each level in the hierarchy
 * individually or for the entire hierarchy.  Basic load balance
 * parameters can be set from an input file, while more complex
 * behavior can be set at run-time via member functions, including
 * dynamic reconfiguration of balance operations.
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b>
 *    - \b bin_pack_method
 *       String value indicating the type of bin packing to use to map patches
 *       to processors.  Currently, two options are supported: "GREEDY" and
 *       "SPATIAL".  The "GREEDY" method simply maps each patch (box) to the
 *       first processor (bin), in ascending tbox::MPI process number, whose
 *       difference between the average workload and its current workload is
 *       less than the workload of the patch in question.  The "SPATIAL" method
 *       first constructs an ordering of the patches (boxes) by passing a
 *       Morton-type curve through the center of each box.  Then, it attempts
 *       to map the patches to processors by assigning patches that are near
 *       each other on the curve to the same processor.  If no input value is
 *       specified, a default value of "SPATIAL" is used.  The input value
 *       will be used for all levels and will below.
 *
 *    - \b max_workload_factor
 *       Double array (length = number of levels) used during the box-splitting
 *       phase to determine which boxes to split.  Specifically, boxes will be
 *       chopped if their estimated workload is greater than
 *       max_workload_factor * A, where A is the average workload
 *       (i.e., A = (total work)/(num processors)).  The default value for this
 *       parameter is 1.0.  It can be set to any value greater than zero,
 *       either in the input file or via the setMaxWorkloadFactor() member
 *       function below.
 *
 *    - \b workload_tolerance
 *       Double array (length = number of levels) used during the box-splitting
 *       phase to determine which boxes to split.  The tolerance value can be
 *       use to prevent splitting of boxes when the computed box workload is
 *       close to the computed ideal workload.  A box is split if:<br>
 *       ( box_workload <= ( (1. + workload_tolerance) * ideal_workload ) )<br>
 *       Tolerance values should be greater than or equal to 0.0 and less
 *       then 1.0.  Large values will probably have undesirable results.
 *       It can be set either in the input file or via the
 *       setWorkloadTolerance() member function below.<br>
 *       NOTE: If a length is less than max levels then finest value
 *       specified is use for finer levels.  If length is greater
 *       than max levels, the values are ignored.
 *
 *    - \b ignore_level_box_union_is_single_box
 *       Boolean flag to control chopping of level boxes when the union of the
 *       input boxes passed to the loadBalanceBoxes() routine is a single box.
 *       The default value is false, which means that the domain will be
 *       chopped to make patch boxes based on the (single box) union of the
 *       boxes describing the level regardless of the input boxes.  When the
 *       value is set to true, either via the setIgnoreLevelDomainIsSingleBox()
 *       function or an input file, the domain will be chopped by chopping each
 *       of the input boxes.
 *
 *    - \b processor_layout
 *       Integer array (length = DIM) indicating the way in which the domain
 *       should be chopped when a level can be described as a single
 *       parallelepiped region (i.e., a box).  If no input value is provided,
 *       or if the product of these entries does not equal the number of
 *       processors, then the processor layout computed will be computed from
 *       the size of the domain box and the number of processors in use if
 *       necessary.<br>
 *       NOTE: The largest patch size constraint specified in the
 *       input for the GriddingAlgorithm object takes precedence
 *       over the processor layout specification.  That is, if the
 *       processor layout indicates that the resulting level patches
 *       would be larger than the largest patch size, the layout will
 *       be ignored and boxes obeying the patch size constrint will
 *       result.
 *
 *   - \b tile_size
 *   Tile size when using tile mode.  Tile mode restricts box cuts
 *   to tile boundaries.
 *
 * <b> Details: </b> <br>
 * <table>
 *   <tr>
 *     <th>parameter</th>
 *     <th>type</th>
 *     <th>default</th>
 *     <th>range</th>
 *     <th>opt/req</th>
 *     <th>behavior on restart</th>
 *   </tr>
 *   <tr>
 *     <td>bin_pack_method</td>
 *     <td>string</td>
 *     <td>"SPATIAL"</td>
 *     <td>"SPATIAL", "GREEDY"</td>
 *     <td>opt</td>
 *     <td>Not written to restart.  Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>max_workload_factor</td>
 *     <td>array of doubles</td>
 *     <td>none</td>
 *     <td>0.0 <= all values</td>
 *     <td>opt</td>
 *     <td>Not written to restart.  Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>workload_tolerance</td>
 *     <td>array of doubles</td>
 *     <td>none</td>
 *     <td>0.0 <= all values < 1.0</td>
 *     <td>opt</td>
 *     <td>Not written to restart.  Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>ignore_level_box_union_is_single_box</td>
 *     <td>bool</td>
 *     <td>FALSE</td>
 *     <td>TRUE, FALSE</td>
 *     <td>opt</td>
 *     <td>Not written to restart.  Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>processor_layout</td>
 *     <td>int[]</td>
 *     <td>N/A</td>
 *     <td>N/A</td>
 *     <td>opt</td>
 *     <td>Not written to restart.  Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>tile_size</td>
 *     <td>IntVector</td>
 *     <td>1</td>
 *     <td>1-</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 * </table>
 *
 * A sample input file entry might look like:
 *
 * @code
 *    processor_layout = 4 , 4 , 4    // number of processors is 64
 *    bin_pack = "GREEDY"
 *    max_workload_factor = 0.9
 *    ignore_level_box_union_is_single_box = TRUE
 * @endcode
 *
 * Performance warning: This class implements a sequential algorithm.
 * The time it takes to this balancer increases with processor count.
 * However, you can probably use this load balancer on up to 1K
 * processors before its performance degrades noticably.
 *
 * @see LoadBalanceStrategy
 */

class ChopAndPackLoadBalancer:
   public LoadBalanceStrategy
{

public:
   /*!
    * Construct load balancer object, including setting default object state
    * and reading input data from the input data base, if required.
    *
    * @param[in] dim
    * @param[in] name   User-defined string identifier used for error
    *                   reporting.  This string must be non-empty.
    * @param[in] input_db (optional) database pointer providing parameters from
    *                   input file.  This pointer may be null indicating no
    *                   input will be read.
    *
    * @pre !name.empty()
    */
   ChopAndPackLoadBalancer(
      const tbox::Dimension& dim,
      const std::string& name,
      const std::shared_ptr<tbox::Database>& input_db =
         std::shared_ptr<tbox::Database>());

   /*!
    * Construct load balancer object, including setting default object state
    * and reading input data from the input data base, if required.  The only
    * difference between this constructor and the previous one is the string
    * identifier input.  If this constructor is used, the default object name
    * "ChopAndPackLoadBalancer" applies.
    *
    * @param[in] dim
    * @param[in] input_db (optional) database pointer providing parameters from
    *                     input file.  This pointer may be null indicating no
    *                     input will be read.
    */
   explicit ChopAndPackLoadBalancer(
      const tbox::Dimension& dim,
      const std::shared_ptr<tbox::Database>& input_db =
         std::shared_ptr<tbox::Database>());

   /*!
    * The virtual destructor releases all internal storage.
    */
   virtual ~ChopAndPackLoadBalancer();

   /*!
    * Set the max workload factor for either the specified level or all
    * hierarchy levels.  See discussion about inputs above for information
    * on how this value is used during load balancing operations.
    *
    * @param factor        Double value of multiplier for average workload
    *                      used in box chopping.  The default value is 1.0.
    * @param level_number  Optional integer number for level to which factor
    *                      is applied. If no value is given, the factor will
    *                      be used for all levels.
    *
    * @pre factor > 0.0
    */
   void
   setMaxWorkloadFactor(
      double factor,
      int level_number = -1);

   /*!
    * Set the workload tolerance for either the specified level or all
    * hierarchy levels.  See discussion about inputs above for information
    * on how this value is used during load balancing operations.
    *
    * @param tolerance     Double value of tolerance. The default value is 0.0;
    *
    * @param level_number  Optional integer number for level to which factor
    *                      is applied. If no value is given, the value will
    *                      be used for all levels.
    *
    * @pre tolerance > 0.0
    */
   void
   setWorkloadTolerance(
      double tolerance,
      int level_number = -1);

   /*!
    * Configure the load balancer to use the data stored in the hierarchy at
    * the specified descriptor index for estimating the workload on each cell.
    *
    * @param data_id       Integer value of patch data identifier for workload
    *                      estimate on each cell.
    * @param level_number  Optional integer number for level on which data id
    *                      is used.  If no value is given, the data will be
    *                      used for all levels.
    *
    * @pre hier::VariableDatabase::getDatabase()->getPatchDescriptor()->getPatchDataFactory(data_id) is actually a  std::shared_ptr<pdat::CellDataFactory<double> >
    */
   void
   setWorkloadPatchDataIndex(
      int data_id,
      int level_number = -1);

   /*!
    * Configure the load balancer to load balance boxes by assuming all cells
    * on the specified level or all hierarchy levels are weighted equally.
    *
    * @param level_number  Optional integer number for level on which uniform
    *                      workload estimate will be used.  If the level
    *                      number is not specified, a uniform workload
    *                      estimate will be used on all levels.
    */
   void
   setUniformWorkload(
      int level_number = -1);

   /*!
    * Configure the load balancer to use the bin-packing procedure for
    * mapping patches to processors indicated by the string.
    *
    * @param method        String value indicating bin-packing method to use.
    *                      See input file description above for valid options.
    *                      The default value is "GREEDY".
    * @param level_number  Optional integer number for level on which
    *                      bin-packing method will be used. If no value is
    *                      given, the prescribed methods will be used on all
    *                      levels.
    *
    * @pre (method == "GREEDY") || (method == "SPATIAL")
    */
   void
   setBinPackMethod(
      const std::string& method,
      int level_number = -1);

   /*!
    * Set the boolean flag to control chopping of level boxes when the union of
    * the input boxes passed to the loadBalanceBoxes() routine is a single box.
    * The default value is false, which means that the domain will be chopped
    * to make patch boxes based on the (single box) union of the boxes describing
    * the level regardless of the input boxes.  When the value is set to true,
    * the domain will be chopped by chopping each of the input boxes.
    *
    * @param flag          Boolean value indicating whether to ignore the set of
    *                      input boxes to the loadBalanceBoxes() routine when the
    *                      union of those boxes is a single box.
    */
   void
   setIgnoreLevelDomainIsSingleBox(
      bool flag)
   {
      d_ignore_level_box_union_is_single_box = flag;
   }

   /*!
    * Return true if load balancing procedure for given level depends on
    * patch data on mesh; otherwise return false.  This can be used to
    * determine whether a level needs to be rebalanced although its box
    * configuration is unchanged.  This function is pure virtual in
    * the LoadBalanceStrategy base class.
    *
    * @return Boolean value indicating whether load balance depends on
    *         patch data.
    *
    * @param level_number  Integer patch level number.
    */
   bool
   getLoadBalanceDependsOnPatchData(
      int level_number) const;

   /*!
    * Given a list of boxes, representing the domain of a level in the AMR
    * hierarchy, generate an array of boxes and an associated processor
    * mapping from which the patches for the level will be generated and
    * assigned.  The resulting boxes and processor mapping will be determined
    * based on parameters set via input or member functions above.  This
    * function is pure virtual in the LoadBalanceStrategy base class.
    *
    * The load balancing algorithm should ignore any periodic image Boxes
    * in the input balance_box_level.
    *
    * @param balance_box_level
    * @param balance_to_anchor
    * @param hierarchy       Input patch hierarchy in which level will reside.
    * @param level_number    Input integer number of level in patch hierarchy.
    *                        This value must be >= 0.
    * @param min_size        Input integer vector of minimum sizes for
    *                        output boxes. All entries must be > 0.
    * @param max_size        Input integer vector of maximum sizes for
    *                        output boxes. All entries must be >= min_size.
    * @param domain_box_level
    * @param bad_interval    Input integer vector used to create boxes near
    *                        physical domain boundary with sufficient number
    *                        of cells.  No box face will be closer to the
    *                        boundary than the corresponding interval of cells
    *                        to the boundary (the corresponding value is given
    *                        by the normal direction of the box face) unless
    *                        the face coincides with the boundary itself.  The
    *                        point of this argument is to have no patch live
    *                        within a certain ghost cell width of the boundary
    *                        if its boundary does not coincide with that
    *                        boundary .  That is, all ghost cells along a face
    *                        will be either in the domain interior or outside
    *                        the domain.  All entries must be >= 0. See
    *                        hier::BoxUtilities documentation for more details.
    * @param cut_factor      Input integer vector used to create boxes with
    *                        correct sizes.  The length of each box
    *                        direction will be an integer multiple of the
    *                        corresponding cut factor vector entry.  All
    *                        vector entries must be > 0.  See hier::BoxUtilities
    *                        documentation for more details.
    * @param rank_group      Needed for compatibility with parent class.
    *                        This argument is ignored.
    *
    * @pre !balance_to_anchor || balance_to_anchor->hasTranspose()
    * @pre (d_dim == balance_box_level.getDim()) &&
    *      (d_dim == min_size.getDim()) && (d_dim == max_size.getDim()) &&
    *      (d_dim == domain_box_level.getDim()) &&
    *      (d_dim == bad_interval.getDim()) && (d_dim == cut_factor.getDim())
    */
   void
   loadBalanceBoxLevel(
      hier::BoxLevel& balance_box_level,
      hier::Connector* balance_to_anchor,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const hier::IntVector& min_size,
      const hier::IntVector& max_size,
      const hier::BoxLevel& domain_box_level,
      const hier::IntVector& bad_interval,
      const hier::IntVector& cut_factor,
      const tbox::RankGroup& rank_group = tbox::RankGroup()) const;

   /*!
    * Print out all members of the class instance to given output stream.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

   void
   printStatistics(
      std::ostream& output_stream = tbox::plog) const;

   /*!
    * Returns the name of this object.
    *
    * @return The name of this object.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

private:
   // The following are not implemented, but are provided here for
   // dumb compilers.
   ChopAndPackLoadBalancer(
      const ChopAndPackLoadBalancer&);
   ChopAndPackLoadBalancer&
   operator = (
      const ChopAndPackLoadBalancer&);

   /*
    * Read parameters from input database.
    */
   void
   getFromInput(
      const std::shared_ptr<tbox::Database>& input_db);

   /*!
    * Given a list of boxes, representing the domain of a level in the AMR
    * hierarchy, generate an array of boxes and an associated processor
    * mapping from which the patches for the level will be generated and
    * assigned.  The resulting boxes and processor mapping will be determined
    * based on parameters set via input or member functions above.  This
    * function is pure virtual in the LoadBalanceStrategy base class.
    *
    * @param out_boxes       Output box array for generating patches on level.
    * @param mapping         Output processor mapping for patches on level.
    * @param in_boxes        Input box list representing union of patches on level.
    * @param hierarchy       Input patch hierarchy in which level will reside.
    * @param level_number    Input integer number of level in patch hierarchy.
    *                        This value must be >= 0.
    * @param physical_domain Array of boxes describing the physical extent of
    *                        the problem domain in the index space associated
    *                        with the level.  This box array cannot be empty.
    * @param ratio_to_hierarchy_level_zero  Input integer vector indicating
    *                        ratio between index space of level to load balance
    *                        and hierarchy level 0 (i.e., coarsest hierarchy level).
    * @param min_size        Input integer vector of minimum sizes for
    *                        output boxes. All entries must be > 0.
    * @param max_size        Input integer vector of maximum sizes for
    *                        output boxes. All entries must be >= min_size.
    * @param cut_factor      Input integer vector used to create boxes with
    *                        correct sizes.  The length of each box
    *                        direction will be an integer multiple of the
    *                        corresponding cut factor vector entry.  All
    *                        vector entries must be > 0.  See hier::BoxUtilities
    *                        documentation for more details.
    * @param bad_interval    Input integer vector used to create boxes near
    *                        physical domain boundary with sufficient number
    *                        of cells.  No box face will be closer to the
    *                        boundary than the corresponding interval of cells
    *                        to the boundary (the corresponding value is given
    *                        by the normal direction of the box face) unless
    *                        the face coincides with the boundary itself.  The
    *                        point of this argument is to have no patch live
    *                        within a certain ghost cell width of the boundary
    *                        if its boundary does not coincide with that
    *                        boundary .  That is, all ghost cells along a face
    *                        will be either in the domain interior or outside
    *                        the domain.  All entries must be >= 0. See
    *                        hier::BoxUtilities documentation for more details.
    *
    * @pre (ratio_to_hierarchy_level_zero.getDim() == min_size.getDim()) &&
    *      (ratio_to_hierarchy_level_zero.getDim() == max_size.getDim()) &&
    *      (ratio_to_hierarchy_level_zero.getDim() == cut_factor.getDim()) &&
    *      (ratio_to_hierarchy_level_zero.getDim() == bad_interval.getDim())
    * @pre hierarchy
    * @pre level_number >= 0
    * @pre !physical_domain.empty()
    * @pre min_size > hier::IntVector::getZero(d_dim)
    * @pre max_size >= min_size
    * @pre cut_factor > hier::IntVector::getZero(d_dim)
    * @pre bad_interval >= hier::IntVector::getZero(d_dim)
    */
   void
   loadBalanceBoxes(
      hier::BoxContainer& out_boxes,
      hier::ProcessorMapping& mapping,
      const hier::BoxContainer& in_boxes,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      int level_number,
      const hier::BoxContainer& physical_domain,
      const hier::IntVector& ratio_to_hierarchy_level_zero,
      const hier::IntVector& min_size,
      const hier::IntVector& max_size,
      const hier::IntVector& cut_factor,
      const hier::IntVector& bad_interval) const;

   /*
    * Chop single box using uniform workload estimate.
    */
   void
   chopUniformSingleBox(
      hier::BoxContainer& out_boxes,
      std::vector<double>& out_workloads,
      const hier::Box& in_box,
      const hier::IntVector& min_size,
      const hier::IntVector& max_size,
      const hier::IntVector& cut_factor,
      const hier::IntVector& bad_interval,
      const hier::BoxContainer& physical_domain,
      const tbox::SAMRAI_MPI& mpi) const;

   /*
    * Chop boxes in list using uniform workload estimate.
    */
   void
   chopBoxesWithUniformWorkload(
      hier::BoxContainer& out_boxes,
      std::vector<double>& out_workloads,
      const hier::BoxContainer& in_boxes,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      int level_number,
      const hier::IntVector& min_size,
      const hier::IntVector& max_size,
      const hier::IntVector& cut_factor,
      const hier::IntVector& bad_interval,
      const hier::BoxContainer& physical_domain,
      const tbox::SAMRAI_MPI& mpi) const;

   /*
    * Chop boxes in list using non-uniform workload estimate.
    */
   void
   chopBoxesWithNonuniformWorkload(
      hier::BoxContainer& out_boxes,
      std::vector<double>& out_workloads,
      const hier::BoxContainer& in_boxes,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      int level_number,
      const hier::IntVector& ratio_to_coarsest_hierarchy_level,
      int wrk_indx,
      const hier::IntVector& min_size,
      const hier::IntVector& max_size,
      const hier::IntVector& cut_factor,
      const hier::IntVector& bad_interval,
      const hier::BoxContainer& physical_domain,
      const tbox::SAMRAI_MPI& mpi) const;

   /*
    * Map boxes to processors using chosen bin pack method.
    */
   void
   binPackBoxes(
      hier::BoxContainer& boxes,
      hier::ProcessorMapping& mapping,
      std::vector<double>& workloads,
      const std::string& bin_pack_method) const;

   /*!
    * @brief All-to-all communication of box arrays and associated weights.
    *
    * On invocation, each processor has a (possibly empty) array of
    * 'owned' boxes, and each box has a weight.  On return,
    * each processor has a array that contains all boxes owned
    * by all processors, and their associated weights.
    * If all processors input arrays have zero length, an error
    * is thrown.
    *
    * @pre box_list_in.size() == weights_in.size()
    */
   void
   exchangeBoxContainersAndWeightArrays(
      const hier::BoxContainer& box_list_in,
      std::vector<double>& weights_in,
      hier::BoxContainer& box_list_out,
      std::vector<double>& weights_out,
      const tbox::SAMRAI_MPI& mpi) const;

   /*
    * Utility functions to determine parameter values for level.
    */
   int
   getWorkloadDataId(
      int level_number) const
   {
      TBOX_ASSERT(level_number >= 0);
      return level_number < static_cast<int>(d_workload_data_id.size()) ?
             d_workload_data_id[level_number] :
             d_master_workload_data_id;
   }

   double
   getMaxWorkloadFactor(
      int level_number) const
   {
      TBOX_ASSERT(level_number >= 0);
      return level_number < static_cast<int>(d_max_workload_factor.size()) ?
             d_max_workload_factor[level_number] :
             d_master_max_workload_factor;
   }

   double
   getWorkloadTolerance(
      int level_number) const
   {
      TBOX_ASSERT(level_number >= 0);
      return level_number < static_cast<int>(d_workload_tolerance.size()) ?
             d_workload_tolerance[level_number] :
             d_master_workload_tolerance;
   }

   std::string
   getBinPackMethod(
      int level_number) const
   {
      TBOX_ASSERT(level_number >= 0);
      return level_number < static_cast<int>(d_bin_pack_method.size()) ?
             d_bin_pack_method[level_number] :
             d_master_bin_pack_method;
   }

   /*!
    * @brief Set up timers.
    */
   void
   setupTimers();

   /*
    * Object dimension.
    */
   const tbox::Dimension d_dim;

   /*
    * String identifier for load balancer object.
    */
   std::string d_object_name;

   /*
    * Specification of processor layout.
    */
   bool d_processor_layout_specified;
   hier::IntVector d_processor_layout;

   /*
    * Flag to control domain chopping when union of boxes for level
    * is a single box.
    */
   bool d_ignore_level_box_union_is_single_box;

   /*
    * Values for workload estimate data, workload factor, and bin pack method
    * that will be used for all levels unless specified for individual levels.
    */
   int d_master_workload_data_id;
   double d_master_max_workload_factor;
   double d_master_workload_tolerance;
   std::string d_master_bin_pack_method;

   /*
    * Values for workload estimate data, workload factor, and bin pack method
    * used on individual levels when specified as such.
    */
   std::vector<int> d_workload_data_id;
   std::vector<double> d_max_workload_factor;
   std::vector<double> d_workload_tolerance;
   std::vector<std::string> d_bin_pack_method;

   bool d_opt_for_single_box;

   /*!
    * @brief Tile size, when restricting cuts to tile boundaries,
    * Set to 1 when not restricting.
    */
   hier::IntVector d_tile_size;

   mutable std::vector<double> d_load_stat;

   /*
    * Performance timers.
    */
   std::shared_ptr<tbox::Timer> t_load_balance_box_level;
   std::shared_ptr<tbox::Timer> t_load_balance_boxes;
   std::shared_ptr<tbox::Timer> t_load_balance_boxes_remove_intersection;
   std::shared_ptr<tbox::Timer> t_get_global_boxes;
   std::shared_ptr<tbox::Timer> t_bin_pack_boxes;
   std::shared_ptr<tbox::Timer> t_bin_pack_boxes_sort;
   std::shared_ptr<tbox::Timer> t_bin_pack_boxes_pack;
   std::shared_ptr<tbox::Timer> t_chop_boxes;
};

}
}

#endif
