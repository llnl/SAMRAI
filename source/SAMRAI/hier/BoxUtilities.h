/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Routines for processing boxes within a domain of index space.
 *
 ************************************************************************/

#ifndef included_hier_BoxUtilities
#define included_hier_BoxUtilities

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/IntVector.h"

#include <list>
#include <vector>

namespace SAMRAI {
namespace hier {

/**
 * Class BoxUtilities provides several utility routines for processing
 * boxes or collections of boxes.  Many of these operations require
 * information about the location of the input boxes within some region of
 * index space (domain) or are used to compute this sort of information.
 * Often these routines are used in load balancing and communication routines
 * to determine the relationship of a box or set of boxes to some domain
 * boundary, where the domain is specified by a list of box regions.
 *
 * The following provides an explanation some of the concepts
 * common to many of the functions in this class.
 *
 *     - \b cut point
 *
 *        A cut point for a coordinate direction is an integer that
 *        specifies the index value of the cell immediately to the
 *        right of the boundary along which a box will be cut.  For
 *        instance, if we have a cut point value of 12 in the y-coordinate
 *        direction for box [(0,0,0), (5,15,10)], then the box
 *        will be cut into [(0,0,0),(5,11,10)] and [(0,12,0),(5,15,10)].
 *        That is, the box has been cut along the boundary between
 *        the y = 11 and y = 12 cells.
 *
 *     - \b bad cut point
 *        If a box is cut at a "bad cut point", the resulting boxes
 *        will violate some of the box constraints. For example, a
 *        cut point that does not satify the cut factor restriction
 *        or violates the bad_interval constraint (see below)
 *        is a bad cut point.
 *
 *     - \b irregular boundary
 *        An irregular boundary results when the domain cannot be
 *        described by a single box and the domain is non-convex.
 *
 *
 * The following provides an explanation some of the arguments
 * common to many of the functions in this class.
 *
 *    - \b min_size
 *       min_size is a IntVector that specifies the minimum
 *       allowable box length in each direction.  For example, if
 *       min_size = (10,4,15), then the minimum box length in the
 *       x, y, and z directions are 10, 4, and 15 respectively.
 *
 *    - \b max_size
 *       max_size is a IntVector that specifies the maximum
 *       allowable box length in each direction.  For example, if
 *       max_size = (10,40,50), then the maximum box length in the
 *       x, y, and z directions are 10, 40, and 50 respectively.
 *
 *       It should be noted that the max_size constraint has lower
 *       priority than the other constraints.  In instances where
 *       all constraints cannot be simultaneously satisfied, the max_size
 *       constraint is sacrificed.
 *
 *    - \b cut_factor
 *       cut_factor is a IntVector that constrains the
 *       size of a box to be multiples of the components of the
 *       cut_factor.  For instance, if cut_factor = (2,4,5), then the
 *       x, y, and z directions of a 8 box that satisfies the cut_factor
 *       constraint would be multiples of 2, 4, and 5 respectively.
 *
 *       This constraint is usually enforced with the cut_factor equal
 *       to a multiple of the refinement ratio between levels to ensure
 *       that if the box is coarsened, the resulting box is aligned with
 *       the coarse grid (it is assumed that the boundary of the fine box
 *       is aligned with a coarse cell boundary).
 *
 *    - \b bad_interval
 *       bad_interval is a IntVector that limits the distance
 *       a box can be from the boundary of the box list domain so that
 *       the outer boundaries of the box and the box list domain which are
 *       perpendicular to the i-th direction are no closer than the i-th
 *       component of bad_interval.  Another way to think of this is that
 *       the boundary of the box list domain in the i-th direction is not
 *       allowed to lie strictly within the interior of the box after
 *       it has been grown in the i-th direction by the i-th component
 *       of bad_interval.  For example, if bad_interval = (2,3,4), then
 *       the x, y, and z boundaries of the box must be at least 2, 3, and
 *       4 cells away from the x, y, and z boundaries of the box list
 *       domain.
 *
 *       The bad_interval constraint is enforced to avoid the situation
 *       where the ghost region for a box resides partially inside and
 *       partially outside the box list domain which complicates ghost
 *       cell filling.  In addition, this constraint avoids complicated
 *       issues with respect to the numerical accuracy of the solution.
 *
 *       Typically, bad_interval is based on the maximum ghost cell
 *       width over all patch data objects and some coarsen ratio.
 *
 *
 * Note that all member functions of this class are static.  The main intent
 * of the class is to group the functions into one name space.  Thus, you
 * should never attempt to instantiate a class of type BoxUtilities;
 * simply call the functions as static functions; e.g.,
 * BoxUtilities::function(...).  These routines are placed here rather
 * than in the box, box list, box array classes to avoid circular dependencies
 * among these classes.
 *
 * @see Box
 * @see BoxContainer
 */

struct BoxUtilities {
   /**
    * Check the given box for violation of minimum size, cut factor,
    * and box list domain constraints.  If a patch is generated from a box
    * that violates any of these constraints, then some other routine
    * (e.g., ghost cell filling, or inter-patch communication) may fail.
    * Thus, this routine prints an error message describing the violation and
    * executes a program abort.
    *
    * Arguments:
    *
    *    - \b box (input)
    *       box whose constraints are to be checked
    *
    *    - \b min_size (input)
    *       minimum allowed box size.  See class header for further
    *       description.
    *
    *    - \b cut_factor (input)
    *       See class header for description.
    *
    *    - \b bad_interval (input)
    *       See class header for description.
    *
    *       If there is no constraint on the box location within
    *       the box list domain, pass in an empty box array for the
    *       physical_boxes argument.
    *
    *    - \b physical_boxes (input)
    *       box array representing the index space of box list domain
    *
    * @pre (min_size.getDim() == cut_factor.getDim()) &&
    *      (min_size.getDim() == bad_interval.getDim())
    * @pre min_size > IntVector::getZero(min_size.getDim())
    * @pre cut_factor > IntVector::getZero(min_size.getDim())
    * @pre bad_interval >= IntVector::getZero(min_size.getDim())
    */
   static void
   checkBoxConstraints(
      const Box& box,
      const IntVector& min_size,
      const IntVector& cut_factor,
      const IntVector& bad_interval,
      const BoxContainer& physical_boxes);

   /**
    * Replace each box in the list that is too large with a list of non-
    * overlapping smaller boxes whose union covers the same region of
    * index space as the original box.
    *
    * Arguments:
    *
    *    - \b boxes (input)
    *       list of boxes to be chopped
    *
    *    - \b max_size (input)
    *       maximum allowed box size.  See class header for further
    *       description.
    *
    *    - \b min_size (input)
    *       minimum allowed box size.  See class header for further
    *       description.
    *
    *    - \b cut_factor (input)
    *       See class header for description.
    *
    *    - \b bad_interval (input)
    *       See class header for description.
    *
    *    - \b physical_boxes (input)
    *       box array representing the index space of box list domain
    *
    *
    * Notes:
    *
    *    - The resulting boxes will obey the minimum size and cut factor
    *      restrictions if the each of the original boxes does.
    *
    *    - Any box with side length not equal to a multiple of the
    *      cut factor for that direction, will not be chopped along that
    *      direction.
    *
    *    - The maximum size restriction may be sacrificed if the box
    *      cannot be chopped at appropriate points.  However, this is
    *      generally the case only when the box is adjacent to the
    *      box list domain boundary and an irregular boundary configuration
    *      restricts the cut locations or if the maximum size is not a
    *      multiple of the cut factor.
    *
    * @pre (max_size.getDim() == min_size.getDim()) &&
    *      (max_size.getDim() == cut_factor.getDim()) &&
    *      (max_size.getDim() == bad_interval.getDim())
    * @pre min_size > IntVector::getZero(min_size.getDim())
    * @pre max_size >= min_size
    * @pre cut_factor > IntVector::getZero(min_size.getDim())
    * @pre bad_interval >= IntVector::getZero(min_size.getDim())
    * @pre !physical_boxes.empty()
    * @pre !boxes.isOrdered()
    */
   static void
   chopBoxes(
      BoxContainer& boxes,
      const IntVector& max_size,
      const IntVector& min_size,
      const IntVector& cut_factor,
      const IntVector& bad_interval,
      const BoxContainer& physical_boxes);

   /**
    * Chop the box into a collection of boxes according to the collection
    * of cut points specified along each coordinate direction.  Cut points
    * that do not reside within the range of box indices are ignored.
    *
    * Arguments:
    *
    *    - \b boxes (output)
    *       list of boxes into which the "box" argument was chopped
    *
    *    - \b box (input)
    *       box which is to be chopped
    *
    *    - \b cut_points (input)
    *       cut_points is a vector of integer lists, each of which
    *       indicates the indices where the box will be cut in one of
    *       the coordinate directions
    *
    *
    * Assertion checks:
    *
    *    - The cut points for each direction must be on the list in
    *      increasing order.
    *
    *
    * Notes:
    *
    *    - The "boxes" BoxContainer is cleared before any box
    *      operations are performed.  Thus, any boxes on the list when
    *      the function is called will be lost.
    *
    * @pre cut_points.size() == box.getDim().getValue()
    */
   static void
   chopBox(
      BoxContainer& boxes,
      const Box& box,
      const std::vector<std::list<int> >& cut_points);

   /**
    * Extend the box in the list to domain boundary as needed so that
    * the domain boundary does not intersect the ghost cell region around
    * the box in an inappropriate manner.  Intersections that are
    * disallowed are those in which a portion of the domain boundary is
    * parallel to a box face and lies strictly in the interior of the ghost
    * cell box_level adjacent to that face.  In other words, we eliminate
    * ghost cell regions residing outside of a given domain and which are
    * narrower than the specified ghost width.   The boolean return value
    * is true if the input box was extended to the boundary and thus
    * is changed by the routine.  Otherwise, the return value is false.
    *
    * See description of bad_interval in the class header comments for
    * more details.
    *
    * Arguments:
    *
    *    - \b box (input/ouput)
    *       box to be extended
    *
    *    - \b domain (input)
    *       some domain whose interior is the union of boxes in a list
    *
    *    - \b ext_ghosts (input)
    *       IntVector that specifies the size of the desired
    *       ghost cell region
    *
    *
    * Notes:
    *
    *    - The ext_ghosts argument often corresponds to the bad_interval
    *      argument in many of the other functions in class.
    *
    *    - This operation may produce overlap regions among boxes on the
    *      list.
    *
    *    - There exist some bizarre domain configurations for which it is
    *      impossible to grow a box to the boundary and eliminate bad
    *      ghost region intersections.  This routine will extend each box as
    *      far as it can, but will not remedy these degenerate situations
    *      in general.
    *
    * @pre !domain.empty()
    * @pre ext_ghosts >= IntVector::getZero(ext_ghosts.getDim())
    */
   static bool
   extendBoxToDomainBoundary(
      Box& box,
      const BoxContainer& domain,
      const IntVector& ext_ghosts);

   /**
    * Same function as extendBoxToDomainBoundary() above except that it
    * extends each box in a list of boxes to the domain boundary specified
    * by the box list argument as needed.  The boolean return value
    * is true if any box in the input box list was extended to the boundary
    * and thus is changed by the routine.  Otherwise, the return value
    * is false.
    *
    * @pre !domain.empty()
    * @pre ext_ghosts >= IntVector::getZero(ext_ghosts.getDim())
    */
   static bool
   extendBoxesToDomainBoundary(
      BoxContainer& boxes,
      const BoxContainer& domain,
      const IntVector& ext_ghosts);

   /**
    * Grow each box in the list that is smaller than the specified minimum
    * size.
    *
    * Arguments:
    *
    *    - \b boxes (input/output)
    *       list of boxes to be grown to satisfy the min_size constraint
    *
    *    - \b domain (input)
    *       list of boxes whose union is some domain
    *
    *    - \b min_size (input)
    *       minimum allowed box size.  See class header for further
    *       description.
    *
    *
    * Notes:
    *
    *    - Each box that is grown must remain within the union of the
    *      boxes of the given domain.
    *
    *    - If the specified domain is an empty box list, then each box
    *      will be grown to be as large as the minimum size with no
    *      particular restrictions applied.
    *
    *    - This operation may produce overlap regions among boxes on list
    *
    *    - There exist some bizarre domain configurations for which it is
    *      impossible to grow a box sufficiently within the domain.
    *
    *      For instance if the domain is given by
    *      [(0,0),(2,10)], [(0,3),(1,4)], [(0,5),(10,10)]
    *      and the box is given by [(4,1),(6,2)] with a minimum size
    *      of (4,4), there is no way the box can be grown to the minimum
    *      size without have to "cross" the gap in the box list domain.
    *
    *      This routine will grow each box as far as it can, but will not
    *      remedy these situations, generally.
    *
    * @pre min_size > IntVector::getZero(min_size.getDim())
    */
   static void
   growBoxesWithinDomain(
      BoxContainer& boxes,
      const BoxContainer& domain,
      const IntVector& min_size);

   /**
    * Similar to growBoxesWithinDomain but works on one box at
    * a time and the domain is specified by the complement of local
    * parts of the domain.
    *
    * @pre min_size > IntVector::getZero(min_size.getDim())
    */
   static void
   growBoxWithinDomain(
      Box& box,
      const BoxContainer& local_domain_complement,
      const IntVector& min_size);

   /**
    * Determine whether the box may be chopped according to specified
    * min_size, max_size and cut_factor constraints.  For those
    * directions along which the box may be chopped, the cut points are
    * computed.  The cut points for the j-th coordinate direction are
    * placed into a list of integers corresponding to the j-th component
    * of the cut_point array.
    *
    * Return value:
    *
    *    - true is returned if the box may be chopped along any
    *      coordinate direction.  Otherwise, false is returned.
    *
    * Arguments:
    *
    *    - \b cut_points (output)
    *       vector of list of cut points for the box
    *
    *    - \b box (input)
    *       box to be cut
    *
    *    - \b max_size(input)
    *       minimum allowed box size.  See class header for further
    *       description.
    *
    *    - \b min_size(input)
    *       minimum allowed box size.  See class header for further
    *       description.
    *
    *    - \b cut_factor(input)
    *       See class header for description.
    *
    * @pre (max_size.getDim() == min_size.getDim()) &&
    *      (max_size.getDim() == cut_factor.getDim())
    * @pre min_size > IntVector::getZero(max_size.getDim())
    * @pre min_size <= max_size
    * @pre cut_factor > IntVector::getZero(max_size.getDim())
    */
   static bool
   findBestCutPointsGivenMax(
      std::vector<std::list<int> >& cut_points,
      const Box& box,
      const IntVector& max_size,
      const IntVector& min_size,
      const IntVector& cut_factor);

   /**
    * Determine whether the box may be chopped according to specified
    * min_size, max_size and cut_factor constraints along given
    * coordinate direction.  If the box may be chopped, the cut points
    * are computed and placed into a list of integers.
    *
    * Return value:
    *
    *    - true is returned if the box may be chopped along the specified
    *      coordinate direction.  Otherwise, false is returned.
    *
    * Arguments:
    *
    *    - \b idir (input)
    *       coordinate direction along which cut points will be computed
    *
    *    - \b cut_points (output)
    *       list of cut points for the box along the idir coordinate
    *       direction
    *
    *    - \b box (input)
    *       box to be chopped
    *
    *    - \b max_size (input)
    *       maximum allowed box size in idir coordinate direction.
    *
    *    - \b min_size (input)
    *       minimum allowed box size in idir coordinate direction.
    *
    *    - \b cut_factor (input)
    *       See class header for description.
    *
    * @pre !box.empty()
    * @pre min_size > 0
    * @pre max_size >= min_size
    * @pre cut_factor > 0
    */
   static bool
   findBestCutPointsForDirectionGivenMax(
      const tbox::Dimension::dir_t idir,
      std::list<int>& cut_points,
      const Box& box,
      const int max_size,
      const int min_size,
      const int cut_factor);

   /**
    * Determine whether the box may be chopped into the specified
    * number of boxes along each coordinate direction.  For those
    * directions along which the box may be chopped, the cut points are
    * computed.  The cut points for the j-th coordinate direction are
    * placed into a list of integers corresponding to the j-th component
    * of the cut_point array.
    *
    * Return value:
    *
    *    - true is returned if the box may be chopped along any
    *      coordinate direction.  Otherwise, false is returned.
    *
    * Arguments:
    *
    *    - \b cut_points (output)
    *       vector of list of cut points for the box
    *
    *    - \b box (input)
    *       box to be cut
    *
    *    - \b number_boxes (input)
    *       the i-th component of number_boxes specifies the desired
    *       number of cuts to be made along the i-th coordinate
    *       direction.
    *
    *    - \b min_size (input)
    *       minimum allowed box size. See class header for further
    *       description.
    *
    *    - \b cut_factor (input)
    *       See class header for description.
    *
    *
    * Important note: By convention, each integer cut point that is computed
    *                 corresponds to the cell index to the right of cut point.
    *
    * @pre (number_boxes.getDim() == min_size.getDim()) &&
    *      (number_boxes.getDim() == cut_factor.getDim())
    * @pre !box.empty()
    * @pre min_size > IntVector::getZero(number_boxes.getDim())
    * @pre number_boxes > IntVector::getZero(number_boxes.getDim())
    * @pre cut_factor > IntVector::getZero(number_boxes.getDim())
    */
   static bool
   findBestCutPointsGivenNumber(
      std::vector<std::list<int> >& cut_points,
      const Box& box,
      const IntVector& number_boxes,
      const IntVector& min_size,
      const IntVector& cut_factor);

   /**
    * Determine whether the box may be chopped into the specified
    * number of boxes along along given coordinate direction.  If the
    * box may be chopped, the cut points are computed and placed
    * into a list of integers.
    *
    * Return value:
    *
    *    - true is returned if the box may be chopped along the specified
    *      coordinate direction.  Otherwise, false is returned.
    *
    * Arguments:
    *
    *    - \b idir (input)
    *       coordinate direction along which cut points will be computed
    *
    *    - \b cut_points (output)
    *       list of cut points for the box along the idir coordinate
    *       direction
    *
    *    - \b box (input)
    *       box to be chopped
    *
    *    - \b num_boxes (input)
    *       num_boxes specifies the desired number of cuts to be made
    *       along the idir coordinate direction.
    *
    *    - \b min_size (input)
    *       minimum allowed box size in idir coordinate direction.
    *
    *    - \b cut_factor (input)
    *       See class header for description.
    *
    * @pre min_size > 0
    * @pre num_boxes > 0
    * @pre cut_factor > 0
    */
   static bool
   findBestCutPointsForDirectionGivenNumber(
      const tbox::Dimension::dir_t idir,
      std::list<int>& cut_points,
      const Box& box,
      const int num_boxes,
      const int min_size,
      const int cut_factor);

   /**
    * Determine whether box has any bad cut points based on its
    * position within the box list domain.  Information about the potentially
    * bad directions is returned in the IntVector
    * bad_cut_information.  An entry of zero indicates that there are no
    * bad cut points for the box along that coordinate direction.  An entry
    * of one indicates that there may be a bad cut point along that direction.
    *
    * Return value:
    *
    *    - true is returned if the box may potentially have a bad point
    *      along some coordinate direction. Otherwise false is returned.
    *
    * Arguments:
    *
    *    - \b bad_cut_information (output)
    *       A value of 0 in the i-th component of bad_cut_information
    *       indicates that there are no bad cut points in the i-th
    *       coordinate direction.
    *
    *    - \b box (input)
    *       box to be cut.
    *
    *    - \b physical_boxes (input)
    *       box array that represents some domain
    *
    *    - \b bad_interval (input)
    *       See class header for description.
    *
    * @pre (bad_cut_information.getDim() == box.getDim()) &&
    *      (bad_cut_information.getDim() == bad_interval.getDim())
    * @pre bad_interval >= IntVector::getZero(box.getDim())
    */
   static bool
   checkBoxForBadCutPoints(
      IntVector& bad_cut_information,
      const Box& box,
      const BoxContainer& physical_boxes,
      const IntVector& bad_interval);

   /**
    * Determine whether box may have any bad cut points along the specified
    * coordinate direction based on its position within the box array domain.
    *
    * Return value:
    *
    *    - true is returned if the box may potentially have a bad point;
    *      otherwise false is returned.
    *
    * Arguments:
    *
    *    - \b dir (input)
    *       coordinate direction to be checked for bad cut points
    *
    *    - \b box (input)
    *       box to be cut
    *
    *    - \b physical_boxes (input)
    *       box array that represents some domain
    *
    *    - \b bad_interval (input)
    *       See class header for description.
    *
    * @pre box.getDim() == bad_interval.getDim()
    * @pre !box.empty()
    * @pre bad_interval >= IntVector::getZero(box.getDim())
    */
   static bool
   checkBoxForBadCutPointsInDirection(
      const tbox::Dimension::dir_t dir,
      const Box& box,
      const BoxContainer& physical_boxes,
      const IntVector& bad_interval);

   /**
    * Determine bad cut points for box based on the specified box array domain
    * and bad interval.
    *
    * The cut information is returned as an array (size = dim) of arrays
    * (size = number of cells along edge of box) of boolean values.  A
    * false value indicates a good cut point, a true value indicates that
    * the box should not be cut at that point.
    *
    * Arguments:
    *
    *    - \b bad_cuts (output)
    *        stores an array of boolean arrays that indicates whether
    *        a potential cut point is bad.  A value of false indicates
    *        a good cut point, and a true value indicates a bad cut point.
    *
    *    - \b box (input)
    *        box to be cut
    *
    *    - \b physical_boxes (input)
    *        box array that represents some domain
    *
    *    - \b bad_interval (input)
    *        See class header for description.
    *
    * @pre !box.empty()
    * @pre bad_cuts.size() == box.getDim().getValue()
    * @pre bad_interval >= IntVector::getZero(box.getDim())
    */
   static void
   findBadCutPoints(
      std::vector<std::vector<bool> >& bad_cuts,
      const Box& box,
      const BoxContainer& physical_boxes,
      const IntVector& bad_interval);

   /**
    * Find bad cut points for a box given a single coordinate direction.
    * The cut information is returned as an array of boolean values
    * (size = number of cells along specified edge of box).  A false value
    * indicates a good cut point, a true value indicates that the box should
    * not be cut at that point.
    *
    * Arguments:
    *
    *    - \b dir (input)
    *       coordinate direction to be checked for bad cut points
    *
    *    - \b bad_cuts (output)
    *        boolean vector whose entries indicates whether
    *        a potential cut point is bad.
    *
    *    - \b box (input)
    *       box to be cut
    *
    *    - \b physical_boxes (input)
    *       box array that represents some domain
    *
    *    - \b bad_interval (input)
    *       See class header for description.
    *
    * @pre box.getDim() == bad_interval.getDim()
    * @pre !box.empty()
    * @pre bad_interval >= IntVector::getZero(box.getDim())
    */
   static void
   findBadCutPointsForDirection(
      const tbox::Dimension::dir_t dir,
      std::vector<bool>& bad_cuts,
      const Box& box,
      const BoxContainer& physical_boxes,
      const IntVector& bad_interval);

   /**
    * Given a set of potential cut points and a set of bad cut points for
    * a box, adjust the cut points so that they do not coincide with
    * bad cut points.  Typically, the cuts are generated using either of
    * the findBestCutPoints...() functions, and the bad cut points are
    * generated using the findBadCutPoints() function.
    *
    * Arguments:
    *
    *    - \b cuts (input/output)
    *       array of integer lists each of which holds a list of
    *       cut points for the box.  Each list is adjusted so that
    *       no cut points coincide with bad cut points
    *
    *    - \b bad_cuts (input)
    *       array of boolean arrays each of which stores information
    *       about which offsets from the lower corner of the box
    *       are bad cut points
    *
    *    - \b box (input)
    *       box to be cut
    *
    *    - \b min_size (input)
    *       minimum allowed box size.  See class header for further
    *       details.
    *
    *    - \b cut_factor (input)
    *       See class header for description.
    *
    *
    * Assertion checks:
    *
    *    - The cut points for each direction must be strictly increasing
    *      and all satisfy the cut_factor restriction.
    *
    * @pre (box.getDim() == min_size.getDim()) &&
    *      (box.getDim() == cut_factor.getDim())
    * @pre cuts.size() == box.getDim().getValue()
    * @pre bad_cuts.size() == box.getDim().getValue()
    * @pre !box.empty()
    * @pre min_size > IntVector::getZero(box.getDim())
    * @pre cut_factor > IntVector::getZero(box.getDim())
    * @pre for the ith array in bad_cuts, array.size() == box.numberCells(i)
    */
   static void
   fixBadCutPoints(
      std::vector<std::list<int> >& cuts,
      const std::vector<std::vector<bool> >& bad_cuts,
      const Box& box,
      const IntVector& min_size,
      const IntVector& cut_factor);

   /**
    * Given a set of potential cut points and a set of bad cut points for
    * a box, adjust the cut points in the specified coordinate direction
    * so that they do not coincide with bad cut points.  Typically, the
    * cuts are generated using either of the findBestCutPoints...()
    * functions, and the bad cut points are generated using the
    * findBadCutPoints() function.
    *
    * Arguments:
    *
    *    - \b dir (input)
    *       coordinate direction along which to fix cut points
    *
    *    - \b cuts (input/output)
    *       list of integers which holds a list of cut points for the box.
    *       This list is adjusted so that no cut points coincide with bad
    *       cut points.
    *
    *    - \b bad_cuts (input)
    *       array of booleans which stores information about which
    *       offsets from the lower corner of the box are bad cut points
    *
    *    - \b box (input)
    *       box to be cut
    *
    *    - \b min_size (input)
    *       minimum allowed box size along specified coordinate direction.
    *
    *    - \b cut_factor (input)
    *       See class header for description.
    *
    * @pre bad_cuts.size() == box.numberCells(dir)
    * @pre !box.empty()
    * @pre min_size > 0
    * @pre cut_factor > 0
    */
   static void
   fixBadCutPointsForDirection(
      const tbox::Dimension::dir_t dir,
      std::list<int>& cuts,
      const std::vector<bool>& bad_cuts,
      const Box& box,
      const int min_size,
      const int cut_factor);

   /**
    *
    * This function is called by findBadCutPoints(),
    * and the findBadCutPointsForDirection() member functions.  It sets bad
    * cut points near the lower and upper ends of the border box in the
    * given coordinate direction.
    *
    * @pre box.getDim() == border.getDim()
    * @pre (0 <= id) && (id < box.getDim().getValue())
    * @pre bad_cuts.size() == box.numberCells(id)
    * @pre bad_interval >= 0
    *
    */
   static void
   findBadCutPointsForBorderAndDirection(
      const tbox::Dimension::dir_t id,
      std::vector<bool>& bad_cuts,
      const Box& box,
      const Box& border,
      const int bad_interval);

   /**
    * Construct an array of box lists so that each list contains
    * a non-overlapping set of boxes covering some portion of the
    * box at the same array location in the box array.  The regions
    * of index space formed by composing the union of boxes on each
    * box list are mutually disjoint and the union of all boxes in
    * the box lists exactly covers the union of boxes in the original
    * box array.  In other words, this routine partitions the boxes
    * in the "boxes" argument into a set of non-overlapping box collections.
    * If none of the boxes in this box array overlap then each box list
    * in the resulting array has a single box equal to the corresponding
    * box in the box array.  This routine is especially useful for
    * determining a unique set of index points given an array of boxes
    * in some index space.
    *
    * Arguments:
    *
    *    - \b box_list_array (output)
    *       array of box lists which cover mutually exclusive portions
    *       of the index space covered by the "boxes" argument
    *
    *    - \b boxes (input)
    *       an arbitrary box array
    */
   static void
   makeNonOverlappingBoxContainers(
      std::vector<BoxContainer>& box_list_array,
      const BoxContainer& boxes);

   /*!
    * @brief Grow a box and chop it at block boundaries.
    *
    * If growing a box will cause it to extend across a block boundary, this
    * method will chop it into distinct parts that are stored in the output
    * BoxContainer.
    *
    * The output may be coarsened or refined from the input's index space,
    * controlled by the do_refine and do_coarsen arguments.  These arguments
    * may both be false, but at most one may be true.
    *
    * The boxes in the output container that are intersections with neighboring
    * blocks will be defined in the index space of those neighboring blocks
    *
    * @pre (do_refine != do_coarsen || (!do_refine && !do_coarsen))
    * @pre ratio_to_level_zero.getBlockSize() == grid_geom.getNumberBlocks()
    * @pre refine_coarsen_ratio.getBlockSize() == grid_geom.getNumberBlocks()
    *
    * @param[out] grown_boxes  Container to hold the results
    * @param[in]  box          Input box
    * @param[in]  grid_geom    Grid Geometry
    * @param[in]  ratio_to_level_zero  Ratio from input box to level 0
    * @param[in]  refine_coarsen_ratio Ratio to refine or coarsen output
    * @param[in]  grow_width   Width to grow before chopping
    * @param[in]  do_refine
    * @param[in]  do_coarsen
    */
   static void growAndAdjustAcrossBlockBoundary(
      BoxContainer& grown_boxes,
      const Box& box,
      const std::shared_ptr<const BaseGridGeometry>& grid_geom,
      const IntVector& ratio_to_level_zero,
      const IntVector& refine_coarsen_ratio,
      const IntVector& grow_width,
      bool do_refine,
      bool do_coarsen);


};

}
}

#endif
