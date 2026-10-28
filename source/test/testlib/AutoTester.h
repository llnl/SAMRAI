/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   (c) 1997-2026 Lawrence Livermore National Security, LLC
 *                Description:   Simple class used for autotesting.
 *
 ************************************************************************/

#ifndef included_AutoTesterXD
#define included_AutoTesterXD

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/algs/HyperbolicLevelIntegrator.h"
#include "SAMRAI/algs/MethodOfLinesIntegrator.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "BaselineData.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/algs/TimeRefinementIntegrator.h"

#include <memory>

using namespace SAMRAI;

/**
 *  Class AutoTester sets and verifies certain testing information
 *  used for autotesting the applications codes.
 *
 *  The following input parameters may be specified:
 *
 *
 *
 *     - \b test_fluxes (bool) specifies whether we will do
 *                 Riemann test or test on timesteps.
 *     - \b test_iter_num (int) iteration to carry out test.
 *     - \b correct_result (double array) holds correct result
 *     - \b output_correct (bool) specifies whether we will write
 *                 correct result (useful if changing problems
 *                 and want to set "correct" array).
 *
 *
 *
 */

class AutoTester
{
public:
   /**
    * Default constructor for AutoTester
    */
   AutoTester(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<tbox::Database> input_db);

   /**
    * Virtual destructor for AutoTester
    */
   virtual ~AutoTester();

   /**
    * Checks result for applications using TimeRefinementIntegrator
    * and HyperbolicLevelIntegrator.
    */
   int
   evalTestData(
      int iter,
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      const std::shared_ptr<algs::TimeRefinementIntegrator> tri,
      const std::shared_ptr<algs::HyperbolicLevelIntegrator> hli,
      const std::shared_ptr<mesh::GriddingAlgorithm> ga);

   /**
    * Checks result for applications using MethodOfLinesIntegrator
    */
   int
   evalTestData(
      int iter,
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      const double time,
      const std::shared_ptr<algs::MethodOfLinesIntegrator> mol,
      const std::shared_ptr<mesh::GriddingAlgorithm> ga);

   /**
    * Check boxes on specified hierarchy level against test boxes.
    */
   int
   checkHierarchyBoxes(
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      int ln,
      const hier::BoxLevel& correct_box_level,
      int iter);

   static int
   testHierarchyNeighbors(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy);

   static int
   testFlattenedHierarchy(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy);

private:
   /*
    *  Sets the parameters in the struct, based
    *  on data read from input.
    */
   void
   getFromInput(
      std::shared_ptr<tbox::Database> input_db);

   std::string d_object_name;

   const tbox::Dimension d_dim;

   bool d_test_fluxes;    //  specifies whether to check Riemann problem
   bool d_output_correct;   // output  result?
   int d_test_iter_num;     // iteration number to check result.

   std::vector<double> d_correct_result;  // array to hold correct values

   //!@brief Time steps at which to checkHierarchyBoxes().
   std::vector<int> d_test_patch_boxes_at_steps;
   //!@brief checkHierarchyBoxes() at d_test_patch_boxes_at_steps[d_test_patch_boxes_step_count].
   int d_test_patch_boxes_step_count;
   //!@brief Base name of files used in the run and regression test.
   std::string d_base_name;
   //!@brief Whether to write file of boxes for regression check.
   bool d_write_patch_boxes;
   //!@brief Whether to read file of boxes for regression check.
   bool d_read_patch_boxes;

   /*!
    * @brief Database containing correct results for the
    * BoxLevels generated, for comparison check.
    *
    * Database /step_number_SN/level_number_LN is the correct result
    * for the BoxLevel level number LN at the sequence number
    * SN.
    */
   testlib::BaselineData d_baseline_data;

};

#endif
