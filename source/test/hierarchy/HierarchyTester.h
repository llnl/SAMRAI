/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Manager class for patch hierarchy refine/coarsen tests.
 *
 ************************************************************************/

#ifndef included_HierarchyTester
#define included_HierarchyTester

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"

#ifndef included_tbox_String
#include <string>

#define included_String
#endif

#include <memory>

using namespace SAMRAI;
using namespace tbox;
using namespace hier;
using namespace mesh;

namespace SAMRAI {
/**
 * Class HierarchyTester tests patch hierarchy coarsen/refine operations.
 * It sets up an initial patch hierarchy based on input file specifications.
 * Then, depending on how the test is specified in the input file, a
 * new patch hierarchy will be created which is either a coarsened version
 * or a refined version of the initial hierarchy.  Then, all aspects of
 * of the new hierarchy are compared against the initial hierarchy to
 * check whether the hierarchy coarsenin or refining were performed
 * correctly.
 *
 * The functions in this class called from main() are:
 * \begin{enumerate}
 *    - [setupInitialHierarchy(...)]
 *
 *    - [runHierarchyTestAndVerify(...)]
 * \end{enumerate}
 */

class HierarchyTester:public StandardTagAndInitStrategy
{
public:
   /**
    * Constructor initailizes test operations based on input database.
    */
   HierarchyTester(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<Database> hier_test_db);

   /**
    * Destructor deallocates internal storage.
    */
   ~HierarchyTester();

   /**
    * Set up initial hierarchy used to test coarsen/refine operations.
    */
   void
   setupInitialHierarchy(
      std::shared_ptr<Database> main_input_db);

   /**
    * After hierarchy refine/coarsen operations are performed, check results
    * and return integer number of test failures.
    */
   int
   runHierarchyTestAndVerify();

   /**
    * The following two functions are declared pure virtual in the
    * StandardTagAndInitStrategy base class.  Although they do nothing,
    * they must be defined here for proper operation of the
    * GriddingAlgorithm class which creates the initial patch hierarchy.
    */

   virtual void initializeLevelData(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      const bool initial_time,
      const std::shared_ptr<hier::PatchLevel>& old_level,
      const bool allocate_data)
   {
      NULL_USE(hierarchy);
      NULL_USE(level_number);
      NULL_USE(init_data_time);
      NULL_USE(can_be_refined);
      NULL_USE(initial_time);
      NULL_USE(old_level);
      NULL_USE(allocate_data);
   }

   void resetHierarchyConfiguration(
      const std::shared_ptr<PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level)
   {
      NULL_USE(hierarchy);
      NULL_USE(coarsest_level);
      NULL_USE(finest_level);
   }

private:
   /*
    * Object name for error reporting.
    */
   std::string d_object_name;

   const tbox::Dimension d_dim;

   /*
    * Booleans to indicate whether refine or coarsen is operation to test.
    * They are initialized to false and set via input values.
    */
   bool d_do_refine_test;
   bool d_do_coarsen_test;

   /*
    * Ratio to coarsen/refine hierarchy during test.
    */
   IntVector d_ratio;

   /*
    * Initial patch hierarchy set up based on input data and second hierarchy
    * generated during coarsen/refine operations.
    */
   std::shared_ptr<PatchHierarchy> d_initial_patch_hierarchy;
   std::shared_ptr<PatchHierarchy> d_test_patch_hierarchy;

   /*
    * Pointers to gridding algorithm object is cached in test object
    * so that test operates properly.  If this object goes out
    * of scope (i.e., is deleted) before test runs, boundary box
    * calculations will be incorrect since the test assumes
    * internal variables in gridding algorithm exist.
    */
   std::shared_ptr<GriddingAlgorithm> d_gridding_algorithm;

};

}
#endif
