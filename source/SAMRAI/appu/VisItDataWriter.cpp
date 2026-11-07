/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Writes data files for visualization by VisIt
 *
 ************************************************************************/
#include "SAMRAI/appu/VisItDataWriter.h"

#ifdef HAVE_HDF5

#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/pdat/NodeDataFactory.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"


#include <cstring>
#include <ctime>
#include <vector>

extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

void SAMRAI_F77_FUNC(cpfdat2buf3d, CPFDAT2BUF3D) (
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const float *, double *, const int&);
void SAMRAI_F77_FUNC(cpddat2buf3d, CPDDAT2BUF3D) (
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double *, double *, const int&);
void SAMRAI_F77_FUNC(cpidat2buf3d, CPIDAT2BUF3D) (
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, double *, const int&);
}
extern "C" {
void SAMRAI_F77_FUNC(cpfdat2buf2d, CPFDAT2BUF2D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const float *, double *, const int&);
void SAMRAI_F77_FUNC(cpddat2buf2d, CPDDAT2BUF2D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double *, double *, const int&);
void SAMRAI_F77_FUNC(cpidat2buf2d, CPIDAT2BUF2D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, double *, const int&);
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace appu {

const float VisItDataWriter::VISIT_DATAWRITER_VERSION_NUMBER = 2.0;
const int VisItDataWriter::VISIT_NAME_BUFSIZE = 128;
const int VisItDataWriter::VISIT_UNDEFINED_INDEX = -1;
const int VisItDataWriter::VISIT_MASTER = 0;
const int VisItDataWriter::VISIT_FILE_CLUSTER_WRITE_BATON = 117;

bool VisItDataWriter::s_summary_file_opened = false;

tbox::StartupShutdownManager::Handler
VisItDataWriter::s_initialize_handler(
   VisItDataWriter::initializeCallback,
   0,
   0,
   VisItDataWriter::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

std::shared_ptr<tbox::Timer> VisItDataWriter::t_write_plot_data;

/*
 *************************************************************************
 *
 * The constructor --- sets default object state.
 *
 *************************************************************************
 */

VisItDataWriter::VisItDataWriter(
   const tbox::Dimension& dim,
   const std::string& object_name,
   const std::string& dump_directory_name,
   int number_procs_per_file,
   bool is_multiblock):
   d_dim(dim),
   d_mpi(MPI_COMM_NULL)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(number_procs_per_file == 1);

   if ((d_dim < tbox::Dimension(2)) || (d_dim > tbox::Dimension(3))) {
      TBOX_ERROR(
         "VisItDataWriter::VisItDataWriter"
         << "\n          VisItDataWriter only works for"
         << "\n          2D or 3D data" << std::endl);
   }

   d_object_name = object_name;

   d_mpi.dupCommunicator(tbox::SAMRAI_MPI::getSAMRAIWorld());

   d_default_derived_writer = 0;
   d_materials_writer = 0;

   d_number_working_slaves = VISIT_UNDEFINED_INDEX;
   d_file_cluster_size = number_procs_per_file;
   d_number_file_clusters = VISIT_UNDEFINED_INDEX;
   d_my_file_cluster_number = VISIT_UNDEFINED_INDEX;
   d_file_cluster_leader = false;
   d_my_rank_in_file_cluster = VISIT_UNDEFINED_INDEX;
   d_number_files_this_file_cluster = VISIT_UNDEFINED_INDEX;

   d_scaling_ratios.resize(1, hier::IntVector::getOne(dim));

   d_number_visit_variables = 0;
   d_number_visit_variables_plus_depth = 0;
   d_number_species = 0;

   d_time_step_number = VISIT_UNDEFINED_INDEX;
   d_grid_type = VISIT_CARTESIAN;
   d_top_level_directory_name = dump_directory_name;
   d_summary_filename = "summary.samrai";
   d_number_levels = 1;

   d_worker_min_max = 0;

   d_is_multiblock = is_multiblock;
   d_write_ghosts = false;
}

/*
 *************************************************************************
 *
 * The destructor implicitly deallocates the list of plot data items.
 *
 *************************************************************************
 */

VisItDataWriter::~VisItDataWriter()
{
   /*
    * De-allocate min/max structs for each variable.
    */
   if (d_worker_min_max != 0)
      delete[] d_worker_min_max;

   for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
        ipi != d_plot_items.end(); ++ipi) {
      for (int comp = 0; comp < VISIT_MAX_NUMBER_COMPONENTS; ++comp) {
         if (ipi->d_master_min_max[comp] != 0)
            delete[] ipi->d_master_min_max[comp];
      }
   }
   d_mpi.freeCommunicator();
}

/*
 *************************************************************************
 *
 * Register (non-derived) plot quantities: scalar, vector or tensor.
 *
 *************************************************************************
 */
void
VisItDataWriter::registerPlotQuantity(
   const std::string& variable_name,
   const std::string& variable_type,
   const int patch_data_index,
   const int start_depth_index,
   const double scale_factor,
   const std::string& variable_centering)
{
   TBOX_ASSERT(!variable_name.empty());
   TBOX_ASSERT(!variable_type.empty());
   TBOX_ASSERT(patch_data_index >= -1);
   TBOX_ASSERT(start_depth_index >= 0);

   /*
    * Check for name conflicts with existing registered variables.
    */
   for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
        ipi != d_plot_items.end(); ++ipi) {
      if (ipi->d_var_name == variable_name) {
         TBOX_ERROR("VisItDataWriter::registerPlotQuantity()"
            << "\n    Attempting to register variable with name "
            << variable_name << "\n    more than once." << std::endl);
      }
   }

   hier::IntVector ghost_width(hier::IntVector::getZero(d_dim));
   if (d_write_ghosts) {
      ghost_width =
         hier::VariableDatabase::getDatabase()->
         getPatchDescriptor()->
         getPatchDataFactory(patch_data_index)->
         getGhostCellWidth();
   }

   /*
    * Create a plot item and initialize its characteristics.
    */
   VisItItem plotitem;

   initializePlotItem(plotitem,
      variable_name,
      variable_type,
      patch_data_index,
      start_depth_index,
      scale_factor,
      variable_centering,
      ghost_width);

   ++d_number_visit_variables;
   d_number_visit_variables_plus_depth += plotitem.d_depth;
   d_plot_items.push_back(plotitem);

}

/*
 *************************************************************************
 *
 * Register derived plot quantities: scalar, vector, or tensor.  If no
 * derived data strategy is specified and no default derived data
 * strategy is set, an error will result.
 *
 *************************************************************************
 */

void
VisItDataWriter::registerDerivedPlotQuantity(
   const std::string& variable_name,
   const std::string& variable_type,
   VisDerivedDataStrategy* derived_writer,
   double scale_factor,
   const std::string& variable_centering,
   const std::string& variable_mix_type)
{
   TBOX_ASSERT(!variable_name.empty());
   TBOX_ASSERT(!variable_type.empty());

   /*
    * Check for name conflicts with existing registered variables.
    */
   for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
        ipi != d_plot_items.end(); ++ipi) {
      if (ipi->d_var_name == variable_name) {
         TBOX_ERROR("VisItDataWriter::registerDerivedPlotQuantity()"
            << "\n    Attempting to register variable with name "
            << variable_name << "\n    more than once." << std::endl);
      }
   }

   if (variable_name == "Coords") {
      TBOX_ASSERT(variable_type == "VECTOR");
      TBOX_ASSERT(variable_centering == "NODE");
      d_grid_type = VISIT_DEFORMED;
   }

   hier::IntVector ghost_width(hier::IntVector::getZero(d_dim));
   if (d_write_ghosts) {
      ghost_width =
         hier::VariableDatabase::getDatabase()->
         getPatchDescriptor()->getMaxGhostWidth(d_dim);
   }

   /*
    * Create a plot item and initialize its characteristics.
    */
   VisItItem plotitem;

   /*
    * Derived data is packed by the user.  Thus, we specify a dummy
    * patch data id here.
    */
   int patch_data_index = VISIT_UNDEFINED_INDEX;
   int start_depth_index = 0;
   initializePlotItem(plotitem,
      variable_name,
      variable_type,
      patch_data_index,
      start_depth_index,
      scale_factor,
      variable_centering,
      ghost_width);

   if (variable_name == "Coords") {
      plotitem.d_is_deformed_coords = true;

      /*
       * We need to reset the variable name, because it has to be written with
       * a special form to the VisIt readible HDF file.
       */
      char temp_buf[VISIT_NAME_BUFSIZE];
      for (int i = 0; i < plotitem.d_depth; ++i) {
         snprintf(temp_buf, VISIT_NAME_BUFSIZE, ".%02d", i);
         plotitem.d_visit_var_name[i] = variable_name + temp_buf;
      }

      /*
       * In this method, we assume the scale factor is always 1.0.  If a
       * user would like to choose a scale factor, use the
       * "registerSingleNodeCoordinate()" method.
       */
      plotitem.d_coord_scale_factor.resize(d_dim.getValue(), 1.0);
   }

   if (variable_mix_type == "MIXED") {
      plotitem.d_is_material_state_variable = true;
   }

   /*
    * Set characteristics for derived variable.
    */
   plotitem.d_is_derived = true;

   if (derived_writer == 0) {
      if (d_default_derived_writer == 0) {
         TBOX_ERROR("VisItDataWriter::registerDerivedPlotQuantity"
            << "\n    no derived data writer specified for variable:"
            << variable_name
            << "\n    and no default derived data writer set."
            << std::endl);
      } else {
         plotitem.d_derived_writer = d_default_derived_writer;
      }
   } else {
      plotitem.d_derived_writer = derived_writer;
   }

   ++d_number_visit_variables;
   d_number_visit_variables_plus_depth += plotitem.d_depth;
   d_plot_items.push_back(plotitem);
}

/*
 *************************************************************************
 *
 * Reset previously-registered scalar/vector variable to new data id
 * and depth index on the given level.  This allows the use of different
 * patch data ids for the same quantity on different hierarchy levels.
 * We check to make sure that the factory at the given index is
 * defined and consistent with the original registration.
 *
 *************************************************************************
 */

void
VisItDataWriter::resetLevelPlotQuantity(
   const std::string& variable_name,
   const int level_number,
   const int patch_data_index,
   const int start_depth_index)
{
   TBOX_ASSERT(!variable_name.empty());
   TBOX_ASSERT(level_number >= 0);
   TBOX_ASSERT(patch_data_index >= -1);
   TBOX_ASSERT(start_depth_index >= 0);

   /*
    * Verify the supplied patch data index has the same type and centering
    * as the plot item its replacing.
    */
   std::shared_ptr<hier::PatchDataFactory> factory(
      hier::VariableDatabase::getDatabase()->
      getPatchDescriptor()->
      getPatchDataFactory(patch_data_index));

   bool found_type = false;
   variable_data_type vdt = VISIT_DATA_TYPE_BAD;
   variable_centering vc = VISIT_CENTERING_BAD;

   if (!found_type) {
      std::shared_ptr<pdat::CellDataFactory<float> > ffactory(
         std::dynamic_pointer_cast<pdat::CellDataFactory<float>, hier::PatchDataFactory>(
            factory));
      if (ffactory) {
         vdt = VISIT_FLOAT;
         vc = VISIT_CELL;
         found_type = true;
      }
   }
   if (!found_type) {
      std::shared_ptr<pdat::NodeDataFactory<float> > ffactory(
         std::dynamic_pointer_cast<pdat::NodeDataFactory<float>, hier::PatchDataFactory>(
            factory));
      if (ffactory) {
         vdt = VISIT_FLOAT;
         vc = VISIT_NODE;
         found_type = true;
      }
   }

   if (!found_type) {
      std::shared_ptr<pdat::CellDataFactory<double> > dfactory(
         std::dynamic_pointer_cast<pdat::CellDataFactory<double>, hier::PatchDataFactory>(
            factory));
      if (dfactory) {
         vdt = VISIT_DOUBLE;
         vc = VISIT_CELL;
         found_type = true;
      }
   }
   if (!found_type) {
      std::shared_ptr<pdat::NodeDataFactory<double> > dfactory(
         std::dynamic_pointer_cast<pdat::NodeDataFactory<double>, hier::PatchDataFactory>(
            factory));
      if (dfactory) {
         vdt = VISIT_DOUBLE;
         vc = VISIT_NODE;
         found_type = true;
      }
   }
   if (!found_type) {
      std::shared_ptr<pdat::CellDataFactory<int> > ifactory(
         std::dynamic_pointer_cast<pdat::CellDataFactory<int>, hier::PatchDataFactory>(
            factory));
      if (ifactory) {
         vdt = VISIT_INT;
         vc = VISIT_CELL;
         found_type = true;
      }
   }
   if (!found_type) {
      std::shared_ptr<pdat::NodeDataFactory<int> > ifactory(
         std::dynamic_pointer_cast<pdat::NodeDataFactory<int>, hier::PatchDataFactory>(
            factory));
      if (ifactory) {
         vdt = VISIT_INT;
         vc = VISIT_NODE;
         found_type = true;
      }
   }
   if (!found_type) {
      TBOX_ERROR("VisItDataWriter::resetLevelPlotQuantity()"
         << "\n     Unable to determine type and centering"
         << "\n     of supplied patch data index."
         << "\n     ***Exiting" << std::endl);
   }

   /*
    * Find variable in the list of maintained vars.
    */
   bool found_var = false;
   for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
        (!found_var && ipi != d_plot_items.end()); ++ipi) {
      if (ipi->d_var_name == variable_name) {

         /*
          * Check to make sure supplied variable has same type and
          * centering as registered one.
          */
         if ((ipi->d_var_data_type != vdt) ||
             (ipi->d_var_centering != vc)) {
            TBOX_ERROR("VisItDataWriter::resetLevelPlotQuantity()"
               << "\n     The supplied patch data id has a different"
               << "\n     type and centering from the one originally"
               << "\n     registered.  hier::Variable name: "
               << variable_name
               << "\n     ***Exiting" << std::endl);
         }
         if (level_number >=
             static_cast<int>(ipi->d_level_patch_data_index.size())) {
            ipi->d_level_patch_data_index.resize(level_number + 1);
         }
         ipi->d_level_patch_data_index[level_number] = patch_data_index;
         ipi->d_start_depth_index = start_depth_index;
      }
   }

   if (!found_var) {
      TBOX_ERROR("VisItDataWriter::resetLevelPlotQuantity()"
         << "\n     Could not find the variable: "
         << variable_name
         << "\n     in the list of registered plot items."
         << "\n     Be sure registerPlotQuantity() has been"
         << "\n     called for the variable."
         << "\n     ***Exiting" << std::endl);
   }

}

/*
 *************************************************************************
 *
 * Register node coordinates of deformed (moving) grids.
 *
 *************************************************************************
 */

void
VisItDataWriter::registerNodeCoordinates(
   const int patch_data_index,
   const int start_depth_index)
{
   TBOX_ASSERT(patch_data_index >= -1);
   TBOX_ASSERT(start_depth_index >= 0);

   /*
    * Check to make sure "Coords" variable has not already been registered.
    */
   for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
        ipi != d_plot_items.end(); ++ipi) {
      if (ipi->d_var_name == "Coords") {
         TBOX_ERROR("VisItDataWriter::registerNodeCoordinates()"
            << "\n   Coordinates registered more than once." << std::endl);
      }
   }

   /*
    * Set the grid type for the visit data
    */
   d_grid_type = VISIT_DEFORMED;

   /*
    * Verify the supplied patch data index is a valid NODE-centered
    * float or double and has a depth of at least d_dim
    */
   std::shared_ptr<hier::PatchDataFactory> factory(
      hier::VariableDatabase::getDatabase()->
      getPatchDescriptor()->
      getPatchDataFactory(patch_data_index));

   bool found_type = false;
   int var_depth = VISIT_UNDEFINED_INDEX;
   if (!found_type) {

      std::shared_ptr<pdat::NodeDataFactory<float> > ffactory(
         std::dynamic_pointer_cast<pdat::NodeDataFactory<float>, hier::PatchDataFactory>(
            factory));
      if (ffactory) {
         var_depth = ffactory->getDepth();
         found_type = true;
      }
   }
   if (!found_type) {
      std::shared_ptr<pdat::NodeDataFactory<double> > dfactory(
         std::dynamic_pointer_cast<pdat::NodeDataFactory<double>, hier::PatchDataFactory>(
            factory));
      if (dfactory) {
         var_depth = dfactory->getDepth();
         found_type = true;
      }
   }
   if (!found_type) {
      TBOX_ERROR("VisItDataWriter::registerNodeCoordinates"
         << "\n     This variable is NOT a node centered"
         << "\n     float or double type, which is required."
         << "\n     ***Exiting" << std::endl);
   }

   int end_depth = start_depth_index + d_dim.getValue();
   if (var_depth < (end_depth)) {
      TBOX_ERROR("VisItDataWriter::registerNodeCoordinates"
         << "\n     This variable has depth: " << var_depth
         << "\n     It must be a VECTOR type and therefore"
         << "\n     have depth at least d_dim + start_depth_index = "
         << end_depth
         << "\n     ***Exiting" << std::endl);
   }

   hier::IntVector ghost_width(hier::IntVector::getZero(d_dim));
   if (d_write_ghosts) {
      ghost_width = factory->getGhostCellWidth();
   }

   /*
    * Create the coords plot item.
    */
   VisItItem plotitem;

   std::string var_name = "Coords";
   std::string var_type = "VECTOR";
   double scale_factor = 1.0;
   std::string var_cent = "NODE";

   initializePlotItem(plotitem,
      var_name,
      var_type,
      patch_data_index,
      start_depth_index,
      scale_factor,
      var_cent,
      ghost_width);

   plotitem.d_is_deformed_coords = true;

   /*
    * We need to reset the variable name, because it has to be written with
    * a special form to the VisIt readible HDF file.
    */
   char temp_buf[VISIT_NAME_BUFSIZE];
   for (int i = 0; i < plotitem.d_depth; ++i) {
      snprintf(temp_buf, VISIT_NAME_BUFSIZE, ".%02d", i);
      plotitem.d_visit_var_name[i] = var_name + temp_buf;
   }

   /*
    * In this method, we assume the scale factor is always 1.0.  If a
    * user would like to choose a scale factor, use the
    * "registerSingleNodeCoordinate()" method.
    */
   plotitem.d_coord_scale_factor.resize(d_dim.getValue(), 1.0);

   ++d_number_visit_variables;
   d_number_visit_variables_plus_depth += plotitem.d_depth;
   d_plot_items.push_back(plotitem);

}

/*
 *************************************************************************
 *
 * Register node coordinates of deformed (moving) grids.
 *
 *************************************************************************
 */

void
VisItDataWriter::registerSingleNodeCoordinate(
   const int coordinate_number,
   const int patch_data_index,
   const int depth_index,
   const double scale_factor)
{
   TBOX_ASSERT(coordinate_number >= 0 && coordinate_number < d_dim.getValue());
   TBOX_ASSERT(patch_data_index >= -1);
   TBOX_ASSERT(depth_index >= 0);

   /*
    * Set the grid type for the visit data
    */
   d_grid_type = VISIT_DEFORMED;

   /*
    * Verify the supplied patch data index is a valid NODE-centered
    * float or double
    */
   std::shared_ptr<hier::PatchDataFactory> factory(
      hier::VariableDatabase::getDatabase()->
      getPatchDescriptor()->
      getPatchDataFactory(patch_data_index));

   bool found_type = false;
   if (!found_type) {
      std::shared_ptr<pdat::NodeDataFactory<float> > ffactory(
         std::dynamic_pointer_cast<pdat::NodeDataFactory<float>, hier::PatchDataFactory>(
            factory));
      if (ffactory) {
         found_type = true;
      }
   }
   if (!found_type) {
      std::shared_ptr<pdat::NodeDataFactory<double> > dfactory(
         std::dynamic_pointer_cast<pdat::NodeDataFactory<double>, hier::PatchDataFactory>(
            factory));
      if (dfactory) {
         found_type = true;
      }
   }
   if (!found_type) {
      TBOX_ERROR("VisItDataWriter::registerSingleNodeCoordinate"
         << "\n     This variable is NOT a node centered"
         << "\n     float or double type, which is required."
         << "\n     ***Exiting" << std::endl);
   }

   /*
    * Create the coords plot item.  If its the first time this method
    * is called, initialize the "Coords" plot variable.  If it has
    * already been called before (i.e. coordinate_number > 0) then just
    * reset plot variable parameters as necessary.
    */
   if (coordinate_number == 0) {

      /*
       * Check to make sure "Coords" variable has not already been registered.
       */
      for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
           ipi != d_plot_items.end(); ++ipi) {
         if (ipi->d_var_name == "Coords") {
            TBOX_ERROR("VisItDataWriter::registerSingleNodeCoordinate()"
               << "\n   Coordinate registered more than once."
               << std::endl);
         }
      }

      hier::IntVector ghost_width(hier::IntVector::getZero(d_dim));
      if (d_write_ghosts) {
         ghost_width =
            hier::VariableDatabase::getDatabase()->
            getPatchDescriptor()->
            getPatchDataFactory(patch_data_index)->
            getGhostCellWidth();
      }

      VisItItem plotitem;

      std::string var_name = "Coords";
      std::string var_type = "SCALAR";
      std::string var_cent = "NODE";

      initializePlotItem(plotitem,
         var_name,
         var_type,
         patch_data_index,
         depth_index,
         scale_factor,
         var_cent,
         ghost_width);

      plotitem.d_is_deformed_coords = true;

      /*
       * We need to reset the variable name, because it has to be written with
       * a special form to the VisIt readible HDF file.
       */
      char temp_buf[VISIT_NAME_BUFSIZE];
      for (int i = 0; i < plotitem.d_depth; ++i) {
         snprintf(temp_buf, VISIT_NAME_BUFSIZE, ".%02d", i);
         plotitem.d_visit_var_name[i] = var_name + temp_buf;
      }

      plotitem.d_coord_scale_factor.resize(d_dim.getValue());
      plotitem.d_coord_scale_factor[coordinate_number] = scale_factor;
      ++d_number_visit_variables;
      d_number_visit_variables_plus_depth += plotitem.d_depth;

      d_plot_items.push_back(plotitem);

   } else {

      for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
           ipi != d_plot_items.end(); ++ipi) {

         if (ipi->d_is_deformed_coords) {

            ipi->d_var_type = VISIT_VECTOR;
            ipi->d_depth = d_dim.getValue();
            ipi->d_visit_var_name.resize(d_dim.getValue());

            std::string var_name = "Coords";
            char temp_buf[VISIT_NAME_BUFSIZE];
            snprintf(temp_buf, VISIT_NAME_BUFSIZE, ".%02d", coordinate_number);
            ipi->d_visit_var_name[coordinate_number] = var_name + temp_buf;

            ipi->d_coord_scale_factor[coordinate_number] = scale_factor;
            d_number_visit_variables_plus_depth += 1;
         }
      }
   }

}

/*
 *************************************************************************
 *
 * Register material names -- names of all the materials (not species)
 * being used in the application.
 *
 *************************************************************************
 */

void
VisItDataWriter::registerMaterialNames(
   const std::vector<std::string>& material_names)
{
   TBOX_ASSERT(material_names.size() > 0);

   /*
    * Check if we have already tried to register materials.
    */
   if (d_materials_names.size() > 0) {
      TBOX_ERROR("VisItDataWriter::registerMaterialNames"
         << "\n    This method has been called more than once."
         << "\n    The material names may not change during the"
         << "\n    simulation.  ***Exiting" << std::endl);
   }

   /*
    * Register each of the material names as a plot item with material
    * characteristics.
    */
   int num_materials = static_cast<int>(material_names.size());
   d_materials_names.resize(num_materials);
   for (int i = 0; i < num_materials; ++i) {
      if (material_names[i].empty()) {
         TBOX_ERROR("VisItDataWriter::registerMaterialNames"
            << "\n    Material: " << i
            << "\n    has an empty name.  A name must be supplied."
            << "\n    ***Exiting" << std::endl);
      }

      hier::IntVector ghost_width(hier::IntVector::getZero(d_dim));
      if (d_write_ghosts) {
         ghost_width =
            hier::VariableDatabase::getDatabase()->
            getPatchDescriptor()->getMaxGhostWidth(d_dim);
      }

      d_materials_names[i] = material_names[i];

      VisItItem plotitem;

      /*
       * We impose the criteria that materials must be scalar cell centered
       * double quantities.  The user must supply a method to pack
       * this information.  See header for more info.
       */
      std::string var_type = "SCALAR";
      int patch_data_index = VISIT_UNDEFINED_INDEX;
      int start_depth_index = 0;
      double scale_factor = 1.0;
      std::string var_cent = "CELL";

      initializePlotItem(plotitem,
         material_names[i],
         var_type,
         patch_data_index,
         start_depth_index,
         scale_factor,
         var_cent,
         ghost_width);

      plotitem.d_isa_material = true;
      plotitem.d_material_name = material_names[i];

      d_plot_items.push_back(plotitem);
   }
}

/*
 *************************************************************************
 *
 * Register material names -- names of all the materials (not species)
 * being used in the application. Volume fractions (and optionally state
 * variables) will be written in sparse arrays
 *
 *************************************************************************
 */

void
VisItDataWriter::registerSparseMaterialNames(
   const std::vector<std::string>& material_names)
{
   TBOX_ASSERT(material_names.size() > 0);

   /*
    * Check if we have already tried to register materials.
    */
   if (d_materials_names.size() > 0) {
      TBOX_ERROR("VisItDataWriter::registerSparseMaterialNames"
         << "\n    This method has been called more than once."
         << "\n    The material names may not change during the"
         << "\n    simulation.  ***Exiting" << std::endl);
   }

   /*
    * Register each of the material names as a plot item with material
    * characteristics.
    */
   int num_materials = static_cast<int>(material_names.size());
   d_materials_names.resize(num_materials);
   for (int i = 0; i < num_materials; ++i) {
      if (material_names[i].empty()) {
         TBOX_ERROR("VisItDataWriter::registerMaterialNames"
            << "\n    Material: " << i
            << "\n    has an empty name.  A name must be supplied."
            << "\n    ***Exiting" << std::endl);
      }
      d_materials_names[i] = material_names[i];
   }

   hier::IntVector ghost_width(hier::IntVector::getZero(d_dim));
   if (d_write_ghosts) {
      ghost_width =
         hier::VariableDatabase::getDatabase()->
         getPatchDescriptor()->getMaxGhostWidth(d_dim);
   }

   // Sparse Structure
   VisItItem plotitem;

   std::string var_type = "SCALAR";
   int patch_data_index = VISIT_UNDEFINED_INDEX;
   int start_depth_index = 0;
   double scale_factor = 1.0;
   std::string var_cent = "CELL";

   /*
    * this plotitem will write out the mat_list and material_packing_type
    * (=1) the value of mat_list for any zone will be either the material
    * number if the zone is clean, or a negative value that is the negative
    * index into mix_mat, vol_frac and next_mat (assuming the indexing of
    * mix_mat, etc. begin at 1) These arrays will only be written if there are
    * mixed zones.
    */
   initializePlotItem(plotitem,
      "materials",
      var_type,
      patch_data_index,
      start_depth_index,
      scale_factor,
      var_cent,
      ghost_width);
   plotitem.d_isa_material = true;
   plotitem.d_material_name = "sparse_material_list";

   // Use d_is_material_state_variable to mark this as the material list for
   //   the sparse data writing when volume fractions are written out
   plotitem.d_is_material_state_variable = true;

   d_plot_items.push_back(plotitem);
}

/*
 *************************************************************************
 *
 * Register species names for a given material.
 *
 *************************************************************************
 */

void
VisItDataWriter::registerSpeciesNames(
   const std::string& material_name,
   const std::vector<std::string>& species_names)
{
   TBOX_ASSERT(!material_name.empty());
   TBOX_ASSERT(species_names.size() > 0);

   /*
    * Be sure we have already registered materials.
    */
   if (d_materials_names.size() == 0) {
      TBOX_ERROR("VisItDataWriter::registerSpeciesNames"
         << "\n    No materials have yet been registered."
         << "\n    Be sure the 'registerMaterialNames()'"
         << "\n    is called before this method.  ***Exiting"
         << std::endl);
   }

   bool found_material = false;
   for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
        ipi != d_plot_items.end(); ++ipi) {
      if (ipi->d_material_name == material_name) {
         found_material = true;
      }
   }
   if (!found_material) {
      TBOX_ERROR("VisItDataWriter::registerSpeciesNames"
         << "\n    material name = " << material_name
         << "\n    has not been registered." << std::endl);
   }

   /*
    * Find the material in the list of plot items
    */
   VisItItem* material_item = 0;
   for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
        ipi != d_plot_items.end(); ++ipi) {
      if ((ipi->d_material_name == material_name) &&
          ipi->d_isa_material) {
         if (ipi->d_species_names.size() > 0) {
            TBOX_ERROR("VisItDataWriter::registerSpeciesNames"
               << "\n    material name = " << material_name
               << "\n    registerSpeciesNames has already been"
               << "\n    called for this material.  It is not"
               << "\n    possible to reset species during the"
               << "\n    simulation.  ***Exiting" << std::endl);
         }
         material_item = &(*ipi);
      }
   }

   TBOX_ASSERT(material_item != 0);

   d_number_species += static_cast<int>(species_names.size());
   material_item->d_species_names = species_names;

   hier::IntVector ghost_width(hier::IntVector::getZero(d_dim));
   if (d_write_ghosts) {
      ghost_width =
         hier::VariableDatabase::getDatabase()->
         getPatchDescriptor()->getMaxGhostWidth(d_dim);
   }

   /*
    * Create a plot variable for each species of the material.
    */
   for (int i = 0; i < static_cast<int>(species_names.size()); ++i) {

      if (species_names[i].empty()) {
         TBOX_ERROR("VisItDataWriter::registerSpeciesNames"
            << "\n    species i = " << i
            << "\n    for material: " << material_item->d_material_name
            << "\n    is empty.  ***Exiting" << std::endl);
      }

      /*
       * Species are associated with a material.  Hence, we give them
       * the same characteristics (centering, type, etc) as material
       * variables.
       */
      VisItItem plotitem;

      std::string var_type = "SCALAR";
      int patch_data_index = VISIT_UNDEFINED_INDEX;
      int start_depth_index = 0;
      double scale_factor = 1.0;
      std::string var_cent = "CELL";

      initializePlotItem(plotitem,
         species_names[i],
         var_type,
         patch_data_index,
         start_depth_index,
         scale_factor,
         var_cent,
         ghost_width);

      plotitem.d_isa_species = true;

      plotitem.d_species_name = species_names[i];
      plotitem.d_parent_material_pointer = material_item;

      d_plot_items.push_back(plotitem);

   }

}

/*
 *************************************************************************
 *
 * Register VisIt expressions to be embedded in datafile summary
 *   This method may be called multiple times to add more expressions
 *   as needed.
 *
 *************************************************************************
 */
void
VisItDataWriter::registerVisItExpressions(
   const std::vector<std::string>& expression_keys,
   const std::vector<std::string>& expressions,
   const std::vector<std::string>& expression_types)
{
   if ((expressions.size() > 0) &&
       (expressions.size() == expression_keys.size()) &&
       (expressions.size() == expression_types.size())) {
      int num_current_exp = static_cast<int>(d_visit_expressions.size());
      d_visit_expressions.resize(num_current_exp + expressions.size());
      d_visit_expression_keys.resize(num_current_exp + expressions.size());
      d_visit_expression_types.resize(num_current_exp + expressions.size());
      for (int i = 0; i < static_cast<int>(expressions.size()); ++i) {
         d_visit_expressions[num_current_exp + i] = expressions[i];
         d_visit_expression_keys[num_current_exp + i] = expression_keys[i];
         d_visit_expression_types[num_current_exp + i] = expression_types[i];
      }
   }
}

/*
 *************************************************************************
 *
 * Private method which initializes a VisIt variable based on user
 * input.
 *
 *************************************************************************
 */
void
VisItDataWriter::initializePlotItem(
   VisItItem& plotitem,
   const std::string& variable_name,
   const std::string& variable_type,
   const int patch_data_index,
   const int start_depth_index,
   const double scale_factor,
   const std::string& variable_centering,
   const hier::IntVector& ghost_width)
{
   TBOX_ASSERT(!variable_name.empty());
   TBOX_ASSERT(!variable_type.empty());
   TBOX_ASSERT(patch_data_index >= -1);
   TBOX_ASSERT(start_depth_index >= 0);

   plotitem.d_var_name = variable_name;

   /*
    * Set variable type.
    */
   if (variable_type == "SCALAR") {
      plotitem.d_var_type = VISIT_SCALAR;
      plotitem.d_depth = 1;
   } else if (variable_type == "VECTOR") {
      plotitem.d_var_type = VISIT_VECTOR;
      plotitem.d_depth = d_dim.getValue();
   } else if (variable_type == "TENSOR") {
      plotitem.d_var_type = VISIT_TENSOR;
      plotitem.d_depth = d_dim.getValue() * d_dim.getValue();
   } else {
      TBOX_ERROR("VisItDataWriter::registerPlotItem"
         << "\n    variable_type " << variable_type
         << "\n    is unsupported.  You must use SCALAR, VECTOR, or"
         << "\n    TENSOR.  Exiting***" << std::endl);
   }

   /*
    * Check to make sure we have not exceeded max allowed components.
    */
   int num_old_components = d_number_visit_variables_plus_depth
      + static_cast<int>(d_materials_names.size()) // number materials
      + d_number_species;

   int new_num_components = num_old_components + plotitem.d_depth;
   if (new_num_components > VISIT_MAX_NUMBER_COMPONENTS) {
      TBOX_ERROR("VisItDataWriter::registerPlotItem"
         << "\n     Unable to register this quantity because it"
         << "\n     the maximum number of variables allowed in"
         << "\n     the VisItWriter was reached:"
         << "\n       current num variables:"
         << num_old_components
         << "\n       variable depth: "
         << plotitem.d_depth
         << "\n     MAX_NUMBER_COMPONENTS: "
         << VISIT_MAX_NUMBER_COMPONENTS
         << "\n     Contact SAMRAI team for assistance." << std::endl);
   }

   /*
    * Set variable centering.  If the variable supplied a valid patch
    * index, we check the factory of the variable first to try to
    * determine centering from that.  If we cannot determine
    * the type, we use the "variable_centering" provided by the user.
    * Note that we only set the centering for variables that are not
    * derived, materials, or species, since these variables do not have
    * a supplied patch data index because they are packed by the user.
    */
   bool found_type = false;
   int var_depth = 0;
   if (patch_data_index >= 0) {

      std::shared_ptr<hier::PatchDataFactory> factory(
         hier::VariableDatabase::getDatabase()->
         getPatchDescriptor()->
         getPatchDataFactory(patch_data_index));
      if (!factory) {
         TBOX_ERROR("VisItDataWriter::registerPlotItem"
            << "\n    patch data array index = " << patch_data_index
            << "\n    for variable = " << variable_name
            << "\n    is invalid" << std::endl);
      } else {

         if (!found_type) {
            std::shared_ptr<pdat::CellDataFactory<float> > ffactory(
               std::dynamic_pointer_cast<pdat::CellDataFactory<float>,
                                           hier::PatchDataFactory>(factory));
            if (ffactory) {
               plotitem.d_var_centering = VISIT_CELL;
               plotitem.d_var_data_type = VISIT_FLOAT;
               var_depth = ffactory->getDepth();
               found_type = true;
            }
         }
         if (!found_type) {
            std::shared_ptr<pdat::CellDataFactory<double> > dfactory(
               std::dynamic_pointer_cast<pdat::CellDataFactory<double>,
                                           hier::PatchDataFactory>(factory));
            if (dfactory) {
               plotitem.d_var_centering = VISIT_CELL;
               plotitem.d_var_data_type = VISIT_DOUBLE;
               var_depth = dfactory->getDepth();
               found_type = true;
            }
         }
         if (!found_type) {
            std::shared_ptr<pdat::CellDataFactory<int> > ifactory(
               std::dynamic_pointer_cast<pdat::CellDataFactory<int>,
                                           hier::PatchDataFactory>(factory));
            if (ifactory) {
               plotitem.d_var_centering = VISIT_CELL;
               plotitem.d_var_data_type = VISIT_INT;
               var_depth = ifactory->getDepth();
               found_type = true;
            }
         }
         if (!found_type) {
            std::shared_ptr<pdat::NodeDataFactory<float> > ffactory(
               std::dynamic_pointer_cast<pdat::NodeDataFactory<float>,
                                           hier::PatchDataFactory>(factory));
            if (ffactory) {
               plotitem.d_var_centering = VISIT_NODE;
               plotitem.d_var_data_type = VISIT_FLOAT;
               var_depth = ffactory->getDepth();
               found_type = true;
            }
         }
         if (!found_type) {
            std::shared_ptr<pdat::NodeDataFactory<double> > dfactory(
               std::dynamic_pointer_cast<pdat::NodeDataFactory<double>,
                                           hier::PatchDataFactory>(factory));
            if (dfactory) {
               plotitem.d_var_centering = VISIT_NODE;
               plotitem.d_var_data_type = VISIT_DOUBLE;
               var_depth = dfactory->getDepth();
               found_type = true;
            }
         }
         if (!found_type) {
            std::shared_ptr<pdat::NodeDataFactory<int> > ifactory(
               std::dynamic_pointer_cast<pdat::NodeDataFactory<int>,
                                           hier::PatchDataFactory>(factory));
            if (ifactory) {
               plotitem.d_var_centering = VISIT_NODE;
               plotitem.d_var_data_type = VISIT_INT;
               var_depth = ifactory->getDepth();
               found_type = true;
            }
         }

         /*
          * Make sure variable depth is sufficient for the specified type
          * (SCALAR/VECTOR/TENSOR) with the start depth index.
          */
         int end_depth = start_depth_index + plotitem.d_depth;
         if (var_depth < end_depth) {
            TBOX_ERROR("VisItDataWriter::registerPlotItem"
               << "\n    The variable: " << variable_name
               << "\n    has insufficient depth for the type"
               << "\n    and start depth index registered."
               << "\n      var_type:           " << variable_type
               << "\n      start_depth_index:  " << start_depth_index
               << "\n      required min depth: " << end_depth
               << std::endl);
         }

      } // factory is Null

   } // valid patch data index

   if (!found_type) {
      if (variable_centering == "CELL") {
         plotitem.d_var_centering = VISIT_UNKNOWN_CELL;
      } else if (variable_centering == "NODE") {
         plotitem.d_var_centering = VISIT_UNKNOWN_NODE;
      } else {
         TBOX_ERROR("VisItDataWriter::registerPlotItem"
            << "\n     Unable to determine the centering for"
            << "\n     this variable; it must be supplied."
            << "\n     The variable_centering argument is "
            << variable_centering
            << "\n     Possible entries are CELL or NODE"
            << "\n     ***Exiting" << std::endl);
      }

      /*
       * When the var centering is specified by the user, we assume
       * it is of type double.  This implies the user should NOT try
       * to register their undefined variables with type float or int.
       */
      plotitem.d_var_data_type = VISIT_DOUBLE;
   }

   /*
    * Set the patch data index.
    */
   plotitem.d_patch_data_index = patch_data_index;
   plotitem.d_level_patch_data_index.resize(d_number_levels, patch_data_index);

   plotitem.d_visit_var_name.resize(plotitem.d_depth);
   char temp_buf[VISIT_NAME_BUFSIZE];
   for (int i = 0; i < plotitem.d_depth; ++i) {
      if (plotitem.d_depth == 1) {
         plotitem.d_visit_var_name[i] = variable_name;
      } else {
         snprintf(temp_buf, VISIT_NAME_BUFSIZE, ".%02d", i);
         plotitem.d_visit_var_name[i] = variable_name + temp_buf;
      }
   }

   plotitem.d_scale_factor = scale_factor;
   plotitem.d_start_depth_index = start_depth_index;

   /*
    * Initialize min/max information.
    */
   for (int i = 0; i < VISIT_MAX_NUMBER_COMPONENTS; ++i) {
      plotitem.d_master_min_max[i] = 0;
   }

   /*
    * Set derived, coords, material, and species information NULL.
    * If the variable is any of these types, this information
    * should be set by the appropriate registration functions.
    */
   plotitem.d_is_derived = false;
   plotitem.d_derived_writer = 0;
   plotitem.d_is_deformed_coords = false;

   plotitem.d_isa_material = false;
   plotitem.d_materials_writer = 0;

   plotitem.d_isa_species = false;
   plotitem.d_parent_material_pointer = 0;

   // default to CLEAN (not mixed data)
   plotitem.d_is_material_state_variable = false;

   plotitem.d_ghost_width.resize(d_dim.getValue());
   for (int d = 0; d < d_dim.getValue(); ++d) {
      plotitem.d_ghost_width[d] = tbox::MathUtilities<int>::Min(1,ghost_width[d]);
   }
}

/*
 *************************************************************************
 *
 * Private functions for parallel runs which serve as barriers to enable
 * orderly writing of cluster files by passing a baton from the current
 * writer to the next proc in cluster, and so on.
 *
 *************************************************************************
 */

void
VisItDataWriter::dumpWriteBarrierBegin()
{
   int x[1], proc_before_me, len = 1;

   if (d_file_cluster_leader) {
      return;
   } else {
      proc_before_me = (d_my_file_cluster_number * d_file_cluster_size)
         + d_my_rank_in_file_cluster - 1;

      if (d_mpi.getSize() > 1) {
         tbox::SAMRAI_MPI::Status status;
         d_mpi.Recv(x,
            len,
            MPI_INT,
            proc_before_me,
            VISIT_FILE_CLUSTER_WRITE_BATON,
            &status);
      }
   }
}

void
VisItDataWriter::dumpWriteBarrierEnd()
{
   int x[1], proc_after_me;
   int num_procs = d_mpi.getSize();
   proc_after_me = (d_my_file_cluster_number * d_file_cluster_size)
      + d_my_rank_in_file_cluster + 1;
   x[0] = 0;
   if (proc_after_me < num_procs) {
      if (d_mpi.getSize() > 1) {
         d_mpi.Send(x,
            1,
            MPI_INT,
            proc_after_me,
            VISIT_FILE_CLUSTER_WRITE_BATON);
      }
   }
}

/*
 *************************************************************************
 *
 * Write plot data from given hierarchy to HDF file
 *
 *************************************************************************
 */

void
VisItDataWriter::writePlotData(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int time_step_number,
   double simulation_time)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(time_step_number >= 0);
   TBOX_ASSERT(!d_top_level_directory_name.empty());

   /*
    * Currently, this class does not work unless the nodes
    * have globally sequentialized indices.  Check for these.
    */
   hier::BoxLevelConnectorUtils dlbg_edge_utils;

   for (int ln = 0; ln < hierarchy->getNumberOfLevels(); ++ln) {
      const hier::BoxLevel& unsorted_box_level =
         *hierarchy->getPatchLevel(ln)->getBoxLevel();
      std::shared_ptr<hier::BoxLevel> sorted_box_level;
      std::shared_ptr<hier::MappingConnector> unused_sorting_map;
      dlbg_edge_utils.makeSortingMap(
         sorted_box_level,
         unused_sorting_map,
         unsorted_box_level,
         false,
         true);
      if (!d_is_multiblock && *sorted_box_level != unsorted_box_level) {
         TBOX_ERROR(
            "VisItDataWriter: Encountered existing limitation of VisItDataWriter\n"
            << "This class cannot write files unless all patch levels have\n"
            << "globally sequentialized nodes.  This can be accomplished\n"
            << "by the sequentialize_patch_indices = TRUE input flag in\n"
            << "GriddingAlgorithm.  This problem can (and should\n"
            << "be fixed soon.");

      }
   }

   t_write_plot_data->start();

   if (time_step_number <= d_time_step_number) {
      TBOX_ERROR("VisItDataWriter::writePlotData"
         << "\n    data writer with name " << d_object_name
         << "\n    time step number: " << time_step_number
         << " is <= last time step number: " << d_time_step_number
         << std::endl);
   }
   d_time_step_number = time_step_number;

   if ((d_materials_names.size() > 0) && (d_materials_writer == 0)) {
      TBOX_ERROR("VisItDataWriter::writePlotData"
         << "\n    data writer with name " << d_object_name
         << "\n    setMaterialsDataWriter() has not been called,"
         << "\n    this method must be called when using materials."
         << std::endl);
   }

   d_number_levels = hierarchy->getNumberOfLevels();

   d_scaling_ratios.resize(d_number_levels,
      hier::IntVector(d_dim, 1, hierarchy->getNumberBlocks()));

   for (int ln = 1; ln <= hierarchy->getFinestLevelNumber(); ++ln) {
      d_scaling_ratios[ln] =
         hierarchy->getPatchLevel(ln)->getRatioToCoarserLevel();
   }

   if (d_top_level_directory_name.empty()) {
      TBOX_ERROR("VisItDataWriter::writePlotData"
         << "\n    data writer with name " << d_object_name
         << "\n     Dump Directory Name is not set" << std::endl);
   }

   int num_items_to_plot = d_number_visit_variables_plus_depth
      + static_cast<int>(d_materials_names.size()) // number of materials
      + d_number_species;

   if (num_items_to_plot == 0) {
      TBOX_ERROR("VisItDataWriter::writePlotData"
         << "\n    No VisIt variables have been registered."
         << std::endl);
   }

   initializePlotVariableMinMaxInfo(hierarchy);

   writeHDFFiles(hierarchy, simulation_time);

   t_write_plot_data->stop();
}

/*
 *************************************************************************
 *
 * Private function to initialize min/max information for the plot
 * components.  This method will allocate space for the d_mm array on
 * the VISIT_MASTER processsor (which holds min/max info for all plot
 * variables on all patches) and will allocate the d_worker_min_max array
 * on all processors except the VISIT_MASTER.  This latter array is used
 * to store data to be sent to the master when summary information is
 * written.
 *
 *************************************************************************
 */

void
VisItDataWriter::initializePlotVariableMinMaxInfo(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
   TBOX_ASSERT(hierarchy);

   /*
    * Compute max number of patches on this processor.
    */
   unsigned int number_local_patches = 0;
   unsigned int tot_number_of_patches = 0;

   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln) {
      std::shared_ptr<hier::PatchLevel> patch_level(
         hierarchy->getPatchLevel(ln));
      tot_number_of_patches += patch_level->getGlobalNumberOfPatches();
      for (hier::PatchLevel::iterator ip(patch_level->begin());
           ip != patch_level->end(); ++ip) {
         ++number_local_patches;
      }
   }

   int max_number_local_patches = number_local_patches;
   if (d_mpi.getSize() > 1) {
      d_mpi.AllReduce(&max_number_local_patches, 1, MPI_MAX);
   }

   /*
    * Determine number of worker processors.  NOTE: subract one because we
    * don't want to count processor zero.
    */
   int count_me_in = 0;
   if (number_local_patches > 0 || d_mpi.getRank() == VISIT_MASTER)
      count_me_in = 1;
   d_number_working_slaves = count_me_in;
   if (d_mpi.getSize() > 1) {
      d_mpi.AllReduce(&d_number_working_slaves, 1, MPI_SUM);
   }
   d_number_working_slaves -= 1;

   if (d_mpi.getRank() != VISIT_MASTER) {

      /*
       * Worker processor:  Allocate an array large enough to hold patch
       * min max information.  Pack array by var_item, component number,
       * level, and local patch number.
       */
      int num_items_to_plot = d_number_visit_variables_plus_depth
         + static_cast<int>(d_materials_names.size()) // number materials
         + d_number_species;

      int num_components = max_number_local_patches * num_items_to_plot;
      if (d_worker_min_max != 0) {
         delete[] d_worker_min_max;
      }
      if (num_components > 0) {
         d_worker_min_max = new patchMinMaxStruct[num_components];
      }
      memset((char *)d_worker_min_max, 0,
         num_components * sizeof(patchMinMaxStruct));
      for (int i = 0; i < num_components; ++i) {
         d_worker_min_max[i].patch_data_on_disk = false;
         d_worker_min_max[i].min = tbox::MathUtilities<double>::getMax();
         d_worker_min_max[i].max = tbox::MathUtilities<double>::getMin();
         d_worker_min_max[i].material_composition_code =
            VisMaterialsDataStrategy::VISIT_MIXED;
         d_worker_min_max[i].species_composition_code =
            VisMaterialsDataStrategy::VISIT_MIXED;
      }

   } else {   // (d_mpi.getRank() == VISIT_MASTER)

      /*
       * Master processor:  allocate array for each plot item to hold
       * min/max information for ALL patches, on all levels.
       */
      for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
           ipi != d_plot_items.end(); ++ipi) {

         for (int comp = 0; comp < ipi->d_depth; ++comp) {

            /*
             * Create space for master min/max struct, if it doesn't
             * already exist.
             */
            if (ipi->d_master_min_max[comp] != 0) {
               delete[] ipi->d_master_min_max[comp];
            }
            patchMinMaxStruct* mm =
               new patchMinMaxStruct[tot_number_of_patches];
            memset((char *)mm, 0,
               tot_number_of_patches * sizeof(patchMinMaxStruct));
            ipi->d_master_min_max[comp] = mm;

            for (unsigned int pn = 0; pn < number_local_patches; ++pn) {
               ipi->d_master_min_max[comp][pn].patch_data_on_disk = false;
               ipi->d_master_min_max[comp][pn].min =
                  tbox::MathUtilities<double>::getMax();
               ipi->d_master_min_max[comp][pn].max =
                  tbox::MathUtilities<double>::getMin();
               ipi->d_master_min_max[comp][pn].material_composition_code =
                  VisMaterialsDataStrategy::VISIT_MIXED;
               ipi->d_master_min_max[comp][pn].species_composition_code =
                  VisMaterialsDataStrategy::VISIT_MIXED;
            }
         }
      }
   } // proc == VISIT_MASTER

}

/*
 *************************************************************************
 *
 * Private function to coordinate writing HDF plot files.
 *
 *************************************************************************
 */

void
VisItDataWriter::writeHDFFiles(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   double simulation_time)
{

   TBOX_ASSERT(hierarchy);

// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif

   char temp_buf[VISIT_NAME_BUFSIZE];
   std::string dump_dirname;
   tbox::Database* visit_HDFFilePointer;

   int num_procs = d_mpi.getSize();
   int my_proc = d_mpi.getRank();

   if (d_file_cluster_size > num_procs) {
      d_file_cluster_size = num_procs;
   }
   d_my_file_cluster_number = my_proc / d_file_cluster_size;
   d_my_rank_in_file_cluster = my_proc % d_file_cluster_size;

   if (d_my_rank_in_file_cluster == 0) {
      d_file_cluster_leader = true;
   } else {
      d_file_cluster_leader = false;
   }

   d_number_file_clusters =
      static_cast<int>(ceil(static_cast<double>(num_procs)
                          / static_cast<double>(d_file_cluster_size)));
   d_number_files_this_file_cluster = d_file_cluster_size;

   if (d_my_file_cluster_number == (d_number_file_clusters - 1)) {
      // set d_number_files_this_file_cluster for last cluster
      d_number_files_this_file_cluster =
         num_procs - (d_file_cluster_size
                      * (static_cast<int>(
                            floor(static_cast<double>(num_procs)
                               / static_cast<double>(d_file_cluster_size)))));
      if (d_number_files_this_file_cluster == 0) {
         d_number_files_this_file_cluster = d_file_cluster_size;
      }
   }

   d_processor_in_file_cluster_number.resize(num_procs);
   for (int i = 0; i < num_procs; ++i) {
      d_processor_in_file_cluster_number[i] = i / d_file_cluster_size;
   }

   snprintf(temp_buf, VISIT_NAME_BUFSIZE, "%05d", d_time_step_number);
   d_current_dump_directory_name = "visit_dump.";
   d_current_dump_directory_name += temp_buf;
   if (!d_top_level_directory_name.empty() &&
       d_top_level_directory_name[d_top_level_directory_name.length() - 1] == '/') {
      dump_dirname = d_top_level_directory_name;
   } else {
      dump_dirname = d_top_level_directory_name + "/";
   }
   dump_dirname = dump_dirname + d_current_dump_directory_name;
   tbox::Utilities::recursiveMkdir(dump_dirname);

// The baton barrier implementation seems to be buggy, and can affect
// other mpi writes and reads.

//#define USE_BATON_BARRIERS

#ifdef USE_BATON_BARRIERS
   dumpWriteBarrierBegin();
#endif
   {
      // cluster_leader guaranteed to enter this section before anyone else
      snprintf(temp_buf, VISIT_NAME_BUFSIZE, "/processor_cluster.%05d.samrai",
         d_my_file_cluster_number);
      std::string database_name(temp_buf);
      std::string visit_HDFFilename = dump_dirname + database_name;
      visit_HDFFilePointer = new tbox::HDFDatabase(database_name);
      if (d_file_cluster_leader) {

         // creates the HDF file:
         //      dirname/visit_dump.000n/processor_cluster.000m.samrai
         //      where n is timestep #, m is processor number
         visit_HDFFilePointer->create(visit_HDFFilename);

      } else {
         // file already created other procs just need to open it
         const bool read_write_mode(true);
         if (!visit_HDFFilePointer->open(visit_HDFFilename, read_write_mode)) {
            TBOX_ERROR("VisItDataWriter::writeHDFFiles"
               << "\n    data writer with name " << d_object_name
               << "\n    Error attempting to open visit file "
               << visit_HDFFilename << std::endl);
         }
      }

      // create group for this proc
      snprintf(temp_buf, VISIT_NAME_BUFSIZE, "processor.%05d", my_proc);
      std::shared_ptr<tbox::Database> processor_HDFGroup(
         visit_HDFFilePointer->putDatabase(std::string(temp_buf)));
      writeVisItVariablesToHDFFile(processor_HDFGroup,
         hierarchy,
         0,
         hierarchy->getFinestLevelNumber(),
         simulation_time);
      visit_HDFFilePointer->close(); // invokes H5FClose
      delete visit_HDFFilePointer; // deletes tbox::HDFDatabase object
   }

#ifdef USE_BATON_BARRIERS
   dumpWriteBarrierEnd();
#endif

   /*
    * When using DLBG, the globalized data is not saved by default,
    * so it must be generated, requiring communication.
    * To avoid mixing these communications and those required
    * in writing plot data, execute a function to compute and
    * cache the globalized data.
    */
   for (int ln = 0; ln < hierarchy->getNumberOfLevels(); ++ln) {
      hierarchy->getPatchLevel(ln)->getBoxes();
   }

   tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier();

   writeSummaryToHDFFile(dump_dirname,
      hierarchy,
      0,
      hierarchy->getFinestLevelNumber(),
      simulation_time);
}

/*
 *************************************************************************
 *
 * Private function to find global patch number, given level number and
 * local patch number. This is needed because SAMRAI maintains the patch
 * number on each level, but VisIt needs a unique number for each patch
 * written.
 *
 *************************************************************************
 */

int
VisItDataWriter::getGlobalPatchNumber(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const int patch_number)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(level_number >= 0);
   TBOX_ASSERT(patch_number >= 0);

   int global_patch_id = 0;

   for (int i = 0; i < level_number; ++i) {
      global_patch_id +=
         hierarchy->getPatchLevel(i)->getGlobalNumberOfPatches();
   }

   global_patch_id += patch_number;
   return global_patch_id;
}

/*
 *************************************************************************
 *
 * Private function to write variables & materials data to an HDF File.
 *
 *************************************************************************
 */

void
VisItDataWriter::writeVisItVariablesToHDFFile(
   const std::shared_ptr<tbox::Database>& processor_HDFGroup,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int coarsest_level,
   int finest_level,
   double simulation_time)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(coarsest_level >= 0);
   TBOX_ASSERT(finest_level >= 0);

   /*
    * Reset the var_id_ctr - this is used to record min/max summary
    * information for every plotted variable on the patch.  It is incremented
    * for each component (i.e. depth), of each variable, of each patch.
    */
   d_var_id_ctr = 0;

   char temp_buf[VISIT_NAME_BUFSIZE];
   std::shared_ptr<tbox::Database> level_HDFGroup, patch_HDFGroup;

   for (int ln = coarsest_level; ln <= finest_level; ++ln) {

      /*
       * create new HDFGroup for this level
       */
      snprintf(temp_buf, VISIT_NAME_BUFSIZE, "level.%05d", ln);
      level_HDFGroup = processor_HDFGroup->putDatabase(std::string(temp_buf));

      std::shared_ptr<hier::PatchLevel> patch_level(
         hierarchy->getPatchLevel(ln));

      hier::IntVector coarsen_ratio(patch_level->getRatioToCoarserLevel());

      for (hier::PatchLevel::iterator ip(patch_level->begin());
           ip != patch_level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& patch = *ip;

         /*
          * create new HDFGroup for this patch
          */
         int pn = patch->getLocalId().getValue();
         snprintf(temp_buf, VISIT_NAME_BUFSIZE, "patch.%05d", pn);
         patch_HDFGroup = level_HDFGroup->putDatabase(std::string(temp_buf));

         int curr_var_id_ctr = d_var_id_ctr;
         packRegularAndDerivedData(
            patch_HDFGroup, hierarchy, ln, *patch, simulation_time);

         if (d_materials_names.size() > 0) {
            d_var_id_ctr = curr_var_id_ctr;
            packMaterialsData(patch_HDFGroup,
               hierarchy,
               ln,
               *patch);

            d_var_id_ctr = curr_var_id_ctr;
            packSpeciesData(hierarchy,
               ln,
               *patch);
         }

      }
   }

   /*
    * Clean up from packing operations.  If this is not done there is a dangling smart
    * pointer reference to HDF5 groups and the file may not be written/closed.
    */
   for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
        ipi != d_plot_items.end(); ++ipi) {
      ipi->d_species_HDFGroup.reset();
      ipi->d_extents_species_HDFGroup.reset();
   }
}

/*
 *************************************************************************
 *
 * Private function to pack regular and derived VisIt variables into
 * specified HDF database.
 *
 *************************************************************************
 */

void
VisItDataWriter::packRegularAndDerivedData(
   const std::shared_ptr<tbox::Database>& patch_HDFGroup,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   hier::Patch& patch,
   double simulation_time)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(level_number >= 0);

   /*
    * Loop over variables and write out those that are NOT
    * material or species variables.
    */
   for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
        ipi != d_plot_items.end(); ++ipi) {

      /*
       * Only write regular (non-derived) and derived vars, not
       * material or spacies vars.
       */
      if (!(ipi->d_isa_material || ipi->d_isa_species)) {

         /*
          * create buffer to hold patch data
          */
         int buf_size = getBufferSize(patch.getBox(),
               hier::IntVector(ipi->d_ghost_width),
               ipi->d_var_centering);

         hier::Box pack_box(patch.getBox());
         if (d_write_ghosts) {
            pack_box.grow(hier::IntVector(ipi->d_ghost_width));
         }

         double* dbuffer = new double[buf_size]; // used to pack var
         float* fbuffer = new float[buf_size]; // copy to float for writing

         // Check for mixed/clean state variables
         if (!(ipi->d_is_material_state_variable)) { // Conventional variable
            for (int depth_id = 0; depth_id < ipi->d_depth; ++depth_id) {

               /*
                * If its derived data, pack via the derived writer.
                * Otherwise, pack with local private method.
                */
               bool data_exists_on_patch = false;
               int patch_data_id = VISIT_UNDEFINED_INDEX;
               if (ipi->d_is_derived) {

                  // derived data
                  data_exists_on_patch =
                     ipi->d_derived_writer->
                     packDerivedDataIntoDoubleBuffer(
                        dbuffer,
                        patch,
                        pack_box,
                        ipi->d_var_name,
                        depth_id,
                        simulation_time);

               } else {

                  /*
                   * Check if patch data id has been reset on the level.  If
                   * not, just use the original registered data id.
                   */
                  patch_data_id = ipi->d_patch_data_index;
                  if (static_cast<int>(ipi->d_level_patch_data_index.size()) >
                      level_number) {
                     patch_data_id =
                        ipi->d_level_patch_data_index[level_number];
                  }

                  data_exists_on_patch =
                     patch.checkAllocated(patch_data_id);

                  if (data_exists_on_patch) {
                     int new_depth_id = ipi->d_start_depth_index + depth_id;

                     // regular (non-derived) data
                     packPatchDataIntoDoubleBuffer(
                        patch.getPatchData(patch_data_id),
                        new_depth_id,
                        ipi->d_var_data_type,
                        patch.getBox(),
                        dbuffer,
                        ipi->d_var_centering,
                        hier::IntVector(ipi->d_ghost_width));
                  }
               }

               double dmax = -tbox::MathUtilities<double>::getMax();
               double dmin = tbox::MathUtilities<double>::getMax();

               if (data_exists_on_patch) {

                  /*
                   * Scale data (while still double)
                   */
                  const double scale = ipi->d_scale_factor;
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
                  if (scale != 1.0) {
                     for (int i = 0; i < buf_size; ++i) {
                        dbuffer[i] *= scale;
                     }
                  }

                  /*
                   * Determine patch min/max.
                   */

                  for (int i = 0; i < buf_size; ++i) {
                     if (dbuffer[i] > dmax) {
                        dmax = dbuffer[i];
                     }
                     if (dbuffer[i] < dmin) {
                        dmin = dbuffer[i];
                     }
                  }

                  checkFloatMinMax(ipi->d_visit_var_name[depth_id],
                     dmin,
                     dmax,
                     level_number,
                     patch.getLocalId().getValue(),
                     patch_data_id);

                  /*
                   * Convert buffer from double to float
                   */
                  for (int i = 0; i < buf_size; ++i) {
                     fbuffer[i] = static_cast<float>(dbuffer[i]);
                  }

                  /*
                   * Write to disk
                   */
                  std::string vname = ipi->d_visit_var_name[depth_id];
                  patch_HDFGroup->putFloatArray(vname,
                     fbuffer,
                     buf_size);

               } else { // data does not exist on patch

                  dmax = 0.;
                  dmin = 0.;

               }

               /*
                * Write min/max summary info
                */
               if (d_mpi.getRank() == VISIT_MASTER) {
                  int gpn = getGlobalPatchNumber(hierarchy,
                        level_number,
                        patch.getLocalId().getValue());
                  ipi->d_master_min_max[depth_id][gpn].patch_data_on_disk =
                     data_exists_on_patch;
                  ipi->d_master_min_max[depth_id][gpn].min = dmin;
                  ipi->d_master_min_max[depth_id][gpn].max = dmax;
               } else {
                  d_worker_min_max[d_var_id_ctr].patch_data_on_disk =
                     data_exists_on_patch;
                  d_worker_min_max[d_var_id_ctr].min = dmin;
                  d_worker_min_max[d_var_id_ctr].max = dmax;
               }

               /*
                * Increment local var_id counter used for d_mm array.
                */
               ++d_var_id_ctr;

            } // loop over var depths
         } // conventional (not material state) variable
         else { // Data is mixed
            for (int depth_id = 0; depth_id < ipi->d_depth; ++depth_id) {
               std::vector<double> dmix_data;

               /*
                * If its derived data, pack via the derived writer.
                * Otherwise, pack with local private method.
                */
               bool data_exists_on_patch = false;
               int patch_data_id = VISIT_UNDEFINED_INDEX;
               if (ipi->d_is_derived) {

                  // Single function packs clean and mixed data
                  data_exists_on_patch =
                     ipi->d_derived_writer->
                     packMixedDerivedDataIntoDoubleBuffer(
                        dbuffer,
                        dmix_data,
                        patch,
                        pack_box,
                        ipi->d_var_name,
                        depth_id);

               } else {
                  TBOX_ERROR("Mixed Data must be treated as Derived");
               }

               double dmax = -tbox::MathUtilities<double>::getMax();
               double dmin = tbox::MathUtilities<double>::getMax();

               if (data_exists_on_patch) {

                  /*
                   * Scale data (while still double)
                   */
                  const double scale = ipi->d_scale_factor;
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
                  if (scale != 1.0) {
                     for (int i = 0; i < buf_size; ++i) {
                        dbuffer[i] *= scale;
                     }
                     // Scale Mixed data
                     for (std::vector<double>::iterator mi = dmix_data.begin();
                          mi != dmix_data.end(); ++mi) {
                        (*mi) *= scale;
                     }
                  }

                  /*
                   * Determine patch min/max.
                   */

                  int i;
                  for (i = 0; i < buf_size; ++i) {
                     if (dbuffer[i] > dmax) {
                        dmax = dbuffer[i];
                     }
                     if (dbuffer[i] < dmin) {
                        dmin = dbuffer[i];
                     }
                  }
                  // Include mix data in patch min/max
                  for (std::vector<double>::iterator mi = dmix_data.begin();
                       mi != dmix_data.end(); ++mi) {
                     if ((*mi) > dmax) {
                        dmax = (*mi);
                     }
                     if ((*mi) < dmin) {
                        dmin = (*mi);
                     }
                  }

                  checkFloatMinMax(ipi->d_visit_var_name[depth_id],
                     dmin,
                     dmax,
                     level_number,
                     patch.getLocalId().getValue(),
                     patch_data_id);

                  /*
                   * Convert buffer from double to float
                   */
                  for (i = 0; i < buf_size; ++i) {
                     fbuffer[i] = static_cast<float>(dbuffer[i]);
                  }

                  /*
                   * Write to disk
                   */
                  std::string vname = ipi->d_visit_var_name[depth_id];
                  patch_HDFGroup->putFloatArray(vname,
                     fbuffer,
                     buf_size);

                  // If there are no mixed zones in this patch do not write
                  //   mix_zone, mix_mat, vol_fracs, and next_mat
                  int mix_buf_size = static_cast<int>(dmix_data.size());
                  if (mix_buf_size > 0) {
                     // copy mixdata to float for writing
                     // If we had a putFloatVector() this copy could be avoided
                     float* fmixbuffer = new float[mix_buf_size];
                     //std::copy(dmix_data.begin(),dmix_data.end(),fmixbuffer);
                     for (i = 0; i < mix_buf_size; ++i) {
                        fmixbuffer[i] = static_cast<float>(dmix_data[i]);
                     }

                     /*
                      * Write Mixed State Variable to disk
                      * We need to know where to write this
                      */

                     // Use an HDF Group for all material state data
                     std::shared_ptr<tbox::Database> mat_state_HDFGroup;
                     if (patch_HDFGroup->isDatabase("material_state")) {
                        mat_state_HDFGroup =
                           patch_HDFGroup->getDatabase("material_state");
                     } else {
                        mat_state_HDFGroup =
                           patch_HDFGroup->putDatabase("material_state");
                     }
                     mat_state_HDFGroup->putFloatArray(vname,
                        fmixbuffer,
                        mix_buf_size);

                     // For now assume that this buffer cannot be reused.
                     delete[] fmixbuffer;
                  }

               } else { // data does not exist on patch

                  dmax = 0.;
                  dmin = 0.;

               }

               /*
                * Write min/max summary info
                */
               if (d_mpi.getRank() == VISIT_MASTER) {
                  int gpn = getGlobalPatchNumber(hierarchy,
                        level_number,
                        patch.getLocalId().getValue());
                  ipi->d_master_min_max[depth_id][gpn].patch_data_on_disk =
                     data_exists_on_patch;
                  ipi->d_master_min_max[depth_id][gpn].min = dmin;
                  ipi->d_master_min_max[depth_id][gpn].max = dmax;
               } else {
                  d_worker_min_max[d_var_id_ctr].patch_data_on_disk =
                     data_exists_on_patch;
                  d_worker_min_max[d_var_id_ctr].min = dmin;
                  d_worker_min_max[d_var_id_ctr].max = dmax;
               }

               /*
                * Increment local var_id counter used for d_mm array.
                */
               ++d_var_id_ctr;

            } // loop over var depths

         } // material_state variable

         delete[] dbuffer;
         delete[] fbuffer;

      } // var is not species or material
      else {
         for (int depth_id = 0; depth_id < ipi->d_depth; ++depth_id) {
            ++d_var_id_ctr;
         }
      }

   } // iterate over vars

}

/*
 *************************************************************************
 *
 * Private function to pack Material variables into specified HDF
 * database.
 *
 *************************************************************************
 */

void
VisItDataWriter::packMaterialsData(
   const std::shared_ptr<tbox::Database>& patch_HDFGroup,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   hier::Patch& patch)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(level_number >= 0);

   /*
    * Loop over variables and pull out those that are material variables.
    */
   std::shared_ptr<tbox::Database> materials_HDFGroup;
   std::shared_ptr<tbox::Database> material_name_HDFGroup;
   for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
        ipi != d_plot_items.end(); ++ipi) {

      if (ipi->d_isa_material) {

         /*
          * create buffer to hold patch data
          */
         int buf_size = getBufferSize(patch.getBox(),
               hier::IntVector(ipi->d_ghost_width),
               ipi->d_var_centering);

         // Pointers to buffers for dense packing format
         double* dbuffer = 0;  // used to pack var
         float* fbuffer = 0;   // copy to float for writing
         // Pointer to buffer for sparse packing format
         int* ibuffer = 0;     // used to pack mat_list

         // Allocate appropriate memory for packing format
         if (!(ipi->d_is_material_state_variable)) {
            dbuffer = new double[buf_size];
            fbuffer = new float[buf_size];
         } else {
            ibuffer = new int[buf_size];
         }

         for (int depth_id = 0; depth_id < ipi->d_depth; ++depth_id) {

            /*
             * Create the "materials" HDF database entry:
             * <patch HDF group>
             *    "materials"  (only if there are materials)
             *       material_name
             *          "species"    (only if matl has species)
             */
            if (!materials_HDFGroup) {
               materials_HDFGroup =
                  patch_HDFGroup->putDatabase("materials");
            }

            double dmax = -tbox::MathUtilities<double>::getMax();
            double dmin = tbox::MathUtilities<double>::getMax();
            bool data_on_disk = false;

            int return_code;

            /*
             * Method for Sparse representation of materials (and volume
             *   fractions)
             */
            if ((ipi->d_is_material_state_variable)) {
               /*
                * Sparse packing method
                * This requires packing an array of integers (mat_list) with
                *   either the material number of the material occupying the
                *   current cell, or a (negative) index to a set of auxilliary>                 *   vectors.
                *   mat_list:  Material occupying clean zone or (negative)
                *              index to mix_mat, vol_fracs, and next_mat
                *   mix_zones: Cell with which the mixed data is associated
                *   mix_mat:   Material number of partial volumes
                *   vol_fracs: Volume fractions of materials specified by
                *              mix_mat
                *   next_mat:  Next index if there are more materials in mixed>                 *              zone, or zero to indicate mixed zone is complete
                * If the patch is occupied by a single material this can be
                *   indicated by returning
                *   VisMaterialsDataStrategy::VISIT_ALLONE from
                *   packMaterialFractionsIntoSparseBuffers(). This will further
                *   reduce the file size by only writing one entry for mat_list
                *   which indicates the material for the current patch.
                */

               std::vector<int> mix_zones;
               std::vector<int> mix_mat;
               std::vector<double> vol_fracs;
               std::vector<int> next_mat;

               hier::Box pack_box(patch.getBox());
               if (d_write_ghosts) {
                  pack_box.grow(hier::IntVector(ipi->d_ghost_width));
               }

               return_code = d_materials_writer->
                  packMaterialFractionsIntoSparseBuffers(
                     ibuffer, mix_zones, mix_mat, vol_fracs, next_mat,
                     patch, pack_box);

               /*
                * Write to disk
                */
               // Mark material storage type as dense
               std::string mtype = "material_packing_type";
               materials_HDFGroup->putInteger(mtype, 1);
               std::string vname = "mat_list";

               // Is the patch "clean" (ALL_ONE)
               if (return_code == VisMaterialsDataStrategy::VISIT_ALLONE) {
                  materials_HDFGroup->putIntegerArray(vname,
                     ibuffer,
                     1);
               } else { // Otherwise write full mat_list
                  materials_HDFGroup->putIntegerArray(vname,
                     ibuffer,
                     buf_size);
               }
               // The limits should always be 0.0 to 1.0 for volume_fractions
               //   Do this as a check (?)
               dmax = 1.0;
               dmin = 0.0;
               data_on_disk = true;
               // Are there mixed cells?
               if (mix_zones.size() > 0) {
                  TBOX_ASSERT((mix_zones.size() == mix_mat.size()) &&
                     (mix_zones.size() == vol_fracs.size()) &&
                     (mix_zones.size() == next_mat.size()));

                  /*
                   * Determine patch min/max.
                   */
                  for (std::vector<double>::iterator vf = vol_fracs.begin();
                       vf != vol_fracs.end(); ++vf) {
                     if ((*vf) > dmax) {
                        dmax = (*vf);
                     }
                     if ((*vf) < dmin) {
                        dmin = (*vf);
                     }
                  }

                  int dummy_pdata_id = VISIT_UNDEFINED_INDEX;
                  checkFloatMinMax(ipi->d_visit_var_name[depth_id],
                     dmin,
                     dmax,
                     level_number,
                     patch.getLocalId().getValue(),
                     dummy_pdata_id);

                  vname = "mix_zones";
                  materials_HDFGroup->putIntegerVector(vname, mix_zones);

                  vname = "mix_mat";
                  materials_HDFGroup->putIntegerVector(vname, mix_mat);

                  // allocate buffer for volume fraction data
                  float* fmix_data_buffer;
                  fmix_data_buffer = new float[mix_mat.size()];
                  /*
                   * Convert buffer from double to float
                   */
                  for (unsigned int i = 0; i < vol_fracs.size(); ++i) {
                     fmix_data_buffer[i] = static_cast<float>(vol_fracs[i]);
                  }

                  vname = "vol_fracs";
                  materials_HDFGroup->putFloatArray(vname,
                     fmix_data_buffer,
                     static_cast<int>(mix_mat.size()));

                  vname = "next_mat";
                  materials_HDFGroup->putIntegerVector(vname, next_mat);
                  // cleanup buffer
                  delete[] fmix_data_buffer;
               }
            } else { // Legacy material format
                     // Mark material storage type as dense (only once per patch)
               std::string mtype = "material_packing_type";
               if (!(materials_HDFGroup->isInteger(mtype))) {
                  materials_HDFGroup->putInteger(mtype, 0);
               }

               // create materials_name HDF database
               std::string mname = ipi->d_material_name;
               material_name_HDFGroup =
                  materials_HDFGroup->putDatabase(mname);

               // create "species" HDF database for material name
               if (ipi->d_species_names.size() > 0) {
                  ipi->d_species_HDFGroup =
                     material_name_HDFGroup->putDatabase("species");
               }

               hier::Box pack_box(patch.getBox());
               if (d_write_ghosts) {
                  pack_box.grow(hier::IntVector(ipi->d_ghost_width));
               }

               // pack the buffer with material data
               return_code = d_materials_writer->
                  packMaterialFractionsIntoDoubleBuffer(
                     dbuffer,
                     patch,
                     pack_box,
                     ipi->d_material_name);

               // check return code
               if ((return_code != VisMaterialsDataStrategy::VISIT_MIXED)
                   && (return_code !=
                       VisMaterialsDataStrategy::VISIT_ALLONE)
                   && (return_code !=
                       VisMaterialsDataStrategy::VISIT_ALLZERO)) {
                  TBOX_ERROR(
                     "VisItDataWriter::packMaterialsData()"
                     << "\n    Invalid return value from "
                     << "packMaterialFractionsIntoDoubleBuffer()\n");
               }

               if (return_code == VisMaterialsDataStrategy::VISIT_MIXED) {

                  data_on_disk = true;

                  /*
                   * Scale data (while still double)
                   */
                  const double scale = ipi->d_scale_factor;
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
                  if (scale != 1.0) {
                     for (int i = 0; i < buf_size; ++i) {
                        dbuffer[i] *= scale;
                     }
                  }

                  /*
                   * Determine patch min/max.
                   */
                  int i;
                  for (i = 0; i < buf_size; ++i) {
                     if (dbuffer[i] > dmax) {
                        dmax = dbuffer[i];
                     }
                     if (dbuffer[i] < dmin) {
                        dmin = dbuffer[i];
                     }
                  }

                  int dummy_pdata_id = VISIT_UNDEFINED_INDEX;
                  checkFloatMinMax(ipi->d_visit_var_name[depth_id],
                     dmin,
                     dmax,
                     level_number,
                     patch.getLocalId().getValue(),
                     dummy_pdata_id);

                  /*
                   * Convert buffer from double to float
                   */
                  for (i = 0; i < buf_size; ++i) {
                     fbuffer[i] = static_cast<float>(dbuffer[i]);
                  }

                  /*
                   * Write to disk
                   */
                  std::string vname = ipi->d_material_name;
                  vname = vname + "-fractions";
                  material_name_HDFGroup->putFloatArray(vname,
                     fbuffer,
                     buf_size);

               } else if (return_code ==
                          VisMaterialsDataStrategy::VISIT_ALLONE) {

                  data_on_disk = false;
                  dmin = 1.0;
                  dmax = 1.0;

               } else { // return code == VISIT_ALLZERO

                  data_on_disk = false;
                  dmin = 0.0;
                  dmax = 0.0;

               }
            }

            /*
             * Write min/max summary info
             */
            if (d_mpi.getRank() == VISIT_MASTER) {
               int gpn = getGlobalPatchNumber(hierarchy,
                     level_number,
                     patch.getLocalId().getValue());
               ipi->d_master_min_max[depth_id][gpn].patch_data_on_disk =
                  data_on_disk;
               ipi->d_master_min_max[depth_id][gpn].min = dmin;
               ipi->d_master_min_max[depth_id][gpn].max = dmax;
               ipi->d_master_min_max[depth_id][gpn].
               material_composition_code = return_code;
            } else {
               d_worker_min_max[d_var_id_ctr].patch_data_on_disk =
                  data_on_disk;
               d_worker_min_max[d_var_id_ctr].min = dmin;
               d_worker_min_max[d_var_id_ctr].max = dmax;
               d_worker_min_max[d_var_id_ctr].material_composition_code =
                  return_code;
            }

            /*
             * Increment local var_id counter used for d_mm array.
             */
            ++d_var_id_ctr;

         } // loop over var depths

         // Dense packing format delete buffers
         if (!(ipi->d_is_material_state_variable)) {
            delete[] fbuffer;
            delete[] dbuffer;
         } else { // Delete buffer used for sparse format
            delete[] ibuffer;
         }

      } // var is a material
      else {
         for (int depth_id = 0; depth_id < ipi->d_depth; ++depth_id) {
            ++d_var_id_ctr;
         }
      }

   } // iterate over vars

}

/*
 *************************************************************************
 *
 * Private function to pack Species variables into specified HDF
 * database.
 *
 *************************************************************************
 */

void
VisItDataWriter::packSpeciesData(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   hier::Patch& patch)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(level_number >= 0);

   /*
    * Loop over variables and pull out those that are material variables.
    */
   for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
        ipi != d_plot_items.end(); ++ipi) {

      if (ipi->d_isa_species) {

         /*
          * create buffer to hold patch data
          */
         int buf_size = getBufferSize(patch.getBox(),
               hier::IntVector(ipi->d_ghost_width),
               ipi->d_var_centering);

         double* dbuffer = new double[buf_size]; // used to pack var
         float* fbuffer = new float[buf_size]; // copy to float for writing

         hier::Box pack_box(patch.getBox());
         if (d_write_ghosts) {
            pack_box.grow(hier::IntVector(ipi->d_ghost_width));
         }

         for (int depth_id = 0; depth_id < ipi->d_depth; ++depth_id) {

            // pack the buffer with species data
            int return_code = d_materials_writer->
               packSpeciesFractionsIntoDoubleBuffer(
                  dbuffer,
                  patch,
                  pack_box,
                  ipi->d_material_name,
                  ipi->d_species_name);

            // check return code
            if ((return_code != VisMaterialsDataStrategy::VISIT_MIXED)
                && (return_code !=
                    VisMaterialsDataStrategy::VISIT_ALLONE)
                && (return_code ==
                    VisMaterialsDataStrategy::VISIT_ALLZERO)) {
               TBOX_ERROR("VisItDataWriter::packSpeciesData()"
                  << "\n    Invalid return value from "
                  << "packSpeciesFractionsIntoDoubleBuffer()\n"
                  << std::endl);
            }

            double dmax = -tbox::MathUtilities<double>::getMax();
            double dmin = tbox::MathUtilities<double>::getMax();
            bool data_on_disk = false;

            if (return_code == VisMaterialsDataStrategy::VISIT_MIXED) {

               /*
                * Scale data (while still double)
                */
               const double scale = ipi->d_scale_factor;
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
               if (scale != 1.0) {
                  for (int i = 0; i < buf_size; ++i) {
                     dbuffer[i] *= scale;
                  }
               }

               /*
                * Determine patch min/max.
                */
               for (int i = 0; i < buf_size; ++i) {
                  if (dbuffer[i] > dmax) {
                     dmax = dbuffer[i];
                  }
                  if (dbuffer[i] < dmin) {
                     dmin = dbuffer[i];
                  }
               }

               int dummy_pdata_id = VISIT_UNDEFINED_INDEX;
               checkFloatMinMax(ipi->d_visit_var_name[depth_id],
                  dmin,
                  dmax,
                  level_number,
                  patch.getLocalId().getValue(),
                  dummy_pdata_id);

               /*
                * Convert buffer from double to float
                */
               for (int i = 0; i < buf_size; ++i) {
                  fbuffer[i] = static_cast<float>(dbuffer[i]);
               }

               /*
                * Write to disk
                */
               std::string sname = ipi->d_species_name;
               ipi->d_parent_material_pointer->
               d_species_HDFGroup->putFloatArray(sname,
                  fbuffer,
                  buf_size);

            } else if (return_code ==
                       VisMaterialsDataStrategy::VISIT_ALLONE) {

               data_on_disk = false;
               dmin = 1.0;
               dmax = 1.0;

            } else { // return code == VISIT_ALLZERO

               data_on_disk = false;
               dmin = 0.0;
               dmax = 0.0;

            }

            /*
             * Write min/max summary info
             */
            if (d_mpi.getRank() == VISIT_MASTER) {
               int gpn = getGlobalPatchNumber(hierarchy,
                     level_number,
                     patch.getLocalId().getValue());
               ipi->d_master_min_max[depth_id][gpn].patch_data_on_disk =
                  data_on_disk;
               ipi->d_master_min_max[depth_id][gpn].min = dmin;
               ipi->d_master_min_max[depth_id][gpn].max = dmax;
               ipi->d_master_min_max[depth_id][gpn].species_composition_code =
                  return_code;
            } else {
               d_worker_min_max[d_var_id_ctr].patch_data_on_disk =
                  data_on_disk;
               d_worker_min_max[d_var_id_ctr].min = dmin;
               d_worker_min_max[d_var_id_ctr].max = dmax;
               d_worker_min_max[d_var_id_ctr].species_composition_code =
                  return_code;
            }

            /*
             * Increment local var_id counter used for d_mm array.
             */
            ++d_var_id_ctr;

         } // loop over var depths

         delete[] fbuffer;
         delete[] dbuffer;

      } // var is a species
      else {
         for (int depth_id = 0; depth_id < ipi->d_depth; ++depth_id) {
            ++d_var_id_ctr;
         }
      }

   } // iterate over vars

}

/*
 *************************************************************************
 *
 * Private function to check float min/max values.
 *
 *************************************************************************
 */

void
VisItDataWriter::checkFloatMinMax(
   const std::string& var_name,
   const double dmin,
   const double dmax,
   const int level_number,
   const int patch_number,
   const int patch_data_id)
{
   TBOX_ASSERT(level_number >= 0);
   TBOX_ASSERT(patch_number >= 0);
   TBOX_ASSERT(patch_data_id >= -1);

   double fmin = -(tbox::MathUtilities<double>::getMax());
   double fmax = tbox::MathUtilities<double>::getMax();

   if (dmin < fmin) {
      TBOX_ERROR("VisItDataWriter:"
         << "\n    hier::Patch data "
         << var_name
         << " is less than FLT_MIN "
         << "\n    level: " << level_number
         << "  patch: " << patch_number
         << "  patch_data_id: " << patch_data_id
         << "  value: " << dmin
         << "\n    It cannot be read by VisIt."
         << "\n    Make sure data is properly initialized or"
         << "\n    use scale factor to increase its size.");
   }
   if (dmax > fmax) {
      TBOX_ERROR("VisItDataWriter:"
         << "\n    hier::Patch data "
         << var_name
         << " is greater than FLT_MAX "
         << "\n    level: " << level_number
         << "  patch: " << patch_number
         << "  patch_data_id: " << patch_data_id
         << "  value: " << dmax
         << "\n    It cannot be interpreted by VisIt."
         << "\n    Make sure data is properly initialized or"
         << "\n    use scale factor to decrease its size.");
   }
}

/*
 *************************************************************************
 *
 * Private function to write one summary HDF file covering data from all
 * processors for use by VisIt.
 *
 *************************************************************************
 */

void
VisItDataWriter::writeSummaryToHDFFile(
   std::string dump_dirname,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int coarsest_plot_level,
   int finest_plot_level,
   double simulation_time)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(coarsest_plot_level >= 0);
   TBOX_ASSERT(finest_plot_level >= 0);

   int i, ln;

   /*
    * Pack patch min/max information
    */
   exchangeMinMaxPatchInformation(hierarchy,
      coarsest_plot_level,
      finest_plot_level);

   /*
    * The "VISIT_MASTER" writes a set of summary information to
    * the summary file that describes data contained in the visit
    * files written by each MPI process.
    *
    * Although the VISIT_MASTER processor needs the global mesh
    * data, global communication is required, so access the
    * global data to make sure it is built.
    */
   for (ln = 0; ln < hierarchy->getNumberOfLevels(); ++ln) {
      hierarchy->getPatchLevel(ln)->getBoxes();
   }
   int my_proc = d_mpi.getRank();
   if (my_proc == VISIT_MASTER) {
      char temp_buf[VISIT_NAME_BUFSIZE];
      //snprintf(temp_buf, VISIT_NAME_BUFSIZE, "/summary.samrai");
      //string summary_HDFFilename = dump_dirname + temp_buf;
      std::string summary_HDFFilename = dump_dirname + "/" + d_summary_filename;
      std::shared_ptr<tbox::Database> summary_HDFFilePointer(
         std::make_shared<tbox::HDFDatabase>("root"));
      summary_HDFFilePointer->create(summary_HDFFilename);

      /*
       * Create BASIC information HDF Group and provide it the following
       * information:
       *   - VisItWriter version number
       *   - grid_type (CARTESIAN or DEFORMED)
       *   - simulation time
       *   - time step number
       *   - number of processors
       *   - number of file clusters
       *   - DIM
       *   - number of levels
       *   - number of patches at each level (array - int[nlevels])
       *   - total number of patches (on all levels)
       *   - ratio to coarser level (array - int[nlevels][ndim])
       */

      snprintf(temp_buf, VISIT_NAME_BUFSIZE, "BASIC_INFO");
      std::shared_ptr<tbox::Database> basic_HDFGroup(
         summary_HDFFilePointer->putDatabase(std::string(temp_buf)));

      std::shared_ptr<tbox::HDFDatabase> hdf_database(
         SAMRAI_SHARED_PTR_CAST<tbox::HDFDatabase, tbox::Database>(basic_HDFGroup));
      TBOX_ASSERT(hdf_database);
      hid_t basic_group_id = hdf_database->getGroupId();

      std::string key_string = "VDR_version_number";
      basic_HDFGroup->putFloat(key_string, VISIT_DATAWRITER_VERSION_NUMBER);

      key_string = "grid_type";
      std::string data_string;
      if (d_grid_type == VISIT_CARTESIAN) {
         data_string = "CARTESIAN";
      } else if (d_grid_type == VISIT_DEFORMED) {
         data_string = "DEFORMED";
      } else {
         TBOX_ERROR("VisItDataWriter::writeSummaryToHDFFile"
            << "\n    data writer with name " << d_object_name
            << "\n    Illegal grid type: " << d_grid_type << std::endl);
      }
      basic_HDFGroup->putString(key_string, data_string);

      key_string = "time";
      basic_HDFGroup->putDouble(key_string, simulation_time);

      key_string = "time_step_number";
      basic_HDFGroup->putInteger(key_string, d_time_step_number);

      key_string = "number_processors";
      basic_HDFGroup->putInteger(key_string, d_mpi.getSize());

      key_string = "number_file_clusters";
      basic_HDFGroup->putInteger(key_string, d_number_file_clusters);

      key_string = "number_dimensions_of_problem";
      basic_HDFGroup->putInteger(key_string, d_dim.getValue());

      int num_levels;
      unsigned int tot_number_of_patches = 0;

      key_string = "number_levels";
      num_levels = hierarchy->getNumberOfLevels();
      basic_HDFGroup->putInteger(key_string, num_levels);

      key_string = "number_patches_at_level";
      std::vector<int> num_patches_per_level(num_levels);
      for (ln = coarsest_plot_level; ln <= finest_plot_level; ++ln) {
         num_patches_per_level[ln] =
            hierarchy->getPatchLevel(ln)->getGlobalNumberOfPatches();
         tot_number_of_patches += num_patches_per_level[ln];
      }
      basic_HDFGroup->putIntegerVector(key_string, num_patches_per_level);

      key_string = "number_global_patches";
      basic_HDFGroup->putInteger(key_string, tot_number_of_patches);

      /*
       * When writing VisIt data, it expects to see 3D data for
       * xlo, dx, ratios_to_coarser, and number of ghosts.  The
       * VISIT_FIXED_DIM is set to 3, and the third element
       * is zero if we have 2D data.
       */
      key_string = "ratios_to_coarser_levels";
      int idx = 0;
      int* rtcl = new int[num_levels * VISIT_FIXED_DIM];
      for (i = 0; i < num_levels * VISIT_FIXED_DIM; ++i) rtcl[i] = 0;
      for (ln = 0; ln <= finest_plot_level; ++ln) {
         for (i = 0; i < d_dim.getValue(); ++i) {
            idx = ln * VISIT_FIXED_DIM + i;
            if (ln == 0) {
               rtcl[idx] = VISIT_UNDEFINED_INDEX;
            } else {
               rtcl[idx] = d_scaling_ratios[ln](0,i);
            }
         }
      }
      HDFputIntegerArray2D(key_string,
         rtcl,
         num_levels,
         VISIT_FIXED_DIM,
         basic_group_id);
      delete[] rtcl;

      /*
       * Write to BASIC HDF group information about the
       * VisIt plot variables:
       *   - number of visit variables
       *
       *   for each variable {
       *      - name
       *      - centering (1 = CELL, 0 = NODE)
       *      - scale factor
       *      - depth
       *      - ghosts
       *   }
       */
      key_string = "number_visit_variables";
      basic_HDFGroup->putInteger(key_string, d_number_visit_variables);

      std::vector<std::string> var_names(d_number_visit_variables);
      std::vector<int> var_centering(d_number_visit_variables);
      std::vector<double> var_scale_factors(d_number_visit_variables);
      std::vector<int> var_depths(d_number_visit_variables);

      // SGS propose adding array indicating clean/mixed
      std::vector<int> var_material_state_variable(d_number_visit_variables);

      int* var_ghosts = new int[d_number_visit_variables * VISIT_FIXED_DIM];
      for (i = 0; i < d_number_visit_variables * VISIT_FIXED_DIM; ++i) {
         var_ghosts[i] = 0;
      }

      i = 0;
      for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
           ipi != d_plot_items.end(); ++ipi) {

         if (!(ipi->d_isa_material || ipi->d_isa_species)) {
            var_names[i] = ipi->d_var_name;
            if ((ipi->d_var_centering == VISIT_CELL) ||
                (ipi->d_var_centering == VISIT_UNKNOWN_CELL)) {
               var_centering[i] = 1;
            } else {
               var_centering[i] = 0;
            }
            if (ipi->d_is_deformed_coords) {
               var_scale_factors[i] = 0.0;
            } else {
               var_scale_factors[i] = ipi->d_scale_factor;
            }
            var_depths[i] = ipi->d_depth;
            if (!d_write_ghosts) {
               for (int dim = 0; dim < d_dim.getValue(); ++dim) {
                  var_ghosts[i * VISIT_FIXED_DIM + dim] =
                     0;
               }
            } else {
               for (int dim = 0; dim < d_dim.getValue(); ++dim) {
                  var_ghosts[i * VISIT_FIXED_DIM + dim] =
                     ipi->d_ghost_width[dim];
               }
            }
            if (ipi->d_is_material_state_variable) {
               // var_material_state_variable[i] = VISIT_MATERIAL;
               var_material_state_variable[i] = 1;
            } else {
               // var_material_state_variable[i] = VISIT_CLEAN;
               var_material_state_variable[i] = 0;
            }
            ++i;
         }
      }

      key_string = "var_names";
      basic_HDFGroup->putStringVector(key_string, var_names);

      key_string = "var_cell_centered";
      basic_HDFGroup->putIntegerVector(key_string, var_centering);

      key_string = "scaling";
      basic_HDFGroup->putDoubleVector(key_string, var_scale_factors);

      // VCHANGE VisIt needs to read this array so it will
      //         know which variables have mixed material data.
      // SGS
      key_string = "material_state_variable";
      basic_HDFGroup->putIntegerVector(key_string,
         var_material_state_variable);

      if (d_grid_type == VISIT_DEFORMED) {
         std::vector<double> coord_scaling(VISIT_FIXED_DIM);
         for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
              ipi != d_plot_items.end(); ++ipi) {
            if (ipi->d_is_deformed_coords) {
               for (i = 0; i < VISIT_FIXED_DIM; ++i) {
                  coord_scaling[i] = 0.0;
                  if (i < d_dim.getValue()) {
                     coord_scaling[i] = ipi->d_coord_scale_factor[i];
                  }
               }
            }
         }
         key_string = "deformed_coordinate_scaling";
         basic_HDFGroup->putDoubleVector(key_string, coord_scaling);
      }

      key_string = "var_number_components";
      basic_HDFGroup->putIntegerVector(key_string, var_depths);

      key_string = "var_number_ghosts";
      HDFputIntegerArray2D(key_string,
         var_ghosts,
         d_number_visit_variables,
         VISIT_FIXED_DIM,
         basic_group_id);
      delete[] var_ghosts;

      // Embed VisIt expressions
      if (d_visit_expressions.size() > 0) {
         TBOX_ASSERT((d_visit_expressions.size() ==
                      d_visit_expression_keys.size()) &&
            (d_visit_expressions.size() == d_visit_expression_types.size()));

         std::string expdbname = "visit_expressions";
         std::shared_ptr<tbox::Database> expression_HDFGroup(
            summary_HDFFilePointer->putDatabase(expdbname));
         std::string expression_keys("expression_keys");
         std::string expressions("expressions");
         std::string expression_types("expression_types");
         expression_HDFGroup->putStringVector(expression_keys,
            d_visit_expression_keys);
         expression_HDFGroup->putStringVector(expressions,
            d_visit_expressions);
         expression_HDFGroup->putStringVector(expression_types,
            d_visit_expression_types);
      }

      /*
       * Write time and data info to BASIC HDF group
       */
      key_string = "time_of_dump";
/*
 *    std::string temp = dump_dirname + "/temp";
 *    std::ofstream dfile(temp.c_str());
 *    char cmdstr[100];
 *    snprintf(cmdstr, 100, "date > %s",temp.c_str());
 *    system(cmdstr);
 *    std::ifstream date_file(temp.c_str());
 *    std::string ds1, ds2, ds3, ds4, ds5, ds6;
 *    date_file >> ds1 >> ds2 >> ds3 >> ds4 >> ds5 >> ds6;
 *    std::string date = ds1 + " " + ds2 + " " + ds3 + " " + ds4
 + " " + ds5 + " " + ds6;
 *    basic_HDFGroup->putString(key_string, date);
 *    dfile.close();
 *    date_file.close();
 *    snprintf(cmdstr, 100, "rm %s",temp.c_str());
 *    system(cmdstr);
 */
      const int MAXLEN = 256;
      char s[MAXLEN];
      time_t t = time(0);
      strftime(s, MAXLEN, "%a %b %d %H:%M:%S %Z %Y", localtime(&t));
      basic_HDFGroup->putString(key_string, std::string(s));

      /*
       * Write general materials and species information to the
       * materials and species HDF groups:
       *
       *  - material names
       *  - species names
       */
      if (d_materials_names.size() > 0) {
         snprintf(temp_buf, VISIT_NAME_BUFSIZE, "materials");
         std::shared_ptr<tbox::Database> materials_HDFGroup(
            summary_HDFFilePointer->putDatabase(std::string(temp_buf)));

         key_string = "material_names";
         materials_HDFGroup->putStringVector(key_string, d_materials_names);

         int mat_ghosts[VISIT_FIXED_DIM];
         for (i = 0; i < VISIT_FIXED_DIM; ++i) mat_ghosts[i] = 0;
         for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
              ipi != d_plot_items.end(); ++ipi) {
            if (ipi->d_isa_material) {
               if (!d_write_ghosts) {
                  for (int dim = 0; dim < d_dim.getValue(); ++dim) {
                     mat_ghosts[dim] = 0;
                  }
               } else {
                  for (int dim = 0; dim < d_dim.getValue(); ++dim) {
                     mat_ghosts[dim] = ipi->d_ghost_width[dim];
                  }
               }
            }
         }
         key_string = "material_number_ghosts";
         materials_HDFGroup->putIntegerArray(key_string,
            mat_ghosts,
            VISIT_FIXED_DIM);

         snprintf(temp_buf, VISIT_NAME_BUFSIZE, "species");
         std::shared_ptr<tbox::Database> species_HDFGroup(
            materials_HDFGroup->putDatabase(std::string(temp_buf)));

         for (i = 0; i < static_cast<int>(d_materials_names.size()); ++i) {
            key_string = d_materials_names[i];
            for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
                 ipi != d_plot_items.end(); ++ipi) {
               if ((ipi->d_material_name == d_materials_names[i]) &&
                   ipi->d_isa_material) {
                  if (ipi->d_species_names.size() > 0) {
                     species_HDFGroup->putStringVector(key_string,
                        ipi->d_species_names);
                  }
               }
            }
         }
      }

      /*
       * When writing VisIt data, it expects to see 3D data for
       * xlo, dx, ratios_to_coarser, and number of ghosts.  The
       * VISIT_FIXED_DIM is set to 3, and the third element
       * is zero if we have 2D data.
       */
      double geom_lo[VISIT_FIXED_DIM] = { 0., 0., 0. };

      double dx_curr_lev[SAMRAI::MAX_DIM_VAL];
      double patch_xlo, patch_xhi;

      for (i = 0; i < d_dim.getValue(); ++i) {
         dx_curr_lev[i] = 0.0;
      }

      /*
       * Add mesh dx information to BASIC group
       */

      double* dx = new double[VISIT_FIXED_DIM * num_levels];
      double* xlo = new double[VISIT_FIXED_DIM * num_levels];
      for (i = 0; i < VISIT_FIXED_DIM * num_levels; ++i) {
         dx[i] = 0.0;
         xlo[i] = 0.0;
      }
      if (d_grid_type != VISIT_DEFORMED) {
         //This is never entered in multiblock case
         const std::shared_ptr<geom::CartesianGridGeometry> ggeom(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
               hierarchy->getGridGeometry()));
         TBOX_ASSERT(ggeom);
         int next = 0;
         for (ln = coarsest_plot_level; ln <= finest_plot_level; ++ln) {
            for (i = 0; i < VISIT_FIXED_DIM; ++i) {
               if (i < d_dim.getValue()) {
                  if (ln == 0) {
                     xlo[i] = ggeom->getXLower()[i];
                     dx_curr_lev[i] = ggeom->getDx()[i]; // coarsest level dx
                     dx[next] = dx_curr_lev[i];
                  } else {
                     double scale_ratio =
                        static_cast<double>(
                           d_scaling_ratios[ln](0,i));
                     dx_curr_lev[i] = dx_curr_lev[i] / scale_ratio;
                     dx[next] = dx_curr_lev[i];
                  }
               }
               ++next;
            }
         }
      }

      key_string = "dx";
      HDFputDoubleArray2D(key_string,
         dx,
         num_levels,
         VISIT_FIXED_DIM,
         basic_group_id);
      delete[] dx;

      /*
       * Add mesh xlo information to BASIC group
       */
      key_string = "XLO";
      basic_HDFGroup->putDoubleArray(key_string,
         xlo,
         VISIT_FIXED_DIM);
      delete[] xlo;

      /*
       * Write parent/child information
       */
      writeParentChildInfoToSummaryHDFFile(hierarchy, basic_HDFGroup);

      /*
       * Write processor mapping information and domain extents
       * for each patch.
       */

      snprintf(temp_buf, VISIT_NAME_BUFSIZE, "extents");
      std::shared_ptr<tbox::Database> extents_HDFGroup(
         summary_HDFFilePointer->putDatabase(std::string(temp_buf)));
      hdf_database =
         SAMRAI_SHARED_PTR_CAST<tbox::HDFDatabase, tbox::Database>(extents_HDFGroup);
      TBOX_ASSERT(hdf_database);
      hid_t extents_group_id = hdf_database->getGroupId();

      /*
       * Create "patch_map" sub-database of extents group.
       */
      patchMapStruct* pms = new patchMapStruct[tot_number_of_patches];

      for (ln = coarsest_plot_level; ln <= finest_plot_level; ++ln) {
         std::shared_ptr<hier::PatchLevel> patch_level(
            hierarchy->getPatchLevel(ln));
         const std::vector<int>& proc_mapping =
            patch_level->getProcessorMapping().getProcessorMapping();

         for (int p = 0; p < patch_level->getGlobalNumberOfPatches(); ++p) {
            int proc_num = proc_mapping[p];
            int global_patch_id = getGlobalPatchNumber(hierarchy, ln, p);
            pms[global_patch_id].processor_number = proc_num;
            pms[global_patch_id].file_cluster_number =
               d_processor_in_file_cluster_number[proc_num];
            pms[global_patch_id].level_number = ln;
            pms[global_patch_id].patch_number = p;
         }
      }

      key_string = "patch_map";
      HDFputPatchMapStructArray(key_string,
         pms,
         tot_number_of_patches,
         extents_group_id);

      delete[] pms;

      /*
       * Printf-style string on how to name patches:
       * %T --- local patch number
       * %L --- level number
       * %B --- block number
       * %P --- processor number
       * %F --- file number
       */
      key_string = "patch_names_printf";
      std::string patch_names_printf = "f%F R%P B%B L%02L P%05T";
      basic_HDFGroup->putString(key_string, patch_names_printf);

      /*
       * Create "patch_extents" sub-database of extents group.
       */
      patchExtentsStruct* pes = new patchExtentsStruct[tot_number_of_patches];

      for (unsigned int pn = 0; pn < tot_number_of_patches; ++pn) {
         for (i = 0; i < VISIT_FIXED_DIM; ++i) {
            pes[pn].lower[i] = 0;
            pes[pn].upper[i] = 0;
            pes[pn].xlo[i] = 0.;
            pes[pn].xhi[i] = 0.;
         }
      }

      int bdry_type_length = 2*tot_number_of_patches*VISIT_FIXED_DIM;
      std::vector<int> bdry_type(bdry_type_length);
      for (i = 0; i < bdry_type_length; ++i) {
         bdry_type[i] = 0;
      }

      /*
       * Set patch extents
       */
      if (d_grid_type != VISIT_DEFORMED) {
         //This is never entered in multiblock case
         const std::shared_ptr<geom::CartesianGridGeometry> ggeom(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
               hierarchy->getGridGeometry()));
         TBOX_ASSERT(ggeom);
         for (i = 0; i < d_dim.getValue(); ++i) {
            geom_lo[i] = ggeom->getXLower()[i];
            dx_curr_lev[i] = ggeom->getDx()[i]; // coarsest level dx
         }

      } else {

         /*
          * Deformed grid - set extents to 0.
          */
         for (i = 0; i < d_dim.getValue(); ++i) {
            geom_lo[i] = 0.;
            dx_curr_lev[i] = 0.;
         }
      }

      for (ln = coarsest_plot_level; ln <= finest_plot_level; ++ln) {
         const hier::BoxContainer& boxes =
            hierarchy->getPatchLevel(ln)->getBoxes();

         /*
          * Set the dx for the next level
          */
         for (i = 0; i < d_dim.getValue(); ++i) {
            double scale_ratio =
               static_cast<double>(
                  d_scaling_ratios[ln](0,i));
            dx_curr_lev[i] = dx_curr_lev[i] / scale_ratio;
         }

         hier::BoxContainer phys_domain;
         hier::Box phys_domain_box(d_dim); 
         if (d_write_ghosts && hierarchy->getNumberBlocks() == 1) {
             
            if (hierarchy->getGridGeometry()->getDomainIsSingleBox(hier::BlockId(0))) {
               hierarchy->getGridGeometry()->computePhysicalDomain(
                  phys_domain,
                  hierarchy->getPatchLevel(ln)->getRatioToLevelZero(),
                  hier::BlockId(0));
               phys_domain.unorder();
               phys_domain.coalesce();
               TBOX_ASSERT(phys_domain.size() == 1);
               phys_domain_box = phys_domain.front(); 
            }
         }
         int pnum = 0;
         for (hier::BoxContainer::const_iterator itr = boxes.begin();
              itr != boxes.end(); ++itr, ++pnum) {
            int global_patch_id = getGlobalPatchNumber(hierarchy, ln, pnum);
            const hier::Box& box = *itr;
            const int* lower = &box.lower()[0];
            const int* upper = &box.upper()[0];

            for (i = 0; i < d_dim.getValue(); ++i) {
               pes[global_patch_id].lower[i] = lower[i];
               pes[global_patch_id].upper[i] = upper[i];

               if (d_grid_type != VISIT_DEFORMED) {
                  patch_xlo = geom_lo[i]
                     + dx_curr_lev[i] * (double)lower[i];
                  patch_xhi = geom_lo[i]
                     + dx_curr_lev[i] * (double)(upper[i] + 1);
               } else {
                  patch_xlo = 0.;
                  patch_xhi = 0.;
               }
               pes[global_patch_id].xlo[i] = patch_xlo;
               pes[global_patch_id].xhi[i] = patch_xhi;

            }

            if (!phys_domain_box.empty()) {

               int bdry_indx = global_patch_id*2*VISIT_FIXED_DIM;
               for (unsigned short gdim = 0; gdim < d_dim.getValue(); ++gdim) {
                  if (box.lower(gdim) == phys_domain_box.lower(gdim)) {
                     bdry_type[bdry_indx] = 1;
                  }
                  if (box.upper(gdim) == phys_domain_box.upper(gdim)) {
                     bdry_type[bdry_indx+1] = 1;
                  }
                  bdry_indx += 2;
               }
            }

         } // loop over patch boxes
      } // loop over levels

      /*
       * Write patch min/max for each variable.
       */
      std::shared_ptr<tbox::Database> extents_materials_HDFGroup;
      for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
           ipi != d_plot_items.end(); ++ipi) {
         for (int comp = 0; comp < ipi->d_depth; ++comp) {

            /*
             * Regular (i.e. not materials or species) variables
             */
            if (!(ipi->d_isa_material) && !(ipi->d_isa_species)) {

               key_string = ipi->d_visit_var_name[comp] + "-Extents";
               HDFputPatchMinMaxStructArray(
                  key_string,
                  ipi->d_master_min_max[comp],
                  tot_number_of_patches,
                  extents_group_id);

            } else if (ipi->d_isa_material) {

               /*
                * Create materials HDF group
                * set up HDF extents structure for materials:
                *  <extents group>
                *      materials group  (i.e. one constructed below)
                *         material_name group
                *            species group
                */
               if (!extents_materials_HDFGroup) {
                  extents_materials_HDFGroup =
                     extents_HDFGroup->putDatabase("materials");
               }

               key_string = ipi->d_material_name;
               std::shared_ptr<tbox::Database>
               extents_material_name_HDFGroup;
               if (!(ipi->d_is_material_state_variable)) {
                  std::string mname = ipi->d_material_name;
                  // material_name group
                  extents_material_name_HDFGroup =
                     extents_materials_HDFGroup->putDatabase(mname);
                  key_string += "-Fractions";
               } else {
                  // Sparse Format does not need additional group
                  extents_material_name_HDFGroup = extents_materials_HDFGroup;
               }
               std::shared_ptr<tbox::HDFDatabase> extents_database(
                  SAMRAI_SHARED_PTR_CAST<tbox::HDFDatabase, tbox::Database>(
                     extents_material_name_HDFGroup));
               TBOX_ASSERT(extents_database);
               hid_t extents_material_name_group_id =
                  extents_database->getGroupId();

               HDFputPatchMinMaxStructArray(
                  key_string,
                  ipi->d_master_min_max[comp],
                  tot_number_of_patches,
                  extents_material_name_group_id);

               // species group
               if (ipi->d_species_names.size() > 0) {
                  ipi->d_extents_species_HDFGroup =
                     extents_material_name_HDFGroup->putDatabase("species");
               }

            } else if (ipi->d_isa_species) {

               // species
               key_string = ipi->d_species_name;
               std::shared_ptr<tbox::HDFDatabase> extents_database(
                  SAMRAI_SHARED_PTR_CAST<tbox::HDFDatabase, tbox::Database>(
                     ipi->d_parent_material_pointer->d_extents_species_HDFGroup));
               TBOX_ASSERT(extents_database);

               /*
                * species group
                */
               hid_t species_group_id = extents_database->getGroupId();

               HDFputPatchMinMaxStructArray(
                  key_string,
                  ipi->d_master_min_max[comp],
                  tot_number_of_patches,
                  species_group_id);
            }

         } // loop over components
      } // loop over variables

      if (d_worker_min_max) {
         delete[] d_worker_min_max;
      }

      key_string = "patch_extents";
      HDFputPatchExtentsStructArray(key_string,
         pes,
         tot_number_of_patches,
         extents_group_id);

      delete[] pes;

      key_string = "bdry_type";
      HDFputBoundaryTypeArray(key_string,
         bdry_type,
         tot_number_of_patches,
         extents_group_id);

      summary_HDFFilePointer->close();

   } // if VISIT_MASTER

   tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier();

   if (my_proc == VISIT_MASTER) {

      /*
       * Add this dump entry to dumps.visit file
       */
      if (d_time_step_number == 0) s_summary_file_opened = false;
      std::string path;
      if (!d_top_level_directory_name.empty() &&
          d_top_level_directory_name[d_top_level_directory_name.length() - 1] == '/') {
         path = d_top_level_directory_name + "dumps.visit";
      } else {
         path = d_top_level_directory_name + "/dumps.visit";
      }
      std::string file = d_current_dump_directory_name + "/"
         + d_summary_filename;

      /*
       * If summary file has not yet been opened, open and write file.
       * If it has been opened, append this to file.
       */
      if (!s_summary_file_opened) {
         s_summary_file_opened = true;
         std::ofstream sfile(path.c_str(), std::ios::out);
         sfile << file << "\n";
         sfile.close();
      } else {
         std::ofstream sfile(path.c_str(), std::ios::app);
         sfile << file << "\n";
         sfile.close();
      }
   }

   tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier();
}

/*
 *************************************************************************
 *
 * Private function to store min/max information on each patch for each
 * variable.  The "master" processor allocates an array that will hold
 * the global data.  The "worker" processors send this information to
 * the master which in turn unpacks and stores the data.
 *
 *************************************************************************
 */
void
VisItDataWriter::exchangeMinMaxPatchInformation(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarsest_plot_level,
   const int finest_plot_level)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(coarsest_plot_level >= 0);
   TBOX_ASSERT(finest_plot_level >= 0);

   /*
    * Compute max number of patches on any processor, and the total number of
    * patches in the problem.
    */
   int ln, pn, comp, item_ctr;
   int number_local_patches = 0;

   for (ln = coarsest_plot_level; ln <= finest_plot_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> patch_level(
         hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(patch_level->begin());
           ip != patch_level->end(); ++ip) {
         ++number_local_patches;
      }
   }

   int max_number_local_patches = number_local_patches;
   if (d_mpi.getSize() > 1) {
      d_mpi.AllReduce(&max_number_local_patches, 1, MPI_MAX);
   }

   int num_items_to_plot = d_number_visit_variables_plus_depth
      + static_cast<int>(d_materials_names.size()) // number materials
      + d_number_species;

   int message_size = max_number_local_patches
      * static_cast<int>(sizeof(patchMinMaxStruct)) * num_items_to_plot;

   if (d_mpi.getRank() != VISIT_MASTER) {

      /*
       * Worker processor:  send contents of d_worker_min_max array that
       * was setup in "initializePlotVariableMinMaxInfo()" to the
       * master processor.
       */
      if (number_local_patches > 0) {
         if (d_mpi.getSize() > 1) {
            d_mpi.Send(d_worker_min_max,
               message_size,
               MPI_BYTE,
               VISIT_MASTER,
               0);
         }
      }

   } else { // (d_mpi.getRank() == VISIT_MASTER)

      /*
       * Master processor:  Receive the min/max information sent by the
       * worker processors (above).  Unpack into the local d_mm array
       * for each plot variable.
       */

      // recv buffer large enough to receive info from any processor.
      patchMinMaxStruct* buf = 0;
      if (d_number_working_slaves > 0) {
         buf = new patchMinMaxStruct[max_number_local_patches
                                     * num_items_to_plot];
         memset((char *)buf, 0,
            max_number_local_patches * num_items_to_plot
            * sizeof(patchMinMaxStruct));
      }

      /*
       * Receive information sent by "sending_proc".
       */
      int number_msgs_recvd = 0;
      while (number_msgs_recvd < d_number_working_slaves) {
         int sending_proc = -1;
         if (d_mpi.getSize() > 1) {
            tbox::SAMRAI_MPI::Status status;
            d_mpi.Recv(buf,
               message_size,
               MPI_BYTE,
               MPI_ANY_SOURCE,
               MPI_ANY_TAG,
               &status);
            sending_proc = status.MPI_SOURCE;
         }
         ++number_msgs_recvd;

         /*
          * Unpack the information from buf and fill d_worker_min_max array.
          */
         item_ctr = 0;

         for (ln = coarsest_plot_level; ln <= finest_plot_level; ++ln) {
            const std::vector<int>& proc_mapping =
               hierarchy->getPatchLevel(ln)->getProcessorMapping().getProcessorMapping();

            int npatches_on_level = static_cast<int>(proc_mapping.size());
            for (pn = 0; pn < npatches_on_level; ++pn) {
               if (proc_mapping[pn] == sending_proc) {
                  int global_patch_id =
                     getGlobalPatchNumber(hierarchy, ln, pn);
                  for (std::list<VisItItem>::iterator ipi(d_plot_items.begin());
                       ipi != d_plot_items.end(); ++ipi) {
                     for (comp = 0; comp < ipi->d_depth; ++comp) {
                        ipi->d_master_min_max[comp][global_patch_id] =
                           buf[item_ctr];
                        ++item_ctr;
                     }
                  }  // variables
               } // patch from sending proc?
            } // patches
         } // levels
      } // while msgs_recvd < working procs

      if (d_number_working_slaves > 0) {
         delete[] buf;
      }

   } // proc == VISIT_MASTER?

}

/*
 *************************************************************************
 *
 * Private function to find and write to summary file parent & child
 * info.  For each global patch number, find its children using the
 * box_tree.  Record the children, as well as each child's parent, in a
 * child_parent array. A child_ptrs array records for each global patch
 * number, the number of children that patch has, as well as the offset
 * into the child_parent array where the patch numbers of those children
 * are stored.  If a patch has no children, offset = -1. Next, the child
 * info from the child_parent array is copied into the final child array.
 * Then the child_parent array is sorted by child number.  Now all
 * parents of a given patch are grouped together in this sorted array.
 * The parents are then stored in a parent array, and a parent_ptrs
 * array is created similar to the child_ptrs array.
 *
 *************************************************************************
 */
void
VisItDataWriter::writeParentChildInfoToSummaryHDFFile(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const std::shared_ptr<tbox::Database>& basic_HDFGroup)
{
   TBOX_ASSERT(hierarchy);

   struct cpPointerStruct {  // auxiliary info for child/parent data struct
      int offset;
      union {
         int number_children;
         int number_parents;
      } u;
   };

   /*
    * Find child patches of each global patch number
    */
   int tot_number_of_patches = 0;
   int finest_level = hierarchy->getFinestLevelNumber();
   int ln;

   for (ln = 0; ln <= finest_level; ++ln) {
      tot_number_of_patches += hierarchy->getPatchLevel(ln)->getGlobalNumberOfPatches();
   }
   int chunk_size = 2 * tot_number_of_patches;
   int current_child_parent_max_size = chunk_size;
   struct cpPointerStruct* child_ptrs =
      new struct cpPointerStruct[tot_number_of_patches];
   struct childParentStruct* child_parent =
      new struct childParentStruct[chunk_size];
   int child_parent_idx = 0;
   int child_ptrs_idx = 0;

   for (ln = 0; ln <= finest_level; ++ln) {
      const hier::BoxContainer& coarser_boxes =
         hierarchy->getPatchLevel(ln)->getBoxLevel()->getGlobalizedVersion().getGlobalBoxes();

      std::shared_ptr<hier::BoxContainer> child_box_tree;
      hier::IntVector ratio(hier::IntVector::getZero(d_dim));

      if (ln != finest_level) {
         std::shared_ptr<hier::PatchLevel> child_patch_level(
            hierarchy->getPatchLevel(ln + 1));
         ratio = child_patch_level->getRatioToCoarserLevel();

         const hier::BoxContainer& global_child_boxes =
            child_patch_level->getBoxLevel()->
            getGlobalizedVersion().getGlobalBoxes();

         /*
          * We need to strip out periodic Boxes.  This is only
          * necessary in single block case, as multiblock case can never
          * have periodic conditions.
          */
         if (hierarchy->getGridGeometry()->getNumberBlocks() == 1) {

            hier::BoxContainer non_per_child_boxes;
            for (hier::RealBoxConstIterator gi(global_child_boxes.realBegin());
                 gi != global_child_boxes.realEnd(); ++gi) {
               non_per_child_boxes.insert(*gi);
            }

            child_box_tree.reset(
               new hier::BoxContainer(
                  non_per_child_boxes));
            child_box_tree->makeTree(hierarchy->getGridGeometry().get());

         } else {

            child_box_tree.reset(
               new hier::BoxContainer(
                  global_child_boxes));
            child_box_tree->makeTree(hierarchy->getGridGeometry().get());

         }
      }

      for (hier::RealBoxConstIterator bi(coarser_boxes.realBegin());
           bi != coarser_boxes.realEnd(); ++bi) {

         if (ln == finest_level) {
            child_ptrs[child_ptrs_idx].u.number_children = 0;
            child_ptrs[child_ptrs_idx++].offset = VISIT_UNDEFINED_INDEX;
         } else {
            hier::BlockId block_id(bi->getBlockId());
            hier::Box compare_box(*bi);
            compare_box.refine(ratio);

            hier::BoxContainer overlap_boxes;

            child_box_tree->findOverlapBoxes(
               overlap_boxes,
               compare_box,
               ratio);

            int num_kids = static_cast<int>(overlap_boxes.size());
            child_ptrs[child_ptrs_idx].u.number_children = num_kids;

            if (num_kids == 0) {
               child_ptrs[child_ptrs_idx++].offset = VISIT_UNDEFINED_INDEX;
            } else {
               child_ptrs[child_ptrs_idx++].offset = child_parent_idx;
               if ((child_parent_idx + num_kids) >
                   current_child_parent_max_size) {
                  current_child_parent_max_size += chunk_size;
                  struct childParentStruct* temp = child_parent;
                  child_parent =
                     new struct
                     childParentStruct[current_child_parent_max_size];
                  for (int idx = 0; idx < child_parent_idx; ++idx) {
                     child_parent[idx].child = temp[idx].child;
                     child_parent[idx].parent = temp[idx].parent;
                  }
                  delete[] temp;
               }

               for (hier::BoxContainer::iterator ob_itr = overlap_boxes.begin();
                    ob_itr != overlap_boxes.end(); ++ob_itr) {
                  child_parent[child_parent_idx].child =
                     getGlobalPatchNumber(hierarchy, ln + 1,
                        ob_itr->getLocalId().getValue());
                  child_parent[child_parent_idx++].parent =
                     getGlobalPatchNumber(hierarchy, ln,
                        bi->getLocalId().getValue());
               }
            }
         }
      }
   }

   int* parent_array = 0;
   int* child_array = 0;
   int parent_array_length = 0;
   int child_array_length = child_parent_idx;

   struct cpPointerStruct* parent_ptrs = 0;

   // copy child info to child array
   if (child_array_length > 0) {

      child_array = new int[child_array_length];
      for (int idx = 0; idx < child_array_length; ++idx) {
         child_array[idx] = child_parent[idx].child;
      }

      qsort((char *)child_parent,
         child_array_length,
         sizeof(struct childParentStruct),
         &childParentCompareFunc);

      // now record parents in the parent array
      parent_array = new int[chunk_size];
      int parent_size = chunk_size;
      int cp_idx = 0;
      int next_parent = 0;
      parent_ptrs =
         new struct cpPointerStruct[tot_number_of_patches];

      for (int gpn = 0; gpn < tot_number_of_patches; ++gpn) {
         if (gpn < child_parent[cp_idx].child) {
            parent_ptrs[gpn].offset = VISIT_UNDEFINED_INDEX;
            parent_ptrs[gpn].u.number_parents = 0;
         } else {
            int num_pars = 0;
            while (cp_idx < current_child_parent_max_size &&
                   child_parent[cp_idx].child == gpn) {
               if (num_pars == 0) {
                  parent_ptrs[gpn].offset = next_parent;
               }
               ++num_pars;
               if (next_parent >= parent_size) {
                  // increase size of parent_array
                  int old_parent_size = parent_size;
                  int* temp = new int[old_parent_size];
                  for (int i = 0; i < old_parent_size; ++i) {
                     temp[i] = parent_array[i];
                  }
                  delete[] parent_array;
                  parent_size += chunk_size;
                  parent_array = new int[parent_size];
                  for (int i = 0; i < old_parent_size; ++i) {
                     parent_array[i] = temp[i];
                  }
                  delete[] temp;
               }
               parent_array[next_parent++] = child_parent[cp_idx++].parent;
            }
            parent_ptrs[gpn].u.number_parents = num_pars;
         }
      }
      parent_array_length = next_parent;
   }
   delete[] child_parent;

   /*
    * Write parent & child info to summary file
    */
   std::string key_string = "child_array_length";
   basic_HDFGroup->putInteger(key_string, child_array_length);
   key_string = "parent_array_length";
   basic_HDFGroup->putInteger(key_string, parent_array_length);
   std::shared_ptr<tbox::HDFDatabase> hdf_database(
      SAMRAI_SHARED_PTR_CAST<tbox::HDFDatabase, tbox::Database>(basic_HDFGroup));
   TBOX_ASSERT(hdf_database);
   hid_t basic_group_id = hdf_database->getGroupId();
   if (child_array_length > 0) {
      key_string = "child_array";
      basic_HDFGroup->putIntegerArray(key_string,
         child_array,
         child_array_length);

      key_string = "child_pointer_array";
      HDFputChildParentStructArray(
         key_string,
         (void *)child_ptrs,
         tot_number_of_patches,
         basic_group_id,
         sizeof(cpPointerStruct),
         "number_children");
   }
   if (child_array) delete[] child_array;
   if (child_ptrs) delete[] child_ptrs;

   if (parent_array_length > 0) {
      key_string = "parent_array";
      basic_HDFGroup->putIntegerArray(key_string,
         parent_array,
         parent_array_length);
      key_string = "parent_pointer_array";
      HDFputChildParentStructArray(
         key_string,
         (void *)parent_ptrs,
         tot_number_of_patches,
         basic_group_id,
         sizeof(cpPointerStruct),
         "number_parents");
   }

   if (parent_array) {
      delete[] parent_array;
   }
   if (parent_ptrs) {
      delete[] parent_ptrs;
   }

}

/*
 *************************************************************************
 *
 *    childParentCompareFunc() used by qsort to sort child_parent array
 *        by child patch num to find all parents of a given child
 *
 *************************************************************************
 */

int
VisItDataWriter::childParentCompareFunc(
   const void* s1,
   const void* s2)
{
   struct childParentStruct* a = (struct childParentStruct *)s1;
   struct childParentStruct* b = (struct childParentStruct *)s2;
   if (((struct childParentStruct *)a)->child >
       ((struct childParentStruct *)b)->child) {
      return 1;
   } else if (((struct childParentStruct *)b)->child >
              ((struct childParentStruct *)a)->child) {
      return VISIT_UNDEFINED_INDEX;
   } else {
      return 0;
   }
}

/*
 *************************************************************************
 *
 * Private utility function to pack DIM patch data into 1D double
 * precision buffer, omitting ghost data if necessary. Data is packed in
 * column major order, i.e. x0,y0,z0, x1,y0,z0, x2,y0,z0, ...
 *
 *************************************************************************
 */

void
VisItDataWriter::packPatchDataIntoDoubleBuffer(
   const std::shared_ptr<hier::PatchData>& pdata,
   const int depth_index,
   const variable_data_type data_type,
   const hier::Box patch_box,
   double* buffer,
   const variable_centering centering,
   const hier::IntVector& ghost_width)
{
   TBOX_ASSERT(depth_index >= 0);
   TBOX_ASSERT((patch_box * pdata->getGhostBox()).isSpatiallyEqual(patch_box));

   /*
    * Currently, Fortran calls in this method limit it to two or
    * three dimensions only.  It will not work for DIM < 2 or
    * DIM > 3.  Give the user a run-time error asserting this.
    */
   if (d_dim < tbox::Dimension(2) || d_dim > tbox::Dimension(3)) {
      TBOX_ERROR(
         d_object_name << "VisItDataWriter::packPatchDataIntoDoubleBuffer()"
                       << "\n  This case has DIM = " << d_dim
                       << "\n  Dimensions < 2 or > 3 are not supported at"
                       << "\n  this time." << std::endl);
   }

   int buf_size = getBufferSize(patch_box,
         ghost_width,
         centering);

   hier::Index databox_lower = pdata->getGhostBox().lower();
   hier::Index databox_upper = pdata->getGhostBox().upper();

   hier::Box plot_box(patch_box);
   plot_box.grow(ghost_width);

   hier::Index plolower = plot_box.lower();
   hier::Index ploupper = plot_box.upper();

   /*
    * Extend index boundaries by one point if NODE data is used.
    */
   if (centering == VISIT_NODE || centering == VISIT_UNKNOWN_NODE) {
      databox_upper += 1;
      ploupper += 1;
   }

   switch (data_type) {

      case VISIT_FLOAT: {
         const float* dat_ptr = 0;
         if (centering == VISIT_CELL) {
            std::shared_ptr<const pdat::CellData<float> > fpdata(
               SAMRAI_SHARED_PTR_CAST<const pdat::CellData<float>, hier::PatchData>(pdata));

            TBOX_ASSERT(fpdata);

            dat_ptr = fpdata->getPointer(depth_index);
         } else if (centering == VISIT_UNKNOWN_CELL) {
            pdat::CellData<float> cell_copy(plot_box,
                                            1,
                                            hier::IntVector::getZero(d_dim));
            pdata->copy2(cell_copy);
            dat_ptr = cell_copy.getPointer();
         } else if (centering == VISIT_NODE) {
            std::shared_ptr<const pdat::NodeData<float> > fpdata(
               SAMRAI_SHARED_PTR_CAST<const pdat::NodeData<float>, hier::PatchData>(pdata));

            TBOX_ASSERT(fpdata);

            dat_ptr = fpdata->getPointer(depth_index);
         } else if (centering == VISIT_UNKNOWN_NODE) {
            pdat::NodeData<float> node_copy(plot_box,
                                            1,
                                            hier::IntVector::getZero(d_dim));
            pdata->copy2(node_copy);
            dat_ptr = node_copy.getPointer();
         }
         if (d_dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(cpfdat2buf2d, CPFDAT2BUF2D) (databox_lower(0),
               databox_lower(1),
               plolower(0), plolower(1),
               ploupper(0), ploupper(1),
               databox_upper(0), databox_upper(1),
               dat_ptr, buffer, buf_size);
         } else if (d_dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(cpfdat2buf3d, CPFDAT2BUF3D) (databox_lower(0),
               databox_lower(1), databox_lower(2),
               plolower(0), plolower(1), plolower(2),
               ploupper(0), ploupper(1), ploupper(2),
               databox_upper(0), databox_upper(1), databox_upper(2),
               dat_ptr, buffer, buf_size);
         }
         break;
      }

      case VISIT_DOUBLE: {
         const double* dat_ptr = 0;
         if (centering == VISIT_CELL) {
            std::shared_ptr<const pdat::CellData<double> > dpdata(
               SAMRAI_SHARED_PTR_CAST<const pdat::CellData<double>, hier::PatchData>(pdata));
            TBOX_ASSERT(dpdata);

            dat_ptr = dpdata->getPointer(depth_index);
         } else if (centering == VISIT_UNKNOWN_CELL) {
            pdat::CellData<double> cell_copy(plot_box,
                                             1,
                                             hier::IntVector::getZero(d_dim));
            pdata->copy2(cell_copy);
            dat_ptr = cell_copy.getPointer();
         } else if (centering == VISIT_NODE) {
            std::shared_ptr<const pdat::NodeData<double> > dpdata(
               SAMRAI_SHARED_PTR_CAST<const pdat::NodeData<double>, hier::PatchData>(pdata));
            TBOX_ASSERT(dpdata);

            dat_ptr = dpdata->getPointer(depth_index);
         } else if (centering == VISIT_UNKNOWN_NODE) {
            pdat::NodeData<double> node_copy(plot_box,
                                             1,
                                             hier::IntVector::getZero(d_dim));
            pdata->copy2(node_copy);
            dat_ptr = node_copy.getPointer();
         }
         if (d_dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(cpddat2buf2d, CPDDAT2BUF2D) (databox_lower(0),
               databox_lower(1),
               plolower(0), plolower(1),
               ploupper(0), ploupper(1),
               databox_upper(0), databox_upper(1),
               dat_ptr, buffer, buf_size);
         } else if (d_dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(cpddat2buf3d, CPDDAT2BUF3D) (databox_lower(0),
               databox_lower(1), databox_lower(2),
               plolower(0), plolower(1), plolower(2),
               ploupper(0), ploupper(1), ploupper(2),
               databox_upper(0), databox_upper(1), databox_upper(2),
               dat_ptr, buffer, buf_size);
         }
         break;
      }

      case VISIT_INT: {
         const int* dat_ptr = 0;
         if (centering == VISIT_CELL) {
            std::shared_ptr<const pdat::CellData<int> > ipdata(
               SAMRAI_SHARED_PTR_CAST<const pdat::CellData<int>, hier::PatchData>(pdata));

            TBOX_ASSERT(ipdata);

            dat_ptr = ipdata->getPointer(depth_index);
         } else if (centering == VISIT_UNKNOWN_CELL) {
            pdat::CellData<int> cell_copy(plot_box, 1, hier::IntVector::getZero(
                                             d_dim));
            pdata->copy2(cell_copy);
            dat_ptr = cell_copy.getPointer();
         } else if (centering == VISIT_NODE) {
            std::shared_ptr<const pdat::NodeData<int> > ipdata(
               SAMRAI_SHARED_PTR_CAST<const pdat::NodeData<int>, hier::PatchData>(pdata));

            TBOX_ASSERT(ipdata);

            dat_ptr = ipdata->getPointer(depth_index);
         } else if (centering == VISIT_UNKNOWN_NODE) {
            pdat::NodeData<int> node_copy(plot_box, 1, hier::IntVector::getZero(
                                             d_dim));
            pdata->copy2(node_copy);
            dat_ptr = node_copy.getPointer();
         }
         if (d_dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(cpidat2buf2d, CPIDAT2BUF2D) (databox_lower(0),
               databox_lower(1),
               plolower(0), plolower(1),
               ploupper(0), ploupper(1),
               databox_upper(0), databox_upper(1),
               dat_ptr, buffer, buf_size);
         } else if (d_dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(cpidat2buf3d, CPIDAT2BUF3D) (databox_lower(0),
               databox_lower(1), databox_lower(2),
               plolower(0), plolower(1), plolower(2),
               ploupper(0), ploupper(1), ploupper(2),
               databox_upper(0), databox_upper(1), databox_upper(2),
               dat_ptr, buffer, buf_size);
         }
         break;
      }

      default: {
         TBOX_ERROR(
            d_object_name << ":packPatchDataIntoDoubleBuffer()"
                          << "\n  Unknown type.  ***Exiting."
                          << std::endl);
      }
   }
}

/*
 *************************************************************************
 *
 * Create a 2D integer array entry in an HDF database with the specified
 * key name.  The array type is based on the hdf type H5T_NATIVE_INT.
 *
 *************************************************************************
 */

void
VisItDataWriter::HDFputIntegerArray2D(
   const std::string& key,
   const int* data,
   const int nelements0,
   const int nelements1,
   const hid_t group_id)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);
   TBOX_ASSERT((nelements0 > 0) && (nelements1 > 0));

   herr_t errf = 0;
   if ((nelements0 > 0) && (nelements1 > 0)) {
      hsize_t dim[] = { static_cast<hsize_t>(nelements0),
                        static_cast<hsize_t>(nelements1) };
      hid_t space = H5Screate_simple(2, dim, 0);

      TBOX_ASSERT(space >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(group_id,
            key.c_str(),
            H5T_NATIVE_INT,
            space,
            H5P_DEFAULT,
            H5P_DEFAULT,
            H5P_DEFAULT);
#else
      hid_t dataset = H5Dcreate(group_id,
            key.c_str(),
            H5T_NATIVE_INT,
            space,
            H5P_DEFAULT);
#endif

      TBOX_ASSERT(dataset >= 0);

      errf = H5Dwrite(dataset,
            H5T_NATIVE_INT,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            data);
      TBOX_ASSERT(errf >= 0);

      errf = H5Sclose(space);

      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(dataset);
      TBOX_ASSERT(errf >= 0);

   } else {
      TBOX_ERROR("VisItDataWriter::HDFputIntegerArray2D()"
         << "\n    data writer with name " << d_object_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
   NULL_USE(errf);
}

/*
 *************************************************************************
 *
 * Create a 2D double array entry in an HDF database with the specified
 * key name.  The array type is based on the hdf type H5T_NATIVE_DOUBLE.
 *
 *************************************************************************
 */

void
VisItDataWriter::HDFputDoubleArray2D(
   const std::string& key,
   const double* data,
   const int nelements0,
   const int nelements1,
   const hid_t group_id)
{

   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);
   TBOX_ASSERT((nelements0 > 0) && (nelements1 > 0));

   herr_t errf = 0;
   if ((nelements0 > 0) && (nelements1 > 0)) {
      hsize_t dim[] = { static_cast<hsize_t>(nelements0),
                        static_cast<hsize_t>(nelements1) };
      hid_t space = H5Screate_simple(2, dim, 0);

      TBOX_ASSERT(space >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(group_id,
            key.c_str(),
            H5T_NATIVE_DOUBLE,
            space,
            H5P_DEFAULT,
            H5P_DEFAULT,
            H5P_DEFAULT);
#else
      hid_t dataset = H5Dcreate(group_id,
            key.c_str(),
            H5T_NATIVE_DOUBLE,
            space,
            H5P_DEFAULT);
#endif

      TBOX_ASSERT(dataset >= 0);

      errf = H5Dwrite(dataset,
            H5T_NATIVE_DOUBLE,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            data);

      TBOX_ASSERT(errf >= 0);

      errf = H5Sclose(space);
      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(dataset);
      TBOX_ASSERT(errf >= 0);

   } else {
      TBOX_ERROR("VisItDataWriter::HDFputDoubleArray2D()"
         << "\n    data writer with name " << d_object_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
   NULL_USE(errf);
}

/*
 *************************************************************************
 *
 * Create an array of patch extent (pe) structs in an HDF database
 * with the specified key name.
 *
 *************************************************************************
 */

void
VisItDataWriter::HDFputPatchExtentsStructArray(
   const std::string& key,
   const patchExtentsStruct* data,
   const int nelements,
   const hid_t group_id)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);
   TBOX_ASSERT(nelements > 0);

   herr_t errf = 0;
   if (nelements > 0) {
      hid_t space;
      hsize_t dim[1];
      dim[0] = nelements;
      space = H5Screate_simple(1, dim, 0);
      TBOX_ASSERT(space >= 0);

      hid_t pe_id = H5Tcreate(H5T_COMPOUND, sizeof(patchExtentsStruct));
      TBOX_ASSERT(pe_id >= 0);

      hsize_t dim1[1];
      dim1[0] = VISIT_FIXED_DIM;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t intXdType = H5Tarray_create(H5T_NATIVE_INT, 1, dim1);
#else
      hid_t intXdType = H5Tarray_create(H5T_NATIVE_INT, 1, dim1, 0);
#endif
      TBOX_ASSERT(intXdType >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t doubleXdType = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dim1);
#else
      hid_t doubleXdType = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dim1, 0);
#endif
      TBOX_ASSERT(doubleXdType >= 0);

      errf = H5Tinsert(pe_id,
            "lower",
            HOFFSET(patchExtentsStruct, lower),
            intXdType);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tinsert(pe_id,
            "upper",
            HOFFSET(patchExtentsStruct, upper),
            intXdType);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tinsert(pe_id,
            "xlo",
            HOFFSET(patchExtentsStruct, xlo),
            doubleXdType);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tinsert(pe_id,
            "xup",
            HOFFSET(patchExtentsStruct, xhi),
            doubleXdType);
      TBOX_ASSERT(errf >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(group_id,
            key.c_str(),
            pe_id,
            space,
            H5P_DEFAULT,
            H5P_DEFAULT,
            H5P_DEFAULT);
#else
      hid_t dataset = H5Dcreate(group_id,
            key.c_str(),
            pe_id,
            space,
            H5P_DEFAULT);
#endif
      TBOX_ASSERT(dataset >= 0);

      errf = H5Dwrite(dataset,
            pe_id,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            data);
      TBOX_ASSERT(errf >= 0);

      errf = H5Sclose(space);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tclose(pe_id);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tclose(intXdType);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tclose(doubleXdType);
      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(dataset);
      TBOX_ASSERT(errf >= 0);

   } else {
      TBOX_ERROR("VisItDataWriter::HDFputPatchExtentsStructArray()"
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
   NULL_USE(errf);
}


void
VisItDataWriter::HDFputBoundaryTypeArray(
   const std::string& key,
   const std::vector<int>& data,
   const int num_patches,
   const hid_t group_id)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(num_patches > 0);
   TBOX_ASSERT(static_cast<size_t>(2*num_patches*VISIT_FIXED_DIM) == data.size());

   herr_t errf = 0;
   if (num_patches > 0) {
      hid_t space;
      hsize_t dim[1];
      dim[0] = 2*num_patches*VISIT_FIXED_DIM;
      space = H5Screate_simple(1, dim, 0);
      TBOX_ASSERT(space >= 0);

      hid_t int_id = H5Tcopy(H5T_NATIVE_INT);
      TBOX_ASSERT(int_id >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(group_id,
            key.c_str(),
            int_id,
            space,
            H5P_DEFAULT,
            H5P_DEFAULT,
            H5P_DEFAULT);
#else
      hid_t dataset = H5Dcreate(group_id,
            key.c_str(),
            int_id,
            space,
            H5P_DEFAULT);
#endif
      TBOX_ASSERT(dataset >= 0);

      errf = H5Dwrite(dataset,
            int_id,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            &data[0]);
      TBOX_ASSERT(errf >= 0);

      errf = H5Sclose(space);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tclose(int_id);
      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(dataset);
      TBOX_ASSERT(errf >= 0);

   } else {
      TBOX_ERROR("VisItDataWriter::HDFputBoundaryTypeArray()"
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
   NULL_USE(errf);
}

/*
 *************************************************************************
 *
 * Create patch map for each of the patch extents HDF entries
 * with the specified key name.
 *
 *************************************************************************
 */

void
VisItDataWriter::HDFputPatchMapStructArray(
   const std::string& key,
   const patchMapStruct* data,
   const int nelements,
   const hid_t group_id)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);
   TBOX_ASSERT(nelements > 0);

   herr_t errf = 0;
   if (nelements > 0) {
      hid_t space;
      hsize_t dim[1];
      dim[0] = nelements;
      space = H5Screate_simple(1, dim, 0);
      TBOX_ASSERT(space >= 0);

      hid_t pm_id = H5Tcreate(H5T_COMPOUND, sizeof(patchMapStruct));
      TBOX_ASSERT(pm_id >= 0);

      errf = H5Tinsert(pm_id,
            "processor_number",
            HOFFSET(patchMapStruct, processor_number),
            H5T_NATIVE_INT);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tinsert(pm_id,
            "file_cluster_number",
            HOFFSET(patchMapStruct, file_cluster_number),
            H5T_NATIVE_INT);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tinsert(pm_id,
            "level_number",
            HOFFSET(patchMapStruct, level_number),
            H5T_NATIVE_INT);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tinsert(pm_id,
            "patch_number",
            HOFFSET(patchMapStruct, patch_number),
            H5T_NATIVE_INT);
      TBOX_ASSERT(errf >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(group_id,
            key.c_str(),
            pm_id,
            space,
            H5P_DEFAULT,
            H5P_DEFAULT,
            H5P_DEFAULT);
#else
      hid_t dataset = H5Dcreate(group_id,
            key.c_str(),
            pm_id,
            space,
            H5P_DEFAULT);
#endif
      TBOX_ASSERT(dataset >= 0);
      errf = H5Dwrite(dataset,
            pm_id,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            data);

      TBOX_ASSERT(errf >= 0);

      errf = H5Sclose(space);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tclose(pm_id);
      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(dataset);
      TBOX_ASSERT(errf >= 0);

   } else {
      TBOX_ERROR("VisItDataWriter::HDFputPatchMapStructArray()"
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
   NULL_USE(errf);
}

/*
 *************************************************************************
 *
 * Create an array of max-min double (mm) structs an HDF database with
 * the specified key name.
 *
 *************************************************************************
 */
void
VisItDataWriter::HDFputPatchMinMaxStructArray(
   const std::string& key,
   const patchMinMaxStruct* data,
   const int nelements,
   const hid_t group_id)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);
   TBOX_ASSERT(nelements > 0);

   herr_t errf = 0;
   if (nelements > 0) {
      hid_t space;
      hsize_t dim[1];
      dim[0] = nelements;
      space = H5Screate_simple(1, dim, 0);
      TBOX_ASSERT(space >= 0);

      hid_t s1_tid = H5Tcreate(H5T_COMPOUND, sizeof(patchMinMaxStruct));
      TBOX_ASSERT(s1_tid >= 0);

      errf = H5Tinsert(s1_tid,
            "data_is_defined",
            HOFFSET(patchMinMaxStruct, patch_data_on_disk),
            H5T_NATIVE_CHAR);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tinsert(s1_tid,
            "material_composition_flag",
            HOFFSET(patchMinMaxStruct, material_composition_code),
            H5T_NATIVE_INT);

      TBOX_ASSERT(errf >= 0);

      errf = H5Tinsert(s1_tid,
            "species_composition_flag",
            HOFFSET(patchMinMaxStruct, species_composition_code),
            H5T_NATIVE_INT);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tinsert(s1_tid,
            "min",
            HOFFSET(patchMinMaxStruct, min),
            H5T_NATIVE_DOUBLE);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tinsert(s1_tid,
            "max",
            HOFFSET(patchMinMaxStruct, max),
            H5T_NATIVE_DOUBLE);
      TBOX_ASSERT(errf >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(group_id,
            key.c_str(),
            s1_tid,
            space,
            H5P_DEFAULT,
            H5P_DEFAULT,
            H5P_DEFAULT);
#else
      hid_t dataset = H5Dcreate(group_id,
            key.c_str(),
            s1_tid,
            space,
            H5P_DEFAULT);
#endif
      TBOX_ASSERT(dataset >= 0);

      errf = H5Dwrite(dataset,
            s1_tid,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            data);
      TBOX_ASSERT(errf >= 0);

      errf = H5Sclose(space);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tclose(s1_tid);
      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(dataset);
      TBOX_ASSERT(errf >= 0);

   } else {
      TBOX_ERROR("VisItDataWriter::HDFputPatchMinMaxStructArray()"
         << "\n    data writer with name " << d_object_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
   NULL_USE(errf);
}

/*
 *************************************************************************
 *
 * Create an array of child/parent pointer (cpp) structs an HDF
 * database with the specified key name.
 *
 *************************************************************************
 */

void
VisItDataWriter::HDFputChildParentStructArray(
   const std::string& key,
   const void* data,
   const int nelements,
   hid_t group_id,
   const int sizeOfStruct,
   const std::string& field_name)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);
   TBOX_ASSERT(nelements > 0);

   herr_t errf = 0;
   if (nelements > 0) {
      hid_t space;
      hsize_t dim[1];
      dim[0] = nelements;
      space = H5Screate_simple(1, dim, 0);
      TBOX_ASSERT(space >= 0);

      hid_t s1_tid = H5Tcreate(H5T_COMPOUND, sizeOfStruct);
      TBOX_ASSERT(s1_tid >= 0);

      errf = H5Tinsert(s1_tid,
            "offset",
            0,
            H5T_NATIVE_INT);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tinsert(s1_tid,
            field_name.c_str(),
            sizeof(int),
            H5T_NATIVE_INT);
      TBOX_ASSERT(errf >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(group_id,
            key.c_str(),
            s1_tid,
            space,
            H5P_DEFAULT,
            H5P_DEFAULT,
            H5P_DEFAULT);
#else
      hid_t dataset = H5Dcreate(group_id,
            key.c_str(),
            s1_tid,
            space,
            H5P_DEFAULT);
#endif
      TBOX_ASSERT(dataset >= 0);

      errf = H5Dwrite(dataset, s1_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      TBOX_ASSERT(errf >= 0);

      errf = H5Sclose(space);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tclose(s1_tid);
      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(dataset);
      TBOX_ASSERT(errf >= 0);

   } else {
      TBOX_ERROR("VisItDataWriter::HDFputChildParentStructArray()"
         << "\n    data writer with name " << d_object_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
   NULL_USE(errf);
}

/*
 *************************************************************************
 *
 * Private function to calculate required size for a buffer
 *
 *************************************************************************
 */

int
VisItDataWriter::getBufferSize(
   const hier::Box& patch_box,
   const hier::IntVector& ghost_cell_width,
   const variable_centering centering)
{
   int buf_size = 1;
   int cen = 0;
   if (centering == VISIT_CELL || centering == VISIT_UNKNOWN_CELL) {
      cen = 1;
   } else if (centering == VISIT_NODE || centering == VISIT_UNKNOWN_NODE) {
      cen = 2;
   }
   const int* lower = &patch_box.lower()[0];
   const int* upper = &patch_box.upper()[0];

   for (int i = 0; i < d_dim.getValue(); ++i) {
      buf_size *= upper[i] - lower[i] + cen + (2 * ghost_cell_width(i));
   }
   return buf_size;
}

/*
 *************************************************************************
 *
 * Dump plotitem fields for debugging purposes
 *
 *************************************************************************
 */
void
VisItDataWriter::dumpItem(
   VisItItem& plotitem,
   std::ostream& os) const
{
   os << "d_var_name: " << plotitem.d_var_name << "\n";
   std::string type;
   if (plotitem.d_var_type == VISIT_SCALAR) {
      type = "SCALAR";
   } else if (plotitem.d_var_type == VISIT_VECTOR) {
      type = "VECTOR";
   } else if (plotitem.d_var_type == VISIT_TENSOR) {
      type = "TENSOR";
   }
   os << "d_var_type: " << type << "\n";
   std::string data_type;
   if (plotitem.d_var_data_type == VISIT_DOUBLE) {
      data_type = "DOUBLE";
   } else if (plotitem.d_var_data_type == VISIT_FLOAT) {
      data_type = "FLOAT";
   } else if (plotitem.d_var_data_type == VISIT_INT) {
      data_type = "INT";
   }
   os << "d_var_data_type: " << data_type << "\n";
   std::string cent;
   if (plotitem.d_var_centering == VISIT_CELL) {
      cent = "CELL";
   } else if (plotitem.d_var_centering == VISIT_UNKNOWN_CELL) {
      cent = "UNKNOWN_CELL";
   } else if (plotitem.d_var_centering == VISIT_NODE) {
      cent = "NODE";
   } else if (plotitem.d_var_centering == VISIT_UNKNOWN_NODE) {
      cent = "VISIT_UNKNOWN_NODE";
   }
   os << "d_var_centering: " << cent << "\n";
   os << "d_patch_data_index: " << plotitem.d_patch_data_index << "\n";
   os << "d_depth: " << plotitem.d_depth << "\n";
   os << "d_start_depth_index: " << plotitem.d_start_depth_index << "\n";
   os << "d_scale_factor: " << plotitem.d_scale_factor << "\n";
   os << "d_is_material_state_variable: "
      << plotitem.d_is_material_state_variable << "\n";
   os << "d_derived_writer ptr: " << plotitem.d_derived_writer << "\n";

   int i;
   for (i = 0; i < plotitem.d_depth; ++i) {
      os << "   comp_name[" << i << "]: "
         << plotitem.d_visit_var_name[i] << "\n";
   }

   os << "d_isa_material: " << plotitem.d_isa_material << "\n";
   os << "d_material_name: " << plotitem.d_material_name << "\n";
   os << "d_materials_writer ptr: " << plotitem.d_materials_writer << "\n";
   if (plotitem.d_isa_material) {
      os << "  Species for this material: " << "\n";
      for (i = 0; i < static_cast<int>(plotitem.d_species_names.size()); ++i) {
         os << "   species[" << i << "]: "
            << plotitem.d_species_names[i] << "\n";
      }
   }

   os << "d_isa_species: " << plotitem.d_isa_species << "\n";
   if (plotitem.d_parent_material_pointer) {
      os << "parent_material: "
         << plotitem.d_parent_material_pointer->d_material_name << "\n";
   } else {
      os << "parent_material: Not Applicable\n";
   }
}

VisItDataWriter::childParentStruct::childParentStruct():
   child(-1),
   parent(-1)
{
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
