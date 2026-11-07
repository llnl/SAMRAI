/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Writes data files for visualization by VisIt
 *
 ************************************************************************/

#ifndef included_appu_VisItDataWriter
#define included_appu_VisItDataWriter

#include "SAMRAI/SAMRAI_config.h"

/*
 ************************************************************************
 *  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT HDF5
 ************************************************************************
 */
#ifdef HAVE_HDF5

#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/appu/VisMaterialsDataStrategy.h"

#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/HDFDatabase.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

#include <string>
#include <list>
#include <vector>
#include <memory>

namespace SAMRAI {
namespace appu {

/*!
 * @brief Class VisItDataWriter is used by SAMRAI-based
 * application codes to generate VisIt data files.  VisIt provides a
 * wide range of visualization and post-processing capabilities.  This
 * class supports both cell-centered and node-centered 2D and 3D AMR
 * data where the underlying data type is either double, float, or
 * int.  Scalar, vector and 2nd-order tensor variables are
 * supported. This class may be used when the mesh geometry is managed
 * by a geom::CartesianGridGeometry object, or when the mesh itself is
 * stored in a state variable to allow moving deformed grids.
 *
 * The VisIt data writer (VDW) supports CELL or NODE centered data
 * of data type double, float, and integer.  There is support to
 * convert an alternative user-defined type into standard cell or node
 * centered data for visualization purposes.  See method descriptions
 * below for more information.  The data actually written by the VDW
 * is of type float so data should not exceed the defined MAX_FLT
 * and/or MIN_FLT.  An optional scale factor may be supplied for
 * data that exceeds these values.
 *
 * After creating a VDW object, the data items to be dumped for
 * plotting must be registered with that object. Thereafter, the
 * object can be used to generate a series of visualization dump files
 * during the execution of an application.  The dumps will include all
 * data items registered with the VDW.
 *
 * This class supports the dumping of data that
 * resides on an AMR patch hierarchy at the time the dump file is
 * written, and derived data, data that does not live on the patch
 * hierarchy, but that can be computed from data which does live on
 * the patch hierarchy.  Derived data requires the user to implement a
 * method for a concrete class derived from the
 * VisDerivedDataStrategy abstract class.
 *
 * The class also supports writing material and species data.  Material
 * is data which indicates for each cell, the fractional amount of
 * a material compound, e.g. "steel", "gas", "copper", etc., contained
 * in that cell. Species are subcomponents of a material, e.g. if the material
 * is "gas", it may have subcomponents such as "nitrogen", "oxygen",
 * etc.  To use material data it is necessary for a user to implement
 * the relevant packing methods for a concrete class derived from the
 * VisMaterialsDataStrategy abstract class.
 *
 * The VDW requires the compilation of SAMRAI with the HDF5 library.
 *
 *
 * A brief summary of the basic steps typically involved in using the
 * VDW are:
 *
 *    -  Create a VisIt data writer object, specifying a name for the
 *       object and the name of the directory to contain the visit dump
 *       files.  An optional argument number_procs_per_file, applicable
 *       to parallel runs, sets the number of processors that share a
 *       common dump file.  This can reduce parallel I/O contention.
 *       The default value of this arg is 1.  If the value specified
 *       is greater than the number of processors, then all processors
 *       share a single dump file.
 *
 *    - Register hierarchy variable data fields using
 *      registerPlotQuantity(). The variables registered may be scalar,
 *      vector, or tensor (depth 1, dim, and dim*dim, respectively)
 *      which is specified in the argument list.  All variables require a
 *      string identifier and an index into the patch data array on the
 *      AMR hierarchy. Optionally, a start depth index may be supplied
 *      to plot a single component of a multi-component field (i.e.
 *      a scalar quantity within a vector, or a vector within a tensor).
 *      A scale factor may also optionally be specified. Lastly, for
 *      user-defined data that is not cell or node type, the proper
 *      centering (cell or node) for the data may be supplied.
 *
 *    - If using derived data, set a default user-defined
 *      derived data writer using setDefaultDerivedDataWriter(). A derived
 *      data writer may also optionally be specified in the call to
 *      registerDerivedPlotQuantity(); if so, this derived data writer
 *      will be used for the variable being registered.  If no derived
 *      writer is specified in the registration call, the default
 *      derived data writer will be used.
 *
 *    - Register derived variable data fields using
 *      registerDerivedPlotQuantity().  All derived variables require a string
 *      identifier. Optionally, a VisDerivedDataStrategy object,
 *      a scale factor, and centering and may be specified.  By default
 *      By default, the data is assumed to be cell-centered.
 *
 *    - The resetLevelPlotQuantity() method is provided for cases when a
 *      variable lives at different patch data indices on different
 *      levels.  Invoking this method redefines the patch data array
 *      index for a given variable on a given level that will be
 *      written to a visit dump file.  Before this method is called,
 *      the variable must be registered using registerPlotQuantity() method.
 *
 *    - If using deformed structured AMR grids (moving grids), register
 *      the coordinates of the nodes using the registerNodeCoordinates()
 *      method.
 *
 *    - If using materials, set the user-defined data writer for material
 *      data using setMaterialsDataWriter().  Then register the names of
 *      the materials being used with registerMaterialNames().
 *
 *    - If using species of the materials, register the names of
 *      the species of each material using the registerSpeciesNames() method.
 *
 *    - The writer will generate VisIt dump files when the method
 *      writePlotData() is called.  Minimally, only a hierarchy and the
 *      time step number is needed.  A simulation time can also be
 *      specified which will be included as part of the file
 *      information in the dump
 *
 *    - The document "Generating VisIt Visualization Data Files in
 *      SAMRAI" in the SAMRAI documentation directory
 *      (docs/userdocs/VisIt-writer.pdf) gives in-depth details on the
 *      use of the VDW and materials, as well as example code fragments
 *      showing how the various VisItDataWriter methods can be embedded
 *      in an application.
 *
 */

class VisItDataWriter
{
public:
   /*!
    * @brief The constructor initializes the VisIt data writer to a
    * default state.
    *
    * The object_name argument is used primarily for error reporting.
    * The dump_directory_name argument is the name of the directory
    * which will contain the visit dump files.  The directory name may
    * include a path.  If the dump directory or any intermediate
    * directories in the path do not exist, they will be created. The
    * optional number_procs_per_file argument is applicable to
    * parallel runs and specifies the number of processors that share
    * a common dump file; the default value is 1. If the specified
    * number_procs_per_file is greater than the number of processors,
    * then all processors share a single vis dump file.  Reducing the
    * number of files written may reduce parallel I/O contention and
    * thus improve I/O efficiency.  The optional argument is_multiblock
    * defaults to false.  It must be set to true for problems on multiblock
    * domains, and left false in all other cases.
    *
    * Before the data writer object can be used for dumping VisIt
    * data, the variables and material-related data (if any) must be
    * registered.
    *
    * An error results and the program will halt if:
    *   - the data is not 2D nor 3D, i.e. dim != 2 and dim != 3
    *   - when assertion checking is active, the object name string is
    *     empty, or the number_procs_per_file is <= 0.
    *
    * @param dim
    * @param object_name String name for data writer object
    * @param dump_directory_name String name for dump directory, which
    *    may include a path.
    * @param number_procs_per_file Optional integer number processors
    *    (>= 1) to share a common dump file; default is 1.
    * @param is_multiblock Optional argument should be set to true only
    *    only for problems on a multiblock domain.
    *
    * @pre !object_name.empty()
    * @pre number_procs_per_file > 0
    * @pre (dim.getValue() == 2) || (dim.getValue() == 3)
    */
   VisItDataWriter(
      const tbox::Dimension& dim,
      const std::string& object_name,
      const std::string& dump_directory_name,
      int number_procs_per_file = 1,
      bool is_multiblock = false);

   /*!
    * @brief The destructor for a VisItDataWriter object.
    *
    * The destructor for the writer does nothing interesting.
    */
   ~VisItDataWriter();

   /*!
    * @brief This method sets the default data writer to use for derived
    * data.
    *
    * If a non-null derived data writer is supplied by the member function
    * registerDerivedPlotQuantity() it will be used in place of the one
    * given here.
    *
    * An error results and the program will halt if:
    *   - assertion checking is active and the default_derived_writer
    *     pointer is null.
    *
    * @param default_derived_writer Pointer to a VisDerivedDataStrategy
    *    object.
    *
    * @pre default_derived_writer != 0
    */
   void
   setDefaultDerivedDataWriter(
      VisDerivedDataStrategy* default_derived_writer)
   {
      TBOX_ASSERT(default_derived_writer != 0);
      d_default_derived_writer = default_derived_writer;
   }

   /*!
    * @brief This method sets the data writer to use for materials.
    *
    * An error results and the program will halt if:
    *   - assertion checking is active and the materials_data__writer
    *     pointer is null.
    *
    * @param materials_data_writer Pointer to a VisMaterialsDataStrategy
    *    object.
    *
    * @pre materials_data_writer != 0
    */
   void
   setMaterialsDataWriter(
      VisMaterialsDataStrategy* materials_data_writer)
   {
      TBOX_ASSERT(materials_data_writer != 0);
      d_materials_writer = materials_data_writer;
   }

   /*!
    * @brief This method registers a variable with the VisIt data writer.
    *
    * Each plot quantity requires a variable name, which is what VisIt
    * will label the plotted quantity.  The variable type is a string
    * specifying either "SCALAR", "VECTOR", or "TENSOR".  By default, the
    * dimension of a scalar variable is 1, vector is dim, and tensor is
    * dim*dim. The integer patch data array index and optional depth
    * index indicate where the data may be found on patches in the hierarchy.
    *
    * A number of optional parameters may be used to further specify
    * characteristics of the plotted variable. The start depth index
    * allows subsets of variables with depth greater than than the supplied
    * type (scalar, vector, tensor) to be specified. For example, a single
    * depth index of a variable with depth greater than one may be registered
    * as a scalar.  A scale factor may be specified such that each data value
    * is multiplied by this factor before being written to the file. Finally,
    * Finally, the variable centering may be specified for data that is not
    * standard CELL or NODE centered types. By default, the writer will set
    * the centering according to the type of data in the supplied patch data
    * index. It will revert to the supplied type only if it is unable to
    * determine the type from the index.
    *
    * Data does not need to exist on all patches or all levels.
    *
    * An error results and the program will halt if:
    *   - a variable was previously registered with the same name.
    *   - the variable type is not "SCALAR", "VECTOR", or "TENSOR"
    *   - the patch data factory referred to by the patch data array
    *     index is null.
    *   - the start depth index is invalid.
    *   - the supplied variable centering is not "CELL" or "NODE".
    *
    * @param variable_name name of variable.
    * @param variable_type "SCALAR", "VECTOR", "TENSOR"
    * @param patch_data_index patch data descriptor id
    * @param start_depth_index (optional) zero by default; may specify
    *    starting index if patch_data_id has depth greater than
    *    the supplied variable_type
    * @param scale_factor (optional) scale factor for data
    * @param variable_centering (optional) "CELL" or "NODE" - used
    *    only when data being registered is not standard cell or
    *    node type.
    *
    * @pre !variable_name.empty()
    * @pre !variable_type.empty()
    * @pre patch_data_index >= -1
    * @pre start_depth_index >= 0
    */
   void
   registerPlotQuantity(
      const std::string& variable_name,
      const std::string& variable_type,
      const int patch_data_index,
      const int start_depth_index = 0,
      const double scale_factor = 1.0,
      const std::string& variable_centering = "UNKNOWN");

   /*!
    * @brief This method registers a derived variable with the VisIt data
    * writer.
    *
    * Each derived variable requires a variable name, which is what VisIt
    * will label the plotted quantity.  The variable type is a string
    * specifying either "SCALAR", "VECTOR", or "TENSOR".  By default, the
    * dimension of a scalar variable is 1, vector is dim, and tensor is
    * dim*dim.  The derived writer should implement methods defined in
    * the derived data strategy which compute the derived data.
    *
    * Optional parameters may be used to further define characteristics
    * of the derived variable. A scale factor may be specified such that
    * each data value is multiplied by this factor before being written
    * to the file.  The variable centering should specify the variable as
    * "CELL" or "NODE" type (if unspecified, "CELL" is used by default).
    *
    * An error results and the program will halt if:
    *   - a variable was previously registered with the same name.
    *   - the variable type is not "SCALAR", "VECTOR", or "TENSOR"
    *   - the derived_writer is NULL and there is no default derived writer.
    *   - the supplied variable centering is not "CELL" or "NODE".
    *
    * @param variable_name name of variable.
    * @param variable_type "SCALAR", "VECTOR", "TENSOR"
    * @param derived_writer (optional) derived data strategy
    *    object to use to calculate this derived data - will use default
    *    derived data object if not supplied
    * @param scale_factor (optional) scale factor with which to multiply
    *    all data values
    * @param variable_centering (optional) centering of derived data - "CELL"
    *    or "NODE"
    * @param variable_mix_type (optional) indicate whether or not the mixed
    *    material state will be stored, "MIXED", or the default of using cell
    *    averages "CLEAN". If "MIXED" then
    *    packMixedDerivedDataIntoDoubleBuffer() must be provided.
    *
    * @pre !variable_name.empty()
    * @pre !variable_type.empty()
    * @pre (variable_name != "Coords") ||
    *      ((variable_type == "VECTOR") && (variable_centering == "NODE"))
    */
   void
   registerDerivedPlotQuantity(
      const std::string& variable_name,
      const std::string& variable_type,
      VisDerivedDataStrategy* derived_writer = 0,
      const double scale_factor = 1.0,
      const std::string& variable_centering = "CELL",
      const std::string& variable_mix_type = "CLEAN");

   /*!
    * @brief This method resets the patch_data_index, and/or
    * the depth_index, at a specific level, of a previously registered
    * plot variable.
    *
    * The change redefines the patch data object written to the plot
    * file on the specified level to the data at the new patch data
    * array index / depth index.  This method is used when a
    * particular variable lives at different patch data slots
    * on different hierarchy levels.  For example, suppose a
    * variable lives at a patch data array index on every level except
    * the finest hierarchy level, where it lives at a different index.
    * First, the variable must be registered using
    * registerPlotQuantity().  Second, the patch data index for
    * the finest hierarchy level is reset using this method. When the
    * data is plotted, it will appear on all levels in the
    * hierarchy. The patch data array index must refer to data with
    * the same type (SCALAR, VECTOR, or TENSOR), centering (NODE or CELL),
    * and data type (int, float, double) as the data for which the
    * variable was originally registered.
    *
    * An error results and the program will halt if:
    *   - this variable name was not previously registered.
    *   - the patch data referred to by the patch data array index
    *     is null, or the data type is not the same type as the data
    *     which was originally registered.
    *   - the depth index is invalid.
    *
    * @param variable_name name of variable.
    * @param level_number level number on which data index is being reset.
    * @param patch_data_index new patch data array index.
    * @param start_depth_index (optional) argument indicating the new depth
    *    index.
    *
    * @pre !variable_name.empty()
    * @pre level_number >= 0
    * @pre patch_data_index >= -1
    * @pre start_depth_index >= 0
    */
   void
   resetLevelPlotQuantity(
      const std::string& variable_name,
      const int level_number,
      const int patch_data_index,
      const int start_depth_index = 0);

   /*!
    * @brief This method is used to register node coordinates for
    * deformed structured AMR grids (moving grids).
    *
    * The patch data index must correspond to an dim-dimensional vector
    * that defines the coordinate location ([X,Y] in 2D, [X,Y,Z] in 3D).
    * The data defining the node locations must be node centered.
    *
    * An error results and the program will halt if:
    *   - the patch data array index is invalid.
    *   - the depth of the patch data index is less than dim.
    *
    * If the nodal coordinates are not in a NodeData object on the
    * hierarchy, you can use registerDerivedPlotQuantity() with the
    * variable name of "Coords", the type of "VECTOR" and the
    * variable_centering of "NODE".
    *
    * @param patch_data_index patch data index of the coordinate data.
    * @param start_depth_index (optional) start index for case where
    *    coordinate data is a subcomponent of a larger patch data vector
    *
    * @pre patch_data_index >= -1
    * @pre start_depth_index >= 0
    */
   void
   registerNodeCoordinates(
      const int patch_data_index,
      const int start_depth_index = 0);

   /*!
    * @brief Same as above method, but allows registration of single
    * coordinate for deformed structured AMR grids (moving grids).
    *
    * The coordinate number should be 0, 1, or 2 for X, Y, and Z directions,
    * respectively.  The patch data index must either be a scalar,  or a
    * vector with an appropriate depth index.  A scale factor may be used to
    * scale grid data.
    *
    * If the nodal coordinates are not in a NodeData object on the
    * hierarchy, you can use registerDerivedPlotQuantity() with the
    * variable name of "Coords", the type of "VECTOR" and the
    * variable_centering of "NODE".
    *
    * @param coordinate_number must be 0 or 1 for 2D, or 0,1,2 for 3D.
    * @param patch_data_index patch data index of the coordinate data.
    * @param depth_index (optional) index for case where
    *    coordinate data is a subcomponent of a larger patch data vector
    * @param scale_factor scale factor with which to multiply
    *    coordinate data values
    *
    * @pre (coordinate_number >= 0) && (coordinate_number < d_dim.getValue())
    * @pre patch_data_index >= -1
    * @pre depth_index >= 0
    */
   void
   registerSingleNodeCoordinate(
      const int coordinate_number,
      const int patch_data_index,
      const int depth_index = 0,
      const double scale_factor = 1.0);

   /*!
    * @brief This method registers with the VisIt data writer the
    * names of materials being used in the simulation.
    *
    * When a mesh has materials defined over it, every cell will
    * contain a fractional amount f (0 <= f <= 1.0) of every material,
    * called a material fraction or volume fraction.  The sum of all
    * material fractions for every cell must be 1.0.  Since materials
    * are defined over a cell, each materials variable is assumed to
    * be CELL centered.
    *
    * In order to use materials with VisIt, the application class
    * needs to inherit from VisMaterialsDataStrategy, an
    * abstract base class which defines an interface for writing out
    * various materials related fields.  A concrete object of this
    * strategy must be registered with the VisIt data writer.  That
    * concrete object is responsible for providing a concrete
    * implementation of the method packMaterialFractionsIntoDoubleBuffer()
    * which writes out the material fractions for each registered material
    * over a given patch.
    *
    * VisIt uses the material fractions to calculate interpolated
    * material boundaries in cells with multiple materials.  Therefore
    * VisIt can display a material (or materials) as multiple colored
    * contiguous regions.  In addition, VisIt can intersect this
    * volume field with a plane and get accurate 2D boundaries within
    * cells.  If desired, the material fractions for a material can be
    * treated as a scalar field and the usual scalar field plot tools,
    * such as pseudocolor and contour, can be applied.
    *
    * Because species are a subset of materials, it is required that the
    * material names be registered before any species names are registered.
    * New materials may not be added during the simulation; that is, this
    * method should only be called once.
    *
    * An error results and the program will halt if:
    *   - this method is called more than once.
    *   - the new registerSparseMaterialNames() method is called
    *   - when assertion checking is active, the number of names = 0,
    *     or any material name string is empty.
    *
    * @param material_names tbox::Array of strings: the names of the materials.
    *
    * @pre material_names.size() > 0
    * @pre for each member of material_names, mn, !mn.empty()
    * @pre d_materials_names.size() == 0
    */
   void
   registerMaterialNames(
      const std::vector<std::string>& material_names);

   /*!
    * @brief This method registers with the VisIt data writer the
    * names of materials being used in the simulation.
    *
    * When a mesh has materials defined over it, every cell will
    * contain a fractional amount f (0 <= f <= 1.0) of every material,
    * called a material fraction or volume fraction.  The sum of all
    * material fractions for every cell must be 1.0.  Since materials
    * are defined over a cell, each materials variable is assumed to
    * be CELL centered.
    *
    * In order to use materials with VisIt, the application class
    * needs to inherit from VisMaterialsDataStrategy, an
    * abstract base class which defines an interface for writing out
    * various materials related fields.  A concrete object of this
    * strategy must be registered with the VisIt data writer.  That
    * concrete object is responsible for providing a concrete
    * implementation of the method packMaterialFractionsIntoSparseBuffers()
    * which writes out the material fractions for each registered material
    * over a given patch.
    *
    * VisIt uses the material fractions to calculate interpolated
    * material boundaries in cells with multiple materials.  Therefore
    * VisIt can display a material (or materials) as multiple colored
    * contiguous regions.  In addition, VisIt can intersect this
    * volume field with a plane and get accurate 2D boundaries within
    * cells.  If desired, the material fractions for a material can be
    * treated as a scalar field and the usual scalar field plot tools,
    * such as pseudocolor and contour, can be applied.
    *
    * An error results and the program will halt if:
    *   - this method is called more than once.
    *   - the legacy registerMaterialNames() method is called
    *   - when assertion checking is active, the number of names = 0,
    *     or any material name string is empty.
    *
    * @param material_names tbox::Array of strings: the names of the materials.
    *
    * @pre material_names.size() > 0
    * @pre for each member of material_names, mn, !mn.empty()
    * @pre d_materials_names.size() == 0
    */
   void
   registerSparseMaterialNames(
      const std::vector<std::string>& material_names);

   /*!
    * @brief This method registers the names of the species for a
    * material_name with the VisIt data writer.
    *
    * Species are subcomponents of a material. For example, a simulation
    * with 4 materials (e.g. "copper", "gold", "liquid" and "gas") may have
    * one of the materials (e.g. "gas") that is composed of multiple
    * species (e.g. "nitrogen" and "helium"). Each species is associated
    * with a particular material, so it is necessary that the
    * "registerMaterialNames()" method is called before this method.
    *
    * For each species name registered, there is an associated species
    * fraction,  which must be between 0 and 1.0 on each cell. The sum of
    * the species fractions for all species of a given material must
    * equal 1.0 in every cell in which that material appears.  In VisIt,
    * the species fractions for a species can be treated as
    * a scalar field and the usual scalar plot operations applied to
    * the field to show "concentrations".  So, in the example described
    * above, the percentage of nitrogen in the gas can be displayed in a
    * pseudocolor plot.
    *
    * In order to use species with VisIt, the application class
    * needs to inherit from VisMaterialsDataStrategy, an
    * abstract base class which defines an interface for writing out
    * species related fields.  A concrete object of the materials data
    * strategy class is responsible for providing an implementation of the
    * method packSpeciesFractionsIntoDoubleBuffer() to write out
    * the species fractions field for each registered species over a
    * given patch.
    *
    * An error results and the program will halt if:
    *   - this method is called before registerMaterialNames() is
    *     called.
    *   - the supplied material name is not in the list of materials
    *     supplied in the registerMaterialNames() method.
    *   - the method is called more than once for a given material
    *   - when assertion checking is active, the number of names = 0,
    *     or any name string is empty.
    *
    * @param material_name String name of the material whose species
    *    names are being registered.
    * @param species_names tbox::Array of strings: the names of the species
    *    for material_name.
    *
    * @pre !material_name.empty()
    * @pre species_names.size() > 0
    * @pre d_materials_names.size() > 0
    */
   void
   registerSpeciesNames(
      const std::string& material_name,
      const std::vector<std::string>& species_names);

   /*!
    * @brief This method registers expressions that will be embedded in the
    * VisIt datafiles.
    *
    * The three Arrays of strings define the names (or keys) of the
    * expressions, the expressions, and the types of the expressions
    * (scalar, vector, tensor). For more information on defining VisIt
    * expressions see the VisItUsersManual.
    */
   void
   registerVisItExpressions(
      const std::vector<std::string>& expression_keys,
      const std::vector<std::string>& expressions,
      const std::vector<std::string>& expression_types);

   /*!
    * @brief This method causes the VisIt data writer to dump all
    * registered data. The appropriate packing methods will be invoked
    * for material-related data and derived data.
    *
    * The time step number is used as a file name extension for the
    * dump files. It must be non-negative and greater than the
    * previous time step, if any. A simulation time may be provided as
    * an optional argument.  If this time is not specified, a default
    * value of zero is used.  The simulation time can be accessed in
    * VisIt's "File Information" dialog box.
    *
    * An error results and the program will halt if:
    *   - materials have been registered, but setMaterialsDataWriter() has
    *     not been called.
    *   - the time step number is <= the previous time step number.
    *   - when assertion checking is active, the hierarchy pointer is null,
    *     the time step is < 0, or the dump directory name string is empty,
    *
    * @param hierarchy A pointer to the patch hierarchy on which the data
    *    to be plotted is defined.
    * @param time_step Non-negative integer value specifying the current
    *    time step number.
    * @param simulation_time Optional argument specifying the double
    *    precision simulation time. Default is 0.0.
    *
    * @pre hierarchy
    * @pre time_step_number >= 0
    * @pre !d_top_level_directory_name.empty()
    */
   void
   writePlotData(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      int time_step,
      double simulation_time = 0.0);

   /*!
    * @brief Set the name of the summary file.
    *
    * This sets the summary file written at each step of the simulation
    * which describes the data contained in the visit files written by
    * each MPI process.  The supplied string is appended with ".samrai"
    * so the actual name of the file will be "<filename>.samrai".  If no
    * alternative name is supplied, by default the summary file used is
    * "summary.samrai".
    *
    * @pre !filename.empty()
    */
   void
   setSummaryFilename(
      std::string& filename)
   {
      TBOX_ASSERT(!filename.empty());
      d_summary_filename = filename + ".samrai";
   }

   /*!
    * @brief Returns the object name.
    *
    * @return The object name.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

   /*!
    * @brief Set flag for writing ghosts
    *
    * @param write ghosts   True to write ghosts, false to not write ghosts
    */
   void
   setWriteGhosts(bool write_ghosts)
   {
      d_write_ghosts = write_ghosts; 
   }

private:
   /*
    * Static integer constant describing version of VisIt Data Writer.
    */
   static const float VISIT_DATAWRITER_VERSION_NUMBER;

   /*
    * Static integer constant describing the maximum number of components
    * ever written.
    */
   static const int VISIT_MAX_NUMBER_COMPONENTS = 100;

   /*
    * Static integer constant describing the largest base space dimension
    * ever written.
    */
   static const int VISIT_FIXED_DIM = 3;

   /*
    * Static integer constant describing the maximum size of a C char string.
    */
   static const int VISIT_NAME_BUFSIZE;

   /*
    * Static integer constant describing undefined index.
    */
   static const int VISIT_UNDEFINED_INDEX;

   /*
    * Static integer constant describing process which writes single summary
    * file with information from all processors for parallel runs
    */
   static const int VISIT_MASTER;

   /*
    * Static integer constant describing MPI message tag.
    */
   static const int VISIT_FILE_CLUSTER_WRITE_BATON;

   /*
    * Static boolean that specifies if the summary file (d_summary_filename)
    * has been opened.
    */
   static bool s_summary_file_opened;

   /*
    * Struct used to gather min/max information, and
    * to track floating point overflows of data.
    */
   struct patchMinMaxStruct {
      double min;
      double max;
      int material_composition_code;
      int species_composition_code;
      char patch_data_on_disk;
   };

   /*
    * Struct to hold patch extents.
    */
   struct patchExtentsStruct {
      int lower[VISIT_FIXED_DIM];
      int upper[VISIT_FIXED_DIM];
      double xlo[VISIT_FIXED_DIM];
      double xhi[VISIT_FIXED_DIM];
   };

   /*
    * Struct to hold patch processor mapping info.
    */
   struct patchMapStruct {
      int file_cluster_number;
      int processor_number;
      int level_number;
      int patch_number;
   };

   /*
    * Struct used to hold parent and child information for
    * writing data from multiple tbox::MPI processes.
    */
   struct childParentStruct {
      childParentStruct();
      int child;
      int parent;
   };

   /*
    * hier::Variable type:
    *   SCALAR - scalar plot variable (depth = 1)
    *   VECTOR - vector plot variable (depth = dim)
    *   TENSOR - tensor plot variable (depth = dim*dim)
    */
   enum variable_type { VISIT_SCALAR = 0,
                        VISIT_VECTOR = 1,
                        VISIT_TENSOR = 2 };

   /*
    * hier::Variable data type  - float, double, integer
    */
   enum variable_data_type { VISIT_INT = 3,
                             VISIT_FLOAT = 4,
                             VISIT_DOUBLE = 5,
                             VISIT_DATA_TYPE_BAD = 990 };

   /*
    * hier::Variable centering:
    *   CELL         - standard cell centered
    *   NODE         - standard node centered
    *   UNKNOWN_CELL - unknown type, cast to cell centered
    *   UNKNOWN_NODE - unknown type, cast to node centered
    */
   enum variable_centering { VISIT_CELL = 6,
                             VISIT_NODE = 7,
                             VISIT_UNKNOWN_CELL = 8,
                             VISIT_UNKNOWN_NODE = 9,
                             VISIT_CENTERING_BAD = 991 };

   /*
    * Grid type:
    *   CARTESIAN - standard cartesian grid
    *   DEFORMED  - node centered grid where nodes may be deformed
    *               (e.g. sometimes called curvilinear)
    */
   enum grid_type { VISIT_CARTESIAN = 10,
                    VISIT_DEFORMED = 11 };

   /*
    * The following structure is used to store data about each item
    * to be written to a plot file.
    *
    * Standard information (user supplied):
    *   d_var_name - string variable name:
    *   d_var_type - SCALAR, VECTOR, TENSOR
    *   d_var_data_type - INT, FLOAT, DOUBLE
    *   d_var_centering - CELL, NODE, UNKNOWN_CELL, UNKNOWN_NODE
    *   d_patch_data_index - int patch data id
    *   d_depth - int depth of patch data
    *   d_start_depth_index - int starting depth for vector data
    *   d_scale_factor - dbl scaling factor
    *   d_derived_writer - ptr to derived data writer (NULL if not DERIVED)
    *   d_is_material_state_variable - bool, true if mixed state var, in which
    *       case d_derived_writer must be provided and provide
    *       packMixedDerivedDataIntoDoubleBuffer
    *   d_is_species_state_variable - bool, true if mixed state var, in which
    *       case d_derived_writer must be provided and provide
    *       packMixedDerivedDataIntoDoubleBuffer
    *
    * Standard information (writer internal):
    *   d_data_type - INT, FLOAT, DOUBLE
    *   d_visit_var_name - internally maintained visit var name:
    *       for scalar: name[0]  = d_variable_name
    *       for vector: name[0]  = d_variable_name.00,
    *                   name[1]  = d_variable_name.01,
    *                   ..
    *                   name[nn] = d_variable_name.nn
    *   d_master_min_max - ptr to min/max struct on master for each
    *      var.  This is used to gather the min max info for all patches.
    *   d_deformed_coord_id - id of vector defining deformed coordinates
    *   d_coord_scale_factor - scale factor of the different deformed coords
    *   d_level_start_depth_index - int array specifying start depth on
    *                               each level
    *
    * Material information
    *   d_isa_material - boolean specifying if variable is a material
    *   d_material_name - string name of material
    *   d_species_names - string array names of the species
    *   d_materials_writer - ptr to user-supplied materials writer
    *   d_material_name_HDFGroup - hdf group for material
    *
    * Species information
    *   d_isa_species - boolean specifying if variable is a species
    *   d_species_name - string name of species
    *   d_parent_material_pointer - VisItItem ptr to parent material
    *   d_species_HDFGroup - hdf group for species
    *   d_extents_species_HDFGroup - hdf group for species extents
    */
   struct VisItItem {

      /*
       * Standard information (user supplied)
       */
      std::string d_var_name;
      variable_type d_var_type;
      variable_centering d_var_centering;
      int d_patch_data_index;
      int d_start_depth_index;
      double d_scale_factor;
      bool d_is_derived;
      VisDerivedDataStrategy* d_derived_writer;
      bool d_is_deformed_coords;
      // new flag for mixed/clean state variables
      bool d_is_material_state_variable;
      // Do we want to extend this for species, or can that fit into the
      //   material state variable treatment?
      //bool d_is_species_state_variable;
      std::vector<int> d_ghost_width;

      /*
       * Standard information (writer generated)
       */
      int d_depth;
      variable_data_type d_var_data_type;
      std::vector<std::string> d_visit_var_name;
      struct patchMinMaxStruct *
      d_master_min_max[VISIT_MAX_NUMBER_COMPONENTS];
      std::vector<int> d_level_patch_data_index;
      std::vector<double> d_coord_scale_factor;

      /*
       * Material information
       */
      bool d_isa_material;
      std::string d_material_name;
      std::vector<std::string> d_species_names;
      VisMaterialsDataStrategy* d_materials_writer;
      /*
       * Species information
       */
      bool d_isa_species;
      std::string d_species_name;
      VisItItem* d_parent_material_pointer;
      std::shared_ptr<tbox::Database> d_species_HDFGroup;
      std::shared_ptr<tbox::Database> d_extents_species_HDFGroup;
   };

   /*
    * Utility routine to initialize a standard variable for
    * plotting based on user input.  Derived, coordinate, material,
    * and species data all use this method to initialize the variable
    * and then set specific characteristics in their appropriate
    * register methods.
    */
   void
   initializePlotItem(
      VisItItem& plotitem,
      const std::string& variable_name,
      const std::string& variable_type,
      const int patch_data_index,
      const int start_depth_index,
      const double scale_factor,
      const std::string& variable_centering,
      const hier::IntVector& ghost_width);

   /*
    * Utility routine to reset a variable by level for plotting
    * (either vector or scalar variable).
    */
   void
   resetLevelPlotItem(
      const std::string& variable_name,
      int level_number,
      int patch_data_array_index,
      int start_depth_index,
      std::string method_name);

   /*
    * Coordinate writing HDF plot files, both data and summary.
    */
   void
   writeHDFFiles(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      double simulation_time);

   /*
    * Allocate and initialize the min/max structs that hold
    * summary information about each plotted variable.
    */
   void
   initializePlotVariableMinMaxInfo(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy);

   /*
    * Write variable data to HDF file.
    */
   void
   writeVisItVariablesToHDFFile(
      const std::shared_ptr<tbox::Database>& processor_HDFGroup,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      int coarsest_level,
      int finest_level,
      double simulation_time);

   /*
    * Pack regular (i.e. NOT materials or species) and derived data into
    * the supplied HDF database for output.
    */
   void
   packRegularAndDerivedData(
      const std::shared_ptr<tbox::Database>& patch_HDFGroup,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      hier::Patch& patch,
      double simulation_time);

   /*
    * Pack the materials data into the supplied database for output.
    */
   void
   packMaterialsData(
      const std::shared_ptr<tbox::Database>& patch_HDFGroup,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      hier::Patch& patch);

   /*
    * Pack the species data.  The correct HDF database is determined
    * from the parent material.
    */
   void
   packSpeciesData(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      hier::Patch& patch);

   /*
    * Check min/max to make exit cleanly if users data exceeds float
    * min/max values. Otherwise, the writer will dump core when it
    * tries to convert double data to float for writing the vis file.
    */
   void
   checkFloatMinMax(
      const std::string& var_name,
      const double dmin,
      const double dmax,
      const int level_number,
      const int patch_number,
      const int patch_data_id);

   /*
    * Convert level number, patch number, to global patch number.
    */
   int
   getGlobalPatchNumber(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const int patch_number);

   /*
    * Calculate and then write patch parent and
    * child info to summary HDF file.
    */
   void
   writeParentChildInfoToSummaryHDFFile(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const std::shared_ptr<tbox::Database>& basic_HDFGroup);

   /*
    *    Sort function for use by qsort to sort child_parent array
    *    by child patch number so can find all parents of a given child.
    */
   static int
   childParentCompareFunc(
      const void* s1,
      const void* s2);

   /*
    * Barrier functions to enable orderly writing of cluster files by
    * passing a baton from current writer to next proc in cluster.
    */
   void
   dumpWriteBarrierBegin();
   void
   dumpWriteBarrierEnd();

   /*
    * Write summary data for VisIt to HDF file.
    */
   void
   writeSummaryToHDFFile(
      std::string dump_dir_name,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      int coarsest_plot_level,
      int finest_plot_level,
      double simulation_time);

   /*
    * Helper method for writeSummaryToHDFFile() method above.
    * Performs tbox::MPI communications to send min/max information for all
    * variables on all patches to the VISIT_MASTER.
    */
   void
   exchangeMinMaxPatchInformation(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_plot_level,
      const int finest_plot_level);

   /*
    * Pack dim patch data into 1D double precision buffer,
    * eliminating ghost data if necessary
    */
   void
   packPatchDataIntoDoubleBuffer(
      const std::shared_ptr<hier::PatchData>& pdata,
      const int depth_index,
      const variable_data_type type_of_data,
      const hier::Box patch_box,
      double* buffer,
      const variable_centering centering,
      const hier::IntVector& ghost_width);

   /*
    * Create a 2D integer array entry in the database with the specified
    * key name.
    */
   void
   HDFputIntegerArray2D(
      const std::string& key,
      const int* data,
      const int nelements0,
      const int nelements1,
      const hid_t group_id);

   /*
    * Create a 2D double array entry in the database with the specified
    * key name.
    */
   void
   HDFputDoubleArray2D(
      const std::string& key,
      const double* data,
      const int nelements0,
      const int nelements1,
      const hid_t group_id);

   /*
    * Create an array of patch extent structs in the database
    * with the specified key name.
    */
   void
   HDFputPatchExtentsStructArray(
      const std::string& key,
      const patchExtentsStruct* data,
      const int nelements,
      const hid_t group_id);

   /*
    * Create an array that identifies the type of the boundary.
    */
   void
   HDFputBoundaryTypeArray(
      const std::string& key,
      const std::vector<int>& data,
      const int nelements,
      const hid_t group_id);

   /*
    * Create an array of patch map structs in the database
    * with the specified key name.
    */
   void
   HDFputPatchMapStructArray(
      const std::string& key,
      const patchMapStruct* data,
      const int nelements,
      const hid_t group_id);

   /*
    * Create an array of min/max structs in the database with
    * the specified key name.
    */
   void
   HDFputPatchMinMaxStructArray(
      const std::string& key,
      const patchMinMaxStruct* data,
      const int nelements,
      const hid_t group_id);

   /*
    * Create an array of child/parent pointer (CPP) structs in
    * the database with the specified key name.
    */
   void
   HDFputChildParentStructArray(
      const std::string& key,
      const void* data,
      const int nelements,
      const hid_t group_id,
      const int sizeOfStruct,
      const std::string& fieldName);

   /*
    * Calculate buffer size to hold data on patch with given ghost cell
    * width and centering.
    */
   int
   getBufferSize(
      const hier::Box& patch_box,
      const hier::IntVector& ghost_cell_width,
      const variable_centering centering);

   /*
    * Dump item fields for debugging purposes.
    */
   void
   dumpItem(
      VisItItem& plotitem,
      std::ostream& os) const;

   /*
    * Dimension of object
    */
   const tbox::Dimension d_dim;

   /*!
    * @brief Exclusive SAMRAI_MPI duplicated for this object.
    */
   tbox::SAMRAI_MPI d_mpi;

   /*
    * Name of this VisIt data writer object
    */
   std::string d_object_name;

   /*
    * Hierarchy level information to write to a file.
    */
   int d_number_levels;

   /*
    * tbox::Array of mesh-scaling ratios from each level to reference level
    * (i.e., coarsest level).
    */
   std::vector<hier::IntVector> d_scaling_ratios;

   /*
    * Default data writer for user defined data.
    */
   VisDerivedDataStrategy* d_default_derived_writer;

   /*
    * Data writer for materials data.
    */
   VisMaterialsDataStrategy* d_materials_writer;

   /*
    * Directory into which VisIt files will be written.
    */
   std::string d_top_level_directory_name;
   std::string d_current_dump_directory_name;
   std::string d_summary_filename;

   /*
    * Grid type - CARTESIAN or DEFORMED.
    */
   grid_type d_grid_type;

   /*
    * Time step number (passed in by user).
    */
   int d_time_step_number;

   /*
    * Number of worker processors that pass info to VISIT_MASTER
    */
   int d_number_working_slaves;

   /*
    * Cluster information for parallel runs. A file_cluster
    * is a set of processors that all write VisIt data to
    * a single disk file.  d_processor_on_file_cluster[processorNumber]
    * returns the file_clusterNumber of processorNumber.
    * d_file_cluster_leader is controller of file_cluster.
    */
   int d_number_file_clusters;
   int d_my_file_cluster_number;
   std::vector<int> d_processor_in_file_cluster_number;
   bool d_file_cluster_leader;
   int d_file_cluster_size;
   int d_my_rank_in_file_cluster;
   int d_number_files_this_file_cluster;

   /*
    * Number of registered VisIt variables, materials, and species.
    * Each regular and derived and variable (i.e. variables registered with
    * a patch data index) are considered visit variables.  If the grid is
    * deformed, the node coordinates registered are also considered
    * visit variables.  Material and species data are NOT counted in
    * d_number_visit_variables because these are supplied by the user. The
    * d_number_visit_variables does not consider whether the plotted variable
    * is scalar, vector, or tensor.  The d_number_visit_variables_plus_depth
    * does consider this, and is the sume of the depths of all the visit
    * variables, again not counting material and species data.
    *
    * The material names are stored in d_materials_names.  Any or all of the
    * materials may have an associated species, and the TOTAL of these is
    * counted in d_number_species.
    */
   int d_number_visit_variables;
   int d_number_visit_variables_plus_depth;
   std::vector<std::string> d_materials_names;
   int d_number_species;

   /*
    * For parallel runs, this array of min/max structs holds the summary
    * information all local patches on workers prior to being sent to
    * master.  It is filled in this order: by level, then by
    * local_patch_number, then by var_item, then by component_number.
    */
   patchMinMaxStruct* d_worker_min_max;
   int d_var_id_ctr;

   /*
    *  tbox::List of scalar and vector variables registered with
    *  VisItDataWriter.
    */
   std::list<VisItItem> d_plot_items;

   /*
    * Boolean that is set to true only in multiblock problems.
    */
   bool d_is_multiblock;

   /*
    * Boolean to control writing of ghosts
    */
   bool d_write_ghosts;

   /*
    * brief Storage for strings defining VisIt expressions to be embedded in
    * the plot dump.
    */
   std::vector<std::string> d_visit_expression_keys;
   std::vector<std::string> d_visit_expressions;
   std::vector<std::string> d_visit_expression_types;

   //! @brief Timer for writePlotData().
   static std::shared_ptr<tbox::Timer> t_write_plot_data;

   /*!
    * @brief Initialize static objects and register shutdown routine.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback()
   {
      t_write_plot_data = tbox::TimerManager::getManager()->getTimer(
            "appu:VisItDataWriter::writePlotData()");
   }

   /*!
    * @brief Method registered with ShutdownRegister to cleanup statics.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback()
   {
      t_write_plot_data.reset();
   }

   /*
    * Static initialization and cleanup handler.
    */
   static tbox::StartupShutdownManager::Handler
      s_initialize_handler;

};
}
}

#endif
#endif
