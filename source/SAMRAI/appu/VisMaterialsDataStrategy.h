/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface for writing material related data to a VisIt
 *                dump file.
 *
 ************************************************************************/

#ifndef included_appu_VisMaterialsDataStrategy
#define included_appu_VisMaterialsDataStrategy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Utilities.h"

#include <string>
#include <vector>

namespace SAMRAI {
namespace appu {

/*!
 * @brief Class VisMaterialsDataStrategy is an abstract base
 * class that defines an interface allowing an VisItDataWriter
 * object to generate plot files that contain material and species
 * fractions, as well as state variables for individual materials.  A
 * concrete object of this type must be registered with the data
 * writer in order to use materials or species with the data writer.
 * The registration of the concrete object is done using the method
 * setMaterialsDataWriter() from the VisItDataWriter class.  VisIt
 * requires that material fractions, species fractions, and material
 * state variables be cell-centered.  If they are not cell-centered in
 * the simulation, it is the job of the relevant packing method to
 * convert them to a cell-centered basis before packing them into the
 * buffer.
 *
 * The concrete strategy object is responsible for supplying an
 * implementation of packMaterialFractionsIntoDoubleBuffer().  If
 * species are used, packSpeciesFractionsIntoDoubleBuffer() must also
 * be implemented.  If material state variables are used,
 * packMaterialStateVariableIntoDoubleBuffer() must be implemented.
 *
 * @see VisItDataWriter
 */

class VisMaterialsDataStrategy
{
public:
   /*!
    * @brief Default constructor for VisMaterialsDataStrategy.
    */
   VisMaterialsDataStrategy();

   /*!
    * @brief Destructor for VisMaterialsDataStrategy.
    */
   virtual ~VisMaterialsDataStrategy();

   /*!
    * @brief Enumerated type for the allowable return values
    * for packMaterialFractionsIntoDoubleBuffer() and
    * packSpeciesFractionsIntoDoubleBuffer().
    * - \b ALL_ZERO   (Fractions are 0.0 for every cell in patch.)
    * - \b ALL_ONE    (Fractions are 1.0 for every cell in patch.)
    * - \b MIXTURE    (There is some of this species/material in one or
    *                    more cells of this patch, but the above two cases
    *                   do not apply.)
    */

   enum PACK_RETURN_TYPE { VISIT_ALLZERO = 0,
                           VISIT_ALLONE = 1,
                           VISIT_MIXED = 2 };

   /*!
    * @brief This function, which must be implemented whenever
    * materials are used, packs cell-centered material fractions for
    * the given material, patch, and region into the given 1D double
    * precision buffer which will already have been allocated.  If a
    * non-zero ghost cell vector was specified when
    * registerMaterialNames() was invoked, then ghost data
    * corresponding to this ghost cell vector must be packed into this
    * double buffer.  The data must be packed into the buffer in
    * column major order, i.e. (f(x_0,y_0,z_0), f(x_1,y_0,z_0),
    * f(x_2,y_0,z_0), ...).
    *
    * This method will be called once for each material for each patch.
    *
    * A enumerated PACK_RETURN_TYPE is used for a return value.  To
    * save space in the visit data file, you may choose to set the return
    * value to indicate if a material does not exist at all on the patch,
    * or if the material exists fully on the patch. A return of ALL_ZERO
    * indicates there is 0% of the material in each of the cells of the
    * patch, while a return type of ALL_ONE indicates the material consumes
    * 100% on each of the cells of the patch.  If the patch has a mixture
    * of the material (i.e. between 0% and 100%) then return MIXTURE.
    *
    * @param buffer Double precision array into which cell-centered
    *  material fractions are packed.
    * @param patch hier::Patch on which fractions are defined.
    * @param region hier::Box region over which to pack fractions.
    * @param material_name String identifier for the material.
    * @return The enumeration constant
    *    VisMaterialsDataStrategy::ALL_ZERO,
    *    VisMaterialsDataStrategy::ALL_ONE,
    *    or VisMaterialsDataStrategy::MIXTURE.
    */
   virtual int
   packMaterialFractionsIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& material_name) const;

   /*!
    * @brief Pack sparse volume fraction data
    *
    * This function, which must be implemented whenever materials are
    * used if the sparse packing format is to be used. Packing traverses the
    * patch in a column major order (i.e. (f(x_0,y_0,z_0), f(x_1,y_0,z_0),
    * f(x_2,y_0,z_0), ...) ). If the current cell is clean, then the material
    * number of the material occupying the cell should be packed in to the
    * buffer. If the cell is mixed, then a negative index (i.e., with numbering
    * begining at -1) into the auxilliary mix_zones, mix_mat, vol_fracs, and
    * next_mat vectors which will contain the sparse representation of the
    * volume fractions (VisIt will correct the negative index internally and
    * offset to a zero based representation).  For each component of a mixed
    * cell packMaterialFractionsIntoDoubleBuffer() should pack an entry in the
    * following vectors: mix_zones: the cell number with which the fraction is
    * associated mix_mat: the material number of the fraction vol_fracs: the
    * volume fraction of the current material next_mat: either the (positive
    * but still offset by one) index to the next mixed component within the
    * auxilliary vectors or a 0, indicating the end of the mixed components for
    * this cell.
    *
    * If a non-zero ghost cell vector was specified when
    * registerMaterialNames() was invoked, then ghost data
    * corresponding to this ghost cell vector must be packed into this
    * double buffer.
    *
    * This method will be called once for each patch.
    *
    * A enumerated PACK_RETURN_TYPE is used for a return value.  To
    * save space in the visit data file, you may choose to set the return
    * value to ALL_ONE to indicate that a single material occupies 100% of
    * each cell on the patch. Otherwise, MIXTURE should be returned even if
    * individual cells are not mixed, as this indicates that the patch contains
    * multiple materials.
    *
    * @param mat_list Double precision array into which the material
    *  numbers (or negative indices) are packed.
    * @param mix_zones std::vector<int> into which the cell number
    *  associated with the mixed components are packed.
    * @param mix_mat std::vector<int> into which the material numbers
    *  of the mixed components are packed.
    * @param vol_fracs std::vector<double> into which the volume fractions
    *  (between 0.0 and 1.0) of the mixed components are packed.
    * @param next_mat std::vector<int> into which the (positive) index of
    *  the next mixed component or a terminating 0 are packed.
    * @param patch hier::Patch on which fractions are defined.
    * @param region hier::Box region over which to pack fractions.
    * @return The enumeration constant
    *    VisMaterialsDataStrategy::ALL_ONE,
    *    or VisMaterialsDataStrategy::MIXTURE.
    */
   virtual int
   packMaterialFractionsIntoSparseBuffers(
      int* mat_list,
      std::vector<int>& mix_zones,
      std::vector<int>& mix_mat,
      std::vector<double>& vol_fracs,
      std::vector<int>& next_mat,
      const hier::Patch& patch,
      const hier::Box& region) const;

   /*!
    * @brief This function packs cell-centered species fractions for
    * the given species.
    *
    * This user supplied function packs species fractions of the
    * given material, patch, and region into the supplied 1D double
    * precision buffer.  If a non-zero ghost cell vector was specified when
    * registerSpeciesNames() was invoked, then ghost data
    * corresponding to this ghost cell vector must be packed into this
    * double buffer.  The data must be packed into the buffer in
    * column major order.
    *
    * This method will be called once for each species for each patch.
    *
    * The method must return a PACK_RETURN_TYPE of ALL_ONE, ALL_ZERO,
    * or MIXED.  See the discussion above for the
    * "packMaterialFractionsIntoDoubleBuffer()" method for an explanation
    * of correct return values.
    *
    * @param buffer Double precision array into which  cell-centered
    *   species fractions are packed.
    * @param patch hier::Patch on which fractions are defined.
    * @param region hier::Box region over which to pack fractions.
    * @param material_name String identifier for the material to
    *  which the species belongs.
    * @param species_name String identifier for the species.
    * @return The enumeration constant
    *    VisMaterialsDataStrategy::ALL_ZERO,
    *    VisMaterialsDataStrategy::ALL_ONE,
    *    or VisMaterialsDataStrategy::MIXTURE.
    */
   virtual int
   packSpeciesFractionsIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& material_name,
      const std::string& species_name) const;
};
}
}
#endif
