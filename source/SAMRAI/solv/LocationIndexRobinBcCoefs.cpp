/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Robin boundary condition support on cartesian grids.
 *
 ************************************************************************/
#include <stdlib.h>

#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include IOMANIP_HEADER_FILE

namespace SAMRAI {
namespace solv {

/*
 ************************************************************************
 * Constructor using database
 ************************************************************************
 */

LocationIndexRobinBcCoefs::LocationIndexRobinBcCoefs(
   const tbox::Dimension& dim,
   const std::string& object_name,
   const std::shared_ptr<tbox::Database>& input_db):
   d_dim(dim),
   d_object_name(object_name)
{
   TBOX_ASSERT(input_db);

   for (int i = 0; i < 2 * d_dim.getValue(); ++i) {
      d_a_map[i] = tbox::MathUtilities<double>::getSignalingNaN();
      d_b_map[i] = tbox::MathUtilities<double>::getSignalingNaN();
      d_g_map[i] = tbox::MathUtilities<double>::getSignalingNaN();
   }
   getFromInput(input_db);
}

/*
 ************************************************************************
 * Destructor
 ************************************************************************
 */

LocationIndexRobinBcCoefs::~LocationIndexRobinBcCoefs()
{
}

/*
 ********************************************************************
 * Set state from input database
 ********************************************************************
 */

void
LocationIndexRobinBcCoefs::getFromInput(
   const std::shared_ptr<tbox::Database>& input_db)
{
   if (!input_db) {
      return;
   }

   for (int i = 0; i < 2 * d_dim.getValue(); ++i) {
      std::string name = "boundary_" + tbox::Utilities::intToString(i);
      if (input_db->isString(name)) {
         d_a_map[i] = 1.0;
         d_g_map[i] = 0.0;
         std::vector<std::string> specs = input_db->getStringVector(name);
         if (specs[0] == "value") {
            d_a_map[i] = 1.0;
            d_b_map[i] = 0.0;
            if (specs.size() != 2) {
               TBOX_ERROR("LocationIndexRobinBcCoefs::getFromInput error...\n"
                  << "exactly 1 value needed with \"value\" boundary specifier"
                  << std::endl);
            } else {
               d_g_map[i] = atof(specs[1].c_str());
            }
         } else if (specs[0] == "slope") {
            d_a_map[i] = 0.0;
            d_b_map[i] = 1.0;
            if (specs.size() != 2) {
               TBOX_ERROR("LocationIndexRobinBcCoefs::getFromInput error...\n"
                  << "exactly 1 value needed with \"slope\" boundary specifier"
                  << std::endl);
            } else {
               d_g_map[i] = atof(specs[1].c_str());
            }
         } else if (specs[0] == "coefficients") {
            if (specs.size() != 3) {
               TBOX_ERROR("LocationIndexRobinBcCoefs::getFromInput error...\n"
                  << "exactly 2 values needed with \"coefficients\" boundary specifier"
                  << std::endl);
            } else {
               d_a_map[i] = atof(specs[1].c_str());
               d_b_map[i] = atof(specs[2].c_str());
            }
         } else {
            TBOX_ERROR(d_object_name << ": Bad boundary specifier\n"
                                     << "'" << specs[0] << "'.  Use either 'value'\n"
                                     << "'slope' or 'coefficients'.\n");
         }
      } else {
         TBOX_ERROR(d_object_name << ": Missing boundary specifier.\n");
      }
   }
}

/*
 ************************************************************************
 * Set the bc coefficients to their mapped values.
 ************************************************************************
 */

void
LocationIndexRobinBcCoefs::setBcCoefs(
   const std::shared_ptr<pdat::ArrayData<double> >& acoef_data,
   const std::shared_ptr<pdat::ArrayData<double> >& bcoef_data,
   const std::shared_ptr<pdat::ArrayData<double> >& gcoef_data,
   const std::shared_ptr<hier::Variable>& variable,
   const hier::Patch& patch,
   const hier::BoundaryBox& bdry_box,
   double fill_time) const
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY2(d_dim, patch, bdry_box);

   NULL_USE(variable);
   NULL_USE(patch);
   NULL_USE(fill_time);

   int location = bdry_box.getLocationIndex();
   TBOX_ASSERT(location >= 0 && location < 2 * d_dim.getValue());
   if (acoef_data) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *acoef_data);

      acoef_data->fill(d_a_map[location]);
   }
   if (bcoef_data) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *bcoef_data);

      bcoef_data->fill(d_b_map[location]);
   }
   if (gcoef_data) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *gcoef_data);

      gcoef_data->fill(d_g_map[location]);
   }
}

hier::IntVector
LocationIndexRobinBcCoefs::numberOfExtensionsFillable() const
{
   /*
    * Return some really big number.  We have no limits.
    */
   return hier::IntVector(d_dim, 1 << (sizeof(int) - 1));
}

/*
 ************************************************************************
 * Assignment operator
 ************************************************************************
 */

LocationIndexRobinBcCoefs&
LocationIndexRobinBcCoefs::operator = (
   const LocationIndexRobinBcCoefs& r)
{
   d_object_name = r.d_object_name;
   for (int i = 0; i < 2 * d_dim.getValue(); ++i) {
      d_a_map[i] = r.d_a_map[i];
      d_b_map[i] = r.d_b_map[i];
      d_g_map[i] = r.d_g_map[i];
   }
   return *this;
}

}
}
