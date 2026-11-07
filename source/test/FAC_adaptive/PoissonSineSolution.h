/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   PoissonSineSolution class declaration
 *
 ************************************************************************/
#ifndef included_PoissonSineSolution
#define included_PoissonSineSolution

#include "SAMRAI/tbox/Dimension.h"
#include "SinusoidFcn.h"

#include <string>
#include "SAMRAI/tbox/Database.h"

/*
 * SAMRAI classes
 */
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/solv/RobinBcCoefStrategy.h"

using namespace SAMRAI;

/*!
 * @brief Specialized class to provide sine-solution-specific stuff.
 */
class PoissonSineSolution:
   public solv::RobinBcCoefStrategy
{

public:
   PoissonSineSolution(
      const tbox::Dimension& dim);

   PoissonSineSolution(
      /*! Ojbect name */
      const std::string& object_name
      ,
      const tbox::Dimension& dim
      , /*! Input database */
      tbox::Database& database
      , /*! Standard output stream */
      std::ostream * out_stream = 0
      , /*! Log output stream */
      std::ostream * log_stream = 0);

   virtual ~PoissonSineSolution();

   void
   setFromDatabase(
      tbox::Database& database);

   void
   setNeumannLocation(
      int location_index,
      bool flag = true);

   void
   setPoissonSpecifications(
      /*! Object to set */ solv::PoissonSpecifications& sps,
      /*! C id, if used */ int C_patch_data_id,
      /*! D id, if used */ int D_patch_data_id) const;

   /*!
    * @brief Set parameters living on grid.
    *
    * Ignored data are: diffcoef_data and ccoef_data
    * because they are constant.
    */
   void
   setGridData(
      hier::Patch& patch,
      pdat::CellData<double>& exact_data,
      pdat::CellData<double>& source_data);

   virtual void
   setBcCoefs(
      const std::shared_ptr<pdat::ArrayData<double> >& acoef_data,
      const std::shared_ptr<pdat::ArrayData<double> >& bcoef_data,
      const std::shared_ptr<pdat::ArrayData<double> >& gcoef_data,
      const std::shared_ptr<hier::Variable>& variable,
      const hier::Patch& patch,
      const hier::BoundaryBox& bdry_box,
      const double fill_time = 0.0) const;

   hier::IntVector
   numberOfExtensionsFillable() const;

   friend std::ostream&
   operator << (
      std::ostream& os,
      const PoissonSineSolution& r);

private:
   const tbox::Dimension d_dim;

   bool d_neumann_location[2 * SAMRAI::MAX_DIM_VAL];
   double d_linear_coef;
   SinusoidFcn d_exact;

};

#endif  // included_PoissonSineSolution
