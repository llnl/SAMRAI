/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   PoissonGaussianSolution class declaration
 *
 ************************************************************************/
#ifndef included_PoissonGaussianSolution
#define included_PoissonGaussianSolution

#include <string>

#include "SAMRAI/tbox/Dimension.h"
#include "GaussianFcn.h"

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
 * @brief Specialized class to provide Gaussian
 * solution-specific stuff.
 *
 * The exact solution is
 * @f[ u = e^{-\lambda |r-r_0|^2} @f]
 * where @f$ r_0 @f$ is the center of the Gaussian.
 *
 * The diffusion coefficients are 1.
 *
 * Plugging these into the Poisson equation, we get
 * the following source function
 * @f[ 2 \lambda e^{\lambda |r-r_0|^2} ( 3 + 2 \lambda |r-r0|^2 ) @f]
 */
class PoissonGaussianSolution:
   public solv::RobinBcCoefStrategy
{

public:
   PoissonGaussianSolution(
      const tbox::Dimension& dim);

   PoissonGaussianSolution(
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

   virtual ~PoissonGaussianSolution();

   void
   setFromDatabase(
      tbox::Database& database);

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

   //! Compute exact solution for a given coordinate.
   double
   exactFcn(
      double x,
      double y) const;
   double
   exactFcn(
      double x,
      double y,
      double z) const;
   //! Compute source for a given coordinate.
   double
   sourceFcn(
      double x,
      double y) const;
   double
   sourceFcn(
      double x,
      double y,
      double z) const;

   friend std::ostream&
   operator << (
      std::ostream& os,
      const PoissonGaussianSolution& r);

private:
   const tbox::Dimension d_dim;

   //! @brief Gaussian component of solution and source.
   GaussianFcn d_gauss;

};

#endif  // included_PoissonGaussianSolution
