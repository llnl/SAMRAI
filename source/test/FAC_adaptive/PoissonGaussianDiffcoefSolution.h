/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   PoissonGaussianDiffcoefSolution class declaration
 *
 ************************************************************************/
#ifndef included_PoissonGaussianDiffcoefSolution
#define included_PoissonGaussianDiffcoefSolution

#include <string>

#include "SAMRAI/tbox/Dimension.h"
#include "SinusoidFcn.h"
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
 * @brief Specialized class to provide Gaussian-diffcoef
 * solution-specific stuff.
 *
 * The following variable coefficient problem is contrived
 * to test the code on variable coefficients.
 *
 * The diffusion coefficients are
 * @f[ D_x = D_y = e^{-\lambda |r-r_0|^2} @f]
 * where @f$ r_0 @f$ is the center of the Gaussian.
 * The exact solution
 * @f[ u = sin(k_x x + \phi_x) sin(k_y y + \phi_y) sin(k_z z + \phi_z) @f]
 * Source term (derived by substituting the exact solution into the PDE) is
 * @f[
 *  e^{\lambda |r-r_0|^2}
 *  \left\{
 *    k^2 s_x s_y s_z
 *    - 2\lambda \left[ k_x(x-x_0) c_x s_y s_z
 + k_y(y-y_0) s_x c_y s_z
 + k_z(z-z_0) s_x s_y c_z
 *               \right ]
 *    \right\}
 * @f]
 * where @f$ k^2 = k_x^2 + k_y^2 @f$,
 * @f$ s_x = sin(k_x x + \phi_x) @f$,
 * @f$ s_y = sin(k_y y + \phi_y) @f$,
 * @f$ s_z = sin(k_z z + \phi_z) @f$,
 * @f$ c_x = cos(k_x x + \phi_x) @f$,
 * @f$ c_y = cos(k_y y + \phi_y) @f$ and
 * @f$ c_z = cos(k_z z + \phi_z) @f$.
 */
class PoissonGaussianDiffcoefSolution:
   public solv::RobinBcCoefStrategy
{

public:
   PoissonGaussianDiffcoefSolution(
      const tbox::Dimension& dim);

   PoissonGaussianDiffcoefSolution(
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

   virtual ~PoissonGaussianDiffcoefSolution();

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
    * Ignored data are: ccoef_data
    * because it is constant.
    */
   void
   setGridData(
      hier::Patch& patch,
      pdat::SideData<double>& diffcoef_data,
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
   //! Compute diffusion coefficient for a given coordinate.
   double
   diffcoefFcn(
      double x,
      double y) const;
   double
   diffcoefFcn(
      double x,
      double y,
      double z) const;

   friend std::ostream&
   operator << (
      std::ostream& os,
      const PoissonGaussianDiffcoefSolution& r);

private:
   const tbox::Dimension d_dim;

   //! @brief Gaussian component of solution and source.
   GaussianFcn d_gcomp;
   //! @brief Sine-Sine component of solution and source.
   SinusoidFcn d_sscomp;
   //! @brief Cosine-Sine component of solution and source.
   SinusoidFcn d_cscomp;
   //! @brief Sine-Cosine component of solution and source.
   SinusoidFcn d_sccomp;
   //@{
   double d_lambda;
   double d_k[SAMRAI::MAX_DIM_VAL];
   double d_p[SAMRAI::MAX_DIM_VAL];
   double d_k2;
   //@}

};

#endif  // included_PoissonGaussianDiffcoefSolution
