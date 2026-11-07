/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   PoissonMultigaussianSolution class declaration
 *
 ************************************************************************/
#ifndef included_PoissonMultigaussianSolution
#define included_PoissonMultigaussianSolution

#include <string>

#include "SAMRAI/tbox/Dimension.h"
#include "GaussianFcn.h"

/*
 * Explanation of using PACKAGE macro to change code
 * when compiling under SAMRAI:
 * - SAMRAI does not define PACKAGE while this file's
 *  native package does.
 * - We don't use the vector template under SAMRAI
 *  because implicit instantiation is disabled there.
 */

#ifdef PACKAGE
#include <vector>
#endif
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
 * @brief Specialized class to provide solution-specific stuff
 * for a solution with multiple Gaussians.
 *
 * A single Gaussian exact solution is
 * @f[ u = e^{-\lambda |r-r_0|^2} @f]
 * where @f$ r_0 @f$ is the center of the Gaussian.
 *
 * The diffusion coefficients are 1.
 *
 * Plugging these into the Poisson equation, we get
 * the following source function
 * @f[ 2 \lambda e^{\lambda |r-r_0|^2} ( 3 + 2 \lambda |r-r0|^2 ) @f]
 *
 * The multi-Gaussian is a sum of Gaussian solutions.
 */
class PoissonMultigaussianSolution:
   public solv::RobinBcCoefStrategy
{

public:
   PoissonMultigaussianSolution(
      const tbox::Dimension& dim);

   PoissonMultigaussianSolution(
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

   virtual ~PoissonMultigaussianSolution();

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

   void
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
      const PoissonMultigaussianSolution& r);

private:
   const tbox::Dimension d_dim;

   //! @brief Gaussian component of solution and source.
   std::vector<GaussianFcn> d_gauss;
#define d_gauss_size d_gauss.size()
#define d_gauss_begin d_gauss.begin()
#define d_gauss_end d_gauss.end()
#define d_gauss_append(ITEM) d_gauss.insert(d_gauss.end(), ITEM)
#define d_gauss_const_iterator std::vector<GaussianFcn>::const_iterator

};

#endif  // included_PoissonMultigaussianSolution
