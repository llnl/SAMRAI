/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   PoissonPolynomialSolution class declaration
 *
 ************************************************************************/
#ifndef included_PoissonPolynomialSolution
#define included_PoissonPolynomialSolution

#include "SAMRAI/tbox/Dimension.h"
#include "QuarticFcn.h"

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
 * @brief Specialized class to provide polynomial-solution-specific stuff.
 */
class PoissonPolynomialSolution:
   public solv::RobinBcCoefStrategy
{

public:
   PoissonPolynomialSolution(
      const tbox::Dimension& dim);

   PoissonPolynomialSolution(
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

   virtual ~PoissonPolynomialSolution();

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

   friend std::ostream&
   operator << (
      std::ostream& os,
      const PoissonPolynomialSolution& r);

private:
   const tbox::Dimension d_dim;

   QuarticFcn d_exact;
   QuarticFcn d_source;

};

#endif  // included_PoissonPolynomialSolution
