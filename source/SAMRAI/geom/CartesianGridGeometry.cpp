/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Simple Cartesian grid geometry for an AMR hierarchy.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianGridGeometry.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

// Cell data coarsen operators
#include "SAMRAI/geom/CartesianCellComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianCellDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianCellFloatWeightedAverage.h"

// Cell data refine operators
#include "SAMRAI/geom/CartesianCellComplexConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellComplexLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleLinearRefine.h"
#include "SAMRAI/geom/CartesianCellFloatConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellFloatLinearRefine.h"
#include "SAMRAI/geom/CartesianCellConservativeLinearRefine.h"

// Edge data coarsen operators
#include "SAMRAI/geom/CartesianEdgeComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianEdgeDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianEdgeFloatWeightedAverage.h"

// Edge data refine operators
#include "SAMRAI/geom/CartesianEdgeDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianEdgeFloatConservativeLinearRefine.h"

// Face data coarsen operators
#include "SAMRAI/geom/CartesianFaceComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianFaceDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianFaceFloatWeightedAverage.h"

// Face data refine operators
#include "SAMRAI/geom/CartesianFaceDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianFaceFloatConservativeLinearRefine.h"

// Node data refine operators
#include "SAMRAI/geom/CartesianNodeComplexLinearRefine.h"
#include "SAMRAI/geom/CartesianNodeDoubleLinearRefine.h"
#include "SAMRAI/geom/CartesianNodeFloatLinearRefine.h"

// Outerface data coarsen operators
#include "SAMRAI/geom/CartesianOuterfaceComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianOuterfaceDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianOuterfaceFloatWeightedAverage.h"

// Outerside data coarsen operators
#include "SAMRAI/geom/CartesianOutersideDoubleWeightedAverage.h"

// Side data coarsen operators
#include "SAMRAI/geom/CartesianSideComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianSideDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianSideFloatWeightedAverage.h"

// Side data refine operators
#include "SAMRAI/geom/CartesianSideDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianSideFloatConservativeLinearRefine.h"

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/pdat/SideVariable.h"

#include "SAMRAI/hier/BoundaryLookupTable.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <cstdlib>
#include <fstream>
#include <typeinfo>

namespace SAMRAI {
namespace geom {

const int CartesianGridGeometry::GEOM_CARTESIAN_GRID_GEOMETRY_VERSION = 2;

/*
 *************************************************************************
 *
 * Constructors for CartesianGridGeometry.  Both set up operator
 * handlers and register the geometry object with the RestartManager.
 * However, one initializes data members based on arguments.
 * The other initializes the object based on input file information.
 *
 *************************************************************************
 */
CartesianGridGeometry::CartesianGridGeometry(
   const tbox::Dimension& dim,
   const std::string& object_name,
   const std::shared_ptr<tbox::Database>& input_db):
   GridGeometry(dim, object_name, input_db, false),
   d_domain_box(dim)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(input_db);

   buildOperators();

   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }

   getFromInput(input_db, is_from_restart);
}

CartesianGridGeometry::CartesianGridGeometry(
   const std::string& object_name,
   const double* x_lo,
   const double* x_up,
   hier::BoxContainer& domain):
   GridGeometry(object_name, domain),
   d_domain_box(domain.front().getDim())
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(domain.size() > 0);
   TBOX_ASSERT(x_lo != 0);
   TBOX_ASSERT(x_up != 0);

   buildOperators();

   setGeometryData(x_lo, x_up, domain);
}

CartesianGridGeometry::CartesianGridGeometry(
   const std::string& object_name,
   const double* x_lo,
   const double* x_up,
   hier::BoxContainer& domain,
   const std::shared_ptr<hier::TransferOperatorRegistry>& op_reg):
   GridGeometry(object_name, domain, op_reg),
   d_domain_box(domain.front().getDim())
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(domain.size() > 0);
   TBOX_ASSERT(x_lo != 0);
   TBOX_ASSERT(x_up != 0);

   buildOperators();

   setGeometryData(x_lo, x_up, domain);
}

/*
 *************************************************************************
 *
 * Destructor for CartesianGridGeometry deallocates grid storage.
 *
 *************************************************************************
 */

CartesianGridGeometry::~CartesianGridGeometry()
{
}

/*
 *************************************************************************
 *
 * Create and return pointer to refined version of this Cartesian
 * grid geometry object refined by the given ratio.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BaseGridGeometry>
CartesianGridGeometry::makeRefinedGridGeometry(
   const std::string& fine_geom_name,
   const hier::IntVector& refine_ratio) const
{
   const tbox::Dimension dim(getDim());

   TBOX_ASSERT(!fine_geom_name.empty());
   TBOX_ASSERT(fine_geom_name != getObjectName());
   TBOX_ASSERT(refine_ratio > hier::IntVector::getZero(dim));

   hier::BoxContainer fine_domain(getPhysicalDomain());
   fine_domain.refine(refine_ratio);

   std::shared_ptr<hier::BaseGridGeometry> fine_geometry(
      new CartesianGridGeometry(fine_geom_name,
         d_x_lo,
         d_x_up,
         fine_domain,
         d_transfer_operator_registry));

   fine_geometry->initializePeriodicShift(getPeriodicShift(refine_ratio));

   return fine_geometry;
}

/*
 *************************************************************************
 *
 * Create and return pointer to coarsened version of this Cartesian
 * grid geometry object coarsened by the given ratio.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BaseGridGeometry>
CartesianGridGeometry::makeCoarsenedGridGeometry(
   const std::string& coarse_geom_name,
   const hier::IntVector& coarsen_ratio) const
{
   TBOX_ASSERT(!coarse_geom_name.empty());
   TBOX_ASSERT(coarse_geom_name != getObjectName());
   TBOX_ASSERT(coarsen_ratio > hier::IntVector::getZero(getDim()));

   hier::BoxContainer coarse_domain(getPhysicalDomain());
   coarse_domain.coarsen(coarsen_ratio);

   /*
    * Need to check that domain can be coarsened by given ratio.
    */
   const hier::BoxContainer& fine_domain = getPhysicalDomain();
   const int nboxes = fine_domain.size();
   hier::BoxContainer::const_iterator fine_domain_itr = fine_domain.begin();
   hier::BoxContainer::iterator coarse_domain_itr = coarse_domain.begin();
   for (int ib = 0; ib < nboxes; ++ib, ++fine_domain_itr, ++coarse_domain_itr) {
      hier::Box testbox =
         hier::Box::refine(*coarse_domain_itr,
                           coarsen_ratio);
      if (!testbox.isSpatiallyEqual(*fine_domain_itr)) {
#ifdef DEBUG_CHECK_ASSERTIONS
         tbox::plog
         << "CartesianGridGeometry::makeCoarsenedGridGeometry : Box # "
         << ib << std::endl;
         tbox::plog << "      fine box = " << *fine_domain_itr << std::endl;
         tbox::plog << "      coarse box = " << *coarse_domain_itr << std::endl;
         tbox::plog << "      refined coarse box = " << testbox << std::endl;
#endif
         TBOX_ERROR(
            "geom::CartesianGridGeometry::makeCoarsenedGridGeometry() error...\n"
            << "    geometry object with name = " << getObjectName()
            << "\n    Cannot be coarsened by ratio " << coarsen_ratio
            << std::endl);
      }
   }

   std::shared_ptr<hier::BaseGridGeometry> coarse_geometry(
      new CartesianGridGeometry(coarse_geom_name,
         d_x_lo,
         d_x_up,
         coarse_domain,
         d_transfer_operator_registry));

   coarse_geometry->initializePeriodicShift(getPeriodicShift(-coarsen_ratio));

   return coarse_geometry;
}

/*
 *************************************************************************
 *
 * Set data members for this geometry object based on arguments.
 *
 *************************************************************************
 */

void
CartesianGridGeometry::setGeometryData(
   const double* x_lo,
   const double* x_up,
   const hier::BoxContainer& domain)
{
   const tbox::Dimension& dim(getDim());

   TBOX_ASSERT(x_lo != 0);
   TBOX_ASSERT(x_up != 0);

   for (int id = 0; id < dim.getValue(); ++id) {
      d_x_lo[id] = x_lo[id];
      d_x_up[id] = x_up[id];
   }

   if (getPhysicalDomain().empty()) {
      setPhysicalDomain(domain, 1);
   }

   hier::Box bigbox(dim);
   const hier::BoxContainer& block_domain = getPhysicalDomain();
   for (hier::BoxContainer::const_iterator k = block_domain.begin();
        k != block_domain.end(); ++k) {
      bigbox += *k;
   }

   d_domain_box = bigbox;

   hier::IntVector ncells = d_domain_box.numberCells();
   for (int id2 = 0; id2 < dim.getValue(); ++id2) {
      double length = d_x_up[id2] - d_x_lo[id2];
      d_dx[id2] = length / ((double)ncells(id2));
   }
}

/*
 *************************************************************************
 *
 * Create CartesianPatchGeometry geometry object, initializing its
 * boundary and grid information and assign it to the given patch.
 *
 *************************************************************************
 */

void
CartesianGridGeometry::setGeometryDataOnPatch(
   hier::Patch& patch,
   const hier::IntVector& ratio_to_level_zero,
   const TwoDimBool& touches_regular_bdry) const
{
   const tbox::Dimension& dim(getDim());

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY2(dim, patch, ratio_to_level_zero);

   const hier::BlockId& block_id = patch.getBox().getBlockId();
   hier::BlockId::block_t blk = block_id.getBlockValue();

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio must be nonzero.  Additionally,
    * all components not equal to 1 must have the same sign.
    */
   TBOX_ASSERT(ratio_to_level_zero != 0);

   if (dim > tbox::Dimension(1)) {
      for (unsigned int i = 0; i < dim.getValue(); ++i) {
         bool pos0 = ratio_to_level_zero(blk,i) > 0;
         bool pos1 = ratio_to_level_zero(blk,(i + 1) % d_dim.getValue()) > 0;
         TBOX_ASSERT(pos0 == pos1
            || (ratio_to_level_zero(blk,i) == 1)
            || (ratio_to_level_zero(blk,(i + 1) % dim.getValue()) == 1));
      }
   }
#endif

   double dx[SAMRAI::MAX_DIM_VAL];
   double x_lo[SAMRAI::MAX_DIM_VAL];
   double x_up[SAMRAI::MAX_DIM_VAL];

   bool coarsen = false;
   if (ratio_to_level_zero(blk,0) < 0) coarsen = true;
   hier::IntVector tmp_rat = ratio_to_level_zero;
   for (int id2 = 0; id2 < dim.getValue(); ++id2) {
      tmp_rat(blk,id2) = abs(ratio_to_level_zero(blk,id2));
   }

   hier::Box index_box = d_domain_box;
   hier::Box box = patch.getBox();

   if (coarsen) {
      index_box.coarsen(tmp_rat);
      for (tbox::Dimension::dir_t id3 = 0; id3 < dim.getValue(); ++id3) {
         dx[id3] = d_dx[id3] * ((double)tmp_rat(id3));
      }
   } else {
      index_box.refine(tmp_rat);
      for (tbox::Dimension::dir_t id4 = 0; id4 < dim.getValue(); ++id4) {
         dx[id4] = d_dx[id4] / ((double)tmp_rat(id4));
      }
   }

   for (tbox::Dimension::dir_t id5 = 0; id5 < dim.getValue(); ++id5) {
      x_lo[id5] = d_x_lo[id5]
         + ((double)(box.lower(id5) - index_box.lower(id5))) * dx[id5];
      x_up[id5] = x_lo[id5] + ((double)box.numberCells(id5)) * dx[id5];
   }

   std::shared_ptr<CartesianPatchGeometry> geom(
      std::make_shared<CartesianPatchGeometry>(ratio_to_level_zero,
         touches_regular_bdry,
         block_id,
         dx, x_lo, x_up));

   patch.setPatchGeometry(geom);

}

void
CartesianGridGeometry::buildOperators()
{
   // CartesianGridGeometry specific Coarsening Operators
   addCoarsenOperator(
      typeid(pdat::CellVariable<dcomplex>).name(),
      std::make_shared<CartesianCellComplexWeightedAverage>());
   addCoarsenOperator(
      typeid(pdat::CellVariable<double>).name(),
      std::make_shared<CartesianCellDoubleWeightedAverage>());
   addCoarsenOperator(
      typeid(pdat::CellVariable<float>).name(),
      std::make_shared<CartesianCellFloatWeightedAverage>());
   addCoarsenOperator(
      typeid(pdat::EdgeVariable<dcomplex>).name(),
      std::make_shared<CartesianEdgeComplexWeightedAverage>());
   addCoarsenOperator(
      typeid(pdat::EdgeVariable<double>).name(),
      std::make_shared<CartesianEdgeDoubleWeightedAverage>());
   addCoarsenOperator(
      typeid(pdat::EdgeVariable<float>).name(),
      std::make_shared<CartesianEdgeFloatWeightedAverage>());
   addCoarsenOperator(
      typeid(pdat::FaceVariable<dcomplex>).name(),
      std::make_shared<CartesianFaceComplexWeightedAverage>());
   addCoarsenOperator(
      typeid(pdat::FaceVariable<double>).name(),
      std::make_shared<CartesianFaceDoubleWeightedAverage>());
   addCoarsenOperator(
      typeid(pdat::FaceVariable<float>).name(),
      std::make_shared<CartesianFaceFloatWeightedAverage>());
   addCoarsenOperator(
      typeid(pdat::OuterfaceVariable<dcomplex>).name(),
      std::make_shared<CartesianOuterfaceComplexWeightedAverage>());
   addCoarsenOperator(
      typeid(pdat::OuterfaceVariable<double>).name(),
      std::make_shared<CartesianOuterfaceDoubleWeightedAverage>());
   addCoarsenOperator(
      typeid(pdat::OuterfaceVariable<float>).name(),
      std::make_shared<CartesianOuterfaceFloatWeightedAverage>());
   addCoarsenOperator(
      typeid(pdat::OutersideVariable<double>).name(),
      std::make_shared<CartesianOutersideDoubleWeightedAverage>());
   addCoarsenOperator(
      typeid(pdat::SideVariable<dcomplex>).name(),
      std::make_shared<CartesianSideComplexWeightedAverage>());
   addCoarsenOperator(
      typeid(pdat::SideVariable<double>).name(),
      std::make_shared<CartesianSideDoubleWeightedAverage>());
   addCoarsenOperator(
      typeid(pdat::SideVariable<float>).name(),
      std::make_shared<CartesianSideFloatWeightedAverage>());

   // CartesianGridGeometry specific Refinement Operators
   addRefineOperator(
      //typeid(pdat::CellVariable<dcomplex>).name(),
      //std::make_shared<CartesianCellComplexConservativeLinearRefine>());
      typeid(pdat::CellVariable<dcomplex>).name(),
      std::make_shared<CartesianCellConservativeLinearRefine<dcomplex>>());
   addRefineOperator(
      //typeid(pdat::CellVariable<double>).name(),
      //std::make_shared<CartesianCellDoubleConservativeLinearRefine>());
      typeid(pdat::CellVariable<double>).name(),
      std::make_shared<CartesianCellConservativeLinearRefine<double>>());
   addRefineOperator(
      //typeid(pdat::CellVariable<float>).name(),
      //std::make_shared<CartesianCellFloatConservativeLinearRefine>());
      typeid(pdat::CellVariable<float>).name(),
      std::make_shared<CartesianCellConservativeLinearRefine<float>>());
   addRefineOperator(
      typeid(pdat::EdgeVariable<double>).name(),
      std::make_shared<CartesianEdgeDoubleConservativeLinearRefine>());
   addRefineOperator(
      typeid(pdat::EdgeVariable<float>).name(),
      std::make_shared<CartesianEdgeFloatConservativeLinearRefine>());
   addRefineOperator(
      typeid(pdat::FaceVariable<double>).name(),
      std::make_shared<CartesianFaceDoubleConservativeLinearRefine>());
   addRefineOperator(
      typeid(pdat::FaceVariable<float>).name(),
      std::make_shared<CartesianFaceFloatConservativeLinearRefine>());
   addRefineOperator(
      typeid(pdat::SideVariable<double>).name(),
      std::make_shared<CartesianSideDoubleConservativeLinearRefine>());
   addRefineOperator(
      typeid(pdat::SideVariable<float>).name(),
      std::make_shared<CartesianSideFloatConservativeLinearRefine>());
   addRefineOperator(
      typeid(pdat::CellVariable<dcomplex>).name(),
      std::make_shared<CartesianCellComplexLinearRefine>());
   addRefineOperator(
      typeid(pdat::CellVariable<double>).name(),
      std::make_shared<CartesianCellDoubleLinearRefine>());
   addRefineOperator(
      typeid(pdat::CellVariable<float>).name(),
      std::make_shared<CartesianCellFloatLinearRefine>());
   addRefineOperator(
      typeid(pdat::NodeVariable<dcomplex>).name(),
      std::make_shared<CartesianNodeComplexLinearRefine>());
   addRefineOperator(
      typeid(pdat::NodeVariable<double>).name(),
      std::make_shared<CartesianNodeDoubleLinearRefine>());
   addRefineOperator(
      typeid(pdat::NodeVariable<float>).name(),
      std::make_shared<CartesianNodeFloatLinearRefine>());
}

/*
 *************************************************************************
 *
 * Print CartesianGridGeometry class data.
 *
 *************************************************************************
 */

void
CartesianGridGeometry::printClassData(
   std::ostream& os) const
{
   const tbox::Dimension& dim(getDim());

   os << "Printing CartesianGridGeometry data: this = "
      << (CartesianGridGeometry *)this << std::endl;
   os << "d_object_name = " << getObjectName() << std::endl;

   int id;
   os << "d_x_lo = ";
   for (id = 0; id < dim.getValue(); ++id) {
      os << d_x_lo[id] << "   ";
   }
   os << std::endl;
   os << "d_x_up = ";
   for (id = 0; id < dim.getValue(); ++id) {
      os << d_x_up[id] << "   ";
   }
   os << std::endl;
   os << "d_dx = ";
   for (id = 0; id < dim.getValue(); ++id) {
      os << d_dx[id] << "   ";
   }
   os << std::endl;

   os << "d_domain_box = " << d_domain_box << std::endl;

   GridGeometry::printClassData(os);
}

/*
 *************************************************************************
 *
 * Write class version number and object state to restart database.
 *
 *************************************************************************
 */

void
CartesianGridGeometry::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   hier::BaseGridGeometry::putToRestart(restart_db);

   const tbox::Dimension& dim(getDim());

   restart_db->putInteger("GEOM_CARTESIAN_GRID_GEOMETRY_VERSION",
      GEOM_CARTESIAN_GRID_GEOMETRY_VERSION);

   restart_db->putDoubleArray("x_lo", d_x_lo, dim.getValue());
   restart_db->putDoubleArray("x_up", d_x_up, dim.getValue());

}

/*
 *************************************************************************
 *
 * Data is read from input only if the simulation is not from restart.
 * Otherwise, all values specifed in the input database are ignored.
 * In this method data from the database are read to local
 * variables and the setGeometryData() method is called to
 * initialize the data members.
 *
 *************************************************************************
 */

void
CartesianGridGeometry::getFromInput(
   const std::shared_ptr<tbox::Database>& input_db,
   bool is_from_restart)
{
   if (!is_from_restart && !input_db) {
      TBOX_ERROR(": CartesianGridGeometry::getFromInput()\n"
         << "no input database supplied" << std::endl);
   }

   const tbox::Dimension& dim(getDim());

   if (!is_from_restart) {

      double x_lo[SAMRAI::MAX_DIM_VAL],
             x_up[SAMRAI::MAX_DIM_VAL];

      input_db->getDoubleArray("x_lo", x_lo, dim.getValue());
      input_db->getDoubleArray("x_up", x_up, dim.getValue());
      for (int i = 0; i < dim.getValue(); ++i) {
         if (!(x_lo[i] < x_up[i])) {
            INPUT_RANGE_ERROR("x_lo and x_up");
         }
      }

      setGeometryData(x_lo, x_up, getPhysicalDomain());

   } else if (input_db) {
      bool read_on_restart =
         input_db->getBoolWithDefault("read_on_restart", false);
      int num_keys = static_cast<int>(input_db->getAllKeys().size());
      if (num_keys > 0 && read_on_restart) {
         TBOX_WARNING(
            "CartesianGridGeometry::getFromInput() warning...\n"
            << "You want to override restart data with values from\n"
            << "an input database which is not allowed." << std::endl);
      }
   }
}

/*
 *************************************************************************
 *
 * Checks to see if the version number for the class is the same as
 * as the version number of the restart file.
 * If they are equal, then the data from the database are read to local
 * variables and the setGeometryData() method is called to
 * initialize the data members.
 *
 *************************************************************************
 */
void
CartesianGridGeometry::getFromRestart()
{
   std::shared_ptr<tbox::Database> restart_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!restart_db->isDatabase(getObjectName())) {
      TBOX_ERROR("CartesianGridGeometry::getFromRestart() error...\n"
         << "    database with name " << getObjectName()
         << " not found in the restart file" << std::endl);
   }
   std::shared_ptr<tbox::Database> db(
      restart_db->getDatabase(getObjectName()));

   const tbox::Dimension& dim(getDim());

   int ver = db->getInteger("GEOM_CARTESIAN_GRID_GEOMETRY_VERSION");
   if (ver != GEOM_CARTESIAN_GRID_GEOMETRY_VERSION) {
      TBOX_ERROR("CartesianGridGeometry::getFromRestart() error...\n"
         << "    geometry object with name = " << getObjectName()
         << "Restart file version is different than class version"
         << std::endl);
   }

   double x_lo[SAMRAI::MAX_DIM_VAL],
          x_up[SAMRAI::MAX_DIM_VAL];
   db->getDoubleArray("x_lo", x_lo, dim.getValue());
   db->getDoubleArray("x_up", x_up, dim.getValue());

   setGeometryData(x_lo, x_up, getPhysicalDomain());

}

}
}
