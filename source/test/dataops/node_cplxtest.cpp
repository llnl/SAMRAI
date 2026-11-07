/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Main program to test node-centered complex patch data ops
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <memory>

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"

#include "SAMRAI/tbox/SAMRAIManager.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/math/HierarchyDataOpsComplex.h"
#include "SAMRAI/math/HierarchyNodeDataOpsComplex.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/math/HierarchyNodeDataOpsReal.h"
#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/pdat/NodeIterator.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/VariableContext.h"


using namespace SAMRAI;

/* Helper function declarations */
static bool
complexDataSameAsValue(
   int desc_id,
   dcomplex value,
   std::shared_ptr<hier::PatchHierarchy> hierarchy);
static bool
doubleDataSameAsValue(
   int desc_id,
   double value,
   std::shared_ptr<hier::PatchHierarchy> hierarchy);

#define NVARS 4

int main(
   int argc,
   char* argv[]) {

   int num_failures = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   if (argc < 2) {
      TBOX_ERROR("Usage: " << argv[0] << " [dimension]");
   }

   const unsigned short d = static_cast<unsigned short>(atoi(argv[1]));
   TBOX_ASSERT(d > 0);
   TBOX_ASSERT(d <= SAMRAI::MAX_DIM_VAL);
   const tbox::Dimension dim(d);

   const std::string log_fn = std::string("node_cplxtest.")
      + tbox::Utilities::intToString(dim.getValue(), 1) + "d.log";
   tbox::PIO::logAllNodes(log_fn);

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {
      int ln, iv;

      // Make a dummy hierarchy domain
      double lo[SAMRAI::MAX_DIM_VAL];
      double hi[SAMRAI::MAX_DIM_VAL];

      hier::Index clo0(dim);
      hier::Index chi0(dim);
      hier::Index clo1(dim);
      hier::Index chi1(dim);
      hier::Index flo0(dim);
      hier::Index fhi0(dim);
      hier::Index flo1(dim);
      hier::Index fhi1(dim);

      for (int i = 0; i < dim.getValue(); ++i) {
         lo[i] = 0.0;
         clo0(i) = 0;
         flo0(i) = 4;
         fhi0(i) = 7;
         if (i == 1) {
            hi[i] = 0.5;
            chi0(i) = 2;
            clo1(i) = 3;
            chi1(i) = 4;
         } else {
            hi[i] = 1.0;
            chi0(i) = 9;
            clo1(i) = 0;
            chi1(i) = 9;
         }
         if (i == 0) {
            flo1(i) = 8;
            fhi1(i) = 13;
         } else {
            flo1(i) = flo0(i);
            fhi1(i) = fhi0(i);
         }
      }

      hier::Box coarse0(clo0, chi0, hier::BlockId(0));
      hier::Box coarse1(clo1, chi1, hier::BlockId(0));
      hier::Box fine0(flo0, fhi0, hier::BlockId(0));
      hier::Box fine1(flo1, fhi1, hier::BlockId(0));
      hier::IntVector ratio(dim, 2);

      hier::BoxContainer coarse_domain;
      hier::BoxContainer fine_boxes;
      coarse_domain.pushBack(coarse0);
      coarse_domain.pushBack(coarse1);
      fine_boxes.pushBack(fine0);
      fine_boxes.pushBack(fine1);

      std::shared_ptr<geom::CartesianGridGeometry> geometry(
         new geom::CartesianGridGeometry(
            "CartesianGeometry",
            lo,
            hi,
            coarse_domain));

      std::shared_ptr<hier::PatchHierarchy> hierarchy(
         new hier::PatchHierarchy("PatchHierarchy", geometry));

      hierarchy->setMaxNumberOfLevels(2);
      hierarchy->setRatioToCoarserLevel(ratio, 1);

      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
      const int nproc = mpi.getSize();

      const int n_coarse_boxes = coarse_domain.size();
      const int n_fine_boxes = fine_boxes.size();

      std::shared_ptr<hier::BoxLevel> layer0(
         std::make_shared<hier::BoxLevel>(
            hier::IntVector(dim, 1), geometry));
      std::shared_ptr<hier::BoxLevel> layer1(
         std::make_shared<hier::BoxLevel>(ratio, geometry));

      hier::BoxContainer::iterator coarse_itr = coarse_domain.begin();
      for (int ib = 0; ib < n_coarse_boxes; ++ib, ++coarse_itr) {
         if (nproc > 1) {
            if (ib == layer0->getMPI().getRank()) {
               layer0->addBox(hier::Box(*coarse_itr, hier::LocalId(ib),
                     layer0->getMPI().getRank()));
            }
         } else {
            layer0->addBox(hier::Box(*coarse_itr, hier::LocalId(ib), 0));
         }
      }

      hier::BoxContainer::iterator fine_itr = fine_boxes.begin();
      for (int ib = 0; ib < n_fine_boxes; ++ib, ++fine_itr) {
         if (nproc > 1) {
            if (ib == layer1->getMPI().getRank()) {
               layer1->addBox(hier::Box(*fine_itr, hier::LocalId(ib),
                     layer1->getMPI().getRank()));
            }
         } else {
            layer1->addBox(hier::Box(*fine_itr, hier::LocalId(ib), 0));
         }
      }

      hierarchy->makeNewPatchLevel(0, layer0);
      hierarchy->makeNewPatchLevel(1, layer1);

      // Create instance of hier::Variable database
      hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
      std::shared_ptr<hier::VariableContext> dummy(
         variable_db->getContext("dummy"));
      const hier::IntVector no_ghosts(dim, 0);

      // Make some dummy variables and data on the hierarchy
      std::shared_ptr<pdat::NodeVariable<dcomplex> > nvar[NVARS];
      int nvindx[NVARS];
      nvar[0].reset(new pdat::NodeVariable<dcomplex>(dim, "nvar0", 1));
      nvindx[0] = variable_db->registerVariableAndContext(
            nvar[0], dummy, no_ghosts);
      nvar[1].reset(new pdat::NodeVariable<dcomplex>(dim, "nvar1", 1));
      nvindx[1] = variable_db->registerVariableAndContext(
            nvar[1], dummy, no_ghosts);
      nvar[2].reset(new pdat::NodeVariable<dcomplex>(dim, "nvar2", 1));
      nvindx[2] = variable_db->registerVariableAndContext(
            nvar[2], dummy, no_ghosts);
      nvar[3].reset(new pdat::NodeVariable<dcomplex>(dim, "nvar3", 1));
      nvindx[3] = variable_db->registerVariableAndContext(
            nvar[3], dummy, no_ghosts);

      std::shared_ptr<pdat::NodeVariable<double> > nwgt(
         new pdat::NodeVariable<double>(dim, "nwgt", 1));
      int nwgt_id = variable_db->registerVariableAndContext(
            nwgt, dummy, no_ghosts);

      // allocate data on hierarchy
      for (ln = 0; ln < 2; ++ln) {
         hierarchy->getPatchLevel(ln)->allocatePatchData(nwgt_id);
         for (iv = 0; iv < NVARS; ++iv) {
            hierarchy->getPatchLevel(ln)->allocatePatchData(nvindx[iv]);
         }
      }

      std::shared_ptr<math::HierarchyDataOpsComplex> node_ops(
         new math::HierarchyNodeDataOpsComplex(
            hierarchy,
            0,
            1));
      TBOX_ASSERT(node_ops);

      std::shared_ptr<math::HierarchyDataOpsReal<double> > nwgt_ops(
         new math::HierarchyNodeDataOpsReal<double>(
            hierarchy,
            0,
            1));

      std::shared_ptr<hier::Patch> patch;

      // Initialize control volume data for node-centered components
      hier::Box coarse_fine = fine0 + fine1;
      coarse_fine.coarsen(ratio);
      for (ln = 0; ln < 2; ++ln) {
         std::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));
         for (hier::PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            patch = *ip;
            std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
               SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                  patch->getPatchGeometry()));
            TBOX_ASSERT(pgeom);
            const double* dx = pgeom->getDx();
            double node_vol = dx[0];
            for (int i = 1; i < dim.getValue(); ++i) {
               node_vol *= dx[i];
            }
            std::shared_ptr<pdat::NodeData<double> > data(
               SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
                  patch->getPatchData(nwgt_id)));
            TBOX_ASSERT(data);
            data->fillAll(node_vol);
#if defined(HAVE_RAJA)
            tbox::parallel_synchronize();
#endif
            pdat::NodeIndex ni(dim);

            if (dim.getValue() == 2) {
               int plo0 = patch->getBox().lower(0);
               int phi0 = patch->getBox().upper(0);
               int plo1 = patch->getBox().lower(1);
               int phi1 = patch->getBox().upper(1);
               int ic;

               if (ln == 0) {
                  data->fillAll(0.0, (coarse_fine * patch->getBox()));
#if defined(HAVE_RAJA)
                  tbox::parallel_synchronize();
#endif

                  if (patch->getLocalId() == 0) {
                     // bottom face boundaries
                     for (ic = plo0; ic < phi0; ++ic) {
                        ni = pdat::NodeIndex(hier::Index(ic, plo1),
                              pdat::NodeIndex::LowerRight);
                        (*data)(ni) *= 0.5;
                     }
                     // left and right face boundaries
                     for (ic = plo1; ic <= phi1; ++ic) {
                        ni = pdat::NodeIndex(hier::Index(plo0, ic),
                              pdat::NodeIndex::UpperLeft);
                        (*data)(ni) *= 0.5;
                        ni = pdat::NodeIndex(hier::Index(phi0, ic),
                              pdat::NodeIndex::UpperRight);
                        (*data)(ni) *= 0.5;
                     }
                     // corner boundaries
                     (*data)(pdat::NodeIndex(hier::Index(plo0, plo1),
                                pdat::NodeIndex::LowerLeft)) *= 0.25;
                     (*data)(pdat::NodeIndex(hier::Index(phi0, plo1),
                                pdat::NodeIndex::LowerRight)) *= 0.25;
                  } else {
                     // top and bottom face boundaries
                     for (ic = plo0; ic < phi0; ++ic) {
                        ni = pdat::NodeIndex(hier::Index(ic, phi1),
                              pdat::NodeIndex::UpperRight);
                        (*data)(ni) *= 0.5;
                        ni = pdat::NodeIndex(hier::Index(ic, plo1),
                              pdat::NodeIndex::LowerRight);
                        (*data)(ni) = 0.0;
                     }
                     // left and right face boundaries
                     for (ic = plo1; ic < phi1; ++ic) {
                        ni = pdat::NodeIndex(hier::Index(plo0, ic),
                              pdat::NodeIndex::UpperLeft);
                        (*data)(ni) *= 0.5;
                        ni = pdat::NodeIndex(hier::Index(phi0, ic),
                              pdat::NodeIndex::UpperRight);
                        (*data)(ni) *= 0.5;
                     }
                     // corner boundaries
                     (*data)(pdat::NodeIndex(hier::Index(plo0, plo1),
                                pdat::NodeIndex::LowerLeft)) = 0.0;
                     (*data)(pdat::NodeIndex(hier::Index(plo0, phi1),
                                pdat::NodeIndex::UpperLeft)) *= 0.25;
                     (*data)(pdat::NodeIndex(hier::Index(phi0, plo1),
                                pdat::NodeIndex::LowerRight)) = 0.0;
                     (*data)(pdat::NodeIndex(hier::Index(phi0, phi1),
                                pdat::NodeIndex::UpperRight)) *= 0.25;
                  }
               } else {
                  if (patch->getLocalId() == 0) {
                     // top and bottom coarse-fine face boundaries
                     for (ic = plo0; ic <= phi0; ++ic) {
                        ni = pdat::NodeIndex(hier::Index(ic, plo1),
                              pdat::NodeIndex::LowerRight);
                        (*data)(ni) *= 1.5;
                        ni = pdat::NodeIndex(hier::Index(ic, phi1),
                              pdat::NodeIndex::UpperRight);
                        (*data)(ni) *= 1.5;
                     }
                     // left coarse-fine face boundaries
                     for (ic = plo1; ic < phi1; ++ic) {
                        ni = pdat::NodeIndex(hier::Index(plo0, ic),
                              pdat::NodeIndex::UpperLeft);
                        (*data)(ni) *= 1.5;
                     }
                     // coarse-fine corner boundaries
                     (*data)(pdat::NodeIndex(hier::Index(plo0, plo1),
                                pdat::NodeIndex::LowerLeft)) *= 2.25;
                     (*data)(pdat::NodeIndex(hier::Index(plo0, phi1),
                                pdat::NodeIndex::UpperLeft)) *= 2.25;
                  } else {
                     // top and bottom coarse-fine face boundaries
                     for (ic = plo0; ic < phi0; ++ic) {
                        ni = pdat::NodeIndex(hier::Index(ic, plo1),
                              pdat::NodeIndex::LowerRight);
                        (*data)(ni) *= 1.5;
                        ni = pdat::NodeIndex(hier::Index(ic, phi1),
                              pdat::NodeIndex::UpperRight);
                        (*data)(ni) *= 1.5;
                     }
                     // right coarse-fine face boundaries
                     for (ic = plo1; ic < phi1; ++ic) {
                        ni = pdat::NodeIndex(hier::Index(phi0, ic),
                              pdat::NodeIndex::UpperRight);
                        (*data)(ni) *= 1.5;
                     }
                     // coarse-fine corner boundaries
                     (*data)(pdat::NodeIndex(hier::Index(phi0, plo1),
                                pdat::NodeIndex::LowerRight)) *= 2.25;
                     (*data)(pdat::NodeIndex(hier::Index(phi0, phi1),
                                pdat::NodeIndex::UpperRight)) *= 2.25;
                     // shared left boundaries
                     for (ic = plo1; ic <= phi1 + 1; ++ic) {
                        ni = pdat::NodeIndex(hier::Index(plo0, ic),
                              pdat::NodeIndex::LowerLeft);
                        (*data)(ni) = 0;
                     }
                  }
               }
            } else {
               int plo0 = patch->getBox().lower(0);
               int phi0 = patch->getBox().upper(0);
               int plo1 = patch->getBox().lower(1);
               int phi1 = patch->getBox().upper(1);
               int plo2 = patch->getBox().lower(2);
               int phi2 = patch->getBox().upper(2);
               int ic0, ic1, ic2;

               if (ln == 0) {
                  data->fillAll(0.0, (coarse_fine * patch->getBox()));
#if defined(HAVE_RAJA)
                  tbox::parallel_synchronize();
#endif

                  if (patch->getLocalId() == 0) {
                     // front and back face boundaries
                     for (ic0 = plo0; ic0 < phi0; ++ic0) {
                        for (ic1 = plo1; ic1 <= phi1; ++ic1) {
                           ni = pdat::NodeIndex(hier::Index(ic0, ic1, plo2),
                                 pdat::NodeIndex::UUL);
                           (*data)(ni) *= 0.5;
                           ni = pdat::NodeIndex(hier::Index(ic0, ic1, phi2),
                                 pdat::NodeIndex::UUU);
                           (*data)(ni) *= 0.5;
                        }
                     }
                     // bottom face boundary
                     for (ic0 = plo0; ic0 < phi0; ++ic0) {
                        for (ic2 = plo2; ic2 < phi2; ++ic2) {
                           ni = pdat::NodeIndex(hier::Index(ic0, plo1, ic2),
                                 pdat::NodeIndex::ULU);
                           (*data)(ni) *= 0.5;
                        }
                     }
                     // left and right face boundaries
                     for (ic1 = plo1; ic1 <= phi1; ++ic1) {
                        for (ic2 = plo2; ic2 < phi2; ++ic2) {
                           ni = pdat::NodeIndex(hier::Index(plo0, ic1, ic2),
                                 pdat::NodeIndex::LUU);
                           (*data)(ni) *= 0.5;
                           ni = pdat::NodeIndex(hier::Index(phi0, ic1, ic2),
                                 pdat::NodeIndex::UUU);
                           (*data)(ni) *= 0.5;
                        }
                     }
                     // lower front and back edge boundaries
                     for (ic0 = plo0; ic0 < phi0; ++ic0) {
                        ni = pdat::NodeIndex(hier::Index(ic0, plo1, plo2),
                              pdat::NodeIndex::ULL);
                        (*data)(ni) *= 0.25;
                        ni = pdat::NodeIndex(hier::Index(ic0, plo1, phi2),
                              pdat::NodeIndex::ULU);
                        (*data)(ni) *= 0.25;
                     }
                     // lower left and right edge boundaries
                     for (ic2 = plo2; ic2 < phi2; ++ic2) {
                        ni = pdat::NodeIndex(hier::Index(plo0, plo1, ic2),
                              pdat::NodeIndex::LLU);
                        (*data)(ni) *= 0.25;
                        ni = pdat::NodeIndex(hier::Index(phi0, plo1, ic2),
                              pdat::NodeIndex::ULU);
                        (*data)(ni) *= 0.25;
                     }
                     // left and right front and back edge boundaries
                     for (ic1 = plo1; ic1 <= phi1; ++ic1) {
                        ni = pdat::NodeIndex(hier::Index(plo0, ic1, plo2),
                              pdat::NodeIndex::LUL);
                        (*data)(ni) *= 0.25;
                        ni = pdat::NodeIndex(hier::Index(phi0, ic1, plo2),
                              pdat::NodeIndex::UUL);
                        (*data)(ni) *= 0.25;
                        ni = pdat::NodeIndex(hier::Index(plo0, ic1, phi2),
                              pdat::NodeIndex::LUU);
                        (*data)(ni) *= 0.25;
                        ni = pdat::NodeIndex(hier::Index(phi0, ic1, phi2),
                              pdat::NodeIndex::UUU);
                        (*data)(ni) *= 0.25;
                     }
                     // corner boundaries
                     (*data)(pdat::NodeIndex(hier::Index(plo0, plo1, plo2),
                                pdat::NodeIndex::LLL)) *= 0.125;
                     (*data)(pdat::NodeIndex(hier::Index(phi0, plo1, plo2),
                                pdat::NodeIndex::ULL)) *= 0.125;
                     (*data)(pdat::NodeIndex(hier::Index(plo0, plo1, phi2),
                                pdat::NodeIndex::LLU)) *= 0.125;
                     (*data)(pdat::NodeIndex(hier::Index(phi0, plo1, phi2),
                                pdat::NodeIndex::ULU)) *= 0.125;
                  } else {
                     // front and back face boundaries
                     for (ic0 = plo0; ic0 < phi0; ++ic0) {
                        for (ic1 = plo1; ic1 < phi1; ++ic1) {
                           ni = pdat::NodeIndex(hier::Index(ic0, ic1, plo2),
                                 pdat::NodeIndex::UUL);
                           (*data)(ni) *= 0.5;
                           ni = pdat::NodeIndex(hier::Index(ic0, ic1, phi2),
                                 pdat::NodeIndex::UUU);
                           (*data)(ni) *= 0.5;
                        }
                     }
                     // top face boundary
                     for (ic0 = plo0; ic0 < phi0; ++ic0) {
                        for (ic2 = plo2; ic2 < phi2; ++ic2) {
                           ni = pdat::NodeIndex(hier::Index(ic0, phi1, ic2),
                                 pdat::NodeIndex::UUU);
                           (*data)(ni) *= 0.5;
                        }
                     }
                     // left and right face boundaries
                     for (ic1 = plo1; ic1 < phi1; ++ic1) {
                        for (ic2 = plo2; ic2 < phi2; ++ic2) {
                           ni = pdat::NodeIndex(hier::Index(plo0, ic1, ic2),
                                 pdat::NodeIndex::LUU);
                           (*data)(ni) *= 0.5;
                           ni = pdat::NodeIndex(hier::Index(phi0, ic1, ic2),
                                 pdat::NodeIndex::UUU);
                           (*data)(ni) *= 0.5;
                        }
                     }
                     // upper and lower front and back edge boundaries
                     for (ic0 = plo0; ic0 < phi0; ++ic0) {
                        ni = pdat::NodeIndex(hier::Index(ic0, phi1, plo2),
                              pdat::NodeIndex::UUL);
                        (*data)(ni) *= 0.25;
                        ni = pdat::NodeIndex(hier::Index(ic0, phi1, phi2),
                              pdat::NodeIndex::UUU);
                        (*data)(ni) *= 0.25;
                        ni = pdat::NodeIndex(hier::Index(ic0, plo1, plo2),
                              pdat::NodeIndex::ULL);
                        (*data)(ni) = 0.0;
                        ni = pdat::NodeIndex(hier::Index(ic0, plo1, phi2),
                              pdat::NodeIndex::ULU);
                        (*data)(ni) = 0.0;
                     }
                     // upper and lower left and right edge boundaries
                     for (ic2 = plo2; ic2 < phi2; ++ic2) {
                        ni = pdat::NodeIndex(hier::Index(plo0, phi1, ic2),
                              pdat::NodeIndex::LUU);
                        (*data)(ni) *= 0.25;
                        ni = pdat::NodeIndex(hier::Index(phi0, phi1, ic2),
                              pdat::NodeIndex::UUU);
                        (*data)(ni) *= 0.25;
                        ni = pdat::NodeIndex(hier::Index(plo0, plo1, ic2),
                              pdat::NodeIndex::LLU);
                        (*data)(ni) = 0.0;
                        ni = pdat::NodeIndex(hier::Index(phi0, plo1, ic2),
                              pdat::NodeIndex::ULU);
                        (*data)(ni) = 0.0;
                     }
                     // front and back left and right edge boundaries
                     for (ic1 = plo1; ic1 < phi1; ++ic1) {
                        ni = pdat::NodeIndex(hier::Index(plo0, ic1, plo2),
                              pdat::NodeIndex::LUL);
                        (*data)(ni) *= 0.25;
                        ni = pdat::NodeIndex(hier::Index(phi0, ic1, plo2),
                              pdat::NodeIndex::UUL);
                        (*data)(ni) *= 0.25;
                        ni = pdat::NodeIndex(hier::Index(plo0, ic1, phi2),
                              pdat::NodeIndex::LUU);
                        (*data)(ni) *= 0.25;
                        ni = pdat::NodeIndex(hier::Index(phi0, ic1, phi2),
                              pdat::NodeIndex::UUU);
                        (*data)(ni) *= 0.25;
                     }
                     // corner boundaries
                     (*data)(pdat::NodeIndex(hier::Index(plo0, plo1, plo2),
                                pdat::NodeIndex::LLL)) = 0.0;
                     (*data)(pdat::NodeIndex(hier::Index(phi0, plo1, plo2),
                                pdat::NodeIndex::ULL)) = 0.0;
                     (*data)(pdat::NodeIndex(hier::Index(plo0, plo1, phi2),
                                pdat::NodeIndex::LLU)) = 0.0;
                     (*data)(pdat::NodeIndex(hier::Index(phi0, plo1, phi2),
                                pdat::NodeIndex::ULU)) = 0.0;
                     (*data)(pdat::NodeIndex(hier::Index(plo0, phi1, plo2),
                                pdat::NodeIndex::LUL)) *= 0.125;
                     (*data)(pdat::NodeIndex(hier::Index(phi0, phi1, plo2),
                                pdat::NodeIndex::UUL)) *= 0.125;
                     (*data)(pdat::NodeIndex(hier::Index(plo0, phi1, phi2),
                                pdat::NodeIndex::LUU)) *= 0.125;
                     (*data)(pdat::NodeIndex(hier::Index(phi0, phi1, phi2),
                                pdat::NodeIndex::UUU)) *= 0.125;
                     // bottom face boundary
                     for (ic0 = plo0; ic0 < phi0; ++ic0) {
                        for (ic2 = plo2; ic2 < phi2; ++ic2) {
                           ni = pdat::NodeIndex(hier::Index(ic0, plo1, ic2),
                                 pdat::NodeIndex::ULU);
                           (*data)(ni) = 0.0;
                        }
                     }
                  }
               } else {
                  if (patch->getLocalId() == 0) {
                     // front and back coarse-fine face boundaries
                     for (ic0 = plo0; ic0 <= phi0; ++ic0) {
                        for (ic1 = plo1; ic1 < phi1; ++ic1) {
                           ni = pdat::NodeIndex(hier::Index(ic0, ic1, plo2),
                                 pdat::NodeIndex::UUL);
                           (*data)(ni) *= 1.5;
                           ni = pdat::NodeIndex(hier::Index(ic0, ic1, phi2),
                                 pdat::NodeIndex::UUU);
                           (*data)(ni) *= 1.5;
                        }
                     }
                     // top and bottom coarse-fine face boundaries
                     for (ic0 = plo0; ic0 <= phi0; ++ic0) {
                        for (ic2 = plo2; ic2 < phi2; ++ic2) {
                           ni = pdat::NodeIndex(hier::Index(ic0, plo1, ic2),
                                 pdat::NodeIndex::ULU);
                           (*data)(ni) *= 1.5;
                           ni = pdat::NodeIndex(hier::Index(ic0, phi1, ic2),
                                 pdat::NodeIndex::UUU);
                           (*data)(ni) *= 1.5;
                        }
                     }
                     // left coarse-fine face boundary
                     for (ic1 = plo1; ic1 < phi1; ++ic1) {
                        for (ic2 = plo2; ic2 < phi2; ++ic2) {
                           ni = pdat::NodeIndex(hier::Index(plo0, ic1, ic2),
                                 pdat::NodeIndex::LUU);
                           (*data)(ni) *= 1.5;
                        }
                     }
                     // upper and lower front and back edge boundaries
                     for (ic0 = plo0; ic0 <= phi0; ++ic0) {
                        ni = pdat::NodeIndex(hier::Index(ic0, phi1, plo2),
                              pdat::NodeIndex::UUL);
                        (*data)(ni) *= 2.25;
                        ni = pdat::NodeIndex(hier::Index(ic0, phi1, phi2),
                              pdat::NodeIndex::UUU);
                        (*data)(ni) *= 2.25;
                        ni = pdat::NodeIndex(hier::Index(ic0, plo1, plo2),
                              pdat::NodeIndex::ULL);
                        (*data)(ni) *= 2.25;
                        ni = pdat::NodeIndex(hier::Index(ic0, plo1, phi2),
                              pdat::NodeIndex::ULU);
                        (*data)(ni) *= 2.25;
                     }
                     // upper and lower left edge boundaries
                     for (ic2 = plo2; ic2 < phi2; ++ic2) {
                        ni = pdat::NodeIndex(hier::Index(plo0, phi1, ic2),
                              pdat::NodeIndex::LUU);
                        (*data)(ni) *= 2.25;
                        ni = pdat::NodeIndex(hier::Index(plo0, plo1, ic2),
                              pdat::NodeIndex::LLU);
                        (*data)(ni) *= 2.25;
                     }
                     // front and back left edge boundaries
                     for (ic1 = plo1; ic1 < phi1; ++ic1) {
                        ni = pdat::NodeIndex(hier::Index(plo0, ic1, plo2),
                              pdat::NodeIndex::LUL);
                        (*data)(ni) *= 2.25;
                        ni = pdat::NodeIndex(hier::Index(plo0, ic1, phi2),
                              pdat::NodeIndex::LUU);
                        (*data)(ni) *= 2.25;
                     }
                     // coarse-fine corner boundaries
                     (*data)(pdat::NodeIndex(hier::Index(plo0, plo1, plo2),
                                pdat::NodeIndex::LLL)) *= 3.375;
                     (*data)(pdat::NodeIndex(hier::Index(plo0, phi1, plo2),
                                pdat::NodeIndex::LUL)) *= 3.375;
                     (*data)(pdat::NodeIndex(hier::Index(plo0, plo1, phi2),
                                pdat::NodeIndex::LLU)) *= 3.375;
                     (*data)(pdat::NodeIndex(hier::Index(plo0, phi1, phi2),
                                pdat::NodeIndex::LUU)) *= 3.375;
                  } else {
                     // front and back coarse-fine face boundaries
                     for (ic0 = plo0; ic0 < phi0; ++ic0) {
                        for (ic1 = plo1; ic1 < phi1; ++ic1) {
                           ni = pdat::NodeIndex(hier::Index(ic0, ic1, plo2),
                                 pdat::NodeIndex::UUL);
                           (*data)(ni) *= 1.5;
                           ni = pdat::NodeIndex(hier::Index(ic0, ic1, phi2),
                                 pdat::NodeIndex::UUU);
                           (*data)(ni) *= 1.5;
                        }
                     }
                     // top and bottom coarse-fine face boundaries
                     for (ic0 = plo0; ic0 < phi0; ++ic0) {
                        for (ic2 = plo2; ic2 < phi2; ++ic2) {
                           ni = pdat::NodeIndex(hier::Index(ic0, plo1, ic2),
                                 pdat::NodeIndex::ULU);
                           (*data)(ni) *= 1.5;
                           ni = pdat::NodeIndex(hier::Index(ic0, phi1, ic2),
                                 pdat::NodeIndex::UUU);
                           (*data)(ni) *= 1.5;
                        }
                     }
                     // right coarse-fine face boundaries
                     for (ic1 = plo1; ic1 < phi1; ++ic1) {
                        for (ic2 = plo2; ic2 < phi2; ++ic2) {
                           ni = pdat::NodeIndex(hier::Index(phi0, ic1, ic2),
                                 pdat::NodeIndex::LUU);
                           (*data)(ni) *= 1.5;
                        }
                     }
                     // upper and lower front and back edge boundaries
                     for (ic0 = plo0; ic0 < phi0; ++ic0) {
                        ni = pdat::NodeIndex(hier::Index(ic0, phi1, plo2),
                              pdat::NodeIndex::UUL);
                        (*data)(ni) *= 2.25;
                        ni = pdat::NodeIndex(hier::Index(ic0, phi1, phi2),
                              pdat::NodeIndex::UUU);
                        (*data)(ni) *= 2.25;
                        ni = pdat::NodeIndex(hier::Index(ic0, plo1, plo2),
                              pdat::NodeIndex::ULL);
                        (*data)(ni) *= 2.25;
                        ni = pdat::NodeIndex(hier::Index(ic0, plo1, phi2),
                              pdat::NodeIndex::ULU);
                        (*data)(ni) *= 2.25;
                     }
                     // upper and lower right edge boundaries
                     for (ic2 = plo2; ic2 < phi2; ++ic2) {
                        ni = pdat::NodeIndex(hier::Index(phi0, phi1, ic2),
                              pdat::NodeIndex::UUU);
                        (*data)(ni) *= 2.25;
                        ni = pdat::NodeIndex(hier::Index(phi0, plo1, ic2),
                              pdat::NodeIndex::ULU);
                        (*data)(ni) *= 2.25;
                     }
                     // front and back right edge boundaries
                     for (ic1 = plo1; ic1 < phi1; ++ic1) {
                        ni = pdat::NodeIndex(hier::Index(phi0, ic1, plo2),
                              pdat::NodeIndex::UUL);
                        (*data)(ni) *= 2.25;
                        ni = pdat::NodeIndex(hier::Index(phi0, ic1, phi2),
                              pdat::NodeIndex::UUU);
                        (*data)(ni) *= 2.25;
                     }
                     // coarse-fine corner boundaries
                     (*data)(pdat::NodeIndex(hier::Index(phi0, plo1, plo2),
                                pdat::NodeIndex::ULL)) *= 3.375;
                     (*data)(pdat::NodeIndex(hier::Index(phi0, phi1, plo2),
                                pdat::NodeIndex::UUL)) *= 3.375;
                     (*data)(pdat::NodeIndex(hier::Index(phi0, plo1, phi2),
                                pdat::NodeIndex::ULU)) *= 3.375;
                     (*data)(pdat::NodeIndex(hier::Index(phi0, phi1, phi2),
                                pdat::NodeIndex::UUU)) *= 3.375;
                     // shared left boundaries
                     for (ic1 = plo1; ic1 <= phi1 + 1; ++ic1) {
                        for (ic2 = plo2; ic2 <= phi2 + 1; ++ic2) {
                           ni = pdat::NodeIndex(hier::Index(plo0, ic1, ic2),
                                 pdat::NodeIndex::LLL);
                           (*data)(ni) = 0;
                        }
                     }
                  }
               }
            }
         }
      }

      // Test #1: Print out control volume data and compute its integral

      // Test #1a: Check control volume data set properly
      // Expected: cwgt = 0.01 on coarse (except where finer patch exists) and
      // 0.0025 on fine level
/*   bool vol_test_passed = true;
 *   for (ln = 0; ln < 2; ++ln) {
 *   for (hier::PatchLevel::iterator ip(hierarchy->getPatchLevel(ln)->begin()); ip != hierarchy->getPatchLevel(ln)->end(); ++ip) {
 *   patch = hierarchy->getPatchLevel(ln)->getPatch(ip());
 *   std::shared_ptr< pdat::NodeData<double> > cvdata = patch->getPatchData(cwgt_id);
 *
 *   pdat::NodeIterator cend(pdat::NodeGeometry::end(cvdata->getBox()));
 *   for (pdat::NodeIterator c(pdat::NodeGeometry::begin(cvdata->getBox())); c != cend && vol_test_passed; ++c) {
 *   pdat::NodeIndex cell_index = *c;
 *
 *   if (ln == 0) {
 *   if ((coarse_fine * patch->getBox()).contains(cell_index)) {
 *   if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),0.0) ) {
 *   vol_test_passed = false;
 *   }
 *   } else {
 *   if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),0.01) ) {
 *   vol_test_passed = false;
 *   }
 *   }
 *   }
 *
 *   if (ln == 1) {
 *   if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),0.0025) ) {
 *   vol_test_passed = false;
 *   }
 *   }
 *   }
 *   }
 *   }
 *   if (!vol_test_passed) {
 *   ++num_failures;
 *   tbox::perr << "FAILED: - Test #1a: Check control volume data set properly" << std::endl;
 *   cwgt_ops->printData(cwgt_id, tbox::plog);
 *   }
 */
/*   tbox::plog << "node control volume data" << std::endl;
 *   nwgt_ops->printData(nwgt_id, tbox::plog);
 */

      // Test #1b: math::HierarchyNodeDataOpsComplex::sumControlVolumes()
      // Expected: norm = 0.5
      double norm =
         node_ops->sumControlVolumes(nvindx[0], nwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(norm, 0.5)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #1: math::HierarchyNodeDataOpsComplex::sumControlVolumes()\n"
         << "Expected value = 0.5 , Computed value = "
         << norm << std::endl;
      }

      // Test #2: math::HierarchyNodeDataOpsComplex::numberOfEntries()
      // Expected: num_data_points = 121
      size_t num_data_points = node_ops->numberOfEntries(nvindx[0]);
      {
         size_t compare;
         if (dim.getValue() == 2) {
            compare = 121;
         } else {
            compare = 1001;
         }
         if (num_data_points != compare) {
            ++num_failures;
            tbox::perr
            << "FAILED: - Test #2: math::HierarchyNodeDataOpsComplex::numberOfEntries()\n"
            << "Expected value = " << compare << ", Computed value = "
            << num_data_points << std::endl;
         }
      }

      // Test #3a: math::HierarchyNodeDataOpsComplex::setToScalar()
      // Expected: v0 = (2.0,1.5)
      dcomplex val0 = dcomplex(2.0, 1.5);
      node_ops->setToScalar(nvindx[0], val0);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!complexDataSameAsValue(nvindx[0], val0, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #3a: math::HierarchyNodeDataOpsComplex::setToScalar()\n"
         << "Expected: v0 = " << val0 << std::endl;
         node_ops->printData(nvindx[0], tbox::plog);
      }

      // Test #3b: math::HierarchyNodeDataOpsComplex::setToScalar()
      // Expected:  v1 = (4.0, 3.0)
      dcomplex val1(4.0, 3.0);
      node_ops->setToScalar(nvindx[1], dcomplex(4.0, 3.0));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!complexDataSameAsValue(nvindx[1], val1, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #3b: math::HierarchyNodeDataOpsComplex::setToScalar()\n"
         << "Expected: v1 = " << val1 << std::endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #4: math::HierarchyNodeDataOpsComplex::copyData()
      // Expected:   v2 = v1 = (4.0, 3.0)
      node_ops->copyData(nvindx[2], nvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!complexDataSameAsValue(nvindx[2], val1, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #4: math::HierarchyNodeDataOpsComplex::copyData()\n"
         << "Expected: v2 = " << val1 << std::endl;
         node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test #5: math::HierarchyNodeDataOpsComplex::swapData()
      // Expected:  v0 = (4.0, 3.0), v1 = (2.0,1.5)
      node_ops->swapData(nvindx[0], nvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!complexDataSameAsValue(nvindx[0], val1, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #5a: math::HierarchyNodeDataOpsComplex::swapData()\n"
         << "Expected: v0 = " << val1 << std::endl;
         node_ops->printData(nvindx[0], tbox::plog);
      }
      if (!complexDataSameAsValue(nvindx[1], val0, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #5b: math::HierarchyNodeDataOpsComplex::swapData()\n"
         << "Expected: v1 = " << val0 << std::endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #6: math::HierarchyNodeDataOpsComplex::scale()
      // Expected:  v2 = 0.25 * v2 = (1.0,0.75)
      node_ops->scale(nvindx[2], 0.25, nvindx[2]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_scale(1.0, 0.75);
      if (!complexDataSameAsValue(nvindx[2], val_scale, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #6: math::HierarchyNodeDataOpsComplex::scale()\n"
         << "Expected: v2 = " << val_scale << std::endl;
         node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test #7: math::HierarchyNodeDataOpsComplex::add()
      // Expected: v3 = v0 + v1 = (6.0, 4.5)
      node_ops->add(nvindx[3], nvindx[0], nvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_add(6.0, 4.5);
      if (!complexDataSameAsValue(nvindx[3], val_add, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #7: math::HierarchyNodeDataOpsComplex::add()\n"
         << "Expected: v3 = " << val_add << std::endl;
         node_ops->printData(nvindx[3], tbox::plog);
      }

      // Reset v0: v0 = (0.0,4.5)
      node_ops->setToScalar(nvindx[0], dcomplex(0.0, 4.5));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      // Test #8: math::HierarchyNodeDataOpsComplex::subtract()
      // Expected: v1 = v3 - v0 = (6.0,0.0)
      node_ops->subtract(nvindx[1], nvindx[3], nvindx[0]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_sub(6.0, 0.0);
      if (!complexDataSameAsValue(nvindx[1], val_sub, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #8: math::HierarchyNodeDataOpsComplex::subtract()\n"
         << "Expected: v1 = " << val_sub << std::endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #9a: math::HierarchyNodeDataOpsComplex::addScalar()
      // Expected: v1 = v1 + (0.0,-4.0) = (6.0,-4.0)
      node_ops->addScalar(nvindx[1], nvindx[1], dcomplex(0.0, -4.0));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_addScalar(6.0, -4.0);
      if (!complexDataSameAsValue(nvindx[1], val_addScalar, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #9a: math::HierarchyNodeDataOpsComplex::addScalar()\n"
         << "Expected: v1 = " << val_addScalar << std::endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #9b: math::HierarchyNodeDataOpsComplex::addScalar()
      // Expected: v2 = v2 + (0.0,0.25) = (1.0,1.0)
      node_ops->addScalar(nvindx[2], nvindx[2], dcomplex(0.0, 0.25));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      val_addScalar = dcomplex(1.0, 1.0);
      if (!complexDataSameAsValue(nvindx[2], val_addScalar, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #9b: math::HierarchyNodeDataOpsComplex::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << std::endl;
         node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test #9c: math::HierarchyNodeDataOpsComplex::addScalar()
      // Expected: v2 = v2 + (3.0,-4.0) = (4.0,-3.0)
      node_ops->addScalar(nvindx[2], nvindx[2], dcomplex(3.0, -4.0));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      val_addScalar = dcomplex(4.0, -3.0);
      if (!complexDataSameAsValue(nvindx[2], val_addScalar, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #9c: math::HierarchyNodeDataOpsComplex::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << std::endl;
         node_ops->printData(nvindx[2], tbox::plog);
      }

      // Reset v3: v3 = (0.5, 0.0)
      node_ops->setToScalar(nvindx[3], dcomplex(0.5, 0.0));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      // Test #10: math::HierarchyNodeDataOpsComplex::multiply()
      // Expected: v1 = v3 * v1 = (3.0,-2.0)
      node_ops->multiply(nvindx[1], nvindx[3], nvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_mult(3.0, -2.0);
      if (!complexDataSameAsValue(nvindx[1], val_mult, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #10: math::HierarchyNodeDataOpsComplex::multiply()\n"
         << "Expected: v1 = " << val_mult << std::endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #11: math::HierarchyNodeDataOpsComplex::divide()
      // Expected: v0 = v2 / v1 = (1.3846153846154,-0.076923076923077)
      node_ops->divide(nvindx[0], nvindx[2], nvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_div(1.3846153846154, -0.076923076923077);
      if (!complexDataSameAsValue(nvindx[0], val_div, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #11: math::HierarchyNodeDataOpsComplex::divide()\n"
         << "Expected: v0 = " << val_div << std::endl;
         node_ops->printData(nvindx[0], tbox::plog);
      }

      // Test #12: math::HierarchyNodeDataOpsComplex::reciprocal()
      // Expected: v1 = 1 / v1 = (0.23076923076923, 0.15384615384615)
      node_ops->reciprocal(nvindx[1], nvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_rec(0.23076923076923, 0.15384615384615);
      if (!complexDataSameAsValue(nvindx[1], val_rec, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #12: math::HierarchyNodeDataOpsComplex::reciprocal()\n"
         << "Expected: v1 = " << val_rec << std::endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }
      // Test #13: Place some bogus values on coarse level
      std::shared_ptr<pdat::NodeData<dcomplex> > ndata;

      // set values
      std::shared_ptr<hier::PatchLevel> level_zero(
         hierarchy->getPatchLevel(0));
      for (hier::PatchLevel::iterator ip(level_zero->begin());
           ip != level_zero->end(); ++ip) {
         patch = *ip;
         ndata = SAMRAI_SHARED_PTR_CAST<pdat::NodeData<dcomplex>,
                            hier::PatchData>(patch->getPatchData(nvindx[2]));
         TBOX_ASSERT(ndata);
         hier::Index index0(dim, 2);
         hier::Index index1(dim, 3);
         index1(0) = 5;
         if (dim.getValue() == 2) {
            if (patch->getBox().contains(index0)) {
               (*ndata)(pdat::NodeIndex(index0, pdat::NodeIndex::LowerLeft), 0) =
                  dcomplex(100.0, -50.0);
            }
            if (patch->getBox().contains(index1)) {
               (*ndata)(pdat::NodeIndex(index1, pdat::NodeIndex::UpperRight), 0) =
                  dcomplex(-1000.0, 20.0);
            }
         } else {
            if (patch->getBox().contains(index0)) {
               (*ndata)(pdat::NodeIndex(index0, pdat::NodeIndex::LLL), 0) =
                  dcomplex(100.0, -50.0);
            }
            if (patch->getBox().contains(index1)) {
               (*ndata)(pdat::NodeIndex(index1, pdat::NodeIndex::UUU), 0) =
                  dcomplex(-1000.0, 20.0);
            }
         }
      }

      // check values
      bool bogus_value_test_passed = true;
      for (hier::PatchLevel::iterator ipp(level_zero->begin());
           ipp != level_zero->end(); ++ipp) {
         patch = *ipp;
         ndata = SAMRAI_SHARED_PTR_CAST<pdat::NodeData<dcomplex>,
                            hier::PatchData>(patch->getPatchData(nvindx[2]));
         TBOX_ASSERT(ndata);
         hier::Index idx0(dim, 2);
         hier::Index idx1(dim, 3);
         idx1(0) = 5;
         pdat::NodeIndex::Corner corner0;
         pdat::NodeIndex::Corner corner1;
         if (dim.getValue() == 2) {
            corner0 = pdat::NodeIndex::LowerLeft;
            corner1 = pdat::NodeIndex::UpperRight;
         } else {
            corner0 = pdat::NodeIndex::LLL;
            corner1 = pdat::NodeIndex::UUU;
         }
         pdat::NodeIndex index0(idx0, corner0);
         pdat::NodeIndex index1(idx1, corner1);

         pdat::NodeIterator cend(pdat::NodeGeometry::end(ndata->getBox()));
         for (pdat::NodeIterator c(pdat::NodeGeometry::begin(ndata->getBox()));
              c != cend && bogus_value_test_passed; ++c) {
            pdat::NodeIndex node_index = *c;

            if (node_index == index0) {
               if (!tbox::MathUtilities<dcomplex>::equalEps((*ndata)(node_index),
                      dcomplex(100.0, -50.0))) {
                  bogus_value_test_passed = false;
               }
            } else {
               if (node_index == index1) {
                  if (!tbox::MathUtilities<dcomplex>::equalEps((*ndata)(node_index),
                         dcomplex(-1000.0, 20.0))) {
                     bogus_value_test_passed = false;
                  }
               } else {
                  if (!tbox::MathUtilities<dcomplex>::equalEps((*ndata)(node_index),
                         dcomplex(4.0, -3.0))) {
                     bogus_value_test_passed = false;
                  }
               }
            }
         }
      }
      if (!bogus_value_test_passed) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #13:  Place some bogus values on coarse level"
         << std::endl;
         node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test norms on patch data with nvindx[2] on hierarchy with bogus values

      // Test #14: math::HierarchyNodeDataOpsComplex::L1Norm() - w/o control weight
      // Expected:  bogus_l1_norm = 1787.0034
      double bogus_l1_norm = node_ops->L1Norm(nvindx[2]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      {
         double compare;
         if (dim.getValue() == 2) {
            compare = 1787.0034;
         } else {
            compare = 6837.00337888;
         }
         if (!tbox::MathUtilities<double>::equalEps(bogus_l1_norm, compare)) {
            ++num_failures;
            tbox::perr
            << "FAILED: - Test #14: math::HierarchyCellDataOpsComplex::L1Norm()"
            << "Expected value = " << compare << ", Computed value = "
            << std::setprecision(12) << bogus_l1_norm << std::endl;
         }
      }

      // Test #15: math::HierarchyNodeDataOpsComplex::L1Norm() - w/ control weight
      // Expected: l1_norm = 2.5
      double correct_l1_norm = node_ops->L1Norm(nvindx[2], nwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(correct_l1_norm, 2.5)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #15: math::HierarchyNodeDataOpsComplex::L1Norm()"
         << " - w/ control weight\n"
         << "Expected value = 2.5, Computed value = "
         << correct_l1_norm << std::endl;
      }

      // Test #16: math::HierarchyNodeDataOpsComplex::L2Norm()
      // Expected: l2_norm = 3.53553390593
      double l2_norm = node_ops->L2Norm(nvindx[2], nwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(l2_norm, 3.53553390593)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #16: math::HierarchyNodeDataOpsComplex::L2Norm()\n"
         << "Expected value = 3.53553390593, Computed value = "
         << l2_norm << std::endl;
      }

      // Test #17: math::HierarchyNodeDataOpsComplex::maxNorm() -w/o control weight
      // Expected: bogus_max_norm = 1000.19998
      double bogus_max_norm = node_ops->maxNorm(nvindx[2]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(bogus_max_norm, 1000.19998)) {
         ++num_failures;
         tbox::perr << "FAILED: - Test #17: Bogus maxNorm of v2\n"
                    << "Expected value = 1000.19998, Computed value = "
                    << bogus_max_norm << std::endl;
      }

      // Test #18: math::HierarchyNodeDataOpsComplex::maxNorm() -w/control weight
      // Expected: max_norm = 5.0
      double max_norm = node_ops->maxNorm(nvindx[2], nwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(max_norm, 5.0)) {
         ++num_failures;
         tbox::perr << "FAILED: - Test #18: maxNorm of v2\n"
                    << "Expected value = 5.0, Computed value = "
                    << max_norm << std::endl;
      }

      // Reset data and test sums, axpy's
      node_ops->setToScalar(nvindx[0], dcomplex(1.0, -3.0));
      node_ops->setToScalar(nvindx[1], dcomplex(2.5, 3.0));
      node_ops->setToScalar(nvindx[2], dcomplex(7.0, 0.0));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      // Test #19: math::HierarchyNodeDataOpsComplex::linearSum()
      // Expected:  v3 = (2.0,5.0)
      node_ops->linearSum(nvindx[3],
         dcomplex(2.0, 0.0), nvindx[1], dcomplex(0.0, -1.0), nvindx[0]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_linearSum(2.0, 5.0);
      if (!complexDataSameAsValue(nvindx[3], val_linearSum, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #19: math::HierarchyNodeDataOpsComplex::linearSum()\n"
         << "Expected: v3 = " << val_linearSum << std::endl;
         node_ops->printData(nvindx[3], tbox::plog);
      }

      // Test #20: math::HierarchyNodeDataOpsComplex::axmy()
      // Expected:  v3 = (6.5,12.0)
      node_ops->axmy(nvindx[3], 3.0, nvindx[1], nvindx[0]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_axmy(6.5, 12.0);
      if (!complexDataSameAsValue(nvindx[3], val_axmy, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #20: math::HierarchyNodeDataOpsComplex::axmy()\n"
         << "Expected: v3 = " << val_axmy << std::endl;
         node_ops->printData(nvindx[3], tbox::plog);
      }

      // Test #21a: math::HierarchyNodeDataOpsComplex::dot()
      // Expected:  cdot = (8.75,-10.5)
      dcomplex cdot = node_ops->dot(nvindx[2], nvindx[1], nwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex ans_2_dot_1(8.75, -10.5);
      if (!tbox::MathUtilities<dcomplex>::equalEps(cdot, ans_2_dot_1)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #21a: math::HierarchyNodeDataOpsComplex::dot()\n"
         << "Expected value = (8.75,-10.5), Computed value = "
         << cdot << std::endl;
      }

      // Test #21b: math::HierarchyNodeDataOpsComplex::dot()
      // Expected:  cdot = (8.75,10.5)
      dcomplex cdot2 = node_ops->dot(nvindx[1], nvindx[2], nwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex ans_1_dot_2(8.75, 10.5);
      if (!tbox::MathUtilities<dcomplex>::equalEps(cdot2, ans_1_dot_2)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #21b: math::HierarchyNodeDataOpsComplex::dot()\n"
         << "Expected value = (8.75,10.5), Computed value = "
         << cdot2 << std::endl;
      }

      // Test #22: math::HierarchyNodeDataOpsComplex::abs()
      // Expected: abs(v0) = 5.0
      node_ops->setToScalar(nvindx[0], dcomplex(4.0, -3.0));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      node_ops->abs(nwgt_id, nvindx[0]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!doubleDataSameAsValue(nwgt_id, 5.0, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #22: math::HierarchyNodeDataOpsComplex::abs()\n"
         << "Expected: abs(v0) = 5.0" << std::endl;
         nwgt_ops->printData(nwgt_id, tbox::plog);
      }

      // deallocate data on hierarchy
      for (ln = 0; ln < 2; ++ln) {
         hierarchy->getPatchLevel(ln)->deallocatePatchData(nwgt_id);
         for (iv = 0; iv < NVARS; ++iv) {
            hierarchy->getPatchLevel(ln)->deallocatePatchData(nvindx[iv]);
         }
      }

      for (iv = 0; iv < NVARS; ++iv) {
         nvar[iv].reset();
      }
      nwgt.reset();

      geometry.reset();
      hierarchy.reset();
      node_ops.reset();
      nwgt_ops.reset();

      if (num_failures == 0) {
         tbox::pout << "\nPASSED:  node cplxtest" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return num_failures;
}

/*
 * Returns true if all the data in the hierarchy is equal to the specified
 * value.  Returns false otherwise.
 */
static bool
complexDataSameAsValue(
   int desc_id,
   dcomplex value,
   std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   bool test_passed = true;

   int ln;
   std::shared_ptr<hier::Patch> patch;
   for (ln = 0; ln < 2; ++ln) {

      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         patch = *ip;
         std::shared_ptr<pdat::NodeData<dcomplex> > nvdata(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<dcomplex>, hier::PatchData>(
               patch->getPatchData(desc_id)));

         TBOX_ASSERT(nvdata);

         pdat::NodeIterator cend(pdat::NodeGeometry::end(nvdata->getBox()));
         for (pdat::NodeIterator c(pdat::NodeGeometry::begin(nvdata->getBox()));
              c != cend && test_passed; ++c) {
            pdat::NodeIndex node_index = *c;
            if (!tbox::MathUtilities<dcomplex>::equalEps((*nvdata)(node_index),
                   value)) {
               test_passed = false;
            }
         }
      }
   }

   return test_passed;
}

/*
 * Returns true if all the data in the hierarchy is equal to the specified
 * value.  Returns false otherwise.
 */
static bool
doubleDataSameAsValue(
   int desc_id,
   double value,
   std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   bool test_passed = true;

   int ln;
   std::shared_ptr<hier::Patch> patch;
   for (ln = 0; ln < 2; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         patch = *ip;
         std::shared_ptr<pdat::NodeData<double> > nvdata(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
               patch->getPatchData(desc_id)));

         TBOX_ASSERT(nvdata);

         pdat::NodeIterator cend(pdat::NodeGeometry::end(nvdata->getBox()));
         for (pdat::NodeIterator c(pdat::NodeGeometry::begin(nvdata->getBox()));
              c != cend && test_passed; ++c) {
            pdat::NodeIndex node_index = *c;
            if (!tbox::MathUtilities<double>::equalEps((*nvdata)(node_index),
                   value)) {
               test_passed = false;
            }
         }
      }
   }

   return test_passed;
}
