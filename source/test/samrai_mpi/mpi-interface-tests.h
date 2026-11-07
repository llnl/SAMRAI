/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test multiple SAMRAI interfaces.
 *
 ************************************************************************/

/*!
 * @brief Test multiple SAMRAI_MPI interfaces for a given runtime_mpi mode.
 *
 * @param fail_count Increment this count by number of failures.
 * @param runtime_mpi Whether to use run-time MPI.
 */
void
mpiInterfaceTests(
   int& argc,
   char **& argv,
   int& fail_count,
   bool runtime_mpi,
   bool mpi_disabled);

/*!
 * @brief Test SAMRAI_MPI::Bcast.
 *
 * @param fail_count Increment this count by number of failures.
 *
 * @return number of failures found.
 */
int
mpiInterfaceTestBcast(
   int& fail_count);

/*!
 * @brief Test SAMRAI_MPI::Allreduce.
 *
 * @param fail_count Increment this count by number of failures.
 *
 * @return number of failures found.
 */
int
mpiInterfaceTestAllreduce(
   int& fail_count);

/*!
 * @brief Test SAMRAI_MPI::parallelPrefixSum.
 *
 * @param fail_count Increment this count by number of failures.
 *
 * @return number of failures found.
 */
int
mpiInterfaceTestParallelPrefixSum(
   int& fail_count);
