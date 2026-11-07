/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Simple utility class for interfacing with MPI
 *
 ************************************************************************/

#include "SAMRAI/tbox/SAMRAI_MPI.h"

#ifdef SAMRAI_HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif

#ifdef SAMRAI_HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <stdlib.h>

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/Utilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

#ifdef __INSURE__
/*
 * These are defined in mpich mpi.h and break the insure compile.
 * This may impact Globus in some way, at least from the comments
 * in the mpi.h header file.  Why mpich externs something that is
 * not defined in the mpich is confusing and probably just broken.
 */
int MPICHX_TOPOLOGY_DEPTHS;
int MPICHX_TOPOLOGY_COLORS;
int MPICHX_PARALLELSOCKETS_PARAMETERS;
#endif

namespace SAMRAI {
namespace tbox {

const SAMRAI_MPI::Comm SAMRAI_MPI::commWorld = MPI_COMM_WORLD;
const SAMRAI_MPI::Comm SAMRAI_MPI::commNull = MPI_COMM_NULL;

bool SAMRAI_MPI::s_mpi_is_initialized = false;
bool SAMRAI_MPI::s_we_started_mpi(false);
SAMRAI_MPI SAMRAI_MPI::s_samrai_world(MPI_COMM_NULL);

bool SAMRAI_MPI::s_call_abort_in_serial_instead_of_exit = false;
bool SAMRAI_MPI::s_call_abort_in_parallel_instead_of_mpiabort = false;

/*
 **************************************************************************
 * Constructor.
 **************************************************************************
 */
SAMRAI_MPI::SAMRAI_MPI(
   const Comm& comm):
   d_comm(comm),
   d_rank(-1),
   d_size(-1)
{
   if (comm != MPI_COMM_NULL) {
#ifdef HAVE_MPI
      if (s_mpi_is_initialized) {
         MPI_Comm_rank(d_comm, &d_rank);
         MPI_Comm_size(d_comm, &d_size);
      }
#else
      d_rank = 0;
      d_size = 1;
#endif
   }
}

/*
 **************************************************************************
 * Copy constructor.
 **************************************************************************
 */
SAMRAI_MPI::SAMRAI_MPI(
   const SAMRAI_MPI& other):
   d_comm(other.d_comm),
   d_rank(other.d_rank),
   d_size(other.d_size)
{
}

/*
 **************************************************************************
 *
 * Abort the program.
 *
 **************************************************************************
 */

void
SAMRAI_MPI::abort()
{

#ifdef HAVE_MPI
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   if (mpi.getSize() > 1) {
      if (s_call_abort_in_parallel_instead_of_mpiabort) {
         ::abort();
      } else {
         MPI_Abort(mpi.getCommunicator(), -1);
      }
   } else {
      if (s_call_abort_in_serial_instead_of_exit) {
         ::abort();
      } else {
         exit(-1);
      }
   }
#else
   if (s_call_abort_in_serial_instead_of_exit) {
      ::abort();
   } else {
      exit(-1);
   }
#endif

}

/*
 **************************************************************************
 *
 * Wrapper for MPI_Init().
 *
 **************************************************************************
 */
void
SAMRAI_MPI::init(
   int* argc,
   char** argv[])
{
#ifdef HAVE_MPI

   MPI_Init(argc, argv);
   s_mpi_is_initialized = true;
   s_we_started_mpi = true;

   Comm dup_comm;
   MPI_Comm_dup(MPI_COMM_WORLD, &dup_comm);
   s_samrai_world.setCommunicator(dup_comm);

#else
   NULL_USE(argc);
   NULL_USE(argv);
   s_samrai_world.d_comm = MPI_COMM_WORLD;
   s_samrai_world.d_size = 1;
   s_samrai_world.d_rank = 0;
#endif

   if (getenv("SAMRAI_ABORT_ON_ERROR")) {
      SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(true);
      SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort(true);
   }
}

/*
 **************************************************************************
 *
 * Wrapper for MPI_Init().
 *
 **************************************************************************
 */
void
SAMRAI_MPI::init(
   Comm comm)
{
   if (comm == MPI_COMM_NULL) {
      std::cerr << "SAMRAI_MPI::init: invalid initializing Communicator."
                << std::endl;
   }
#ifdef HAVE_MPI

   s_mpi_is_initialized = true;
   s_we_started_mpi = false;

   Comm dup_comm;
   MPI_Comm_dup(comm, &dup_comm);
   s_samrai_world.setCommunicator(dup_comm);

#endif

   if (getenv("SAMRAI_ABORT_ON_ERROR")) {
      SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(true);
      SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort(true);
   }
}

/*
 **************************************************************************
 * Initialize SAMRAI_MPI with MPI disabled.
 **************************************************************************
 */
void
SAMRAI_MPI::initMPIDisabled()
{
   s_mpi_is_initialized = false;
   s_we_started_mpi = false;

   s_samrai_world.d_comm = MPI_COMM_WORLD;
   s_samrai_world.d_size = 1;
   s_samrai_world.d_rank = 0;

   if (getenv("SAMRAI_ABORT_ON_ERROR")) {
      SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(true);
      SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort(true);
   }
}

/*
 **************************************************************************
 *
 * Wrapper for MPI_Finalize().
 *
 **************************************************************************
 */
void
SAMRAI_MPI::finalize()
{
#ifdef HAVE_MPI
   if (s_mpi_is_initialized) {
      MPI_Comm_free(&s_samrai_world.d_comm);
   } else {
      s_samrai_world.d_comm = MPI_COMM_NULL;
   }

   if (s_we_started_mpi) {
      MPI_Finalize();
   }
#endif
}

/*
 *****************************************************************************
 *
 * Methods named like MPI's native interfaces (without the MPI_ prefix)
 * are wrappers for the native interfaces.  The SAMRAI_MPI versions
 * introduce a flag to determine whether MPI is really used at run time.
 * When the run-time flag is on, these wrappers are identical to the MPI
 * versions.  When the flag is off, most of these methods are no-ops
 * (which is not necessarily the same as calling the MPI functions with
 * only 1 process in the communicator).
 *
 *****************************************************************************
 */

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Comm_rank(
   Comm comm,
   int* rank)
{
#ifndef HAVE_MPI
   NULL_USE(comm);
#endif
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Comm_rank is a no-op without run-time MPI!");
   }
   int rval = MPI_SUCCESS;
#ifdef HAVE_MPI
   rval = MPI_Comm_rank(comm, rank);
#else
   *rank = 0;
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Comm_size(
   Comm comm,
   int* size)
{
#ifndef HAVE_MPI
   NULL_USE(comm);
#endif
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Comm_size is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   return MPI_Comm_size(comm, size);

#else
   *size = 1;
   return MPI_SUCCESS;

#endif
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Comm_compare(
   Comm comm1,
   Comm comm2,
   int* result)
{
#ifndef HAVE_MPI
   NULL_USE(comm1);
   NULL_USE(comm2);
   NULL_USE(result);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Comm_compare is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Comm_compare(comm1, comm2, result);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Comm_free(
   Comm* comm)
{
#ifndef HAVE_MPI
   NULL_USE(comm);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Comm_free is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Comm_free(comm);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Finalized(
   int* flag)
{
   int rval = MPI_SUCCESS;
#ifdef HAVE_MPI
   rval = MPI_Finalized(flag);
#else
   *flag = true;
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Get_count(
   Status* status,
   Datatype datatype,
   int* count)
{
#ifndef HAVE_MPI
   NULL_USE(status);
   NULL_USE(datatype);
   NULL_USE(count);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Get_count is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Get_count(status, datatype, count);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Request_free(
   Request* request)
{
#ifndef HAVE_MPI
   NULL_USE(request);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Get_count is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Request_free(request);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Test(
   Request* request,
   int* flag,
   Status* status)
{
#ifndef HAVE_MPI
   NULL_USE(request);
   NULL_USE(flag);
   NULL_USE(status);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Test is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Test(request, flag, status);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Test_cancelled(
   Status* status,
   int* flag)
{
#ifndef HAVE_MPI
   NULL_USE(status);
   NULL_USE(flag);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Test_canceled is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Test_cancelled(status, flag);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Wait(
   Request* request,
   Status* status)
{
#ifndef HAVE_MPI
   NULL_USE(request);
   NULL_USE(status);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Wait is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Wait(request, status);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Waitall(
   int count,
   Request* reqs,
   Status* stats)
{
#ifndef HAVE_MPI
   NULL_USE(count);
   NULL_USE(reqs);
   NULL_USE(stats);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Waitall is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Waitall(count, reqs, stats);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Waitany(
   int count,
   Request* array_of_requests,
   int* index,
   Status* status)
{
#ifndef HAVE_MPI
   NULL_USE(count);
   NULL_USE(array_of_requests);
   NULL_USE(index);
   NULL_USE(status);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Waitany is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Waitany(count, array_of_requests, index, status);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Waitsome(
   int incount,
   Request* array_of_requests,
   int* outcount,
   int* array_of_indices,
   Status* array_of_statuses)
{
#ifndef HAVE_MPI
   NULL_USE(incount);
   NULL_USE(array_of_requests);
   NULL_USE(outcount);
   NULL_USE(array_of_indices);
   NULL_USE(array_of_statuses);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Waitsome is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Waitsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 * If MPI is enabled, use MPI_Wtime.
 * Else if POSIX time is available, use POSIX time.
 * Else, return 0.
 *****************************************************************************
 */
double
SAMRAI_MPI::Wtime()
{
   double rval = 0.0;
   if (!s_mpi_is_initialized) {
#ifdef SAMRAI_HAVE_SYS_TIMES_H
      // Without MPI, use POSIX time.
      struct tms tmp_tms;
      clock_t clock_ticks_since_reference = times(&tmp_tms);
      const double clock_ticks_per_second = double(sysconf(_SC_CLK_TCK));
      rval = double(clock_ticks_since_reference) / clock_ticks_per_second;
#endif
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Wtime();
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *
 * MPI wrappers with the communicator removed from the argument list.
 *
 *****************************************************************************
 */

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Allgather(
   void* sendbuf,
   int sendcount,
   Datatype sendtype,
   void* recvbuf,
   int recvcount,
   Datatype recvtype) const
{
#ifndef HAVE_MPI
   NULL_USE(sendbuf);
   NULL_USE(sendcount);
   NULL_USE(sendtype);
   NULL_USE(recvbuf);
   NULL_USE(recvcount);
   NULL_USE(recvtype);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Allgather is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, d_comm);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Allgatherv(
   void* sendbuf,
   int sendcounts,
   Datatype sendtype,
   void* recvbuf,
   int* recvcounts,
   int* displs,
   Datatype recvtype) const
{
#ifndef HAVE_MPI
   NULL_USE(sendbuf);
   NULL_USE(sendcounts);
   NULL_USE(sendtype);
   NULL_USE(recvbuf);
   NULL_USE(recvcounts);
   NULL_USE(displs);
   NULL_USE(recvtype);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Algatherv is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Allgatherv(sendbuf,
            sendcounts,
            sendtype,
            recvbuf,
            recvcounts,
            displs,
            recvtype,
            d_comm);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Allreduce(
   void* sendbuf,
   void* recvbuf,
   int count,
   Datatype datatype,
   Op op) const
{
#ifndef HAVE_MPI
   NULL_USE(sendbuf);
   NULL_USE(recvbuf);
   NULL_USE(count);
   NULL_USE(datatype);
   NULL_USE(op);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Allreduce is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, d_comm);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Attr_get(
   int keyval,
   void* attribute_val,
   int* flag) const
{
#ifndef HAVE_MPI
   NULL_USE(keyval);
   NULL_USE(attribute_val);
   NULL_USE(flag);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Attr_get is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Attr_get(d_comm, keyval, attribute_val, flag);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Barrier() const
{
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      // A no-op is OK for sequential Barrier.
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Barrier(d_comm);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Bcast(
   void* buffer,
   int count,
   Datatype datatype,
   int root) const
{
#ifndef HAVE_MPI
   NULL_USE(buffer);
   NULL_USE(count);
   NULL_USE(datatype);
   NULL_USE(root);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      // A no-op is OK for sequential Bcast.
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Bcast(buffer, count, datatype, root, d_comm);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Comm_dup(
   Comm* newcomm) const
{
#ifndef HAVE_MPI
   NULL_USE(newcomm);
#endif
   *newcomm = commNull;
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Comm_dup is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Comm_dup(d_comm, newcomm);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Comm_rank(
   int* rank) const
{
#ifndef HAVE_MPI
   NULL_USE(rank);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      *rank = 0;
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Comm_rank(d_comm, rank);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Comm_size(
   int* size) const
{
#ifndef HAVE_MPI
   NULL_USE(size);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      *size = 1;
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Comm_size(d_comm, size);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Gather(
   void* sendbuf,
   int sendcount,
   Datatype sendtype,
   void* recvbuf,
   int recvcount,
   Datatype recvtype,
   int root) const
{
#ifndef HAVE_MPI
   NULL_USE(sendbuf);
   NULL_USE(sendcount);
   NULL_USE(sendtype);
   NULL_USE(recvbuf);
   NULL_USE(recvcount);
   NULL_USE(recvtype);
   NULL_USE(root);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Gather is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, d_comm);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Gatherv(
   void* sendbuf,
   int sendcount,
   Datatype sendtype,
   void* recvbuf,
   int* recvcounts,
   int* displs,
   Datatype recvtype,
   int root) const
{
#ifndef HAVE_MPI
   NULL_USE(sendbuf);
   NULL_USE(sendcount);
   NULL_USE(sendtype);
   NULL_USE(recvbuf);
   NULL_USE(recvcounts);
   NULL_USE(displs);
   NULL_USE(recvtype);
   NULL_USE(root);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Gatherv is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Gatherv(sendbuf,
            sendcount,
            sendtype,
            recvbuf,
            recvcounts,
            displs,
            recvtype,
            root,
            d_comm);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Iprobe(
   int source,
   int tag,
   int* flag,
   Status* status) const
{
#ifndef HAVE_MPI
   NULL_USE(source);
   NULL_USE(tag);
   NULL_USE(flag);
   NULL_USE(status);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Iprobe is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Iprobe(source, tag, d_comm, flag, status);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Isend(
   void* buf,
   int count,
   Datatype datatype,
   int dest,
   int tag,
   Request* req) const
{
#ifndef HAVE_MPI
   NULL_USE(buf);
   NULL_USE(count);
   NULL_USE(datatype);
   NULL_USE(dest);
   NULL_USE(tag);
   NULL_USE(req);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Isend is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Isend(buf, count, datatype, dest, tag, d_comm, req);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Irecv(
   void* buf,
   int count,
   Datatype datatype,
   int source,
   int tag,
   Request* request) const
{
#ifndef HAVE_MPI
   NULL_USE(buf);
   NULL_USE(count);
   NULL_USE(datatype);
   NULL_USE(source);
   NULL_USE(tag);
   NULL_USE(request);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Irecv is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Irecv(buf, count, datatype, source, tag, d_comm, request);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Probe(
   int source,
   int tag,
   Status* status) const
{
#ifndef HAVE_MPI
   NULL_USE(source);
   NULL_USE(tag);
   NULL_USE(status);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Probe is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Probe(source, tag, d_comm, status);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Recv(
   void* buf,
   int count,
   Datatype datatype,
   int source,
   int tag,
   Status* status) const
{
#ifndef HAVE_MPI
   NULL_USE(buf);
   NULL_USE(count);
   NULL_USE(datatype);
   NULL_USE(source);
   NULL_USE(tag);
   NULL_USE(status);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Recv is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Recv(buf, count, datatype, source, tag, d_comm, status);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Reduce(
   void* sendbuf,
   void* recvbuf,
   int count,
   Datatype datatype,
   Op op,
   int root) const
{
#ifndef HAVE_MPI
   NULL_USE(sendbuf);
   NULL_USE(recvbuf);
   NULL_USE(count);
   NULL_USE(datatype);
   NULL_USE(op);
   NULL_USE(root);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Reduce is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, d_comm);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Send(
   void* buf,
   int count,
   Datatype datatype,
   int dest,
   int tag) const
{
#ifndef HAVE_MPI
   NULL_USE(buf);
   NULL_USE(count);
   NULL_USE(datatype);
   NULL_USE(dest);
   NULL_USE(tag);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Send is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Send(buf, count, datatype, dest, tag, d_comm);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Sendrecv(
   void* sendbuf, int sendcount, Datatype sendtype, int dest, int sendtag,
   void* recvbuf, int recvcount, Datatype recvtype, int source, int recvtag,
   Status* status) const
{
#ifndef HAVE_MPI
   NULL_USE(sendbuf);
   NULL_USE(sendcount);
   NULL_USE(sendtype);
   NULL_USE(dest);
   NULL_USE(sendtag);
   NULL_USE(recvbuf);
   NULL_USE(recvcount);
   NULL_USE(recvtype);
   NULL_USE(source);
   NULL_USE(recvtag);
   NULL_USE(status);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Send is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Sendrecv(
            sendbuf, sendcount, sendtype, dest, sendtag,
            recvbuf, recvcount, recvtype, source, recvtag,
            d_comm, status);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *****************************************************************************
 */
int
SAMRAI_MPI::Scan(
   void* sendbuf,
   void* recvbuf,
   int count,
   Datatype datatype,
   Op op) const
{
#ifndef HAVE_MPI
   NULL_USE(sendbuf);
   NULL_USE(recvbuf);
   NULL_USE(count);
   NULL_USE(datatype);
   NULL_USE(op);
#endif
   int rval = MPI_SUCCESS;
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::Scan is a no-op without run-time MPI!");
   }
#ifdef HAVE_MPI
   else {
      rval = MPI_Scan(sendbuf, recvbuf, count, datatype, op, d_comm);
   }
#endif
   return rval;
}

/*
 *****************************************************************************
 *
 * Methods named like MPI's native interfaces (without the MPI_ prefix)
 * are wrappers for the native interfaces.  The SAMRAI_MPI versions
 * introduce a flag to determine whether MPI is really used at run time.
 * When the run-time flag is on, these wrappers are identical to the MPI
 * versions.  When the flag is off, most of these methods are no-ops
 * (which is not necessarily the same as calling the MPI functions with
 * only 1 process in the communicator).
 *
 *****************************************************************************
 */

/*
 **************************************************************************
 * Specialized Allreduce for integers.
 **************************************************************************
 */
int
SAMRAI_MPI::AllReduce(
   int* x,
   int count,
   Op op,
   int* ranks_of_extrema) const
{
#ifndef HAVE_MPI
   NULL_USE(x);
   NULL_USE(count);
   NULL_USE(op);
   NULL_USE(ranks_of_extrema);
#endif
   if ((op == MPI_MINLOC || op == MPI_MAXLOC) &&
       ranks_of_extrema == 0) {
      TBOX_ERROR("SAMRAI_MPI::AllReduce: If you specify reduce\n"
         << "operation MPI_MINLOC or MPI_MAXLOC, you must\n"
         << "provide space for the ranks in the 'ranks_of_extrema'\n"
         << "argument.");
   }
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::AllReduce is a no-op without run-time MPI!");
   }

   int rval = MPI_SUCCESS;
   /*
    * Get ranks of extrema if user operation specified it or user
    * specified min/max operation and provides space for rank.
    */
   bool get_ranks_of_extrema =
      op == MPI_MINLOC ? true :
      op == MPI_MAXLOC ? true :
      ranks_of_extrema != 0 && (op == MPI_MIN || op == MPI_MAX);

   if (!get_ranks_of_extrema) {
      std::vector<int> recv_buf(count);
      rval = Allreduce(x, &recv_buf[0], count, MPI_INT, op);
      for (int c = 0; c < count; ++c) {
         x[c] = recv_buf[c];
      }
   } else {
      Op locop =
         op == MPI_MIN ? MPI_MINLOC :
         op == MPI_MAX ? MPI_MAXLOC :
         op;
      IntIntStruct* send_buf = new IntIntStruct[count];
      IntIntStruct* recv_buf = new IntIntStruct[count];
      for (int c = 0; c < count; ++c) {
         send_buf[c].j = x[c];
         send_buf[c].i = d_rank;
      }
      rval = Allreduce(send_buf, recv_buf, count, MPI_2INT, locop);
      for (int c = 0; c < count; ++c) {
         x[c] = recv_buf[c].j;
         ranks_of_extrema[c] = recv_buf[c].i;
      }

      delete[] send_buf;
      delete[] recv_buf;
   }

   return rval;
}

/*
 **************************************************************************
 * Specialized Allreduce for doubles.
 **************************************************************************
 */
int
SAMRAI_MPI::AllReduce(
   double* x,
   int count,
   Op op,
   int* ranks_of_extrema) const
{
#ifndef HAVE_MPI
   NULL_USE(x);
   NULL_USE(count);
   NULL_USE(op);
   NULL_USE(ranks_of_extrema);
#endif
   if ((op == MPI_MINLOC || op == MPI_MAXLOC) &&
       ranks_of_extrema == 0) {
      TBOX_ERROR("SAMRAI_MPI::AllReduce: If you specify reduce\n"
         << "operation MPI_MINLOC or MPI_MAXLOC, you must\n"
         << "provide space for the ranks in the 'ranks_of_extrema'\n"
         << "argument.");
   }
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::AllReduce is a no-op without run-time MPI!");
   }

   int rval = MPI_SUCCESS;
   /*
    * Get ranks of extrema if user operation specified it or user
    * specified min/max operation and provides space for rank.
    */
   bool get_ranks_of_extrema =
      op == MPI_MINLOC ? true :
      op == MPI_MAXLOC ? true :
      ranks_of_extrema != 0 && (op == MPI_MIN || op == MPI_MAX);

   if (!get_ranks_of_extrema) {
      std::vector<double> recv_buf(count);
      rval = Allreduce(x, &recv_buf[0], count, MPI_DOUBLE, op);
      for (int c = 0; c < count; ++c) {
         x[c] = recv_buf[c];
      }
   } else {
      Op locop =
         op == MPI_MIN ? MPI_MINLOC :
         op == MPI_MAX ? MPI_MAXLOC :
         op;
      DoubleIntStruct* send_buf = new DoubleIntStruct[count];
      DoubleIntStruct* recv_buf = new DoubleIntStruct[count];
      for (int c = 0; c < count; ++c) {
         send_buf[c].d = x[c];
         send_buf[c].i = d_rank;
      }
      rval = Allreduce(send_buf, recv_buf, count, MPI_DOUBLE_INT, locop);
      for (int c = 0; c < count; ++c) {
         x[c] = recv_buf[c].d;
         ranks_of_extrema[c] = recv_buf[c].i;
      }

      delete[] send_buf;
      delete[] recv_buf;
   }

   return rval;
}

/*
 **************************************************************************
 * Specialized Allreduce for doubles.
 **************************************************************************
 */
int
SAMRAI_MPI::AllReduce(
   float* x,
   int count,
   Op op,
   int* ranks_of_extrema) const
{
#ifndef HAVE_MPI
   NULL_USE(x);
   NULL_USE(count);
   NULL_USE(op);
   NULL_USE(ranks_of_extrema);
#endif
   if ((op == MPI_MINLOC || op == MPI_MAXLOC) &&
       ranks_of_extrema == 0) {
      TBOX_ERROR("SAMRAI_MPI::AllReduce: If you specify reduce\n"
         << "operation MPI_MINLOC or MPI_MAXLOC, you must\n"
         << "provide space for the ranks in the 'ranks_of_extrema'\n"
         << "argument.");
   }
   if (!s_mpi_is_initialized) {
      TBOX_ERROR("SAMRAI_MPI::AllReduce is a no-op without run-time MPI!");
   }

   int rval = MPI_SUCCESS;
   /*
    * Get ranks of extrema if user operation specified it or user
    * specified min/max operation and provides space for rank.
    */
   bool get_ranks_of_extrema =
      op == MPI_MINLOC ? true :
      op == MPI_MAXLOC ? true :
      ranks_of_extrema != 0 && (op == MPI_MIN || op == MPI_MAX);

   if (!get_ranks_of_extrema) {
      std::vector<float> recv_buf(count);
      rval = Allreduce(x, &recv_buf[0], count, MPI_FLOAT, op);
      for (int c = 0; c < count; ++c) {
         x[c] = recv_buf[c];
      }
   } else {
      Op locop =
         op == MPI_MIN ? MPI_MINLOC :
         op == MPI_MAX ? MPI_MAXLOC :
         op;
      FloatIntStruct* send_buf = new FloatIntStruct[count];
      FloatIntStruct* recv_buf = new FloatIntStruct[count];
      for (int c = 0; c < count; ++c) {
         send_buf[c].f = x[c];
         send_buf[c].i = d_rank;
      }
      rval = Allreduce(send_buf, recv_buf, count, MPI_FLOAT_INT, locop);
      for (int c = 0; c < count; ++c) {
         x[c] = recv_buf[c].f;
         ranks_of_extrema[c] = recv_buf[c].i;
      }

      delete[] send_buf;
      delete[] recv_buf;
   }

   return rval;
}

/*
 **************************************************************************
 * Parallel prefix sum for ints.
 *
 * This method implements the parallel prefix sum algorithm.
 * The distance loop is expected to execute (ln d_size) times,
 * doing up to 1 send and 1 receive each time.
 *
 * Note: I'm not sure we have to use all non-blocking calls to get
 * good performance, but it probably can't hurt.  --BTNG
 **************************************************************************
 */
int
SAMRAI_MPI::parallelPrefixSum(
   int* x,
   int count,
   int tag) const
{
#ifndef HAVE_MPI
   NULL_USE(x);
   NULL_USE(count);
   NULL_USE(tag);
#endif

   // Scratch data.
   std::vector<int> send_scr(count), recv_scr(count);

   Request send_req, recv_req;
   Status send_stat, recv_stat;
   int mpi_err = MPI_SUCCESS;

   for (int distance = 1; distance < d_size; distance *= 2) {

      const int recv_from = d_rank - distance;
      const int send_to = d_rank + distance;

      if (recv_from >= 0) {
         mpi_err = Irecv(&recv_scr[0], count, MPI_INT, recv_from, tag, &recv_req);
         if (mpi_err != MPI_SUCCESS) {
            return mpi_err;
         }
      }

      if (send_to < d_size) {
         send_scr.clear();
         send_scr.insert(send_scr.end(), x, x + count);
         mpi_err = Isend(&send_scr[0], count, MPI_INT, send_to, tag, &send_req);
         if (mpi_err != MPI_SUCCESS) {
            return mpi_err;
         }
      }

      if (recv_from >= 0) {
         mpi_err = Wait(&recv_req, &recv_stat);
         if (mpi_err != MPI_SUCCESS) {
            return mpi_err;
         }
         for (int i = 0; i < count; ++i) {
            x[i] += recv_scr[i];
         }
      }

      if (send_to < d_size) {
         mpi_err = Wait(&send_req, &send_stat);
         if (mpi_err != MPI_SUCCESS) {
            return mpi_err;
         }
      }

   }

   return MPI_SUCCESS;
}

/*
 **************************************************************************
 * Check whether there is a receivable message, for use in guarding
 * against errant messages (message from an unrelated communication)
 * that may be mistakenly received.  This check is imperfect; it can
 * detect messages that have arrived but it can't detect messages that
 * has not arrived.
 *
 * The barriers prevent processes from starting or finishing the check
 * too early.  Early start may miss recently sent errant messages from
 * slower processes.  Early finishes can allow the process to get ahead
 * and send a valid message that may be mistaken as an errant message
 * by the receiver doing the Iprobe.
 **************************************************************************
 */
bool
SAMRAI_MPI::hasReceivableMessage(
   Status* status,
   int source,
   int tag) const
{
   int flag = false;
   if (s_mpi_is_initialized) {
      SAMRAI_MPI::Status tmp_status;
      Barrier();
      int mpi_err = Iprobe(source, tag, &flag, status ? status : &tmp_status);
      if (mpi_err != MPI_SUCCESS) {
         TBOX_ERROR("SAMRAI_MPI::hasReceivableMessage: Error probing for message." << std::endl);
      }
      Barrier();
   }
   return flag == true;
}

/*
 **************************************************************************
 **************************************************************************
 */
void
SAMRAI_MPI::dupCommunicator(
   const SAMRAI_MPI& r)
{
#ifdef HAVE_MPI
   int rval = r.Comm_dup(&d_comm);
   if (rval != MPI_SUCCESS) {
      TBOX_ERROR("SAMRAI_MPI::dupCommunicator: Error duplicating\n"
         << "communicator.");
   }
   MPI_Comm_rank(d_comm, &d_rank);
   MPI_Comm_size(d_comm, &d_size);
   TBOX_ASSERT(d_rank == r.d_rank);
   TBOX_ASSERT(d_size == r.d_size);
#else
   d_comm = r.d_comm;
   d_rank = r.d_rank;
   d_size = r.d_size;
#endif
}

/*
 **************************************************************************
 **************************************************************************
 */
void
SAMRAI_MPI::freeCommunicator()
{
#ifdef HAVE_MPI
   if (d_comm != MPI_COMM_NULL) {
      TBOX_ASSERT(SAMRAI_MPI::usingMPI());
      Comm_free(&d_comm);
      // d_comm is now set to MPI_COMM_NULL;
   }
#else
   d_comm = MPI_COMM_NULL;
#endif
   d_rank = 0;
   d_size = 1;
}

/*
 **************************************************************************
 **************************************************************************
 */
int
SAMRAI_MPI::compareCommunicator(
   const SAMRAI_MPI& r) const
{
#ifdef HAVE_MPI
   int compare_result;
   int mpi_err = Comm_compare(
         d_comm,
         r.d_comm,
         &compare_result);
   if (mpi_err != MPI_SUCCESS) {
      TBOX_ERROR("SAMRAI_MPI::compareCommunicator: Error comparing two communicators.");
   }
   return compare_result;

#else
   NULL_USE(r);
   return d_comm == r.d_comm ? MPI_IDENT : MPI_CONGRUENT;

#endif
}

/*
 **************************************************************************
 **************************************************************************
 */
void
SAMRAI_MPI::setCommunicator(
   const Comm& comm)
{
   d_comm = comm;
   d_rank = 0;
   d_size = 1;
#ifdef HAVE_MPI
   if (s_mpi_is_initialized) {
      if (d_comm != MPI_COMM_NULL) {
         MPI_Comm_rank(d_comm, &d_rank);
         MPI_Comm_size(d_comm, &d_size);
      }
   }
#endif
}

#ifndef HAVE_MPI
/*
 **************************************************************************
 **************************************************************************
 */
SAMRAI_MPI::Status::Status():
   MPI_SOURCE(-1),
   MPI_TAG(-1),
   MPI_ERROR(-1)
{
}
#endif

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Unsuppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
