/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility functions for error reporting, file manipulation, etc.
 *
 ************************************************************************/

#ifndef included_tbox_Utilities
#define included_tbox_Utilities

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/Logger.h"

#include <string>

#include <sys/types.h>
#include <sys/stat.h>

#if defined(SAMRAI_HAVE_CALIPER)
#include <caliper/cali.h>
#endif

namespace SAMRAI {
namespace tbox {

#ifdef _MSC_VER
#include <direct.h>
typedef int mode_t;
#define  S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
#define S_IRUSR 0
#define S_IWUSR 0
#define S_IXUSR 0
#endif

/*!
 * A statement that does nothing, for insure++ make it something
 * more complex than a simple C null statement to avoid a warning.
 */
#ifdef __INSURE__
#define NULL_STATEMENT if (0) int nullstatement = 0
#else
#define NULL_STATEMENT
#endif

/*!
 * A null use of a parameter to a function, use to avoid compiler warnings
 * about unused parameters in a function.
 */
#define NULL_USE_PARAM(variable)

/*!
 * A null use of a variable, use to avoid GNU compiler
 * warnings about unused variables.
 */
#define NULL_USE(variable)                               \
   static_cast< void >( variable )

/*!
 * Throw an error assertion from within any C++ source code.  The
 * macro argument may be any standard ostream expression.  The file and
 * line number of the abort are also printed.
 */
#define TBOX_ERROR(X)                                           \
   do {                                                         \
      std::ostringstream tboxos;                                \
      tboxos << X;                                              \
      SAMRAI::tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__); \
   } while (0)

/*!
 * Print a warning without exit.  Print file and line number of the warning.
 */
#define TBOX_WARNING(X)                        \
   do {                                        \
      std::ostringstream tboxos;               \
      tboxos << X;                             \
      SAMRAI::tbox::Logger::getInstance()->logWarning( \
         tboxos.str(), __FILE__, __LINE__);    \
   } while (0)

/*!
 * Print a debug without exit.  Print file and line number of the debug.
 */
#define TBOX_DEBUG(X)                        \
   do {                                      \
      std::ostringstream tboxos;             \
      tboxos << X;                           \
      SAMRAI::tbox::Logger::getInstance()->logDebug( \
         tboxos.str(), __FILE__, __LINE__);  \
   } while (0)

/*!
 * Throw an error assertion from within any C++ source code if the
 * given expression is not true.  This is a parallel-friendly version
 * of assert.
 * The file and line number of the abort are also printed.
 */
#ifdef DEBUG_CHECK_ASSERTIONS

#define TBOX_ASSERT(EXP)                                           \
   do {                                                            \
      if (!(EXP)) {                                                \
         std::ostringstream tboxos;                                \
         tboxos << "Failed assertion: " << # EXP;                  \
         SAMRAI::tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__); \
      }                                                            \
   } while (0)
#else

/*
 * No assertion checking
 */
#define TBOX_ASSERT(EXP)

#endif

#ifdef DEBUG_CHECK_ASSERTIONS

#define TBOX_CONSTEXPR_ASSERT(EXP)                                 \
   do {                                                            \
      if (!(EXP)) {                                                \
         printf("Failed assertion: ");                         \
         printf( # EXP );                                                \
         printf("\n");                                                   \
         SAMRAI::tbox::Utilities::abort("", __FILE__, __LINE__); \
      }                                                            \
   } while (0)
#else

/*
 * No assertion checking
 */
#define TBOX_CONSTEXPR_ASSERT(EXP)

#endif

/*!
 * Throw an error assertion from within any C++ source code if the
 * given expression is not true.  This is a parallel-friendly version
 * of assert.
 * The file and line number of the abort are also printed along with the
 * supplied message.
 */
#ifdef DEBUG_CHECK_ASSERTIONS

#define TBOX_ASSERT_MSG(EXP, MSG)                                  \
   do {                                                            \
      if (!(EXP)) {                                                \
         std::ostringstream tboxos;                                \
         tboxos << "Failed assertion: " << # EXP << std::endl << # MSG; \
         SAMRAI::tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__); \
      }                                                            \
   } while (0)
#else

/*
 * No assertion checking
 */
#define TBOX_ASSERT_MSG(EXP, MSG)

#endif

/*!
 */
#ifdef DEBUG_CHECK_ASSERTIONS

#define SAMRAI_SHARED_PTR_CAST std::dynamic_pointer_cast
#define CPP_CAST dynamic_cast

#else

#define SAMRAI_SHARED_PTR_CAST std::static_pointer_cast
#define CPP_CAST static_cast

#endif

/*!
 * Throw an error assertion from within any C++ source code if the
 * given expression is not true.  This version is used for assertions
 * that are useful checking internal SAMRAI state for developers working
 * on SAMRAI.  User level assertions should use TBOX_ASSERT.
 *
 * This is a parallel-friendly version of assert.  The file and line
 * number of the abort are also printed.
 */
#ifdef DEBUG_CHECK_DEV_ASSERTIONS

#define TBOX_DEV_ASSERT(EXP)                                       \
   do {                                                            \
      if (!(EXP)) {                                                \
         std::ostringstream tboxos;                                \
         tboxos << "Failed assertion: " << # EXP;                  \
         SAMRAI::tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__); \
      }                                                            \
   } while (0)
#else

/*
 * No SAMRAI internal developer assertion checking
 */
#define TBOX_DEV_ASSERT(EXP)

#endif

#define TBOX_ASSERT_OBJDIM_EQUALITY2(arg1, arg2) \
   TBOX_DIM_ASSERT(                              \
   ((arg1).getDim() == (arg2).getDim())          \
   )

#define TBOX_ASSERT_OBJDIM_EQUALITY3(arg1, arg2, arg3) \
   TBOX_DIM_ASSERT(                                    \
   ((arg1).getDim() == (arg2).getDim()) &&             \
   ((arg1).getDim() == (arg3).getDim())                \
   )

#define TBOX_ASSERT_OBJDIM_EQUALITY4(arg1, arg2, arg3, arg4) \
   TBOX_DIM_ASSERT(                                          \
   ((arg1).getDim() == (arg2).getDim()) &&                   \
   ((arg1).getDim() == (arg3).getDim()) &&                   \
   ((arg1).getDim() == (arg4).getDim())                      \
   )

#define TBOX_ASSERT_OBJDIM_EQUALITY5(arg1, arg2, arg3, arg4, arg5) \
   TBOX_DIM_ASSERT(                                                \
   ((arg1).getDim() == (arg2).getDim()) &&                         \
   ((arg1).getDim() == (arg3).getDim()) &&                         \
   ((arg1).getDim() == (arg4).getDim()) &&                         \
   ((arg1).getDim() == (arg5).getDim())                            \
   )

#define TTBOX_ASSERT_OBJDIM_EQUALITY6(arg1, arg2, arg3, arg4, arg5, arg6) \
   TBOX_DIM_ASSERT(                                                       \
   ((arg1).getDim() == (arg2).getDim()) &&                                \
   ((arg1).getDim() == (arg3).getDim()) &&                                \
   ((arg1).getDim() == (arg4).getDim()) &&                                \
   ((arg1).getDim() == (arg5).getDim()) &&                                \
   ((arg1).getDim() == (arg6).getDim())                                   \
   )

#define TBOX_ASSERT_OBJDIM_EQUALITY7(arg1, \
                                     arg2, \
                                     arg3, \
                                     arg4, \
                                     arg5, \
                                     arg6, \
                                     arg7) \
   TBOX_DIM_ASSERT(                        \
   ((arg1).getDim() == (arg2).getDim()) && \
   ((arg1).getDim() == (arg3).getDim()) && \
   ((arg1).getDim() == (arg4).getDim()) && \
   ((arg1).getDim() == (arg5).getDim()) && \
   ((arg1).getDim() == (arg6).getDim()) && \
   ((arg1).getDim() == (arg7).getDim())    \
   )

#define TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dim, arg1) \
   TBOX_DIM_ASSERT(                                 \
   ((dim) == (arg1).getDim())                       \
   )                                                \

#define TBOX_ASSERT_DIM_OBJDIM_EQUALITY2(dim, arg1, arg2) \
   TBOX_DIM_ASSERT(                                       \
   ((dim) == (arg1).getDim()) &&                          \
   ((dim) == (arg2).getDim())                             \
   )

#define TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, arg1, arg2, arg3) \
   TBOX_DIM_ASSERT(                                             \
   ((dim) == (arg1).getDim()) &&                                \
   ((dim) == (arg2).getDim()) &&                                \
   ((dim) == (arg3).getDim())                                   \
   )

#define TBOX_ASSERT_DIM_OBJDIM_EQUALITY4(dim, arg1, arg2, arg3, arg4) \
   TBOX_DIM_ASSERT(                                                   \
   ((dim) == (arg1).getDim()) &&                                      \
   ((dim) == (arg2).getDim()) &&                                      \
   ((dim) == (arg3).getDim()) &&                                      \
   ((dim) == (arg4).getDim())                                         \
   )

#define TBOX_ASSERT_DIM_OBJDIM_EQUALITY5(dim, arg1, arg2, arg3, arg4, arg5) \
   TBOX_DIM_ASSERT(                                                         \
   ((dim) == (arg1).getDim()) &&                                            \
   ((dim) == (arg2).getDim()) &&                                            \
   ((dim) == (arg3).getDim()) &&                                            \
   ((dim) == (arg4).getDim()) &&                                            \
   ((dim) == (arg5).getDim())                                               \
   )

#define TBOX_ASSERT_DIM_OBJDIM_EQUALITY6(dim,  \
                                         arg1, \
                                         arg2, \
                                         arg3, \
                                         arg4, \
                                         arg5, \
                                         arg6) \
   TBOX_DIM_ASSERT(                            \
   ((dim) == (arg1).getDim()) &&               \
   ((dim) == (arg2).getDim()) &&               \
   ((dim) == (arg3).getDim()) &&               \
   ((dim) == (arg4).getDim()) &&               \
   ((dim) == (arg5).getDim()) &&               \
   ((dim) == (arg6).getDim())                  \
   )

#define TBOX_ASSERT_DIM_OBJDIM_EQUALITY7(dim,  \
                                         arg1, \
                                         arg2, \
                                         arg3, \
                                         arg4, \
                                         arg5, \
                                         arg6, \
                                         arg7) \
   TBOX_DIM_ASSERT(                            \
   ((dim) == (arg1).getDim()) &&               \
   ((dim) == (arg2).getDim()) &&               \
   ((dim) == (arg3).getDim()) &&               \
   ((dim) == (arg4).getDim()) &&               \
   ((dim) == (arg5).getDim()) &&               \
   ((dim) == (arg6).getDim()) &&               \
   ((dim) == (arg7).getDim())                  \
   )

#define TBOX_ASSERT_DIM_OBJDIM_EQUALITY8(dim,  \
                                         arg1, \
                                         arg2, \
                                         arg3, \
                                         arg4, \
                                         arg5, \
                                         arg6, \
                                         arg7, \
                                         arg8) \
   TBOX_DIM_ASSERT(                            \
   ((dim) == (arg1).getDim()) &&               \
   ((dim) == (arg2).getDim()) &&               \
   ((dim) == (arg3).getDim()) &&               \
   ((dim) == (arg4).getDim()) &&               \
   ((dim) == (arg5).getDim()) &&               \
   ((dim) == (arg6).getDim()) &&               \
   ((dim) == (arg7).getDim()) &&               \
   ((dim) == (arg8).getDim())                  \
   )

#ifdef DEBUG_CHECK_DIM_ASSERTIONS

#define TBOX_DIM_ASSERT(EXP)                                             \
   do {                                                                  \
      if (!(EXP)) {                                                      \
         std::ostringstream tboxos;                                      \
         tboxos << "Failed dimension assertion: " << # EXP;              \
         SAMRAI::tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__);       \
      }                                                                  \
   } while (0)

#else

/*
 * No dimensional assertion checking
 */
#define TBOX_DIM_ASSERT(EXP)

#endif

#ifdef DEBUG_CHECK_DIM_ASSERTIONS

#define TBOX_CONSTEXPR_DIM_ASSERT(EXP)                                   \
   do {                                                                  \
      if (!(EXP)) {                                                      \
         printf("Failed dimension assertion: ");                         \
         printf( # EXP );                                                \
         printf("\n");                                                   \
         SAMRAI::tbox::Utilities::abort("", __FILE__, __LINE__);         \
      }                                                                  \
   } while (0)

#else

/*
 * No dimensional assertion checking
 */
#define TBOX_CONSTEXPR_DIM_ASSERT(EXP)

#endif

/**
 * Throw an error assertion from within any C++ source code.  This is
 * is similar to TBOX_ERROR(), but is designed to be invoked after a
 * call to a PETSc library function.  In other words, it acts similarly
 * to the PETSc CHKERRQ(ierr) macro.
 */
#ifdef HAVE_PETSC

/*
 * In the following, "CHKERRCONTINUE(ierr);" will cause PETSc to print out
 * a stack trace that led to the error; this may be useful for debugging.
 */

#define PETSC_SAMRAI_ERROR(ierr)                                   \
   do {                                                            \
      if (ierr) {                                                  \
         std::ostringstream tboxos;                                \
         SAMRAI::tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__); \
      }                                                            \
   } while (0)
#endif

#if defined(SAMRAI_HAVE_CALIPER)
#define SAMRAI_CALI_CXX_MARK_FUNCTION CALI_CXX_MARK_FUNCTION
#define SAMRAI_CALI_MARK_BEGIN(label) CALI_MARK_BEGIN(label)
#define SAMRAI_CALI_MARK_END(label)   CALI_MARK_END(label)
#else
#define SAMRAI_CALI_CXX_MARK_FUNCTION
#define SAMRAI_CALI_MARK_BEGIN(label)
#define SAMRAI_CALI_MARK_END(label)
#endif

/*
 * Macros defined for host-device compilation.
 */
#define SAMRAI_INLINE inline

#if defined(HAVE_CUDA) && defined(__CUDACC__)
#define SAMRAI_HOST_DEVICE __host__ __device__
#elif defined(HAVE_HIP) 
#define SAMRAI_HOST_DEVICE __host__ __device__
#else
#define SAMRAI_HOST_DEVICE
#endif


/*!
 * Utilities is a Singleton class containing basic routines for error
 * reporting, file manipulations, etc.
 */

struct Utilities {
   /*!
    * Create the directory specified by the path string.  Permissions are set
    * by default to rwx by user.  The intermediate directories in the
    * path are created if they do not already exist.  When
    * only_node_zero_creates is true, only node zero creates the
    * directories.  Otherwise, all nodes create the directories.
    */
   static void
   recursiveMkdir(
      const std::string& path,
      mode_t mode = (S_IRUSR | S_IWUSR | S_IXUSR),
      bool only_node_zero_creates = true);

   /*!
    * Rename a file from old file name to new file name.
    *
    * @pre !old_filename.empty()
    * @pre !new_filename.empty()
    */
   static void
   renameFile(
      const std::string& old_filename,
      const std::string& new_filename)
   {
      TBOX_ASSERT(!old_filename.empty());
      TBOX_ASSERT(!new_filename.empty());
      rename(old_filename.c_str(), new_filename.c_str());
   }

   /*!
    * Convert an integer to a string.
    *
    * The returned string is padded with zeros as needed so that it
    * contains at least the number of characters indicated by the
    * minimum width argument.  When the number is positive, the
    * string is padded on the left. When the number is negative,
    * the '-' sign appears first, followed by the integer value
    * padded on the left with zeros.  For example, the statement
    * intToString(12, 5) returns "00012" and the statement
    * intToString(-12, 5) returns "-0012".
    */
   static std::string
   intToString(
      int num,
      int min_width = 1);

   /*!
    * Convert a size_t to a string.
    *
    * The returned string is padded with zeros as needed so that it
    * contains at least the number of characters indicated by the
    * minimum width argument.  When the number is positive, the
    * string is padded on the left. When the number is negative,
    * the '-' sign appears first, followed by the integer value
    * padded on the left with zeros.  For example, the statement
    * intToString(12, 5) returns "00012" and the statement
    * intToString(-12, 5) returns "-0012".
    */
   static std::string
   sizetToString(
      size_t num,
      int min_width = 1);

   /*!
    * Convert common integer values to strings.
    *
    * These are simply wrappers around intToString that ensure the
    * same width is uniformally used when converting to string
    * representations.
    */
   static std::string
   nodeToString(
      int num)
   {
      return intToString(num, s_node_width);
   }
   static std::string
   processorToString(
      int num)
   {
      return intToString(num, s_processor_width);
   }
   static std::string
   patchToString(
      int num)
   {
      return intToString(num, s_patch_width);
   }
   static std::string
   levelToString(
      int num)
   {
      return intToString(num, s_level_width);
   }
   static std::string
   blockToString(
      int num)
   {
      return intToString(num, s_block_width);
   }

   /*!
    * Aborts the run after printing an error message with file and
    * linenumber information.
    */
   static void
   abort(
      const std::string& message,
      const std::string& filename,
      const int line);

private:
   /*
    * Sizes for converting integers to fixed width strings
    * for things like filenames etc.
    */
   static const int s_node_width = 7;
   static const int s_processor_width = 7;
   static const int s_patch_width = 7;
   static const int s_level_width = 4;
   static const int s_block_width = 7;

};

}
}

#endif
