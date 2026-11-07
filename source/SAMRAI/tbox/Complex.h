/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   dcomplex class for old-style complex and new complex<double>
 *
 ************************************************************************/

#ifndef included_tbox_Complex
#define included_tbox_Complex

#include "SAMRAI/SAMRAI_config.h"

#include <complex>

/*!
 * @page toolbox_complex Toolbox Complex Type
 *
 * @brief dcomplex is a typedef to overcome C++ compiler issues with
 * the std::complex type.
 *
 * The std::complex type should be a template however some older C++ compilers
 * implement complex as a double complex.  dcomplex is used to hide this
 * platform issue behind a typedef.
 *
 * @internal NOTE: This should be removed when no longer required.
 *
 */

#ifndef LACKS_TEMPLATE_COMPLEX
typedef std::complex<double> dcomplex;
#else
typedef std::complex dcomplex;
#endif

#endif
