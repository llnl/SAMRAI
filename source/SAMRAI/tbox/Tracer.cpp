/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A simple call sequence tracking class
 *
 ************************************************************************/

#include "SAMRAI/tbox/Tracer.h"
#include "SAMRAI/tbox/PIO.h"

namespace SAMRAI {
namespace tbox {

std::ostream * Tracer::s_stream = &plog;

Tracer::Tracer(
   const std::string& message)
{
   d_message = message;
   if (s_stream) {
      (*s_stream) << "Entering " << d_message << std::endl << std::flush;
   }
}

Tracer::~Tracer()
{
   if (s_stream) {
      (*s_stream) << "Exiting " << d_message << std::endl << std::flush;
   }
}

}
}
