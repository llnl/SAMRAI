/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A simple call sequence tracking class
 *
 ************************************************************************/

#ifndef included_tbox_Tracer
#define included_tbox_Tracer

#include "SAMRAI/SAMRAI_config.h"

#include <iostream>
#include <string>

namespace SAMRAI {
namespace tbox {

/**
 * Class Tracer allows one to trace entrances and exits of
 * class member functions. An example usage is:
 * \code
 * #include "SAMRAI/tbox/Tracer.h"
 * ....
 * void MyClass::myClassMemberFunction()
 * {
 *    Tracer t("MyClass::myClassMemberFunction");
 *     ....
 * }
 * \endcode
 * When the function `myClassMemberFunction' is called, a tracer
 * object  local to the function scope is created and the message
 * {\em Entering MyClass::myClassMemberFunction} is sent to the
 * tracer output stream. Upon exiting the function, the tracer object
 * is destroyed and {\em Exiting MyClass::myClassMemberFunction}
 * is sent to the tracer output stream.  By default, the tracer
 * class sends data to the parallel log stream, although the default
 * output stream can be changed through a call to static member function
 * setTraceStream(), which will set the tracer output stream for all
 * subsequent calls.
 */

class Tracer
{
public:
   /**
    * The constructor for Tracer prints ``Entering \<message\>''
    * to the tracer output stream.
    */
   explicit Tracer(
      const std::string& message);

   /**
    * The destructor for Tracer prints ``Exiting \<message\>''
    * to the tracer output stream.
    */
   ~Tracer();

   /**
    * Set the tracer output stream for all tracer output.  By default,
    * this is set to the parallel log stream plog.  If this argument is
    * NULL, then all output to trace streams is disabled.
    */
   static void
   setTraceStream(
      std::ostream* stream)
   {
      s_stream = stream;
   }

private:
   // Unimplemented default constructor.
   Tracer();

   // Unimplemented copy constructor.
   Tracer(
      const Tracer& other);

   // Unimplemented assignment operator.
   Tracer&
   operator = (
      const Tracer& rhs);

   std::string d_message;
   static std::ostream* s_stream;
};

}
}

#endif
