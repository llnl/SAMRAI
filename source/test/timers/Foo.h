/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Simple example to demonstrate input/restart of patch data.
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include <string>

#define included_String
#include "SAMRAI/tbox/Serializable.h"

using namespace SAMRAI;

class Foo
{
public:
   Foo();
   ~Foo();

   void
   timerOff();
   void
   timerOn();
   void
   zero(
      int depth);
   void
   one(
      int depth);
   void
   two(
      int depth);
   void
   three(
      int depth);
   void
   four(
      int depth);
   void
   five(
      int depth);
   void
   six(
      int depth);
   void
   seven(
      int depth);
   void
   startAndStop(
      std::string& name);
   void
   setMaxDepth(
      int max_depth);
   void
   start(
      std::string& name);
   void
   stop(
      std::string& name);

private:
   int d_depth;
   int d_max_depth;

};
