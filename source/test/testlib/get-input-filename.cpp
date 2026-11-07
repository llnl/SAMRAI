/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility for getting input file name.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"
#include <string>

#include "test/testlib/get-input-filename.h"

int get_input_filename(
   int* argc,
   char* argv[],
   std::string& input_filename) {

   int rval = 0;
   std::string argv0(argv[0]);
   if (*argc > 1) {
      // Input file is the first argument.  Shift other arguments down.
      input_filename = argv[1];
      --(*argc);
      int i;
      for (i = 1; i < (*argc); ++i) {
         argv[i] = argv[i + 1];
      }
   } else if (*argc == 1 && argv0.rfind("check-") < argv0.size()) {
      /*
       * No argument but input file is implicit in program name
       * which has the form check-<input file without .input>.
       */
      input_filename = argv0.substr(argv0.rfind("check-") + 6) + ".input";
   } else if (*argc == 1) {
      // No argument and not invoked as "check-blah".
      rval = 1;
   }
   return rval;
}
