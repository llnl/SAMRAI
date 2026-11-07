/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility for getting input file name.
 *
 ************************************************************************/
#include <string>

/*!
 * This function gets the input file name from the argument list.
 * For most cases, this is simply the first
 * command-line argument (after the program name).
 * The variable @c input_filename is set to the filename,
 * and all arguments following it are shifed down one.
 * Unless you are a developer, this is all you need to know.
 *
 * For the special case with no command-line argument,
 * if the program is invoked with the name "check-blah",
 * it assumes the input file is "blah.input".
 * The purpose of this setup is so that it is compatible with
 * the automake make-check mechanism, which runs programs without
 * arguments.  At the same time, an input file can be manually
 * specified for operational flexibility.
 *
 * @return 0 if filename is found, 1 otherwise.
 */
int get_input_filename(int* argc, char* argv[], std::string & input_filename);
