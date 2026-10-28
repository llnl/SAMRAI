/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Parser that reads the input database grammar
 *
 ************************************************************************/

#include "SAMRAI/tbox/Parser.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/PIO.h"

#ifdef __INTEL_COMPILER
// Ignore Intel warnings about external declarations
#pragma warning (disable:1419)
#endif

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

extern int
yyparse();
extern void
yyrestart(
   FILE *);

extern void
parser_static_table_initialize();

namespace SAMRAI {
namespace tbox {

Parser * Parser::s_default_parser = 0;
bool Parser::s_static_tables_initialized = 0;

/*
 *************************************************************************
 *
 * The constructor creates an unitialized parser object.  All of the
 * interesting work is done by member function parse().
 *
 *************************************************************************
 */

Parser::Parser()
   :
   d_read_locally(false)
{
   if (!s_static_tables_initialized) {
      parser_static_table_initialize();
      s_static_tables_initialized = 1;
   }
}

/*
 *************************************************************************
 *
 * The destructor automatically deallocates the parser object data.
 *
 *************************************************************************
 */

Parser::~Parser()
{
}

/*
 *************************************************************************
 *
 * Begin parsing the input database file.  Return the number of errors
 * encountered in the parse.
 *
 *************************************************************************
 */

int
Parser::parse(
   const std::string& filename,
   FILE* fstream,
   const std::shared_ptr<Database>& database,
   bool read_locally)
{
   d_errors = 0;
   d_warnings = 0;
   d_read_locally = read_locally;

   // Find the path in the filename, if one exists
   std::string::size_type slash_pos = filename.find_last_of('/');
   if (slash_pos == std::string::npos) {
      d_pathname = "";
   } else {
      d_pathname = filename.substr(0, slash_pos + 1);
   }

   ParseData pd;
   pd.d_filename = filename;
   pd.d_fstream = fstream;
   pd.d_linenumber = 1;
   pd.d_cursor = 1;
   pd.d_nextcursor = 1;
   d_parse_stack.clear();
   d_parse_stack.push_front(pd);

   d_scope_stack.clear();
   d_scope_stack.push_front(database);

   s_default_parser = this;
   yyrestart(0);
   if (yyparse() && (d_errors == 0)) {
      error("Unexpected parse error");
   }
   s_default_parser = 0;

   d_parse_stack.clear();
   d_scope_stack.clear();

   return d_errors;
}

/*
 *************************************************************************
 *
 * Advance the cursor to the next line in the current input file.
 *
 *************************************************************************
 */

void
Parser::advanceLine(
   const int nline)
{
   Parser::ParseData& pd = d_parse_stack.front();
   pd.d_linenumber += nline;
   pd.d_cursor = 1;
   pd.d_nextcursor = 1;
}

/*
 *************************************************************************
 *
 * Advance the cursor position by the token in the specified string.
 * Tabs are expanded assuming tab stops at eight character markers.
 *
 *************************************************************************
 */

void
Parser::advanceCursor(
   const std::string& token)
{
   Parser::ParseData& pd = d_parse_stack.front();
   pd.d_cursor = pd.d_nextcursor;
   for (std::string::const_iterator i = token.begin(); i != token.end(); ++i) {
      if (*i == '\t') {
         pd.d_nextcursor = ((pd.d_nextcursor + 7) & (~7)) + 1;
      } else {
         ++(pd.d_nextcursor);
      }
   }
}

/*
 *************************************************************************
 *
 * Print out errors to pout and track the number of errors.
 *
 *************************************************************************
 */

void
Parser::error(
   const std::string& message)
{
   Parser::ParseData& pd = d_parse_stack.front();

   pout << "Error in " << pd.d_filename << " at line " << pd.d_linenumber
        << " column " << pd.d_cursor
        << " : " << message << std::endl << std::flush;

   pout << pd.d_linebuffer << std::endl << std::flush;

   for (int i = 0; i < pd.d_cursor; ++i)
      pout << " ";
   pout << "^\n";

   ++d_errors;
}

/*
 *************************************************************************
 *
 * Print out warnings to pout and track the number of warnings.
 *
 *************************************************************************
 */

void
Parser::warning(
   const std::string& message)
{
   Parser::ParseData& pd = d_parse_stack.front();

   pout << "Warning in " << pd.d_filename << " at line " << pd.d_linenumber
        << " column " << pd.d_cursor
        << " : " << message << std::endl << std::flush;

   pout << pd.d_linebuffer << std::endl << std::flush;

   for (int i = 0; i < pd.d_cursor; ++i)
      pout << " ";
   pout << "^\n";

   ++d_warnings;
}

/*
 *************************************************************************
 *
 * Iterate through the database scopes, looking for the first match on
 * the key value.
 *
 *************************************************************************
 */

std::shared_ptr<Database>
Parser::getDatabaseWithKey(
   const std::string& key)
{
   std::list<std::shared_ptr<Database> >::iterator i = d_scope_stack.begin();
   for ( ; i != d_scope_stack.end(); ++i) {
      if ((*i)->keyExists(key)) {
         return *i;
      }
   }
   return std::shared_ptr<Database>();
}

/*
 *************************************************************************
 *
 * Create a new parse state on the parse stack and open the specified
 * new file for reading.
 *
 *************************************************************************
 */

bool
Parser::pushIncludeFile(
   const std::string& filename)
{
   FILE* fstream = 0;
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());

   std::string filename_with_path;

   // If this is not a fully qualified pathname use
   // current search path
   std::string::size_type slash_pos;
   slash_pos = filename.find_first_of('/');
   if (slash_pos == 0) {
      filename_with_path = filename;
   } else {
      filename_with_path = d_pathname;
      filename_with_path += filename;
   }

   if (d_read_locally || mpi.getRank() == 0) {
      fstream = fopen(filename_with_path.c_str(), "r");
   }

   int worked = (fstream ? 1 : 0);

   if (!d_read_locally) {
      mpi.Bcast(&worked, 1, MPI_INT, 0);
   }

   if (!worked) {
      error("Could not open include file ``" + filename_with_path + "''");
   } else {
      ParseData pd;
      pd.d_filename = filename_with_path;
      pd.d_fstream = fstream;
      pd.d_linenumber = 1;
      pd.d_cursor = 1;
      pd.d_nextcursor = 1;
      d_parse_stack.push_front(pd);
   }

   return worked ? true : false;
}

/*
 *************************************************************************
 *
 * Close the current input file and pop the parse stack.
 *
 *************************************************************************
 */

void
Parser::popIncludeFile()
{
   Parser::ParseData& pd = d_parse_stack.front();
   if (pd.d_fstream) {
      fclose(pd.d_fstream);
   }
   d_parse_stack.pop_front();
}

/*
 *************************************************************************
 *
 * Manage the input reading for the flex scanner.  If running with MPI,
 * the node zero reads the data and broadcasts the length and the data
 * to all processors.
 *
 *************************************************************************
 */

int
Parser::yyinput(
   char* buffer,
   const int max_size)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   int byte = 0;
   if (d_read_locally || mpi.getRank() == 0) {
      byte = static_cast<int>(fread(buffer,
                                 1,
                                 max_size,
                                 d_parse_stack.front().d_fstream));
   }
   if (!d_read_locally) {
      mpi.Bcast(&byte, 1, MPI_INT, 0);
      if (byte > 0) {
         mpi.Bcast(buffer, byte, MPI_CHAR, 0);
      }
   }
   return byte;
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
