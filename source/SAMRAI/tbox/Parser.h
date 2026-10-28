/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Parser that reads the input database grammar
 *
 ************************************************************************/

#ifndef included_tbox_Parser
#define included_tbox_Parser

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/Database.h"

#include <cstdio>
#include <string>
#include <list>
#include <memory>

namespace SAMRAI {
namespace tbox {

/**
 * Class Parser parses the user input file and places the resulting
 * (key,value) pairs into the input database object.  The parser object
 * controls the overall parsing of the input file, which includes error
 * handing and tracking file, line number, and cursor position.  If
 * running on multiple processors, only node zero reads in data from the
 * specified input file and broadcasts that data to the other processors.
 * The input file argument for the other processors is ignored and may be
 * NULL.
 *
 * The parser class also defines a ``default'' parser that may be accessed
 * via a static member function.  The default parser may only be accessed
 * by the yacc/lex routines called during the input file parsing.  This
 * singleton-like approach provides a clean way to communicate parser
 * information to the yacc/lex files without global variables.
 *
 * This parser (and the associated yacc and lex files) assume the GNU flex
 * and bison utilities.  These utilities are only required if the grammar or
 * scanner is changed, since the SAMRAI distribution already contains the
 * output files from flex and bison.
 */

class Parser
{
public:
   /**
    * The parser constructor simply creates an uninitialized parser object.
    * Member function parse() must be called before any other member function
    * to initialize the object and parse the input data.  Function parse()
    * may be called multiple times to parse multiple input files, but all
    * state values (such as the number of errors or warnings) are reset at
    * the beginning of each new parse pass.
    */
   Parser();

   /**
    * Destroy the parser object and deallocate parser data.
    */
   ~Parser();

   /**
    * Parse the input file from the specified file stream.  The number of
    * syntax errors is returned.  A successful parse will return zero errors.
    * The parse() function takes the initial filename (for informational
    * purposes) and the filestream from which to read the parse data.
    * All (key,value) pairs are placed in the specified database.  If
    * running in parallel, the fstream must be valid on node zero, but
    * is ignored on other nodes and may be set to NULL.  Multiple input
    * files may be parsed by calling parse() for each file, but all variables
    * are reset at the beginning of each parse.
    *
    * When @p read_locally is true, every rank reads its own file stream.
    * This is intended for rank-local auxiliary data; normal application
    * input should retain the default collective behavior.
    *
    * @param filename Name of the input file
    * @param fstream Open input stream.  In collective mode, it is required
    *                only on rank zero.
    * @param database Database that receives the parsed key-value pairs.
    * @param read_locally When true, each rank reads from its own @p fstream
    *                     without broadcasting input from rank zero.
    */
   int
   parse(
      const std::string& filename,
      FILE* fstream,
      const std::shared_ptr<Database>& database,
      bool read_locally = false);

   /**
    * Return the total number of errors resulting from the parse.
    */
   int
   getNumberErrors() const
   {
      return d_errors;
   }

   /**
    * Return the total number of warnings resulting from the parse.
    */
   int
   getNumberWarnings() const
   {
      return d_warnings;
   }

   /**
    * Return the parser object.  This mechanism is useful for communicating
    * with the yacc/lex routines during the input file parse.  The default
    * parser will be NULL outside of the parse call.
    */
   static Parser *
   getParser()
   {
      return s_default_parser;
   }

   /**
    * Return the current database scope.  The current scope is modified
    * through the enterScope() and leaveScope() member functions.
    */
   std::shared_ptr<Database>&
   getScope()
   {
      return d_scope_stack.front();
   }

   /**
    * Create a new database scope with the specified name.  This new scope
    * will be the default scope until leaveScope() is called.
    */
   void
   enterScope(
      const std::string& name)
   {
      d_scope_stack.push_front(d_scope_stack.front()->putDatabase(name));
   }

   /**
    * Leave the current database scope and return to the previous scope.
    * It is an error to leave the outermost scope.
    */
   void
   leaveScope()
   {
      d_scope_stack.pop_front();
   }

   /**
    * Lookup the scope that contains the specified key.  If the scope does
    * not exist, then return a NULL pointer to the database.
    */
   std::shared_ptr<Database>
   getDatabaseWithKey(
      const std::string& name);

   /**
    * Save the current context and switch to the specified input file.
    * This routine returns true if the file exists and the switch was
    * successful and false otherwise.
    */
   bool
   pushIncludeFile(
      const std::string& filename);

   /**
    * Pop the include file context off of the stack and return to the
    * previous include file.
    */
   void
   popIncludeFile();

   /**
    * Report a parsing error with the specified error message.  This routine
    * will only be called from the parser or the scanner.  Errors are printed
    * to pout, since it is assumed that all nodes are parsing the same input
    * file.
    */
   void
   error(
      const std::string& message);

   /**
    * Report a parsing warning with the specified warning message.  This
    * routine will only be called from the parser or the scanner.  Errors
    * are printed to pout, since it is assumed that all nodes are parsing
    * the same input file.
    */
   void
   warning(
      const std::string& message);

   /**
    * Set the input line which is currently being parsed.
    */
   void
   setLine(
      const std::string& line)
   {
      Parser::ParseData& pd = d_parse_stack.front();
      pd.d_linebuffer = line;
   }

   /**
    * Advance the line number by the specified number of lines.  If no
    * argument is given, then the line number is advanced by one line.
    */
   void
   advanceLine(
      const int nline = 1);

   /**
    * Advance the position of the cursor on the line using the values
    * in the specified character string.  Tab characters in the string
    * are assumed to advance the cursor to eight character tab stops.
    * The cursor position is automatically reset to one whenever the
    * line number is changed.
    */
   void
   advanceCursor(
      const std::string& token);

   /**
    * Define the input reading routine used by flex.  Under MPI, node zero
    * reads the input and broadcasts the character data to all processors.
    */
   int
   yyinput(
      char * buffer,
      const int max_size);

private:
   Parser(
      const Parser&);           // not implemented
   Parser&
   operator = (
      const Parser&);           // not implemented

   struct ParseData {
      std::string d_filename;   // filename for description
      FILE* d_fstream;          // input stream to parse
      std::string d_linebuffer; // line being parsed
      int d_linenumber;         // line number in input stream
      int d_cursor;             // cursor position in line
      int d_nextcursor;         // next cursor position in line
   };

   int d_errors;                // total number of parse errors
   int d_warnings;              // total number of warnings

   std::list<Parser::ParseData> d_parse_stack;

   std::list<std::shared_ptr<Database> > d_scope_stack;

   static Parser* s_default_parser;

   static bool s_static_tables_initialized;

   std::string d_pathname;           // path to filename for including
   bool d_read_locally;
};

}
}

#endif
