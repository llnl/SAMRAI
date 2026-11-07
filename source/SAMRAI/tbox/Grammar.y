%{
//
// This file is part of the SAMRAI distribution.  For full copyright
// information, see COPYRIGHT and LICENSE.
//
// Copyright:   (c) 1997-2026 Lawrence Livermore National Security, LLC
// Description: Yacc grammar description for the input database
//

#include "SAMRAI/SAMRAI_config.h"
#include <math.h>

#include STL_SSTREAM_HEADER_FILE


#if !defined(OSTRINGSTREAM_TYPE_IS_BROKEN) && defined(OSTRSTREAM_TYPE_IS_BROKEN)
typedef ostringstream ostrstream;
#endif

#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Parser.h"
#include <string>
#include <memory>

#ifdef __xlC__
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif


using namespace SAMRAI;
using namespace tbox;

extern int yylex();
void yyerror(const char *const error)
{
   Parser::getParser()->error(error);

}

// Do not change the numbering of keys without checking promotion logic

#define KEY_COMPLEX (0)
#define KEY_DOUBLE  (1)
#define KEY_INTEGER (2)
#define KEY_BOOL    (3)
#define KEY_BOX     (4)
#define KEY_CHAR    (5)
#define KEY_STRING  (6)

static std::string type_names[] = {
   "complex", "double", "int", "bool", "box", "char", "string"
};

#define IS_NUMBER(X) (((X) >= 0) && ((X) < KEY_BOOL))
#define PROMOTE(X,Y) ((X) < (Y) ? (X) : (Y))

struct KeyData
{
   int             d_node_type;	// KEYDATA node type (see defines above)
   int             d_array_type;// array type (numbers may be promoted)
   int             d_array_size;// total size of the array if head element
   KeyData*   d_next;	// pointer to next (key,data) pair
   bool            d_bool;	// boolean if node is KEY_BOOL
   DatabaseBox     d_box;	// box if node is KEY_BOX
   char            d_char;	// character if node is KEY_CHAR
   dcomplex        d_complex;	// complex if node is KEY_COMPLEX
   double          d_double;	// double if node is KEY_DOUBLE
   int             d_integer;	// integer if node is KEY_INTEGER
   std::string     d_string;	// string if node is KEY_STRING
};

static void delete_list(KeyData*);
static void to_boolean(KeyData*);
static void to_integer(KeyData*);
static void to_double(KeyData*);
static void to_complex(KeyData*);
static KeyData* binary_op(KeyData*, KeyData*, const int);
static KeyData* compare_op(KeyData*, KeyData*, const int);
static KeyData* eval_function(KeyData*, const std::string&);
static KeyData* lookup_variable(const std::string&, const int, const bool);

%}

%union 
{
  char          u_char;
  double        u_double;
  int           u_integer;
  KeyData* u_keydata;
  std::string*       u_keyword;
  std::string*       u_string;
}

%token             T_AND
%token             T_ASSIGN
%token <u_char>    T_CHAR
%token             T_CLOSE_CURLY
%token             T_CLOSE_PAREN
%token             T_CLOSE_SQBKT
%token             T_COMMA
%token             T_DIV
%token <u_double>  T_DOUBLE
%token             T_ELSE
%token             T_EXP
%token             T_EQUALS
%token             T_GREATER_EQUALS
%token             T_GREATER
%token             T_LESS_EQUALS
%token             T_LESS
%token             T_FALSE
%token <u_integer> T_INTEGER
%token <u_keyword> T_KEYWORD
%token             T_MINUS
%token             T_MULT
%token             T_NOT
%token             T_NOT_EQUALS
%token             T_OR
%token             T_OPEN_CURLY
%token             T_OPEN_PAREN
%token             T_OPEN_SQBKT
%token             T_PLUS
%token             T_QUESTION
%token             T_SEMI
%token <u_string>  T_STRING
%token             T_TRUE

%type  <u_keydata> P_BOX
%type  <u_keydata> P_COMPLEX
%type  <u_keydata> P_EXPRESSION
%type  <u_keydata> P_EXPRESSION_LIST
%type  <u_keydata> P_INTEGER_VECTOR
%type  <u_keydata> P_PRIMITIVE_TYPE

%left  T_OR
%left  T_AND
%left  T_NOT
%left  T_EQUALS T_NOT_EQUALS T_GREATER_EQUALS T_GREATER T_LESS_EQUALS T_LESS
%left  T_MINUS T_PLUS
%left  T_MULT  T_DIV
%right T_EXP
%nonassoc T_NEGATION

%%

/*
 * The specification is the grammar start state.  An input file consists
 * of a list of definitions.
 */

P_SPECIFICATION
 : P_DEFINITION_LIST 
 ;

/*
 * A definition list is zero or more definitions.
 */

P_DEFINITION_LIST
 : /* empty */
 | P_DEFINITION_LIST P_DEFINITION
 ;

/*
 * A definition is either a new database scope or a (key,value) pair where
 * the value is an array (and an array may only have one element).  In the
 * latter case, the array data is entered into the input database.
 */

P_DEFINITION
 : T_KEYWORD T_OPEN_CURLY {

   /* This is a hack to make a warning message go away from
      a symbol flex defines but does not use */
   if(0) {
      goto yyerrlab1;
   }

      if (Parser::getParser()->getScope()->keyExists(*$1)) {
	 std::string tmp("Redefinition of key ``");
         tmp += *$1;
         tmp += "''";
         Parser::getParser()->warning(tmp);
      }
      Parser::getParser()->enterScope(*$1);
   } P_DEFINITION_LIST T_CLOSE_CURLY {
      Parser::getParser()->leaveScope();
      delete $1;
   }
 | T_KEYWORD T_ASSIGN {
      if (Parser::getParser()->getScope()->keyExists(*$1)) {
	 std::string tmp("Redefinition of key ``");
         tmp += *$1;
         tmp += "''";
         Parser::getParser()->warning(tmp);
      }
   } P_EXPRESSION_LIST {
      KeyData* list = $4;
      const int n = list->d_array_size;

      switch (list->d_array_type) {
         case KEY_BOOL: {
            std::vector<bool> data(n);
            for (int i = n-1; i >= 0; i--) {
               data[i] = list->d_bool;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putBoolVector(*$1, data);
            break;
         }
         case KEY_BOX: {
            std::vector<DatabaseBox> data(n);
            for (int i = n-1; i >= 0; i--) {
               data[i] = list->d_box;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putDatabaseBoxVector(*$1, data);
            break;
         }
         case KEY_CHAR: {
            std::vector<char> data(n);
            for (int i = n-1; i >= 0; i--) {
               data[i] = list->d_char;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putCharVector(*$1, data);
            break;
         }
         case KEY_COMPLEX: {
            std::vector<dcomplex> data(n);
            for (int i = n-1; i >= 0; i--) {
               to_complex(list);
               data[i] = list->d_complex;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putComplexVector(*$1, data);
            break;
         }
         case KEY_DOUBLE: {
            std::vector<double> data(n);
            for (int i = n-1; i >= 0; i--) {
               to_double(list);
               data[i] = list->d_double;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putDoubleVector(*$1, data);
            break;
         }
         case KEY_INTEGER: {
            std::vector<int> data(n);
            for (int i = n-1; i >= 0; i--) {
               data[i] = list->d_integer;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putIntegerVector(*$1, data);
            break;
         }
         case KEY_STRING: {
            std::vector<std:string> data(n);
            for (int i = n-1; i >= 0; i--) {
               data[i] = list->d_string;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putStringVector(*$1, data);
            break;
         }
         default:
            Parser::getParser()->error("Internal parser error!");
            break;
      }

      delete_list($4);
      delete $1;
   }
 | T_SEMI {
      Parser::getParser()->warning(
         "Semicolon found in keyword phrase (ignored)");
   }
 ;

/*
 * Parse a list of expressions.  As each new expression is added to the
 * list, verify that the type of the new expression is compatible with the
 * type of the list.  If not, then print an error message and do not add
 * the expression to the list.  The list size is incremented as each new
 * expression is added to the list.  The final expression list is returned
 * in reverse order since this is a left-recursive grammar.
 */

P_EXPRESSION_LIST
 : P_EXPRESSION {
      $$ = $1;
   }
 | P_EXPRESSION_LIST T_COMMA P_EXPRESSION {
      switch($1->d_array_type) {
         case KEY_BOOL:
         case KEY_CHAR:
         case KEY_STRING:
            if ($3->d_node_type != $1->d_array_type) {
               Parser::getParser()->error("Type mismatch in array");
               delete $3;
               $$ = $1;
            } else {
               $3->d_array_size = $1->d_array_size + 1;
               $3->d_next       = $1;
               $$               = $3;
            }
            break;
         case KEY_BOX:
            if ($3->d_node_type != KEY_BOX) {
               Parser::getParser()->error("Type mismatch in box array");
               delete $3;
               $$ = $1;
            } else if ($3->d_box.getDimVal() != $1->d_box.getDimVal()) {
               Parser::getParser()->error("Box array dimension mismatch");
               delete $3;
               $$ = $1;
            } else {
               $3->d_array_size = $1->d_array_size + 1;
               $3->d_next       = $1;
               $$               = $3;
            }
            break;
         case KEY_COMPLEX:
         case KEY_DOUBLE:
         case KEY_INTEGER:
            if (!IS_NUMBER($3->d_node_type)) {
               Parser::getParser()->error("Type mismatch in number array");
               delete $3;
               $$ = $1;
            } else {
               $3->d_array_type = PROMOTE($1->d_array_type, $3->d_node_type);
               $3->d_array_size = $1->d_array_size + 1;
               $3->d_next       = $1;
               $$               = $3;
            }
            break;
      }
   }
 ;

/*
 * Parse a simple expression grammar for the input files.  Expressions are
 * one of the following:
 *
 *	- E = ( bool | box | char | complex | double | int | string )
 *	- E = variable
 *	- E = array[n]
 *	- E = ( E1 )
 *	- E = func(E1)
 *	- E = ( E1 ? E2 : E3 )
 *	- E = !E1
 *	- E = E1 || E2
 *	- E = E1 && E2
 *	- E = E1 == E2
 *	- E = E1 != E2
 *	- E = E1 < E2
 *	- E = E1 <= E2
 *	- E = E1 > E2
 *	- E = E1 >= E2
 *	- E = E1 + E2
 *	- E = E1 - E2
 *	- E = E1 * E2
 *	- E = E1 / E2
 *	- E = E1 ^ E2
 *	- E = -E1
 */

P_EXPRESSION
 : P_PRIMITIVE_TYPE {
      $$ = $1;
   }
 | T_KEYWORD {
      $$ = lookup_variable(*$1, 0, false);
      delete $1;
   }
 | T_KEYWORD T_OPEN_SQBKT P_EXPRESSION T_CLOSE_SQBKT {
      to_integer($3);
      $$ = lookup_variable(*$1, $3->d_integer, true);
      delete $1;
      delete $3;
   }
 | T_OPEN_PAREN P_EXPRESSION T_CLOSE_PAREN {
      $$ = $2;
   }
 | T_KEYWORD T_OPEN_PAREN P_EXPRESSION T_CLOSE_PAREN {
      $$ = eval_function($3, *$1);
      delete $1;
   }
 | T_OPEN_PAREN P_EXPRESSION T_QUESTION {
      if ($2->d_node_type != KEY_BOOL) {
         Parser::getParser()->error("X in (X ? Y : Z) is not a boolean");
      }
   } P_EXPRESSION T_ELSE P_EXPRESSION T_CLOSE_PAREN {
      if (($2->d_node_type == KEY_BOOL) && ($2->d_bool)) {
         $$ = $5;
         delete $7;
      } else {
         $$ = $7;
         delete $5;
      }
      delete $2;
   }
 | T_NOT P_EXPRESSION {
      to_boolean($2);
      $2->d_bool = !$2->d_bool;
      $$ = $2;
   }
 | P_EXPRESSION T_OR P_EXPRESSION {
      to_boolean($1);
      to_boolean($3);
      $1->d_bool = $1->d_bool || $3->d_bool;
      delete $3;
      $$ = $1;
   }
 | P_EXPRESSION T_AND P_EXPRESSION {
      to_boolean($1);
      to_boolean($3);
      $1->d_bool = $1->d_bool && $3->d_bool;
      delete $3;
      $$ = $1;
   }
 | P_EXPRESSION T_EQUALS P_EXPRESSION {
      $$ = compare_op($1, $3, T_EQUALS);
   }
 | P_EXPRESSION T_NOT_EQUALS P_EXPRESSION {
      $$ = compare_op($1, $3, T_EQUALS);
      $$->d_bool = !($$->d_bool);
   }
 | P_EXPRESSION T_GREATER_EQUALS P_EXPRESSION {
      $$ = compare_op($1, $3, T_LESS);
      $$->d_bool = !($$->d_bool);
   }
 | P_EXPRESSION T_GREATER P_EXPRESSION {
      $$ = compare_op($1, $3, T_GREATER);
   }
 | P_EXPRESSION T_LESS_EQUALS P_EXPRESSION {
      $$ = compare_op($1, $3, T_GREATER);
      $$->d_bool = !($$->d_bool);
   }
 | P_EXPRESSION T_LESS P_EXPRESSION {
      $$ = compare_op($1, $3, T_LESS);
   }
 | P_EXPRESSION T_PLUS P_EXPRESSION {
      if (($1->d_node_type == KEY_STRING) && ($3->d_node_type == KEY_STRING)) {
	 std::string tmp($1->d_string);
	 tmp += $3->d_string;
         $1->d_string = tmp;
         delete $3;
         $$ = $1;
      } else {
         $$ = binary_op($1, $3, T_PLUS);
      }
   }
 | P_EXPRESSION T_MINUS P_EXPRESSION {
      $$ = binary_op($1, $3, T_MINUS);
   }
 | P_EXPRESSION T_MULT P_EXPRESSION {
      $$ = binary_op($1, $3, T_MULT);
   }
 | P_EXPRESSION T_DIV P_EXPRESSION {
      $$ = binary_op($1, $3, T_DIV);
   }
 | P_EXPRESSION T_EXP P_EXPRESSION {
      $$ = binary_op($1, $3, T_EXP);
   }
 | T_MINUS P_EXPRESSION %prec T_NEGATION {
      switch ($2->d_node_type) {
         case KEY_INTEGER:
            $2->d_integer = -($2->d_integer);
            break;
         case KEY_DOUBLE:
            $2->d_double = -($2->d_double);
            break;
         case KEY_COMPLEX:
            $2->d_complex = -($2->d_complex);
            break;
         default:
            Parser::getParser()->error("X in -X is not a number");
            break;
      }
      $$ = $2;
   }
 ;

/*
 * Parse a primitive type: bool, box, char, complex, double, int, string.
 */

P_PRIMITIVE_TYPE
 : T_TRUE {
      $$ = new KeyData;
      $$->d_node_type  = KEY_BOOL;
      $$->d_array_type = KEY_BOOL;
      $$->d_array_size = 1;
      $$->d_next       = 0;
      $$->d_bool       = true;
   }
 | T_FALSE {
      $$ = new KeyData;
      $$->d_node_type  = KEY_BOOL;
      $$->d_array_type = KEY_BOOL;
      $$->d_array_size = 1;
      $$->d_next       = 0;
      $$->d_bool       = false;
   }
 | P_BOX {
      $$ = $1;
   }
 | T_CHAR {
      $$ = new KeyData;
      $$->d_node_type  = KEY_CHAR;
      $$->d_array_type = KEY_CHAR;
      $$->d_array_size = 1;
      $$->d_next       = 0;
      $$->d_char       = $1;
   }
 | P_COMPLEX {
      $$ = $1;
   }
 | T_DOUBLE {
      $$ = new KeyData;
      $$->d_node_type  = KEY_DOUBLE;
      $$->d_array_type = KEY_DOUBLE;
      $$->d_array_size = 1;
      $$->d_next       = 0;
      $$->d_double     = $1;
   }
 | T_INTEGER {
      $$ = new KeyData;
      $$->d_node_type  = KEY_INTEGER;
      $$->d_array_type = KEY_INTEGER;
      $$->d_array_size = 1;
      $$->d_next       = 0;
      $$->d_integer    = $1;
   }
 | T_STRING {
      $$ = new KeyData;
      $$->d_node_type  = KEY_STRING;
      $$->d_array_type = KEY_STRING;
      $$->d_array_size = 1;
      $$->d_next       = 0;
      $$->d_string     = *$1;
      delete $1;
   }
 ;

/*
 * Parse a complex number as (x,y), where x and y are double expressions.
 */

P_COMPLEX
 : T_OPEN_PAREN P_EXPRESSION T_COMMA P_EXPRESSION T_CLOSE_PAREN {
      to_double($2);
      to_double($4);
      $2->d_complex    = dcomplex($2->d_double, $4->d_double);
      $2->d_node_type  = KEY_COMPLEX;
      $2->d_array_type = KEY_COMPLEX;
      delete $4;
      $$ = $2;
   }
 ;

/*
 * Parse a box description of the form [(l0,l1,...,ln),(u0,u1,...,un)].
 * The two integer vectors must have the same length, and the length
 * must be between one and three, inclusive.
 */

P_BOX
 : T_OPEN_SQBKT P_INTEGER_VECTOR T_COMMA P_INTEGER_VECTOR T_CLOSE_SQBKT {
      $$ = new KeyData;
      $$->d_node_type  = KEY_BOX;
      $$->d_array_type = KEY_BOX;
      $$->d_array_size = 1;
      $$->d_next       = 0;

      if ($2->d_array_size != $4->d_array_size) {
         Parser::getParser()->error("Box lower/upper dimension mismatch");
      } else if ($2->d_array_size > SAMRAI::MAX_DIM_VAL) {
         Parser::getParser()->error("Box dimension too large (> SAMRAI::MAX_DIM_VAL)");
      } else {
         const int n = $2->d_array_size;
	 const tbox::Dimension dim(static_cast<unsigned short>(n));
         $$->d_box.setDim(dim);

         KeyData* list_lower = $2;
         KeyData* list_upper = $4;
         for (int i = n-1; i >= 0; i--) {
            $$->d_box.lower(i) = list_lower->d_integer;
            $$->d_box.upper(i) = list_upper->d_integer;
            list_lower = list_lower->d_next;
            list_upper = list_upper->d_next;
         }

         delete_list($2);
         delete_list($4);
      }
   }
 ;

/*
 * Parse an integer vector of the form (i0, i1, ..., in).
 */

P_INTEGER_VECTOR
 : T_OPEN_PAREN P_EXPRESSION_LIST T_CLOSE_PAREN {
      KeyData* list = $2;
      while (list) {
         to_integer(list);
         list = list->d_next;
      }
      $$ = $2;
   }
 ;

%%

/*
 * Delete all elements in a keyword list.
 */

static void delete_list(KeyData* list)
{
   while (list) {
      KeyData* byebye = list;
      list = list->d_next;
      delete byebye;
   }
}

/*
 * Verify that the number is a boolean; otherwise, report an error and
 * convert the argument into a boolean false.
 */

static void to_boolean(KeyData* keydata)
{
   if (keydata->d_node_type != KEY_BOOL) {
      Parser::getParser()->error("Cannot convert type into boolean");
      keydata->d_bool       = false;
      keydata->d_node_type  = KEY_BOOL;
      keydata->d_array_type = KEY_BOOL;
   }
}

/*
 * Convert the number into an integer.  If the conversion cannot be
 * performed, then print an error and return an integer zero.
 */

static void to_integer(KeyData* keydata)
{
   switch (keydata->d_node_type) {
      case KEY_INTEGER:
         break;
      case KEY_DOUBLE:
         Parser::getParser()->warning("Double truncated to integer");
         keydata->d_integer = static_cast<int>(keydata->d_double);
         break;
      case KEY_COMPLEX:
         Parser::getParser()->warning("Complex truncated to integer");
         keydata->d_integer = static_cast<int>(keydata->d_complex.real());
         break;
      default:
         Parser::getParser()->error("Cannot convert type into integer");
         keydata->d_integer = 0;
         break;
   }
   keydata->d_node_type  = KEY_INTEGER;
   keydata->d_array_type = KEY_INTEGER;
}

/*
 * Convert the number in the keydata structure to a double.  If the
 * conversion cannot be performed, then print an error and return zero.
 */

static void to_double(KeyData* keydata)
{
   switch (keydata->d_node_type) {
      case KEY_INTEGER:
         keydata->d_double = (double)(keydata->d_integer);
         break;
      case KEY_DOUBLE:
         break;
      case KEY_COMPLEX:
         Parser::getParser()->warning("Complex truncated to double");
         keydata->d_double = keydata->d_complex.real();
         break;
      default:
         Parser::getParser()->error("Cannot convert type into double");
         keydata->d_double = 0.0;
         break;
   }
   keydata->d_node_type  = KEY_DOUBLE;
   keydata->d_array_type = KEY_DOUBLE;
}

/*
 * Convert the number in the keydata structure to a complex.  If the
 * conversion cannot be performed, then print an error and return zero.
 */

static void to_complex(KeyData* keydata)
{
   switch (keydata->d_node_type) {
      case KEY_INTEGER:
         keydata->d_complex = dcomplex((double) keydata->d_integer, 0.0);
         break;
      case KEY_DOUBLE:
         keydata->d_complex = dcomplex(keydata->d_double, 0.0);
         break;
      case KEY_COMPLEX:
         break;
      default:
         Parser::getParser()->error("Cannot convert type into complex");
         keydata->d_complex = dcomplex(0.0, 0.0);
         break;
   }
   keydata->d_node_type  = KEY_COMPLEX;
   keydata->d_array_type = KEY_COMPLEX;
}

/*
 * Perform one of the standard binary operations +, -, *, /, or ^ on numeric
 * types.  Concatenation for strings is implemented above.  Return an integer
 * zero if there is a type mismatch.
 */

static KeyData* binary_op(KeyData* a, KeyData* b, const int op)
{
   if (!IS_NUMBER(a->d_node_type) || !IS_NUMBER(b->d_node_type)) {
      Parser::getParser()->error(
         "Cannot perform numerical operations on non-numeric types");
      a->d_integer    = 0;
      a->d_node_type  = KEY_INTEGER;
      a->d_array_type = KEY_INTEGER;
   } else {
      const int result_type = PROMOTE(a->d_node_type, b->d_node_type);
      switch (result_type) {
         case KEY_INTEGER:
            switch (op) {
               case T_DIV:
                  a->d_integer = a->d_integer / b->d_integer;
                  break;
               case T_EXP:
                  a->d_integer =
                     static_cast<int>(pow((double) a->d_integer, (double) b->d_integer));
                  break;
               case T_MINUS:
                  a->d_integer = a->d_integer - b->d_integer;
                  break;
               case T_MULT:
                  a->d_integer = a->d_integer * b->d_integer;
                  break;
               case T_PLUS:
                  a->d_integer = a->d_integer + b->d_integer;
                  break;
            }
            break;
         case KEY_DOUBLE:
            to_double(a);
            to_double(b);
            switch (op) {
               case T_DIV:
                  a->d_double = a->d_double / b->d_double;
                  break;
               case T_EXP:
                  a->d_double = pow(a->d_double, b->d_double);
                  break;
               case T_MINUS:
                  a->d_double = a->d_double - b->d_double;
                  break;
               case T_MULT:
                  a->d_double = a->d_double * b->d_double;
                  break;
               case T_PLUS:
                  a->d_double = a->d_double + b->d_double;
                  break;
            }
            break;
         case KEY_COMPLEX:
            to_complex(a);
            to_complex(b);
            switch (op) {
               case T_DIV:
                  a->d_complex = a->d_complex / b->d_complex;
                  break;
               case T_EXP:
                  /*
		   * SGS this is broken for insure++ and gcc 3.3.2
		   * a->d_complex = pow(a->d_complex, b->d_complex);
		   * replaced with the defn from the header file.
		   */
		  a->d_complex = exp(a->d_complex * log(b->d_complex));
                  break;
               case T_MINUS:
                  a->d_complex = a->d_complex - b->d_complex;
                  break;
               case T_MULT:
                  a->d_complex = a->d_complex * b->d_complex;
                  break;
               case T_PLUS:
                  a->d_complex = a->d_complex + b->d_complex;
                  break;
            }
            break;
      }
   }
   delete b;
   return(a);
}

/*
 * Perform one of the standard comparison operations ==, <, or >.  The other
 * operators !=, >=, and <= are computed above by using one of the first three
 * and then negating the result.  Return a boolean false if there is a type
 * mismatch problem.
 */

static KeyData* compare_op(KeyData* a, KeyData* b, const int op)
{
   if (!IS_NUMBER(a->d_node_type) || !IS_NUMBER(b->d_node_type)) {
      if (a->d_node_type != b->d_node_type) {
         Parser::getParser()->error(
            "Cannot compare different non-numeric types");
         a->d_bool = false;
      } else if (op != T_EQUALS) {
         Parser::getParser()->error(
            "Cannot apply <, >, <=, or >= to non-numeric types");
         a->d_bool = false;
      } else {
         switch(a->d_node_type) {
            case KEY_BOOL:
               a->d_bool = (a->d_bool == b->d_bool);
               break;
            case KEY_BOX:
               a->d_bool = (a->d_box == b->d_box);
               break;
            case KEY_CHAR:
               a->d_bool = (a->d_char == b->d_char);
               break;
            case KEY_STRING:
               a->d_bool = (a->d_string == b->d_string);
               break;
         }
      }
   } else {
      const int promoted = PROMOTE(a->d_node_type, b->d_node_type);
      switch (promoted) {
         case KEY_INTEGER:
            switch (op) {
               case T_EQUALS:
                  a->d_bool = (a->d_integer == b->d_integer);
                  break;
               case T_LESS:
                  a->d_bool = (a->d_integer < b->d_integer);
                  break;
               case T_GREATER:
                  a->d_bool = (a->d_integer > b->d_integer);
                  break;
            }
            break;
         case KEY_DOUBLE:
// Intel warns about comparison of floating point numbers
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
            to_double(a);
            to_double(b);
            switch (op) {
               case T_EQUALS:
                  a->d_bool = (a->d_double == b->d_double);
                  break;
               case T_LESS:
                  a->d_bool = (a->d_double < b->d_double);
                  break;
               case T_GREATER:
                  a->d_bool = (a->d_double > b->d_double);
                  break;
            }
            break;
         case KEY_COMPLEX:
            to_complex(a);
            to_complex(b);
            switch (op) {
               case T_EQUALS:
                  a->d_bool = (a->d_complex == b->d_complex);
                  break;
               case T_LESS:
               case T_GREATER:
                  Parser::getParser()->error(
                     "Operators <, >, <=, and >= are not defined for complex");
                  a->d_bool = false;
                  break;
            }
            break;
      }
   }
   a->d_node_type  = KEY_BOOL;
   a->d_array_type = KEY_BOOL;
   delete b;
   return(a);
}

/*
 * Perform a function evaluation on the specified argument.
 */

struct arith_functions {
   std::string     d_name;
   double   (*d_r2r_func)(double);
   dcomplex (*d_c2c_func)(const dcomplex&);
   double   (*d_c2r_func)(const dcomplex&);
};



#if 0
// Static initialization was not working with SGI
// compiler; so use an initialization function to 
// create the table
static arith_functions af[] = {
   { "abs"  , fabs , 0   , abs  },
   { "acos" , acos , 0   , 0    },
   { "asin" , asin , 0   , 0    },
   { "atan" , atan , 0   , 0    },
   { "ceil" , ceil , 0   , 0    },
   { "conj" , 0    , conj, 0    },
   { "cos"  , cos  , cos , 0    },
   { "cosh" , cosh , cosh, 0    },
   { "exp"  , exp  , exp , 0    },
   { "fabs" , fabs , 0   , 0    },
   { "floor", floor, 0   , 0    },
   { "imag" , 0    , 0   , imag },
   { "log10", log10, 0   , 0    },
   { "log"  , log  , log , 0    },
   { "real" , 0    , 0   , real },
   { "sin"  , sin  , sin , 0    },
   { "sinh" , sinh , sinh, 0    },
   { "sqrt" , sqrt , sqrt, 0    },
   { "tan"  , tan  , 0   , 0    },
   { ""     , 0    , 0   , 0    }
};
#endif

static arith_functions af[20];

// These are needed to deal with imag/real returning a reference
// under GCC 3.4.x
static double imag_thunk(const dcomplex &a)
{
   return imag(a);
}

static double real_thunk(const dcomplex &a)
{
   return real(a);
}

void parser_static_table_initialize()
{
   af[0].d_name =    "abs";
   af[0].d_r2r_func = fabs;
   af[0].d_c2c_func = 0;
   af[0].d_c2r_func = std::abs;


   af[1].d_name =    "acos";
   af[1].d_r2r_func = acos;
   af[1].d_c2c_func = 0;
   af[1].d_c2r_func = 0;
   
   af[2].d_name =    "asin";
   af[2].d_r2r_func = asin;
   af[2].d_c2c_func = 0;
   af[2].d_c2r_func = 0;
   
   af[3].d_name =    "atan";
   af[3].d_r2r_func = atan;
   af[3].d_c2c_func = 0;
   af[3].d_c2r_func = 0;
   
   af[4].d_name =    "ceil";
   af[4].d_r2r_func = ceil;
   af[4].d_c2c_func = 0;
   af[4].d_c2r_func = 0;

   af[5].d_name =    "conj";
   af[5].d_r2r_func = 0;
   af[5].d_c2c_func = conj;
   af[5].d_c2r_func = 0;


   af[6].d_name =    "cos";
   af[6].d_r2r_func = ::cos;
   af[6].d_c2c_func = std::cos;
   af[6].d_c2r_func = 0;

   af[7].d_name =    "cosh";
   af[7].d_r2r_func = ::cosh;
   af[7].d_c2c_func = std::cosh;
   af[7].d_c2r_func = 0;

   af[8].d_name =    "exp";
   af[8].d_r2r_func = ::exp;
   af[8].d_c2c_func = std::exp;
   af[8].d_c2r_func = 0;

   af[9].d_name =    "fabs";
   af[9].d_r2r_func = fabs;
   af[9].d_c2c_func = 0;
   af[9].d_c2r_func = 0;

   af[10].d_name =    "floor";
   af[10].d_r2r_func = floor;
   af[10].d_c2c_func = 0;
   af[10].d_c2r_func = 0;

   af[11].d_name =    "imag";
   af[11].d_r2r_func = 0;
   af[11].d_c2c_func = 0;
   af[11].d_c2r_func = imag_thunk;

   af[12].d_name =    "log10";
   af[12].d_r2r_func = ::log10;
   af[12].d_c2c_func = 0;
   af[12].d_c2r_func = 0;

   af[13].d_name =    "log";
   af[13].d_r2r_func = ::log;
   af[13].d_c2c_func = std::log;
   af[13].d_c2r_func = 0;

   af[14].d_name =    "real";
   af[14].d_r2r_func = 0;
   af[14].d_c2c_func = 0;
   af[14].d_c2r_func = real_thunk;

   af[15].d_name =    "sin";
   af[15].d_r2r_func = ::sin;
   af[15].d_c2c_func = std::sin;
   af[15].d_c2r_func = 0;

   af[16].d_name =    "sinh";
   af[16].d_r2r_func = ::sinh;
   af[16].d_c2c_func = std::sinh;
   af[16].d_c2r_func = 0;

   af[17].d_name =    "sqrt";
   af[17].d_r2r_func = ::sqrt;
   af[17].d_c2c_func = std::sqrt;
   af[17].d_c2r_func = 0;

   af[18].d_name =    "tan";
   af[18].d_r2r_func = tan;
   af[18].d_c2c_func = 0;
   af[18].d_c2r_func = 0;

   af[19].d_name =    "";
   af[19].d_r2r_func = 0;
   af[19].d_c2c_func = 0;
   af[19].d_c2r_func = 0;
}

static KeyData* eval_function(KeyData* arg, const std::string& func)
{
   if (!IS_NUMBER(arg->d_node_type)) {
      std::string tmp("Unknown function ");
      tmp += func;
      tmp += "(";
      tmp += type_names[arg->d_node_type];
      tmp += ")";
      Parser::getParser()->error(tmp);
   } else if (func == "int") {
      to_double(arg);
      arg->d_integer    = static_cast<int>(arg->d_double);
      arg->d_node_type  = KEY_INTEGER;
      arg->d_array_type = KEY_INTEGER;
   } else {
      for (int f = 0; af[f].d_name.length() > 0; f++) {
         if (af[f].d_name == func) {
            if (arg->d_node_type == KEY_COMPLEX) {
               if (af[f].d_c2c_func) {
                  arg->d_complex = (*af[f].d_c2c_func)(arg->d_complex);
               } else if (af[f].d_c2r_func) {
                  arg->d_double     = (*af[f].d_c2r_func)(arg->d_complex);
                  arg->d_node_type  = KEY_DOUBLE;
                  arg->d_array_type = KEY_DOUBLE;
               } else {
                  to_double(arg);
                  arg->d_double = (*af[f].d_r2r_func)(arg->d_double);
               }
            } else {
               if (af[f].d_r2r_func) {
                  to_double(arg);
                  arg->d_double = (*af[f].d_r2r_func)(arg->d_double);
               } else if (af[f].d_c2r_func) {
                  to_complex(arg);
                  arg->d_double     = (*af[f].d_c2r_func)(arg->d_complex);
                  arg->d_node_type  = KEY_DOUBLE;
                  arg->d_array_type = KEY_DOUBLE;
               } else {
                  to_complex(arg);
                  arg->d_complex = (*af[f].d_c2c_func)(arg->d_complex);
               }
            }
            return(arg);
         }
      }

      std::string tmp("Unknown function ");
      tmp += func;
      tmp += "(";
      tmp += type_names[arg->d_node_type];
      tmp += ")";
      Parser::getParser()->error(tmp);
   }
   return(arg);
}

/*
 * Fetch a variable in the database.  If there is an error, then print
 * an error message and return an integer zero as result.
 */

static KeyData* lookup_variable(
   const std::string& key, const int index, const bool is_array)
{
   KeyData* result = new KeyData;
   result->d_node_type  = KEY_INTEGER;
   result->d_array_type = KEY_INTEGER;
   result->d_array_size = 1;
   result->d_next       = 0;
   result->d_integer    = 0;

   Parser *parser = Parser::getParser();
   std::shared_ptr<Database> db(parser->getDatabaseWithKey(key));

   if (!db) {
      std::string tmp("Variable ``");
      tmp += key;
      tmp += "'' not found in database";
      parser->error(tmp);
   } else if (!is_array && (db->getArraySize(key) > 1)) {
      std::string tmp("Variable ``");
      tmp += key;
      tmp += "'' is not a scalar value";
      parser->error(tmp);
   } else if ((index < 0) || (index >= db->getArraySize(key))) {
      ostrstream oss;
      oss << index;
      std::string tmp("Variable ``");
      tmp += key;
      tmp += "[";
      tmp += oss.str();
      tmp += "]'' out of range";
      parser->error(tmp);
   } else if (db->isInteger(key)) {
      result->d_integer    = db->getIntegerVector(key)[index];
      result->d_node_type  = KEY_INTEGER;
      result->d_array_type = KEY_INTEGER;

   } else if (db->isDouble(key)) {
      result->d_double     = db->getDoubleVector(key)[index];
      result->d_node_type  = KEY_DOUBLE;
      result->d_array_type = KEY_DOUBLE;

   } else if (db->isComplex(key)) {
      result->d_complex    = db->getComplexVector(key)[index];
      result->d_node_type  = KEY_COMPLEX;
      result->d_array_type = KEY_COMPLEX;

   } else if (db->isBool(key)) {
      result->d_bool       = db->getBoolVector(key)[index];
      result->d_node_type  = KEY_BOOL;
      result->d_array_type = KEY_BOOL;

   } else if (db->isDatabaseBox(key)) {
      result->d_box        = db->getDatabaseBoxVector(key)[index];
      result->d_node_type  = KEY_BOX;
      result->d_array_type = KEY_BOX;

   } else if (db->isChar(key)) {
      result->d_char       = db->getCharVector(key)[index];
      result->d_node_type  = KEY_CHAR;
      result->d_array_type = KEY_CHAR;

   } else if (db->isString(key)) {
      result->d_string     = db->getStringVector(key)[index];
      result->d_node_type  = KEY_STRING;
      result->d_array_type = KEY_STRING;

   } else {
      parser->error("Unknown type for variable="+key);
   }

   return(result);
}

#ifdef __xlC__
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
