#ifdef __GNUC__
#ifndef __INTEL_COMPILER
#if __GNUC__ > 4 ||               (__GNUC__ == 4 && (__GNUC_MINOR__ > 2 ||                                  (__GNUC_MINOR__ == 2 &&                                   __GNUC_PATCHLEVEL__ > 0)))
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wconversion"
#endif
#endif
#endif

#ifdef __INTEL_COMPILER
// Ignore Intel warnings about unreachable statements
#pragma warning (disable:177)
// Ignore Intel warnings about external declarations
#pragma warning (disable:1419)
// Ignore Intel warnings about type conversions
#pragma warning (disable:810)
// Ignore Intel remarks about non-pointer conversions
#pragma warning (disable:2259)
// Ignore Intel remarks about zero used for undefined preprocessor syms
#pragma warning (disable:193)
// Ignore Intel remarks about unreachable code
#pragma warning (disable:111)
#endif

#ifdef __xlC__
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif
/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     T_AND = 258,
     T_ASSIGN = 259,
     T_CHAR = 260,
     T_CLOSE_CURLY = 261,
     T_CLOSE_PAREN = 262,
     T_CLOSE_SQBKT = 263,
     T_COMMA = 264,
     T_DIV = 265,
     T_DOUBLE = 266,
     T_ELSE = 267,
     T_EXP = 268,
     T_EQUALS = 269,
     T_GREATER_EQUALS = 270,
     T_GREATER = 271,
     T_LESS_EQUALS = 272,
     T_LESS = 273,
     T_FALSE = 274,
     T_INTEGER = 275,
     T_KEYWORD = 276,
     T_MINUS = 277,
     T_MULT = 278,
     T_NOT = 279,
     T_NOT_EQUALS = 280,
     T_OR = 281,
     T_OPEN_CURLY = 282,
     T_OPEN_PAREN = 283,
     T_OPEN_SQBKT = 284,
     T_PLUS = 285,
     T_QUESTION = 286,
     T_SEMI = 287,
     T_STRING = 288,
     T_TRUE = 289,
     T_NEGATION = 290
   };
#endif
/* Tokens.  */
#define T_AND 258
#define T_ASSIGN 259
#define T_CHAR 260
#define T_CLOSE_CURLY 261
#define T_CLOSE_PAREN 262
#define T_CLOSE_SQBKT 263
#define T_COMMA 264
#define T_DIV 265
#define T_DOUBLE 266
#define T_ELSE 267
#define T_EXP 268
#define T_EQUALS 269
#define T_GREATER_EQUALS 270
#define T_GREATER 271
#define T_LESS_EQUALS 272
#define T_LESS 273
#define T_FALSE 274
#define T_INTEGER 275
#define T_KEYWORD 276
#define T_MINUS 277
#define T_MULT 278
#define T_NOT 279
#define T_NOT_EQUALS 280
#define T_OR 281
#define T_OPEN_CURLY 282
#define T_OPEN_PAREN 283
#define T_OPEN_SQBKT 284
#define T_PLUS 285
#define T_QUESTION 286
#define T_SEMI 287
#define T_STRING 288
#define T_TRUE 289
#define T_NEGATION 290




/* Copy the first part of user declarations.  */


//
// This file is part of the SAMRAI distribution.  For full copyright
// information, see COPYRIGHT and LICENSE.
//
// Copyright:	(c) 1997-2026 Lawrence Livermore National Security, LLC
// Description:	Yacc grammar description for the input database
//

#include "SAMRAI/SAMRAI_config.h"
#include <math.h>

#include STL_SSTREAM_HEADER_FILE


#if !defined(OSTRINGSTREAM_TYPE_IS_BROKEN) && defined(OSTRSTREAM_TYPE_IS_BROKEN)
typedef std::ostringstream ostrstream;
#endif

#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Parser.h"
#include <string>

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
static KeyData* lookup_variable(const std::string&, const size_t, const bool);



/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE

{
  char          u_char;
  double        u_double;
  int           u_integer;
  KeyData* u_keydata;
  std::string*  u_keyword;
  std::string*  u_string;
}
/* Line 193 of yacc.c.  */

	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */


#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do {  } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  3
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   249

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  36
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  13
/* YYNRULES -- Number of rules.  */
#define YYNRULES  44
/* YYNRULES -- Number of states.  */
#define YYNSTATES  83

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   290

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     5,     6,     9,    10,    16,    17,    22,
      24,    26,    30,    32,    34,    39,    43,    48,    49,    58,
      61,    65,    69,    73,    77,    81,    85,    89,    93,    97,
     101,   105,   109,   113,   116,   118,   120,   122,   124,   126,
     128,   130,   132,   138,   144
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      37,     0,    -1,    38,    -1,    -1,    38,    39,    -1,    -1,
      21,    27,    40,    38,     6,    -1,    -1,    21,     4,    41,
      42,    -1,    32,    -1,    43,    -1,    42,     9,    43,    -1,
      45,    -1,    21,    -1,    21,    29,    43,     8,    -1,    28,
      43,     7,    -1,    21,    28,    43,     7,    -1,    -1,    28,
      43,    31,    44,    43,    12,    43,     7,    -1,    24,    43,
      -1,    43,    26,    43,    -1,    43,     3,    43,    -1,    43,
      14,    43,    -1,    43,    25,    43,    -1,    43,    15,    43,
      -1,    43,    16,    43,    -1,    43,    17,    43,    -1,    43,
      18,    43,    -1,    43,    30,    43,    -1,    43,    22,    43,
      -1,    43,    23,    43,    -1,    43,    10,    43,    -1,    43,
      13,    43,    -1,    22,    43,    -1,    34,    -1,    19,    -1,
      47,    -1,     5,    -1,    46,    -1,    11,    -1,    20,    -1,
      33,    -1,    28,    43,     9,    43,     7,    -1,    29,    48,
       9,    48,     8,    -1,    28,    42,     7,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   158,   158,   165,   167,   177,   177,   196,   196,   281,
     297,   300,   376,   379,   383,   389,   392,   396,   396,   410,
     415,   422,   429,   432,   436,   440,   443,   447,   450,   461,
     464,   467,   470,   473,   497,   505,   513,   516,   524,   527,
     535,   543,   559,   577,   613
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "T_AND", "T_ASSIGN", "T_CHAR",
  "T_CLOSE_CURLY", "T_CLOSE_PAREN", "T_CLOSE_SQBKT", "T_COMMA", "T_DIV",
  "T_DOUBLE", "T_ELSE", "T_EXP", "T_EQUALS", "T_GREATER_EQUALS",
  "T_GREATER", "T_LESS_EQUALS", "T_LESS", "T_FALSE", "T_INTEGER",
  "T_KEYWORD", "T_MINUS", "T_MULT", "T_NOT", "T_NOT_EQUALS", "T_OR",
  "T_OPEN_CURLY", "T_OPEN_PAREN", "T_OPEN_SQBKT", "T_PLUS", "T_QUESTION",
  "T_SEMI", "T_STRING", "T_TRUE", "T_NEGATION", "$accept",
  "P_SPECIFICATION", "P_DEFINITION_LIST", "P_DEFINITION", "@1", "@2",
  "P_EXPRESSION_LIST", "P_EXPRESSION", "@3", "P_PRIMITIVE_TYPE",
  "P_COMPLEX", "P_BOX", "P_INTEGER_VECTOR", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    36,    37,    38,    38,    40,    39,    41,    39,    39,
      42,    42,    43,    43,    43,    43,    43,    44,    43,    43,
      43,    43,    43,    43,    43,    43,    43,    43,    43,    43,
      43,    43,    43,    43,    45,    45,    45,    45,    45,    45,
      45,    45,    46,    47,    48
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     0,     2,     0,     5,     0,     4,     1,
       1,     3,     1,     1,     4,     3,     4,     0,     8,     2,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     5,     5,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       3,     0,     2,     1,     0,     9,     4,     7,     5,     0,
       3,    37,    39,    35,    40,    13,     0,     0,     0,     0,
      41,    34,     8,    10,    12,    38,    36,     0,     0,     0,
      33,    19,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     6,
       0,     0,    15,     0,    17,     0,     0,    11,    21,    31,
      32,    22,    24,    25,    26,    27,    29,    30,    23,    20,
      28,    16,    14,     0,     0,    44,     0,    42,     0,    43,
       0,     0,    18
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     1,     2,     6,    10,     9,    22,    23,    74,    24,
      25,    26,    34
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -21
static const yytype_int16 yypact[] =
{
     -21,    10,   -17,   -21,     7,   -21,   -21,   -21,   -21,    28,
     -21,   -21,   -21,   -21,   -21,   -20,    28,    28,    28,   -14,
     -21,   -21,     8,   183,   -21,   -21,   -21,    69,    28,    28,
     -21,   219,    56,    28,     9,    28,    28,    28,    28,    28,
      28,    28,    28,    28,    28,    28,    28,    28,    28,   -21,
      81,   102,   -21,    28,   -21,    -2,   -14,   183,   219,    22,
      22,    70,    70,    70,    70,    70,    -7,    22,    70,   201,
      -7,   -21,   -21,   123,    28,   -21,    30,   -21,   144,   -21,
      28,   165,   -21
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -21,   -21,    26,   -21,   -21,   -21,    11,   -16,   -21,   -21,
     -21,   -21,   -15
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      30,    31,    32,    37,     4,    75,    38,    35,    28,    29,
       3,     7,    50,    51,    33,     5,    45,    35,    56,    57,
      58,    59,    60,    61,    62,    63,    64,    65,    66,    67,
      68,    69,    70,    11,     8,    38,    27,    73,    79,    12,
       0,    76,     0,     0,    55,     0,     0,    13,    14,    15,
      16,     0,    17,     0,     0,     0,    18,    19,    78,    36,
       0,    20,    21,    52,    81,    53,    37,     0,     0,    38,
      39,    40,    41,    42,    43,    49,     0,     0,    44,    45,
      37,    46,    47,    38,    36,     0,    48,    54,    71,     0,
       4,    37,    44,    45,    38,    39,    40,    41,    42,    43,
      48,     5,     0,    44,    45,    36,    46,    47,     0,     0,
      72,    48,    37,     0,     0,    38,    39,    40,    41,    42,
      43,     0,     0,     0,    44,    45,    36,    46,    47,     0,
      77,     0,    48,    37,     0,     0,    38,    39,    40,    41,
      42,    43,     0,     0,     0,    44,    45,    36,    46,    47,
       0,     0,     0,    48,    37,     0,    80,    38,    39,    40,
      41,    42,    43,     0,     0,     0,    44,    45,    36,    46,
      47,     0,    82,     0,    48,    37,     0,     0,    38,    39,
      40,    41,    42,    43,     0,     0,    36,    44,    45,     0,
      46,    47,     0,    37,     0,    48,    38,    39,    40,    41,
      42,    43,     0,     0,    36,    44,    45,     0,    46,    47,
       0,    37,     0,    48,    38,    39,    40,    41,    42,    43,
       0,     0,     0,    44,    45,     0,    46,     0,     0,    37,
       0,    48,    38,    39,    40,    41,    42,    43,     0,     0,
       0,    44,    45,     0,    46,     0,     0,     0,     0,    48
};

static const yytype_int8 yycheck[] =
{
      16,    17,    18,    10,    21,     7,    13,     9,    28,    29,
       0,     4,    28,    29,    28,    32,    23,     9,     9,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,     5,    27,    13,    10,    53,     8,    11,
      -1,    56,    -1,    -1,    33,    -1,    -1,    19,    20,    21,
      22,    -1,    24,    -1,    -1,    -1,    28,    29,    74,     3,
      -1,    33,    34,     7,    80,     9,    10,    -1,    -1,    13,
      14,    15,    16,    17,    18,     6,    -1,    -1,    22,    23,
      10,    25,    26,    13,     3,    -1,    30,    31,     7,    -1,
      21,    10,    22,    23,    13,    14,    15,    16,    17,    18,
      30,    32,    -1,    22,    23,     3,    25,    26,    -1,    -1,
       8,    30,    10,    -1,    -1,    13,    14,    15,    16,    17,
      18,    -1,    -1,    -1,    22,    23,     3,    25,    26,    -1,
       7,    -1,    30,    10,    -1,    -1,    13,    14,    15,    16,
      17,    18,    -1,    -1,    -1,    22,    23,     3,    25,    26,
      -1,    -1,    -1,    30,    10,    -1,    12,    13,    14,    15,
      16,    17,    18,    -1,    -1,    -1,    22,    23,     3,    25,
      26,    -1,     7,    -1,    30,    10,    -1,    -1,    13,    14,
      15,    16,    17,    18,    -1,    -1,     3,    22,    23,    -1,
      25,    26,    -1,    10,    -1,    30,    13,    14,    15,    16,
      17,    18,    -1,    -1,     3,    22,    23,    -1,    25,    26,
      -1,    10,    -1,    30,    13,    14,    15,    16,    17,    18,
      -1,    -1,    -1,    22,    23,    -1,    25,    -1,    -1,    10,
      -1,    30,    13,    14,    15,    16,    17,    18,    -1,    -1,
      -1,    22,    23,    -1,    25,    -1,    -1,    -1,    -1,    30
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    37,    38,     0,    21,    32,    39,     4,    27,    41,
      40,     5,    11,    19,    20,    21,    22,    24,    28,    29,
      33,    34,    42,    43,    45,    46,    47,    38,    28,    29,
      43,    43,    43,    28,    48,     9,     3,    10,    13,    14,
      15,    16,    17,    18,    22,    23,    25,    26,    30,     6,
      43,    43,     7,     9,    31,    42,     9,    43,    43,    43,
      43,    43,    43,    43,    43,    43,    43,    43,    43,    43,
      43,     7,     8,    43,    44,     7,    48,     7,    43,     8,
      12,    43,     7
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(SAMRAI_yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (SAMRAI_yychar == YYEMPTY && yylen == 1)				\
    {								\
      SAMRAI_yychar = (Token);						\
      SAMRAI_yylval = (Value);						\
      yytoken = YYTRANSLATE (SAMRAI_yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args) do {} while (0)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top) do {} while (0)
# define YY_REDUCE_PRINT(Rule) do {} while (0)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int SAMRAI_yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (SAMRAI_yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int SAMRAI_yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE SAMRAI_yylval;

/* Number of syntax errors so far.  */
int SAMRAI_yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  SAMRAI_yynerrs = 0;
  SAMRAI_yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (SAMRAI_yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      SAMRAI_yychar = YYLEX;
    }

  if (SAMRAI_yychar <= YYEOF)
    {
      SAMRAI_yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (SAMRAI_yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &SAMRAI_yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &SAMRAI_yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (SAMRAI_yychar != YYEOF)
    SAMRAI_yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = SAMRAI_yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 5:

    {

   /* This is a hack to make a warning message go away from
      a symbol flex defines but does not use */
   if(0) {
      goto yyerrlab1;
   }

      if (Parser::getParser()->getScope()->keyExists(*(yyvsp[(1) - (2)].u_keyword))) {
	 std::string tmp("Redefinition of key ``");
         tmp += *(yyvsp[(1) - (2)].u_keyword);
         tmp += "''";
         Parser::getParser()->warning(tmp);
      }
      Parser::getParser()->enterScope(*(yyvsp[(1) - (2)].u_keyword));
   }
    break;

  case 6:

    {
      Parser::getParser()->leaveScope();
      delete (yyvsp[(1) - (5)].u_keyword);
   }
    break;

  case 7:

    {
      if (Parser::getParser()->getScope()->keyExists(*(yyvsp[(1) - (2)].u_keyword))) {
	 std::string tmp("Redefinition of key ``");
         tmp += *(yyvsp[(1) - (2)].u_keyword);
         tmp += "''";
         Parser::getParser()->warning(tmp);
      }
   }
    break;

  case 8:

    {
      KeyData* list = (yyvsp[(4) - (4)].u_keydata);
      const int n = list->d_array_size;

      switch (list->d_array_type) {
         case KEY_BOOL: {
            std::vector<bool> data(n);
            for (int i = n-1; i >= 0; i--) {
               data[i] = list->d_bool;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putBoolVector(*(yyvsp[(1) - (4)].u_keyword), data);
            break;
         }
         case KEY_BOX: {
            std::vector<DatabaseBox> data(n);
            for (int i = n-1; i >= 0; i--) {
               data[i] = list->d_box;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putDatabaseBoxVector(*(yyvsp[(1) - (4)].u_keyword), data);
            break;
         }
         case KEY_CHAR: {
            std::vector<char> data(n);
            for (int i = n-1; i >= 0; i--) {
               data[i] = list->d_char;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putCharVector(*(yyvsp[(1) - (4)].u_keyword), data);
            break;
         }
         case KEY_COMPLEX: {
            std::vector<dcomplex> data(n);
            for (int i = n-1; i >= 0; i--) {
               to_complex(list);
               data[i] = list->d_complex;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putComplexVector(*(yyvsp[(1) - (4)].u_keyword), data);
            break;
         }
         case KEY_DOUBLE: {
            std::vector<double> data(n);
            for (int i = n-1; i >= 0; i--) {
               to_double(list);
               data[i] = list->d_double;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putDoubleVector(*(yyvsp[(1) - (4)].u_keyword), data);
            break;
         }
         case KEY_INTEGER: {
            std::vector<int> data(n);
            for (int i = n-1; i >= 0; i--) {
               data[i] = list->d_integer;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putIntegerVector(*(yyvsp[(1) - (4)].u_keyword), data);
            break;
         }
         case KEY_STRING: {
            std::vector<std::string> data(n);
            for (int i = n-1; i >= 0; i--) {
               data[i] = list->d_string;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putStringVector(*(yyvsp[(1) - (4)].u_keyword), data);
            break;
         }
         default:
            Parser::getParser()->error("Internal parser error!");
            break;
      }

      delete_list((yyvsp[(4) - (4)].u_keydata));
      delete (yyvsp[(1) - (4)].u_keyword);
   }
    break;

  case 9:

    {
      Parser::getParser()->warning(
         "Semicolon found in keyword phrase (ignored)");
   }
    break;

  case 10:

    {
      (yyval.u_keydata) = (yyvsp[(1) - (1)].u_keydata);
   }
    break;

  case 11:

    {
      switch((yyvsp[(1) - (3)].u_keydata)->d_array_type) {
         case KEY_BOOL:
         case KEY_CHAR:
         case KEY_STRING:
            if ((yyvsp[(3) - (3)].u_keydata)->d_node_type != (yyvsp[(1) - (3)].u_keydata)->d_array_type) {
               Parser::getParser()->error("Type mismatch in array");
               delete (yyvsp[(3) - (3)].u_keydata);
               (yyval.u_keydata) = (yyvsp[(1) - (3)].u_keydata);
            } else {
               (yyvsp[(3) - (3)].u_keydata)->d_array_size = (yyvsp[(1) - (3)].u_keydata)->d_array_size + 1;
               (yyvsp[(3) - (3)].u_keydata)->d_next       = (yyvsp[(1) - (3)].u_keydata);
               (yyval.u_keydata)               = (yyvsp[(3) - (3)].u_keydata);
            }
            break;
         case KEY_BOX:
            if ((yyvsp[(3) - (3)].u_keydata)->d_node_type != KEY_BOX) {
               Parser::getParser()->error("Type mismatch in box array");
               delete (yyvsp[(3) - (3)].u_keydata);
               (yyval.u_keydata) = (yyvsp[(1) - (3)].u_keydata);
            } else if ((yyvsp[(3) - (3)].u_keydata)->d_box.getDimVal() != (yyvsp[(1) - (3)].u_keydata)->d_box.getDimVal()) {
               Parser::getParser()->error("Box array dimension mismatch");
               delete (yyvsp[(3) - (3)].u_keydata);
               (yyval.u_keydata) = (yyvsp[(1) - (3)].u_keydata);
            } else {
               (yyvsp[(3) - (3)].u_keydata)->d_array_size = (yyvsp[(1) - (3)].u_keydata)->d_array_size + 1;
               (yyvsp[(3) - (3)].u_keydata)->d_next       = (yyvsp[(1) - (3)].u_keydata);
               (yyval.u_keydata)               = (yyvsp[(3) - (3)].u_keydata);
            }
            break;
         case KEY_COMPLEX:
         case KEY_DOUBLE:
         case KEY_INTEGER:
            if (!IS_NUMBER((yyvsp[(3) - (3)].u_keydata)->d_node_type)) {
               Parser::getParser()->error("Type mismatch in number array");
               delete (yyvsp[(3) - (3)].u_keydata);
               (yyval.u_keydata) = (yyvsp[(1) - (3)].u_keydata);
            } else {
               (yyvsp[(3) - (3)].u_keydata)->d_array_type = PROMOTE((yyvsp[(1) - (3)].u_keydata)->d_array_type, (yyvsp[(3) - (3)].u_keydata)->d_node_type);
               (yyvsp[(3) - (3)].u_keydata)->d_array_size = (yyvsp[(1) - (3)].u_keydata)->d_array_size + 1;
               (yyvsp[(3) - (3)].u_keydata)->d_next       = (yyvsp[(1) - (3)].u_keydata);
               (yyval.u_keydata)               = (yyvsp[(3) - (3)].u_keydata);
            }
            break;
      }
   }
    break;

  case 12:

    {
      (yyval.u_keydata) = (yyvsp[(1) - (1)].u_keydata);
   }
    break;

  case 13:

    {
      (yyval.u_keydata) = lookup_variable(*(yyvsp[(1) - (1)].u_keyword), 0, false);
      delete (yyvsp[(1) - (1)].u_keyword);
   }
    break;

  case 14:

    {
      to_integer((yyvsp[(3) - (4)].u_keydata));
      (yyval.u_keydata) = lookup_variable(*(yyvsp[(1) - (4)].u_keyword), (yyvsp[(3) - (4)].u_keydata)->d_integer, true);
      delete (yyvsp[(1) - (4)].u_keyword);
      delete (yyvsp[(3) - (4)].u_keydata);
   }
    break;

  case 15:

    {
      (yyval.u_keydata) = (yyvsp[(2) - (3)].u_keydata);
   }
    break;

  case 16:

    {
      (yyval.u_keydata) = eval_function((yyvsp[(3) - (4)].u_keydata), *(yyvsp[(1) - (4)].u_keyword));
      delete (yyvsp[(1) - (4)].u_keyword);
   }
    break;

  case 17:

    {
      if ((yyvsp[(2) - (3)].u_keydata)->d_node_type != KEY_BOOL) {
         Parser::getParser()->error("X in (X ? Y : Z) is not a boolean");
      }
   }
    break;

  case 18:

    {
      if (((yyvsp[(2) - (8)].u_keydata)->d_node_type == KEY_BOOL) && ((yyvsp[(2) - (8)].u_keydata)->d_bool)) {
         (yyval.u_keydata) = (yyvsp[(5) - (8)].u_keydata);
         delete (yyvsp[(7) - (8)].u_keydata);
      } else {
         (yyval.u_keydata) = (yyvsp[(7) - (8)].u_keydata);
         delete (yyvsp[(5) - (8)].u_keydata);
      }
      delete (yyvsp[(2) - (8)].u_keydata);
   }
    break;

  case 19:

    {
      to_boolean((yyvsp[(2) - (2)].u_keydata));
      (yyvsp[(2) - (2)].u_keydata)->d_bool = !(yyvsp[(2) - (2)].u_keydata)->d_bool;
      (yyval.u_keydata) = (yyvsp[(2) - (2)].u_keydata);
   }
    break;

  case 20:

    {
      to_boolean((yyvsp[(1) - (3)].u_keydata));
      to_boolean((yyvsp[(3) - (3)].u_keydata));
      (yyvsp[(1) - (3)].u_keydata)->d_bool = (yyvsp[(1) - (3)].u_keydata)->d_bool || (yyvsp[(3) - (3)].u_keydata)->d_bool;
      delete (yyvsp[(3) - (3)].u_keydata);
      (yyval.u_keydata) = (yyvsp[(1) - (3)].u_keydata);
   }
    break;

  case 21:

    {
      to_boolean((yyvsp[(1) - (3)].u_keydata));
      to_boolean((yyvsp[(3) - (3)].u_keydata));
      (yyvsp[(1) - (3)].u_keydata)->d_bool = (yyvsp[(1) - (3)].u_keydata)->d_bool && (yyvsp[(3) - (3)].u_keydata)->d_bool;
      delete (yyvsp[(3) - (3)].u_keydata);
      (yyval.u_keydata) = (yyvsp[(1) - (3)].u_keydata);
   }
    break;

  case 22:

    {
      (yyval.u_keydata) = compare_op((yyvsp[(1) - (3)].u_keydata), (yyvsp[(3) - (3)].u_keydata), T_EQUALS);
   }
    break;

  case 23:

    {
      (yyval.u_keydata) = compare_op((yyvsp[(1) - (3)].u_keydata), (yyvsp[(3) - (3)].u_keydata), T_EQUALS);
      (yyval.u_keydata)->d_bool = !((yyval.u_keydata)->d_bool);
   }
    break;

  case 24:

    {
      (yyval.u_keydata) = compare_op((yyvsp[(1) - (3)].u_keydata), (yyvsp[(3) - (3)].u_keydata), T_LESS);
      (yyval.u_keydata)->d_bool = !((yyval.u_keydata)->d_bool);
   }
    break;

  case 25:

    {
      (yyval.u_keydata) = compare_op((yyvsp[(1) - (3)].u_keydata), (yyvsp[(3) - (3)].u_keydata), T_GREATER);
   }
    break;

  case 26:

    {
      (yyval.u_keydata) = compare_op((yyvsp[(1) - (3)].u_keydata), (yyvsp[(3) - (3)].u_keydata), T_GREATER);
      (yyval.u_keydata)->d_bool = !((yyval.u_keydata)->d_bool);
   }
    break;

  case 27:

    {
      (yyval.u_keydata) = compare_op((yyvsp[(1) - (3)].u_keydata), (yyvsp[(3) - (3)].u_keydata), T_LESS);
   }
    break;

  case 28:

    {
      if (((yyvsp[(1) - (3)].u_keydata)->d_node_type == KEY_STRING) && ((yyvsp[(3) - (3)].u_keydata)->d_node_type == KEY_STRING)) {
	 std::string tmp((yyvsp[(1) - (3)].u_keydata)->d_string);
	 tmp += (yyvsp[(3) - (3)].u_keydata)->d_string;
         (yyvsp[(1) - (3)].u_keydata)->d_string = tmp;
         delete (yyvsp[(3) - (3)].u_keydata);
         (yyval.u_keydata) = (yyvsp[(1) - (3)].u_keydata);
      } else {
         (yyval.u_keydata) = binary_op((yyvsp[(1) - (3)].u_keydata), (yyvsp[(3) - (3)].u_keydata), T_PLUS);
      }
   }
    break;

  case 29:

    {
      (yyval.u_keydata) = binary_op((yyvsp[(1) - (3)].u_keydata), (yyvsp[(3) - (3)].u_keydata), T_MINUS);
   }
    break;

  case 30:

    {
      (yyval.u_keydata) = binary_op((yyvsp[(1) - (3)].u_keydata), (yyvsp[(3) - (3)].u_keydata), T_MULT);
   }
    break;

  case 31:

    {
      (yyval.u_keydata) = binary_op((yyvsp[(1) - (3)].u_keydata), (yyvsp[(3) - (3)].u_keydata), T_DIV);
   }
    break;

  case 32:

    {
      (yyval.u_keydata) = binary_op((yyvsp[(1) - (3)].u_keydata), (yyvsp[(3) - (3)].u_keydata), T_EXP);
   }
    break;

  case 33:

    {
      switch ((yyvsp[(2) - (2)].u_keydata)->d_node_type) {
         case KEY_INTEGER:
            (yyvsp[(2) - (2)].u_keydata)->d_integer = -((yyvsp[(2) - (2)].u_keydata)->d_integer);
            break;
         case KEY_DOUBLE:
            (yyvsp[(2) - (2)].u_keydata)->d_double = -((yyvsp[(2) - (2)].u_keydata)->d_double);
            break;
         case KEY_COMPLEX:
            (yyvsp[(2) - (2)].u_keydata)->d_complex = -((yyvsp[(2) - (2)].u_keydata)->d_complex);
            break;
         default:
            Parser::getParser()->error("X in -X is not a number");
            break;
      }
      (yyval.u_keydata) = (yyvsp[(2) - (2)].u_keydata);
   }
    break;

  case 34:

    {
      (yyval.u_keydata) = new KeyData;
      (yyval.u_keydata)->d_node_type  = KEY_BOOL;
      (yyval.u_keydata)->d_array_type = KEY_BOOL;
      (yyval.u_keydata)->d_array_size = 1;
      (yyval.u_keydata)->d_next       = 0;
      (yyval.u_keydata)->d_bool       = true;
   }
    break;

  case 35:

    {
      (yyval.u_keydata) = new KeyData;
      (yyval.u_keydata)->d_node_type  = KEY_BOOL;
      (yyval.u_keydata)->d_array_type = KEY_BOOL;
      (yyval.u_keydata)->d_array_size = 1;
      (yyval.u_keydata)->d_next       = 0;
      (yyval.u_keydata)->d_bool       = false;
   }
    break;

  case 36:

    {
      (yyval.u_keydata) = (yyvsp[(1) - (1)].u_keydata);
   }
    break;

  case 37:

    {
      (yyval.u_keydata) = new KeyData;
      (yyval.u_keydata)->d_node_type  = KEY_CHAR;
      (yyval.u_keydata)->d_array_type = KEY_CHAR;
      (yyval.u_keydata)->d_array_size = 1;
      (yyval.u_keydata)->d_next       = 0;
      (yyval.u_keydata)->d_char       = (yyvsp[(1) - (1)].u_char);
   }
    break;

  case 38:

    {
      (yyval.u_keydata) = (yyvsp[(1) - (1)].u_keydata);
   }
    break;

  case 39:

    {
      (yyval.u_keydata) = new KeyData;
      (yyval.u_keydata)->d_node_type  = KEY_DOUBLE;
      (yyval.u_keydata)->d_array_type = KEY_DOUBLE;
      (yyval.u_keydata)->d_array_size = 1;
      (yyval.u_keydata)->d_next       = 0;
      (yyval.u_keydata)->d_double     = (yyvsp[(1) - (1)].u_double);
   }
    break;

  case 40:

    {
      (yyval.u_keydata) = new KeyData;
      (yyval.u_keydata)->d_node_type  = KEY_INTEGER;
      (yyval.u_keydata)->d_array_type = KEY_INTEGER;
      (yyval.u_keydata)->d_array_size = 1;
      (yyval.u_keydata)->d_next       = 0;
      (yyval.u_keydata)->d_integer    = (yyvsp[(1) - (1)].u_integer);
   }
    break;

  case 41:

    {
      (yyval.u_keydata) = new KeyData;
      (yyval.u_keydata)->d_node_type  = KEY_STRING;
      (yyval.u_keydata)->d_array_type = KEY_STRING;
      (yyval.u_keydata)->d_array_size = 1;
      (yyval.u_keydata)->d_next       = 0;
      (yyval.u_keydata)->d_string     = *(yyvsp[(1) - (1)].u_string);
      delete (yyvsp[(1) - (1)].u_string);
   }
    break;

  case 42:

    {
      to_double((yyvsp[(2) - (5)].u_keydata));
      to_double((yyvsp[(4) - (5)].u_keydata));
      (yyvsp[(2) - (5)].u_keydata)->d_complex    = dcomplex((yyvsp[(2) - (5)].u_keydata)->d_double, (yyvsp[(4) - (5)].u_keydata)->d_double);
      (yyvsp[(2) - (5)].u_keydata)->d_node_type  = KEY_COMPLEX;
      (yyvsp[(2) - (5)].u_keydata)->d_array_type = KEY_COMPLEX;
      delete (yyvsp[(4) - (5)].u_keydata);
      (yyval.u_keydata) = (yyvsp[(2) - (5)].u_keydata);
   }
    break;

  case 43:

    {
      (yyval.u_keydata) = new KeyData;
      (yyval.u_keydata)->d_node_type  = KEY_BOX;
      (yyval.u_keydata)->d_array_type = KEY_BOX;
      (yyval.u_keydata)->d_array_size = 1;
      (yyval.u_keydata)->d_next       = 0;

      if ((yyvsp[(2) - (5)].u_keydata)->d_array_size != (yyvsp[(4) - (5)].u_keydata)->d_array_size) {
         Parser::getParser()->error("Box lower/upper dimension mismatch");
      } else if ((yyvsp[(2) - (5)].u_keydata)->d_array_size > SAMRAI::MAX_DIM_VAL) {
         Parser::getParser()->error("Box dimension too large (> SAMRAI::MAX_DIM_VAL)");
      } else {
         const int n = (yyvsp[(2) - (5)].u_keydata)->d_array_size;
	 const Dimension dim(static_cast<unsigned short>(n));
         (yyval.u_keydata)->d_box.setDim(dim);

         KeyData* list_lower = (yyvsp[(2) - (5)].u_keydata);
         KeyData* list_upper = (yyvsp[(4) - (5)].u_keydata);
         for (int i = n-1; i >= 0; i--) {
            (yyval.u_keydata)->d_box.lower(i) = list_lower->d_integer;
            (yyval.u_keydata)->d_box.upper(i) = list_upper->d_integer;
            list_lower = list_lower->d_next;
            list_upper = list_upper->d_next;
         }

         delete_list((yyvsp[(2) - (5)].u_keydata));
         delete_list((yyvsp[(4) - (5)].u_keydata));
      }
   }
    break;

  case 44:

    {
      KeyData* list = (yyvsp[(2) - (3)].u_keydata);
      while (list) {
         to_integer(list);
         list = list->d_next;
      }
      (yyval.u_keydata) = (yyvsp[(2) - (3)].u_keydata);
   }
    break;


/* Line 1267 of yacc.c.  */

      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++SAMRAI_yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, SAMRAI_yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, SAMRAI_yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (SAMRAI_yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (SAMRAI_yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &SAMRAI_yylval);
	  SAMRAI_yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = SAMRAI_yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (SAMRAI_yychar != YYEOF && SAMRAI_yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &SAMRAI_yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}





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
   af[5].d_c2c_func = std::conj;
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
   const std::string& key, const size_t index, const bool is_array)
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
   } else if (index >= db->getArraySize(key)) {
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

