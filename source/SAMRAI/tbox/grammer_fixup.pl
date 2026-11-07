#!/usr/bin/perl

#########################################################################
##
## This file is part of the SAMRAI distribution.  For full copyright 
## information, see COPYRIGHT and LICENSE. 
##
## Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
## Description:   $Description 
##
#########################################################################

## File:        $URL$
## Package:     SAMRAI tests
## Copyright:   (c) 1997-2026 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision$
## Description: Code-generating script in inputdb package.

while(<>) {
    s/^#line.*//;
    s/.*Revision:.*//;
    s/.*Date:.*//;
    
    # substitution to replace [yynerrs,yychar,yylval] with SAMRAI_[yynerrs,yychar,yylval]
    s/yynerrs/SAMRAI_yynerrs/g;
    s/yychar/SAMRAI_yychar/g;
    s/yylval/SAMRAI_yylval/g;
    
    # These fixup some warning messages coming from insure++

    s/^# define YYDPRINTF\(Args\)$/# define YYDPRINTF(Args) do {} while (0)/;
    s/^# define YYDSYMPRINT\(Args\)$/# define YYDSYMPRINT(Args) do {} while (0)/;
    s/^# define YYDSYMPRINTF\(Title, Token, Value, Location\)$/# define YYDSYMPRINTF(Title, Token, Value, Location) do {} while (0)/;
    s/^# define YY_STACK_PRINT\(Bottom, Top\)$/# define YY_STACK_PRINT(Bottom, Top) do {} while (0)/;
    s/^# define YY_REDUCE_PRINT\(Rule\)$/# define YY_REDUCE_PRINT(Rule) do {} while (0)/;

    s/(\s+);}/$1\}/;

    # a more silent null use
    s/\(void\) yyvaluep;/if(0) {char *temp = (char *)&yyvaluep; temp++;}/;

    s/\/\* empty \*\/;//;

    print $_;
}

