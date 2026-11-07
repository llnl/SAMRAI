#!/usr/bin/perl
#########################################################################
##
## This file is part of the SAMRAI distribution.  For full copyright 
## information, see COPYRIGHT and LICENSE. 
##
## Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
## Description:   Script in input database package. 
##
#########################################################################

while(<>) {
    s/^#line.*//;
    s/.*Revision:.*//;
    s/.*Date:.*//;
    s/.*Header:.*//;

    s/#include <unistd.h>/#ifdef SAMRAI_HAVE_UNISTD_H\n#include <unistd.h>\n#endif/;

    # substitution to replace [yylval] with SAMRAI_[yylval]
    s/yylval/SAMRAI_yylval/g;

    s/YY_DO_BEFORE_ACTION;/YY_DO_BEFORE_ACTION/;
    s/^(\s)+;$/$1do {} while(0);/;

    s/#if YY_STACK_USED/#ifdef YY_STACK_USED/;
    s/#if YY_ALWAYS_INTERACTIVE/#ifdef YY_ALWAYS_INTERACTIVE/;
    s/#if YY_MAIN/#ifdef YY_MAIN/;

    print $_;
}
