#########################################################################
##
## This file is part of the SAMRAI distribution.  For full copyright 
## information, see COPYRIGHT and LICENSE. 
##
## Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
## Description:   simple shell script to generate flex and bison files 
##
#########################################################################

dir_name=`echo ${0} | sed -e 's:^\([^/]*\)$:./\1:' -e 's:/[^/]*$::'`;
cd $dir_name

#
# Use yacc since ASCI red does not support alloca() function used by bison
#

bison -d -p yy Grammar.y
perl grammer_fixup.pl Grammar.tab.c > Grammar.cpp
perl grammer_fixup.pl Grammar.tab.h > Grammar.h
rm Grammar.tab.c
rm Grammar.tab.h

#
# Scanner requires flex due to input reading method with MPI
#

flex -Pyy -otemp.$$ Scanner.l
perl scanner_fixup.pl temp.$$ > Scanner.cpp 
rm temp.$$

# Add some pragma's to ignore warnings from compilers for
# machine generated code.
cat >> temp.$$ <<-EOF 
#ifdef __GNUC__
#ifndef __INTEL_COMPILER
#if __GNUC__ > 4 || \
              (__GNUC__ == 4 && (__GNUC_MINOR__ > 2 || \
                                 (__GNUC_MINOR__ == 2 && \
                                  __GNUC_PATCHLEVEL__ > 0)))
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
EOF
cat temp.$$ Scanner.cpp > temp2.$$
mv temp2.$$ Scanner.cpp
cat temp.$$ Grammar.cpp > temp2.$$
mv temp2.$$ Grammar.cpp
rm temp.$$

