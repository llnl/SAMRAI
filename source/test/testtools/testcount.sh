#!/bin/sh
#########################################################################
##
## This file is part of the SAMRAI distribution.  For full copyright 
## information, see COPYRIGHT and LICENSE. 
##
## Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
## Description:   script for testing the number of tests passed 
##
#########################################################################

num_tests=`echo $1 |  sed -e 's/,/ /g' | wc -w`
num_tests=`expr $num_tests \* $2`
num_tests=`expr $num_tests \+ $3`
num_passed=`grep PASSED $4 | wc -l`
if test $num_tests -ne $num_passed; 
   then echo "FAILED:  Count of passed tests in directory " `pwd` " passed was $num_passed expected $num_tests"
fi
