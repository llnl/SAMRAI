#!/usr/apps/bin/perl

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
## Package:     SAMRAI application
## Copyright:   (c) 1997-2026 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision$
## Description: Script in source-manipulation utility.


# This script removes the <DIM> from constructor method names

require "find.pl";
use Cwd;
use File::Basename;
use Fcntl ':mode';

&find('.');

exit;

sub wanted {
    if ($File::Find::dir =~ m/.svn/o ) {
	$File::Find::prune = true;
    } elsif ( m/.svn/ ) {
    } elsif ( -f ) {
	doit($name);
    }
}

sub doit {
    my ($filename) = shift;
    print("$filename\n") if $debug;

    # The doit is invoked in the subdir so need to strip off 
    # the directory
    $filename=basename $filename;

    $ofilename=$filename . ".tmp";

    $mode = S_IMODE((stat($filename))[2]);

    stat($filename);
    if( -x _ ) {
	$executable=1;
    } else{
	$executable=0;
    }

    open(IN, $filename) || die "Cannot open input file $filename : $!";
    open OUT, "> $ofilename" || die "Cannot open output file $ofilename";

     while (<IN>) {
	 } else {
	     print OUT;
	 }
     }

     close IN  || die "Cannot close file $filename";
     close OUT || die "Cannot close file $ofilename";

    unlink($filename);
    rename($ofilename, $filename);

    chmod $mode, $filename;

}
