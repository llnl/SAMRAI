#!/usr/bin/perl
#########################################################################
##
## This file is part of the SAMRAI distribution.  For full copyright 
## information, see COPYRIGHT and LICENSE. 
##
## Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
## Description:   perl script to indent SAMRAI 
##
#########################################################################

use strict;

use File::Basename;
use File::Find;
use File::Path;
use File::Compare;
use Cwd;
use Text::Wrap;
use File::Copy;

my $pwd = cwd;

my $debug = 1;

#
# File pattern to look for
#
my $filePattern;
my $fileExcludePattern;
my $excludeDirPattern;

my @allfiles = ();
sub selectFile {
    if ( $File::Find::dir =~ m!$excludeDirPattern! ) {
	$File::Find::prune = 1;
    }
    elsif ( -f && m/$filePattern/ && ! m/$fileExcludePattern/ ) {
	push @allfiles, $File::Find::name;
	$allfiles[$#allfiles] =~ s|^\./||;
    }
}

# Flush I/O on write to avoid buffering
$|=1;


my $end_of_line = $/;

@allfiles = ();
$filePattern = q!(.*\.(([ChI]))$)!;
$fileExcludePattern = q!(Grammar.[Ch]|Scanner.cpp)!;
$excludeDirPattern=q!/(.svn|CVS|automatic|include)$!;
find( \&selectFile, "." );
print "files=@allfiles" if ($debug > 2);
foreach my $file (@allfiles) {
    print "Working on $file\n";
    my $directory = dirname $file;

    my $filebasename = basename $file;

    system "uncrustify -l CPP -c source/scripts/uncrustify.cfg --replace $file";
}
