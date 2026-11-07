#! /usr/bin/perl
#########################################################################
##
## This file is part of the SAMRAI distribution.  For full copyright 
## information, see COPYRIGHT and LICENSE. 
##
## Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
## Description:   perl script to update Xd sed files to templates on DIM 
##
#########################################################################

use File::Basename;
use File::Find;
use Cwd;

my $pwd = cwd;
#die basename($0) . " should not be run from your current directory"
#    if $pwd =~ m:\b(examples|source/test|source/scripts)(/|$):;


my $debug = 0;

#
# Read in sed.data to get the substitution patterns for X strings.
#

#
# Get prefix pattern to look for, or use default prefix pattern.
#

my $prepat;
$prepat = q|(.*\.h(.sed)?$)|;
print "prepat: $prepat\n" if ( $debug > 0 );


#
# Find the h files to search
#

@allfiles = ();
sub selectXfile {
    # This subroutine selects the X files in a find command.
    # print "-$File::Find::dir- -$File::Find::name-\n";
#    if ( $File::Find::dir =~ m!/(examples|source/test|source/scripts|CVS|[123]d)$!o ) {

    if ( $File::Find::dir =~ m!/(.svn|CVS|[123]d)$!o ) {
	# print "pruned\n";
	$File::Find::prune = true;
    }
    elsif ( -f && m/$prepat/o ) {
	push @allfiles, $File::Find::name;
	$allfiles[$#allfiles] =~ s|^\./||o;
    }
}
print "Scanning...\n" if ( $debug > 0 );
find( \&selectXfile, '.' );
print "Done.\n" if ( $debug > 0 );

for $xfile (@allfiles) {
#    print "Working on $xfile\n";
    $xdir = dirname $xfile;
    print "xdir: $xdir\n" if ( $debug > 0 );

    print "File Suffix " . substr($dfile, -1, 1) . "\n" if ($debug > 0);
    print "dfile: $dfile\n" if ( $debug > 0 );
    
    open XF, "< $xfile" || die "Cannot open file $xfile";
    while ( <XF> ) {
	# This strips off comments 
	$/ = undef;
	s#/\*[^*]*\*+([^/*][^*]*\*+)*/|//[^\n]*|("(\\.|[^"\\])*"|'(\\.|[^'\\])*'|.[^/"'\\]*)#$2#gs;  

	if ( m/.*class(\s+)(\w+)(\s+).*/g ) {
	    print "$2\n";
	}

    }

    close XF || die "Cannot close file $xfile";
}





