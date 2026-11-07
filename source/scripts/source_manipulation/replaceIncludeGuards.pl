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

#
# Disallow running from certain directories.
#

my $pwd = cwd;
#die basename($0) . " should not be run from your current directory"
#    if $pwd =~ m:\b(examples|source/test|source/scripts)(/|$):;

my $debug = 0;

#
# Read in sed.data to get the substitution patterns for X strings.
#

my $datfile = (dirname $0) . "/headers.data";
print "datfile: $datfile\n" if ( $debug > 0 );
open(DATFILE, "$datfile") || die "Cannot open input sed file $datfile";
while (<DATFILE>) {
    if ( m/^\#/o ) {
    } else {
	if ( m/(\w*)\s+(\w*)/o ){ 
	    $header_to_package{$2} = $1;
	    push @hpattern, $2;
	}
    }
    
}
close(DATFILE);

my $headerpatterns = join '|', @hpattern;
print "hpatterns = $headerpatterns\n" if $debug;

#
# Get prefix pattern to look for, or use default prefix pattern.
#

my $prepat;
$prepat = q|(.*\.[ChI](.sed)?$)|;
print "prepat: $prepat\n" if ( $debug > 0 );


#
# Find the X files to convert.
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
    print "Working on $xfile\n" if $debug;
    $xdir = dirname $xfile;
    print "xdir: $xdir\n" if ( $debug > 0 );

    ( $dfile = basename $xfile ) =~ s/(.*)X([-\.].*)\.sed$/\1X\2/o;
    $dfile .= ".tmp";

    print "File Suffix " . substr($dfile, -1, 1) . "\n" if ($debug > 0);
    print "dfile: $dfile\n" if ( $debug > 0 );
    
    open XF, "< $xfile" || die "Cannot open file $xfile";
    open TF, "> $dfile" || die "Cannot open file $tpath";
    while ( <XF> ) {

	if ( m/included\_($headerpatterns)/o ) {
	    $pack=$header_to_package{$1};
	    s/included\_($allSimplepatterns)/included_$pack\_$1/go;
	    print TF;
	} else {
	    print TF;
	}
    }

    close XF || die "Cannot close file $xfile";
    close TF || die "Cannot close file $dfile";

    printf "unlink $xfile\n"  if ( $debug > 0 );
    unlink($xfile);
    printf "rename $dfile $xfile\n"  if ( $debug > 0 );
    rename( $dfile, $xfile);
}





