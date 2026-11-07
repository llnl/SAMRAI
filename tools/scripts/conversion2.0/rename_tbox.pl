#! /usr/bin/perl
##
## File:        $URL$
## Package:     SAMRAI scripts
## Copyright:   (c) 1997-2026 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision$
## Description: perl script to update Xd sed files to templates on DIM
##

use File::Basename;
use File::Find;
use Cwd;

#
# Disallow running from certain directories.
#

$prefix=@ARGV[0];

my $pwd = cwd;
#die basename($0) . " should not be run from your current directory"
#    if $pwd =~ m:\b(examples|source/test|source/scripts)(/|$):;


my $debug = 0;


#
# Read in sed.data to get the substitution patterns for X strings.
#

my $datfile = (dirname $0) . "/$prefix.data";
print "datfile: $datfile\n" if ( $debug > 0 );
open(DATFILE, "$datfile") || die "Cannot open input sed file $datfile";
while (<DATFILE>) {
    if ( m/^([^#][^ ]+)\n/o ) {
      push @Xpattern, ${1};
    }
}
close(DATFILE);
my $allXpatterns = join '|', @Xpattern;

my $allSimplepatterns = join '|', @Xpattern;
$allSimplepatterns =~ s/($prefix)_//go;

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

    if ( $File::Find::name =~ m!/(.svn|CVS|[123]d|\{arch\})$!o ) {
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
    print "Working on $xfile\n";
    $xdir = dirname $xfile;
    print "xdir: $xdir\n" if ( $debug > 0 );

    $xfilebasename = basename $xfile;
    print "<$xfilebasename>\n" if ( $debug > 0);

    ( $dfile = basename $xfile ) =~ s/(.*)X([-\.].*)\.sed$/\1X\2/o;
    $dfile .= ".tmp";

    ( $newclassname = basename $xfile ) =~ s/(.*)X([-\.].*)\.sed$/\1/o;
    ( $oldclassname = basename $xfile ) =~ s/(.*)X([-\.].*)\.sed$/\1X/o;

    print "File Suffix " . substr($dfile, -1, 1) . "\n" if ($debug > 0);
    print "dfile: $dfile\n" if ( $debug > 0 );
    print "oldclassname: $oldclassname\n" if ( $debug > 0 );
    print "newclassname: $newclassname\n" if ( $debug > 0 );
    
    open XF, "< $xfile" || die "Cannot open file $xfile";
    open TF, "> $dfile" || die "Cannot open file $tpath";
    while ( $str = <XF> ) {

	# replace filenames for includes
#        $str =~ s/(\#include.*)\"($allSimplepatterns)\.([hIC])\"/\1\"sgs\2\.\3/go;
#        $str =~ s/(\#include)(\s*)\"(.*)\"/\1\2\"tbox\/\3"/go;
        $str =~ s/\#include(\s*)\"($allSimplepatterns)\.([hIC])\"/#include\1\"tbox\/\2\.\3"/go;
	print TF $str;
    }

    close XF || die "Cannot close file $xfile";
    close TF || die "Cannot close file $dfile";

    printf "unlink $xfile\n" if ( $debug > 0 );
    unlink($xfile);
    printf "rename $dfile $xfile\n" if ( $debug > 0 );
    rename( $dfile, $xfile);
}





