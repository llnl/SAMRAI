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

	# Special class name changes
	$str =~ s/tbox_Box/tbox::DatabaseBox/go;

	# Special method changes
#	$str =~ s/getBox/getDatabaseBox/go;
#	$str =~ s/putBox/putDatabaseBox/go;
#	$str =~ s/getBoxArray/getDatabaseBoxArray/go;
#	$str =~ s/putBoxArray/putDatabaseBoxArray/go;

	# Included lines are the include guards which should retain underscores

	$included=0;
        if ( $str =~ s/included_($allXpatterns)[X123]/included_\1/go ) {
	    $included = 1;
	}

        if ( $str =~ m/included_$prefix/go ) {
	    $included = 1;
	}

	# File lines are statements that are including a file
	$fileline=0;
	# replace filenames for includes
        if ( $str =~ s/($allXpatterns)[X123]\.([hI])/\1\.\2/go ) {
	    $fileline = 1;
	}

	if ( $str =~ m/^\#include/ ) {
	    $fileline = 1;
	}

	if ($fileline) {
	    $str =~ s/($prefix)_//go;

	    $str =~ s/($allSimplepatterns)1/\1/go;
	    $str =~ s/($allSimplepatterns)2/\1/go;
	    $str =~ s/($allSimplepatterns)3/\1/go;
	}
	
	# replace class forward declares
#	$str =~ s/class ($allXpatterns)X;/template<int DIM> class \1;/go;

	# replace classnames

	if(!$fileline && !$included) {

	    # replace documentation references
	    $str =~ s/see $prefix_($allSimplepatterns)X/see $prefix\_\1/go;

	    if ($included) {
		print "this is bogus\n";
	    }
	    # If the class is a template, need <DIM,TYPE>
	    $str =~ s/($allXpatterns)X<(.*)>/\1<NDIM,\2>/go;
	    $str =~ s/($allXpatterns)1<(.*)>/\1<1,\2>/go;
	    $str =~ s/($allXpatterns)2<(.*)>/\1<2,\2>/go;
	    $str =~ s/($allXpatterns)3<(.*)>/\1<3,\2>/go;
	    
	    # Class is not already templated, just add <DIM>
	    $str =~ s/($allXpatterns)X/\1<NDIM>/go;
	    $str =~ s/($allXpatterns)1/\1<1>/go;
	    $str =~ s/($allXpatterns)2/\1<2>/go;
	    $str =~ s/($allXpatterns)3/\1<3>/go;
	    
	    
	    $replaced = 0;

 	    # If the class is a template, need <DIM,TYPE>
 	    if ($str =~ s/(\W)($allSimplepatterns)<(.*)>/\1$prefix::\2<NDIM,\3>/go )
	    {
		$replaced = 1;
	    }
	    $str =~ s/(\W)($allSimplepatterns)X<(.*)>/\1$prefix::\2<NDIM,\3>/go;
	    $str =~ s/(\W)($allSimplepatterns)1<(.*)>/\1$prefix::\2<1,\3>/go;
	    $str =~ s/(\W)($allSimplepatterns)2<(.*)>/\1$prefix::\2<2,\3>/go;
	    $str =~ s/(\W)($allSimplepatterns)3<(.*)>/\1$prefix::\2<3,\3>/go;
	    
	    # Class is not already templated, just add <DIM>
	    if (!$replaced) {
		$str =~ s/(\W)($allSimplepatterns)(\W)/\1$prefix::\2<NDIM>\3/go;
	    }
	    $str =~ s/(\W)($allSimplepatterns)X/\1$prefix::\2<NDIM>/go;
	    $str =~ s/(\W)($allSimplepatterns)1/\1$prefix::\2<1>/go;
	    $str =~ s/(\W)($allSimplepatterns)2/\1$prefix::\2<2>/go;
	    $str =~ s/(\W)($allSimplepatterns)3/\1$prefix::\2<3>/go;
	    
	    # fix up nested templates 
	    $str =~ s/\<NDIM\>\>/\<NDIM\> \>/o;
	    $str =~ s/\<1\>\>/\<1\> \>/o;
	    $str =~ s/\<2\>\>/\<2\> \>/o;
	    $str =~ s/\<3\>\>/\<3\> \>/o;
	    
	    # strip off package names

	    $str =~ s/($prefix)_([A-Z])/\1::\2/go;

	    # fix up plog with package name

	    $str =~ s/(\s)pout/\1tbox::pout/go;
	    $str =~ s/(\s)perr/\1tbox::perr/go;
	    $str =~ s/(\s)plog/\1tbox::plog/go;
	}

	print TF $str;
    }

    close XF || die "Cannot close file $xfile";
    close TF || die "Cannot close file $dfile";

    printf "unlink $xfile\n"  if ( $debug > 0 );
    unlink($xfile);
    printf "rename $dfile $xfile\n"  if ( $debug > 0 );
    rename( $dfile, $xfile);
}





