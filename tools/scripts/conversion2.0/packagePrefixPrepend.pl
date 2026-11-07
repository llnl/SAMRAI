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

my $datfile = (dirname $0) . "/$prefix.classnames.data";
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
$allSimplepatterns =~ s/($prefix):://go;

print "$allSimplepatterns\n" if $debug;

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

    ( $dfile = basename $xfile ) =~ s/(.*)X([-\.].*)\.sed$/\1X\2/o;
    $dfile .= ".tmp";

    open XF, "< $xfile" || die "Cannot open file $xfile";
    open TF, "> $dfile" || die "Cannot open file $tpath";
    while ( $str = <XF> ) {


	# Do not prepend when inside a set of quotes.  Prepend a funky
	# pattern and strip it off after the generic prepend.
	$str =~ s/\"(.*)($allSimplepatterns)(.*)\"/\"\1ZZZXXXYYY_\2\3\"/g;

	# Prepend package prefix
	# SGS NOTE:  Why does this need to be repeated?  
	$str =~ s/([^a-zA-Z_0-9\:])($allSimplepatterns)(\W)/\1$prefix\:\:\2\3/g;
	$str =~ s/([^a-zA-Z_0-9\:])($allSimplepatterns)(\W)/\1$prefix\:\:\2\3/g;
	$str =~ s/([^a-zA-Z_0-9\:])($allSimplepatterns)(\W)/\1$prefix\:\:\2\3/g;
	$str =~ s/([^a-zA-Z_0-9\:])($allSimplepatterns)(\W)/\1$prefix\:\:\2\3/g;
	$str =~ s/([^a-zA-Z_0-9\:])($allSimplepatterns)(\W)/\1$prefix\:\:\2\3/g;

	$str =~ s/ZZZXXXYYY_($allSimplepatterns)/\1/g;

	print TF $str;
    }

    close XF || die "Cannot close file $xfile";
    close TF || die "Cannot close file $dfile";

    printf "unlink $xfile\n"  if ( $debug > 0 );
    unlink($xfile);
    printf "rename $dfile $xfile\n"  if ( $debug > 0 );
    rename( $dfile, $xfile);
}





