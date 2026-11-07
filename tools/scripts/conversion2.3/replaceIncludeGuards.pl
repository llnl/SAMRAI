#! /usr/bin/perl
##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/scripts/source_manipulation/replaceIncludeGuards.pl $
## Package:     SAMRAI scripts
## Copyright:   (c) 1997-2026 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1917 $
## Description: perl script to update Xd sed files to templates on DIM
##

use File::Basename;
use File::Find;
use Cwd;

#
# Disallow running from certain directories.
#

my $pwd = cwd;

my $debug = 1;

#
# Get prefix pattern to look for, or use default prefix pattern.
#

my $prepat;
$prepat = q|(.*\.[ChI]$)|;
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

if(1) {
    $xdir = dirname $xfile;
    print "xdir: $xdir\n" if ( $debug > 0 );

    ( $dfile = basename $xfile ) =~ s/(.*)X([-\.].*)\.sed$/\1X\2/o;
    $dfile .= ".tmp";

    print "File Suffix " . substr($dfile, -1, 1) . "\n" if ($debug > 0);
    print "dfile: $dfile\n" if ( $debug > 0 );
    
    open XF, "< $xfile" || die "Cannot open file $xfile";
    open TF, "> $dfile" || die "Cannot open file $tpath";

    # read in more whole paragraph, not just one line
#    $/ = '';
    undef $/;

    while ( <XF> ) {
#	s/#ifndef (.*)$//mgi;
#	s/^#ifndef (\w*)(\s*)$#include (.*)//smgi;
	s/^#ifndef DEBUG_NO_INLINE/#FIXUP/smgi;
	s/^#ifndef include(\w*)(\s*)#include \"([\w\.\/]*)\"(\s*)#endif/#include \"$3\"/smgi;
	s/^#FIXUP/#ifndef DEBUG_NO_INLINE/smgi;
	print TF;
    }

    close XF || die "Cannot close file $xfile";
    close TF || die "Cannot close file $dfile";

    printf "unlink $xfile\n"  if ( $debug > 0 );
    unlink($xfile);
    printf "rename $dfile $xfile\n"  if ( $debug > 0 );
    rename( $dfile, $xfile);
}
}





