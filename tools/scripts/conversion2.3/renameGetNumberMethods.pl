#! /usr/bin/perl
##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/tools/scripts/conversion2.0/renameXd.pl $
## Package:     SAMRAI scripts
## Copyright:   (c) 1997-2026 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1917 $
## Description: perl script to rename getNumber methods to be getNumberOf 
##

use File::Basename;
use File::Find;
use Cwd;

#
# Disallow running from certain directories.
#

$patternFile="getNumberMethods.data";

my $pwd = cwd;

#
# Read in datafile to get the names for replacement
#
my $patternFile = (dirname $0) . "/$patternFile";
open(PATTERNFILE, "$patternFile") || die "Cannot open input sed file $patternFile";
while (<PATTERNFILE>) {
    if ( m/^([^#][^ ]+)\n/o ) {
      push @replacePatterns, ${1};
    }
}
close(PATTERNFILE);

my $fileExtensionPattern;
$fileExtensionPattern = q/(.*\.[ChI]$)|(.*\.CPP$)|(.*\.cpp$)|(.*\.cxx$)|(.*\.CXX$)|(.*\.H$)|(.*\.hxx$)|(.*\.Hxx$)|(.*\.HXX$)/;

#
# Find the files to convert.
#
# Excludes files that are in internal source code control directories.
#
@filesToProcess = ();
sub selectFiles {
    if ( $File::Find::name =~ m!/(.svn|CVS|\{arch\})$!o ) {
	$File::Find::prune = true;
    }
    elsif ( -f && m/$fileExtensionPattern/o ) {
	push @filesToProcess, $File::Find::name;
	$filesToProcess[$#filesToProcess] =~ s|^\./||o;
    }
}
find( \&selectFiles, '.' );

for $file (@filesToProcess) {
    print "Working on $file\n";
    $directory = dirname $file;

    $filebasename = basename $file;

    $tempFile = $filebasename . ".samrai.tmp";

    open FILE, "< $file" || die "Cannot open file $file";
    open TEMPFILE, "> $tempFile" || die "Cannot open temporary work file $tempFile";
    while ( $str = <FILE> ) {

	for $pattern (@replacePatterns) {
	    $new_name = ${pattern};
	    $new_name =~ s/getNumber/getNumberOf/g;
	    $str =~ s/${pattern}/${new_name}/g;
	}

	print TEMPFILE $str;
    }

    close FILE || die "Cannot close file $file";
    close TEMPFILE || die "Cannot close file $tempFile";

    unlink($file);
    rename( $tempFile, $file);
}





