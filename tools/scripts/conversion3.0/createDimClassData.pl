#! /usr/bin/perl
##
## File:        $URL$
## Package:     SAMRAI scripts
## Copyright:   (c) 1997-2026 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1917 $
## Description: Create list of classes templated on DIM.  Run in v2 
##              directory to generate files.
##

use strict;

use File::Basename;
use File::Find;
use Cwd;

# Flush I/O on write to avoid buffering
$|=1;

my $debug=0;

my $end_of_line = $/;

#
# Remove duplicated values
#
sub unique {
    foreach my $test (@_){
	my $i = -1;
	my @indexes = map {$i++;$_ eq $test ? $i : ()} @_;
	shift @indexes;
	foreach my $index (@indexes){
	    splice(@_,$index,1);
	}
    }
    return @_;
}

my $pwd = cwd;


#=============================================================================
# Determine classes that are templated on dimension.
#=============================================================================
my $fileExtensionPattern;
$fileExtensionPattern = q/.*\.h$/;

#
# Find the files to convert.
my @filesToProcess = ();
print @filesToProcess if $debug > 2;

#
# Build list of files in which to look for templates.
#
sub selectHFiles {
    if ( $File::Find::name =~ m!/(.svn|CVS|include|scripts|\{arch\})$!o ) {
	$File::Find::prune = 1;
    }
    elsif ( -f && m/.*\.h$/o ) {
	push @filesToProcess, $File::Find::name;
    }
}
find( \&selectHFiles, '.' );

my $templatesClassNamesPattern = "";

#
# For classes with templating on DIM and other classes.
#
my %templatesInFile = ();
my @templatesOnDIM = ();
my @fileWithTemplatesOnDIM = ();

#
# For classes with templating only DIM.
#
my %templatesOnDimOnlyInFile = ();
my @templatesOnDIMOnly = ();
my @fileWithTemplatesOnDIMOnly = ();

#
# Find list of classes that are templated on dimension and build
# list of those classes, files with those classes and
# hash of files to classes.
#
for my $file (@filesToProcess) {
    print "Looking for templates on DIM in $file\n" if $debug > 1;
    my $directory = dirname $file;
    my $filebasename = basename $file;

    open FILE, "< $file" || die "Cannot open file $file";
    
    # read in whole paragraph, not just one line
    undef $/;
    while ( <FILE> ) {
	# This pattern could use improvement to better find class
        # declarations vs foward declaration.  Currently assumes a
        # space after classname will suffice to distinquish which is
        # not correct.
	while ( $_ =~ m/template\s*<\s*int\s*DIM\s*>\s*class\s+(\w+)\s+/sgm ) {
	    print "Found template on DIM only $1\n" if $debug > 1;
	    $templatesClassNamesPattern .= "$1|";

	    push @templatesOnDIMOnly, $1;
	    push @fileWithTemplatesOnDIMOnly, $file;
	    push @{$templatesOnDimOnlyInFile{$file}} ,$1;
	}
    }
    close FILE;

    open FILE, "< $file" || die "Cannot open file $file";
    
    while ( <FILE> ) {
	# This pattern could use improvement to better find class
	# declarations vs foward declaration.  Currently assumes a
	# space after classname will suffice to distinquish which is
	# not correct.
	while ( $_ =~ m/template\s*<\s*int\s*DIM,\s*.*?>\s*class\s+(\w+)\s+/sgm ) {
	    print "Found template on DIM with other $1\n" if $debug > 1;
	    $templatesClassNamesPattern .= "$1";
	    push @templatesOnDIM, $1;
	    push @fileWithTemplatesOnDIM, $file;
	    push @{$templatesInFile{$file}} ,$1;
	}
    }
    close FILE;
}

my $filename="SAMRAI_classes_templated_on_DIM_only.txt";
open FILE, '>', $filename or die "Can't open $filename : $!";
foreach ( @templatesOnDIMOnly )
{
    print FILE "$_\n";
}
close FILE;


my $filename="SAMRAI_classes_templated_on_DIM_and_other.txt";
open FILE, '>', $filename or die "Can't open $filename : $!";
foreach ( @templatesOnDIM )
{
    print FILE "$_\n";
}
close FILE;


