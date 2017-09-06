#!/usr/bin/perl  

## Version 1.00

## Introduction.
## Split the maf.gz file into maf file by target name.
## Store in the directory named as the input file name.

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 1 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

  This script is to create individual maf files named after the target names.
A gzip maf file should be provided.
  Output goes to the directory named after the input gzip maf file.	

  --help or no arguments - show this message, ignore undefined inputs

  NOTE1: The record in the gzip maf file should be ordered by target name !!!!
  NOTE2: This script can only process paiwise maf file (a targe and a query) !
  NOTE3: No single line alignment is allowed !
  NOTE4: The prexisting files under the same directory will be overwritten !!! 
   
  USAGE: mafSplitByTarget.pl input.maf.gz  
	 
USAGES

print "========== Start at "; system("date"); print "\n\n";
print "$0 $Arg_list\n";
my $timer = time();

#set the output_field_separator to " " and autoflushing
$,=' ';
$|=1;

# set the over-write flag, to over-write any existing files
my $Force=0;
if ($Arg_list =~ m/--Force/) {
	$Force = 1;
	print "Set to OVER-WRITING mode!\n";
}

my ($out_dir) = $ARGV[0] =~ m/(.+).gz/;
die "The name of $ARGV[0] is not correct ! Die !!!\n" unless defined($out_dir);
my $in_fileH;
open($in_fileH, "gunzip -c $ARGV[0] |") or die "Can't open file $ARGV[0] ! Die !!!\n";
unless (-d $out_dir) { mkdir $out_dir; }

my $target_name = '';
my $out_fileH;
my ($aline, $tline, $qline);
while (<$in_fileH>){
	
	next unless m/^a\s+score=/;
	
	$aline = $_;
	$tline = <$in_fileH>;
	$qline = <$in_fileH>;
	
	$tline =~ m/^s\s+[-\w]+\.([-\w]+)/;
	if ($target_name ne $1){
		$target_name=$1;
		open($out_fileH, ">>$out_dir\/$target_name.maf") or die "Can't open $out_dir\/$target_name.maf ! Die !!!\n";
		print $out_fileH "##maf version=1 scoring=blastz\n";
	}
	
	print $out_fileH $aline;
	print $out_fileH $tline;
	print $out_fileH $qline;
	print $out_fileH "\n";		

}



print "\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";

