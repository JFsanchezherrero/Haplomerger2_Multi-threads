#!/usr/bin/perl  

## Version 1.00

## Introduction.
## This script runs score matrix inferring by invocate lasz_D.
## This script uses optimized parameters to speed up the inferring process. 
## File "scoreInference.ctl" is not a must, if ignored, default setting will be used.

## Dependences.
## lastz_D, gzip

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 1 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

   This script is simply to run lastz_D to infer score matrix for two genome 
(target and query)	sequences. The obtained score file (target.query.raw.q)
requires further editing for special use.

   Executables lastz_D and gzip are required.
   
Options: (case sensitive!!!)

   --help or no arguments, show this message. Undefined inputs are ignored.
   --identity=<min>[..<max>], The default is '25..75'.
   --target=<target_fa>, target sequence fa file, mandatory.
   --query=<query_fa>, query sequence fa file, mandatory.  
   --Force,      Over-writing existing files is allowed. the default is not.

Notes:
   1) This script is for close species pairs. For divergent species pair, the 
      parameters used in this script may not be sensitive enough.
   2) Fasta files (with suffix .fa) are required.
   3) Gzip files are also accepted.
   4) The current directory should be the working directory, and two .fa files
      should be placed in the current directory.
   5) The target and query names used for "target.query.raw.q" are stripped
      from the .fa file names.
      
Example:
   lastz_D_Wrapper.pl --target=human.fa.gz --query=chimp.fa.gz
   And the result file is "human.chimp.raw.time.xx_xx.q".        
	 
USAGES

print "\n\n========== Start at "; system("date");
print "$0 $Arg_list\n\n";
my $timer = time();

#set the output_field_separator to " "
$,=' ';

# set the over-write flag, to over-write any existing files
my $Force=0;
if ($Arg_list =~ m/--Force/) {
	$Force = 1;
	print "Set to OVER-WRITING mode!\n";
}

# Store species names and to check if all prerequesite data are in proper conditions.
my ($target_file, $target, $target_is_gz);
my ($query_file,  $query,  $query_is_gz );
if ($Arg_list =~ m/--target=([-\.\w]+)/){
	$target_file=$1;
	$target_is_gz = $target_file =~ m/.gz$/ ? 1 : 0;
	$target_file =~ m/^([-\w]+)/;
	$target = $1;
}else{
	die "Target file is not provided !\n";
}
if ($Arg_list =~ m/--query=([-\.\w]+)/){
	$query_file=$1;
	$query_is_gz = $query_file =~ m/.gz$/ ? 1 : 0;
	$query_file =~ m/^([-\w]+)/;
	$query = $1;
}else{
	die "Query file is not provided !\n";
}

unless (-f $target_file and -f $query_file){
	die "Target or query files are not found !\n";
}

if (-f "$target.$query.raw.q" and $Force == 0){
	die "$target.$query.raw.q is already existed !\n";
}

if ($target_is_gz == 1){
	print "Uncompress $target_file to $target.$timer first ...\n";
	system("gunzip -c $target_file >$target.$timer");
	$target_file = "$target.$timer";
}

if ($query_is_gz == 1){
	print "Uncompress $query_file to $query.$timer first ...\n";
	system("gunzip -c $query_file >$query.$timer");
	$query_file = "$query.$timer";
}

my $cmd;
$cmd  = "lastz_D $target_file";
$cmd .= "[multiple] $query_file "; #be careful about the $target_file[multiple]	
$cmd .= "--inferonly --step=20 --seed=12of19 --notransition --notrivial --verbosity=0 ";

my $id_range='25_75';
if ($Arg_list =~ m/--identity=(\d+\.\.\d+)/){
	$cmd .= "--identity=$1 ";
	$id_range = $1;
	$id_range =~ s/\.\./_/;	
  $cmd .= "--infscores=$target.$query.$timer.$id_range.q ";
}else{
	$cmd .= '--identity=25..75 '; 
  $cmd .= "--infscores=$target.$query.$timer.$id_range.q ";	
}

print "$cmd\n";
system($cmd);
print "Score matrix inference is finished!\n";

if (-f "$target.$timer"){
	print "Deleting $target.$timer ...\n";	
	unlink "$target.$timer";
}
if (-f "$query.$timer"){
	print "Deleting $query.$timer ...\n";	
	unlink "$query.$timer";
}

{
	open(my $qInFH,"<$target.$query.$timer.$id_range.q") or die " Can not read $target.$query.$timer.$id_range.q !\n";
	my $qMatrix;
	while(<$qInFH>){
		if(m/^\s*A\s+C\s+G\s+T/ or m/^[ACGT]\s+[0-9-]+/){
			s/^.//;
			$qMatrix .= $_;
		}
	}
	close $qInFH;
	
	open(my $qOutFH,">>$target.$query.$timer.$id_range.q") or die " Can not read $target.$query.$timer.$id_range.q !\n";
	print $qOutFH "\n\n==== For HaploMerger users, you may copy the following score Matrix to scoreMatrix.q ====\n";
	print $qOutFH $qMatrix; 
	close $qOutFH;
}

print "$target.$query.raw.$timer.xx_xx.q is needed to be modified\n";
print "and has its name changed to $target.$query.q in order to feed to axtChain!\n";

print "\n\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";
