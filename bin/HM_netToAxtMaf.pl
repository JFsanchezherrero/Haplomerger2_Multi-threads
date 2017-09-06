#!/usr/bin/perl  

## Version 1.00

## Introduction.
## Provide a net file and its corresponding chain(.gz) file.
## create maf from Net.
## Each time only one pair of species can be processed.

## Recommand log file: species1.species2.result/_synNetMaf.log

## Assumptions.
## 1)The current/home directory is the working directory.
## 2)size files.
## 3)lastz and axtChainRecipBestNet.pl have been run.
## 4)This script is in its very preliminary form, so users of this script take their own risk. 
## 5)the script does not check for argument syntax, the users need to ensure correct syntax. 

## Dependences.
## axtChain, chainAntiRepeat, chainMergeSort, chainStitchId, chainNet, NetSyntenic, chainSwap ...

## Update
## V1.00. 
## creation of  this script.


use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 1 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

	This script is to create maf from net file.
	
  Required excutables http://genome.ucsc.edu/admin/jksrc.zip,
including netToAxt, axtSort, axtFilter, and axtToMaf.	

   --help or no arguments - show this message. Undefined inputs are ignored.

   --Species -   Mandatory, to provide species names which are wanted to be 
                 processed. NOTE THAT the species name should be given 
                 in ORDER, because the species name comes first will be used
                 as the target species.
                 Note: each run only allow two species.
                 Example: --Species human chimp
           
   --netFile=<file_name> - specify the net file name,
                 do not include the path! The default is 'all.rbest.Net.gz'.
                 chain file name will be deduced from the net file name,
                 chain file can be .gz or not.
   --chainFile=<file_name> - normally the chain file is under same prefix
                 with the net file
                 This option can force to used the chain
                 file specified.
   --recreateChain - recreate chain.gz corresponding for the provided net file
                 Using provided net and chain files
                 
   --createMirrorAxt - to create mirror axt
   							 This function requires script HM_axtToPerfectMirrorAxt.pl.
   --minScore=N - filter out the alignment whose score is lower 
                 than N (default=0)             
   --axt -       output axt file instead of maf   
                                   
   --Force -     over-writing existing files (default=no)
   
Notes:
1)The current/home directory is the working directory.
2)lastz and axtChainRecipBestNet.pl have been run.
3)The results maf files are put into ./species1.species2.result.  
	 
USAGES

print "\n\n========== Start at "; system("date");
print "$0 $Arg_list\n\n";
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

# Check and store species names
my @Species;
unless ($Arg_list =~ m/--Species\s+([^-]*)/) {die "No --Species argument or species_names found! Die!\n" }
unless (@Species = $1 =~ m/(\w+)/g)  {die "No species names found! Die!\n" };
unless (scalar(@Species) == 2) { die "only two species are allowed! Die!\n"; } 
print "Species included: ", @Species, "\n";

sub netToMaf;

netToMaf;

print "\n\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";

######################### subroutine defination #####################################

sub netToMaf {
	
	print "=== netToMaf start ===\n";
	
	my ($target,$query) = @Species[0,1];
	my $path = "$target.$query.result";
	
	# Store and check net file and chain.gz file
	my ($net_file, $chain_file);	#source 
	if ($Arg_list =~ m/--netFile=([-\.\w]+)\.net\.gz/) {
		if(-f "$path/$1.net.gz"){
			system("gunzip -c $path/$1.net.gz >$path/$1.net");
		}else{
			die "Can not open $path/$1.net.gz!\n"; 
		}
		$net_file =  "$path/$1.net";
		$chain_file = "$path/$1.chain.gz";
	}elsif($Arg_list =~ m/--netFile=([-\.\w]+)\.net/) {
		$net_file =  "$path/$1.net";
		$chain_file = "$path/$1.chain.gz";
	}else{
		$net_file = "$path/all.rbest.net";
		$chain_file = "$path/all.rbest.chain.gz";
	}
	
	if ($Arg_list =~ m/--chainFile=([-\.\w]+\.chain\.gz)/) { 
		my $temp = $1;  
		if ($Arg_list =~ m/--recreateChain/ and $Arg_list =~ m/--netFile/){
			if (-f $chain_file and $Force == 0) { die "No chain file recreation needed! $chain_file is already existed! die!\n"; }
			my $cmd  = "gunzip -c $path/$temp ";
			$cmd .= "| netChainSubset $net_file /dev/stdin /dev/stdout 2>>$path/_netChainSubset.log ";
			$cmd .= "| chainStitchId /dev/stdin /dev/stdout 2>>$path/_chainStitchId.log ";
			$cmd .= "| gzip -c >$chain_file";
			print "Recreate the corresponding chain file ($chain_file) for the provided net file ($net_file) ...\n";
			print $cmd,"\n";
			system($cmd);			
		}else{
			$chain_file = "$path/$temp";
		}
	}
	
	print "Checking the following two files ...\n$chain_file\n$net_file\n";
	unless (-f $chain_file and -f $net_file) { print "Files are not found! Die!\n"; die "\n"; }
	
	# checking sizes files
	print "Checking missing sizes files ...\n";
	my $missing_sizes =0;
	unless (-f "$target.sizes") { $missing_sizes =1; print "$target.sizes\n"; }
	unless (-f "$query.sizes") { $missing_sizes =1;  print "$query.sizes\n"; }
	die "File missing! Die! \n" if ($missing_sizes == 1);
	
	# checking nib files
	my $missing_nib = 0;
	my $sizeFH;
	foreach my $temp (@Species) {
		print "checking nib files in $temp.seq/ ...\n";
		open($sizeFH, "<$temp.sizes") or die "Can not open file $temp.sizes! Die!\n";
		while (<$sizeFH>) {
			next unless m/([-\.\w]+)\t/;
			unless (-f "$temp.seq/$1.nib") { $missing_nib=1; print "$1.nib\n"; }
		}
		close $sizeFH;
	}
	die "Some nib files are missing! Die!\n" if ($missing_nib == 1);
	
	my ($axt_file,$maf_file);
	# checking old net.maf.gz files
	if ($Arg_list =~ m/--netFile=([-\.\w]+)\.net/) {
		$maf_file = "$path/$1.net.maf.gz";
		$axt_file = "$path/$1.net.axt.gz";
	}else{
		$maf_file = "$path/all.rbest.net.maf.gz";
		$axt_file = "$path/all.rbest.net.axt.gz";
	}
	print "checking net.maf/axt file $maf_file and $axt_file ...\n";
	if (-f $maf_file or -f $axt_file) {
		$Force ? print "File $maf_file or $axt_file existed, will be over-written!\n" : die "File existed! Die!\n";
	}

	
	# creating the maf/axt file
	my $cmd;
	if($Arg_list !~ m/--createMirrorAxt/){
		$cmd  = "gunzip -c $chain_file ";
		$cmd .= "| netToAxt -maxGap=300 $net_file /dev/stdin $target.seq $query.seq /dev/stdout ";
		$cmd .= "| axtSort /dev/stdin /dev/stdout ";
	
		if($Arg_list =~ m/--minScore=(\d+)/){
			$cmd .= "|axtFilter -minScore=$1 /dev/stdin ";
		}
	
		if($Arg_list =~ m/--axt/){
			$cmd .= "| gzip -c > $axt_file";
		}else{
			$cmd .= "| axtToMaf -tPrefix=$target. -qPrefix=$query. /dev/stdin $target.sizes $query.sizes /dev/stdout";
			$cmd .= "| gzip -c > $maf_file";
		}
	
		print $cmd,"\n";
		system($cmd);
	}else{
		print "Creating mirrored axt and maf.\n";
		$cmd  = "gunzip -c $chain_file ";
		$cmd .= "| netToAxt -maxGap=300 $net_file /dev/stdin $target.seq $query.seq /dev/stdout ";
	
		if($Arg_list =~ m/--minScore=(\d+)/){
			$cmd .= "|axtFilter -minScore=$1 /dev/stdin ";
		}
		
		if($Arg_list =~ m/--axt/){
			$cmd .= "| HM_axtToPerfectMirrorAxt.pl $target.sizes $query.sizes /dev/stdin | axtSort /dev/stdin /dev/stdout | gzip -c > $axt_file";
		}else{
			$cmd .= "| HM_axtToPerfectMirrorAxt.pl $target.sizes $query.sizes /dev/stdin | axtSort /dev/stdin /dev/stdout ";
			$cmd .= "| axtToMaf -tPrefix=$target. -qPrefix=$query. /dev/stdin $target.sizes $query.sizes /dev/stdout";
			$cmd .= "| gzip -c > $maf_file";
		}
		
		print $cmd,"\n";
		system($cmd);
	}	
	
	if ($Arg_list =~ m/--netFile=([-\.\w]+)\.net\.gz/){
		unlink "$path/$1.net";
	}	 
	print "=== netToMaf done! ===\n\n";
}
