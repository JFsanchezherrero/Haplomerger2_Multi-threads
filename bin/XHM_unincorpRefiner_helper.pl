#!/usr/bin/perl  

## Version 2.30

## Introduction.


use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 1 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit; }

  This script is to recreate the hm.unpaired_updated for a given file. 
        
	the recreated file is put in the current directory.
	
	?/XHM_unincorpRefiner_helper.pl sequence.fa.gz
	
USAGES

##0:new_scaffold_name	1:old_scaffold_name	2:size	3:start	4:len	5:full_scaffold_or_part_of_the_scaffold	6:N_counts	7:lowcase_count
my @records;

my $inFH;
if($ARGV[0]=~/\.gz$/){
	open($inFH, "gunzip -c $ARGV[0] |") or die "Can not open gunzip -c $ARGV[0] | .\n";
}else{
  open($inFH, "< $ARGV[0]") or die "Can not open $ARGV[0] .\n";
}


my ($line,$name,$desc,$seq,$is_the_end)=("","","","",0);

while($line=<$inFH>){
	last if $line =~ m/^>/;
}	
die "Problems with the fasta file; no symbol > found !\n" unless defined($line);

while($is_the_end==0){
	
	($name,$desc,$seq)=('','','');
	
	if($line =~ m/>\s*(\S+)/){
		$name=$1;
	}
	if($line =~ m/>\s*(\S+)\s+(.+)/){
		$desc=$2;
		chomp $desc;
	}	

	while($line=<$inFH>){
		last if $line =~ m/^>/;
		chomp $line;
		$seq .= $line;
	}
	$is_the_end=1 unless defined($line);
	
	my $tt;
	($tt->[0],$tt->[1])=($name,$name);
	$tt->[2]=length($seq);
	$tt->[3]=0;
	$tt->[4]=$tt->[2];
	$tt->[5]='full';
	$tt->[6]= $seq =~ tr/Nn//;
	$tt->[7]= $seq =~ tr/atcg//;
	
	push @records,$tt;

}

open(my $outFH, ">hm.unpaired_updated") or die "Can not open hm.unpaired_updated.\n";
print $outFH "# This table is created by XHM_unincorpRefiner_helper.pl .\n";
print $outFH "# This table is not intended to be manually modified.\n\n";
print $outFH "#new_scaffold_name	old_scaffold_name	size	start	len	full_scaffold_or_part_of_the_scaffold	N_counts	lowcase_count\n";
$"="\t"; #"
foreach (@records){
	print $outFH "@$_\n";
}