#!/usr/bin/perl  


## notes
## The agp file's coordinates are 0-base; the cabog6.1 uses agp specification v1.1.
## Only gaps between contig are considered.

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 0 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

   This script is to counting contig and scaffold N50 from an agp(v1.1)
   files (tested for ca6.1).
   Only gaps between contigs are count.
   Connected contigs are counted as one.
   Scaffold N50 with or without gaps.
   
   The count is just approximate!
   
   using stdin and stdout.
    	 
USAGES

my (@ctgs,@scfs,@scfs_gaps);
my ($pre_scf);
my ($len0,$len1,$len2)=(0,0,0);
	
while(<STDIN>){
	next if m/^#|^\s+/;
	chomp;
	my @tt=split /\t/,$_;
	if(!defined($pre_scf) or $pre_scf ne $tt[0]){
		push @ctgs,$len0 if defined($pre_scf);
		push @scfs,$len1 if defined($pre_scf);
		push @scfs_gaps,$len2 if defined($pre_scf);
		($len0,$len1,$len2)=(0,0,0);
		$pre_scf=$tt[0];
	}
	$len2+=$tt[2]-$tt[1]+1;
	$len1+=$tt[2]-$tt[1]+1 if $tt[4] eq 'W';
	if($tt[4] eq 'N'){
		push @ctgs,$len0;
		$len0=0;
	}else{
		$len0+=$tt[2]-$tt[1]+1;
	}
}

@ctgs = sort { $b <=> $a } @ctgs;
@scfs = sort { $b <=> $a } @scfs;
@scfs_gaps = sort { $b <=> $a } @scfs_gaps;

my ($type,$tt,$total,$N50);
$type = "contig N50:";
$tt=\@ctgs;
($total,$N50)=(0,0);
for(my $i=0;$i<scalar(@$tt);$i++){
	$total+=$tt->[$i];
}
for(my $i=0;$i<scalar(@$tt);$i++){
	$N50+=$tt->[$i];
	if($N50>=$total/2){
		print "total $total; $type $tt->[$i].\n\n";
		last;
	}
}

$type = "scaffold N50:";
$tt=\@scfs;
($total,$N50)=(0,0);
for(my $i=0;$i<scalar(@$tt);$i++){
	$total+=$tt->[$i];
}
for(my $i=0;$i<scalar(@$tt);$i++){
	$N50+=$tt->[$i];
	if($N50>=$total/2){
		print "total $total; $type $tt->[$i].\n\n";
		last;
	}
}

$type = "scaffold with gaps N50:";
$tt=\@scfs_gaps;
($total,$N50)=(0,0);
for(my $i=0;$i<scalar(@$tt);$i++){
	$total+=$tt->[$i];
}
for(my $i=0;$i<scalar(@$tt);$i++){
	$N50+=$tt->[$i];
	if($N50>=$total/2){
		print "total $total; $type $tt->[$i].\n\n";
		last;
	}
}
