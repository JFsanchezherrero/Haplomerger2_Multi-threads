#!/usr/bin/perl  

## notes
## to break the misjoins

## bug:
## not found 

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 3 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

   This script is to use the info from the misjoin table to break
   the potential misjoins and output the new scf.
   
   --scf=<scf file> - gzipped fasta file containing the original scaffolds
   
   --misjoin=<misjion info> - hm.assembly_errors
   --aliFilter=N - minimum score required to trigger a misjoin breakage; 
                   default=5000000
   --overhangFilter=N - minimum length of the overhang tails to trigger
                        a misjoin breakage; default=50000bp                

   NOTE: 
   1. In default, for each misjion, this script break both scaffolds involved
      if the misjoin meets the standard set by --aliFilter/overhangFilter.
   2. On the other hand, users can modify the column (*breakpoint_setting)
      in the misjoin file (hm.assembly_errors) to overwrite this default
      behavior: 
      0=default behavior; 
      1=break the target scaffold regardless of --aliFilter/overhangFilter;
      2=break the query scaffold regardless of --aliFilter/overhangFilter;
      3=break both scaffolds regardless of --aliFilter/overhangFilter.               

   
   STDOUT outputs the new scaffold sequences;
   STDERR outputs the logs.
    	 
USAGES


#### ################################
#### reading arguments


my ($in_scf);
if ($Arg_list =~ m/--scf=(\S+)/) {
	$in_scf=$1;
}else{
	die "no --scf arg.\n";
}

my ($in_mis);
if ($Arg_list =~ m/--misjoin=(\S+)/) {
	$in_mis=$1;
}else{
	die "no --misjoin arg.\n";
}

my $ali_filter=5000000;
if ($Arg_list =~ m/--aliFilter=(\S+)/) {
	$ali_filter=$1;
}
print STDERR "--aliFilter=$ali_filter.\n";

my $overhang_filter=50000;
if ($Arg_list =~ m/--overhangFilter=(\S+)/) {
	$overhang_filter=$1;
}
print STDERR "--overhangFilter=$overhang_filter.\n";


#### ################################
#### reading scf

my %scfs;

{
	my $inFH_name=$in_scf;
	my $inFH;
	if($inFH_name=~m/\.gz$/){
		open($inFH,"gunzip -c $inFH_name |") or die "Can not open the file $inFH_name.\n";
	}else{
		open($inFH,"<$inFH_name") or die "Can not open the file $inFH_name.\n"; 
	}
	
	my ($line,$name,$seq,$is_the_end)=("","","",0);
	
	while($line=<$inFH>){
		last if $line =~ m/^>/;
	}	
	die "Problems with the fasta file; no symbol > found !\n" unless defined($line);
	
	while($is_the_end==0){
		
		($name,$seq)=('','');
		
		if($line =~ m/>\s*(\S+)/){
			$name=$1;
		}
	
		while($line=<$inFH>){
			last if $line =~ m/^>/;
			chomp $line;
			$seq .= $line;
		}
		$is_the_end=1 unless defined($line);
		
		#$name=~s/^scf//; ################
		$scfs{$name}=$seq;
	}
	close $inFH;
}

print STDERR "Finished scf fasta reading.\n";


#### ################################
#### reading misjoin

#0:node_id	1:tsc_ali_len	2:ali_len_or_score	3:direction(5or3)	
#4:tsc_name	5:tsc_size	6:tsc_strand	7:ali_start_of_tsc	8:ali_end_of_tsc	9:joining_point_of_tsc 10:overhang_len_of_tsc	
#11:qsc_name	12:qsc_size	13:qsc_strand	14:ali_start_of_qsc	15:ali_end_of_qsc	16:joining_point_of_qsc	17:overhang_len_of_qsc	
#18:*tsc_breakpoint	19:*qsc_breakpoint	20:*breakpoint_setting
#21:tsc_num 22:tsc_ovt 23:tsc_space 24:tsc_N, 25:qsc_num 26:qsc_ovt 27:qsc_space 28:qsc_N 29:tsc_it 30:qsc_it	
my %breaks; #breaks{scf_id}=array of breakpoints

{
	my $inFH_name=$in_mis;
	my $inFH;
	if($inFH_name=~m/\.gz$/){
		open($inFH,"gunzip -c $inFH_name |") or die "Can not open the file $inFH_name.\n";
	}else{
		open($inFH,"<$inFH_name") or die "Can not open the file $inFH_name.\n"; 
	}
	
	while(<$inFH>){
		next if m/^#|^\s+/;
		chomp;
		my @tt =  split /\t/,$_;
		if($tt[20]==1 or $tt[20]==3){
			push @{$breaks{$tt[4]}},$tt[9]; 
		}
		if($tt[20]==2 or $tt[20]==3){
			push @{$breaks{$tt[11]}},$tt[16];
		}
		if($tt[20]==0 and $tt[2]>=$ali_filter and $tt[10]>=$overhang_filter and $tt[17]>=$overhang_filter){
			push @{$breaks{$tt[4]}},$tt[9];
			push @{$breaks{$tt[11]}},$tt[16];
		}
	}
	close $inFH;
	
	my ($total_bp_count,$total_scf_with_bp_count)=(0,0);
	foreach my $tt (keys %breaks){
		push @{$breaks{$tt}},length($scfs{$tt});
		push @{$breaks{$tt}},0;		
		@{$breaks{$tt}} = sort { $a <=> $b } @{$breaks{$tt}};
		$total_scf_with_bp_count++;
		$total_bp_count+=scalar(@{$breaks{$tt}})-2;
	}
	print STDERR "Total scaffolds that suffer from misjoin: $total_scf_with_bp_count .\n";
	print STDERR "Total breakpoints for misjoined scaffolds: $total_bp_count .\n";	
}

print STDERR "Finished misjoin/breaks reading.\n";

#### ################################
#### output scaffolds 
#### 

{
	print STDERR "\n\n========== The coordinates (0-based) of the breakpoints ==========\n";
	
	my $seq;
	my $tt;
	foreach my $scf_id (keys %scfs){
		print STDERR "processing breakpoints of $scf_id: ";
		if(!exists($breaks{$scf_id}) or scalar(@{$breaks{$scf_id}})<=2){
			$seq=$scfs{$scf_id};
			$tt=length($seq);
			print ">$scf_id\n";
			while(length($seq)>0){
				print substr($seq,0,60,'')."\n";
			}
			print STDERR "intact\n";			
		}else{
			print STDERR "0";
			for(my $i=1;$i<scalar(@{$breaks{$scf_id}});$i++){
				$tt=$breaks{$scf_id}->[$i];
				$seq=substr($scfs{$scf_id},$breaks{$scf_id}->[$i-1],$breaks{$scf_id}->[$i]-$breaks{$scf_id}->[$i-1]);
				print ">$scf_id"."_"."$i\n";
				while(length($seq)>0){
					print substr($seq,0,60,'')."\n";
				}
				print STDERR ":$tt";
			}
			print STDERR "\n";
		}
	}
}

print STDERR "\n====== Finished outputing new scf. ========\n\n";