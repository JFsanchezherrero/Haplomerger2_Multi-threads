#!/usr/bin/perl  


## notes
## to break the misjoins
## for testing

# bug: 

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 3 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

   This script is to use the info from the misjoin table to create
   new scf and agp.
   
   *based on the result of XHM_break_misjoin*.pl
   *breaking between contigs!!!!!!!!! having bug for dividing agp items   
   
   Gzip files are supported.
   
   NOTE: 1-break the 1; 2-break the 2; 3-break both; other-break nothing
   
   --agp=<agp file> - input agp v1.1 for scf
   --scf=<scf file> - input scaffold sequences
   --misjoin=<misjion info> - hm.assembly_errors

   NOTE that this script is for ca6.1 original file, which
   means: 1) the scf.fa file have a "scf" prefix for every scf seq;
   2) the agp file have no "scf" prefix for every scf item.   
   
   
   Outputs include (placed in the same dir of the source file):
   *.new.agp.gz; *.new.scf.fa.gz
    	 
USAGES


#### ################################
#### reading arguments

my ($in_agp,$out_agp,$out_scf);
if ($Arg_list =~ m/--agp=(\S+)/) {
	$in_agp=$1;
	$in_agp=~m/(.+)\.agp(\.gz)*$/;
	$out_agp=$1.'.new.agp.gz';
	$out_scf=$1.'.new.scf.fa.gz';
}else{
	die "no --agp arg.\n";
}

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
		
		$name=~s/^scf//;
		$scfs{$name}=$seq;
	}
	close $inFH;
}

print STDERR "Finished scf fasta reading.\n";

#### ################################
#### reading agp

## note that in the ca6 ctg.fasta, ctg name includes a leading 'ctg', but not in the agp. 
# for ctg:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:W, 5:ctg_id, 6:ctg_start, 7:ctg_end, 8:strand, 9:seq]
# for gap:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:N, 5:gap_size, 6:fragment, 7:yes, 8:null, 9:seq]
# need to convert from 1-base to 0 base
my %agp; 

{
	#### reading agp
	my $inFH_name=$in_agp;
	my $inFH;
	if($inFH_name=~m/\.gz$/){
		open($inFH,"gunzip -c $inFH_name |") or die "Can not open the file $inFH_name.\n";
	}else{
		open($inFH,"<$inFH_name") or die "Can not open the file $inFH_name.\n"; 
	}
	
	while(<$inFH>){
		next if m/^#|^\s+/;
		chomp;
		my @tt= split /\t/,$_;
		$tt[1]--;
		if($tt[4] eq 'N'){ push @tt,0; $tt[9]='N' x $tt[5]; } #add + and Nseq for gaps
		if($tt[4] eq 'W'){ #add seq accupacy for contigs; count the time of ctg usage
			$tt[6]--;
			push @tt,substr($scfs{$tt[0]},$tt[1],$tt[2]-$tt[1]); 
		}
		push @{$agp{$tt[0]}},\@tt;
	}
	foreach my $scf_id (keys %agp){
		my $tt=$agp{$scf_id};
		@$tt=sort {$a->[3] <=>$b->[3]} @$tt;
	}	
	close $inFH;	
}

print STDERR "Finished agp reading.\n";

#### ################################
#### reading misjoin

#0:node_id	1:tsc_ali_len	2:ali_len_or_score	3:direction(5or3)	
#4:tsc_name	5:tsc_size	6:tsc_strand	7:ali_start_of_tsc	8:ali_end_of_tsc	9:joining_point_of_tsc 10:overhang_len_of_tsc	
#11:qsc_name	12:qsc_size	13:qsc_strand	14:ali_start_of_qsc	15:ali_end_of_qsc	16:joining_point_of_qsc	17:overhang_len_of_qsc	
#18:*tsc_breakpoint	19:*qsc_breakpoint	20:*breakpoint_setting
#21:tsc_num 22:tsc_ovt 23:tsc_space 24:tsc_N, 25:qsc_num 26:qsc_ovt 27:qsc_space 28:qsc_N 29:tsc_it 30:qsc_it	
my %breaks; #breaks{scf_id}=array of end positions of the agp items

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
		$tt[4]=~s/scf//; $tt[11]=~s/scf//; ### ca6.1
		if($tt[20]==1 or $tt[20]==3){
			push @{$breaks{$tt[4]}},$tt[29]-1 if $tt[3]==3; 
			push @{$breaks{$tt[4]}},$tt[29]-2 if $tt[3]==5; 
		}
		if($tt[20]==2 or $tt[20]==3){
			if($tt[13] eq '+'){
			push @{$breaks{$tt[11]}},$tt[30]-1 if $tt[3]==3; 
			push @{$breaks{$tt[11]}},$tt[30]-2 if $tt[3]==5;				
			}else{
			push @{$breaks{$tt[11]}},$tt[30]-1 if $tt[3]==5; 
			push @{$breaks{$tt[11]}},$tt[30]-2 if $tt[3]==3;			
			}
		}
	}
	close $inFH;
	
	foreach my $tt (keys %breaks){
		@{$breaks{$tt}} = sort { $a <=> $b } @{$breaks{$tt}};
	}
}

print STDERR "Finished misjoin/breaks reading.\n";

#### ################################
#### processing agps

## note that in the ca6 ctg.fasta, ctg name includes a leading 'ctg', but not in the agp. 
# for ctg:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:W, 5:ctg_id, 6:ctg_start, 7:ctg_end, 8:strand, 9:seq]
# for gap:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:N, 5:gap_size, 6:fragment, 7:yes, 8:split_flag, 9:seq ]
# need to convert from 1-base to 0 base
#$clk{ctg1@ctg2}=[ovt,num]; sta:B/C is not counted

{
	foreach my $scf_id (keys %breaks){
		my $tt=$breaks{$scf_id};
		my ($start,$end)=(0,0);
		my $i=0;
		for(;$i<scalar(@$tt);$i++){
			$end=$tt->[$i];
			push @{$agp{$scf_id.'_'.$i}},@{$agp{$scf_id}}[$start .. $end];
			$start=$end+1;
		}
		$end=scalar(@{$agp{$scf_id}})-1;
		push @{$agp{$scf_id.'_'.$i}}, @{$agp{$scf_id}}[$start .. $end];		
		delete $agp{$scf_id};	
	}
	
	foreach my $scf_id (keys %agp){
		my $tt=$agp{$scf_id};
		shift @$tt if $tt->[0][4] eq 'N';
		pop   @$tt if scalar(@$tt)>0 and $tt->[scalar(@$tt)-1][4] eq 'N';
		delete $agp{$scf_id} if scalar(@$tt)==0;
	}
	
	foreach my $scf_id (keys %agp){
		my $tt=$agp{$scf_id};		
		my ($start,$end)=(0,0);
		for(my $i=0;$i<scalar(@$tt);$i++){
			$end=$start+length($tt->[$i][9]);
			$tt->[$i][0]=$scf_id; $tt->[$i][1]=$start; $tt->[$i][2]=$end; $tt->[$i][3]=$i+1;
			$start=$end;
		}
	}	
}

print STDERR "Finished processing agp.\n";

#### ################################
#### output agp
#### 
{
	$"="\t"; #"
	open(my $outFH,"|gzip -c >$out_agp") or die "can not create agp |gzip -c >$out_agp.\n";
	foreach my $scf_id (keys %agp){
		my $tt=$agp{$scf_id};
		for(my $i=0;$i<scalar(@$tt);$i++){
			$tt->[$i][1]++; 
			if($tt->[$i][4] eq 'W'){
				print $outFH "@{$tt->[$i]}[0..8]\n";
			}else{
				print $outFH "@{$tt->[$i]}[0..7]\n";
			}
			$tt->[$i][1]--;
		}
	}
	close $outFH;
	$"=" "; #"
}
print STDERR "Finished outputing agp.\n";

#### ################################
#### output scaffolds 
#### 

{
	open(my $outFH,"| gzip -c >$out_scf") or die "can not create | gzip -c >$out_scf.\n";
	foreach my $scf_id (keys %agp){
		my $tt=$agp{$scf_id};
		my $seq;
		for(my $i=0;$i<scalar(@$tt);$i++){
			$seq .= $tt->[$i][9];
		}
		print $outFH ">scf$scf_id\n";
		while(length($seq)>0){
			print $outFH substr($seq,0,70,'')."\n";
		}
	}
	close $outFH;
}

print STDERR "Finished outputing scf.\n";