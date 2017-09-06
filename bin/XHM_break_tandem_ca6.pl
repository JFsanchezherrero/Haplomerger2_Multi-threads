#!/usr/bin/perl  


## notes
## to break the misjoins
## for testing

## assumption
## 1. alignment是准确的, tandem之间必定缺少overlap
## 2. contig是准确的

# bug: 

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 3 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

   This script is to break the tandem by finding the Ngap between 
   two tandem arranged contigs and then break it there.
   
   Gzip files are supported.
   
   --agp=<agp file>
   --scf=<scf file>
   --tandem=<tandem info>
   
   
   Outputs include (placed in the same dir of the source file):
   *.new.tandem.
    	 
USAGES


#### ################################
#### reading arguments

my ($in_agp,$out_agp,$out_scf,$out_tandem);
if ($Arg_list =~ m/--agp=(\S+)/) {
	$in_agp=$1;
	$in_agp=~m/(.+)\.agp(\.gz)*$/;
	$out_agp=$1.'.new.agp.gz';
	$out_scf=$1.'.new.scf.fa.gz';
	$out_tandem=$1.'.new.tandem';	
}else{
	die "no --agp arg.\n";
}

my ($in_scf);
if ($Arg_list =~ m/--scf=(\S+)/) {
	$in_scf=$1;
}else{
	die "no --scf arg.\n";
}

my ($in_tandem);
if ($Arg_list =~ m/--tandem=(\S+)/) {
	$in_tandem=$1;
}else{
	die "no --tandem arg.\n";
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
# for ctg:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:W, 5:ctg_id, 6:ctg_start, 7:ctg_end, 8:strand, 9:seq, 10:del_flag, 11:5'cut_len, 12:3'cut-len]
# for gap:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:N, 5:gap_size, 6:fragment, 7:yes, 8:null, 9:seq, 10:del_flag, 11:5'cut_len, 12:3'cut-len ]
# 10:del_flag(1 for del, 0 for nothing), 11/12:(0-for cut nothing, len for cut)
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
		if($tt[4] eq 'N'){ push @tt,0; $tt[9]='N' x $tt[5]; push @tt,(0,0,0); } #add + and Nseq for gaps
		if($tt[4] eq 'W'){ #add seq accupacy for contigs; count the time of ctg usage
			$tt[6]--;
			push @tt,(substr($scfs{$tt[0]},$tt[1],$tt[2]-$tt[1]),0,0,0); 
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
#### reading tandem

#0:sc_id	1:1st_start	2:1st_len	3:2nd_start	4:2nd_len	5:interval	6:1st_Ns	7:1st_LCs	8:2nd_Ns	9:2nd_LCs	10:deletion_flag	11:is_delete
# only process those tandem with gap or GF_ contig in between two tandem segments
my @tandems;

{
	my $inFH_name=$in_tandem;
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
		$tt[0]=~s/scf//;
		$tt[2]+=$tt[1]; $tt[4]+=$tt[3];
		push @tt,0;
		push @tandems,\@tt;
	}
	#foreach my $ts (keys %tandems){
	#	@{$tandems{$tt[0]}} = sort {$b->[4] <=> $a->[4]} @{$tandems{$tt[0]}};
	#}
	close $inFH;
}

print STDERR "Finished tandems reading.\n";

#### ################################
#### processing tandem

#0:sc_id	1:1st_start	2:1st_len	3:2nd_start	4:2nd_len	5:interval	6:1st_Ns	7:1st_LCs	8:2nd_Ns	9:2nd_LCs	10:deletion_flag	11:is_deleted
# only process those tandem with gap or GF_ contig in between two tandem segments	

## note that in the ca6 ctg.fasta, ctg name includes a leading 'ctg', but not in the agp. 
# for ctg:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:W, 5:ctg_id, 6:ctg_start, 7:ctg_end, 8:strand, 9:seq, 10:del_flag, 11:5'cut_len, 12:3'cut-len]
# for gap:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:N, 5:gap_size, 6:fragment, 7:yes, 8:null, 9:seq, 10:del_flag, 11:5'cut_len, 12:3'cut-len ]
# 10:del_flag(1 for del, 0 for nothing), 11/12:(0-for cut nothing, len for cut)

{
	my ($p1,$p2,$p3,$p4); # for four position from 5' to 3'
	foreach my $ts (@tandems){
		next if $ts->[10]<1;
		my $ta=$agp{$ts->[0]};
		($p1,$p2,$p3,$p4)=(-1,-1,-1,-1);
		for(my $i=0;$i<scalar(@$ta);$i++){
			$p1=$i if $ta->[$i][1]<=$ts->[1] and $ta->[$i][2]>=$ts->[1];
			$p2=$i if $ta->[$i][1]<=$ts->[2] and $ta->[$i][2]>=$ts->[2];	
			$p3=$i if $ta->[$i][1]<=$ts->[3] and $ta->[$i][2]>=$ts->[3];
			$p4=$i if $ta->[$i][1]<=$ts->[4] and $ta->[$i][2]>=$ts->[4];								
		}
		
		#the tandem should be deleted or not?
		for(my $i=$p2;$i<=$p3;$i++){
			$ts->[11]=1 if $ta->[$i][4] eq 'N';
		}
		next if $ts->[11]<1; # should not be deleted due to obvious overlap
		
		#delete the portion in the interval
		#for(my $i=$p2+1;$i<$p3;$i++){
		#	$ta->[$i][10]=1;
		#}
		
		#set deletion for the 5'/'3 tandems
		if($ts->[10]==1){
			$ta->[$p1][12]=$ta->[$p1][2]-$ts->[1];
			$ta->[$p3][11]=$ts->[3]-$ta->[$p3][1];
			#delete the portion in the interval
			for(my $i=$p1+1;$i<$p3;$i++){
				$ta->[$i][10]=1;
			}			
		}else{
			$ta->[$p4][11]=$ts->[4]-$ta->[$p4][1];
			$ta->[$p2][12]=$ta->[$p2][2]-$ts->[2];
			#delete the portion in the interval
			for(my $i=$p2+1;$i<$p4;$i++){
				$ta->[$i][10]=1;
			}					
		}
	}
	
	## remove tandem
	foreach my $scf_id (keys %agp){
		my $ta=$agp{$scf_id};
		my $count=scalar(@$ta);
		for(my $i=0;$i<$count;$i++){
			if($ta->[$i][10]==1){
				splice @$ta,$i,1; $i--; $count--;
			}elsif($ta->[$i][11]+$ta->[$i][12]>=$ta->[$i][2]-$ta->[$i][1]){
				splice @$ta,$i,1; $i--;  $count--;
			}elsif($ta->[$i][11]>0 or $ta->[$i][12]>0){
				substr($ta->[$i][9],0,$ta->[$i][11],'');
				substr($ta->[$i][9],length($ta->[$i][9])-$ta->[$i][12],$ta->[$i][12],'');
			}
		}
	}
	
	## adjust the agp files
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
			if($tt->[$i][4] ne 'N'){ $tt->[$i][6]=1; $tt->[$i][7]=length($tt->[$i][9]); }
			$start=$end;
		}
	}		
}

#### ################################
#### outputing tandem

{
	my $outFH_name=$out_tandem;
	open(my $outFH,">$outFH_name") or die "can not create >$outFH_name.\n";
	
	#print out the misjoin file
	print $outFH "#by XHM_break_tandem_ca6.pl\n";
	print $outFH "#0:sc_id	1:1st_start	2:1st_len	3:2nd_start	4:2nd_len	5:interval	6:1st_Ns	7:1st_LCs	8:2nd_Ns	9:2nd_LCs	10:deletion_flag	11:is_deleted\n";
	
	$"="\t"; #"
	foreach my $tt (@tandems){
		$tt->[0]="scf".$tt->[0];
		$tt->[2]=$tt->[2]-$tt->[1];
		$tt->[4]=$tt->[4]-$tt->[3];
		print $outFH "@$tt\n";
	}
	$"=" "; #"
}
print STDERR "Finished outputing tandem.\n";

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
