#!/usr/bin/perl  


## notes
## to break the misjoins
## for testing

## assumption
## 1. alignment是准确的, break只能发生在alignment外
## 2. contig的准确的

# bug: 

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 4 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

   This script is to break the misjoins.
   
   Gzip files are supported.
   
   --agp=<agp file> - agp v1.1 for the scf
   --scf=<scf file> - scaffold sequences
   --clk=<clk file> - mate pairs info
   --misjoin=<misjion info> - hm.assembly_errors
   
   NOTE that this script is NOT for ca6.1 original file, which
   means: 1) the scf_id in scf.fa file is the same as
   that for the agp file.
   
   Outputs include (placed in the same dir of the source file):
   *.new.misjoin.
    	 
USAGES


#### ################################
#### reading arguments

my ($in_agp,$out_agp,$out_scf,$out_mis);
if ($Arg_list =~ m/--agp=(\S+)/) {
	$in_agp=$1;
	$in_agp=~m/(.+)\.agp(\.gz)*$/;
	$out_agp=$1.'.new.agp.gz';
	$out_scf=$1.'.new.scf.fa.gz';
	$out_mis=$1.'.new.misjoin';	
}else{
	die "no --agp arg.\n";
}

my ($in_scf);
if ($Arg_list =~ m/--scf=(\S+)/) {
	$in_scf=$1;
}else{
	die "no --scf arg.\n";
}

my ($in_clk);
if ($Arg_list =~ m/--clk=(\S+)/) {
	$in_clk=$1;
}else{
	die "no --clk arg.\n";
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
		
		#$name=~s/^scf//; #######################
		$scfs{$name}=$seq;
	}
	close $inFH;
}

print STDERR "Finished scf fasta reading.\n";

#### ################################
#### reading agp

## note that in the ca6 ctg.fasta, ctg name includes a leading 'ctg', but not in the agp. 
# for ctg:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:W, 5:ctg_id, 6:ctg_start, 7:ctg_end, 8:strand, 9:seq, 10:split_flag(1after,0not,-1before)]
# for gap:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:N, 5:gap_size, 6:fragment, 7:yes, 8:null, 9:seq, 10:null ]
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
		if($tt[4] eq 'N'){ push @tt,0; $tt[9]='N' x $tt[5]; push @tt,0; } #add + and Nseq for gaps
		if($tt[4] eq 'W'){ #add seq accupacy for contigs; count the time of ctg usage
			$tt[6]--;
			push @tt,(substr($scfs{$tt[0]},$tt[1],$tt[2]-$tt[1]),0); 
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
#### reading clk
my %clk; #$clk{ctg1@ctg2}=[ovt,num]; sta:B/C is not counted

{
	my $inFH_name=$in_clk;
	my $inFH;
	if($inFH_name=~m/\.gz$/){
		open($inFH,"gunzip -c $inFH_name |") or die "Can not open the file $inFH_name.\n";
	}else{
		open($inFH,"<$inFH_name") or die "Can not open the file $inFH_name.\n"; 
	}

	$/="}\n";
	while(<$inFH>){
		next if $_ !~ m/^.CLK\n/; #.=}
		if($_ =~ m/co1:(\S+)\nco2:(\S+)\nori:(\S+)\novt:(\S+)\n.+num:(\S+)\nsta:(\S+)\n/s){
			next if $6 eq 'B' or $5 eq 'C'; # could be Bad, Chimera, Assembly, Unknown
			my $ori= $3;
			my $ovt= $4 eq 'O' ? 1 : 0;
			my $num= $5;
			#die "$1 e $2\n" if $1 eq $2; #testing#
			if(exists($clk{$1.'@'.$2.'@'.$3})){
				$clk{$1.'@'.$2.'@'.$3}->[0]+=$ovt;
				$clk{$1.'@'.$2.'@'.$3}->[1]+=$num;
				#$clk{$1.'@'.$2}->[3].=$5;
			}else{
				$clk{$1.'@'.$2.'@'.$3}=[$ovt,$num];
			}
			#testing# print "$1 $2 $ovt $num $5\n" if $5 ne 'U'; #testing# 
		}	
	}
	close $inFH;	
	$/="\n";	
}

print STDERR "Finished clk reading.\n";

#### ################################
#### reading misjoin

#0:node_id	1:tsc_ali_len	2:ali_len_or_score	3:direction(5or3)	
#4:tsc_name	5:tsc_size	6:tsc_strand	7:ali_start_of_tsc	8:ali_end_of_tsc	9:joining_point_of_tsc 10:overhang_len_of_tsc	
#11:qsc_name	12:qsc_size	13:qsc_strand	14:ali_start_of_qsc	15:ali_end_of_qsc	16:joining_point_of_qsc	17:overhang_len_of_qsc	
#18:*tsc_breakpoint	19:*qsc_breakpoint	20:*breakpoint_setting
#21:tsc_num 22:tsc_space_ali 23:tsc_space_hang 24:tsc_N, 25:qsc_num 26:qsc_space_ali 27:qsc_space_hang 28:qsc_N 29:tsc_it 30:qsc_it	
my @misjoin;

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
		#$tt[4]=~s/scf//; $tt[11]=~s/scf//; #################################
		push @tt,(-1,-1, -1,-1, -1,-1, -1,-1, -1,-1);
		push @misjoin,\@tt;
	}
	close $inFH;
}

print STDERR "Finished misjoin reading.\n";

#### ################################
#### processing misjoin

#0:node_id	1:tsc_ali_len	2:ali_len_or_score	3:direction(5or3)	
#4:tsc_name	5:tsc_size	6:tsc_strand	7:ali_start_of_tsc	8:ali_end_of_tsc	9:joining_point_of_tsc 10:overhang_len_of_tsc	
#11:qsc_name	12:qsc_size	13:qsc_strand	14:ali_start_of_qsc	15:ali_end_of_qsc	16:joining_point_of_qsc	17:overhang_len_of_qsc	
#18:*tsc_breakpoint	19:*qsc_breakpoint	20:*breakpoint_setting
#21:tsc_num 22:tsc_space_ali 23:tsc_space_hang 24:tsc_N, 25:qsc_num 26:qsc_space_ali 27:qsc_space_hang 28:qsc_N 29:tsc_it 30:qsc_it	

## note that in the ca6 ctg.fasta, ctg name includes a leading 'ctg', but not in the agp. 
# for ctg:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:W, 5:ctg_id, 6:ctg_start, 7:ctg_end, 8:strand, 9:seq]
# for gap:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:N, 5:gap_size, 6:fragment, 7:yes, 8:split_flag, 9:seq ]
# need to convert from 1-base to 0 base
#$clk{ctg1@ctg2}=[ovt,num]; sta:B/C is not counted

{
	my $range=50000; # 50kb range
	foreach my $tt (@misjoin){
	
		my ($tsc,$qsc, $tsc_point,$qsc_point)=($tt->[4],$tt->[11], $tt->[9], $tt->[16]);
		my ($tsc_it, $qsc_it); #the current item order
		my $direction= $tt->[3]; 
		my $strand = $tt->[13] eq '+' ? 1 : -1;
		my ($tsc_count, $qsc_count)=(scalar(@{$agp{$tsc}}), scalar(@{$agp{$qsc}})); # total items of a scaffold
		my (@ctg5,@ctg3,$ctgrr);
		my ($ovt,$num);
		
		##################

		print STDERR "\n=======$tsc vs $qsc ====\n";
		
		# the tsc ctg containing the breakpoint
		for(my $i=0;$i<$tsc_count;$i++){ # locate the tsc contig for having the break point
			next if $agp{$tsc}->[$i][4] ne 'W';
			if($agp{$tsc}->[$i][1]<=$tsc_point and $agp{$tsc}->[$i][2]>=$tsc_point){
				$tsc_it=$i;
				$tt->[29]=$i+1;
				last;
			}
		}

	# for 22:tsc_space_ali 23:tsc_space_hang 24:tsc_N	
		if($direction==3){
			$tt->[22]=$tsc_point-$agp{$tsc}->[$tsc_it][1];
			$tt->[23]=$agp{$tsc}->[$tsc_it][2]-$tsc_point;
			if(exists($agp{$tsc}->[$tsc_it-1]) and $agp{$tsc}->[$tsc_it-1][4] eq 'N'){			
				$tt->[24]=$agp{$tsc}->[$tsc_it-1][5].'-';			
			}else{
				$tt->[24]="0-";
			}
			if(exists($agp{$tsc}->[$tsc_it+1]) and $agp{$tsc}->[$tsc_it+1][4] eq 'N'){
				$tt->[24].=$agp{$tsc}->[$tsc_it+1][5];
			}else{
				$tt->[24].=0;
			}
		}else{
			$tt->[22]=$agp{$tsc}->[$tsc_it][2]-$tsc_point;			
			$tt->[23]=$tsc_point-$agp{$tsc}->[$tsc_it][1];
			if(exists($agp{$tsc}->[$tsc_it-1]) and $agp{$tsc}->[$tsc_it-1][4] eq 'N'){			
				$tt->[24]=$agp{$tsc}->[$tsc_it-1][5].'-';			
			}else{
				$tt->[24]="0-";
			}
			if(exists($agp{$tsc}->[$tsc_it+1]) and $agp{$tsc}->[$tsc_it+1][4] eq 'N'){
				$tt->[24].=$agp{$tsc}->[$tsc_it+1][5];
			}else{
				$tt->[24].=0;
			}			
		}
	
		#count mate numbers for tsc
		undef @ctg5; undef @ctg3;
		($ovt,$num)=(0,0);
		$ctgrr = $direction==3 ? $tsc_it : $tsc_it-1;
		for(my $k=0;$k<=$ctgrr;$k++){
			push @ctg5,[$agp{$tsc}->[$k][5],$agp{$tsc}->[$k][8]] if $agp{$tsc}->[$k][4] eq 'W' and $agp{$tsc}->[$k][2]>=$tsc_point-$range;
		}
		for(my $k=$ctgrr+1;$k<$tsc_count;$k++){
			push @ctg3,[$agp{$tsc}->[$k][5],$agp{$tsc}->[$k][8]] if $agp{$tsc}->[$k][4] eq 'W' and $agp{$tsc}->[$k][1]<=$tsc_point+$range;
		}
		print STDERR "target: @ctg5 \ntarget:@ctg3\n";
		for(my $i=0;$i<scalar(@ctg5);$i++){
			for(my $k=0;$k<scalar(@ctg3);$k++){
				my ($id1,$id2);
				if($ctg5[$i]->[1] eq '+' and $ctg3[$k]->[1] eq '+'){
					$id1=$ctg5[$i]->[0].'@'.$ctg3[$k]->[0].'@N';
					$id2=$ctg3[$k]->[0].'@'.$ctg5[$i]->[0].'@R';
				}elsif($ctg5[$i]->[1] eq '+' and $ctg3[$k]->[1] eq '-'){
					$id1=$ctg5[$i]->[0].'@'.$ctg3[$k]->[0].'@I';
					$id2=$ctg3[$k]->[0].'@'.$ctg5[$i]->[0].'@I';		
				}elsif($ctg5[$i]->[1] eq '-' and $ctg3[$k]->[1] eq '-'){
					$id1=$ctg5[$i]->[0].'@'.$ctg3[$k]->[0].'@R';
					$id2=$ctg3[$k]->[0].'@'.$ctg5[$i]->[0].'@N';			
				}elsif($ctg5[$i]->[1] eq '-' and $ctg3[$k]->[1] eq '+'){
					$id1=$ctg5[$i]->[0].'@'.$ctg3[$k]->[0].'@O';
					$id2=$ctg3[$k]->[0].'@'.$ctg5[$i]->[0].'@O';			
				}				
				if(exists $clk{$id1}){ 
					$num+=$clk{$id1}->[1]; 
					print STDERR "  target $tsc $tsc_point at $agp{$tsc}->[$tsc_it][5] : $ctg5[$i]->[0] - $ctg3[$k]->[0] : ".$clk{$id1}->[1]."\n"; #testing#
				}
				if(exists $clk{$id2}){ 
					$num+=$clk{$id2}->[1];
					print STDERR "  target $tsc $tsc_point at $agp{$tsc}->[$tsc_it][5] : $ctg3[$k]->[0] - $ctg5[$i]->[0] : ".$clk{$id2}->[1]."\n"; #testing#
				}
			}
		}
		$tt->[21]=$num;

		#####################
		
		# the qsc ctg containing the breakpoint	
		for(my $i=0;$i<$qsc_count;$i++){ # locate the qsc contig for having the break point
			next if $agp{$qsc}->[$i][4] ne 'W';
			if($agp{$qsc}->[$i][1]<=$qsc_point and $agp{$qsc}->[$i][2]>=$qsc_point){
				$qsc_it=$i;
				$tt->[30]=$i+1;
				last;
			}
		}
	
		# for 24:qsc_space 25:qsc_N	
		if(($direction==3 and $strand==1) or ($direction==5 and $strand==-1)){
			$tt->[26]=$qsc_point-$agp{$qsc}->[$qsc_it][1];			
			$tt->[27]=$agp{$qsc}->[$qsc_it][2]-$qsc_point;
			if(exists($agp{$qsc}->[$qsc_it-1]) and $agp{$qsc}->[$qsc_it-1][4] eq 'N'){			
				$tt->[28]=$agp{$qsc}->[$qsc_it-1][5].'-';			
			}else{
				$tt->[28]="0-";
			}
			if(exists($agp{$qsc}->[$qsc_it+1]) and $agp{$qsc}->[$qsc_it+1][4] eq 'N'){
				$tt->[28].=$agp{$qsc}->[$qsc_it+1][5];
			}else{
				$tt->[28].=0;
			}
		}else{
			$tt->[26]=$agp{$qsc}->[$qsc_it][2]-$qsc_point;			
			$tt->[27]=$qsc_point-$agp{$qsc}->[$qsc_it][1];
			if(exists($agp{$qsc}->[$qsc_it+1]) and $agp{$qsc}->[$qsc_it+1][4] eq 'N'){			
				$tt->[28]=$agp{$qsc}->[$qsc_it+1][5].'-';			
			}else{
				$tt->[28]="0-";
			}
			if(exists($agp{$qsc}->[$qsc_it-1]) and $agp{$qsc}->[$qsc_it-1][4] eq 'N'){
				$tt->[28].=$agp{$qsc}->[$qsc_it-1][5];
			}else{
				$tt->[28].=0;
			}				
		}
	
		#count mate numbers for qsc
		undef @ctg5; undef @ctg3;
		($ovt,$num)=(0,0);
		$ctgrr = (($direction==3 and $strand==1) or ($direction==5 and $strand==-1))  ? $qsc_it : $qsc_it-1; 
		for(my $k=0;$k<=$ctgrr;$k++){
			push @ctg5,[$agp{$qsc}->[$k][5],$agp{$qsc}->[$k][8]] if $agp{$qsc}->[$k][4] eq 'W' and $agp{$qsc}->[$k][2]>=$qsc_point-$range;
		}
		for(my $k=$ctgrr+1;$k<$qsc_count;$k++){
			push @ctg3,[$agp{$qsc}->[$k][5],$agp{$qsc}->[$k][8]] if $agp{$qsc}->[$k][4] eq 'W' and $agp{$qsc}->[$k][1]<=$qsc_point+$range;
		}
		my @ttt=(@ctg5,@ctg3);
		if($strand>0){  #####
			print STDERR "query: @ctg5 \nquery:@ctg3\n";
		}else{
			@ctg5=reverse @ctg5; @ctg3 =reverse @ctg3;
			print STDERR "query: @ctg3 \nquery:@ctg5\n";	
		} ##########
		for(my $i=0;$i<scalar(@ctg5);$i++){
			for(my $k=0;$k<scalar(@ctg3);$k++){
				my ($id1,$id2);
				if($ctg5[$i]->[1] eq '+' and $ctg3[$k]->[1] eq '+'){
					$id1=$ctg5[$i]->[0].'@'.$ctg3[$k]->[0].'@N';
					$id2=$ctg3[$k]->[0].'@'.$ctg5[$i]->[0].'@R';
				}elsif($ctg5[$i]->[1] eq '+' and $ctg3[$k]->[1] eq '-'){
					$id1=$ctg5[$i]->[0].'@'.$ctg3[$k]->[0].'@I';
					$id2=$ctg3[$k]->[0].'@'.$ctg5[$i]->[0].'@I';		
				}elsif($ctg5[$i]->[1] eq '-' and $ctg3[$k]->[1] eq '-'){
					$id1=$ctg5[$i]->[0].'@'.$ctg3[$k]->[0].'@R';
					$id2=$ctg3[$k]->[0].'@'.$ctg5[$i]->[0].'@N';			
				}elsif($ctg5[$i]->[1] eq '-' and $ctg3[$k]->[1] eq '+'){
					$id1=$ctg5[$i]->[0].'@'.$ctg3[$k]->[0].'@O';
					$id2=$ctg3[$k]->[0].'@'.$ctg5[$i]->[0].'@O';			
				}						
				if(exists $clk{$id1}){ 
					$num+=$clk{$id1}->[1];
					print STDERR "  query $qsc $qsc_point at $agp{$qsc}->[$qsc_it][5] : $ctg5[$i]->[0] - $ctg3[$k]->[0] : ".$clk{$id1}->[1]."\n"; #testing#
				}
				if(exists $clk{$id2}){
					$num+=$clk{$id2}->[1];
					print STDERR "  query $qsc $qsc_point at $agp{$qsc}->[$qsc_it][5] : $ctg3[$k]->[0] - $ctg5[$i]->[0] : ".$clk{$id2}->[1]."\n"; #testing#				
				}
			}
		}
		$tt->[25]=$num;	
		
		##################	决定对tsc还是qsc进行切断 ################	
		if($tt->[21]>2*$tt->[25] and $tt->[27]-$tt->[23]<1000){
			$tt->[20]=2;
		}elsif($tt->[25]>2*$tt->[21] and $tt->[23]-$tt->[27]<1000){
			$tt->[20]=1;
		}elsif($tt->[21]>$tt->[25] and $tt->[27]==0 and $tt->[23]>5000){
			$tt->[20]=2;
		}elsif($tt->[25]>$tt->[21] and $tt->[23]==0 and $tt->[27]>5000){
			$tt->[20]=1;
		}else{
			$tt->[20]=3;
		}
								
	}
}

#### ################################
#### outputing misjoin

{
	my $outFH_name=$out_mis;
	open(my $outFH,">$outFH_name") or die "can not create >$outFH_name.\n";
	
	#print out the misjoin file
	print $outFH "#by XHM_break_misjoin.pl\n";
	print $outFH "#0:node_id	1:tsc_ali_len	2:ali_len_or_score	3:direction(5or3)	";
	print $outFH "4:tsc_name	5:tsc_size	6:tsc_strand	7:ali_start_of_tsc	8:ali_end_of_tsc	9:joining_point_of_tsc	10:overhang_len_of_tsc	";
	print $outFH "11:qsc_name	12:qsc_size	13:qsc_strand	14:ali_start_of_qsc	15:ali_end_of_qsc	16:joining_point_of_qsc	17:overhang_len_of_qsc	";
	print $outFH "18:*tsc_breakpoint	19:*qsc_breakpoint	20:*breakpoint_setting	";
	print $outFH "21:tsc_num	22:tsc_space_ali	23:tsc_space_hang	24:tsc_N	25:qsc_num	26:qsc_space_ali	27:qsc_space_hang	28:qsc_N	29:tsc_it	30:qsc_it\n";
	
	$"="\t"; #"
	foreach my $tt (@misjoin){
		#$tt->[4]="scf".$tt->[4]; $tt->[11]="scf".$tt->[11]; ############# old
		print $outFH "@$tt\n";
	}
	$"=" "; #"
}
