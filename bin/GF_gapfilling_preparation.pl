#!/usr/bin/perl  


## notes
## The agp file's coordinates are 0-base; the cabog6.1 uses agp specification v1.1.
## Only gaps between contig are considered.

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 2 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

   This script is to output a gap4fill file to aid gap filling process.
   
   Gzip files are supported.
   
   --scf=<scf_fasta> - scf file, should be *.scf.fasta(.gz)
   --agp=<agp file> - the agp v1.1 file, should be *.agp(.gz)
   
   --trimN - trim the low quality tails of contigs; def=not trim 
   
   Outputs include (placed in the same dir of the source file):
   *.scf.ad.fasta,
   *.ad.agp,
   *.gap4fill.
    	 
USAGES


print STDERR "\n\n========== Start at " . qx(date);
print STDERR "$0 $Arg_list\n\n";
my $timer = time();

#### ################################
#### reading arguments

my ($in_scf_name);
my ($out_scf_name); 
if ($Arg_list =~ m/--scf=(\S+)/) {
	$in_scf_name=$1;
	$in_scf_name=~m/(.+)(\.fa|\.fasta)+(\.gz)*$/;
	$out_scf_name=$1.'.ad.fa.gz';
}else{
	die "no --scf arg.\n";
}

my ($in_agp_name,$out_agp_name); 
my $out_g4f_name; 
if ($Arg_list =~ m/--agp=(\S+)/) {
	$in_agp_name=$1;
	$in_agp_name=~m/(.+)\.agp(\.gz)*$/;
	$out_agp_name=$1.'.ad.agp';
	$out_g4f_name=$1.'.ad.gap4fill';
}else{
	die "no --agp arg.\n";
}

my $trimN=0;
#if ($Arg_list =~ m/--trimN/) {
#	$trimN=1;
#}

#### ################################
#### reading scfs  

my %scfs;

{
	my $fasta_name=$in_scf_name;
	my $inFH;
	if($fasta_name=~m/\.gz$/){
		open($inFH,"gunzip -c $fasta_name |") or die "Can not open gunzip -c $fasta_name |.\n";
	}else{
		open($inFH,"<$fasta_name") or die "Can not open $fasta_name.\n"; 
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
		
		$scfs{$name}=uc $seq;
	}
	close $inFH;	
}

#### ################################
#### reading and process the info from the agp file 
#### (assume coordinate data of the agp file is in accending order)

## note that in the ca6 ctg.fasta, ctg name includes a leading 'ctg', but not in the agp. 

# for ctg:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:W, 5:ctg_id, 6:ctg_start, 7:ctg_end, 8:strand, 9:seq]
# for gap:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:N, 5:gap_size, 6:fragment, 7:yes, 8:+, 9:seq]
# need to convert from 1-base to 0 base
my %agp; 

my ($rm_scf_count,$cut_seq_count,$cut_seq_len)=(0,0,0);
{
	#### reading agp
	my $agpFH;
	if($in_agp_name=~m/\.gz$/){
		open($agpFH,"gunzip -c $in_agp_name |") or die "Can not open the in agp file $in_agp_name.\n";
	}else{
		open($agpFH,"<$in_agp_name") or die "Can not open the in agp file $in_agp_name.\n"; 
	}
	
	while(<$agpFH>){
		if( m/^#|^\s/){ next; }
		chomp;
		my @tt= split /\t/,$_;
		$tt[1]--;
		if($tt[4] eq 'N'){ push @tt,'+'; $tt[9]='N' x $tt[5]; } #add + and Nseq for gaps
		if($tt[4] eq 'W'){  #add seq 
			die "the scf $tt[0] does not exist in $in_scf_name \n" if !exists($scfs{$tt[0]});
			$tt[9]=substr($scfs{$tt[0]},$tt[1],$tt[2]-$tt[1]);
		}
		push @{$agp{$tt[0]}},\@tt;
	}
	foreach my $scf_id (keys %agp){
		my $tt=$agp{$scf_id};
		@$tt=sort {$a->[3] <=>$b->[3]} @$tt;
	}	
	close $agpFH;
	
	#### reporting double used contigs 
	#foreach my $cc (keys %ctgs_count){
	#	die "ctg $cc has been used >1 time.\n" if $ctgs_count{$cc}>1;
	#}
	
	#### remove leading or ending Ngaps
	#foreach my $scf_id (keys %agp){
	#	my $agp_sc=$agp{$scf_id};
	#	while($agp_sc->[0][4] eq 'N'){ print STDERR "$agp_sc->[0][0] $agp_sc->[0][3]\n"; shift @$agp_sc; }
	#	while($agp_sc->[scalar(@$agp_sc)-1][4] eq 'N'){ print STDERR "$agp_sc->[scalar(@$agp_sc)-1][0] $agp_sc->[scalar(@$agp_sc)-1][3]\n"; pop @$agp_sc; }
	#}
	
	#### cut Ns at the ctg border
	# require that the leading and ending 42bp should not containing any gaps 
	# for ctg:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:W, 5:ctg_id, 6:ctg_start, 7:ctg_end, 8:strand, 9:seq]
	# for gap:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:N, 5:gap_size, 6:fragment, 7:yes, 8:+, 9:seq]
	if($trimN>0){
	my ($count,$seq,$leading,$ending);
	foreach my $scf_id (keys %agp){
		my $tt=$agp{$scf_id};
		$count=scalar(@$tt);
		for(my $i=0;$i<$count;$i++){
			if($tt->[$i][4] eq 'W'){
				$seq=$tt->[$i][9];
				# processing the leading seq
				if(($i==0 or $tt->[$i-1][4] eq 'N') and $seq=~m/[ATCG]{42}/){
					$leading=$-[0];
					if($leading>0){
						substr($seq,0,$leading,'');
						if($i!=0){ $tt->[$i-1][9] .= 'N' x $leading; }
						$cut_seq_count++; $cut_seq_len+=$leading;
					}
				}elsif($seq=~m/^N+/){
					$leading=$+[0];
					substr($seq,0,$leading,'');
					if($i!=0){ $tt->[$i-1][9] .= 'N' x $leading; }
					$cut_seq_count++;	$cut_seq_len+=$leading;				
				}
				
				#elsif(($i==0 or $tt->[$i-1][4] eq 'N') and $seq!~m/[ATCG]{50} and ($i==$count-1 or $tt->[$i-1][4] eq 'N')){
				#	# to remove the contig
				#	if($i-1>=0){ $tt->[$i-1][9] .= 'N' x length($seq); }
				#	if($i+1<$count){ 
				#		$tt->[$i-1][9] .= $tt->[$i+1][9]; splice(@$tt,$i,2); $count-=2; 
				#	}else{
				#		splice(@$tt,$i,1); $count-=1;
				#	}
				#	$rm_ctg_count++;
				#	next;
				#}
				
				# processing the ending seq
				$seq=reverse $seq;
				if(($i==$count-1 or $tt->[$i+1][4] eq 'N') and $seq=~m/[ATCG]{42}/){
					$ending=$-[0];
					if($ending>0){
						substr($seq,0,$ending,'');
						if($i!=$count-1){ $tt->[$i+1][9] .= 'N' x $ending; }
						$cut_seq_count++; $cut_seq_len+=$ending;
					}
				}elsif($seq=~m/^N+/){
					$ending=$+[0];
					substr($seq,0,$ending,'');
					if($i!=$count-1){ $tt->[$i+1][9] .= 'N' x $ending; }
					$cut_seq_count++;	$cut_seq_len+=$ending;				
				}
				
				# retake the processed seq
				$tt->[$i][9]=reverse $seq;
			}
		}
	}
  }

	#### remove scaffolds containing only Ns
	#foreach my $scf_id (keys %agp){
	#	my $tt=$agp{$scf_id};
	#	my $have_ctg=0;
	#	for(my $i=0;$i<scalar(@$tt);$i++){
	#		$have_ctg = 1 if $tt->[$i][4] eq 'W';
	#	}
	#	if($have_ctg<1){ delete $agp{$scf_id}; $rm_scf_count++;  }
	#}
	
	#### adjust the coordinates
	# require that the leading and ending 50bp should not containing any gaps 
	# for ctg:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:W, 5:ctg_id, 6:ctg_start, 7:ctg_end, 8:strand, 9:seq]
	# for gap:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:N, 5:gap_size, 6:fragment, 7:yes, 8:+, 9:seq]
	#foreach my $scf_id (keys %agp){
	#	my $tt=$agp{$scf_id};
	#	my ($scf_start, $scf_end, $len)=(0,0,0);
	#	for(my $i=0;$i<scalar(@$tt);$i++){
	#		$tt->[$i][3]=$i+1;
	#		$len=length($tt->[$i][9]);
	#		$tt->[$i][1]=$scf_start;
	#		$tt->[$i][2]=$scf_start+$len;
	#		$scf_start+=$len;
	#		if($tt->[$i][4] eq 'W'){
	#			$tt->[$i][6]=0;
	#			$tt->[$i][7]=$len;
	#		}else{
	#			$tt->[$i][5]=$len;
	#		}
	#	}
	#}
	
	print STDERR "deleted scf count: $rm_scf_count; cut seq count: $cut_seq_count; cut seq len: $cut_seq_len.\n";			
}

#### ################################
#### output contigs 
#### 
my $outFH;

#open($outFH,"| gzip -c >$out_ctg_name") or die "can not create | gzip -c >$out_ctg_name.\n";
#foreach my $scf_id (keys %agp){
#	my $tt=$agp{$scf_id};
#	my $seq;
#	for(my $i=0;$i<scalar(@$tt);$i++){
#		if($tt->[$i][4] eq 'W'){
#			print $outFH ">ctg$tt->[$i][5]\n";
#			$seq=$tt->[$i][9];
#			if($tt->[$i][8] eq '-'){
#				$seq=~tr/ACGT/TGCA/;
#				$seq=reverse $seq;
#			}
#			while(length($seq)>0){
#				print $outFH substr($seq,0,70,'')."\n";
#			}
#		}
#	}
#}
#close $outFH;


#### ################################
#### output agp
#### 

$"="\t"; #"
open($outFH,">$out_agp_name") or die "can not create agp $out_agp_name.\n";
foreach my $scf_id (keys %agp){
	my $tt=$agp{$scf_id};
	for(my $i=0;$i<scalar(@$tt);$i++){
		$tt->[$i][1]++; 
		if($tt->[$i][4] eq 'W'){
			$tt->[$i][6]++;
			print $outFH "@{$tt->[$i]}[0..8]\n";
			$tt->[$i][6]--;
		}else{
			print $outFH "@{$tt->[$i]}[0..7]\n";
		}
		$tt->[$i][1]--;
	}
}

close $outFH;
$"=" "; #"

#### ################################
#### output scaffolds 
#### 

open($outFH,"| gzip -c >$out_scf_name") or die "can not create | gzip -c >$out_scf_name.\n";
foreach my $scf_id (keys %agp){
	my $tt=$agp{$scf_id};
	my $seq;
	for(my $i=0;$i<scalar(@$tt);$i++){
		$seq .= $tt->[$i][9];
	}
	print $outFH ">$scf_id\n";
	while(length($seq)>0){
		print $outFH substr($seq,0,60,'')."\n";
	}
}
close $outFH;

print STDERR "finished ctg/scf/agp adjustment.\n";

#### ################################
#### creating gap4fill tables 
#### 

#### break contigs at Ngaps
# for ctg:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:W, 5:ctg_id, 6:ctg_start, 7:ctg_end, 8:strand, 9:seq]
# for gap:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:N, 5:gap_size, 6:fragment, 7:yes, 8:+, 9:seq]
# for gap:: 5:WN 
foreach my $scf_id (keys %agp){
	my $tt=$agp{$scf_id};
	my ($seq,$count,$start,$end)=('',scalar(@$tt),0,0);
	for(my $i=0;$i<$count;$i++){
		if($tt->[$i][4] eq 'W'){
			$seq=$tt->[$i][9];
			if($seq=~m/N+/){
				($start,$end)=($-[0],$+[0]);
				my $ttt1=[$tt->[$i][0],$tt->[$i][1],$tt->[$i][1]+$start,0,'W',$tt->[$i][5],$tt->[$i][6],$tt->[$i][6]+$start,$tt->[$i][8],''];
				my $ttt2=[$tt->[$i][0],$tt->[$i][1]+$start,$tt->[$i][1]+$end,0,'WN',$tt->[$i][5],$tt->[$i][6]+$start,$tt->[$i][6]+$end,$tt->[$i][8],''];
				my $ttt3=[$tt->[$i][0],$tt->[$i][1]+$end,$tt->[$i][2],0,'W',$tt->[$i][5],$tt->[$i][6]+$end,$tt->[$i][7],$tt->[$i][8],''];
				splice(@$tt,$i,1,$ttt1,$ttt2,$ttt3);
				$i+=1; $count+=2;
			}
		}
	}
}

open($outFH,">$out_g4f_name") or die "can not open gap4fill table $out_g4f_name.\n";
print $outFH "#agpv1.1; 0-based. ctg_start/end refers to the strand used in the scaffold, not the original strand.\n";
print $outFH "#for gaps between cotig, 4N, 5len, 6fragment, 7yes, 8null\n";
print $outFH "#0scf_id\t1scf_start\t2scf_end\t3count\t4type\t5ctg_id\t6ctg_start\t7ctg_end\t8ctg_strand\t9gap_size\t10gap_count\t11filling_seq\n";
$"="\t"; #"
foreach my $scf_id (keys %agp){
	my $tt=$agp{$scf_id};
	my ($count,$gap_size,$gap_count)=(0,0,0);
	for(my $i=0;$i<scalar(@$tt);$i++){
		if($tt->[$i][4] eq 'WN'){
			$gap_count++; $tt->[$i][10]=$gap_count;
			$count++; $tt->[$i][3]=$count;
			$gap_size =	$tt->[$i][7]-$tt->[$i][6];
			print $outFH "@{$tt->[$i]}[0..8]\t$gap_size\t$gap_count\t\n";
		}elsif($tt->[$i][4] eq 'N'){
			$gap_count++; $tt->[$i][10]=$gap_count;
			$count++; $tt->[$i][3]=$count;
			print $outFH "@{$tt->[$i]}[0..5]\tfragment\tyes\t\t$tt->[$i][5]\t$gap_count\t\n";
		}else{
			$count++; $tt->[$i][3]=$count;
			print $outFH "@{$tt->[$i]}[0..8]\t-1\t-1\t\n";
		}
	}
}


print STDERR "\n\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";
