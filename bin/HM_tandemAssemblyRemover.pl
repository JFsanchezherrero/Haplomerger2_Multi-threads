#!/usr/bin/perl  

## Version 1.00

## Introduction.

## notes
## 1) (solved by cut off large gaps) occasionally, (<0.2%) incorrect recognition of tandems, 
## these tandems are from one net of many/large gaps
## parameters affecting detection:
## ali_len, ali_coverage%, minimum gap size to cut, --minLen, --maxInterval 

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 1 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit; }

  This script is to remove the tandem assembly errors from the old assemlby.
A file containing a new assembly free of tandem errors is provided instead.
It requires chain file named self.rbest.chain.gz!
  
  --help or no arguments - show this message, ignore undefined inputs
  
  --Species -   mandatory, to provide species names which are wanted to be
                processed
                NOTE THAT the species name should be given
                in ORDER, because the name comes first will be used
                as the target.
                Example: --Species bbev1a bbev1b
                In this script, both names refer to the same genome assembly,
                the name come first is used as the target genome, the second
                name is used as the query genome.
                
  --runChainNet - default(=1) to run chainNet first, 
                set to 0 start from hm.tandem_assemblies 
                
  --filter=N - filter out alignment block whose ali_len low than N 
               (def=2500 bp of ali_len, minimum=2500 bp)                
  --minLen=N - minimum alignment length (on target sequence) for consideration
               (default=10000 bp, minimum=5000 bp)
  --maxInterval=N - maximum interval between tandom assembly allowed
               (default=50(%) of the alignment length; a large
                value like 10,000,000 leads to the removal of translocation-like
                mis-assemblies)
  --minCoverage=N - ignore sequences whose alignment-coverage < N (def=66%)                                   
 
  --verbose - print more information        
  --Force - force to overwrite existing result files
  --Delete - delete intermediate and temp files (no effect on logfiles) 
  
  Output files include: 
  ~result/species[0]_xt.fa.gz,
  ~result/hm.tandem_assemblies.  
	 
USAGES

print "========== Start at "; system("date"); print "\n\n";
print "$0 $Arg_list\n";
my $timer = time();

#### set the output_field_separator to " " and autoflushing
$,='';
$|=1;

#### special parameters
# Store species names and result directory name
my (@Species, $out_dir); 
unless ($Arg_list =~ m/--Species\s+([^-]*)/) {die "No --Species argument or species_names found!\n"; }
unless (@Species = $1 =~ m/(\w+)/g)  {die "No species names found!\n"; }
$out_dir = "$Species[0].$Species[1].result";
print "Species included: ", @Species, "\n";

#### parse the arguments
my %para = ("--filter" => 2500, "--minLen" => 10000, "--maxInterval" => 50, "--runChainNet"=>0, 
						"--minCoverage"=>66, "--Force" => 0, "--Delete" => 0, "--verbose" => 0  );
# parameters with values
while($Arg_list =~ m/(--\w+)=(\w+)/g){ 
	if(exists($para{$1})){
		$para{$1}=$2;
	}else{
		die "Parameter $1 is unknown!\n";
	}
}
# parameters of toggling type
while($Arg_list =~ m/(--\w+)\s+-/g){ 
	if(exists($para{$1})){
		$para{$1}=1;
	}else{
		die "Parameter $1 is unknown!\n";
	}
}
while($Arg_list =~ m/(--\w+)\s*\z/g){ 
	if(exists($para{$1})){
		$para{$1}=1;
	}else{
		die "Parameter $1 is unknown!\n";
	}
}

$para{'--minLen'} =5000 if $para{'--minLen'}<5000;
$para{'--maxInterval'} = $para{'--maxInterval'}/100;
$para{'--minCoverage'} =$para{'--minCoverage'}/100;

#### set the over-write flag and deletion flag
print "Set to OVER-WRITING mode!\n" if $para{"--Force"}>0;
print "Set to DELETING mode!\n" if $para{"--Delete"}>0;
print "Set to VERBOSE mode!\n" if $para{"--verbose"}>0;

#########################################
#### read the filtered.rbest.net file and produce nodes' data structure
# @nets, store each pairwise alignment(net) block, ordered by target position
# 	$nets[.]->[0..16]=[0:reserved,
#                       1:tsc_id, 2:strand(1,-1), 3:tstart, 4:tend, 
#	                      5:qsc_id, 6:strand(1,-1), 7:qstart, 8:qend, 
#                       9:ali_len, 10:ali_score, 11: chain_id, 12:reserved
####
# 	$gaps[.]->[0..16]=[0:reserved,
#                       1:tsc_id, 2:strand(1,-1), 3:tstart, 4:tend, 
#	                      5:qsc_id, 6:strand(1,-1), 7:qstart, 8:qend, 
####
# 	$blocks[.]->[0..16]=[0:reserved,
#                       1:tsc_id, 2:strand(1,-1), 3:tstart, 4:tend, 
#	                      5:qsc_id, 6:strand(1,-1), 7:qstart, 8:qend,
my (@nets,@gaps,@blocks,@blocks_tmp,%sc); # nets and gaps are actually stack

#########################################
#### read the filtered.rbest.net file and produce nodes' data structure
# @nets, store each pairwise alignment(net) block, ordered by target position
# 	$nets[.]->[0..16]=[0:reserved,
#                       1:tsc_id, 2:strand(1,-1), 3:tstart, 4:tend, 
#	                      5:qsc_id, 6:strand(1,-1), 7:qstart, 8:qend, 
#                       9:ali_len, 10:ali_score, 11: chain_id, 12:reserved
####
# 	$gaps[.]->[0..16]=[0:reserved,
#                       1:tsc_id, 2:strand(1,-1), 3:tstart, 4:tend, 
#	                      5:qsc_id, 6:strand(1,-1), 7:qstart, 8:qend, 
####
# 	$blocks[.]->[0..16]=[0:reserved(ali_len),
#                       1:tsc_id, 2:strand(1,-1), 3:tstart, 4:tend, 
#	                      5:qsc_id, 6:strand(1,-1), 7:qstart, 8:qend, 
####  
# $tandems[.]->[0..16]=[0:sc_id,
#                       1:1st_start, 2:1st_len 3:2nd_start, 4:2nd_len, 
#	                      5:intervals, 
#                      	6:1st_N_count, 7:1st_lowcase_count, 8:2nd_N_count, 9:2nd_lowcase_count,
#                      10:delete or not (0:not delete, 1:delete 1st_portion, 3:delete 2nd_portion)
my @tandems;

#### temp variables
my ($temp, $temp1, $temp2);
my $cmd;

####
sub find_tandems();
sub remove_tandems();

find_tandems() if $para{'--runChainNet'}>0;
remove_tandems();

if($para{'--Delete'}>0){
	print "Force to delete temp files in $out_dir/selfOnly.nets/.\n";
	unlink glob "$out_dir/selfOnly.nets/*.*";
	rmdir "$out_dir/selfOnly.nets";
}

print "\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";

##################################     subs     ##################################################

##############################
#### 
##############################
sub find_tandems(){
	
	####
	#### split net and filter them
	my @net_files;	

	if(-d "$out_dir/selfOnly.nets" and $para{'--Force'}==0){
		print "$out_dir/selfOnly.nets already exists, die!\n";
		die;
	}
	unlink glob "$out_dir/selfOnly.nets/*.*";
	rmdir "$out_dir/selfOnly.nets";
	if(! -f "$out_dir/all.rbest.chain.gz"){ print "$out_dir/all.rbest.chain.gz is not found!\n"; }
	$cmd="chainSplit $out_dir/selfOnly.nets $out_dir/all.rbest.chain.gz";
	print $cmd,"\n";
	system($cmd);	
		
	@net_files = glob "$out_dir/selfOnly.nets/*.chain";
	print "ChainNet and filter nets ...\n";

	for(my $i=0;$i<scalar(@net_files);$i++){
		$net_files[$i] =~ m/nets\/(.+)\.chain/;
		$net_files[$i] = $1;
		
		$cmd  = "chainFilter -q=$net_files[$i] $out_dir/selfOnly.nets/$net_files[$i].chain 2>/dev/null";
		$cmd .= "| chainPreNet /dev/stdin $Species[0].sizes $Species[1].sizes /dev/stdout 2>/dev/null";
		$cmd .= "| chainNet -minSpace=25 -minScore=2000 /dev/stdin $Species[0].sizes $Species[1].sizes /dev/stdout /dev/null 2>/dev/null";
		$cmd .= "| netSyntenic /dev/stdin /dev/stdout 2>/dev/null >$out_dir/selfOnly.nets/$net_files[$i].net 2>/dev/null";	
		print $cmd,"\n";
		system($cmd);
	}
	
	####
	####               
	print "Read the net files from $out_dir/selfOnly.nets and parse them ...\n";
	my ($tsc, $qsc, $step, $strand, $delete_step);
	my $block_pt=-1;
	
	foreach $tsc (@net_files){
	
	open(my $netFH, "<$out_dir/selfOnly.nets/$tsc.net") or die "Can not read $out_dir/selfOnly.nets/$tsc.net !\n";
	$_=<$netFH>;
	unless(defined($_) and $_=~m/^net/){ close $netFH; next; } # no net found!
	
	($qsc, $step, $strand, $delete_step)=($tsc, -1, 1, -1);
	
	while(<$netFH>){
		##
		if(m/^(\s+)fill\s+(\d+)\s+(\d+)\s+([-\.\w]+)\s+([+-])\s+(\d+)\s+(\d+)\s+id\s+\d+\s+score\s+\d+\s+ali\s+(\d+)/){
			$temp=(length($1)-1)/2;
			if($delete_step>=0 and $temp>=$delete_step){
				next;
			}
			if($3<$para{'--minLen'}*0.75 or $7<$para{'--minLen'}*0.75 or $8<$para{'--filter'}){
				$delete_step=$temp+1;
				next;
			}
			if($temp==$step){
				$block_pt++;
				$blocks[$block_pt]=[ @{$nets[$step]} ];
			}elsif($temp<$step){
				for(my $i=$step;$i>=$temp;$i--){
					$block_pt++;
					$blocks[$block_pt]=[ @{$nets[$i]} ];			
				}
			}
			$delete_step=-1;
			$step=$temp;
			$strand = $5 eq '-' ? -1 : 1;
			$nets[$step] = [-1, $tsc,1,$2,$2+$3, $qsc,$strand,$6,$6+$7];							  
		}elsif(m/^(\s+)gap\s+(\d+)\s+(\d+)\s+([-\.\w]+)\s+([+-])\s+(\d+)\s+(\d+)/){
			$temp=(length($1)-2)/2;
			if($delete_step>=0 and $temp>=$delete_step-1){
				next;
			}
			if($temp<$step){
				for(my $i=$step;$i>$temp;$i--){
					$block_pt++;
					$blocks[$block_pt]=[ @{$nets[$i]} ];			
				}					
			}
			if($3>=$para{'--minLen'}*0.75 or $7>=$para{'--minLen'}*0.75){ # minimum gap size use to cut the alignment, should be smaller than alignment len
				$block_pt++;
				if($5 eq '+'){
					$blocks[$block_pt]=[ @{$nets[$temp]}[0,1,2,3], $2, @{$nets[$temp]}[5,6,7], $6 ];
					$nets[$temp]->[3]=$2+$3;
					$nets[$temp]->[7]=$6+$7;
				}else{
					$blocks[$block_pt]=[ @{$nets[$temp]}[0,1,2,3], $2, @{$nets[$temp]}[5,6], $6+$7, $nets[$temp]->[8] ];
					$nets[$temp]->[3]=$2+$3;
					$nets[$temp]->[8]=$6;				
				}
			}			
			$step=$temp;
			$strand = $5 eq '-' ? -1 : 1;
			$gaps[$step] = [-1, $tsc,1,$2,$2+$3, $qsc,$strand,$6,$6+$7];
		}
	}
	if($step>-1){
		for(my $i=$step;$i>=0;$i--){
			$block_pt++;
			$blocks[$block_pt]=[ @{$nets[$i]} ];			
		}
	}
	close $netFH;
	}
	
	####
	#### %sc
	%sc=();
	@blocks = sort { $a->[1] cmp $b->[1] or $a->[3] <=> $b->[3] or $a->[4] <=> $b->[4] } @blocks;
	$temp=scalar(@blocks);
	for(my $i=$temp-1;$i>=0;$i--){
		splice(@blocks,$i,1) if $blocks[$i]->[4]-$blocks[$i]->[3]<150 or $blocks[$i]->[8]-$blocks[$i]->[7]<150;
	}
	for(my $i=0;$i<scalar(@blocks);$i++){
		$blocks[$i]->[0]=$i;
	}
	$temp=$blocks[0]->[1];
	$sc{$temp}->[0]=0;
	for(my $i=0;$i<scalar(@blocks);$i++){
		if($blocks[$i]->[1] ne $temp){
			$sc{$temp}->[1]=$i-1;
			$temp=$blocks[$i]->[1];
			$sc{$temp}->[0]=$i;
		}
	}
	$sc{$temp}->[1]=scalar(@blocks)-1;

	####
	#### output info
	$,="\t";
	print "\nShowing original blocks ... \n";
	for(my $i=0;$i<scalar(@blocks);$i++){
		print @{$blocks[$i]},"\n"; 
	}
	print "\n\n";
	$,='';
		
	####
	####
	my @cuts;
	$temp2=-1;
	foreach $temp (keys %sc){
		@cuts=();
		$temp1=-1;
		for(my $i=$sc{$temp}->[0];$i<=$sc{$temp}->[1];$i++){
			$temp1++;
			($cuts[$temp1]->[1],$cuts[$temp1]->[0])=($blocks[$i]->[3],$blocks[$i]->[7]);
			$temp1++;
			($cuts[$temp1]->[1],$cuts[$temp1]->[0])=($blocks[$i]->[4],$blocks[$i]->[8]);
		}
		@cuts = sort {$a->[0] <=> $b->[0]} @cuts;
		
		my ($ucut,$dcut);
		for(my $i=$sc{$temp}->[0];$i<=$sc{$temp}->[1];$i++){
			for(my $j=0;$j<scalar(@cuts);$j++){
				if( ($cuts[$j]->[0]>=$blocks[$i]->[3] and $cuts[$j]->[0]<=$blocks[$i]->[4])
						 and ($cuts[$j]->[1]>=$blocks[$i]->[7] and $cuts[$j]->[1]<=$blocks[$i]->[8]) ){
					##
					$ucut=$cuts[$j]->[0];
					$dcut=$cuts[$j]->[1];
					##
					$temp2++;
					$blocks_tmp[$temp2]=[-1, @{$blocks[$i]}[1 .. 3],	$ucut, @{$blocks[$i]}[5 .. 7], $dcut ];
					($blocks[$i]->[3],$blocks[$i]->[7])=($ucut,$dcut);			
				}
			}
			$temp2++;
			$blocks_tmp[$temp2]=[-1, @{$blocks[$i]}[1 .. 4],	@{$blocks[$i]}[5 .. 8] ];
		}	
	}
	
	####
	#### filtering out small blocks, and fill %sc
	@blocks = @blocks_tmp;
	@blocks_tmp=();
	$temp=scalar(@blocks);
	for(my $i=$temp-1;$i>=0;$i--){
		splice(@blocks,$i,1) if $blocks[$i]->[4]-$blocks[$i]->[3]<=0 or $blocks[$i]->[8]-$blocks[$i]->[7]<=0;
	}
	@blocks = sort { $a->[1] cmp $b->[1] or $a->[3] <=> $b->[3] or $a->[4] <=> $b->[4] } @blocks;
	for(my $i=0;$i<scalar(@blocks);$i++){
		$blocks[$i]->[0]=0; # 0:not visit -1:visit
	}
	%sc=();
	$temp=$blocks[0]->[1];
	$sc{$temp}->[0]=0;
	for(my $i=0;$i<scalar(@blocks);$i++){
		if($blocks[$i]->[1] ne $temp){
			$sc{$temp}->[1]=$i-1;
			$temp=$blocks[$i]->[1];
			$sc{$temp}->[0]=$i;
		}
	}
	$sc{$temp}->[1]=scalar(@blocks)-1;
	
	####
	#### output info
	$,="\t";
	print "\nShowing blocks after cutting ... \n";
	for(my $i=0;$i<scalar(@blocks);$i++){
		print @{$blocks[$i]},"\n"; 
	}
	$,='';
	
	my ($cn,$nn, $ct5,$ct3,$ctl,  $cq5,$cq3,$cql,  $nt5,$nt3,$ntl,  $nq5,$nq3,$nql );
	foreach my $s (keys(%sc)){
		for(my $i=$sc{$s}->[0];$i<=$sc{$s}->[1]-1;$i++){
			$cn=$blocks[$i];
			for(my $k=$i+1;$k<=$sc{$s}->[1];$k++){
				$nn=$blocks[$k];
				##
				($ct5,$ct3,$ctl, $cq5,$cq3,$cql)=($cn->[3],$cn->[4],$cn->[4]-$cn->[3],   $cn->[7],$cn->[8],$cn->[8]-$cn->[7]);
				($nt5,$nt3,$ntl, $nq5,$nq3,$nql)=($nn->[3],$nn->[4],$nn->[4]-$nn->[3],   $nn->[7],$nn->[8],$nn->[8]-$nn->[7]);
				##
				if( $ct5>$nq5-$nql*0.1 and $ct5<$nq5+$nql*0.1 and $cq5>$nt5-$ntl*0.1 and $cq5<$nt5+$ntl*0.1
						and $nt5>$cq5-$cql*0.1 and $nt5<$cq5+$cql*0.1 and $nq5>$ct5-$ctl*0.1 and $nq5<$ct5+$ctl*0.1 ){
					my ($xx,$yy);
					my ($xend,$yend)=($i,$k);
					my ($tcov,$qcov);
					my $found_flag=0;
					for($xx=$k-1;$xx>=$i and $found_flag==0;$xx--){
						for($yy=$k;$yy<=$sc{$s}->[1];$yy++){
							if( $blocks[$xx]->[4]>$blocks[$yy]->[8]-($blocks[$yy]->[8]-$blocks[$yy]->[7])*0.1 and $blocks[$xx]->[4]<$blocks[$yy]->[8]+($blocks[$yy]->[8]-$blocks[$yy]->[7])*0.1
									and $blocks[$xx]->[8]>$blocks[$yy]->[4]-($blocks[$yy]->[4]-$blocks[$yy]->[3])*0.1 and $blocks[$xx]->[8]<$blocks[$yy]->[4]+($blocks[$yy]->[4]-$blocks[$yy]->[3])*0.1
									and $blocks[$yy]->[4]>$blocks[$xx]->[8]-($blocks[$xx]->[8]-$blocks[$xx]->[7])*0.1 and $blocks[$yy]->[4]<$blocks[$xx]->[8]+($blocks[$xx]->[8]-$blocks[$xx]->[7])*0.1
									and $blocks[$yy]->[8]>$blocks[$xx]->[4]-($blocks[$xx]->[4]-$blocks[$xx]->[3])*0.1 and $blocks[$yy]->[8]<$blocks[$xx]->[4]+($blocks[$xx]->[4]-$blocks[$xx]->[3])*0.1 
									and ($i==$xx or $blocks[$xx]->[7]>=$blocks[$i]->[8]) 
									and ($k==$yy or $blocks[$yy]->[7]>=$blocks[$k]->[8]) ){
									
									# tesing coverage and contiguity
									$found_flag=1;
									if($i!=$xx){
										$tcov=$blocks[$i]->[4]-$blocks[$i]->[3]+$blocks[$xx]->[4]-$blocks[$xx]->[3];
										$qcov=$blocks[$i]->[8]-$blocks[$i]->[7]+$blocks[$xx]->[8]-$blocks[$xx]->[7];
										$temp=$i;
										for(my $t1=$i+1;$t1<=$xx-1;$t1++){
											if($blocks[$t1]->[7]>=$blocks[$temp]->[8] and $blocks[$t1]->[8]<=$blocks[$xx]->[7]){
												$tcov+=$blocks[$t1]->[4]-$blocks[$t1]->[3];
												$qcov+=$blocks[$t1]->[8]-$blocks[$t1]->[7];
												$temp=$t1;
											}  
										}
										$found_flag=0 if $tcov/($blocks[$xx]->[4]-$blocks[$i]->[3])<$para{'--minCoverage'} and $qcov/($blocks[$xx]->[8]-$blocks[$i]->[7])<$para{'--minCoverage'};
	
										$tcov=$blocks[$k]->[4]-$blocks[$k]->[3]+$blocks[$yy]->[4]-$blocks[$yy]->[3];
										$qcov=$blocks[$k]->[8]-$blocks[$k]->[7]+$blocks[$yy]->[8]-$blocks[$yy]->[7];										
										$temp=$k;
										for(my $t1=$k+1;$t1<=$yy-1 and $found_flag>0;$t1++){
											if($blocks[$t1]->[7]>=$blocks[$temp]->[8] and $blocks[$t1]->[8]<=$blocks[$yy]->[7]){
												$tcov+=$blocks[$t1]->[4]-$blocks[$t1]->[3];
												$qcov+=$blocks[$t1]->[8]-$blocks[$t1]->[7];
												$temp=$t1;
											}  
										}
										$found_flag=0 if $tcov/($blocks[$yy]->[4]-$blocks[$k]->[3])<$para{'--minCoverage'} and $qcov/($blocks[$yy]->[8]-$blocks[$k]->[7])<$para{'--minCoverage'};									
									}
									
									if($found_flag>0){
										($ct3,$cq3)=($blocks[$xx]->[4],$blocks[$xx]->[8]);
										($nt3,$nq3)=($blocks[$yy]->[4],$blocks[$yy]->[8]);
										($xend,$yend)=($xx,$yy);
										last;
									}
							}
						}
					}
					if($found_flag>0){
						($ctl,$cql,$ntl,$nql)=($ct3-$ct5,$cq3-$cq5,$nt3-$nt5,$nq3-$nq5);
						if( $nt5-$ct3<$ctl*$para{'--maxInterval'} and $nt5-$ct3<$cql*$para{'--maxInterval'}
								and $nt5-$ct3<$ntl*$para{'--maxInterval'} and $nt5-$ct3<$nql*$para{'--maxInterval'} ){
							push @tandems,[$s, $ct5,$ctl,$nt5,$ntl, $nt5-$ct3, -1,-1,-1,-1, 0];
							$i=$yend;
							last; # next round from [$i]
						}
					}
				}
			}
		}
	}
	
	####
	####
	%sc=();
	@tandems = sort { $a->[0] cmp $b->[0] or $a->[1] <=> $b->[1] } @tandems;
	$temp=$tandems[0]->[0];
	$sc{$temp}->[0]=0;
	for(my $i=0;$i<scalar(@tandems);$i++){
		$temp1=$tandems[$i]->[2]-int($tandems[$i]->[2]/3)*3;
		$tandems[$i]->[2] = $temp1==1 ? $tandems[$i]->[2]-1 : $temp1==2 ? $tandems[$i]->[2]-2 : $tandems[$i]->[2];
		$temp1=$tandems[$i]->[4]-int($tandems[$i]->[4]/3)*3;
		$tandems[$i]->[4] = $temp1==1 ? $tandems[$i]->[4]-1 : $temp1==2 ? $tandems[$i]->[4]-2 : $tandems[$i]->[4];	
		if($tandems[$i]->[0] ne $temp){
			$sc{$temp}->[1]=$i-1;
			$temp=$tandems[$i]->[0];
			$sc{$temp}->[0]=$i;
		}
	}
	$sc{$temp}->[1]=scalar(@tandems)-1;
	
}

##############################
#### 
##############################
sub remove_tandems(){
	
	####
	#### read in files
	if($para{'--runChainNet'}==0){
		@tandems=();
		open(my $tandemFH,"<$out_dir/hm.tandem_assemblies") 
			or die "Can not read $out_dir/hm.tandem_assemblies,\n to recreate this file please run the script with option --runChainNet=0 !\n";
		while(<$tandemFH>){
			next if m/^#|^\s/;
			chomp;
			push @tandems,[ split /\t/ ];
		}
		close $tandemFH;
		##
		%sc=();
		@tandems = sort { $a->[0] cmp $b->[0] or $a->[1] <=> $b->[1] } @tandems;
		$temp=$tandems[0]->[0];
		$sc{$temp}->[0]=0;
		for(my $i=0;$i<scalar(@tandems);$i++){
			if($tandems[$i]->[0] ne $temp){
				$sc{$temp}->[1]=$i-1;
				$temp=$tandems[$i]->[0];
				$sc{$temp}->[0]=$i;
			}
		}
		$sc{$temp}->[1]=scalar(@tandems)-1;	
	}
	
	####
	#### processing
	my ($line,$is_the_end, $name, $ss,  $seq,@new_seq, $id, $Ns,$lcs,  $start, $end, $len);
	
	open(my $inFH, "gunzip -c $Species[0].fa.gz |") or die "Can not read: gunzip -c $Species[0].fa.gz | !\n";
	open(my $outFH,"| gzip -c >$out_dir/$Species[0]_xt.fa.gz") or die "Can not create: | gzip -c >$out_dir/$Species[0]_xt.fa.gz !\n";
	
	#### look for the first seq
	while($line=<$inFH>){
		next if $line =~ m/^\s+/;
		last if $line =~ m/^>/;
	}	
	die "Problems with $Species[0].fa.gz on line ... $line \n" if !defined($line);
	
	#### store the seq
	$is_the_end=0;
	while($is_the_end==0){
		chomp $line;
		$line =~ /^>(\S+)/;
		$name=$1;
		$seq ='';
		
		while($line=<$inFH>){
			next if $line =~ m/^\s+/;
			last if $line =~ m/^>/;
			chomp $line;
			$seq .= $line;
		}
		$is_the_end=1 if !defined($line);
		
		####
		@new_seq=(); $new_seq[0]='';
		if(! exists($sc{$name})){
			$new_seq[0]=$seq; 
		}else{
			for(my $i=$sc{$name}->[0];$i<=$sc{$name}->[1];$i++){
				####
				($Ns,$lcs)=(0,0);
				$ss=substr($seq,$tandems[$i]->[1],$tandems[$i]->[2]);
				while($ss=~/(.)/g){
					$Ns++  if $1 eq 'N' or $1 eq 'n';
					$lcs++ if $1 =~ m/[acgt]/;
				}
				($tandems[$i]->[6],$tandems[$i]->[7])=($Ns,$lcs);
				####
				($Ns,$lcs)=(0,0);
				$ss=substr($seq,$tandems[$i]->[3],$tandems[$i]->[4]);
				while($ss=~/(.)/g){
					$Ns++  if $1 eq 'N' or $1 eq 'n';
					$lcs++ if $1 =~ m/[acgt]/;
				}
				($tandems[$i]->[8],$tandems[$i]->[9])=($Ns,$lcs);	
				#### seting cut off critera
				if($tandems[$i]->[2]<$para{'--minLen'} or $tandems[$i]->[4]<$para{'--minLen'}){
					$tandems[$i]->[10]=0;
				}elsif($tandems[$i]->[6]+$tandems[$i]->[7]>$tandems[$i]->[8]+$tandems[$i]->[9]){
					$tandems[$i]->[10]=1;
				}else{
					$tandems[$i]->[10]=3;
				} 
			}
			$start=0; $end=0; $id=0;
			for(my $i=$sc{$name}->[0];$i<=$sc{$name}->[1];$i++){
				$temp=$tandems[$i]->[10];
				next if $temp<1;
				$end=$tandems[$i]->[$temp];
				$new_seq[0].=substr($seq,$start,$end-$start);
				$id++;
				$new_seq[$id]=substr($seq,$end,$tandems[$i]->[$temp+1]);
				$start=$end+$tandems[$i]->[$temp+1];
			}
			$new_seq[0].=substr($seq,$start,length($seq)-$start);
		}
		####
		for(my $i=0;$i<scalar(@new_seq);$i++){
			if($i==0){
				print $outFH ">$name\n"; 
			}else{
				print $outFH ">$name"."_$i\n";
			}
			$len=length($new_seq[$i]);
			while($len>0){ print $outFH substr($new_seq[$i],0,100,''),"\n"; $len-=100;} 
		}
	}
	close $inFH;
	close $outFH;
	
	####
	#### output info
	print "\n#*****************************************************************#\n";
	print "# Finished removing of tandem assembly errors.\n\n";
	print "# The minimum portion length for deletion is set to $para{'--minLen'}.\n\n";
		
	####
	#### output info
	if(-f "$out_dir/hm.tandem_assemblies"){
		unlink "$out_dir/hm.tandem_assemblies.bak";
		rename "$out_dir/hm.tandem_assemblies","$out_dir/hm.tandem_assemblies.bak";
	}
	open(my $tandemFH,">$out_dir/hm.tandem_assemblies") or die "Can not create $out_dir/hm.tandem_assemblies!\n";
	
	print $tandemFH "# Finished removing of tandem assembly errors.\n\n";
	print $tandemFH "# The minimum portion length for deletion is set to $para{'--minLen'},\n\n";
	print $tandemFH "# This table shows the deleted and retained scaffold portions.\n\n";
	print $tandemFH "# deletion_flag==1 means that the first portion has been deleted.\n";
	print $tandemFH "# deletion_flag==3 means that the second portion has been deleted.\n"; 
	print $tandemFH "# deletion_flag==0 or -2 means that no portion has been deleted.\n\n";	
	print $tandemFH "#sc_id\t1st_start\t1st_len\t2nd_start\t2nd_len\tinterval";
	print $tandemFH "\t1st_Ns\t1st_LCs\t2nd_Ns\t2nd_LCs\tdeletion_flag\n";
	$,="\t";
	for(my $i=0;$i<scalar(@tandems);$i++){
		print $tandemFH @{$tandems[$i]};
		print $tandemFH "\n";
	}
	close $tandemFH;
	$,='';
}




