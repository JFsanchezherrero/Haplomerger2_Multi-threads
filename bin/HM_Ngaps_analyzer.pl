#!/usr/bin/perl  

## Version 1.40

## Introduction.


use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 1 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit; }

  This script is to find out Ngaps in an assembly and the haplotype sequences
that can used to minimize the size of Ngaps.
  
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
                
  --minFlankingAli=N - consider Ngaps flanked by gap-free and N-free
                       alignments of length >=N (N>=1,def=30bp)           
 
  --verbose - print more information        
  --Force - force to overwrite existing result files
  --Delete - delete intermediate and temp files (no effect on logfiles) 
  
  Required files: 
   ~result/hm.scaffolds,
   ~result/hm.nodes,
   ~result/zeroMinSpace.rbest.net.gz,
   species.fa.gz (original genome sequences).
   
  Output files include:
    ~result/hm.n_gaps.  
	 
USAGES

print "========== Start at "; system("date"); print "\n\n";
print "$0 $Arg_list\n";
my $timer = time();

#### set the output_field_separator to " " and autoflushing
$,=' ';
$|=1;

#### special parameters
# Store species names and result directory name
my (@Species, $out_dir); 
unless ($Arg_list =~ m/--Species\s+([^-]*)/) {die "No --Species argument or species_names found!\n"; }
unless (@Species = $1 =~ m/(\w+)/g)  {die "No species names found!\n"; }
$out_dir = "$Species[0].$Species[1].result";
print "Species included: ", @Species, "\n";

#### parse the arguments
my %para = ("--minFlankingAli"=>30, "--Force" => 0, "--Delete" => 0, "--verbose" => 0  );
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

#### temp variables
my ($temp, $temp1, $temp2);


#########################################
#### ngaps info
# @ngaps, [id]->[0:reserved, 1:scaffold_name, 2:scaffold_id, 3:scaffold_size, 4:Ngap_start, 5:Ngap_end, 6:Ngap_len,
#                7:node_id, 8:chain_id, 9:ali_len, 10:is_mirror, 11:strand_sign, 
#                12:tsc_id, 13:tsc_start, 14:tsc_end, 15:qsc_id, 16:qsc_start, 17:qsc_end,
#                18:tsc_replace_start, 19:tsc_replace_end, 20:qsc_replace_start, 21:qsc_replace_end, 22:new_Ngap_vs_old_Ngap_% ]
#
####
# @sc_ids, stores the name and the size of a (target) scaffold, ordered by size, descendent
#		$sc_ids[sc_id]->[0:name, 1:size, 2:start_of_node_id, 3:end_of_node_id, 4:start_of_gap_id, 5:end_of_gap_id]
#
####
# %nets{$node_id}->[id][0:reserve, 1:tsc_ali_start, 2:tsc_al_end, 3:qsc_ali_start, 4:qsc_ali_end]
#
####
# @node_ids[$id]->[0;reserve, 1:start_of_gap_id, 2:end_of_gap_id] ->@ngaps_by_node
#  
my (@ngaps, @sc_ids, %sc_names, @node_ids, %nets);

#########################################
#### read the filtered.rbest.net file and produce nodes' data structure
# @nodes, store each pairwise alignment(net) block, ordered by target position
# 	$nodes[.]->[0..16]=[0:node_id,
#                       1:tsc_id, 2:strand(1,-1), 3:tstart, 4:tlen, 
#	                      5:qsc_id, 6:strand(1,-1), 7:qstart, 8:qlen, 
#                       9:ali_score, 10:ali_len,   11:predecessor,  12:successor, 
#                      13:tsc_portion_ref,      14:qsc_portion_ref,
#											 15:tsc_direction_predecessor, 16:tsc_direction_successor,
#											 17:qsc_direction_predecessor, 18:tsc_direction_successor,
#                      19:deletion_visit_flag
#											 ( +1=visited during perfect path pair, +2=visited in mis-join breaking, +4=visit in path finding, 3=both 1 and 2, 5=1 and 4, 6=2 and 4, 7=1,2,4
#                      ( >=0:no deletion,  -1: manual deletion, -2:score filtered, -3:lowcases , -4: deletion due too many Ns_lowcases,
#                       -5:strand switich, -6:self_loop, -7:loop, -8:switch_loop -9:deleted_during_anti_mirror, 
#                      -10:, -11:trivial node (stringent), -12:trivial node(adj_node_count=0), 
#											 -13:trival node(count=1), -14:trivial node(count=2), -15:trivial node(count=3),-16:trivial node(count=4),-17:self-self nontrivial nodes
#                      20:chain_id, 21:net_id (internal), 22:mirror_or_not(0:no, 1:yes)
#                      23:[new_scaffold_id],
#                      24:[order_in_new_scaffold],
#											 25:tsc_N_count, 26:tsc_lowcase_count, 27:qsc_N_count, 28:qsc_lowcase_count,
#                      29:[manual_deletion]
#                      30:special for this script: pointer to mirror nodes
# 
my (@nodes, @nodes_by_net_id);

#### ############################################
#### read in @sc_ids and %sc_names and node table
open(my $scaffoldFH, "<$out_dir/hm.scaffolds") or die "Can not open $out_dir/hm.scaffolds !\n";
while(<$scaffoldFH>){
	next if m/^#|^\s/;
	chomp;
	my $tt = [ split /\t/ ];
	my $temp = shift @$tt;
	($tt->[4],$tt->[5])=(-1,-1);
	$sc_ids[$temp] = $tt;
}
close $scaffoldFH;
##
for(my $i=0;$i<scalar(@sc_ids);$i++) {
	$sc_names{"$sc_ids[$i]->[0]"} = $i;
}
##
open(my $nodeFH,"<$out_dir/hm.nodes") or die "Can not open $out_dir/hm.nodes!\n";
while(<$nodeFH>){
	next if m/^#|^\s/;
	chomp;
	my $tt = [ split /\t/ ];
	shift @$tt; #delete tsc_name
	shift @$tt; #delete qsc_name
	splice @$tt,10,1;
	splice @$tt,5,1;
	push @$tt,undef; # 30: pointer to mirror node
	$nodes[$tt->[0]] = $tt;
}
close $nodeFH;
##
@nodes_by_net_id = sort { $a->[21] <=> $b->[21] } @nodes;
for(my $i=0;$i<scalar(@nodes_by_net_id);$i+=2){
	$nodes_by_net_id[$i]->[30]=$nodes_by_net_id[$i+1];
	$nodes_by_net_id[$i+1]->[30]=$nodes_by_net_id[$i];
}
@nodes_by_net_id=();

#### ############################################
#### read genome sequence and store Ngaps
my ($line, $is_the_end, $name, $seq);
open(my $gzFH,"gunzip -c $Species[0].fa.gz |") or die "Can't create pipe gunzip -c $Species[0].fa.gz !\n";
## look for the first seq
while($line=<$gzFH>){
	next if $line =~ m/^\s+/;
	last if $line =~ m/^>/;
}	
die "Problems with $Species[0].fa.gz on line ...\n $line \n" if !defined($line);
## store the seq
$is_the_end=0;
$temp=-1;
while($is_the_end==0){
	chomp $line;
	$line =~ /^>(\S+)/;
	$name=$1;
	$seq ='';
	
	while($line=<$gzFH>){
		next if $line =~ m/^\s+/;
		last if $line =~ m/^>/;
		chomp $line;
		$seq .= $line;
	}
	$is_the_end=1 if !defined($line);
	
	####
	while($seq =~ m/(N+)/gi){
		$temp1=length($1);
		$temp2=pos($seq);
		push @ngaps,[-1, $name, $sc_names{$name}, $sc_ids[$sc_names{$name}]->[1], $temp2-$temp1, $temp2, $temp1,
									-1, -1, -1, -1, -1,   -1,-1,-1,-1,-1,-1,      -1,-1, -1,-1, 100 ]; 
	}
}
close $gzFH;
## order them
@ngaps=sort { $a->[2] <=> $b->[2] or $a->[4] <=> $b->[4] } @ngaps;
$temp = $ngaps[0]->[2]; # the first one
$sc_ids[$temp]->[4]=0;
for(my $i=1;$i<scalar(@ngaps);$i++){ 
	if($temp != $ngaps[$i]->[2]){
		$sc_ids[$temp]->[5] = $i-1;
		$temp = $ngaps[$i]->[2];
		$sc_ids[$temp]->[4] = $i;
	}
}
$temp = $ngaps[scalar(@ngaps)-1]->[2]; # the last one
$sc_ids[$temp]->[5] = scalar(@ngaps)-1; 
## append node info to @ngaps
for(my $i=0;$i<scalar(@ngaps);$i++){
	$temp=$ngaps[$i]->[2];
	for(my $k=$sc_ids[$temp]->[2];$k<=$sc_ids[$temp]->[3] and $k>-1;$k++){
		if($ngaps[$i]->[4]>=$nodes[$k]->[3] and $ngaps[$i]->[5]<=$nodes[$k]->[3]+$nodes[$k]->[4]){
			@{$ngaps[$i]}[7 .. 17] = ($k, $nodes[$k]->[20], $nodes[$k]->[10], $nodes[$k]->[22], $nodes[$k]->[6],
																$nodes[$k]->[1],$nodes[$k]->[3],$nodes[$k]->[3]+$nodes[$k]->[4],
																$nodes[$k]->[5],$nodes[$k]->[7],$nodes[$k]->[7]+$nodes[$k]->[8]);
			last;
		}
	}
}
## order @ngaps by nodes, @node_ids
for(my $i=0;$i<scalar(@nodes);$i++){
	$node_ids[$i]=[0,-1,-1];
}
#@ngaps=sort { $b->[7] <=> $a->[7] or $a->[4] <=> $b->[4] } @ngaps;
#$temp=scalar(@ngaps)-1;
#for(my $i=$temp;$i>=0 and $ngaps[$i]->[7]<0;$i--){
#	splice @ngaps,$i,1;
#}
@ngaps=sort { $a->[7] <=> $b->[7] or $a->[4] <=> $b->[4] } @ngaps;
$temp = $ngaps[0]->[7]; # the first one
$node_ids[$temp]->[1]=0;
for(my $i=1;$i<scalar(@ngaps);$i++){ 
	if($temp != $ngaps[$i]->[7]){
		$node_ids[$temp]->[2] = $i-1;
		$temp = $ngaps[$i]->[7];
		$node_ids[$temp]->[1] = $i;
	}
}
$temp = $ngaps[scalar(@ngaps)-1]->[7]; # the last one
$node_ids[$temp]->[2] = scalar(@ngaps)-1; 

#### ############################################
#### read nets and parse it into N- and gap-free portions
{
	my ($tsc_id, $qsc_id, $nid,$strand, $found_flag, $finish_flag, $tstart, $tend, $qstart, $qend);
	
	##
	my $netFH;
	open($netFH,"<$out_dir/zeroMinSpace.rbest.net") 
		or open($netFH,"gunzip -c $out_dir/zeroMinSpace.rbest.net.gz | ") or die "Can not open $out_dir/zeroMinSpace.rbest.net(.gz)!\n";
	$found_flag=0;
	$finish_flag=1;
	## extract gap-free alignments for each nets
	while(<$netFH>){
		if(m/^net\s+([-\.\w]+)\s+/){
			if($finish_flag<1){
				push @{$nets{$nid}},[$strand, $tstart, $tend, $qstart, $qend];
				$finish_flag=1;
			}				
			$tsc_id=$sc_names{$1};
		}elsif(m/^\sfill\s+(\d+)\s+(\d+)\s+([-\.\w]+)\s+([+-])\s+(\d+)\s+(\d+)\s+id\s+(\d+)/){
			if($finish_flag<1){
				push @{$nets{$nid}},[$strand, $tstart, $tend, $qstart, $qend];
				$finish_flag=1;
			}
			$qsc_id=$sc_names{$3};
			$found_flag=0;
			for(my $i=$sc_ids[$tsc_id]->[2];$i<=$sc_ids[$tsc_id]->[3] and $i>-1 and $found_flag==0;$i++){
				if($nodes[$i]->[20] == $7 and $nodes[$i]->[22]<1 and $nodes[$i]->[5]==$qsc_id and $1<=$nodes[$i]->[3] and $nodes[$i]->[4]+$nodes[$i]->[3]<=$1+$2){
					$found_flag=1;
					$nid=$i;
					$finish_flag=0;
					$strand= $4 eq '+' ? 1 : -1;
					$tstart= $nodes[$i]->[3];
					$qstart= $nodes[$i]->[7];
					$tend  = $nodes[$i]->[3]+$nodes[$i]->[4];	
					$qend  = $nodes[$i]->[7]+$nodes[$i]->[8];			
				}
			}
		}elsif($found_flag>0 and m/^\s\sgap\s+(\d+)\s+(\d+)\s+[-\.\w]+\s+[+-]\s+(\d+)\s+(\d+)/){
			if(	$strand>0 and ($1>=$tstart and $1+$2<=$tend) ){
				push @{$nets{$nid}},[$strand, $tstart, $1, $qstart, $3];
				$tstart=$1+$2;
				$qstart=$3+$4;
			}elsif(	$strand<0 and ($1>=$tstart and $1+$2<=$tend) ){
				push @{$nets{$nid}},[$strand, $tstart, $1, $3+$4, $qend];
				$tstart=$1+$2;
				$qend  =$3;				
			}
		}
	}
	if($finish_flag<1){
		push @{$nets{$nid}},[$strand, $tstart, $tend, $qstart, $qend];
	}	
	close $netFH;

	## remove N from gap-free alignments (divide alignments by Ns)
	my (@qe1,@qe2);
	foreach $nid (keys %nets){
		## processing the target side
		$temp=$nid;
		@qe1=@{$nets{$temp}};
		@qe2=();
		$strand=$qe1[0]->[0];
		for(my $i=0;$i<scalar(@qe1);$i++){
			($tstart,$tend,$qstart,$qend)=@{$qe1[$i]}[1,2,3,4];
			$finish_flag=0;
			for(my $k=$node_ids[$temp]->[1];$k<=$node_ids[$temp]->[2] and $k>-1;$k++){
				if( $ngaps[$k]->[4]>$tstart and $ngaps[$k]->[4]<$tend ){
					if($strand>0){
						push @qe2,[ $strand, $tstart, $ngaps[$k]->[4], $qstart, $qstart+($ngaps[$k]->[4]-$tstart) ];
						if($ngaps[$k]->[5]>=$tend){ $finish_flag=1; last; }
						$tstart=$ngaps[$k]->[5];
						$qstart=$qend-($tend-$tstart);
					}else{
						push @qe2,[ $strand, $tstart, $ngaps[$k]->[4], $qend-($ngaps[$k]->[4]-$tstart),$qend ];
						if($ngaps[$k]->[5]>=$tend){ $finish_flag=1; last; }
						$tstart=$ngaps[$k]->[5];
						$qend=$qstart+($tend-$tstart);					
					}
				}elsif( $ngaps[$k]->[5]>$tstart and $ngaps[$k]->[5]<$tend ){
					if($strand>0){
						$tstart=$ngaps[$k]->[5];
						$qstart=$qend-($tend-$tstart);
					}else{
						$tstart=$ngaps[$k]->[5];
						$qend=$qstart+($tend-$tstart);					
					}					
				}elsif(	$ngaps[$k]->[4]<=$tstart and $ngaps[$k]->[5]>=$tend	){
					$finish_flag=1; last; 
				}
			}
			if($finish_flag<1){
				push @qe2,[ $strand, $tstart, $tend, $qstart, $qend ];
			}
		}
		$temp1=scalar(@qe2);
		for(my $i=$temp1-1;$i>=0;$i--){
			splice @qe2,$i,1 if $qe2[$i]->[2]==$qe2[$i]->[1];
		}
		if(scalar(@qe2)<1){
			delete $nets{$nid};
			next;
		}
		
		## processing the query side
		$temp=$nodes[$nid]->[30]->[0];
		@qe1=();
		$strand=$qe2[0]->[0];
		for(my $i=0;$i<scalar(@qe2);$i++){
			($tstart,$tend,$qstart,$qend)=@{$qe2[$i]}[1,2,3,4];
			$finish_flag=0;
			for(my $k=$node_ids[$temp]->[1];$k<=$node_ids[$temp]->[2] and $k>-1;$k++){
				if( $ngaps[$k]->[4]>$qstart and $ngaps[$k]->[4]<$qend ){
					if($strand>0){
						push @qe1,[ $strand, $tstart, $tstart+($ngaps[$k]->[4]-$qstart), $qstart, $ngaps[$k]->[4] ];
						if($ngaps[$k]->[5]>=$qend){ $finish_flag=1; last; }
						$qstart=$ngaps[$k]->[5];
						$tstart=$tend-($qend-$qstart);
					}else{
						push @qe1,[ $strand, $tend-($ngaps[$k]->[4]-$qstart), $tend, $qstart, $ngaps[$k]->[4] ];
						if($ngaps[$k]->[5]>=$qend){ $finish_flag=1; last; }
						$qstart=$ngaps[$k]->[5];
						$tend=$tstart+($qend-$qstart);					
					}
				}elsif( $ngaps[$k]->[5]>$qstart and $ngaps[$k]->[5]<$qend ){
					if($strand>0){
						$qstart=$ngaps[$k]->[5];
						$tstart=$tend-($qend-$qstart);
					}else{
						$qstart=$ngaps[$k]->[5];
						$tend=$tstart+($qend-$qstart);					
					}					
				}elsif(	$ngaps[$k]->[4]<=$qstart and $ngaps[$k]->[5]>=$qend	){
					$finish_flag=1; last; 
				}
			}
			if($finish_flag<1){
				push @qe1,[ $strand, $tstart, $tend, $qstart, $qend ];
			}
		}
		$temp1=scalar(@qe1);
		for(my $i=$temp1-1;$i>=0;$i--){
			splice @qe1,$i,1 if $qe1[$i]->[2]==$qe1[$i]->[1];
		}
				
		## save the new @nets	portions
		if(scalar(@qe1)<1){
			delete $nets{$nid};
			next;
		}else{
			@{$nets{$nid}} = sort { $a->[1]<=>$b->[1] } @qe1;
		} 	
	}
	
	## create mirror nets-portions
	foreach $temp (keys %nets){
		$nid=$nodes[$temp]->[30]->[0];
		$temp1=scalar(@{$nets{$temp}})-1;
		if($nodes[$temp]->[6]>0){
			for(my $i=0;$i<=$temp1;$i++){
				push @{$nets{$nid}},[ @{$nets{$temp}->[$i]}[0, 3,4, 1,2] ];
			}
		}else{
			for(my $i=$temp1;$i>=0;$i--){
				push @{$nets{$nid}},[ @{$nets{$temp}->[$i]}[0, 3,4, 1,2] ];
			}			
		}	
	}
}

#### ############################################
#### add replacement info
{	
	my $ali_len=$para{'--minFlankingAli'}; # minimum perfect alignment (free of gaps and Ns)
	my $offset= $ali_len>5 ? 5 : $ali_len;
	my $gap_size=1; # minimum gap size for considering
	my ($nid, @qe, $strand, $ep, $sp, $tnlen, $qnlen);
	
	for(my $i=0;$i<scalar(@ngaps);$i++){
		next if $ngaps[$i]->[7]<0;
		next if ! exists $nets{$ngaps[$i]->[7]};
		$nid=$ngaps[$i]->[7];
		$strand=$ngaps[$i]->[11];
		
		#### find flanking alignments
		@qe=@{$nets{$nid}};
		($sp,$ep)=(-1,-1);
		for(my $k=0;$k<scalar(@qe);$k++){
			if($qe[$k]->[1]>=$ngaps[$i]->[5] and $qe[$k]->[2]-$qe[$k]->[1]>=$ali_len){
				$ep=$k;
				last;
			}	
		}
		for(my $k=$ep;$k>=0;$k--){
			if($qe[$k]->[2]<=$ngaps[$i]->[4] and $qe[$k]->[2]-$qe[$k]->[1]>=$ali_len){
				$sp=$k;
				last;
			}	
		}
		next if $ep<0 or $sp<0;
		
		#### calculate Ngap size
		($tnlen, $qnlen)=(0,0);
		($temp,$temp1)=($qe[$sp]->[2],$qe[$ep]->[1]);
		for(my $k=$node_ids[$nid]->[1];$k<=$node_ids[$nid]->[2];$k++){
			if($ngaps[$k]->[4]>=$temp and $ngaps[$k]->[5]<=$temp1){
				$tnlen+=$ngaps[$k]->[6];
			}elsif($ngaps[$k]->[4]>$temp1){
				last;
			}
		}
		($temp,$temp1)= $strand>0 ? ($qe[$sp]->[4],$qe[$ep]->[3]) : ($qe[$ep]->[4],$qe[$sp]->[3]);
		for(my $k=$node_ids[$nodes[$nid]->[30]->[0]]->[1];$k<=$node_ids[$nodes[$nid]->[30]->[0]]->[2];$k++){
			if($ngaps[$k]->[4]>=$temp and $ngaps[$k]->[5]<=$temp1){
				$qnlen+=$ngaps[$k]->[6];
			}elsif($ngaps[$k]->[4]>$temp1){
				last;
			}
		}
		
		#### add replacing info
		$ngaps[$i]->[22]=int(100*$qnlen/$tnlen);
		if($strand>0){
			@{$ngaps[$i]}[18,19,20,21]=($qe[$sp]->[2]-$offset,$qe[$ep]->[1]+$offset,$qe[$sp]->[4]-$offset,$qe[$ep]->[3]+$offset);
		}else{
			@{$ngaps[$i]}[18,19,20,21]=($qe[$sp]->[2]-$offset,$qe[$ep]->[1]+$offset,$qe[$ep]->[4]-$offset,$qe[$sp]->[3]+$offset);
		}
		$ngaps[$i]->[23]=$nid;
	}
}

####
#### print out the @ngaps
if(-f "$out_dir/hm.n_gaps"){
	unlink "$out_dir/hm.n_gaps.bak";
	rename "$out_dir/hm.n_gaps","$out_dir/hm.n_gaps.bak";
}
print "Output the Ngaps information to $out_dir/hm.n_gaps.\n";
open(my $ngapFH,">$out_dir/hm.n_gaps") or die "Can not create $out_dir/hm.n_gaps! \n";
print $ngapFH "#scaffold_name\tsc_id\tsc_size\tgap_start\tgap_end\tgap_len";
print $ngapFH "\tnode_id\tchain_id\tali_len\tis_mirror_or_not\tstrand_sign";
print $ngapFH "\ttsc_id\ttsc_start\ttsc_end\tqsc_id\tqsc_start\tqsc_end";
print $ngapFH "\ttsc_start_of_being_replaced_portion\ttsc_end_of_being_replaced_portion";
print $ngapFH "\tqsc_start_of_replacing_portion\tqsc_end_of_replacing_portion\ttsc_Ngap_vs_qsc_Ngap_%\n";
$,="\t";
for(my $i=0;$i<scalar(@ngaps);$i++){
	print $ngapFH @{$ngaps[$i]}[1 .. 22];
	print $ngapFH "\n";
}
$,="";
close $ngapFH;



if($para{"--Delete"}==1){
	#unlink "$out_dir/filtered.net";
	#"$out_dir/filtered.chain.gz";
	#unlink "$out_dir/filtered.rbest.net","$out_dir/filtered.rbest.chain.gz";
}

print "\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";

