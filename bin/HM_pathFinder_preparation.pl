#!/usr/bin/perl  

## Version 2.10

## Introduction.


use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 1 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit; }

  This script is to produce several tables that can be used to facilitate
the selection of represetative haplotypes and the production of a mosaic 
haplome from a diploid genome assembly.
  
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
                
  --scoreScheme=<ali|score> - choose alignment_length or alignment_score
                as the scoring system (default=score) 
  --filter=N - filter out alignment block (nodes) low than N 
               (default=200 for ali, 20000 for score)
               
  --noNsLCs - skip counting Ns and low cases (soft-masked) for each alignment,
              this saving 90-99% running time of the script (default=counting)             
 
  --verbose - print more information        
  --Force - force to overwrite existing result files
  --Delete - delete intermediate and temp files (no effect on logfiles) 
  
  Required files:
    ~result/zeroMinSpace.rbest.net.gz.
  
  Output files include:
    ~result/hm.scaffolds, 
    ~result/hm.nodes,
    ~result/hm.sc_portions.  
	 
USAGES

print "\n\n========== Start at "; system("date");
print "$0 $Arg_list\n\n";
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
my %para = ("--scoreScheme" => "score",  "--filter" => 0, "--noNsLCs" => 0,
						"--Force" => 0, "--Delete" => 0, "--verbose" => 0, "--antiEmbed" => 0  );
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

#### set score scheme
my $scheme=$para{"--scoreScheme"};
die "--scoreScheme=$scheme is not legal.\n" if $scheme ne "ali" and $scheme ne "score";
my $filter = $para{"--filter"}!=0 ? $para{"--filter"} : $scheme eq "score" ? 20000 : 200;
print "Set the scoring scheme to $scheme, and set the filter score/ali_len to $filter .\n";  

#### set the over-write flag and deletion flag
print "Set to OVER-WRITING mode!\n" if $para{"--Force"}>0;
print "Set to DELETING mode!\n" if $para{"--Delete"}>0;
print "Set to VERBOSE mode!\n" if $para{"--verbose"}>0;

#### temp variables
my ($temp, $temp1, $temp2, $tt, $ttt);


#########################################
#### scaffold info reading
# @tsc_ids, stores the name and the size of a (target) scaffold, ordered by size, descendent
#		$tsc_ids[tsc_id]->[name, size, start_of_node_id, end_of_node_id, new_scaffold_flag]
# @qsc_ids, stores the name and the size of a (query) scaffold, ordered by size, descendent
#		$qsc_ids[qsc_id]->[name, size, start_of_node_by_query_id, end_of_node_by_query_id, new_scaffold_flag]
# @sc_names, store the name(hash key) and the id (hash value) of a scaffold
#
# note that scaffold ids of all data structures use the same system !
my (@tsc_ids, @qsc_ids, %sc_names);

open(my $sizeFH,"<$Species[0].sizes") or die "Can not open file $Species[0].sizes !\n";
while(<$sizeFH>){
	next unless m/([-\.\w]+)\s+(\d+)/;
	push @tsc_ids,[$1,$2,-1,-1, 0]; #[name, size, start_of_node_id, end_of_node_id, tsc_id] #a node refers to an alignment/net block
	push @qsc_ids,[$1,$2,-1,-1, 0]; #[name, size, start_of_node_id, end_of_node_id, qsc_id] #a node refers to an alignment/net block	
}
close $sizeFH;

@tsc_ids = sort { $b->[1] <=> $a->[1] } @tsc_ids;
@qsc_ids = sort { $b->[1] <=> $a->[1] } @qsc_ids;
for(my $i=0;$i<scalar(@tsc_ids);$i++) {
	$sc_names{"$tsc_ids[$i]->[0]"} = $i;
}
my $sc_count = scalar(@tsc_ids);
print "Produce tsc/qsc_ids, total id are $sc_count !\n";

#########################################
#### read the filtered.rbest.net file and produce nodes' data structure
# @nodes, store each pairwise alignment(net) block, ordered by target position
# 	$nodes[.]->[0..16]=[0:node_id,
#                       1:tsc_id, 2:strand(1,-1), 3:tstart, 4:tlen, 
#	                      5:qsc_id, 6:strand(1,-1), 7:qstart, 8:qlen, 
#                       9:[ali_len, 10:ali_score,   11:predecessor,  12:successor, 
#                      13:tsc_portion_ref,      14:qsc_portion_ref,
#											 15:tsc_direction_predecessor, 16:tsc_direction_successor,
#											 17:qsc_direction_predecessor, 18:tsc_direction_successor,
#                      19:deletion_visit_flag
#											 ( +1=visited during perfect path pair, +2=visited in path finding, 3=both 1 and 2,
#                      ( >=0:no deletion,  -1: manual deletion, -2:score filtered, -3: low scored connection, -4: deletion due too many Ns_lowcases,
#                       -5:strand switich, -6:self_loop, -7:loop, -9:deleted_during_anti_mirror,
#												-10:trivial node(stringent), -11:trivial node (relaxed), -12:trivial node(adj_node_count=1, 
#											  -13:trival node(count=2), -14:trivial node(count=3), -15:trivial node(count=4), ),
#                      20:chain_id, 21:net_id (internal), 22:mirror_or_not(0:no, 1:yes)
#                      23:[A-reference_pointer to mirrored node at this script||B-new_scaffold_id],
#                      24:[A-tsc_net_start||B-order_in_new_scaffold],
#                      25:[A-tsc_net_end], 
#											 [->25]:-tsc_N_count, [->26]:tsc_lowcase_count, [->27]:qsc_N_count, [->28]:qsc_lowcase_count,
#											 [->29]:manual_deletion   ]
#
# %qsc_counts, count the nodes from the same query for the present target
#   $qsc_counts{"query_sc_id"}->[0,1..n]=[0:the number of occurrence, 1..n:node_ids from the same query scaffold]
my (@nodes, %qsc_counts);

#########################################
#### store nets and gaps
# @nets[net_id]->[0:tsc_id, 1:tsc_start, 2:tend, 3:qsc_id, 4:qstart, 5:qend, 6:gaps_id_start-1, 7:gaps_id_end]
# @gaps[gap_id]->[0:tstart, 1:tend, 2:qstart, 3:qend, 4:included_ali_len] 
my (@nets, @gaps);

####
####               
print "Read the zeroMinSpace.rbest.net.gz file and do some basic node filtering ...\n";
my ($tsc_pt, $qsc_pt, $node_pt, $net_pt)=(undef,undef,-1, -1);
my ($node_deletion_count, $score_addup_count)=(0,0);
my ($net_id, $gap_id)=(-1, -1);

my $netFH;
open($netFH, "<$out_dir/zeroMinSpace.rbest.net") 
	or open($netFH, "gunzip -c $out_dir/zeroMinSpace.rbest.net.gz | ") or die "Can not read $out_dir/zeroMinSpace.rbest.net(.gz) !\n";

while(<$netFH>){
	####
	if(m/^\s\sgap\s+(\d+)\s+(\d+)\s+[-\.\w]+\s+[+-]\s+(\d+)\s+(\d+)/){ ############# @nets, @gaps
		$gap_id++;
		$gaps[$gap_id]=[$1,$1+$2,$3,$3+$4,0];
		next;
	}
	####
	if(m/^net\s([-\.\w]+)\s/){
		$temp = $1;
		## if not the first target, delete the tops from the same query for the current target
		if (defined($tsc_pt)) { 
			my @deletions;
			foreach $temp1 (values(%qsc_counts)){
				next unless $temp1->[0] > 1;
				for(my $i=1;$i<=$temp1->[0];$i++){
					for(my $j=1;$j<=$temp1->[0];$j++){
						next unless $i != $j;
						if($nodes[$temp1->[$i]]->[7]>=$nodes[$temp1->[$j]]->[7] 
								and $nodes[$temp1->[$i]]->[7]+$nodes[$temp1->[$i]]->[8] <= $nodes[$temp1->[$j]]->[7]+$nodes[$temp1->[$j]]->[8]){
							push @deletions,$temp1->[$i];
							last;
						}
					}
				}
			}
			foreach $temp2 ( sort { $b<=>$a } @deletions ){
				splice @nodes,$temp2,1;
			}
			$node_pt -= scalar(@deletions);
			$node_deletion_count += scalar(@deletions); # runing messages
			%qsc_counts = ();			
		}
		## initialize a new target scaffold
		$tsc_pt = $sc_names{$temp}; # the next target scaffold
		next;
	}
	####
	if(m/^(\s+)fill\s+(\d+)\s+(\d+)\s+([-\.\w]+)\s+([+-])\s+(\d+)\s+(\d+)\s+id\s+(\d+)\s+score\s+(\d+)\s+ali\s+(\d+)\s/){
		if(length($1)<2){ ###### processing top
			$node_pt++;
			$net_pt++;
			$qsc_pt = $sc_names{"$4"};
			my $qsc_strand = $5 eq '-' ? -1 : 1;
			##
			$nodes[$node_pt] = [-1, $tsc_pt,1,$2,$3, $qsc_pt,$qsc_strand,$6,$7, $9,$10, -1,-1,-1,-1, -1,-1,-1,-1, 0,$8,$net_pt,0, undef, $2,$2+$3, $6,$6+$7];							  
			if(defined($qsc_counts{$qsc_pt})){
				$qsc_counts{$qsc_pt}->[0]++;
				push @{$qsc_counts{$qsc_pt}},$node_pt;
			}else{
				$qsc_counts{$qsc_pt} = [1,$node_pt];
			}
			
			##
			if($net_pt>0){ $nets[$net_pt-1]->[7]=$gap_id;	} ############# @nets, @gaps
			$nets[$net_pt]=[$tsc_pt, $2,$2+$3, $qsc_pt,$6,$6+$7, $gap_id,$gap_id]; ############# @nets, @gaps
			
		}elsif( $qsc_pt == $sc_names{"$4"} ###### processing non-top 
							and $nodes[$node_pt]->[7]<=$6 
							and $nodes[$node_pt]->[7]+$nodes[$node_pt]->[8]>=$6+$7 ){
			$nodes[$node_pt]->[9] += $9;
			$nodes[$node_pt]->[10] += $10;
			$gaps[$gap_id]->[4]+=$10; ############# @nets, @gaps
			$score_addup_count++; # runing messages
		}
		next;
	}
}
close $netFH;

#### finalizing the processing of the last target scaffold (delete the tops from the same query for the current target)
$nets[$node_pt]->[7]=$gap_id; ############# @nets, @gaps

if (defined($tsc_pt)) { 
	my @deletions;
	foreach $temp1 (values(%qsc_counts)){
		next unless $temp1->[0] > 1;
		for(my $i=1;$i<=$temp1->[0];$i++){
			for(my $j=1;$j<=$temp1->[0];$j++){
				next unless $i != $j;
				if($nodes[$temp1->[$i]]->[7]>=$nodes[$temp1->[$j]]->[7] 
						and $nodes[$temp1->[$i]]->[7]+$nodes[$temp1->[$i]]->[8]<=$nodes[$temp1->[$j]]->[7]+$nodes[$temp1->[$j]]->[8]){
					push @deletions,$temp1->[$i];
					last;
				}
			}
		}
	}
	foreach $temp2 ( sort { $b<=>$a } @deletions ){
		splice @nodes,$temp2,1;
	}
	$node_pt -= scalar(@deletions);
	$node_deletion_count += scalar(@deletions); # running messages
	%qsc_counts = ();			
}
print "Total nodes are ",$node_pt+1,";\ndeleted nodes are $node_deletion_count;\nadd-up score counts are $score_addup_count.\n"; 

#########################################
#### advanced node filtering and producing the final data set of nodes
####
print "To perform advanced node filtering and to produce the final data set of nodes ...\n";
print scalar(@nodes),"left before advanced filtering.\n";
print "Set filter score/ali_len to $filter for node advanced filtering ...\n";

my ($low_scored_node_count, $embeded_node_count) = (0,0); # running messages

#### filter out low scored nodes
$temp = $scheme eq "ali" ? 10 : 9;
for(my $i=scalar(@nodes)-1;$i>=0;$i--){
	#### filter out low scored nodes
	if($nodes[$i]->[$temp]<$filter){ 
		$nodes[$i]->[19]=-2; 
		$low_scored_node_count++; # running message
		next;
	}
}
for(my $i=scalar(@nodes)-1;$i>=0;$i--){
	splice @nodes,$i,1 if $nodes[$i]->[19]<0;
}

#### filter out embeded nodes
if($para{"--antiEmbed"}>0){
for(my $i=scalar(@nodes)-1;$i>=0;$i--){
	#### filter out embeded nodes
	for(my $j=0;$j<scalar(@nodes);$j++){
		next unless $i != $j;
		if( ($nodes[$j]->[1] == $nodes[$i]->[1] 
				 and $nodes[$j]->[3] <= $nodes[$i]->[3] and $nodes[$i]->[3]+$nodes[$i]->[4] <= $nodes[$j]->[3]+$nodes[$j]->[4])
				or 
				($nodes[$j]->[5] == $nodes[$i]->[5]
				 and $nodes[$j]->[7] <= $nodes[$i]->[7] and $nodes[$i]->[7]+$nodes[$i]->[8] <= $nodes[$j]->[7]+$nodes[$j]->[8])){
			$nodes[$i]->[19]=-4;
			$embeded_node_count++; # running message
			last;		
		}
	}
}
for(my $i=scalar(@nodes)-1;$i>=0;$i--){
	splice @nodes,$i,1 if $nodes[$i]->[19]<0;
}
}

print "$low_scored_node_count low-scored nodes have been filtered (deletion status value = -2)!\n"; 
print "$embeded_node_count embeded nodes have been filtered (deletion status value = -4)!\n";
print scalar(@nodes),"nodes left after advanced filtering (no mirror).\n";

#########################################
#### Recreate the complete synmetric mirror of node data	(double it)
#### Note that lastz is run in --noself mode
my (@nodes_by_score);

my $node_count=scalar(@nodes);

print "Recreate the perfect mirror alignments (nodes) ...\n";
for(my $i=0;$i<$node_count;$i++){
	my $tt1=$nodes[$i];
	my $tt2 = [-1,  @{$tt1}[5 .. 8], @{$tt1}[1 .. 4],   @{$tt1}[9 .. 23], @{$tt1}[26 .. 27], @{$tt1}[24 .. 25] ];
	if($tt2->[2]<0){ ($tt2->[2],$tt2->[6]) = ($tt2->[6],$tt2->[2]); }
	$tt2->[23]=$tt1;
	$tt2->[22]=1; 
	$tt1->[23]=$tt2;
	push @nodes,$tt2;
}

#### processing the target side
@nodes = sort { $a->[1] <=> $b->[1] or $a->[3] <=> $b->[3] or $a->[4] <=> $b->[4] or $a->[5] <=> $b->[5] or $a->[7] <=> $b->[7] } @nodes;
for(my $i=0;$i<scalar(@nodes);$i++){
	$nodes[$i]->[0]=$i;
}
#### perfect mirror for by_score
$temp = $scheme eq "ali" ? 10 : 9;
@nodes_by_score  = sort { $b->[$temp] <=> $a->[$temp] or $a->[21] <=> $b->[21] } @nodes;

#### $nodes->[19]=-9 means the node need to be deleted when producing complete mirror set of nodes
my ($cid, $csc);
my ($cstart, $cend);
for(my $i=0;$i<scalar(@nodes_by_score);$i++){
	next if $nodes_by_score[$i]->[19]<0;
	($cid, $csc)    = ($nodes_by_score[$i]->[0],$nodes_by_score[$i]->[1]);
	($cstart, $cend)= ($nodes_by_score[$i]->[24],$nodes_by_score[$i]->[25]);
	
	for(my $j=$cid-1;$j>=0;$j--){
		last if $nodes[$j]->[1] ne $csc;
		if( ($nodes[$j]->[24] >= $cstart and $nodes[$j]->[25] <= $cend) 
				or ($nodes[$j]->[24] <= $cstart and $nodes[$j]->[25] >= $cend)  ){
			$nodes[$j]->[19]=-9;
			$nodes[$j]->[23]->[19]=-9;
		}		
	}
	for(my $j=$cid+1;$j<scalar(@nodes);$j++){
		last if $nodes[$j]->[1] ne $csc;
		if( ($nodes[$j]->[24] >= $cstart and $nodes[$j]->[25] <= $cend) 
				or ($nodes[$j]->[24] <= $cstart and $nodes[$j]->[25] >= $cend)  ){
			$nodes[$j]->[19]=-9;
			$nodes[$j]->[23]->[19]=-9;	
		}	
	}	
}

#### verbose::: printing out @nodes
if($para{"--verbose"}>=1){
	$,="\t";
	print "\nThis new node table (mirror) created right before the filtering of overlapping nodes. \n";
	print "#tsc_name\tqsc_name\tnode_id\ttsc_id\tstrand\ttstart\ttend\tqsc_id\tstrand\tqstart\tqend\tscore\tali_len";
	print "\tpredecessor\tsuccessor\ttsc_portion_ref\tqsc_portion_ref\ttsc_direction_predecessor\ttsc_direction_successor\tqsc_direction_predecessor\tqsc_direction_successor";
	print "\tdeletion_visit_flag\tchain_id\tnet_id\tis_mirror\tpointer_to_mirrored_node\ttsc_start_position\ttsc_end_position\tqsc_start_position\tqsc_end_position\n";
	for(my $i=0;$i<scalar(@nodes);$i++){
		print $tsc_ids[$nodes[$i]->[1]]->[0],$qsc_ids[$nodes[$i]->[5]]->[0],@{$nodes[$i]}[0 .. 27];
		print "\n";
	}				
	print "\n\n";
	$,=" ";
}

#### remove the overlapped portion of two alignments
# REMINDER: this part needs to be simplified and improved
for(my $i=0;$i<scalar(@nodes_by_score);$i++){
	next if $nodes_by_score[$i]->[19]<0;
	($cid, $csc)    = ($nodes_by_score[$i]->[0],$nodes_by_score[$i]->[1]);
	($cstart, $cend)= ($nodes_by_score[$i]->[24],$nodes_by_score[$i]->[25]);

	my ($is_succeed);
	my ($cn); 
	my ($sign, $net_id);
	my ($diff);
	
	######## process the 5-terminal one ########
	for(my $j=$cid-1;$j>=0;$j--){
		last if $nodes[$j]->[1] ne $csc;
		next if $nodes[$j]->[19]<0;
		if($nodes[$j]->[25] > $cstart  ){
			$cn=$nodes[$j];
			($sign,$net_id)=($cn->[6],$cn->[21]);
		#### not mirror
		if($cn->[22]==0){
			$is_succeed=0;
			for(my $i=$nets[$net_id]->[7];$i>$nets[$net_id]->[6];$i--){
				next if $cn->[25] <= $gaps[$i]->[0];
			#### case 1
				if($cstart>=$gaps[$i]->[1]){
					$diff=$cstart-$gaps[$i]->[1];
					## query sign>0
					if($sign>0){
						$cn->[25] = $gaps[$i]->[1]+$diff;
						$cn->[27] = $gaps[$i]->[3]+$diff;
					## query sign<0
					}else{
						$cn->[25] = $gaps[$i]->[1]+$diff;
						$cn->[26] = $gaps[$i]->[2]-$diff;
					}
					$is_succeed=1;
					last;
			#### case 2	
				}elsif($cstart>=$gaps[$i]->[0]){
					## query sign>0
					if($sign>0){
						$cn->[25] = $gaps[$i]->[0];
						$cn->[27] = $gaps[$i]->[2];
					## query sign<0
					}else{
						$cn->[25] = $gaps[$i]->[0];
						$cn->[26] = $gaps[$i]->[3];
					}
					$is_succeed=1;
					last;
				}					
			}
			#### case 3
			if($is_succeed==0){
				$diff=$cstart-$cn->[24];
				## query sign>0
				if($sign>0){
					$cn->[25] = $cn->[24]+$diff;
					$cn->[27] = $cn->[26]+$diff;
				## query sign<0
				}else{
					$cn->[25] = $cn->[24]+$diff;
					$cn->[26] = $cn->[27]-$diff;
				}			
			}
		#### mirrored node	
		}elsif($cn->[22]==1){
			$is_succeed=0;
			if($sign>0){ ###### $sign>0
			for(my $i=$nets[$net_id]->[7];$i>$nets[$net_id]->[6];$i--){
				next if $cn->[25] <= $gaps[$i]->[2];
			#### case 1
				if($cstart>=$gaps[$i]->[3]){
					$diff=$cstart-$gaps[$i]->[3];
					$cn->[25] = $gaps[$i]->[3]+$diff;
					$cn->[27] = $gaps[$i]->[1]+$diff;
					$is_succeed=1;
					last;
			#### case 2	
				}elsif($cstart>=$gaps[$i]->[2]){
					$cn->[25] = $gaps[$i]->[2];
					$cn->[27] = $gaps[$i]->[0];
					$is_succeed=1;
					last;
				}					
			}
			#### case 3
			if($is_succeed==0){
				$diff=$cstart-$cn->[24];
				$cn->[25] = $cn->[24]+$diff;
				$cn->[27] = $cn->[26]+$diff;
			}				
			}else{ ###### $sign<0
			for(my $i=$nets[$net_id]->[6]+1;$i<=$nets[$net_id]->[7];$i++){
				next if $cn->[25] <= $gaps[$i]->[2];
			#### case 1
				if($cstart>=$gaps[$i]->[3]){
					$diff=$cstart-$gaps[$i]->[3];
					$cn->[25] = $gaps[$i]->[3]+$diff;
					$cn->[26] = $gaps[$i]->[0]-$diff;
					$is_succeed=1;
					last;
			#### case 2	
				}elsif($cstart>=$gaps[$i]->[2]){
					$cn->[25] = $gaps[$i]->[2];
					$cn->[26] = $gaps[$i]->[1];
					$is_succeed=1;
					last;
				}					
			}
			#### case 3
			if($is_succeed==0){
				$diff=$cstart-$cn->[24];
				$cn->[25] = $cn->[24]+$diff;
				$cn->[26] = $cn->[27]-$diff;
			}			
			}
		}
		#### update the mirrored node
		@{$cn->[23]}[24,25,26,27]=@{$cn}[26,27,24,25];
		}
	}
		
	######## process the 3-terminal one ########
	for(my $j=$cid+1;$j<scalar(@nodes);$j++){
		last if $nodes[$j]->[1] ne $csc;
		next if $nodes[$j]->[19]<0;
		if($nodes[$j]->[24] < $cend  ){
			$cn=$nodes[$j];
			($sign,$net_id)=($cn->[6],$cn->[21]);
		#### not mirror
		if($cn->[22]==0){
			$is_succeed=0;
			for(my $i=$nets[$net_id]->[6]+1;$i<=$nets[$net_id]->[7];$i++){
				next if $cn->[24] >= $gaps[$i]->[1];
			#### case 1
				if($cend<=$gaps[$i]->[0]){
					$diff=$gaps[$i]->[0]-$cend;
					## query sign>0
					if($sign>0){
						$cn->[24] = $gaps[$i]->[0]-$diff;
						$cn->[26] = $gaps[$i]->[2]-$diff;
					## query sign<0
					}else{
						$cn->[24] = $gaps[$i]->[0]-$diff;
						$cn->[27] = $gaps[$i]->[3]+$diff;
					}
					$is_succeed=1;
					last;
			#### case 2	
				}elsif($cend<=$gaps[$i]->[1]){
					## query sign>0
					if($sign>0){
						$cn->[24] = $gaps[$i]->[1];
						$cn->[26] = $gaps[$i]->[3];
					## query sign<0
					}else{
						$cn->[24] = $gaps[$i]->[1];
						$cn->[27] = $gaps[$i]->[2];
					}
					$is_succeed=1;
					last;
				}					
			}
			#### case 3
			if($is_succeed==0){
				$diff=$cn->[19]-$cend;
				## query sign>0
				if($sign>0){
					$cn->[24] = $cn->[25]-$diff;
					$cn->[26] = $cn->[27]-$diff;
				## query sign<0
				}else{
					$cn->[24] = $cn->[25]-$diff;
					$cn->[27] = $cn->[26]+$diff;
				}			
			}
		#### mirrored node	
		}elsif($cn->[22]==1){
			$is_succeed=0;
			if($sign>0){ ###### $sign>0
			for(my $i=$nets[$net_id]->[6]+1;$i<=$nets[$net_id]->[7];$i++){
				next if $cn->[24] >= $gaps[$i]->[3];
			#### case 1
				if($cend<=$gaps[$i]->[2]){
					$diff=$gaps[$i]->[2]-$cend;
					$cn->[24] = $gaps[$i]->[2]-$diff;
					$cn->[26] = $gaps[$i]->[0]-$diff;
					$is_succeed=1;
					last;
			#### case 2	
				}elsif($cend<=$gaps[$i]->[3]){
					$cn->[24] = $gaps[$i]->[3];
					$cn->[26] = $gaps[$i]->[1];
					$is_succeed=1;
					last;
				}					
			}
			#### case 3
			if($is_succeed==0){
				$diff=$cn->[25]-$cend;
				$cn->[24] = $cn->[25]-$diff;
				$cn->[26] = $cn->[27]-$diff;
			}				
			}else{ ###### $sign<0
			for(my $i=$nets[$net_id]->[7];$i>$nets[$net_id]->[6];$i--){
				next if $cn->[24] >= $gaps[$i]->[3];
			#### case 1
				if($cend<=$gaps[$i]->[2]){
					$diff=$gaps[$i]->[2]-$cend;
					$cn->[24] = $gaps[$i]->[2]-$diff;
					$cn->[27] = $gaps[$i]->[1]+$diff;
					$is_succeed=1;
					last;
			#### case 2	
				}elsif($cend<=$gaps[$i]->[3]){
					$cn->[24] = $gaps[$i]->[3];
					$cn->[27] = $gaps[$i]->[0];
					$is_succeed=1;
					last;
				}					
			}
			#### case 3
			if($is_succeed==0){
				$diff=$cn->[25]-$cend;
				$cn->[24] = $cn->[25]-$diff;
				$cn->[27] = $cn->[26]+$diff;
			}			
			}
		}		
		#### update the mirrored node
		@{$cn->[23]}[24,25,26,27]=@{$cn}[26,27,24,25];
		}
	}				
}	

#### adjust the score and positions
$temp = $scheme eq "ali" ? 10 : 9;
for(my $i=0;$i<scalar(@nodes);$i++){
	next if $nodes[$i]->[19]<0 or $nodes[$i]->[22]>0;
	my ($start,$end)=($nodes[$i]->[24],$nodes[$i]->[25]);
	$temp1=0;
	if($nets[$nodes[$i]->[21]]->[6]==$nets[$nodes[$i]->[21]]->[7]){
		$temp1+=$end-$start;
	}else{
		for(my $j=$nets[$nodes[$i]->[21]]->[6]+1;$j<=$nets[$nodes[$i]->[21]]->[7];$j++){
			last if $end<=$gaps[$j]->[0];
			next if $start>=$gaps[$j]->[1];
			$temp1+=$gaps[$j]->[0]-$start+$gaps[$j]->[4];
			$start=$gaps[$j]->[1];
		}
		$temp1+=$end-$start;
	}
	$nodes[$i]->[9]=int($nodes[$i]->[9]*$temp1/$nodes[$i]->[10]);
	$nodes[$i]->[10]=$temp1;
	@{$nodes[$i]}[3,4,7,8]=($nodes[$i]->[24],$nodes[$i]->[25]-$nodes[$i]->[24], $nodes[$i]->[26],$nodes[$i]->[27]-$nodes[$i]->[26]);
	$nodes[$i]->[19]=-2 if $nodes[$i]->[$temp]<$filter;
	
	@{$nodes[$i]->[23]}[3,4,7,8,9,10,19]=( @{$nodes[$i]}[7,8], @{$nodes[$i]}[3,4], @{$nodes[$i]}[9,10,19] );
}

#### verbose::: printing out @nodes
if($para{"--verbose"}>=1){
	$,="\t";
	print "\nThis new node table (mirror) created right after the filtering of overlapping nodes. \n";
	print "#tsc_name\tqsc_name\tnode_id\ttsc_id\tstrand\ttstart\ttend\tqsc_id\tstrand\tqstart\tqend\tali_score\tali_len";
	print "\tpredecessor\tsuccessor\ttsc_portion_ref\tqsc_portion_ref\ttsc_direction_predecessor\ttsc_direction_successor\tqsc_direction_predecessor\tqsc_direction_successor";
	print "\tdeletion_visit_flag\tchain_id\tnet_id\tis_mirror\tpointer_to_mirrored_node\ttsc_start_position\ttsc_end_position\tqsc_start_position\tqsc_end_position\n";
	for(my $i=0;$i<scalar(@nodes);$i++){
		print $tsc_ids[$nodes[$i]->[1]]->[0],$qsc_ids[$nodes[$i]->[5]]->[0],@{$nodes[$i]}[0 .. 27];
		print "\n";
	}				
	print "\n\n";
	$,=" ";
}

#### deleting useless nodes
for(my $i=scalar(@nodes)-1;$i>=0;$i--){
	splice @nodes,$i,1 if $nodes[$i]->[19]<0;
}

#### add the node id to the node record
@nodes = sort { $a->[1] <=> $b->[1] or $a->[3] <=> $b->[3] or $a->[4] <=> $b->[4] or $a->[5] <=> $b->[5] or $a->[7] <=> $b->[7] } @nodes;
for(my $i=0;$i<scalar(@nodes);$i++){
	$nodes[$i]->[0] = $i;
}
print scalar(@nodes),"nodes left (with mirror).\n";

#########################################
#### Producing the scaffold portion data structure	
# @tsc_portions, store each portions and the node which the portion belonged to, ordered by start position
#		$tsc_portions[tsc_id]->[portion_id]->
#			[0:portion_id, 1:tsc_id, 2:start_position, 3:len, 4:node_id, 5:node_status, 6:portion_status, 7:is_incorporated, 8:broken_due_to_assembly_errors, 9:manual_breaking] 
# @qsc_portions, store each portions and the node which the portion belonged to, ordered by start position
#		$qsc_portions[qsc_id]->[portion_id]->
#			[0:portion_id, 1:qsc_id, 2:start_position, 3:len, 4:node_id, 5:node_status, 6:portion_status, 7:is_incorporated, 8:broken_due_to_assembly_errors, 9:manual_breaking] 
#			#4:node_id=-1 means not a node portion
#			#4:node_id=id point to a node
#			#5:node_status=3 means both 1 and 2
#			#5:node_status=2 means the node has been finally visited
#			#5:node_status=1 means the node has been visited by perfect path finding
#			#5:node_status=0 means the node has not been visited
#			#5:node_status<0 means the node has been deleted(-1..-8), 
#                      ( -1: manual deletion, -2:score filtered, -3: trivial node, -4: embeded node,
#                       -5:strand switich,-6:self_loop,   -7:loop,           -8: switch_loop)
#			#6:portion_status=1 regular
#			#6:portion_status=0 means the portion have no alignment hit
#			#6:portion_status<0 means the portion have been marked as holding a breakpoint
#					#(-1:strand switch, -2:self_loop, -3:loop, -4:switch_loop,  
#					  -5:branch(reserved/not broken side), -6:branch(broken side), -7:manual breaking, -8:broken_due_to_assembly_errors)
my (@tsc_portions, @qsc_portions);

print "To produce scaffold portion data structure ...\n";

my (@nodes_by_target, @nodes_by_query);

#### produce @nodes_by_target
@nodes_by_target = @nodes;
$tsc_pt = $nodes_by_target[0]->[1]; # the first one
$tsc_ids[$tsc_pt]->[2] = 0;
for(my $i=1;$i<scalar(@nodes_by_target);$i++){ #the ones between the first and the last
	if($tsc_pt != $nodes_by_target[$i]->[1]){
		$tsc_ids[$nodes_by_target[$i-1]->[1]]->[3] = $i-1;
		$tsc_ids[$nodes_by_target[$i]->[1]]->[2] = $i;
		$tsc_pt = $nodes_by_target[$i]->[1];
	}
}
$tsc_pt = $nodes_by_target[scalar(@nodes_by_target)-1]->[1]; # the last one
$tsc_ids[$tsc_pt]->[3] = scalar(@nodes_by_target)-1;

#### produce @nodes_by_query
@nodes_by_query = sort { $a->[5] <=> $b->[5] or $a->[7] <=> $b->[7] } @nodes;
$qsc_pt = $nodes_by_query[0]->[5]; # the first one
$qsc_ids[$qsc_pt]->[2] = 0;
for(my $i=1;$i<scalar(@nodes_by_query);$i++){ #the ones between the first and the last
	if($qsc_pt != $nodes_by_query[$i]->[5]){
		$qsc_ids[$nodes_by_query[$i-1]->[5]]->[3] = $i-1;
		$qsc_ids[$nodes_by_query[$i]->[5]]->[2] = $i;
		$qsc_pt = $nodes_by_query[$i]->[5];
	}
}
$qsc_pt = $nodes_by_query[scalar(@nodes_by_query)-1]->[5]; # the last one
$qsc_ids[$qsc_pt]->[3] = scalar(@nodes_by_query)-1;

#### Cut the scaffold into portions to produce the scaffold portion data structure		
print "Cut the scaffold into portions ...\n";

#### processing target scaffolds
my $no_hit_sc_count = 0;
for(my $i=0;$i<scalar(@tsc_ids);$i++){
	if($tsc_ids[$i]->[2]<0){
		$tsc_portions[$i]->[0] = [0, $i, 0, $tsc_ids[$i]->[1], -1, 0, 0, 0];
		$no_hit_sc_count++; # running messages		
	}else{
		my $cnode_id=$tsc_ids[$i]->[2];
		my $cnode=$nodes_by_target[$cnode_id];
		my ($cnode_start,$cnode_end)=($cnode->[3],$cnode->[3]+$cnode->[4]);
		my $pnode_end=0;
		my $portion_id = 0;
		
		$tsc_portions[$i]->[$portion_id] = [$portion_id,$i,$pnode_end,$cnode_start-$pnode_end,-1,0,1,0];

		if(($cnode_start-$pnode_end)>0){ #### a 0-lengthed portion flanking a node portions
			$portion_id++;
			$tsc_portions[$i]->[$portion_id] = [$portion_id,$i,$cnode_start,0,-1,0,1,0];
		}
		
		$portion_id++;
		$tsc_portions[$i]->[$portion_id] = [$portion_id,$i,$cnode_start,$cnode_end-$cnode_start,$cnode_id,0,1,0];
		$cnode->[13]=$portion_id;
		
		for($cnode_id=$tsc_ids[$i]->[2]+1;$cnode_id<=$tsc_ids[$i]->[3];$cnode_id++){
			$pnode_end=$cnode_end;
			$cnode=$nodes_by_target[$cnode_id];
			($cnode_start,$cnode_end)=($cnode->[3],$cnode->[3]+$cnode->[4]);

			if(($cnode_start-$pnode_end)>0){ #### a 0-lengthed portion flanking a node portions
				$portion_id++;
				$tsc_portions[$i]->[$portion_id] = [$portion_id,$i,$pnode_end,0,-1,0,1,0];
			}
			
			$portion_id++;
			$tsc_portions[$i]->[$portion_id] = [$portion_id,$i,$pnode_end,$cnode_start-$pnode_end,-1,0,1,0];

			if(($cnode_start-$pnode_end)>0){ #### a 0-lengthed portion flanking a node portions
				$portion_id++;
				$tsc_portions[$i]->[$portion_id] = [$portion_id,$i,$cnode_start,0,-1,0,1,0];
			}			

			$portion_id++;			
			$tsc_portions[$i]->[$portion_id] = [$portion_id,$i,$cnode_start,$cnode_end-$cnode_start,$cnode_id,0,1,0];
			$cnode->[13]=$portion_id;
		}

		if(($tsc_ids[$i]->[1]-$cnode_end)>0){ #### a 0-lengthed portion flanking a node portions
			$portion_id++;
			$tsc_portions[$i]->[$portion_id] = [$portion_id,$i,$cnode_end,0,-1,0,1,0];
		}
				
		$portion_id++;
		$tsc_portions[$i]->[$portion_id] = [$portion_id,$i,$cnode_end, $tsc_ids[$i]->[1]-$cnode_end, -1,0,1,0];		
	}
}
print "$no_hit_sc_count target scaffolds have no alignment hit !\n";

#### processing query scaffolds
$no_hit_sc_count = 0;
for(my $i=0;$i<scalar(@qsc_ids);$i++){
	if($qsc_ids[$i]->[2]<0){
		$qsc_portions[$i]->[0] = [0, $i, 0, $qsc_ids[$i]->[1], -1, 0, 0,0];
		$no_hit_sc_count++; # running messages		
	}else{
		my $cnode_id=$qsc_ids[$i]->[2];
		my $cnode=$nodes_by_query[$cnode_id];
		my ($cnode_start,$cnode_end)=($cnode->[7],$cnode->[7]+$cnode->[8]);
		my $pnode_end=0;
		my $portion_id = 0;
		
		$qsc_portions[$i]->[$portion_id] = [$portion_id,$i,$pnode_end,$cnode_start-$pnode_end,-1,0,1,0];
		
		if(($cnode_start-$pnode_end)>0){ #### a 0-lengthed portion flanking a node portions
			$portion_id++;
			$qsc_portions[$i]->[$portion_id] = [$portion_id,$i,$cnode_start,0,-1,0,1,0];
		}
		
		$portion_id++;
		$qsc_portions[$i]->[$portion_id] = [$portion_id,$i,$cnode_start,$cnode_end-$cnode_start,$cnode_id,0,1,0];
		$cnode->[14]=$portion_id;
		
		for($cnode_id=$qsc_ids[$i]->[2]+1;$cnode_id<=$qsc_ids[$i]->[3];$cnode_id++){
			$pnode_end=$cnode_end;
			$cnode=$nodes_by_query[$cnode_id];
			($cnode_start,$cnode_end)=($cnode->[7],$cnode->[7]+$cnode->[8]);
			
			if(($cnode_start-$pnode_end)>0){ #### a 0-lengthed portion flanking a node portions
				$portion_id++;
				$qsc_portions[$i]->[$portion_id] = [$portion_id,$i,$pnode_end,0,-1,0,1,0];
			}
			
			$portion_id++;
			$qsc_portions[$i]->[$portion_id] = [$portion_id,$i,$pnode_end,$cnode_start-$pnode_end,-1,0,1,0];
			
			if(($cnode_start-$pnode_end)>0){ #### a 0-lengthed portion flanking a node portions
				$portion_id++;
				$qsc_portions[$i]->[$portion_id] = [$portion_id,$i,$cnode_start,0,-1,0,1,0];
			}			

			$portion_id++;			
			$qsc_portions[$i]->[$portion_id] = [$portion_id,$i,$cnode_start,$cnode_end-$cnode_start,$cnode_id,0,1,0];
			$cnode->[14]=$portion_id;
		}
		
		if(($qsc_ids[$i]->[1]-$cnode_end)>0){ #### a 0-lengthed portion flanking a node portions
			$portion_id++;
			$qsc_portions[$i]->[$portion_id] = [$portion_id,$i,$cnode_end,0,-1,0,1,0];
		}		
		
		$portion_id++;
		$qsc_portions[$i]->[$portion_id] = [$portion_id,$i,$cnode_end, $qsc_ids[$i]->[1]-$cnode_end, -1,0,1,0];		
	}
}
print "$no_hit_sc_count query  scaffolds have no alignment hit !\n";

#############################################
#### add in info about Ns and lowcases count
#### nodes[25 .. 28], [23]-pointer to mirror node, nodes_by_target, tsc_ids[2,3], sc_names{}

print "Appending infomation of Ns and lowcases to nodes ...\n"; 
my ($seq,$line,$name,$is_the_end,$ss,$Ns,$lcs);

for(my $i=0;$i<scalar(@nodes);$i++){
	@{$nodes[$i]}[25,26,27,28]=(-1,-1,-1,-1);
}

unless(exists($para{'--noNsLCs'}) and $para{'--noNsLCs'}==1){ 
	
open(my $gzFH,"gunzip -c $Species[0].fa.gz |") or die "Can't create pipe gunzip -c $Species[0].fa.gz !\n";

#### look for the first seq
while($line=<$gzFH>){
	next if $line =~ m/^\s+/;
	last if $line =~ m/^>/;
}	
die "Problems with $Species[0].fa.gz on line ...\n $line \n" if !defined($line);

#### store the seq
$is_the_end=0;
while($is_the_end==0){
	chomp $line;
	$line =~ /^>(.+)/;
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
	for(my $i=$tsc_ids[$sc_names{$name}]->[2];$i>=0 and $i<=$tsc_ids[$sc_names{$name}]->[3];$i++){
		($Ns,$lcs)=(0,0);
		$ss=substr($seq,$nodes[$i]->[3],$nodes[$i]->[4]);
		$Ns= $ss=~tr/Nn/xx/;
		$lcs= $ss=~tr/acgt/xxxx/;
		($nodes[$i]->[25],$nodes[$i]->[26])=($Ns,$lcs);
	}
}
close $gzFH;

for(my $i=0;$i<scalar(@nodes);$i++){
	($nodes[$i]->[23]->[27],$nodes[$i]->[23]->[28])=($nodes[$i]->[25],$nodes[$i]->[26]);
}

}

#############################################
#### output other infomations
$,="\t";

####
if(-f "$out_dir/hm.scaffolds"){
	unlink "$out_dir/hm.scaffolds.bak";
	rename "$out_dir/hm.scaffolds","$out_dir/hm.scaffolds.bak";
}
print "Output the scaffold infomation to $out_dir/hm.scaffolds.\n";
open(my $scaffoldFH,">$out_dir/hm.scaffolds") or die "Can not create $out_dir/hm.scaffolds! \n";
print $scaffoldFH "#### You may open this file in excel.\n\n";
print $scaffoldFH "#### This table is not intended to be manually edited.\n\n";
print $scaffoldFH "#tsc_id\tscaffold_name\tsize\tstart_of_node_id\tend_of_node_id\tis_incorporated\n";
for(my $i=0;$i<scalar(@tsc_ids);$i++){
	print $scaffoldFH $i,@{$tsc_ids[$i]};
	print $scaffoldFH "\n";
} 
close $scaffoldFH;

####
if(-f "$out_dir/hm.sc_portions"){
	unlink "$out_dir/hm.sc_portions.bak";
	rename "$out_dir/hm.sc_portions","$out_dir/hm.sc_portions.bak";
}
print "Output the scaffold portion infomation to $out_dir/hm.sc_portions.\n";
open(my $tsc_portionFH,">$out_dir/hm.sc_portions") or die "Can not create $out_dir/hm.sc_portions! \n";
print $tsc_portionFH "#### You may open this file in excel.\n";
print $tsc_portionFH "#### In this table, only the *manual_breaking column can be manually modified,\n";
print $tsc_portionFH "#### which can be set to 1 (0:no breaking, 1:to break the scaffold at this portion).\n";
print $tsc_portionFH "#### Note that column broken_by_assembly_error overrides *manual_breaking.\n\n";
print $tsc_portionFH "#### After editing or updating, this file may be saved as hm.sc_portions_edited or other names\n";
print $tsc_portionFH "#### in case being overwritten by another invocation of HM_pathFinder_preparation.pl.\n\n";
print $tsc_portionFH "#tsc_name\tportion_id\ttsc_id\tstart_position\tlen\tnode_id\tnode_status\tportion_status\tis_incorporated\tbroken_due_to_assembly_errors\t*manual_breaking\n";
for(my $i=0;$i<scalar(@tsc_portions);$i++){
	my $portion_count=scalar(@{$tsc_portions[$i]});
	my $portion_ref; 
	for(my $j=0;$j<$portion_count;$j++){
		$portion_ref = $tsc_portions[$i]->[$j];
		print $tsc_portionFH $tsc_ids[$portion_ref->[1]]->[0],@{$portion_ref},0,0;
		print $tsc_portionFH "\n";
	}
}
close $tsc_portionFH;

####
if(-f "$out_dir/hm.nodes"){
	unlink "$out_dir/hm.nodes.bak";
	rename "$out_dir/hm.nodes","$out_dir/hm.nodes.bak";
}
print "Output the node (alignment block) infomation to $out_dir/hm.nodes.\n";
open(my $nodeFH,">$out_dir/hm.nodes") or die "Can not create $out_dir/hm.nodes! \n";
print $nodeFH "#### You may open this file in excel.\n";
print $nodeFH "#### In this table, only the *manual_deletion column can be manually modified,\n";
print $nodeFH "#### which can be manually set to 1 (0:not deleted, 1:manually deleted).\n\n";
print $nodeFH "#### After editing or updating, this file may be saved as hm.nodes_edited or other names\n";
print $nodeFH "#### in case being overwritten by another invocation of  HM_pathFinder_preparation.pl.\n\n";
print $nodeFH "#tsc_name\tqsc_name\tnode_id\ttsc_id\tstrand\ttstart\ttlen\ttend\tqsc_id\tstrand\tqstart\tqlen\tqend\tscore\tali_len\tpredecessor\tsuccessor\ttsc_portion_ref\tqsc_portion_ref";
print $nodeFH "\ttsc_direction_predecessor\ttsc_direction_successor\tqsc_direction_predecessor\tqsc_direction_successor";
print $nodeFH "\tdeletion_visit_flag\tchain_id\tnet_id\tis_mirror\tnew_scaffold_id\torder_in_new_scaffolds\tNs_of_tsc\tlowcases_of_tsc\tNs_of_qsc\tlowcases_of_qsc\t*manual_deletion\n";
for(my $i=0;$i<scalar(@nodes);$i++){
	print $nodeFH $tsc_ids[$nodes[$i]->[1]]->[0],$qsc_ids[$nodes[$i]->[5]]->[0],@{$nodes[$i]}[0 .. 4],$nodes[$i]->[3]+$nodes[$i]->[4],@{$nodes[$i]}[5 .. 8],$nodes[$i]->[7]+$nodes[$i]->[8],@{$nodes[$i]}[9 .. 22],-1,-1,@{$nodes[$i]}[25 .. 28],0;
	print $nodeFH "\n";
}
close $nodeFH;

if($para{"--Delete"}==1){
	#unlink "$out_dir/filtered.net";
	#"$out_dir/filtered.chain.gz";
	#unlink "$out_dir/filtered.rbest.net","$out_dir/filtered.rbest.chain.gz";
}

print "Finished preparation for pathFinder.pl.\n";

print "\n\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";

