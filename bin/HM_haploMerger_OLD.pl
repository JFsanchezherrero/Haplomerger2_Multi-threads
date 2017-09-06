#!/usr/bin/perl  

## Version 2.10

## Introduction.
##

## notes.
## @unincorps contain those unpaired scaffolds or scaffold portions

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 1 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit; }

  This script is to create haplome (assembly) from the table produced by
HM_pathFinder.pl (hm.new_scaffolds). The table file hm.new_scaffolds 
contains all the descriptive information of the relation between two haplotpes
of the assembly. The script can also output unpaired sequences.
  
  --help or no arguments - show this message, undefined inputs are ignored

  --Species -   mandatory, to provide species names which are wanted to be
                processed
                NOTE THAT the species name should be given
                in ORDER, because the name comes first will be used
                as the target.
                Example: --Species bbev1a bbev1b
                In this script, both names refer to the same genome assembly,
                the name come first is used as the target genome, the second
                name is used as the query genome.
                
                The hm.new_scaffolds file should be place under the directory
                ./species0.species1.result.
                  
  --minOverlap=<number> - set the minimum required alignment length which used
                to connect two different scaffolds(default=0)
                NOTE that set to a very large figure (999,999,999) causes the
                output of an assembly without extension from one haplotypes
                to another
  
  --selectLongHaplotype - when two haplotype is avalable, select the longer
                one into the final assembly. This may reduce gene 
                losses caused by micro-assembly errors, but as a side effect, 
                it causes many switches between haplotypes (def=0, not perform)
                Note that this option will not override manual setting.
  --minLen=N - consider portions with length >=N (def=5000bp)
  --relativeLen=N - select the haplotype with N% longer (def=10 %)
  --absoluteLen=N - select the hapotype with N bp longer (def=2500bp)             
                
  --gapFilling - fill N gaps with another haplotype sequence. As a side-effect,
               this will cause many switches between haplotypes (def=0)
  --minAli=N   - min-len of alignments considered for gapFilling (def=4000bp)               
  --maxGapLeft=N - a successful gap filling require that after filling 
                 the Ngap size should be <=N% of the before (def=10%)                 

  --haploMergingMode=<original|updated> -
                ***original***, 
                to use the orginal data in the table hm.new_scaffolds.
                This is the default behavoir.
                ***updated***, 
                Table hm.new_scaffolds can be automatically updated or
                manually edited. 
                This option allows the use the updated or edited data
                in this table.
      
  --updatedNewScaffolds=<new_scaffold_table_file> - 
                required for --haploMergingMode=updated
  --updatedConnections=<connetions_file> - 
                optional for --haploMergingMode=updated                

  --optimalContinuity - output optiNewScaffolds.fa.gz(this is the default setting)
  --unpaired - output unpaired.fa.gz (this is the default setting)
  --twoHaplotypes - output two haplotype assemblies (def=1, yes)
	 
  --verbose - print more information       
  --Force - force to overwrite existed result files
  --Delete - delete intermediate and temp files (no effect on logfiles)
  
  NOTES: Softmasked sequences will be retained in the new assembly.
  
  Required files:
  ~result/hm.n_gaps (for --gapFilling),
  ~result/hm.new_scaffolds,
  ~result/hm.unpaired,
  ~result/hm.new_scaffolds(_edited), (for update mode),
  ~result/hm.connections(_edited), (for update mode),
  species.fa.gz.
   
  
  There can be five output files:
  1) optiNewScaffolds.fa.gz, a gzipped fasta file containing the haplome
     assembly that 1) selects one representative from two haplotypes, 
     2) extends scaffold length by connecting old scaffolds from different
     haplotypes based on overlapping alignments, 3) avoids switchs of different
     haplotypes inside a scaffold.
  2) (not availabe at this version) maxiNewScaffolds.fa.gz,
     a gzipped fasta file containing the haplome
     assembly that 1) selects one representative from two haplotypes, 
     2) extends scaffold length by connecting old scaffolds from different
     haplotypes based on overlapping alignments, 3) maximizes the contig N50
     length. To achieve 3), old scaffolds from different 
     haplotypes are first connected based on overlapping alignments
     in order to obtain a longer scaffold; then for each alignments between
     two haplotypes, the portion of the haplotype with better continuity is
     selected into the new assembly. As a side effect, this aseembly contains
     many switches between different haplotypes.  
  3) unpaired.fa.gz, contains raw unpaired sequences.
  4) assembly_hapA.fa.gz. Haplotype A assembly (of the target assembly).
  5) assembly_hapB.fa.gz. Haplotype B assembly (of the query  assembly).    
  6) hm.unpaired_updated, describe the unpaired sequences.
  7) hm.connections, describe the connections used for optiNewScaffolds.fa.gz.
  
  There will one more output file when invocate --haploMergingMode=updated:
  8) hm.new_scaffolds_updated 
     (containing updated data from hm.new_scaffolds(_edited) 
      and hm.connections(_edited)).    
	 
USAGES

print "========== Start at "; system("date"); print "\n\n";
print "$0 $Arg_list\n";
my $timer = time();

#### set the output_field_separator to " " and autoflushing
$,=' ';
$|=1;

# special parameters
# Store species names and result directory name
my (@Species, $out_dir); 
unless ($Arg_list =~ m/--Species\s+([^-]*)/) {die "No --Species argument or species_names found!\n"; }
unless (@Species = $1 =~ m/(\w+)/g)  {die "No species names found!\n"; }
$out_dir = "$Species[0].$Species[1].result";
print "Species included: ", @Species, "\n";

#### parse the arguments 
my %para = ("--minOverlap" => 0,  "--Force" => 0, "--Delete" => 0, "--verbose" => 0,
						"--selectLongHaplotype" => 0, "--minLen" => 5000, "--relativeLen" => 10, "--absoluteLen" => 2500, 
						"--gapFilling" => 0, "--minAli" => 4000, "--maxGapLeft"=>10,
						"--haploMergingMode" => "original", 
						"--updatedNewScaffolds" => "",
						"--updatedConnections" => "",
						"--optimalContinuity" => 1, "--unpaired" => 1, "--twoHaplotypes" => 1);

# parameters with values
while($Arg_list =~ m/(--\S+)=(\S+)/g){ 
	if(exists($para{$1})){
		$para{$1}=$2;
	}else{
		die "Parameter $1 is unknown!\n";
	}
}

# parameters of toggle type
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

#### set the over-write flag and deletion flag
print "Set to OVER-WRITING mode!\n" if $para{"--Force"}>0;
print "Set to DELETING mode!\n" if $para{"--Delete"}>0;
print "Overlapping alignments with length < ",$para{'--minOverlap'}, " are not used for connecting two old scaffolds.\n";
print "Option --selectLongHaplotype is set to ",$para{'--selectLongHaplotype'}," !\n";
print "Option --gapFilling is set to ",$para{'--gapFilling'}," !\n";
print "Option --haploMergingMode is set to ",$para{'--haploMergingMode'}," !\n";
print "Option --updatedNewScaffolds is set to ",$para{'--updatedNewScaffolds'}," !\n";
print "Option --updatedConnections is set to ",$para{'--updatedConnections'}," !\n";
print "Option --optimalContinuity is set to ",$para{'--optimalContinuity'}," !\n";
print "Option --unpaired is set to ",$para{'--unpaired'}," !\n";

#########################################
#### portion info
# $portions[id]->[ 0:new_scaffold_id, x:new_scaffold_len(deleted), 1:sub_id, 2:portion_id, 3:active_portion,
#                  4:old_scaffold1_name,  5:old_scffold1_size,  6:old_scaffold1_id,  7:start1, x:end1,  8:strand1,  9:len1,
#									10:old_scaffold2_name, 11:old_scffold2_size, 12:old_scaffold2_id, 13:start2, x:end2, 14:strand2, 15:len2,
#									16:connection,         17:score,
#						>>>		18:N_count_1,          19:N_count_2,
#									20:seq1,               21:seq2,
#                 22:optimal_nsc_internal_id,     23:active_portion_for_optimal_nsc,
#									24:minimal_nsc_internal_id,     25:active_portion_for_minimal_nsc,
#									26:maximal_nsc_internal_id,			27:active_portion_for_maximal_nsc,
#									28:optimal_nsc_position_start,	29:optimal_nsc_position_end 
#                 30:tsc_Ns%	31:tsc_LCs%	32:qsc_Ns%	33:qsc_LCs%	
#                 34:active_portion_updated	35:connection_updated	36:*active_portion_manual	37*connection_manual ]
my (@portions, @portions_1st, @portions_2nd);

#########################################
#### old scaffold position info
# $sc_1st/2nd{'old_scaffold_name'}->[ 0:start_portion_id_on_@portions_1st/2nd,  1:start_portion_id_on_@portions_1st/2nd ]; 
my (%osc_1st, %osc_2nd);

#########################################
#### new scaffold info
# $xxx_nsc[id]->[ 0:start_portion_id, 1:end_portion_id, 2:size, 3:refined_size ];
my (@opt_nsc, @min_nsc, @max_nsc, $zeros);

#########################################
#### new scaffold connections
# $connections[$i]->[ 0:new_scaffold_id, 1:sub_id, 2:portion_id,
#                     3:new_scaffold_name, 4:new_scaffold_position_start(0-based), 5:new_scaffold_position_end,
# 										6:target_old_scaffold_name, 7:target_old_scaffold_postion_start, 8:target_old_scaffold_postion_end,
# 										9:query_old_scaffold_name,  10:query_old_scaffold_postion_start, 11:query_old_scaffold_postion_end,
# 										12:connection,              13:score,  14:tsc_NsLCs%, 15:qsc_NsLCs%, 16:connection_updated,  17:*connection_manual              ]
#my (@connections);

#########################################
#### ngaps info
# @ngaps, [id]->[0:reserved, 1:scaffold_name, 2:scaffold_id, 3:scaffold_size, 4:Ngap_start, 5:Ngap_end, 6:Ngap_len,
#                7:node_id, 8:chain_id, 9:ali_len, 10:is_mirror, 11:strand_sign, 
#                12:tsc_id, 13:tsc_start, 14:tsc_end, 15:qsc_id, 16:qsc_start, 17:qsc_end,
#                18:tsc_replace_start, 19:tsc_replace_end, 20:qsc_replace_start, 21:qsc_replace_end, 22:new_Ngap_vs_old_Ngap_% ]
#
my (@ngaps, %sc_ids);

my ($haplome_size, $haplome_scaffold_N50_size,$haplome_scaffold_N50_number)=(0,0,0);
my ($raw_size_of_unpaired_scaffolds, $raw_size_of_unpaired_scaffold_portions)=(0,0);
my ($temp,$temp1,$temp2);

sub read_and_preprocess_portion_data();
sub read_Ngaps_info();
sub read_and_process_sequences(); # also output unpaired sequences
sub output_optimal_assembly();

read_and_preprocess_portion_data();
read_Ngaps_info() if $para{'--gapFilling'}==1;
read_and_process_sequences(); # also output haplotypeA/B and unpaired sequences
output_optimal_assembly() if $para{'--optimalContinuity'}==1;

$,="";

### output size statistic
print "\n";
print "The haplome size is $haplome_size bp (only caculate sequences with length >= 500 bp);\n";
print "the haplome scaffold N50 number is $haplome_scaffold_N50_number;\n";
print "the haplome scaffold N50 size is $haplome_scaffold_N50_size bp;\n";
print "the raw size of unpaired scaffolds is $raw_size_of_unpaired_scaffolds bp (>= 500 bp);\n";
print "the raw size of unpaired scaffold portions is $raw_size_of_unpaired_scaffold_portions bp (>= 500 bp).\n";

###
print "\nPrinting the optimalContinuity assembly data ...\n";
print "#Sc_scaffold_id\tstart_portion_id\tend_portion_id\tsize\trefined_size\n";
for(my $i=0;$i<scalar(@opt_nsc);$i++){
	$zeros='0' x (7-length($i));
	print "Sc",$zeros,$i,"\t",$opt_nsc[$i]->[0],"\t",$opt_nsc[$i]->[1],"\t",$opt_nsc[$i]->[2],"\t",$opt_nsc[$i]->[3],"\n";
}

###
if($para{'--verbose'}>0){
	print "\nVERBOSE: Printing new_scaffold info table ...\n";
	print "0:new_scaffold_id\t1:sub_id\t2:portion_id\t3:active_portion\t";
	print "4:old_scaffold1_name\t5:old_scffold1_size\t6:old_scaffold1_id\t7:start1\t8:strand1\t9:len1\t";
	print "10:old_scaffold2_name\t11:old_scffold2_size\t12:old_scaffold2_id\t13:start2\t14:strand2\t15:len2\t";
	print "16:connection\t17:score\t18:N_count_1\t19:N_count_2\t";
	print "22:optimal_nsc_internal_id\t23:active_portion_for_optimal_nsc\t24:minimal_nsc_internal_id\t25:active_portion_for_minimal_nsc\t26:maximal_nsc_internal_id\t27:active_portion_for_maximal_nsc\t28:optimal_nsc_position_start\t29:optimal_nsc_position_end\ttsc_NsLCs%\tqscNsLCs%\n";
	$,="\t";
	for(my $i=0;$i<scalar(@portions);$i++){
		my $pp=$portions[$i];
		print @$pp[0 .. 19],@$pp[22 .. 31],"\n"; 
	}
	print "\n";
	$,=' '; 
}

###
if($para{'--verbose'}>0){
	print  "\n#filled\tscaffold_name\tsc_id\tsc_size\tgap_start\tgap_end\tgap_len";
	print  "\tnode_id\tchain_id\tali_len\tis_mirror_or_not\tstrand_sign";
	print  "\ttsc_id\ttsc_start\ttsc_end\tqsc_id\tqsc_start\tqsc_end";
	print  "\ttsc_start_of_being_replaced_portion\ttsc_end_of_being_replaced_portion";
	print  "\tqsc_start_of_replacing_portion\tqsc_end_of_replacing_portion\ttsc_Ngap_vs_qsc_Ngap_%\n";
	$,="\t";
	for(my $i=0;$i<scalar(@ngaps);$i++){
		print  @{$ngaps[$i]}[0 .. 22];
		print  "\n";
	}
	$,="";
}

#########################################
if($para{"--Delete"}==1){
}

print "\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";

######################################### subroutines #########################################

#########################################
#
#########################################
sub read_and_preprocess_portion_data(){

	my ($temp,$temp1);
	my $updatedScFH;
	
	#### reading hm.connections(_edited)
	my %conn;
	if($para{"--haploMergingMode"} eq "updated"){
		print "Runing on UPDATE mode ... !\n";
		
		##
		$temp=$para{"--updatedConnections"};
		unless(-f $temp){
			print "Connection file $temp for UPDATE mode is not existed, continune ... \n";
		}else{
			print "Connection file $temp for UPDATE mode is existed, continune ... \n";
			open(my $conFH, "<$temp") or die "Can not read connection file $temp !\n";
			my $tt;
			while(<$conFH>){
				next if m/^#|^\s/;
				chomp;
				$tt = [ split /\t/ ];
				$tt->[18] = $tt->[19]!=0 ? $tt->[19] : $tt->[18];
				if($tt->[18]!=0){	$conn{$tt->[0].'_'.$tt->[1].'_'.$tt->[2]} = $tt->[18]>0 ? 1 : -1; }
			}
			close $conFH;
		}
		
		##
		$temp=$para{"--updatedNewScaffolds"};
		unless(-f $temp){ 
			print "New scaffold file $temp for UPDATE mode is not found, die.\n"; die; 
		}else{
			print "New scaffold file $temp for UPDATE mode is found, contiune ...\n"; 
		}
		
		##
		if(-f "$out_dir/hm.new_scaffolds_updated"){
			unlink "$out_dir/hm.new_scaffolds_updated.bak";
			rename "$out_dir/hm.new_scaffolds_updated","$out_dir/hm.new_scaffolds_updated.bak";		
		}
		open($updatedScFH,">$out_dir/hm.new_scaffolds_updated") or die "Can not create $out_dir/hm.new_scaffolds_updated!\n";
	}	
	
	#### open hm.new_scaffolds files
	$temp= "$out_dir/hm.new_scaffolds" if $para{"--haploMergingMode"} ne "updated";
	print "reading $temp file and creating @portions data structure ... \n";
	open(my $portionsFH, "<$temp") or die "Can not read new_scaffolds file $temp !\n";
	
	#### creating @portions data structure
	while(<$portionsFH>){
		if(m/^#|^\s/){
			print $updatedScFH $_ if $para{"--haploMergingMode"} eq "updated"; 
			next; 
		}
		chomp;
		my $tt = [ split /\t/ ];
		if($para{"--haploMergingMode"} eq "updated"){
			
			# integrating the updated data from hm.connections* and ouput them to hm.new_scaffolds_updated
			$tt->[28]=$conn{$tt->[0].'_'.$tt->[2].'_'.$tt->[3]} if exists($conn{$tt->[0].'_'.$tt->[2].'_'.$tt->[3]}) and $tt->[28]==0;
			$,="\t";
			print $updatedScFH @$tt; print $updatedScFH "\n";
			$,=' ';
			
			# using updated data from hm.new_scaffolds*
			$tt->[4]  = ($tt->[25]==1 or $tt->[25]==2) ? $tt->[25] : $tt->[4]; 
			$tt->[19] = $tt->[19]!=1 ? $tt->[19] : ($tt->[26]==1 or $tt->[26]==-1) ? $tt->[26]*9 : $tt->[19];
			
			# using manually-edited data from hm.new_scaffolds*
			$tt->[4]  = ($tt->[27]==1 or $tt->[27]==2) ? $tt->[27] : $tt->[4];
			$tt->[19] = ($tt->[19]==0 or $tt->[19]==-1 or $tt->[19]==-2) ? $tt->[19] : ($tt->[28]==1 or $tt->[28]==-1) ? $tt->[28]*9 : $tt->[19];

		}
		push @$tt,(-1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1,$tt->[21],$tt->[22],$tt->[23],$tt->[24],$tt->[25],$tt->[26],$tt->[27],$tt->[28]); # move NsLCs/updating_data to the tail
		splice @$tt,21,8; # !!! delete NsLCs% 
		splice @$tt,16,1;
		splice @$tt, 9,1;
		splice @$tt, 1,1;
		push @portions, $tt;
	}
	close $portionsFH;
	close $updatedScFH if $para{"--haploMergingMode"} eq "updated";
	
	@portions     = sort { $a->[0] <=> $b->[0] or $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] } @portions;
	@portions_1st = sort { $a->[4] cmp $b->[4]   } @portions;
	@portions_2nd = sort { $a->[10] cmp $b->[10] } @portions;
	
	#### filling %osc_1st, %osc_2nd
	my ($csc1,$csc2)=('','');
	for(my $i=0;$i<scalar(@portions);$i++){
		if($csc1 ne $portions_1st[$i]->[4]){ 
			$osc_1st{$csc1}->[1] = $i-1 if exists $osc_1st{$csc1};
			$csc1 = $portions_1st[$i]->[4];
			$osc_1st{$csc1}->[0] = $i;
		}
		if($csc2 ne $portions_2nd[$i]->[10]){ 
			$osc_2nd{$csc2}->[1] = $i-1 if exists $osc_2nd{$csc2};
			$csc2 = $portions_2nd[$i]->[10];
			$osc_2nd{$csc2}->[0] = $i;
		}
	}
	$osc_1st{$csc1}->[1] = scalar(@portions_1st)-1;
	$osc_2nd{$csc2}->[1] = scalar(@portions_2nd)-1;
	
	#### processing connecting junctions between old haplotype scaffolds
	print "Processing junctions between old haplotype scaffolds ... \n";
	for(my $i=0;$i<scalar(@portions);$i++){
		if($portions[$i]->[16]==1){
			$portions[$i]->[16] = ($portions[$i]->[9]<$para{'--minOverlap'} or $portions[$i]->[15]<$para{'--minOverlap'}) ? -9 : 9;
		}
	}
	
}

#########################################
#
#########################################
sub read_Ngaps_info(){
	
	open(my $gapFH,"<$out_dir/hm.n_gaps") or die "Can not open $out_dir/hm.n_gaps!\n";
	while(<$gapFH>){	
		next if m/^#|^\s/;
		chomp;
		my $tt = [ split /\t/ ];
		# filtering the ngap records
		unshift @$tt,-1;
		if($tt->[7]>-1 and $tt->[9]>=$para{'--minAli'} and $tt->[18]>-1 and $tt->[22]<=$para{'--maxGapLeft'}){
			push @ngaps,$tt;
		}
	}
	close $gapFH;
	
	## order @ngaps by nodes, @node_ids
	@ngaps = sort {$a->[7]<=>$b->[7] or $a->[18] <=> $b->[18]} @ngaps;
	my %node_ids=();
	$temp = $ngaps[0]->[7]; # the first one
	$node_ids{$temp}->[0]=0;
	for(my $i=0;$i<scalar(@ngaps);$i++){ 
		if($temp != $ngaps[$i]->[7]){
			$node_ids{$temp}->[1] = $i-1;
			$temp = $ngaps[$i]->[7];
			$node_ids{$temp}->[0] = $i;
		}
	}
	$temp = $ngaps[scalar(@ngaps)-1]->[7]; # the last one
	$node_ids{$temp}->[1] = scalar(@ngaps)-1; 	
	
	## merge overlapping portions
	my @ngaps_tmp=@ngaps;
	@ngaps=();
	my ($sp,$ep);
	foreach my $nid (keys %node_ids){
		($sp,$ep)=($node_ids{$nid}->[0],$node_ids{$nid}->[0]);
		for(my $i=$node_ids{$nid}->[0]+1;$i<=$node_ids{$nid}->[1];$i++){
			if($ngaps_tmp[$i]->[18]<=$ngaps_tmp[$i-1]->[19]){
				$ep=$i;
			}else{
				if($ngaps_tmp[$sp]->[11]>0){
					( $ngaps_tmp[$sp]->[19],$ngaps_tmp[$sp]->[21] ) = ( $ngaps_tmp[$ep]->[19],$ngaps_tmp[$ep]->[21] );		
				}else{
					( $ngaps_tmp[$sp]->[19],$ngaps_tmp[$sp]->[20] ) = ( $ngaps_tmp[$ep]->[19],$ngaps_tmp[$ep]->[20] );
				}
				push @ngaps,$ngaps_tmp[$sp];
				($sp,$ep)=($i,$i);			
			}
		}
		if($ngaps_tmp[$sp]->[11]>0){
			( $ngaps_tmp[$sp]->[19],$ngaps_tmp[$sp]->[21] ) = ( $ngaps_tmp[$ep]->[19],$ngaps_tmp[$ep]->[21] );		
		}else{
			( $ngaps_tmp[$sp]->[19],$ngaps_tmp[$sp]->[20] ) = ( $ngaps_tmp[$ep]->[19],$ngaps_tmp[$ep]->[20] );
		}
		push @ngaps,$ngaps_tmp[$sp];		
	}
	@ngaps_tmp=();
	
	## order @ngaps by node id
	%sc_ids=();
	@ngaps = sort {$a->[2]<=>$b->[2] or $a->[18] <=> $b->[18]} @ngaps;
	$temp = $ngaps[0]->[2]; # the first one
	$sc_ids{$temp}->[0]=0;
	for(my $i=0;$i<scalar(@ngaps);$i++){ 
		if($temp != $ngaps[$i]->[2]){
			$sc_ids{$temp}->[1] = $i-1;
			$temp = $ngaps[$i]->[2];
			$sc_ids{$temp}->[0] = $i;
		}
	}
	$temp = $ngaps[scalar(@ngaps)-1]->[2]; # the last one
	$sc_ids{$temp}->[1] = scalar(@ngaps)-1; 	

}

#########################################
#
#########################################
sub read_and_process_sequences(){
	
	my $temp;
	my $line;
	my $is_the_end=0;
	my ($name,$seq);
	my $zero_fill;
	my $sclen;
	my ($pp,$ss);
	
	open(my $gzFH,"gunzip -c $Species[0].fa.gz |") or die "Can't create pipe gunzip -c $Species[0].fa.gz !\n";
	
	#### unincorps[id]->[ 0:scaffold_name, 1:size, 2:start, 3:len, 4:full_scaffold_or_part_of_the_scaffold, 5:N_counts, 6:lowcase_count, 7:seq ]
	my ($un_inFH,$un_outFH,@unincorps, %unindex);
	open($un_inFH,"<$out_dir/hm.unpaired") or die "Can't open $out_dir/hm.unpaired !\n";
	while(<$un_inFH>){
		next if m/^#|^\s/;
		chomp;
		my $tt = [ split /\t/ ];
		push @$tt,(0,0,undef);
		push @unincorps,$tt;
	}
	close $un_inFH;
	if(scalar(@unincorps)>0){ # when no unpaired sequences
		@unincorps = sort { $a->[0] cmp $b->[0] } @unincorps;
		$temp = $unincorps[0]->[0];
		$unindex{$temp}->[0]=0;
		for(my $i=0;$i<scalar(@unincorps);$i++){
			next if $temp eq $unincorps[$i]->[0];
			$unindex{$temp}->[1]=$i-1;
			$temp=$unincorps[$i]->[0];
			$unindex{$temp}->[0]=$i;
		}
		$unindex{$temp}->[1]=scalar(@unincorps)-1;
	}
	
	#### look for the first seq
	while($line=<$gzFH>){
		next if $line =~ m/^\s+/;
		last if $line =~ m/^>/;
	}	
	die "Problems with $Species[0].fa.gz on line ...\n $line \n" if !defined($line);
	
	#### store the seq
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
		if(exists $osc_1st{$name}){
			for(my $i=$osc_1st{$name}->[0];$i<=$osc_1st{$name}->[1];$i++){
				$pp=$portions_1st[$i];
				$ss=substr($seq,$pp->[7],$pp->[9]);
				####if($pp->[8]<0){ $ss=reverse($ss); $ss=~tr/ACGTacgt/TGCAtgca/; }
				$pp->[20]=$ss;
			}
		}
		
		####
		if(exists $osc_2nd{$name}){
			for(my $i=$osc_2nd{$name}->[0];$i<=$osc_2nd{$name}->[1];$i++){
				$pp=$portions_2nd[$i];
				$ss=substr($seq,$pp->[13],$pp->[15]);
				####if($pp->[14]<0){ $ss=reverse($ss); $ss=~tr/ACGTacgt/TGCAtgca/; }
				$pp->[21]=$ss;
			}
		}
		
		#### 
		#  unincorps[id]->[ 0:scaffold_name, 1:size, 2:start, 3:len, 4:full_scaffold_or_part_of_the_scaffold, 5:N_counts, 6:lowcase_count, 7:seq ]
		# my ($un_inFH,$un_outFH,@unincorps, %unindex);
		if(exists $unindex{$name}){
			for(my $i=$unindex{$name}->[0];$i<=$unindex{$name}->[1];$i++){
				$ss=substr($seq,$unincorps[$i]->[2],$unincorps[$i]->[3]);
				while($ss=~/(.)/g){
					$unincorps[$i]->[5]++ if $1 eq 'N' or $1 eq 'n';
					$unincorps[$i]->[6]++ if $1 =~ m/[acgt]/;
				}
				$unincorps[$i]->[7]=$ss;
			}
		}
	
	}
	close $gzFH;
	
	#### ##################################
	#### output haplotype A and B
	if($para{'--twoHaplotypes'}>0){
		my ($name,$size,$id,$sclen,$found_flag,$start,$end,$strand,$hapFH);
		
		#### haplotype A
		open ($hapFH,"|gzip -c >$out_dir/assembly_hapA.fa.gz") or die "Can not create |gzip -c >$out_dir/assembly_hapA.fa.gz\n";
		print "\nOutputing haplotype assembly A ... \n";
		print "#name\toriginal_size\thap_id\tstart\tend\tstrand\n";
		($name,$size,$id,$found_flag,$start,$end,$strand,$temp,$temp1)=(undef,0,-1,0,0,0,0,'',0);
		for(my $i=0;$i<scalar(@portions);$i++){
			next if ! defined($name) and $portions[$i]->[8]==0; 
			if( $portions[$i]->[8]==0 ){
				next if $found_flag<1;
				$id++;
				$sclen=length($temp);
				print $name."\t".$size."\tA".$id."\t".$start."\t".$end."\t".$strand."\n";
				print $hapFH ">".$name."_A".$id."\toriginal_size:".$size.";original_range:".$start."-".$end.";strand:".$strand."\n";
				while($sclen>0){ print $hapFH substr($temp,0,100,''),"\n"; $sclen-=100; }
				$found_flag=0;
				next;
			}elsif( $found_flag==1 and $portions[$i]->[4] ne $name ){
				$id++;
				$sclen=length($temp);
				print $name."\t".$size."\tA".$id."\t".$start."\t".$end."\t".$strand."\n";
				print $hapFH ">".$name."_A".$id."\toriginal_size:".$size.";original_range:".$start."-".$end.";strand:".$strand."\n";
				while($sclen>0){ print $hapFH substr($temp,0,100,''),"\n"; $sclen-=100; }
				$found_flag=0;				
			}
			if($found_flag==0){
				$strand=$portions[$i]->[8];
				$size=$portions[$i]->[5];
				$name=$portions[$i]->[4];			
				($start,$end) =	($portions[$i]->[7],$portions[$i]->[7]+$portions[$i]->[9]);
				$found_flag=1;				
			}
			($start,$end) = $strand>0 ? ($start,$portions[$i]->[7]+$portions[$i]->[9]) : ($portions[$i]->[7],$end); 
			$temp1 = $portions[$i]->[20];
			if($strand<0){ $temp1 = reverse $temp1; $temp1 =~ tr/ACGTacgt/TGCAtgca/; }
			$temp .= $temp1;
		}
		if($found_flag>0){
			$id++;
			$sclen=length($temp);
			print $name."\t".$size."\tA".$id."\t".$start."\t".$end."\t".$strand."\n";
			print $hapFH ">".$name."_A".$id."\toriginal_size:".$size.";original_range:".$start."-".$end.";strand:".$strand."\n";
			while($sclen>0){ print $hapFH substr($temp,0,100,''),"\n"; $sclen-=100; }
		}
		close $hapFH;
		
		#### haplotype B
		open ($hapFH,"|gzip -c >$out_dir/assembly_hapB.fa.gz") or die "Can not create |gzip -c >$out_dir/assembly_hapB.fa.gz\n";
		print "\nOutputing haplotype assembly B ... \n";
		print "#name\toriginal_size\thap_id\tstart\tend\tstrand\n";
		($name,$size,$id,$found_flag,$start,$end,$strand,$temp,$temp1)=(undef,0,-1,0,0,0,0,'',0);
		for(my $i=0;$i<scalar(@portions);$i++){
			next if ! defined($name) and $portions[$i]->[14]==0; 			
			if($portions[$i]->[14] == 0 ){
				next if $found_flag<1;
				$id++;
				$sclen=length($temp);
				print $name."\t".$size."\tB".$id."\t".$start."\t".$end."\t".$strand."\n";
				print $hapFH ">".$name."_B".$id."\toriginal_size:".$size.";original_range:".$start."-".$end.";strand:".$strand."\n";
				while($sclen>0){ print $hapFH substr($temp,0,100,''),"\n"; $sclen-=100; }
				$found_flag=0;
				next;				
			}elsif( $found_flag==1 and $portions[$i]->[10] ne $name ){
				$id++;
				$sclen=length($temp);
				print $name."\t".$size."\tB".$id."\t".$start."\t".$end."\t".$strand."\n";
				print $hapFH ">".$name."_B".$id."\toriginal_size:".$size.";original_range:".$start."-".$end.";strand:".$strand."\n";
				while($sclen>0){ print $hapFH substr($temp,0,100,''),"\n"; $sclen-=100; }
				$found_flag=0;				
			}
			if($found_flag==0){
				$strand=$portions[$i]->[14];
				$size=$portions[$i]->[11];
				$name=$portions[$i]->[10];
				($start,$end) =  ($portions[$i]->[13],$portions[$i]->[13]+$portions[$i]->[15]);
				$found_flag=1;
			}
			($start,$end) = $strand>0 ? ($start,$portions[$i]->[13]+$portions[$i]->[15]) : ($portions[$i]->[13],$end); 
			$temp1 = $portions[$i]->[21];
			if($strand<0){ $temp1 = reverse $temp1; $temp1 =~ tr/ACGTacgt/TGCAtgca/; }
			$temp .= $temp1;
		}
		if($found_flag>0){
			$id++;
			$sclen=length($temp);
			print $name."\t".$size."\tB".$id."\t".$start."\t".$end."\t".$strand."\n";
			print $hapFH ">".$name."_B".$id."\toriginal_size:".$size.";original_range:".$start."-".$end.";strand:".$strand."\n";
			while($sclen>0){ print $hapFH substr($temp,0,100,''),"\n"; $sclen-=100; }
		}
		close $hapFH;	
		
		print "\n\n";
	}	
	
	#### ##################################
	#### select active haplotype, gap filling, and switch to the right strand
	print "Select the the active haplotype for each portion ...\n";
	for(my $i=0;$i<scalar(@portions);$i++){
		if($para{'--selectLongHaplotype'}>0 
				and $portions[$i]->[34]+$portions[$i]->[36]==0
				and $portions[$i]->[9]>=$para{'--minLen'} and $portions[$i]->[15]>=$para{'--minLen'}){
			$temp =$portions[$i]->[9] *(1-$portions[$i]->[30]/100);
			$temp1=$portions[$i]->[15]*(1-$portions[$i]->[32]/100);
			if($temp>$temp1+$para{'--absoluteLen'} or $temp>$temp1+$temp1*$para{'--relativeLen'}/100){
				$portions[$i]->[3]=1;
			}elsif($temp1>$temp+$para{'--absoluteLen'} or $temp1>$temp+$temp*$para{'--relativeLen'}/100){
				$portions[$i]->[3]=2;
			}
		}
		($portions[$i]->[23],$portions[$i]->[25])=($portions[$i]->[3],$portions[$i]->[3]); 
	}
	##
	if($para{'--gapFilling'}==1){
		print "Fill N-gaps with the haplotype sequences ...\n";
		my ($fill_count,$fill_len)=(0,0);
		my ($tid,$tstart,$tend,$tseq,  $qid,$qstart,$qend,$qseq);
		for(my $i=0;$i<scalar(@portions);$i++){
			next if $portions[$i]->[16]<0;
			if($portions[$i]->[23]==1){
				($tid,$tstart,$tend,$tseq,  $qid,$qstart,$qend,$qseq) = (@{$portions[$i]}[6,7,9,20, 12,13,15,21]);
				$tend=$tstart+$tend; $qend=$qstart+$qend;
			}else{
				($tid,$tstart,$tend,$tseq,  $qid,$qstart,$qend,$qseq) = (@{$portions[$i]}[12,13,15,21, 6,7,9,20]);
				$tend=$tstart+$tend; $qend=$qstart+$qend;			
			}
			$temp=''; $temp1=$tstart;
			if(exists $sc_ids{$tid}){
				for(my $i=$sc_ids{$tid}->[0];$i<=$sc_ids{$tid}->[1];$i++){
					last if $ngaps[$i]->[18]>=$tend;
					next if $ngaps[$i]->[19]<=$tstart;
					next if $qid ne $ngaps[$i]->[15];					
					if($ngaps[$i]->[18]>=$tstart and $ngaps[$i]->[19]<=$tend){
						$fill_count++; $fill_len+=$ngaps[$i]->[19]-$ngaps[$i]->[18]; ####
						if($ngaps[$i]->[20]>=$qstart and $ngaps[$i]->[21]<=$qend){
							$ngaps[$i]->[0]=1;
							if($ngaps[$i]->[11]>0){
								$temp .= substr($tseq,$temp1-$tstart,$ngaps[$i]->[18]-$temp1);
								$temp1 = $ngaps[$i]->[19];
								$temp .= substr($qseq,$ngaps[$i]->[20]-$qstart,$ngaps[$i]->[21]-$ngaps[$i]->[20]);
							}else{
								$temp .= substr($tseq,$temp1-$tstart,$ngaps[$i]->[18]-$temp1);
								$temp1 = $ngaps[$i]->[19];
								$temp2 = reverse( substr(	$qseq,$ngaps[$i]->[20]-$qstart,$ngaps[$i]->[21]-$ngaps[$i]->[20]	)	);
								$temp2 =~tr/ACGTacgt/TGCAtgca/;
								$temp .= $temp2;						
							}
						}else{
							#print "xxdebuggingxx tid:$tid-$ngaps[$i]->[12]-$tstart-$tend; qid:$qid-$ngaps[$i]->[15]-$qstart-$qend; $ngaps[$i]->[18]>>$ngaps[$i]->[19]; $ngaps[$i]->[20]>>$ngaps[$i]->[21]\n";
							print "Coordinate problem no.1 !\n"; die;
						}
					}elsif($ngaps[$i]->[18]>=$tstart and $ngaps[$i]->[19]>$tend){
						print "Coordinate problems no.2 !\n"; die;
					}
				}
			}
			$temp .= substr($tseq,$temp1-$tstart,$tend-$temp1);
			if($portions[$i]->[23]==1){
				$portions[$i]->[20]=$temp;
			}else{
				$portions[$i]->[21]=$temp;			
			}			
		}
		print "Filled $fill_count N-gaps, totally length is $fill_len bp.\n";
	}
	##
	print "Adjust the strand sign of the selected haplotype sequence ...\n";
	for(my $i=0;$i<scalar(@portions);$i++){	
		if($portions[$i]->[23]==1){
			if($portions[$i]->[8]<0){
				$portions[$i]->[20]=reverse($portions[$i]->[20]);
				$portions[$i]->[20]=~tr/ACGTacgt/TGCAtgca/;
			}
		}else{
			if($portions[$i]->[14]<0){
				$portions[$i]->[21]=reverse($portions[$i]->[21]);
				$portions[$i]->[21]=~tr/ACGTacgt/TGCAtgca/;
			}			
		}
	}

	#### ##################################
	#### output unpaired information
	@unincorps = sort { $a->[4] cmp $b->[4] or $b->[3] <=> $a->[3] } @unincorps;
	## output hm.info.unpaired_updated
	$,="\t";
	if(-f "$out_dir/hm.unpaired_updated"){
		unlink "$out_dir/hm.unpaired_updated.bak";
		rename "$out_dir/hm.unpaired_updated","$out_dir/hm.unpaired_updated.bak";		
	}	
	open(my $unFH,">$out_dir/hm.unpaired_updated") or die "Can not open out_dir/hm.unpaired_updated !\n";
	print $unFH "# This table contains UPDATED information of the scaffolds and scaffold portions that failed to be incorporated\n";
	print $unFH "# into the new scaffold assebmly in the path_finding process.\n";
	print $unFH "# These sequence may be used for further evaluation to see if they could be incorporated.\n";
	#print $unFH "# Note that only scaffold portions longer than 5% of full scaffold length were output, \n";
	print $unFH "# and only sequences of at least 500bp were output. \n\n";
	print $unFH "# This table is not intended to be manually modified.\n\n";
	print $unFH "#new_scaffold_name\told_scaffold_name\tsize\tstart\tlen\tfull_scaffold_or_part_of_the_scaffold\tN_counts\tlowcase_count\n";
	for(my $i=0;$i<scalar(@unincorps);$i++){
		my $un=$unincorps[$i];
		$zero_fill='0' x (7-length($i));
		print $unFH "x".substr($unincorps[$i]->[4],0,1)."Sc".$zero_fill.$i,@{$un}[0 .. 6],"\n";
	}
	$,=' ';
	close $unFH;	 
	## Output unpaired_raw.fa.gz
	if($para{'--unpaired'}==1){
		print "Output unpaired.fa.gz (only output scaffolds with size >= 500 bp) ...\n";
		open($un_outFH,"| gzip -c >$out_dir/unpaired.fa.gz") or die "Can't create pipe | gzip -c >$out_dir/unpaired.fa.gz !\n";		
		$,='';
		for(my $i=0;$i<scalar(@unincorps);$i++){
			next if  $unincorps[$i]->[3]<500;
			$zero_fill='0' x (7-length($i));
			$sclen=length($unincorps[$i]->[7]);
			print $un_outFH ">x",substr($unincorps[$i]->[4],0,1),"Sc",$zero_fill,$i," old_",$unincorps[$i]->[0],";size:",$unincorps[$i]->[1],";start:",$unincorps[$i]->[2],";length:",$sclen,"\n";
			while($sclen>0){ print $un_outFH substr($unincorps[$i]->[7],0,100,''),"\n"; $sclen-=100; }
		}
		close $un_outFH;
		$,=' ';
	}
	
	## output total unpaired sequence sizes
	foreach $temp (@unincorps){
		next if $temp->[3]<500;	 #### should be < 2000
		$raw_size_of_unpaired_scaffolds += $temp->[3] if $temp->[4] eq "full";
		$raw_size_of_unpaired_scaffold_portions += $temp->[3] if $temp->[4] eq "part";
	}
	
	##
	@unincorps=();	
		
}

#########################################
#
#########################################
sub output_optimal_assembly(){
	
	print "Processing optimal_assembly data ...\n";
	my ($temp);
	my @sc_sizes;
	my $zero_fill;
	my $sclen;
	my $opt_id=0;
	my $size=0;
	my $psize=0;
	my $nsc_id=$portions[0]->[0];
	my $sub_id=$portions[0]->[1];
	$opt_nsc[$opt_id]->[0]=0;
	$opt_nsc[$opt_id]->[3]=0; # refined-size
	for(my $i=0;$i<scalar(@portions);$i++){
		if(    $nsc_id == $portions[$i]->[0] and $sub_id == $portions[$i]->[1] and $portions[$i]->[16] > -9){
			$psize=$size; #################
			$size += $portions[$i]->[23] == 1 ? $portions[$i]->[9] : $portions[$i]->[15];
			($portions[$i]->[28],$portions[$i]->[29])=($psize,$size) if abs($portions[$i]->[16])==9; #################
		}elsif($nsc_id == $portions[$i]->[0] and $sub_id == $portions[$i]->[1] and $portions[$i]->[16]== -9){
			if($size==0){
				$psize=$size; #################
				$size += $portions[$i]->[23] == 1 ? $portions[$i]->[9] : $portions[$i]->[15];
				($portions[$i]->[28],$portions[$i]->[29])=($psize,$size) if abs($portions[$i]->[16])==9; #################
			}else{
				$psize=$size; #################
				$size += $portions[$i]->[23] == 1 ? $portions[$i]->[9] : $portions[$i]->[15];
				($portions[$i]->[28],$portions[$i]->[29])=($psize,$size) if abs($portions[$i]->[16])==9; #################
				if($nsc_id == $portions[$i+1]->[0]){
					$opt_nsc[$opt_id]->[1]=$i;
					$opt_nsc[$opt_id]->[2] = $size;
					$size = 0;
					$opt_id++;
					$opt_nsc[$opt_id]->[0]=$i+1;
					$sub_id = $portions[$i+1]->[1];
				}
			}
		}elsif($nsc_id == $portions[$i]->[0] and $sub_id != $portions[$i]->[1] and $portions[$i]->[16] > -9){
			$psize=$size; #################
			$size += $portions[$i]->[23] == 1 ? $portions[$i]->[9] : $portions[$i]->[15];
			($portions[$i]->[28],$portions[$i]->[29])=($psize,$size) if abs($portions[$i]->[16])==9; #################
			$sub_id = $portions[$i]->[1];					
		}elsif($nsc_id == $portions[$i]->[0] and $sub_id != $portions[$i]->[1] and $portions[$i]->[16]== -9){
			$opt_nsc[$opt_id]->[1]=$i-1;
			$opt_nsc[$opt_id]->[2] = $size;
			$psize=0;     #################
			$size = $portions[$i]->[23] == 1 ? $portions[$i]->[9] : $portions[$i]->[15];
			($portions[$i]->[28],$portions[$i]->[29])=($psize,$size) if abs($portions[$i]->[16])==9; #################
			$opt_id++;
			$opt_nsc[$opt_id]->[0]=$i;
			$sub_id = $portions[$i]->[1];	
		}elsif($nsc_id != $portions[$i]->[0]){
			$opt_nsc[$opt_id]->[1]=$i-1;
			$opt_nsc[$opt_id]->[2]=$size;
			
			$opt_id++;
			$psize=0;
			$size = $portions[$i]->[23] == 1 ? $portions[$i]->[9] : $portions[$i]->[15];
			($portions[$i]->[28],$portions[$i]->[29])=($psize,$size) if abs($portions[$i]->[16])==9; #################
			$opt_nsc[$opt_id]->[0]=$i;
			$nsc_id=$portions[$i]->[0];
			$sub_id=$portions[$i]->[1];					
		}
	}
	$opt_nsc[$opt_id]->[1]=scalar(@portions)-1;
	$opt_nsc[$opt_id]->[2]=$size;
	
	print "Printing out optimal_assembly data (only output scaffold with length >= 500 bp) ...\n";
	$,='';
	@opt_nsc = sort { $b->[2] <=> $a->[2] } @opt_nsc;	
	open(my $optFH,"|gzip -c >$out_dir/optiNewScaffolds.fa.gz") or die "Can't create pipe |gzip -c >$out_dir/optiNewScaffolds.fa.gz !\n";	
	for(my $i=0;$i<scalar(@opt_nsc);$i++){
		my $temp='';
		for(my $j=$opt_nsc[$i]->[0];$j<=$opt_nsc[$i]->[1];$j++){
			$portions[$j]->[22]=$i;
			$temp .= $portions[$j]->[23] == 1 ? $portions[$j]->[20] : $portions[$j]->[21];
		}
		$zero_fill='0' x (7-length($i));
		$sclen=length($temp);
		push @sc_sizes,$sclen;
		$opt_nsc[$i]->[3]=$sclen; # refined-size	
		next if $opt_nsc[$i]->[3]<500;		#### should be < 500			
		print $optFH ">Sc",$zero_fill,$i," size:",$sclen,"bp\n";
		while($sclen>0){ print $optFH substr($temp,0,100,''),"\n"; $sclen-=100;}
	}
	close $optFH;
	
	### output haplome size and scaffold N50 size
	@sc_sizes = sort { $b <=> $a } @sc_sizes;
	foreach $temp (@sc_sizes){
		last if $temp<500;		#### should be < 500
		$haplome_size+=$temp;
	}
	$temp=0;
	for(my $i=0;$i<scalar(@sc_sizes);$i++){
		$temp+=$sc_sizes[$i];
		if($temp>=0.5*$haplome_size){
			$haplome_scaffold_N50_size=$sc_sizes[$i];
			$haplome_scaffold_N50_number=$i+1;
			last;
		}
	}
	
	###
	print "Printing out hm.connections ... \n";
	$,="\t";
	if(-f "$out_dir/hm.connections"){
		unlink "$out_dir/hm.connections.bak";
		rename "$out_dir/hm.connections","$out_dir/hm.connections.bak";		
	}
	open(my $connectionFH,">$out_dir/hm.connections") or die "Can not create connection file $out_dir/hm.connections !\n";
	print $connectionFH "# You may manually change the value of the last column: *connection_manual:\n";
	print $connectionFH "# 0 = no manual setting, >0(<0) = force (not) to use the connection.\n\n";
	print $connectionFH "# After editing or updating, this file may be saved as hm.connections_edited or other names\n";
	print $connectionFH "# in case being overwritten by another invocation of  HM_haploMerger.pl.\n\n";
	print $connectionFH "#new_scaffold_id\tsub_id\tportion_id\tnew_scaffold_name\tnew_scaffold_position_start\tnew_scaffold_position_end\t";
	print $connectionFH "target_old_scaffold_name\ttarget_old_scaffold_postion_start\ttarget_old_scaffold_postion_end\t";
	print $connectionFH "query_old_scaffold_name\tquery_old_scaffold_postion_start\tquery_old_scaffold_postion_end\t";
	print $connectionFH	"connection\tscore\ttsc_Ns%\ttsc_LCs%\tqsc_Ns%\tqsc_LCs%\tconnection_updated\t*connection_manual\n";
	for(my $i=0;$i<scalar(@portions);$i++){
		next if $portions[$i]->[16] !=1 and $portions[$i]->[16] !=9 and $portions[$i]->[16] !=-9;
		$zero_fill='0' x (7-length($portions[$i]->[22]));
		print $connectionFH $portions[$i]->[0], $portions[$i]->[1], $portions[$i]->[2];
		print $connectionFH "\t";
		print $connectionFH "Sc".$zero_fill.$portions[$i]->[22], $portions[$i]->[28], $portions[$i]->[29];
		print $connectionFH "\t";
		print $connectionFH $portions[$i]->[4], $portions[$i]->[7], $portions[$i]->[7]+$portions[$i]->[9];
		print $connectionFH "\t";
		print $connectionFH $portions[$i]->[10], $portions[$i]->[13], $portions[$i]->[13]+$portions[$i]->[15];	
		print $connectionFH "\t";
		print $connectionFH $portions[$i]->[16],$portions[$i]->[17], $portions[$i]->[30],$portions[$i]->[31], $portions[$i]->[32],$portions[$i]->[33], 0, 0;	
		print $connectionFH "\n";							
	}
	close $connectionFH;
		
	$,=' ';
}
