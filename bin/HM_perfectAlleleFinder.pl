#!/usr/bin/perl  

## Version 1.10

## Introduction.


use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 1 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit; }

  This script is to isolated scaffolds which are perfectly covered by
other scaffolds. These scaffolds can be treated as redundant alleles. 
  
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
  
  --startFromRedundantInfo - output redundant alleles based on the file
                ~result/hm.perfect_redundancy_info; you may edit this file
                to include manual adjustment.
  --redundantInfoFile=<file_name> - default as hm.perfect_redundancy_info_edited            
                
  --coverage=N - N% coverage required to claim a redundant allele
                 scaffold (def=85%)
  --aliCoverage=N - N% alignment (gap-excluded) coverage required to claim
                 a redundant allele scaffold (def=85%)                           
               
  --verbose - print more information        
  --Force - force to overwrite existing result files
  --Delete - delete intermediate and temp files (no effect on logfiles) 
  
	Required files:
	  ~result/hm.scaffolds,
		~result/all.tbest.net.gz,
		~/genome.fa.gz.	
	
  Output files include:
    ~result/hm.perfect_redundancy_info, 
    ~result/alleleSet.fa.gz,   # major scaffolds
    ~result/alleleSet_redundant.fa.gz. # scaffolds covered by major scaffolds 
	 
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
my %para = ("--filter" => 0, "--startFromRedundantInfo" => 0, "--redundantInfoFile" => "hm.perfect_redundancy_info_edited",  
						"--Force" => 0, "--Delete" => 0, "--verbose" => 0, "--coverage" => 85, "--aliCoverage" => 85  );
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

#### setting coverage
print "Set minimum coverage to ",$para{"--coverage"},"\n";  
print "Set minimum alignment coverage to ",$para{"--aliCoverage"},"\n";

#### set the over-write flag and deletion flag
print "Set to OVER-WRITING mode!\n" if $para{"--Force"}>0;
print "Set to DELETING mode!\n" if $para{"--Delete"}>0;
print "Set to VERBOSE mode!\n" if $para{"--verbose"}>0;


#########################################
#########################################
#########################################

sub start_from_scratch();
sub start_from_redundant_info();

if($para{"--startFromRedundantInfo"}>0){
	print "\nStart from provided redundant info from ".$para{'--redundantInfoFile'}.".\n";
	start_from_redundant_info();
}else{
	print "\nStart from scratch .\n";
	start_from_scratch();
}	 

sub start_from_scratch(){
	
	#########################################
	#### read in scaffold info
	# @tsc_ids, stores the name and the size of a (target) scaffold, ordered by size, descendent
	#		$tsc_ids[tsc_id]
	#			->[ 0:tname,  1:tsize,
	#				  2:qname,  3:qsize, 4:tstart, 5:tend, 6:qstart, 7:qend, 
	#         8:strand, 9:ali_score, 10:ali_cover, 11:ali_len, 12:NsLCs_5, 13:NsLCs_mid, 14:NsLCs_3, 15:is_redundant(1=yes, 0=no) ]
	my (@tsc_ids, %sc_names);
	
	open(my $scaffoldFH, "<$out_dir/hm.scaffolds") or die "Can not open $out_dir/hm.scaffolds !\n";
	while(<$scaffoldFH>){
		next if m/^#|^\s/;
		chomp;
		my $tt = [ split /\t/ ];
		my $temp = shift @$tt;
		splice @$tt,2,2;
		@{$tt}[2 .. 15]=('',0,  0,0,0,0,  0,0,0,0,  0,0,0,  0);
		$tsc_ids[$temp] = $tt;
	  $sc_names{$tt->[0]}=$temp;
	}
	close $scaffoldFH;
	
	#########################################
	#### read the all.tbest.net file
	
	print "Read in the all.tbest.net file ...\n";
	my $temp;
	my ($is_found); #$is_found=the right top alignment for a scaffold
	my ($tname, $qname, $tsize, $qsize, $tstart, $tend, $qstart, $qend, $strand, $ali_score, $ali_cover, $ali_len);
	
	my $netFH;
	open($netFH, "gunzip -c $out_dir/all.tbest.net.gz |") 
		or open($netFH, "<$out_dir/all.tbest.net") or die "Can not read $out_dir/all.tbest.net(.gz) !\n";
	
	$qname=undef;
	$tname=undef;
	$is_found=0; 
	while(<$netFH>){
	
		### not info line
		next if m/^\s+gap/;
		
		### net line
		if(m/^net\s+([-\.\w]+)\s+/){
			my $temp = $1;
			## if is not the first target
			if (defined $qname) {
				@{$tsc_ids[$sc_names{$tname}]}[2 .. 11]=($qname, $tsc_ids[$sc_names{$qname}]->[1], $tstart, $tend, $qstart, $qend, 
				                                        $strand, $ali_score, $ali_cover, $ali_len );
			}
			
			## initialize a new target scaffold
			$tname = $temp; # the next target scaffold
			($qname, $tstart, $tend, $qstart, $qend, $strand, $ali_score, $ali_cover, $ali_len)=(undef, 0,0,0,0, 0,0,0,0);
			$is_found=0;
			next;
		}
		
		### fill line
		if(m/^(\s+)fill\s+(\d+)\s+(\d+)\s+([-\.\w]+)\s+([+-])\s+(\d+)\s+(\d+)\s+id\s+\d+\s+score\s+(\d+)\s+ali\s+(\d+).+type\s+(\w+)/){
			#    $1           $2      $3        $4        $5        $6     $7                         $8            $9             $10
			if(length($1)<2){ ###### a top fill
				if($3>0.5*$tsc_ids[$sc_names{$tname}]->[1]){ ###### an alignment covers >50% of the scaffold
					$is_found=1;
					$temp = $5 eq '+' ? 1 : -1;
					($qname, $tstart, $tend, $qstart, $qend, $strand, $ali_score, $ali_cover, $ali_len)
						= ($4, $2, $2+$3, $6, $6+$7, $temp, $8, $3, $9);		
				}else{
					$is_found=0;
				}
			}elsif($is_found>0 and ($10 eq 'inv' or $10 eq 'syn')){
				$ali_score+=$8; $ali_len+=$9;
			} 
		}
	}
	# wrap up the last item
	if (defined $qname and defined $tname) {
		@{$tsc_ids[$sc_names{$tname}]}[2 .. 11]=($qname, $tsc_ids[$sc_names{$qname}]->[1], $tstart, $tend, $qstart, $qend,  
		                                        $strand, $ali_score, $ali_cover, $ali_len );
	}
	close $netFH;
	
	#########################################
	#### read sequences and count Ns and LCs
	
	my ($seq,$line,$name,$is_the_end, $ss, $NsLCs);
	open(my $primaryFH,  "| gzip -c >$out_dir/alleleSet.fa.gz") or die "Can't create pipe | gzip -c >$out_dir/alleleSet.fa.gz !\n";
	open(my $redundantFH,"| gzip -c >$out_dir/alleleSet_redundant.fa.gz") or die "Can't create pipe | gzip -c >$out_dir/alleleSet_redundant.fa.gz !\n";
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
		
		## calculate the NsLCs
		$ss=substr($seq, 0, $tsc_ids[$sc_names{$name}]->[4]);
		$NsLCs = $ss=~tr/Nnatgc/xxxxxx/;
		$tsc_ids[$sc_names{$name}]->[12]=$NsLCs;
		
		$ss=substr($seq, $tsc_ids[$sc_names{$name}]->[4], $tsc_ids[$sc_names{$name}]->[5]-$tsc_ids[$sc_names{$name}]->[4]);
		$NsLCs = $ss=~tr/Nnatgc/xxxxxx/;
		$tsc_ids[$sc_names{$name}]->[13]=$NsLCs;		
		
		$ss=substr($seq, $tsc_ids[$sc_names{$name}]->[5], $tsc_ids[$sc_names{$name}]->[1]-$tsc_ids[$sc_names{$name}]->[5]);
		$NsLCs = $ss=~tr/Nnatgc/xxxxxx/;
		$tsc_ids[$sc_names{$name}]->[14]=$NsLCs;
		
		## looking for perfect redundant allele/scaffolds
		## my ($tstart,$tend,$tsize, $qstart,$qend,$qsize, $ali_cover, $ali_len, $strand)
		my $csc=$tsc_ids[$sc_names{$name}];
		($tstart,$tend,$tsize, $qstart,$qend,$qsize, $ali_cover, $ali_len, $strand)=@{$csc}[4,5,1,  6,7,3, 10, 11, 8];
		if($ali_cover/$tsize > $para{"--coverage"}/100 or $ali_len/$tsize > $para{"--aliCoverage"}/100){
			if( $strand>0 and 
					(	($tstart<$qstart and $tsize-$tend <= $qsize-$qend) 
							or ($tstart<=$qstart and $tsize-$tend < $qsize-$qend) 
							or ($tstart==$qstart and $tsize-$tend == $qsize-$qend	and $csc->[0] lt $csc->[2]) ) ){
				$csc->[15]=1;
			}
			#if($strand<0 and $tstart<=$qsize-$qend and $tsize-$tend <= $qstart ){
			if( $strand<0 and 
					(	($tstart < $qsize-$qend and $tsize-$tend <= $qstart) 
							or ($tstart <= $qsize-$qend and $tsize-$tend < $qstart) 
							or ($tstart == $qsize-$qend and $tsize-$tend == $qstart	and $csc->[0] lt $csc->[2]) ) ){			
				$csc->[15]=1;
			}		
		}
				
		## output sequences
		if($tsc_ids[$sc_names{$name}]->[15]>0){
			print $redundantFH ">$name\n";
			$temp=$tsc_ids[$sc_names{$name}]->[1];
			while($temp>0){ print $redundantFH substr($seq,0,100,''),"\n"; $temp-=100; }
		}else{
			print $primaryFH ">$name\n";
			$temp=$tsc_ids[$sc_names{$name}]->[1];
			while($temp>0){ print $primaryFH substr($seq,0,100,''),"\n"; $temp-=100; }	
		}
	}
	close $gzFH;
		
	#########################################
	#### outputing info
	$,="\t";
	
	#### ~result/hm.perfect_redundancy_info
	if(-f "$out_dir/hm.perfect_redundancy_info"){
		unlink "$out_dir/hm.perfect_redundancy_info.bak";
		rename "$out_dir/hm.perfect_redundancy_info","$out_dir/hm.perfect_redundancy_info.bak";
	}
	print "Output the redundancy info (alignment block) to $out_dir/hm.perfect_redundancy_info.\n";
	open(my $rFH,">$out_dir/hm.perfect_redundancy_info") or die "Can not create $out_dir/hm.perfect_redundancy_info! \n";
	print $rFH "#### You may open this file in excel.\n\n";
	print $rFH "#### You may set *redundancy_flag to 1 to label the scaffold as redundant.\n\n";
	print $rFH "#tsc_id\ttsc_name\ttsize\tqname\tqsize\tstart\tend\tqstart\tqend\tstrand\tali_score\tali_cover\tali_len";
	print $rFH "\tNsLCs_terminal5\tNsLCs_mid\tNsLCs_terminal3\t*redundancy_flag(1=redundant, 0=not)\n";

	#########################################
	#### read in scaffold info
	# @tsc_ids, stores the name and the size of a (target) scaffold, ordered by size, descendent
	#		$tsc_ids[tsc_id]
	#			->[ 0:tname,  1:tsize,
	#				  2:qname,  3:qsize, 4:tstart, 5:tend, 6:qstart, 7:qend, 
	#         8:strand, 9:ali_score, 10:ali_cover, 11:ali_len, 12:NsLCs_5, 13:NsLCs_mid, 14:NsLCs_3, 15:is_redundant(1=yes, 0=no) ]


	for(my $i=0;$i<scalar(@tsc_ids);$i++){
		print $rFH "$i\t";
		print $rFH @{$tsc_ids[$i]};
		print $rFH "\n";
	}
	close $rFH;
	
	#### scaffold lists to standard output
	@tsc_ids = sort { $a->[15] <=> $b->[15] or $b->[1] <=> $a->[1] } @tsc_ids;
	print "\n\nPrimary scaffolds listing:\n";
	print "id\tname\tsize\n";
	my $kk=0;
	my $tsc_count=scalar(@tsc_ids);
	for(;$kk<$tsc_count and $tsc_ids[$kk]->[15]<1;$kk++){
		print $kk,$tsc_ids[$kk]->[0],$tsc_ids[$kk]->[1];
		print "\n";
	}
	print "\n\nRedundant scaffolds listing:\n";
	print "id\tname\tsize\n";
	for(;$kk<$tsc_count;$kk++){
		print $kk,$tsc_ids[$kk]->[0],$tsc_ids[$kk]->[1];
		print "\n";
	}
	
	$,="\t";	

}

#########################################
sub	start_from_redundant_info(){

	my $temp;
	
	open(my $rFH,"<",$para{"--redundantInfoFile"}) or die "Can not read redundantInfoFile ! \n";
	my (@tsc_ids,%sc_names);
	while(<$rFH>){
		next if m/^#|^\s/;
		chomp;
		my $tt = [ split /\t/ ];
		my $temp = shift @$tt;
		$tsc_ids[$temp]=$tt;
		$sc_names{$tt->[0]}=$temp;
	}
	close $rFH;
	
	#########################################
	#### read sequences
	
	my ($seq,$line,$name,$is_the_end);
	open(my $primaryFH,  "| gzip -c >$out_dir/alleleSet.fa.gz") or die "Can't create pipe | gzip -c >$out_dir/alleleSet.fa.gz !\n";
	open(my $redundantFH,"| gzip -c >$out_dir/alleleSet_redundant.fa.gz") or die "Can't create pipe | gzip -c >$out_dir/alleleSet_redundant.fa.gz !\n";
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
		
		## output sequences
		if($tsc_ids[$sc_names{$name}]->[15]>0){
			print $redundantFH ">$name\n";
			$temp=$tsc_ids[$sc_names{$name}]->[1];
			while($temp>0){ print $redundantFH substr($seq,0,100,''),"\n"; $temp-=100; }
		}else{
			print $primaryFH ">$name\n";
			$temp=$tsc_ids[$sc_names{$name}]->[1];
			while($temp>0){ print $primaryFH substr($seq,0,100,''),"\n"; $temp-=100; }	
		}
	}
	close $gzFH;
		
	#########################################
	#### outputing info
	$,="\t";
	
	#### scaffold lists to standard output
	@tsc_ids = sort { $a->[15] <=> $b->[15] or $b->[1] <=> $a->[1] } @tsc_ids;
	print "\n\nPrimary scaffolds listing:\n";
	print "id\tname\tsize\n";
	my $kk=0;
	my $tsc_count=scalar(@tsc_ids);
	for(;$kk<$tsc_count and $tsc_ids[$kk]->[15]<1;$kk++){
		print $kk,$tsc_ids[$kk]->[0],$tsc_ids[$kk]->[1];
		print "\n";
	}
	print "\n\nRedundant scaffolds listing:\n";
	print "name\tsize\n";
	for(;$kk<$tsc_count;$kk++){
		print $kk,$tsc_ids[$kk]->[0],$tsc_ids[$kk]->[1];
		print "\n";
	}
	
	$,="\t";	

}



print "Finished preparation for HM_perfectAllelFinder.pl.\n";
print "\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";


