#!/usr/bin/perl  

## Version 2.40

## Introduction.

## bug fixing
## Now explicitly use perl to run each perl script.

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 1 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit; }

  This script is to refine the unpaired.fa.gz by deleting 
    1) heavily soft-masked sequences,
    2) sequences containing too many N/n,
    3) redundant sequences.
  A file called unpaired_refined.fa.gz will be created.  
  
  Require unpaired_raw.fa.gz file and 
  unpaired_raw.optimalContinuity.result/all.tbest.net file.

  --help or no arguments - show this message, ignore undefined inputs
  
  --Species -   mandatory, to provide species names which are wanted to be
                processed. NOTE THAT the species name should be given
                in ORDER, because the name comes first will be used
                as the target.
                Example: --Species bbev1a bbev1b
                In this script, both names refer to the same genome assembly,
                the name come first is used as the target genome, the second
                name is used as the query genome. 
  
  --runLastzChainNet - run lastz and chainNet first (default=no)
  --threads=N - use N cpu/threads for --runLastzChainNet (default=1)
  --identity=N - identity threshold for HM_all_lastz_mThreads.pl (def=50)               

  --maskFilter=<XX>, filter sequences in which XX% are soft-masked and
           Ns (default=85, recommand 80, 85,90)
  --redundantFilter=<XX>, filter sequences in which XX% are redundant,
           convered by alignment_length+Ns (default=85, recommand 80,85,90)          
     
  --Force, force to overwrite existed result files.
  --Delete, to delete intermediate and temp files (no effect on logfiles).        

	 
USAGES

print "\n\n========== Start at "; system("date");
print "$0 $Arg_list\n\n";
my $timer = time();

#### set the output_field_separator to " " and autoflushing
$,=' ';
$|=1;

#### parse the arguments
my %para = ("--maskFilter" => 85,  "--redundantFilter" => 85,	"--Force" => 0, "--Delete" => 0, "--runLastzChainNet" => 0,
						"--threads"=>1, "--identity"=>50  );
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

# Store species names and result directory name
my (@Species, $out_dir); 
unless ($Arg_list =~ m/--Species\s+([^-]*)/) {die "No --Species argument or species_names found!\n"; }
unless (@Species = $1 =~ m/(\w+)/g)  {die "No species names found!\n"; }
$out_dir = "$Species[0].$Species[1].result";
print "Species included: ", @Species, "\n";

#### set the over-write flag and deletion flag
print "Set to OVER-WRITING mode!\n" if $para{"--Force"}>0;
print "Set to DELETING mode!\n" if $para{"--Delete"}>0;
print "--maskFilter is set to ",$para{"--maskFilter"},"% !\n";
print "--redundantFilter is set to ",$para{"--redundantFilter"},"% !\n";

#### definition of cmd (command line string)
my $cmd;

### run lastz and chainNet
if($para{"--runLastzChainNet"}>0){
	
	unlink "unpaired.fa.gz", "optiNewScaffolds.fa.gz";
	system("ln -s $out_dir/unpaired.fa.gz unpaired.fa.gz");
	system("ln -s $out_dir/optiNewScaffolds.fa.gz optiNewScaffolds.fa.gz");
	
	unlink glob "unpaired.seq/*.nib";
	unlink glob "optiNewScaffolds.seq/*.nib";
	
	$cmd = "initiation.pl --faSplit --faToNib --faSize --Species unpaired optiNewScaffolds --Force --Delete";
	$cmd.= " 1>_hm.un_initiation.log 2>>_hm.un_initiation.log";
	system($cmd);
	print "Finish ",$cmd,"\n";
	
	$cmd = "HM_all_lastz_mThreads.pl --Species unpaired optiNewScaffolds  --threads=" . $para{"--threads"};
	$cmd.= " --targetSize=10000000 --identity=" . $para{'--identity'};
	$cmd.= " --Force --Delete 1>_hm.un_all_lastz.log 2>>_hm.un_all_lastz.log";
	system($cmd);
	print "Finish ",$cmd,"\n"; 

	$cmd = "HM_axtChainRecipBestNet.pl --rbestNet --axtChain --tbest --threads=" . $para{"--threads"};
	$cmd.= " --linearGap=medium --minScore=10000 --minSpace=200 --minScore2=20000 --Species unpaired optiNewScaffolds ";
	$cmd.= " --Force --Delete 1>_hm.un_axtChainRecipBestNet.log 2>>_hm.un_axtChainRecipBestNet.log";
	system($cmd);
	print "Finish ",$cmd,"\n";	
	
	system("cp unpaired.optiNewScaffolds.result/all.tbest.net.gz $out_dir/unpaired.tbest.net.gz");
	unlink "unpaired.fa.gz", "optiNewScaffolds.fa.gz";
	unlink "unpaired.sizes", "optiNewScaffolds.sizes";
}

#########################################
my $net_name="$out_dir/unpaired.tbest.net.gz";
my $seq_name="$out_dir/unpaired.fa.gz";
unless (-f $net_name and -f $seq_name)  {die "$net_name or $seq_name files are not found!\n"; }
#$cmd = "netFilter -minScore=20000 $net_name >unpaired.optiNewScaffolds.result/filtered.tbest.net";
#print $cmd,"...\n";
#system($cmd);
#$net_name="$out_dir/unpaired.tbest.net";

#########################################
#### unincorps[id]->[ 0:new_scaffold_name, 1:old_scaffold_name, 2:size, 3:start, 4:len, 
####                  5:full_scaffold_or_part_of_the_scaffold, 6:N_counts, 7:lowcase_count, 8:ali_len, 9:ali_converage, 10:is_retained ]
print "Read in $out_dir/hm.unpaired_updated ...\n";
my $info_name="$out_dir/hm.unpaired_updated";
my (@unincorps,%names);
open(my $un_inFH,"<$info_name") 
	or die "Can't open $info_name, to recreate this file please run the script with --runLastzChainNet=1 !\n";
while(<$un_inFH>){
	next if m/^#|^\s/;
	chomp;
	my $tt = [ split /\t/ ];
	push @$tt,(0,0,1);
	push @unincorps,$tt;
}
close $un_inFH;
@unincorps = sort { $a->[0] cmp $b->[0] } @unincorps;
for(my $i=0;$i<scalar(@unincorps);$i++){
	$names{$unincorps[$i]->[0]}=$i;
}

#########################################
#### read the filtered.tbest.net file and produce nodes' data structure
print "Read in the $net_name file ...\n";
my ($temp, $cname, $is_top, $ali_cover, $ali_len);

my $netFH;
open($netFH, "gunzip -c $net_name |") or open($netFH, "<$net_name") or die "Can not read $net_name !\n";
$is_top=0;
while(<$netFH>){

	### not info line
	next if m/^\s+gap/;
	
	### net line
	if(m/^net\s+([-\.\w]+)\s+/){
		my $temp = $1;
		## if is not the first target
		if (defined($cname)) {
			$unincorps[$names{$cname}]->[8]= $ali_len;
			$unincorps[$names{$cname}]->[9]= $ali_cover;
		}
		
		## initialize a new target scaffold
		$ali_cover=0; $ali_len=0; $is_top=0;
		$cname = $temp; # the next target scaffold
		next;
	}
	
	### fill line
	if(m/^(\s+)fill\s+\d+\s+(\d+)\s+[-\.\w]+\s+[+-]\s+\d+\s+\d+\s+id\s+\d+\s+score\s+\d+\s+ali\s+(\d+).+type\s+(\w+)/){
		if(length($1)<2){ ###### a top fill
			##if($2>50000 or $2>=0.45*$unincorps[$names{$cname}]->[4] 
			##		or $3>25000 or $3>=0.33*($unincorps[$names{$cname}]->[4]-$unincorps[$names{$cname}]->[6])){
			if( $3>40000  
					or ( ($unincorps[$names{$cname}]->[4]-$unincorps[$names{$cname}]->[6])>10000 
									and $3>0.34*($unincorps[$names{$cname}]->[4]-$unincorps[$names{$cname}]->[6]) )
					or ( $3>0.5*($unincorps[$names{$cname}]->[4]-$unincorps[$names{$cname}]->[6]) ) ){						
				$is_top=1; $ali_cover+=$2; $ali_len+=$3;	
			}else{
				$is_top=0;
			}
		}elsif($is_top>0 and ($4 eq 'inv' or $4 eq 'syn')){
			$ali_len+=$3;
		} 
	}
}
$unincorps[$names{$cname}]->[8]= $ali_len;
$unincorps[$names{$cname}]->[9]= $ali_cover;
close $netFH;

#########################################
for(my $i=0;$i<scalar(@unincorps);$i++){
	
	if($unincorps[$i]->[4]<500){ $unincorps[$i]->[10]=0; next; }
	if($unincorps[$i]->[4]-$unincorps[$i]->[6]-$unincorps[$i]->[7] <500){ $unincorps[$i]->[10]=0; next; }
	if($unincorps[$i]->[4]-$unincorps[$i]->[6]-$unincorps[$i]->[8] <500){ $unincorps[$i]->[10]=0; next; }	
	#if($unincorps[$i]->[4]-$unincorps[$i]->[9] <200){ $unincorps[$i]->[10]=0; next; }
	
	if($unincorps[$i]->[6]+$unincorps[$i]->[7] >= $para{'--maskFilter'}*$unincorps[$i]->[4]/100){ $unincorps[$i]->[10]=0; next; }
	if($unincorps[$i]->[6]+$unincorps[$i]->[8] >= $para{'--redundantFilter'}*$unincorps[$i]->[4]/100){ $unincorps[$i]->[10]=0; next; }
	#if($unincorps[$i]->[9] >= $para{'--redundantFilter'}*$unincorps[$i]->[4]/100){ $unincorps[$i]->[10]=0; next; }
}

#########################################
print "Output sequence to unpaired_refined.fa.gz ...\n";
open(my $inFH,"gunzip -c $out_dir/unpaired.fa.gz |") or die "Can't create pipe gunzip -c $out_dir/unpaired.fa.gz | !\n";
open(my $outFH,"| gzip -c > $out_dir/unpaired_refined.fa.gz") or die "Can't create pipe | gzip -c > $out_dir/unpaired_refined.fa.gz !\n";
my ($line,$seq,$name,$desc,$is_the_end);

#### look for the first seq name
while($line=<$inFH>){
	next if $line =~ m/^\s+/;
	last if $line =~ m/^>/;
}	
die "Problems with unpaired.fa.gz on line ... $line \n" if !defined($line);

#### store the seq
$is_the_end=0;
while($is_the_end==0){
	chomp $line;
	$line =~ /^>([-\w\d\.]+)\s+(.+)/;
	($name,$desc)=($1,$2);
	$seq ='';
	
	while($line=<$inFH>){
		next if $line =~ m/^\s+/;
		last if $line =~ m/^>/;
		$seq .= $line;
		chomp $line;
	}
	$is_the_end=1 if !defined($line);
		
	if($unincorps[$names{$name}]->[10]>0){
		print $outFH ">$name $desc\n";
		print $outFH $seq;
	}
}
close $inFH;
close $outFH;
	
#########################################
# calculate the sequence size
my ($size_of_unpaired_scaffolds, $size_of_unpaired_scaffold_portions)=(0,0);
foreach my $temp (@unincorps){
	$size_of_unpaired_scaffolds += $temp->[4] if $temp->[10] >0 and $temp->[5] eq 'full';
	$size_of_unpaired_scaffold_portions += $temp->[4] if $temp->[10] >0 and $temp->[5] eq 'part';
}

# printing info
$,="\t";
print "\nThe size of refined unpaired scaffolds is $size_of_unpaired_scaffolds bp;\n";
print "the size of refined unpaired scaffold portions is $size_of_unpaired_scaffold_portions bp.\n"; 

# printing table
print "\n\n";
print "# This table contains REFINED information of the scaffolds and scaffold portions that failed to be incorporated\n";
print "# into the new scaffold assebmly in the path_finding process.\n";
print "# and only sequences with length >=500bp were output. \n\n";
print "#new_scaffold_name\told_scaffold_name\tsize\tstart\tlen\tfull_scaffold_or_part_of_the_scaffold\tN_counts\tlowcase_count\tali_len\tali_coverage\tis_retained\n";
for(my $i=0;$i<scalar(@unincorps);$i++){
	my $un=$unincorps[$i];
	print @$un,"\n";
}

if($para{"--Delete"}>0){
	print "Delete optiNewScaffolds.seq, unpaired.seq, unincorporated.optiNewScaffolds.result.\n";
	system("rm -f -r unpaired.seq optiNewScaffolds.seq unpaired.optiNewScaffolds.result");
}

print "\n\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";

