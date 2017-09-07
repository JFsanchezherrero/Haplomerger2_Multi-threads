#!/usr/bin/perl  

## Version 1.6
## Changes:
## Threads implemented by JFSanchezHerrero

## Version 1.50
## changes:
## Do not correct illegal characters in the fasta sequences, just report errors and die.
## Now strictly require fa.gz input file.
## Read in all sequences of all species into the memory.
## remove the gzMultiFa because no longer needed.
## rewrite the faSize function so that it whould quit because of too many files in the species.seq directory.

## Version 1.00

## Introduction.
## This script runs faSplit, faToNib, faSize and create tba-specific fasta genome seq file.
## By doing so, create a suitable working environment for the following blastz, axtChainNet and tba procedures. 

## Recommand log file: ./_initiation.log

## Assumptions.
## 1)The current/home directory is the working directory;
## 3)If start from the faSplit phase, a fasta, soft-masked genome seq file (e.g. species_name.fa) containing all chr/contigs
##   should be placed in the home directory.  
## 2)Other other executable programs required for blastz/tba can be found in the $PATH;
## 3)All pertinent files, such as fasta files, are in correct format, in proper place, naming and other condition;
## 4)This script is in its very preliminary form, so users of this script take their own risk. 
## 5)the script does not check for argument syntax, the users need to ensure correct syntax. 

## Dependences.
## faSplit, faToNib, faSize, nibSize, nibFrag

## Update
## V0.99. 
## creation of  this script.
## V1.00. 
## 1)faSplit. replace the "faSplit" executable with a customized splitting procedure,
## in which, a new scaning is incoporated: each input fasta sequence is checked for 
## a) non standard alphabatic characters (other than ATCGN/atcgn), once found, those
## characters would be replaced by 'N'; b) non-alphabetic characters, [^a-zA-Z], once
## found, deleted. 
## 2)chromosome/contig name now allows "-" and "." characters.
## 3)gzMultiFa is added to produce a gz-compressed multi-fa file for blastz/lasz.
## 4)nibFrag is added in a module to convert nib file into fa file
## 5)multifa is added to reconstruct a 'species_name.fa' from nib files. This file holds 
##   all sequences in the directory in species_name.seq/. 

use strict; 
use warnings;
use Data::Dumper;
use FindBin;
use lib $FindBin::Bin."/lib";
use Parallel::ForkManager;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 1 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

	This script is to build up the directories and data files 
required by the subquent runing of lastz, chainNet and TBA.

  Requrie excutables from  http://genome.ucsc.edu/admin/jksrc.zip,
including faToNib, faSize and nibFrag. 

usages: (!!!case sensitive)
   --help or no arguments,   show this message. Undefined inputs are ignored.
   

Main procedures:	
   --faSize - the next step after --faToNib, create a size file
      to contain all fa's sizes
      The sizes file's format is <seq_name>\t<size>\n.
   --faSplit - break down a big multiple-fa file into a directory
      of single-fa file. File name is the sequence name.
   --faToNib - the next step after --faSplit, convert single-fa
      file into nib file, retaining the softmasking
   --tbaFaCr - the next step after the --faSize, create a multiple-fa file
      from all single-fa file
      The multiple-fa file name is the species name.
      This file cotains special headers for sequences which can be recognized
      and utilized by TBA.

Other procedures:	
   --nibFrag - to reconstruct all single-fa files from all nib files 
      NOTES: a correct fa sizes file is required.
   --multiFa - to convert all nib files into a a multiple-fa file.\
      This is used to reconstruct the 'species_name.fa' from nib files.
      NOTES: a correct fa sizes file is required. 

Options:
   --threads=N - specified N threads/CPUs for use (default=1)   
   --Species - mandatory, to provide species names which are wanted to be 
      processed
      Example: --Species human chimp rhesus
      NOTE that all sequences are kept in the memory, so it is
      better to process one species at one time.
            
   --Force - over-writing existing files is allowed (default=not)
   --Delete - use rm -fr to delete the dir species.seq.fa (default=not)

Notes:
1) The current directory is the working directory.
2) If start from the faSplit phase, a fasta, soft-masked genome seq file
   (e.g. species_name.fa.(gz)) containing all chr/contigs should be placed in
   the working directory.  	
3) Five main procedures can be run separately or together, but since each 
   procedure is based on the previous, you should run them step by step.
	 
USAGES

print "\n\n========== Start at "; system("date");
print "$0 $Arg_list\n\n";
my $timer = time();

# set the over-write flag, to over-write any existed files
my $Force=0;
if ($Arg_list =~ m/--Force/) {
	$Force = 1; print "Set to OVER-WRITING mode!\n";
}
#set the output_field_separator to " "
$,=' ';

# Store species names
my @Species;
unless ($Arg_list =~ m/--Species\s+([^-]*)/) {die "No --Species argument or species_names found!\n" }
unless (@Species = $1 =~ m/(\w+)/g)  {die "No species names found!\n" };
print "Species included: ".join(",",@Species)."\n";

#my $file_name;
#unless ($Arg_list =~ m/--Species\s+(\S+)/) {die "No --Species argument or species_names found!\n" }
#$file_name = $1;
#print "Species included: ", $file_name, "\n";

my $threads = 1;
if ($Arg_list =~ m/--threads=(\d+)/){ $threads = $1 > 0 ? $1 : 1; }
print "Thread number is set to ... $threads !\n";

my $tmpDir = "./intermediate_results";
mkdir $tmpDir, 0755;
my %Species_files;
for (my $i=0; $i < scalar @Species; $i++) {
	
	print "Reading: ".$Species[$i]."\n";
	
	## Split given fasta file in as many CPUs as expected and send multiple threads for each
	# Splits fasta file and takes into account to add the whole sequence if it is broken

	my $GZ_file = $Species[$i].".fa.gz";
  	my $fa_file = $Species[$i].".fa";
  	if (-f $GZ_file){ system("gunzip $GZ_file"); } elsif (-f $fa_file) {} else {die "File missing\n";}
  	
	my $file_size = -s $fa_file; #To get only size
	my $block = int($file_size/$threads);
	
	open (FH, "<$fa_file") or die "Could not open source file. $!";
	print "\t- Splitting file into blocks of $block characters aprox ...\n";
	my $j = 0; 
	while (1) {
			my $chunk;
			my $block_file = $tmpDir."/".$Species[$i]."_part-".$j."_tmp.fasta";
			push (@{ $Species_files{$Species[$i]} }, $block_file);
			
			open(OUT, ">$block_file") or die "Could not open destination file";
			if (!eof(FH)) { read(FH, $chunk,$block);  
					if ($j > 0) { $chunk = ">".$chunk; }
					print OUT $chunk;
			} ## Print the amount of chars  
			if (!eof(FH)) { $chunk = <FH>; print OUT $chunk; } ## print the whole line if it is broken      
			if (!eof(FH)) { 
					$/ = ">"; ## Telling perl where a new line starts
					$chunk = <FH>; chop $chunk; print OUT $chunk; 
					$/ = "\n";
			} ## print the sequence if it is broken
			$j++; close(OUT); last if eof(FH);
	}
	close(FH);
}
## print Dumper \%Species_files;

# read in the fasta file
my %fasta; #{species}->{name}=seq
foreach my $keys (keys %Species_files) {
  
  my @files = @{ $Species_files{$keys} };
  for (my $g=0; $g < scalar @files; $g++) {  
	  my $inFile = $files[$g];
	  print $inFile."\n";
	  $fasta{$keys}{ $inFile }={};
	  my $inFH;
	  unless(-f $inFile){ die "Can not find $inFile file!\n"; }
	  open($inFH, "<$inFile") or die "Can not open $inFile file!\n";
	  
	  my ($line,$name,$seq,$is_the_end)=("","","",0);
	  while($line=<$inFH>){ last if $line =~ m/^>/;  }	
	  die "Problems with the fasta file; no symbol > found !\n" unless defined($line);  
	
	  while($is_the_end==0){
		($name,$seq)=('','');
		if($line =~ m/^>(\S+)/){ $name=$1; }
		if(length($name)<1){ die "Species ($keys) - $inFile: fasta name format error.\n"; }
	  
		while($line=<$inFH>){
			last if $line =~ m/^>/;
			chomp $line;
			$seq .= $line;
		}
		$is_the_end=1 unless defined($line);
		if($seq=~m/[^acgtACGTnN]/){ die "File ( $inFile ) has a sequence ( $name ) containing illegal characters !\n"; }
		if(!defined($fasta{$keys}{$inFile}{$name})){
		  $fasta{$keys}{$inFile}{$name}=$seq;
		}else{
		  die "File ( $keys ) has a duplicate seq name ($name  ) !\n";
	  	}
	}
	close $inFH;
}}

sub faSplit;
sub faToNib;
sub faSize;
sub tbaFaCr;	# tbaFasta file creation! 
sub nibFrag;
sub multiFa;

faSize  if ($Arg_list =~ m/--faSize/);
faSplit if ($Arg_list =~ m/--faSplit/);
faToNib if ($Arg_list =~ m/--faToNib/);
tbaFaCr if ($Arg_list =~ m/--tbaFaCr/);
nibFrag if ($Arg_list =~ m/--nibFrag/);
multiFa if ($Arg_list =~ m/--multiFa/);

if ($Arg_list =~ m/--Delete/) { foreach my $temp (@Species) { system("rm -fr $temp.seq.fa");}}
print "\n\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";

############################## All subroutines ####################################

####
sub faSize {
	print "====faSize started!====\n";
	
	##Step1: check for existing sizes file, if existes, die or not depends on the $Force
	print "checking existing sizes files ...\n";
	
	##Step2: create the sizes files	
	## Parallelize
	foreach my $keys (keys %Species_files) {
		my $pm =  new Parallel::ForkManager($threads); ## Number of subprocesses not equal to CPUs. Each subprocesses will have multiple CPUs if available
		$pm->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; print "\n\n** Child process finished with PID $pid and exit code: $exit_code\n\n"; } );
		$pm->run_on_start( sub { my ($pid,$ident)=@_; } );

		my @files2check = @{ $Species_files{$keys} };
		for (my $j = 0; $j < scalar @files2check; $j++) {
			my $pid = $pm->start($j) and next; print "\nSending child command\n\n";
		
			my $tmp_sizeFile = $files2check[$j].".sizes";
			print "Printing into: $tmp_sizeFile\n";
			my $sizeFH;
			open($sizeFH, ">$tmp_sizeFile") or die "Can not open $tmp_sizeFile!\n";	  
    		foreach my $temp2 (keys %{ $fasta{$keys}{$files2check[$j]} } ) {
				print $sizeFH "$temp2\t".length( $fasta{$keys}{$files2check[$j]}{$temp2})."\n";
			} close $sizeFH;

			$pm->finish($j); # pass an exit code to finish
		}
		$pm->wait_all_children; print "\n** All child processes have finished...\n\n";
		system ("cat $tmpDir/$keys*sizes >> $keys.sizes");
		#system ("rm $tmpDir/$keys*sizes");
	}
	
	print "====faSize done!====\n\n";	
}

####
sub faSplit {
	print "====faSplit started!====\n";
	
	#checking the existing directory
	print "Checking existing directories ... \n";
	
	#spliting files
	foreach my $keys (keys %Species_files) {
		mkdir("$keys.seq.fa");
	
		my $pm =  new Parallel::ForkManager($threads); ## Number of subprocesses not equal to CPUs. Each subprocesses will have multiple CPUs if available
		$pm->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; print "\n\n** Child process finished with PID $pid and exit code: $exit_code\n\n"; } );
		$pm->run_on_start( sub { my ($pid,$ident)=@_; } );

		my @files2check = @{ $Species_files{$keys}};
		for (my $j = 0; $j < scalar @files2check; $j++) {
			my $pid = $pm->start($j) and next; print "\nSending child command\n\n";
		
    		foreach my $temp2 (keys %{ $fasta{$keys}{$files2check[$j]} } ) {
				my $sfaFH;
				open($sfaFH, ">$keys.seq.fa/$temp2.fa") or die "Can not open $keys.seq.fa/$temp2.fa!\n";
				my $seq = $fasta{ $keys }{ $files2check[$j]}{$temp2};
				print $sfaFH ">$temp2\n$seq\n";
				close $sfaFH;
			} 
			$pm->finish($j); # pass an exit code to finish
			print "$files2check[$j] has been added to into $keys.seq.fa/\n";
		}
		$pm->wait_all_children; print "\n** All child processes have finished...\n\n";
	}
	# after this, the fasta sequences are no longer needed to be kept in the memory.
	print "====faSplit done!====\n\n";
}

####
sub faToNib {
	my ($temp1, $temp2, $has_nib, $has_fa, $cmd);
	print "====faToNib started!====\n";
	
	#spliting files
	foreach my $keys (keys %Species_files) {
		
		#Step1: checking the existing directory and make it if not existed
		print "Checking existing directories ... \n";
		mkdir("$keys.seq");
	
		##Step3: converting the fa files to nib files
		print "converting the fa files to nib files ...\n";
		print "Converting directory $keys.seq.fa/*.fa to $keys.seq/*.nib ...\n";
	
		my $pm =  new Parallel::ForkManager($threads); ## Number of subprocesses not equal to CPUs. Each subprocesses will have multiple CPUs if available
		$pm->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; print "\n\n** Child process finished with PID $pid and exit code: $exit_code\n\n"; } );
		$pm->run_on_start( sub { my ($pid,$ident)=@_; } );

		my @files2check = @{ $Species_files{$keys}};
		for (my $j = 0; $j < scalar @files2check; $j++) {
			my $pid = $pm->start($j) and next; print "\nSending child command\n\n";
    		foreach my $temp2 (keys %{ $fasta{$keys}{$files2check[$j]} } ) {
				my $rv=system("faToNib -softMask $keys.seq.fa/$temp2.fa $keys.seq/$temp2.nib");
				#print "faToNib -softMask $keys.seq.fa/$temp2.fa $keys.seq/$temp2.nib\n";
				die "Can not find $keys.seq.fa/$temp2.fa !\n" if $rv>0;
			} 
			$pm->finish($j); # pass an exit code to finish
			print "$files2check[$j] has been added to into $keys.seq.fa/\n";
		}
		$pm->wait_all_children; print "\n** All child processes have finished...\n\n";
	}
	print "====faToNib done!====\n\n";
}

=head1
####
sub tbaFaCr {
	my ($temp1, $temp2, $missing_fa, $has_tbaFa, $sizeFH, $faFH, $tbaFH, $size, $cmd);
	print "====tbaFaCr started!====\n";
	
	##Step1: check if the directory exists, if not exist, die!
	print "checking if the directory is existed ... \n";
	foreach $temp1 (@Species) { die "WARNING: ./$temp1.seq.fa/ for $temp1 is not existed! Die!\n" unless (-d "$temp1.seq.fa") }
	
	##Step2: check for existing sizes file, if not exist, die!
	print "checking existing sizes files ... \n";
	foreach $temp1 (@Species) { 
		die "$temp1.sizes does not exist! Die!\n" unless (-f "$temp1.sizes");		
	}
	
	##Step3: check if all seqs from sizes file are present in the directory, if not, die!
	print "check if all seqs from sizes file are present in the directory ...\n";
	$missing_fa=0;
	foreach $temp1 (@Species) {
		print "Checking missing fa file in $temp1.seq.fa ...\n";
		open($sizeFH, "<$temp1.sizes") or die "Cant not open directory $temp1.sizes!\n";
		while (<$sizeFH>) {
			next unless m/([-\.\w]+)\t/;
			$temp2=$1;
			unless (-f "$temp1.seq.fa/$temp2.fa") { $missing_fa=1; print "$temp2.fa\n"; }
		}
		close $sizeFH;
	}
	die "Some fa files are missing! Die!\n" if ($missing_fa == 1);
	
	##Step4: check for existing tba-fasta file, if exists, die or not depends on $Force!
	print "checking tba-sepcific fasta files ... \n";
	$has_tbaFa=0;
	foreach $temp1 (@Species) { 
		if (-f "$temp1") { $has_tbaFa = 1; print "$temp1\n"; }		
	}
	die   "tba-fa files found! Die!\n" if ($has_tbaFa==1 and $Force == 0);
	print "tba-fa files found! forced to going on ... \n" if ($has_tbaFa==1 and $Force == 1);
	
	##Step5: create the tba-specific fasta files for each species
	foreach $temp1 (@Species) {
		open($tbaFH, ">$temp1") or die "Can not open $temp1\n";
		open($sizeFH, "<$temp1.sizes") or die "Can not open $temp1.sizes\n";
		while (<$sizeFH>) {
			next unless m/([-\.\w]+)\t(\d+)/;
			$temp2 = $1; 	#chromosome/contig name
			$size  = $2; 	#seq size/length
			print $tbaFH ">$temp1:$temp2:1:+:$size\n"; 
			open($faFH, "<$temp1.seq.fa/$temp2.fa") or die "Can not open $temp1.seq.fa/$temp2.fa!\n";
			while (<$faFH>) { next if m/>/; print $tbaFH $_; } 	#copying sequence
			close $faFH;
		}
		print "tba-fa file $temp1 done!\n";
		close $sizeFH; close $tbaFH;
	}
	print "====tbaFaCr done!====\n\n";
}

####
sub nibFrag {
	my ($temp1, $temp2, $missing_nib, $has_fa, $cmd, $sizeFH, $size);
	print "====nibFrag started!====\n";
	
	##Step1: to test if the directory is existed, if not existed, die!
	foreach $temp1 (@Species) { die "WARNING: ./$temp1.seq/ for $temp1 is not existed! Die!\n" unless (-d "$temp1.seq") }
	
	##Step2: check for existing sizes file, if not exist, die!
	print "checking existing sizes files ... \n";
	foreach $temp1 (@Species) { die "$temp1.sizes does not exist! Die!\n" unless (-f "$temp1.sizes");	}	
	
	##Step3: check if all seqs from sizes file are present in nib format file, if not, die!
	print "check if every seq from sizes file has a corresponding nib file ...\n";
	$missing_nib=0;
	foreach $temp1 (@Species) {
		print "Checking missing nib file in $temp1.seq ...\n";
		open($sizeFH, "<$temp1.sizes") or die "Cant not open directory $temp1.sizes!\n";
		while (<$sizeFH>) {
			next unless m/([-\.\w]+)\t/;
			$temp2=$1;
			unless (-f "$temp1.seq/$temp2.nib") { $missing_nib=1; print "$temp2.nib\n"; }
		}
		close $sizeFH;
	}
	die "Some nib files are missing! Die!\n" if ($missing_nib == 1);
	

	
	##Step4: check if any seqs from sizes file has single-fa format file, if not, die depending on $Force!
	print "check if any seq from sizes file has corresponding single-fa files ...\n";
	$has_fa=0;
	foreach $temp1 (@Species) {
		print "Checking single-fa files in $temp1.seq ...\n";
		open($sizeFH, "<$temp1.sizes") or die "Cant not open directory $temp1.sizes!\n";
		while (<$sizeFH>) {
			next unless m/([-\.\w]+)\t/;
			$temp2=$1;
			if (-f "$temp1.seq/$temp2.fa") { $has_fa=1; print "$temp2.fa\n"; }
		}
		close $sizeFH;
	}
	die "Some single-fa files found! Die!\n" if ($has_fa == 1 and $Force == 0);
	print "Some single-fa files found, will be over-written!\n" if ($has_fa==1 and $Force == 1);
	
	##Step5: converting the nib files to single-fa files
	print "converting the nib files to singl-fa files ...\n";
	foreach $temp1 (@Species) {
		print "Processing directory $temp1.seq ...\n";
		open($sizeFH, "<$temp1.sizes") or die "Cant not open directory $temp1.sizes!\n";
		unless (-d "$temp1.seq.fa") { mkdir "$temp1.seq.fa";}
		while (<$sizeFH>) {
			next unless m/([-\.\w]+)\t(\d+)/;
			($temp2, $size) = ($1,$2);
			$cmd = "nibFrag -masked -name=$temp2 $temp1.seq/$temp2.nib 0 $size + $temp1.seq.fa/$temp2.fa";
			system($cmd);
			print "$temp1.seq/$temp2.nib has been converted to $temp1.seq.fa/$temp2.fa\n";
		}
		close $sizeFH;
	}
	print "====nibFrag done!====\n\n";	
}

####

sub multiFa {
	my ($temp1, $temp2, $missing_nib, $has_fa, $cmd, $sizeFH, $size, $file_created);
	print "====multiFa started!====\n";
	
	##Step1: to test if the directory is existed, if not existed, die!
	foreach $temp1 (@Species) { die "WARNING: ./$temp1.seq/ for $temp1 is not existed! Die!\n" unless (-d "$temp1.seq") }
	
	##Step2: check for existing sizes file, if not exist, die!
	print "checking existing sizes files ... \n";
	foreach $temp1 (@Species) { die "$temp1.sizes does not exist! Die!\n" unless (-f "$temp1.sizes");	}	
	
	##Step3: check if all seqs from sizes file are present in nib format file, if not, die!
	print "check if every seq from sizes file has a corresponding nib file ...\n";
	$missing_nib=0;
	foreach $temp1 (@Species) {
		print "Checking missing nib file in $temp1.seq ...\n";
		open($sizeFH, "<$temp1.sizes") or die "Cant not open directory $temp1.sizes!\n";
		while (<$sizeFH>) {
			next unless m/([-\.\w]+)\t/;
			$temp2=$1;
			unless (-f "$temp1.seq/$temp2.nib") { $missing_nib=1; print "$temp2.nib\n"; }
		}
		close $sizeFH;
	}
	die "Some nib files are missing! Die!\n" if ($missing_nib == 1);
	
	##Step4: check for existing multi-fa file, if exists, die depending on $Force!
	print "checking existing multi-fa files ... \n";
	$has_fa = 0;
	foreach $temp1 (@Species) { 
		if (-f "$temp1.fa") { $has_fa=1; print "$temp1.fa\n"; 	}
	}
	die "Multiple fa file found! Die!\n" if ($has_fa==1 and $Force==0);
	print "Multiple fa file found, will be over-written!\n" if ($has_fa==1 and $Force==1);	
	
	##Step5: converting the nib files to multiple-fa files
	print "converting the nib files to multiple-fa files ...\n";
	foreach $temp1 (@Species) {
		print "Processing directory $temp1.seq ...\n";
		open($sizeFH, "<$temp1.sizes") or die "Cant not open directory $temp1.sizes!\n";
		$file_created = 0;
		while (<$sizeFH>) {
			next unless m/([-\.\w]+)\t(\d+)/;
			($temp2, $size) = ($1,$2);
			if ($file_created == 0) {
				$file_created = 1;
				$cmd = "nibFrag -masked -name=$temp2 $temp1.seq/$temp2.nib 0 $size + /dev/stdout | cat >$temp1.fa";
			}else{
				$cmd = "nibFrag -masked -name=$temp2 $temp1.seq/$temp2.nib 0 $size + /dev/stdout | cat >>$temp1.fa";
			}
			system($cmd);
			print "$temp1.seq/$temp2.nib has been converted and added to $temp1.fa\n";
		}
		close $sizeFH;
	}
	print "====multiFa done!====\n\n";	
}
=tail

