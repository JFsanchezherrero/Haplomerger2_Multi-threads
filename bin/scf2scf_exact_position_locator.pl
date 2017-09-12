#!/usr/bin/perl  

use strict; 
use warnings;
use FindBin;
use lib $FindBin::Bin."/lib";
use Parallel::ForkManager;
use Data::Dumper;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 0 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

   Given that an original genome assembly has been linked by
   hierarchical scaffolder.
   
   This script is to locate the exact position of the scaffolds of
   the original (query) genome assembly in the target genome assembly. 
  
   ./scf2scf_exact_position_locator.pl new/target_genome old/query_genome
    	 
USAGES

print STDERR "\n\n========== Start at " . qx(date);
print STDERR "$0 $Arg_list\n\n";
my $timer = time();

#### ################################
#### reading scf

my %tscfs; my %qscfs;
my $threads =$ARGV[2];
my $file_name = $ARGV[3];
my @files;

#### %tscfs
{
	my $fasta_name=$ARGV[0];
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
		$tscfs{$name}=uc $seq;
	}
	close $inFH;	
	print STDERR "Finished tscf fasta reading.\n";
}

#### %qscfs
{
	my $fasta_name=$ARGV[1];
	my $tmpDir = "./intermediate_results";
	mkdir $tmpDir, 0755;
	my %files;

	## Split given fasta file in as many CPUs as expected and send multiple threads for each
	# Splits fasta file and takes into account to add the whole sequence if it is broken
	
	my $GZ_file = $fasta_name;
	my @name = split(".fa.gz", $fasta_name);
	my $fa_file = $name[0].".fasta";
	if (-f $GZ_file){ 
		system("gunzip -c $GZ_file > $fa_file"); 
	} else { die "File missing\n"; }
	
	my $file_size = -s $fa_file; #To get only size
	my $block = int($file_size/$threads);
		
	open (FH, "<$fa_file") or die "Could not open source file. $!";
	print "\t- Splitting file into blocks of $block characters aprox ...\n";
	my $j = 0; 
	while (1) {
			my $chunk;
			my $block_file = $tmpDir."/".$fa_file."_part-".$j."_tmp.fasta";
			push (@files, $block_file);
			
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

	for (my $j=0; $j < scalar @files; $j++) {	
		my $inFH;
		open($inFH,"<$files[$j]") or die "Can not open $files[$j].\n";
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
			$qscfs{$files[$j]}{$name}=uc $seq;
		}
		close $inFH;
	}
	print STDERR "Finished qscf fasta reading.\n";		
}

#### scf2scf
my %scf2scf; #0-base; {target_sc_id}->[0:start_position(0-base),1:end_position, 2:sign, 3:query_sc_id];

my $pm =  new Parallel::ForkManager($threads); ## Number of subprocesses not equal to CPUs. Each subprocesses will have multiple CPUs if available
$pm->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; print "\n\n** Child process finished with PID $pid and exit code: $exit_code\n\n"; } );
$pm->run_on_start( sub { my ($pid,$ident)=@_; } );
my $count = 0;
foreach my $files (keys %qscfs){
	$count++;	
	print "Checking file: $files\n";
	my $pid = $pm->start($count) and next; print "\nSending child command\n\n";
	my %scf2scf_tmp; 

	my $out_file = $files."_dump.txt";
	open (DUMP, ">$out_file"); 

	foreach my $qsc_id (keys %{ $qscfs{$files} } ){		
		my $q_id = $qscfs{$files}{$qsc_id};
		my $tt=reverse $q_id;
		$tt=~ tr/ACGTacgt/TGCAtgca/;
		my $flag=0;
		foreach my $tsc_id (keys %tscfs){
			if( $tscfs{$tsc_id} =~ m/$q_id/ ){
				print DUMP $names."\t".$-[0]."\t".$+[0]."\t1\t".$qsc_id."\n";
				$flag+=1;
				print STDERR "position $qsc_id in target genome >>>+ $tsc_id; times: $flag.\n";
				last;
			}
			if( $tscfs{$tsc_id} =~ m/$tt/ ){
				print DUMP $names."\t".$-[0]."\t".$+[0]."\t-1\t".$qsc_id."\n";
				$flag+=1;
				print STDERR "position $qsc_id in target genome >>>- $tsc_id; times: $flag.\n";
				last;
			}
		}
		print STDERR "warning: can not locate $qsc_id in target genome.\n" if $flag==0;
		die "\nError: $qsc_id can be mapped to the target genome twice.\n\n" if $flag>1; 
	}
	close (DUMP);
	
	$pm->finish($count); # pass an exit code to finish
}
$pm->wait_all_children; print "\n** All child processes have finished...\n\n";

## read all hash
my $tmp_file = $tmpDir."/coordinates_tmp.txt"
for (my $j=0; $j < scalar @files; $j++) {	
	my $out_file = $files[$j]."_dump.txt";
	system("cat $out_file >> $tmp_file");
}

#### output
print "#new_scaffold_id\tstart\tend\tstrand\told_scaffold_id\n";
system("sort $tmp_file >> $file_name");
print STDERR "\n\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";


__END__
#### scf2scf
foreach my $qsc_id (keys %qscfs){
  my $tt=reverse $qscfs{$qsc_id};
  $tt=~ tr/ACGTacgt/TGCAtgca/;
  my $flag=0;
  foreach my $tsc_id (keys %tscfs){
    if( $tscfs{$tsc_id} =~ m/$qscfs{$qsc_id}/ ){
      push @{$scf2scf{$tsc_id}},[$-[0],$+[0],1,$qsc_id];
      $flag+=1;
      print STDERR "position $qsc_id in target genome >>>+ $tsc_id; times: $flag.\n";
      last;
    }
    if( $tscfs{$tsc_id} =~ m/$tt/ ){
      push @{$scf2scf{$tsc_id}},[$-[0],$+[0],-1,$qsc_id];
      $flag+=1;
      print STDERR "position $qsc_id in target genome >>>- $tsc_id; times: $flag.\n";
      last;
    }
  }
  print STDERR "warning: can not locate $qsc_id in target genome.\n" if $flag==0;
  die "\nError: $qsc_id can be mapped to the target genome twice.\n\n" if $flag>1;
}
