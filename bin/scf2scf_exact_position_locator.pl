#!/usr/bin/perl  

use strict; 
use warnings;

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

my %tscfs;
my %qscfs;
my %scf2scf; #0-base; {target_sc_id}->[0:start_position(0-base),1:end_position, 2:sign, 3:query_sc_id];

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
		
		$qscfs{$name}=uc $seq;
	}
	close $inFH;
	print STDERR "Finished qscf fasta reading.\n";		
}

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

#### output
print "#new_scaffold_id\tstart\tend\tstrand\told_scaffold_id\n";
foreach my $id (sort(keys %scf2scf)){
  my $tt=$scf2scf{$id};
  @$tt = sort { $a->[0] <=> $b->[0] } @$tt;
  for(my $i=0;$i<@$tt;$i++){
    print "$id\t$tt->[$i][0]\t$tt->[$i][1]\t$tt->[$i][2]\t$tt->[$i][3]\n";
  }
}


print STDERR "\n\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";