#!/usr/bin/perl  

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 2 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

   This script is to recreate an new assembly based on
   the new_scaffold versus old_scaffold relation 
   (target2query_links) information.
  
   ./scf2scf_assembly_recreation.pl old_genome/query new2old_scaffolds.txt
    	 
USAGES

print STDERR "\n\n========== Start at " . qx(date);
print STDERR "$0 $Arg_list\n\n";
my $timer = time();

#### ################################
#### reading scf

my %tscfs;
my %qscfs;
my %scf2scf; #0-base; {target_sc_id}->[0:start_position(0-base),1:end_position, 2:sign, 3:query_sc_id];

#### %qscfs
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
		
		$qscfs{$name}=$seq;
	}
	close $inFH;
	print STDERR "Finished qscf fasta reading.\n";		
}

#### read in scf2scf
{
  open(my $inFH,"<$ARGV[1]") or die "can not open $ARGV[1].\n";
  while(<$inFH>){
  	next if m/^(#|\s)/;
    chomp;
    next unless m/\S+\t\S+\t\S+\t\S+/;
    my $tt=[split /\t/,$_];
    push @{$scf2scf{$tt->[0]}},[$tt->[1],$tt->[2],$tt->[3],$tt->[4]];
  }
  close $inFH; 
  print STDERR "Finished reading scf2scf_links.\n";		
}

#### scf2scf
my ($Ns);
foreach my $tsc_id (keys %scf2scf){
  my $tt=$scf2scf{$tsc_id};
  @$tt = sort { $a->[0] <=> $b->[0] } @$tt;
  my $i=0;
  for($i=0;$i<scalar(@$tt)-1;$i++){
    if($tt->[$i][2]<0){
      $qscfs{$tt->[$i][3]}=reverse $qscfs{$tt->[$i][3]};
      $qscfs{$tt->[$i][3]}=~ tr/ACGTacgt/TGCAtgca/;
    }
    $tscfs{$tsc_id}.=$qscfs{$tt->[$i][3]};
    $Ns= 'n' x ($tt->[$i+1][0]-$tt->[$i][1]);
    $tscfs{$tsc_id}.=$Ns;
  }
  
  # finishing the last query scaffold
  $i=scalar(@$tt)-1;
  if($tt->[$i][2]<0){
    $qscfs{$tt->[$i][3]}=reverse $qscfs{$tt->[$i][3]};
    $qscfs{$tt->[$i][3]}=~ tr/ACGTacgt/TGCAtgca/;
  }
  $tscfs{$tsc_id}.=$qscfs{$tt->[$i][3]}; 
}

#### output
foreach my $id (sort(keys %tscfs)){
  print ">$id\n";
  while(length($tscfs{$id})>0){
    print substr($tscfs{$id},0,60,'')."\n";
  }
}

print STDERR "\n\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";
