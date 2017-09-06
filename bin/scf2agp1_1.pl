#!/usr/bin/perl  


## assumption
## 0. based on agp v1.1
## 1. contig is separated by N
## 2. contig is correct

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 0 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

   This script is to recreate an agp v1.1 for a given set of scaffolds. 
  
   --mono|ngap - def=mono; 
                 mono = one scaffold one contig;
                 ngap = a scaffold will be broken down into 
                 multiple contigs based on the Ngaps.
   
   cat scf.fa |./scf2agp1_1.pl >scf.agp
    	 
USAGES

# set the over-write flag, to over-write any existing files
my $method='--mono';
if ($Arg_list =~ m/--ngap/) {
	$method='--ngap';
}
print STDERR "Set to $method mode!\n";

#### ################################
#### reading scf

my %scfs;

{
	my ($line,$name,$seq,$is_the_end)=("","","",0);
	
	while($line=<STDIN>){
		last if $line =~ m/^>/;
	}	
	die "Problems with the fasta file; no symbol > found !\n" unless defined($line);
	
	while($is_the_end==0){
		
		($name,$seq)=('','');
		
		if($line =~ m/>\s*(\S+)/){
			$name=$1;
		}
	
		while($line=<STDIN>){
			last if $line =~ m/^>/;
			chomp $line;
			$seq .= $line;
		}
		$is_the_end=1 unless defined($line);
		
		$scfs{$name}=$seq;
	}
}

print STDERR "Finished scf fasta reading.\n";
print "## Created by scf2agp1_1.pl; method $method .\n";

#### ################################
#### create agp 1.1

## note that in the ca6 ctg.fasta, ctg name includes a leading 'ctg', but not in the agp. 
# for ctg:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:W, 5:ctg_id, 6:ctg_start, 7:ctg_end, 8:strand]
# for gap:: $agp{scf_id}->[n][0:scf_id, 1:start, 2:end, 3:count, 4:N, 5:gap_size, 6:'fragment', 7:'yes', 8:null]
# need to convert from 1-base to 0 base
my %agp; 

if($method eq '--mono'){
	# one scaffold one contig
  my ($count,$ctg_id,$clen);
  my ($start,$end); #0-base
  
  foreach my $sc_id (keys %scfs){
    my $seq=$scfs{$sc_id};
    $count=1;
    $clen=length($seq);
    ($start,$end)=(1,$clen);
    $ctg_id=$sc_id.'_ctg'.$count;
    print "$sc_id\t$start\t$end\t$count\tW\t$ctg_id\t1\t$clen\t+\n";
  }  	
}else{
	# breaking scaffolds into contigs based on the Ngaps
  my ($count,$ctg_id,$clen);
  my ($start,$end); #0-base
  my ($nlen,$nstart,$nend);
  
  foreach my $sc_id (keys %scfs){
    my $seq=$scfs{$sc_id};
    $count=0;
    ($start,$end)=(1,length($seq));
    ($nstart,$nend,$nlen)=($end+1,$end+1,0);
    while($seq=~/(N+)/gi){
      ($nstart,$nend,$nlen)=($-[1]+1,$+[1],$+[1]-$-[1]);
      $clen=$nstart-$start;
      if($clen>0){ $count++; $ctg_id=$sc_id.'_ctg'.$count; print "$sc_id\t$start\t$-[1]\t$count\tW\t$ctg_id\t1\t$clen\t+\n"; }
      $count++;
      print "$sc_id\t$nstart\t$nend\t$count\tN\t$nlen\tfragment\tyes\n";
      $start=$nend+1;
      ($nstart,$nend,$nlen)=($end+1,$end+1,0);
    }
    if($start<=$end){ $clen=$end-$start+1; $count++; $ctg_id=$sc_id.'_ctg'.$count; print "$sc_id\t$start\t$end\t$count\tW\t$ctg_id\t1\t$clen\t+\n"; }
  }
}	

print STDERR "Finished agp outputing.\n";

