#!/usr/bin/perl  


## notes
## The file's coordinates are 0-base; the cabog6.1 uses agp specification v1.1.
## Only gap-closing contigs are extract, ans stored in low cases.

## NOTE the gap order is agree with the original scf only if the original scf don't have leading N at each sequences.
## NOTE require the gapcloser fasta and fill file are intact!!!!

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 1 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

   This script is to extract the contig sequences that can
   fully close a gap from the gapcloser results (gc v1.12).
   
   Gzip files are supported.
   
   --fasta=N - the filled scf
   --info=N - the fill file
   
   NOTE that fill file allows negative gaps!
   
   Outputs include (#scf_id, gap_id, ctg_seq): STDOUT.
    	 
USAGES


print STDERR "\n\n========== Start at " . qx(date);
print STDERR "$0 $Arg_list\n\n";
my $timer = time();

#### ################################
#### reading arguments

my ($fasta_name);
if ($Arg_list =~ m/--fasta=(\S+)/) {
	$fasta_name=$1;
}else{
	die "no --fasta arg.\n";
}

my ($info_name); 
if ($Arg_list =~ m/--info=(\S+)/) {
	$info_name=$1;
}else{
	die "no --info arg.\n";
}

#my @pairs;
#{
#	open(my $inFH,"<$ARGV[0]") or die "can not open $ARGV[0].\n";
#	while(<$inFH>){
#		next if /^#|^\s+/;
#		chomp;
#		push @pairs,[split /\s+/,$_];
#	}
#	close $inFH;
#}

my %info;
my %info_list;
my %scfs;

#### ################################
#### reading fasta  

  {
  	undef %scfs;
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
  		
  		$scfs{$name}=uc $seq;
  		#print STDERR "read scaffold $name\n";
  	}
  	close $inFH;	
  }
  
#### ################################
#### reading gapcloser's fill file (require intact fill file!!!!) 
  
  {
  
  	my $inFH;
  	if($info_name=~m/\.gz$/){
  		open($inFH,"gunzip -c $info_name |") or die "Can not open gunzip -c $info_name |.\n";
  	}else{
  		open($inFH,"<$info_name") or die "Can not open $info_name.\n"; 
  	}
  		
  	my ($curr_scf,$count)=(undef,0);
  	
  	while(<$inFH>){
  		chomp;
  		if(m/^\>(\S+)/){
  			$curr_scf=$1; $count=0;
  			next;
  		}
  		my @tt=split /\t/,$_;
  		$count++;
  		next if $tt[4]<1;
  		die "can not find the corresponding scf.\n" unless exists $scfs{$curr_scf};
  		my $seq;
  		if($tt[1]-$tt[0]>0){
  			$seq = substr($scfs{$curr_scf},$tt[0],$tt[1]-$tt[0]);
  			die "can not get the filling seq.\n" unless length($seq)==$tt[1]-$tt[0];
  		}elsif($tt[1]-$tt[0]==0){
  		  $seq='';  			
  		}else{
  			$seq=$tt[1]-$tt[0]; # negative value
  		}
  		
  		if(!exists($info_list{$curr_scf.'@'.$count})){
  			$info_list{$curr_scf.'@'.$count}=1;
  			push @{$info{$curr_scf}},[ $curr_scf, $count, $seq ];
  		}
  	}
  	close $inFH;
  }


#### ################################
#### output info file

{
	#my $outFH;
	#open($outFH,"| gzip -c >$gc_fctgs_name") or die "can not open | gzip -c >$gc_fctgs_name.\n";
	print "#0-based.\n#0scf_id\tgap_order\tfilled_seq\n";
	foreach my $scf_id (keys %info){
		my $tt=$info{$scf_id};
		for(my $i=0;$i<scalar(@$tt);$i++){
			print "$tt->[$i][0]\t$tt->[$i][1]\t$tt->[$i][2]\n";
		}
	}
	#close $outFH;
}

print STDERR "\n\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";
