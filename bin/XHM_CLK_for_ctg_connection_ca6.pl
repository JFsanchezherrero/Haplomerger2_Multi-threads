#!/usr/bin/perl  


## notes
## to break the misjoins
## for testing

## assumption
## 1. alignment是准确的
## 2. contig的准确的
## 3. 短的alignment加上mate pair是可靠的

#{CLK
#co1:7180001489473
#co2:7180001639009
#ori:I  #{N=normal=FF, A=anti-normal=RR, O=outie=RF, I=innie=FR}
#ovt:O  #O=overlap+mates, N=all mates
#ipc:0
#mea:-115.982
#std:1.737
#num:3  #link count
#sta:U #B C not count
#jls:
#GB3UXMY02HT82Ua,GB3UXMY02HT82Ub,M
#FPXQG6302TB8FMa,FPXQG6302TB8FMb,M
#}


# bug: 

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 2 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

   This script is to break the misjoins.
   
   Gzip files are supported.
   
   --clk=<clk file>
   --connection=<connection info file>
   
   
   Outputs include (placed in the same dir of the source file):
   *.new.connection.
    	 
USAGES


#### ################################
#### reading arguments

my ($in_clk);
if ($Arg_list =~ m/--clk=(\S+)/) {
	$in_clk=$1;
}else{
	die "no --clk arg.\n";
}

my ($in_connection);
if ($Arg_list =~ m/--connection=(\S+)/) {
	$in_connection=$1;
}else{
	die "no --connection arg.\n";
}

#### ################################
#### reading clk
my %clk; #$clk{ctg1@ctg2}=[ovt,num]; sta:B/C is not counted

{
	my $inFH_name=$in_clk;
	my $inFH;
	if($inFH_name=~m/\.gz$/){
		open($inFH,"gunzip -c $inFH_name |") or die "Can not open the file $inFH_name.\n";
	}else{
		open($inFH,"<$inFH_name") or die "Can not open the file $inFH_name.\n"; 
	}

	$/="}\n";
	while(<$inFH>){
		next if $_ !~ m/^.CLK\n/; #.=}
		if($_ =~ m/co1:(\S+)\nco2:(\S+)\nori:(\S+)\novt:(\S+)\n.+num:(\S+)\nsta:(\S+)\n/s){
			next if $6 eq 'B' or $6 eq 'C'; # could be Bad, Chimera, Assembly, Unknown
			my $ori=$3;
			my $ovt=$4;
			my $num=$5;
			#die "$1 e $2\n" if $1 eq $2; #testing#
			if(exists($clk{'ctg'.$1.'@ctg'.$2.'@'.$3})){
				$clk{'ctg'.$1.'@ctg'.$2.'@'.$3}->[0]= ($ovt eq "O" or $clk{'ctg'.$1.'@ctg'.$2.'@'.$3}->[0] eq 'O') ? 'O' : 'N';
				$clk{'ctg'.$1.'@ctg'.$2.'@'.$3}->[1]+=$num;
				#$clk{$1.'@'.$2}->[3].=$5;
			}else{
				$clk{'ctg'.$1.'@ctg'.$2.'@'.$3}=[$ovt,$num];
			}
			#testing# print "$1 $2 $ovt $num $5\n" if $5 ne 'U'; #testing# 
		}	
	}
	close $inFH;	
	$/="\n";	
}

print STDERR "Finished clk reading.\n";

#### ################################
#### reading misjoin

my $connection_head;
my @connection;

{
	my $inFH_name=$in_connection;
	my $inFH;
	if($inFH_name=~m/\.gz$/){
		open($inFH,"gunzip -c $inFH_name |") or die "Can not open the file $inFH_name.\n";
	}else{
		open($inFH,"<$inFH_name") or die "Can not open the file $inFH_name.\n"; 
	}
	
	while(<$inFH>){
		if(m/^#|^\s+/){ $connection_head.=$_; next; }
		chomp;
		my @tt =  split /\t/,$_;
		push @tt,('','');
		push @connection,\@tt;
	}
	close $inFH;
}

print STDERR "Finished connection reading.\n";

#### ################################
#### processing connection

{
	my ($ov,$num);
	my ($ctg1,$ctg2);
	my ($ctg1_sign,$ctg2_sign);
	my ($id1,$id2);
	foreach my $portion (@connection){
		next if $portion->[19]!=1;
		$ov='N';
		$num=0;
		($ctg1,$ctg2,$ctg1_sign,$ctg2_sign)=@{$portion}[5,12,10,17];
		if($ctg1_sign==1 and $ctg2_sign==1){
			$id1=$ctg1.'@'.$ctg2.'@N';
			$id2=$ctg2.'@'.$ctg1.'@R';
		}elsif($ctg1_sign==1 and $ctg2_sign==-1){
			$id1=$ctg1.'@'.$ctg2.'@I';
			$id2=$ctg2.'@'.$ctg1.'@I';			
		}elsif($ctg1_sign==-1 and $ctg2_sign==-1){
			$id1=$ctg1.'@'.$ctg2.'@R';
			$id2=$ctg2.'@'.$ctg1.'@N';			
		}elsif($ctg1_sign==-1 and $ctg2_sign==1){
			$id1=$ctg1.'@'.$ctg2.'@O';
			$id2=$ctg2.'@'.$ctg1.'@O';			
		}
		
		if(exists($clk{$id1})){
			$ov = ($clk{$id1}->[0] eq 'O' or $ov eq 'O') ? 'O' : 'N';
			$num+=$clk{$id1}->[1];
		}
		if(exists($clk{$id2})){
			$ov = ($clk{$id2}->[0] eq 'O' or $ov eq 'O') ? 'O' : 'N';
			$num+=$clk{$id2}->[1];
		}
		
		$portion->[29]=$ov;	
		$portion->[30]=$num;			
	}
	
}

#### ################################
#### outputing new connection

{
	my $outFH_name="hm.new_scaffolds_with_CLK";
	open(my $outFH,">$outFH_name") or die "can not create >$outFH_name.\n";
	
	print $outFH $connection_head;

	$"="\t"; #"
	foreach my $tt (@connection){
		print $outFH "@$tt\n";
	}
	$"=" "; #"
}
